/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include <petscsys.h> 
#include <petsctime.h>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//Rounding
#define RND_EPS 1e-6
double roundn(double n) {

    //break n into fractional part (fract) and integral part (intp)
    double fract, intp;
    fract = modf(n,&intp);
    
    // case where n approximates zero, set n to "positive" zero
    if (abs(intp)==0) {
      if(abs(fract)<=RND_EPS) {
	n=0.000;
      }
    }

    return n;
}


int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("FIELD_HDIV",Hdiv,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("ELEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("PROBLEM","ELEM_HDIV"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //add ents to field and set app. order
  ierr = mField.add_ents_to_field_by_TETs(0,"FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"FIELD_HDIV",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"FIELD_HDIV",5); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELEM_HDIV",MBTET); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.simple_partition_problem("PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROBLEM"); CHKERRQ(ierr);

  struct ApproxAnaliticalFunction {

    FieldData scalar(ublas::vector<FieldData> coords) {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
    }
  
    ublas::vector<FieldData> dx_field;
    ublas::vector<FieldData>& dx_scalar(ublas::vector<FieldData> coords) {
      dx_field.resize(3);
      dx_field[0] = 2*coords[0];
      dx_field[1] = 2*coords[1];
      dx_field[2] = 2*coords[2];
      return dx_field;
    }

  };

  struct HdivApprox: public FEMethod_UpLevelStudent,ApproxAnaliticalFunction {

    HdivApprox(Interface& _moab): FEMethod_UpLevelStudent(_moab) {}; 

    vector<double> g_NTET;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      fe_ent_ptr = fe_ptr->fe_ptr;
      ierr = InitDataStructures(); CHKERRQ(ierr);
      ierr = GlobIndices(); CHKERRQ(ierr);
      ierr = LocalIndices(); CHKERRQ(ierr);
      ierr = DataOp(); CHKERRQ(ierr);
      ierr = ShapeFunctions_TET(g_NTET); CHKERRQ(ierr);
      ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
      ierr = GetRowNMatrix_at_GaussPoint(); CHKERRQ(ierr);
      ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);



      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  };

  HdivApprox fe(moab);
  ierr = mField.loop_finite_elements("PROBLEM","ELEM_HDIV",fe);  CHKERRQ(ierr);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }


}
