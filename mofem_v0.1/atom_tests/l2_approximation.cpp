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
  ierr = mField.add_field("H1FIELD_SCALAR_L2",L2,1); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("ELEM_L2_SCALAR"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_L2_SCALAR","H1FIELD_SCALAR_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_L2_SCALAR","H1FIELD_SCALAR_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_L2_SCALAR","H1FIELD_SCALAR_L2"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("PROBLEM_SCALAR_L2"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("PROBLEM_SCALAR_L2","ELEM_L2_SCALAR"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //add ents to field and set app. order
  ierr = mField.add_ents_to_field_by_TETs(0,"H1FIELD_SCALAR_L2"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"H1FIELD_SCALAR_L2",2); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELEM_L2_SCALAR",MBTET); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_SCALAR_L2",bit_level0); CHKERRQ(ierr);

  //Only for testing create independet MoFEM interface which works on the same mesh (MOAB) database
  FieldCore core2(moab);
  FieldInterface& mField2 = core;

  //build fields
  ierr = mField2.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField2.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField2.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField2.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField2.simple_partition_problem("PROBLEM_SCALAR_L2"); CHKERRQ(ierr);
  ierr = mField2.partition_finite_elements("PROBLEM_SCALAR_L2"); CHKERRQ(ierr);
  ierr = mField2.partition_ghost_dofs("PROBLEM_SCALAR_L2"); CHKERRQ(ierr);

  struct ApproxAnaliticalFunction {

    FieldData scalar(ublas::vector<FieldData> coords) {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
    }

  };

  struct ScalarApprox: public FEMethod_UpLevelStudent,ApproxAnaliticalFunction {

    ScalarApprox(Interface& _moab): FEMethod_UpLevelStudent(_moab) {}; 

    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      vector<ublas::matrix<FieldData> > RowN;
      vector<ublas::matrix<FieldData> > ColN;
      ierr = GetGaussRowNMatrix("H1FIELD_SCALAR_L2",MBTET,RowN); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("H1FIELD_SCALAR_L2",MBTET,ColN); CHKERRQ(ierr);
      ublas::matrix<FieldData> NTN;
      NTN = ublas::zero_matrix<FieldData>(RowN[0].size2(),ColN[0].size2());
      ublas::vector<FieldData> F;
      F = ublas::zero_vector<FieldData>(RowN[0].size2());
      for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) {
	double w = G_TET_W45[gg];
	assert(w == w);
	NTN += w*V*prod(trans(RowN[gg]), ColN[gg]);
	double val = scalar(coords_at_Gauss_nodes[gg]);
	F += V*w*val*ublas::matrix_row<ublas::matrix<FieldData> >(RowN[gg],0);
      }
      ublas::matrix<FieldData> L(NTN.size1(),NTN.size2());
      ublas::vector<FieldData> x = F;
      cholesky_decompose(NTN,L);
      cholesky_solve(L,x,ublas::lower());

      for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"H1FIELD_SCALAR_L2",MBTET,dof)) {
	DofIdx idx = dof->get_EntDofIdx();
	dof->get_FieldData() = x[idx];
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }


  };

  struct ScalarApprox_Check: public FEMethod_UpLevelStudent,ApproxAnaliticalFunction {

    ofstream myfile;
    ScalarApprox_Check(Interface& _moab): FEMethod_UpLevelStudent(_moab) {
      //Open mesh_file_name.txt for writing
      myfile.open("l2_approximation.txt");
    }; 

    ~ScalarApprox_Check() {
      myfile.close();
    }

    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(1*4);
      ShapeMBTET(&g_NTET[0],G_TET_X1,G_TET_Y1,G_TET_Z1,1);
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      vector<ublas::vector<FieldData> > field_value_in_element_centre;
      ierr = GetGaussDataVector("H1FIELD_SCALAR_L2",field_value_in_element_centre); CHKERRQ(ierr);

      double val = scalar(coords_at_Gauss_nodes[0]);
      cout << boost::format("%.3lf") % roundn(val - (field_value_in_element_centre[0])[0] );
      cout << " " << boost::format("%.3lf") % roundn( (field_value_in_element_centre[0])[0] );
      cout << endl;

      myfile << boost::format("%.3lf") % roundn(val - (field_value_in_element_centre[0])[0] );
      myfile << " " << boost::format("%.3lf") % roundn( (field_value_in_element_centre[0])[0] );
      myfile << endl;


      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }


  };

  ScalarApprox fe(moab);
  ierr = mField2.loop_finite_elements("PROBLEM_SCALAR_L2","ELEM_L2_SCALAR",fe);  CHKERRQ(ierr);

  ScalarApprox_Check fe_check(moab);
  ierr = mField2.loop_finite_elements("PROBLEM_SCALAR_L2","ELEM_L2_SCALAR",fe_check);  CHKERRQ(ierr);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }


}
