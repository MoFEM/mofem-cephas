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

#include <MoFEM.hpp>
#include <Projection10NodeCoordsOnField.hpp>

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

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
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

  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("TET_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add problems 
  //ierr = mField.add_problem("EDGE_PROJECTOR_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.add_problem("TET_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("TET_PROBLEM","TET_ELEM"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,true); CHKERR_PETSC(rval);
  Range edges;
  rval = moab.get_entities_by_type(0,MBEDGE,edges,true); CHKERR_PETSC(rval);
  ierr = mField.seed_finite_elements(edges); CHKERRQ(ierr);

  //add ents to field and set app. order

  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_by_TETs(tets,"TET_ELEM"); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("TET_PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("TET_PROBLEM"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("TET_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("TET_PROBLEM"); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);

    //Open mesh_file_name.txt for writing
    ofstream myfile;
    myfile.open ("10node_sphere.txt");
    
    //Output displacements
    cout << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    myfile << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr))
    {
        //if(dof_ptr->get_ent_type()!=MBEDGE) continue;
        
        if(dof_ptr->get_dof_rank()==0)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData();
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==1)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData();
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==2)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData();
            cout << boost::format("%.3lf") % roundn(fval) << endl;
            myfile << boost::format("%.3lf") % roundn(fval) << endl;
        }
        
    }


  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFinalize();

}
