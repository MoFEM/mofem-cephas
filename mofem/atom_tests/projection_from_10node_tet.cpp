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

using namespace MoFEM;

static char help[] = "...\n\n";

//Rounding
#define RND_EPS 1e-6
double roundn(double n) {

    //break n into fractional part (fract) and integral part (intp)
    double fract, intp;
    fract = modf(n,&intp);

    // case where n approximates zero, set n to "positive" zero
    if (std::abs(intp)==0) {
      if(std::abs(fract)<=RND_EPS) {
	n=0.000;
      }
    }

    return n;
}


int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERRQ_MOAB(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //add filds
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);

  //add finite elements
  ierr = m_field.add_finite_element("TET_ELEM"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TET_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add problems
  //ierr = m_field.add_problem("EDGE_PROJECTOR_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.add_problem("TET_PROBLEM"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = m_field.modify_problem_add_finite_element("TET_PROBLEM","TET_ELEM"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,true); CHKERRQ_MOAB(rval);
  Range edges;
  rval = moab.get_entities_by_type(0,MBEDGE,edges,true); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_finite_elements(edges); CHKERRQ(ierr);

  //add ents to field and set app. order

  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_by_type(tets,MBTET,"TET_ELEM"); CHKERRQ(ierr);

  //set problem level
  ierr = m_field.modify_problem_ref_level_add_bit("TET_PROBLEM",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->buildProblem("TET_PROBLEM",false); CHKERRQ(ierr);

  //partition
  ierr = prb_mng_ptr->partitionProblem("TET_PROBLEM"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TET_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TET_PROBLEM"); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method); CHKERRQ(ierr);

    //Open mesh_file_name.txt for writing
    std::ofstream myfile;
    myfile.open ("10node_sphere.txt");

    //Output displacements
    std::cout << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>" << std::endl;
    myfile << "<<<< Dofs (X-Translation, Y-Translation, Z-Translation) >>>>>" << std::endl;

    const DofEntity_multiIndex *dofs_ptr;
    ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
    DofEntity_multiIndex_uid_view dofs_view;
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",dof_ptr)) {
      dofs_view.insert(*dof_ptr);
    }
    for(
      DofEntity_multiIndex_uid_view::iterator
      dit=dofs_view.begin();dit!=dofs_view.end();dit++
    ) {
        //if(dof_ptr->getEntType()!=MBEDGE) continue;

        if((*dit)->getDofCoeffIdx()==0)
        {
            //Round and truncate to 3 decimal places
            double fval = (*dit)->getFieldData();
            std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if((*dit)->getDofCoeffIdx()==1)
        {
            //Round and truncate to 3 decimal places
            double fval = (*dit)->getFieldData();
            std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if((*dit)->getDofCoeffIdx()==2)
        {
            //Round and truncate to 3 decimal places
            double fval = (*dit)->getFieldData();
            std::cout << boost::format("%.3lf") % roundn(fval) << std::endl;
            myfile << boost::format("%.3lf") % roundn(fval) << std::endl;
        }

    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  } catch (const std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << std::endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFinalize();

}
