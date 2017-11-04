/** \file build_problems.cpp

  \brief Atom test for building problems

  \bug Not verifying what if partitioned mesh is loaded.

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
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {




  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

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

    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    const char *option;
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

    //Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("F1",L2,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
    ierr = m_field.add_field("F2",HDIV,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

    //meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    //add entities to field
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"F1"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"F2"); CHKERRQ(ierr);

    //set app. order
    //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
    int order = 2;
    ierr = m_field.set_field_order(root_set,MBTET,"F1",order-1); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBTET,"F2",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(root_set,MBTRI,"F2",order); CHKERRQ(ierr);

    ierr = m_field.build_fields(); CHKERRQ(ierr);

    //add elements
    ierr = m_field.add_finite_element("E1"); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("E2"); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("E3"); CHKERRQ(ierr);

    ierr = m_field.modify_finite_element_add_field_row("E1","F1"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("E1","F1"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("E2","F2"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("E2","F2"); CHKERRQ(ierr);
    //To build composite problem
    ierr = m_field.modify_finite_element_add_field_row("E3","F1"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("E3","F2"); CHKERRQ(ierr);

    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"E1"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"E2"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"E3"); CHKERRQ(ierr);

    ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
    ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

    //Problems
    ierr = m_field.add_problem("P1"); CHKERRQ(ierr);
    ierr = m_field.add_problem("P2"); CHKERRQ(ierr);

    //set refinement level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("P1",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_add_bit("P2",bit_level0); CHKERRQ(ierr);

    ierr = m_field.modify_problem_add_finite_element("P1","E1"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("P2","E2"); CHKERRQ(ierr);

    //build problems
    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRQ(ierr);
    ierr = prb_mng_ptr->buildProblem("P1",true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->buildProblem("P2",true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionProblem("P1"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionProblem("P2"); CHKERRQ(ierr);

    ierr = prb_mng_ptr->partitionFiniteElements("P1"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("P1"); CHKERRQ(ierr);
    ierr = m_field.partition_check_matrix_fill_in("P1",-1,-1,0); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("P2"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("P2"); CHKERRQ(ierr);
    ierr = m_field.partition_check_matrix_fill_in("P2",-1,-1,0); CHKERRQ(ierr);

    //compose problem
    ierr = m_field.add_problem("P3"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_add_bit("P3",bit_level0); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("P3","E3"); CHKERRQ(ierr);

    ierr = prb_mng_ptr->buildProblem("P3",false); CHKERRQ(ierr);
    ierr = prb_mng_ptr->inheritPartition("P3","P1",false,"P2",true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("P3"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("P3"); CHKERRQ(ierr);
    ierr = m_field.partition_check_matrix_fill_in("P3",-1,-1,0); CHKERRQ(ierr);

    /*Mat m;
    ierr = m_field.MatCreateMPIAIJWithArrays("P3",&m); CHKERRQ(ierr);
    PetscViewer viewer;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"build_composite_problem.txt",&viewer); CHKERRQ(ierr);
    MatView(m,viewer);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = MatDestroy(&m); CHKERRQ(ierr);*/

    

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  //finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
