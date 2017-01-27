/** \file build_composite_problems.cpp

  \brief Atom test for building composite problems

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

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

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

    //read mesh and create moab and mofem datastrutures

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    EntityHandle root_set = moab.get_root_set();
    Range tets;
    rval = moab.get_entities_by_type(root_set,MBTET,tets,false); CHKERRQ(rval);

    ProblemsManager *prb_mng_ptr;
    ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionMesh(tets,3,2,m_field.getCommSize()); CHKERRQ(ierr);

    Tag part_tag = pcomm->part_tag();
    Range tagged_sets;
    rval = m_field.get_moab().get_entities_by_type_and_tag(
      0,MBENTITYSET,&part_tag,NULL,1,tagged_sets,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    for(Range::iterator mit = tagged_sets.begin();mit!=tagged_sets.end();mit++) {
      int part;
      rval = moab.tag_get_data(part_tag,&*mit,1,&part); CHKERRQ_MOAB(rval);
      if(part==m_field.getCommRank()) {
        pcomm->partition_sets().insert(*mit);
      }
    }
    rval = pcomm->resolve_shared_ents(0,3); CHKERRQ_MOAB(rval);

    Range owned_tets;
    rval = pcomm->get_part_entities(owned_tets,3); MB_CHK_ERR(rval);

    std::ostringstream file_owned;
    file_owned << "out_owned_" << m_field.getCommRank() << ".vtk";
    EntityHandle meshset_owned;
    rval = moab.create_meshset(MESHSET_SET,meshset_owned); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(meshset_owned,owned_tets); CHKERRQ_MOAB(rval);
    rval = moab.write_file(file_owned.str().c_str(),"VTK","",&meshset_owned,1); CHKERRQ_MOAB(rval);

    Range shared_ents;
    // Get entities shared with all other processors
    rval = pcomm->get_shared_entities(-1,shared_ents);MB_CHK_ERR(rval);

    std::ostringstream file_shared_owned;
    file_shared_owned << "out_shared_owned_" << m_field.getCommRank() << ".vtk";
    EntityHandle meshset_shared_owned;
    rval = moab.create_meshset(MESHSET_SET,meshset_shared_owned); CHKERRQ_MOAB(rval);
    rval = moab.add_entities(meshset_shared_owned,shared_ents); CHKERRQ_MOAB(rval);
    rval = moab.write_file(file_shared_owned.str().c_str(),"VTK","",&meshset_shared_owned,1); CHKERRQ_MOAB(rval);

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle part_set = pcomm->partition_sets().front();
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("F1",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
    ierr = m_field.add_field("F2",H1,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);
    //add entities to field
    ierr = m_field.add_ents_to_field_by_TETs(part_set,"F1"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_TETs(part_set,"F2"); CHKERRQ(ierr);
    int order = 4;
    ierr = m_field.set_field_order(part_set,MBTET,"F1",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBTRI,"F1",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBEDGE,"F1",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBVERTEX,"F1",1); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBTET,"F2",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBTRI,"F2",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBEDGE,"F2",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(part_set,MBVERTEX,"F2",1); CHKERRQ(ierr);
    ierr = m_field.build_fields(); CHKERRQ(ierr);

    //Elements
    ierr = m_field.add_finite_element("E1"); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("E2"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("E1","F1"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("E1","F1"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("E2","F2"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("E2","F2"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TETs(part_set,"E1"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TETs(part_set,"E2"); CHKERRQ(ierr);
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

    //Build problems
    ierr = prb_mng_ptr->buildProblemOnDistributedMesh("P1",true,1); CHKERRQ(ierr);
    ierr = prb_mng_ptr->buildProblemOnDistributedMesh("P2",true,1); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("P1",true,0,m_field.getCommSize(),1); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("P1",1); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("P2",true,0,m_field.getCommSize(),1); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("P2",1); CHKERRQ(ierr);

    if(0) {
      Mat m;
      ierr = m_field.MatCreateMPIAIJWithArrays("P1",&m); CHKERRQ(ierr);
      MatView(m,PETSC_VIEWER_DRAW_WORLD);
      std::string wait;
      std::cin >> wait;
      ierr = MatDestroy(&m); CHKERRQ(ierr);
    }


    ierr = m_field.partition_check_matrix_fill_in("P1",-1,-1,1); CHKERRQ(ierr);
    ierr = m_field.partition_check_matrix_fill_in("P2",-1,-1,1); CHKERRQ(ierr);



  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  // Finish work cleaning memory, getting statistics, etc.
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
