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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm)
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  BARRIER_RANK_END(pcomm)

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("DISP",H1,3); CHKERRQ(ierr);
  ierr = m_field.add_field("TEMP",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("PROB"); CHKERRQ(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("PROB",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"DISP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order_temp = 2;
  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order_temp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order_temp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order_temp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  int order_disp = 3;
  ierr = m_field.set_field_order(root_set,MBTET,"DISP",order_disp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"DISP",order_disp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"DISP",order_disp); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"DISP",1); CHKERRQ(ierr);

  ThermalStressElement thermal_stress_elem(m_field);
  ierr = thermal_stress_elem.addThermalSterssElement("ELAS","DISP","TEMP"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("PROB","ELAS"); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("PROB"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("PROB"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("PROB"); CHKERRQ(ierr);

  //set temerature at nodes
  for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(m_field,"TEMP",MBVERTEX,dof)) {
    EntityHandle ent = dof->get()->getEnt();
    ublas::vector<double> coords(3);
    rval = moab.get_coords(&ent,1,&coords[0]); CHKERRQ_MOAB(rval);
    dof->get()->getFieldData() = 1;
  }

  Vec F;
  ierr = m_field.VecCreateGhost("PROB",ROW,&F); CHKERRQ(ierr);
  ierr = thermal_stress_elem.setThermalStressRhsOperators("DISP","TEMP",F,1); CHKERRQ(ierr);

  ierr = m_field.loop_finite_elements("PROB","ELAS",thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  // PetscViewer viewer;
  // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"forces_and_sources_thermal_stress_elem.txt",&viewer); CHKERRQ(ierr);
  // ierr = VecChop(F,1e-4); CHKERRQ(ierr);
  // ierr = VecView(F,viewer); CHKERRQ(ierr);
  // ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  double sum = 0;
  ierr = VecSum(F,&sum); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %9.8f\n",sum); CHKERRQ(ierr);
  double fnorm;
  ierr = VecNorm(F,NORM_2,&fnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm  = %9.8e\n",fnorm); CHKERRQ(ierr);
  if(fabs(sum)>1e-7) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
  }
  if(fabs(fnorm-2.64638118e+00)>1e-7) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
  }


  ierr = VecZeroEntries(F); CHKERRQ(ierr);

  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
