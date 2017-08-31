/** \file convective_matrix.cpp

 * \ingroup convective_mass_elem
 * \ingroup nonlinear_elastic_elem
 *
 * Atom test for convective mass element
 *
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

 #include <BasicFiniteElements.hpp>
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
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
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

    //ref meshset ref level 0
    ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
    ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
    ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);

    //FE
    ierr = m_field.add_finite_element("ELASTIC"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);

    //define problems
    ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problems
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

    //set refinement level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = m_field.add_ents_to_field_by_type(0,MBTET,"SPATIAL_POSITION"); CHKERRQ(ierr);

    //add finite elements entities
    ierr = m_field.add_ents_to_finite_element_by_bit_ref(bit_level0,BitRefLevel().set(),"ELASTIC",MBTET); CHKERRQ(ierr);

    //set app. order
    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      order = 1;
    }
    ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

    ierr = m_field.add_finite_element("NEUAMNN_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
      ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI,"NEUAMNN_FE"); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
      ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI,"NEUAMNN_FE"); CHKERRQ(ierr);
    }

    //Velocity
    ierr = m_field.add_field("SPATIAL_VELOCITY",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_type(0,MBTET,"SPATIAL_VELOCITY"); CHKERRQ(ierr);
    int order_velocity = 1;
    ierr = m_field.set_field_order(0,MBTET,"SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_VELOCITY",1); CHKERRQ(ierr);

    ierr = m_field.add_field("DOT_SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_POSITION",1); CHKERRQ(ierr);
    ierr = m_field.add_field("DOT_SPATIAL_VELOCITY",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_VELOCITY",order_velocity); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_VELOCITY",1); CHKERRQ(ierr);

    ConvectiveMassElement inertia(m_field,1);
    ierr = inertia.setBlocks(); CHKERRQ(ierr);
    ierr = inertia.addConvectiveMassElement("MASS_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = inertia.addVelocityElement("VELOCITY_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","MASS_ELEMENT"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","VELOCITY_ELEMENT"); CHKERRQ(ierr);

    //build field
    ierr = m_field.build_fields(); CHKERRQ(ierr);

    double scale_positions = 2;
    {
      EntityHandle node = 0;
      double coords[3];
      for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"SPATIAL_POSITION",dof_ptr)) {
        if(dof_ptr->get()->getEntType()!=MBVERTEX) continue;
        EntityHandle ent = dof_ptr->get()->getEnt();
        int dof_rank = dof_ptr->get()->getDofCoeffIdx();
        double &fval = dof_ptr->get()->getFieldData();
        if(node!=ent) {
          rval = moab.get_coords(&ent,1,coords); CHKERRQ_MOAB(rval);
          node = ent;
        }
        fval = scale_positions*coords[dof_rank];
      }
    }

    double scale_velocities = 4;
    {
      EntityHandle node = 0;
      double coords[3];
      for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"DOT_SPATIAL_POSITION",dof_ptr)) {
        if(dof_ptr->get()->getEntType()!=MBVERTEX) continue;
        EntityHandle ent = dof_ptr->get()->getEnt();
        int dof_rank = dof_ptr->get()->getDofCoeffIdx();
        double &fval = dof_ptr->get()->getFieldData();
        if(node!=ent) {
          rval = moab.get_coords(&ent,1,coords); CHKERRQ_MOAB(rval);
          node = ent;
        }
        fval = scale_velocities*coords[dof_rank];
      }
    }

    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

    //build adjacencies
    ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

    //build problem

    ProblemsManager *prb_mng_ptr;
    ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);

    ierr = prb_mng_ptr->buildProblem("ELASTIC_MECHANICS",true); CHKERRQ(ierr);

    //partition
    ierr = prb_mng_ptr->partitionProblem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionGhostDofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //create matrices
    Vec F;
    ierr = m_field.query_interface<VecManager>()->vecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
    Vec D;
    ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
    Mat Aij;
    ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

    ierr = inertia.setConvectiveMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = inertia.setVelocityOperators("SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);

    inertia.getLoopFeMassRhs().ts_F = F;
    inertia.getLoopFeMassRhs().ts_a = 1;
    inertia.getLoopFeMassLhs().ts_B = Aij;
    inertia.getLoopFeMassLhs().ts_a = 1;

    inertia.getLoopFeVelRhs().ts_F = F;
    inertia.getLoopFeVelRhs().ts_a = 1;
    inertia.getLoopFeVelLhs().ts_B = Aij;
    inertia.getLoopFeVelLhs().ts_a = 1;

    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","MASS_ELEMENT",inertia.getLoopFeMassRhs()); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","VELOCITY_ELEMENT",inertia.getLoopFeVelRhs()); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","MASS_ELEMENT",inertia.getLoopFeMassLhs()); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","VELOCITY_ELEMENT",inertia.getLoopFeVelLhs()); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    // PetscViewer viewer;
    // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"convective_matrix.txt",&viewer); CHKERRQ(ierr);
    //ierr = PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_SYMMODU); CHKERRQ(ierr);

    //ierr = VecChop(F,1e-4); CHKERRQ(ierr);
    // ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // ierr = VecView(F,viewer); CHKERRQ(ierr);

    //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
    // MatChop(Aij,1e-4);
    // MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
    // MatView(Aij,viewer);
    //std::string wait;
    //std::cin >> wait;
    //
    //  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

    double sum = 0;
    ierr = VecSum(F,&sum); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %9.8e\n",sum); CHKERRQ(ierr);
    double fnorm;
    ierr = VecNorm(F,NORM_2,&fnorm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm  = %9.8e\n",fnorm); CHKERRQ(ierr);

    double mnorm;
    ierr = MatNorm(Aij,NORM_1,&mnorm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"mnorm  = %9.8e\n",mnorm); CHKERRQ(ierr);


    if(fabs(sum-6.27285463e+00)>1e-8) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
    }
    if(fabs(fnorm-1.28223353e+00)>1e-6) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
    }
    if(fabs(mnorm-1.31250000e+00)>1e-6) {
      SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
    }




    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = MatDestroy(&Aij); CHKERRQ(ierr);



  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

  return 0;


}
