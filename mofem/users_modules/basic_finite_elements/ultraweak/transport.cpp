/** \fi;e transport.cpp
\brief Example implementation of transport problem using ultra-week formulation

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

  //Create mofem interface
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entities bit level
  ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FLUXES",HDIV,1); CHKERRQ(ierr);
  ierr = m_field.add_field("VALUES",L2,1); CHKERRQ(ierr);
  ierr = m_field.add_field("ERROR",L2,1); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FLUXES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"VALUES"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"ERROR"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }

  ierr = m_field.set_field_order(root_set,MBTET,"FLUXES",order+1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FLUXES",order+1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"VALUES",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"ERROR",0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //finite elements

  /** thefine sources and other stuff
    *
    * UltraWeakTransportElement is a class collecting functons, opertors and
    * data for ultra week implementation of transport element. See there to
    * learn how elements are created or how operators look like.
    *
    * Some methods in UltraWeakTransportElement are abstract, f.e. user need to
    * implement own surce therm.
    *
    */
  struct MyUltraWeakFE: public UltraWeakTransportElement {

    MyUltraWeakFE(MoFEM::Interface &m_field): UltraWeakTransportElement(m_field) {};

    PetscErrorCode getFlux(EntityHandle ent,const double x,const double y,const double z,double &flux) {
      PetscFunctionBegin;
      //double d = sqrt(x*x+y*y+z*z);
      flux = 1;//-pow(d,5./4.);
      PetscFunctionReturn(0);
    }

    PetscErrorCode getBcOnValues(
      const EntityHandle ent,
      const double x,const double y,const double z,
      double &value) {
      PetscFunctionBegin;
      value = 1;
      PetscFunctionReturn(0);
    }

    PetscErrorCode getBcOnFluxes(
      const EntityHandle ent,
      const double x,const double y,const double z,
      double &flux) {
      PetscFunctionBegin;
      flux = 0.;
      PetscFunctionReturn(0);
    }


  };

  MyUltraWeakFE ufe(m_field);
  ierr = ufe.addFiniteElements("FLUXES","VALUES","ERROR"); CHKERRQ(ierr);

  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(BitRefLevel().set(0),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
  Skinner skin(&moab);
  Range skin_faces; // skin faces from 3d ents
  rval = skin.find_skin(0,tets,false,skin_faces); CHKERR_MOAB(rval);

  // note: what is essential (dirichlet) is natural (neumann) for ultra wik comparic to classical FE
  /*Range neumann_tris;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|TEMPERATURESET,it)) {

    Range tris;
    ierr = it->getMeshSetIdEntitiesByDimension(m_field.get_moab(),2,tris,true); CHKERRQ(ierr);
    neumann_tris.insert(tris.begin(),tris.end());

  }
  Range dirichel_tris = subtract(skin_faces,neumann_tris);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(dirichel_tris,"ULTRAWEAK_FLUXDIRICHLET"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(neumann_tris,"ULTRAWEAK_FLUXNEUMANN"); CHKERRQ(ierr);*/
  ierr = m_field.add_ents_to_finite_element_by_TRIs(skin_faces,"ULTRAWEAK_FLUXNEUMANN"); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(BitRefLevel().set(0)); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("ULTRAWEAK"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ULTRAWEAK",BitRefLevel().set(0)); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK"); CHKERRQ(ierr);

  //boundary conditions
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_FLUXDIRICHLET"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK","ULTRAWEAK_FLUXNEUMANN"); CHKERRQ(ierr);

  ierr = m_field.add_problem("ULTRAWEAK_CALCULATE_ERROR"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ULTRAWEAK_CALCULATE_ERROR",BitRefLevel().set(0)); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ULTRAWEAK_CALCULATE_ERROR","ULTRAWEAK_ERROR"); CHKERRQ(ierr);

  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("ULTRAWEAK"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ULTRAWEAK"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("ULTRAWEAK"); CHKERRQ(ierr);

  //partition for problem calculating error
  ierr = m_field.partition_simple_problem("ULTRAWEAK_CALCULATE_ERROR"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ULTRAWEAK_CALCULATE_ERROR"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("ULTRAWEAK_CALCULATE_ERROR"); CHKERRQ(ierr);

  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ULTRAWEAK",&Aij); CHKERRQ(ierr);
  Vec F,D,D0;
  ierr = m_field.VecCreateGhost("ULTRAWEAK",COL,&D); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ULTRAWEAK",COL,&D0); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ULTRAWEAK",ROW,&F); CHKERRQ(ierr);

  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(D0); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(D); CHKERRQ(ierr);

  ufe.feTriFluxValue.getOpPtrVector().clear();
  ufe.feTriFluxValue.getOpPtrVector().push_back(new MyUltraWeakFE::OpEvaluateBcOnFluxes(ufe,"FLUXES",D0));
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_FLUXDIRICHLET",ufe.feTriFluxValue); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(D0); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(D0); CHKERRQ(ierr);

  ierr = m_field.set_global_ghost_vector("ULTRAWEAK",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  IS dirchlet_ids;
  ierr = ufe.getDirichletBCIndices(&dirchlet_ids); CHKERRQ(ierr);

  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpFluxDivergenceAtGaussPts(ufe,"FLUXES"));
  // ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpVDotDivSigma_L2Hdiv(ufe,"VALUES","FLUXES",Aij,F));
  // ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpTauDotSigma_HdivHdiv(ufe,"FLUXES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpValuesAtGaussPts(ufe,"VALUES"));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpDivTauU_HdivL2(ufe,"FLUXES","VALUES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpL2Source(ufe,"VALUES",F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpTauDotSigma_HdivHdiv(ufe,"FLUXES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpVDotDivSigma_L2Hdiv(ufe,"VALUES","FLUXES",Aij,F));
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",ufe.feVol); CHKERRQ(ierr);

  ufe.feTriFluxValue.getOpPtrVector().clear();
  ufe.feTriFluxValue.getOpPtrVector().push_back(new MyUltraWeakFE::OpRhsBcOnValues(ufe,"FLUXES",F));
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_FLUXNEUMANN",ufe.feTriFluxValue); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  // for ksp solver vector is moved into rhs side
  // for snes it is left ond the left
  ierr = VecScale(F,-1); CHKERRQ(ierr);

  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;
  ierr = MatZeroRowsColumnsIS(Aij,dirchlet_ids,1,D0,F); CHKERRQ(ierr);

  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
  //std::cin >> wait;

  //solve
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.set_global_ghost_vector("ULTRAWEAK",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"forces_and_sources_ultra_weak_transport.txt",&viewer); CHKERRQ(ierr);

  const double chop = 1e-4;
  ierr = VecChop(D,chop); CHKERRQ(ierr);
  //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
  VecView(D,viewer);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  //evaluate error
  ufe.feVol.getOpPtrVector().clear();
  ufe.feVol.getOpPtrVector().clear();

  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpValuesGradientAtGaussPts(ufe,"VALUES"));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpFluxDivergenceAtGaussPts(ufe,"FLUXES"));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpError_L2Norm(ufe,"ERROR"));
  ierr = m_field.loop_finite_elements("ULTRAWEAK_CALCULATE_ERROR","ULTRAWEAK_ERROR",ufe.feVol); CHKERRQ(ierr);

  Vec E;
  ierr = m_field.VecCreateGhost("ULTRAWEAK_CALCULATE_ERROR",ROW,&E); CHKERRQ(ierr);
  ierr = m_field.set_local_ghost_vector("ULTRAWEAK_CALCULATE_ERROR",ROW,E,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.set_global_ghost_vector("ULTRAWEAK_CALCULATE_ERROR",ROW,E,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecDestroy(&E); CHKERRQ(ierr);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //calculate residuals
  ufe.feVol.getOpPtrVector().clear();
  ufe.feVol.getOpPtrVector().clear();

  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpFluxDivergenceAtGaussPts(ufe,"FLUXES"));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpVDotDivSigma_L2Hdiv(ufe,"VALUES","FLUXES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpTauDotSigma_HdivHdiv(ufe,"FLUXES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpValuesAtGaussPts(ufe,"VALUES"));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpDivTauU_HdivL2(ufe,"FLUXES","VALUES",Aij,F));
  ufe.feVol.getOpPtrVector().push_back(new MyUltraWeakFE::OpL2Source(ufe,"VALUES",F));

  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",ufe.feVol); CHKERRQ(ierr);

  ufe.feTriFluxValue.getOpPtrVector().clear();
  ufe.feTriFluxValue.getOpPtrVector().push_back(new MyUltraWeakFE::OpRhsBcOnValues(ufe,"FLUXES",F));
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_FLUXNEUMANN",ufe.feTriFluxValue); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  ufe.feTriFluxValue.getOpPtrVector().clear();
  ufe.feTriFluxValue.getOpPtrVector().push_back(new MyUltraWeakFE::OpEvaluateBcOnFluxes(ufe,"FLUXES",F));
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK_FLUXDIRICHLET",ufe.feTriFluxValue); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  VecAXPY(F,-1.,D0);

  double nrm2_F;
  ierr = VecNorm(F,NORM_2,&nrm2_F); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"nrm2_F = %6.4e\n",nrm2_F);

  const double eps = 1e-8;
  if(nrm2_F > eps) {
    //SETERRQ(PETSC_COMM_SELF,MOFEM_ATOM_TEST_INVALID,"problem with residual");
  }

  ierr = ISDestroy(&dirchlet_ids); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&D0); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PostProcVolumeOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);

  ierr = post_proc.addFieldValuesPostProc("VALUES"); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",post_proc);  CHKERRQ(ierr);
  ierr = post_proc.writeFile("out_values.h5m"); CHKERRQ(ierr);
  //rval = post_proc.postProcMesh.write_file("out.vtk","VTK",""); CHKERRQ_MOAB(rval);
  ierr = post_proc.clearOperators(); CHKERRQ(ierr);

  ierr = post_proc.addFieldValuesPostProc("FLUXES"); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ULTRAWEAK","ULTRAWEAK",post_proc);  CHKERRQ(ierr);
  ierr = post_proc.writeFile("out_fluxes.h5m"); CHKERRQ(ierr);
  ierr = post_proc.clearOperators(); CHKERRQ(ierr);

  PostProcVolumeOnRefinedMesh post_proc_error(m_field,false,0);
  ierr = post_proc_error.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc_error.addFieldValuesPostProc("ERROR"); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ULTRAWEAK_CALCULATE_ERROR","ULTRAWEAK_ERROR",post_proc_error);  CHKERRQ(ierr);
  ierr = post_proc_error.writeFile("out_error.h5m"); CHKERRQ(ierr);
  ierr = post_proc_error.clearOperators(); CHKERRQ(ierr);

  //if(pcomm->rank()==0) {
    //EntityHandle fe_meshset = m_field.get_finite_element_meshset("ULTRAWEAK");
    //rval = moab.write_file("error.vtk","VTK","",&fe_meshset,1); CHKERRQ_MOAB(rval); CHKERRQ_MOAB(rval);
  //}

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
