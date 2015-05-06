/** \file arc_length_nonlinear_elasticity.cpp
 * \ingroup nonlinear_elastic_elem
 * \brief nonlinear elasticity (arc-length control)
 *
 * Solves nonlinear elastic problem. Using arc length control.
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

 static char help[] = "\
 -my_file mesh file name\n\
 -my_sr reduction of step size\n\
 -my_ms maximal number of steps\n\n";

#include <MoFEM.hpp>
using namespace MoFEM;

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <DirichletBC.hpp>
#include <ArcLengthTools.hpp>
#include <adolc/adolc.h>
#include <NonLienarElasticElement.hpp>
#include <NeoHookean.hpp>

#include <PotsProcOnRefMesh.hpp>
#include <PostProcStresses.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include <ElasticMaterials.hpp>

#include <SurfacePressureComplexForLazy.hpp>
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 3;
  }

  // use this if your mesh is partotioned and you run code on parts,
  // you can solve very big problems
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);

  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,0); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,1); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,2); CHKERR_PETSC(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  }

  //data stored on mesh for restart
  Tag th_step_size,th_step;
  double def_step_size = 1;
  rval = moab.tag_get_handle("_STEPSIZE",1,MB_TYPE_DOUBLE,th_step_size,MB_TAG_CREAT|MB_TAG_MESH,&def_step_size);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  int def_step = 1;
  rval = moab.tag_get_handle("_STEP",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_MESH,&def_step);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR(rval);
  const void* tag_data_step_size[1];
  EntityHandle root = moab.get_root_set();
  rval = moab.tag_get_by_ptr(th_step_size,&root,1,tag_data_step_size); CHKERR_PETSC(rval);
  double& step_size = *(double *)tag_data_step_size[0];
  const void* tag_data_step[1];
  rval = moab.tag_get_by_ptr(th_step,&root,1,tag_data_step); CHKERR_PETSC(rval);
  int& step = *(int *)tag_data_step[0];
  //end of data stored for restart
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Start step %D and step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));
  BitRefLevel problem_bit_level;

  if(step == 1) {

    problem_bit_level = bit_levels.back();

    Range CubitSideSets_meshsets;
    ierr = m_field.get_cubit_meshsets(SIDESET,CubitSideSets_meshsets); CHKERRQ(ierr);

    //Fields
    ierr = m_field.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);
    ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

    ierr = m_field.add_field("LAMBDA",NOFIELD,1); CHKERRQ(ierr);

    //Field for ArcLength
    ierr = m_field.add_field("X0_SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

    //FE
    ierr = m_field.add_finite_element("ELASTIC"); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr); //this is for parmetis
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

    //Define rows/cols and element data
    ierr = m_field.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
    //elem data
    ierr = m_field.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

    //define problems
    ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //set finite elements for problems
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGTH"); CHKERRQ(ierr);

    //set refinment level for problem
    ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

    //this entity will carray data for this finite element
    EntityHandle meshset_FE_ARC_LENGTH;
    rval = moab.create_meshset(MESHSET_SET,meshset_FE_ARC_LENGTH); CHKERR_PETSC(rval);
    //get LAMBDA field meshset
    EntityHandle meshset_field_LAMBDA = m_field.get_field_meshset("LAMBDA");
    //add LAMBDA field meshset to finite element ARC_LENGTH
    rval = moab.add_entities(meshset_FE_ARC_LENGTH,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
    //add finite element ARC_LENGTH meshset to refinment database (all ref bit leveles)
    ierr = m_field.seed_ref_level_MESHSET(meshset_FE_ARC_LENGTH,BitRefLevel().set()); CHKERRQ(ierr);
    //finally add created meshset to the ARC_LENGTH finite element
    ierr = m_field.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGTH,"ARC_LENGTH",false); CHKERRQ(ierr);

    //set app. order
    ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
    //
    ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

    //add neumman finite elemnets to add static boundary conditions
    ierr = m_field.add_finite_element("NEUAMNN_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      Range tris;
      rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
    }
    //add nodal force element
    ierr = MetaNodalForces::addNodalForceElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);
  }

  PetscBool linear;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-is_linear",&linear,&linear); CHKERRQ(ierr);

  //NeoHookean<adouble> neo_hooke_adouble;
  //NeoHookean<double> neo_hooke_double;
  //NonlinearElasticElement elastic(m_field,2);
  //ierr = elastic.setBlocks(&neo_hooke_double,&neo_hooke_adouble); CHKERRQ(ierr);
  NonlinearElasticElement elastic(m_field,2);
  ElasticMaterials elastic_materials(m_field);
  ierr = elastic_materials.setBlocks(elastic.setOfBlocks); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

  //post_processing
  PostPocOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
  map<int,NonlinearElasticElement::BlockData>::iterator sit = elastic.setOfBlocks.begin();
  for(;sit!=elastic.setOfBlocks.end();sit++) {
    post_proc.getRowOpPtrVector().push_back(
	  new PostPorcStress(
	    post_proc.postProcMesh,
	    post_proc.mapGaussPts,
	    "SPATIAL_POSITION",
	    sit->second,
	    post_proc.commonData));
  }

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  if(step==1) {
    //10 node tets
    Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material,0); CHKERRQ(ierr);
    ierr = m_field.set_field(0,MBVERTEX,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.set_field(0,MBEDGE,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.field_axpy(1.,"MESH_NODE_POSITIONS","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.set_field(0,MBTRI,"SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.set_field(0,MBTET,"SPATIAL_POSITION"); CHKERRQ(ierr);
  }

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = m_field.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //build database
  if(is_partitioned) {
    SETERRQ(PETSC_COMM_SELF,1,"Not implemented, problem with arc-length force multiplayer");
    //ierr = m_field.build_partitioned_problems(PETSC_COMM_WORLD,1); CHKERRQ(ierr);
    //ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS",true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    ierr = m_field.build_problems(); CHKERRQ(ierr);
    ierr = m_field.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  }
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_pressure_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_force_set(); CHKERRQ(ierr);

  //print block sets with materials
  ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  ArcLengthCtx* arc_ctx = new ArcLengthCtx(m_field,"ELASTIC_MECHANICS");

  PetscInt M,N;
  ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(Aij,&m,&n);
  ArcLengthMatShell* mat_ctx = new ArcLengthMatShell(Aij,arc_ctx,"ELASTIC_MECHANICS");
  Mat ShellAij;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)mat_ctx,&ShellAij); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))ArcLengthMatMultShellOp); CHKERRQ(ierr);

  ArcLengthSnesCtx snes_ctx(m_field,"ELASTIC_MECHANICS",arc_ctx);

  Range node_set;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"LoadPath",cit)) {
    EntityHandle meshset = cit->get_meshset();
    Range nodes;
    rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); CHKERR_THROW(rval);
    node_set.merge(nodes);
  }
  PetscPrintf(PETSC_COMM_WORLD,"Nb. nodes in load path: %u\n",node_set.size());

  SphericalArcLengthControl* arc_method_ptr = new SphericalArcLengthControl(arc_ctx);
  SphericalArcLengthControl& arc_method = *arc_method_ptr;

  double scaled_reference_load = 1;
  double *scale_lhs = &(arc_ctx->getFieldData());
  double *scale_rhs = &(scaled_reference_load);
  NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,Aij,arc_ctx->F_lambda,scale_lhs,scale_rhs);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_neumann = neumann_forces.getLoopSpatialFe();
  if(linear) {
    fe_neumann.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
  }
  fe_neumann.uSeF = true;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = fe_neumann.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = fe_neumann.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);
  ierr = m_field.get_problem("ELASTIC_MECHANICS",&my_dirichlet_bc.problemPtr); CHKERRQ(ierr);
  ierr = my_dirichlet_bc.iNitalize(); CHKERRQ(ierr);

  struct MyPrePostProcessFEMethod: public FEMethod {

    FieldInterface& m_field;
    ArcLengthCtx *arc_ptr;

    SpatialPositionsBCFEMethodPreAndPostProc *bC;
    Range &nodeSet;

    MyPrePostProcessFEMethod(FieldInterface& _m_field,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc,Range &node_set):
      m_field(_m_field),arc_ptr(_arc_ptr),bC(bc),nodeSet(node_set) {}

    PetscErrorCode ierr;

      PetscErrorCode preProcess() {
        PetscFunctionBegin;

	//PetscAttachDebugger();
        switch(snes_ctx) {
          case CTX_SNESSETFUNCTION: {
            ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          }
          break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }

        PetscFunctionReturn(0);
      }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
	  //snes_f
          ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode potsProcessLoadPath() {
      PetscFunctionBegin;
      NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problemPtr->numered_dofs_rows);
      Range::iterator nit = nodeSet.begin();
      for(;nit!=nodeSet.end();nit++) {
	NumeredDofMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator dit,hi_dit;
	dit = numered_dofs_rows.get<Ent_mi_tag>().lower_bound(*nit);
	hi_dit = numered_dofs_rows.get<Ent_mi_tag>().upper_bound(*nit);
	for(;dit!=hi_dit;dit++) {
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ","LAMBDA",0,arc_ptr->getFieldData());
	  PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e\n",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
	}
      }
      PetscFunctionReturn(0);
    }

  };

  struct AssembleLambdaFEMethod: public FEMethod {

    FieldInterface& m_field;
    ArcLengthCtx *arc_ptr;

    SpatialPositionsBCFEMethodPreAndPostProc *bC;

    AssembleLambdaFEMethod(FieldInterface& _m_field,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc):
      m_field(_m_field),arc_ptr(_arc_ptr),bC(bc) {}

    PetscErrorCode ierr;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
	  //F_lambda
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  for(vector<int>::iterator vit = bC->dofsIndices.begin();vit!=bC->dofsIndices.end();vit++) {
	    ierr = VecSetValue(arc_ptr->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
	  }
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
	  //add F_lambda
	  ierr = VecAXPY(snes_f,arc_ptr->getFieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->getFieldData());
	  double fnorm;
	  ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);
	}
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      PetscFunctionReturn(0);
    }

  };

  MyPrePostProcessFEMethod pre_post_method(m_field,arc_ctx,&my_dirichlet_bc,node_set);
  AssembleLambdaFEMethod assemble_F_lambda(m_field,arc_ctx,&my_dirichlet_bc);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  PetscReal my_tol;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal atol,rtol,stol;
    PetscInt maxit,maxf;
    ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);
    atol = my_tol;
    rtol = atol*1e2;
    ierr = SNESSetTolerances(snes,atol,rtol,stol,maxit,maxf); CHKERRQ(ierr);
  }

  //
  /*ierr = SNESSetType(snes,SNESSHELL); CHKERRQ(ierr);
  ierr = SNESShellSetContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESShellSetSolve(snes,snes_apply_arc_length); CHKERRQ(ierr);*/
  //

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  PCArcLengthCtx* pc_ctx = new PCArcLengthCtx(Aij,ShellAij,arc_ctx);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,pc_ctx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,PCApplyArcLength); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,PCSetupArcLength); CHKERRQ(ierr);

  if(flg == PETSC_TRUE) {
    PetscReal rtol,atol,dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits); CHKERRQ(ierr);
    atol = my_tol*1e-2;
    rtol = atol*1e-2;
    ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr);
  }


  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_method);
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("ELASTIC",&elastic.getLoopFeRhs()));
  //surface focres and preassures
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",arc_ctx->F_lambda,it->get_msId());  CHKERRQ(ierr);
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //arc length
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NONE",&assemble_F_lambda));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_method);
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("ELASTIC",&elastic.getLoopFeLhs()));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_neumann));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_method));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

  ierr = m_field.set_local_ghost_vector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  PetscScalar step_size_reduction;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_sr",&step_size_reduction,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    step_size_reduction = 1.;
  }

  PetscInt max_steps;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ms",&max_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    max_steps = 5;
  }

  int its_d;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_its_d",&its_d,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    its_d = 4;
  }
  PetscScalar max_reudction = 10,min_reduction = 0.1;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_max_step_reduction",&max_reudction,&flg); CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_min_step_reduction",&min_reduction,&flg); CHKERRQ(ierr);

  double gamma = 0.5,reduction = 1;
  //step = 1;
  if(step == 1) {
    step_size = step_size_reduction;
  } else {
    reduction = step_size_reduction;
    step++;
  }
  double step_size0 = step_size;

  if(step>1) {
    ierr = m_field.set_other_global_ghost_vector(
      "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double x0_nrm;
    ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dlambda);
    ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);
  } else {
    ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
    ierr = arc_ctx->setAlphaBeta(0,1); CHKERRQ(ierr);
  }
  ierr = SnesRhs(snes,D,F,&snes_ctx); CHKERRQ(ierr);

  Vec D0,x00;
  ierr = VecDuplicate(D,&D0); CHKERRQ(ierr);
  ierr = VecDuplicate(arc_ctx->x0,&x00); CHKERRQ(ierr);
  bool converged_state = false;

  for(int jj = 0;step<max_steps;step++,jj++) {

    ierr = VecCopy(D,D0); CHKERRQ(ierr);
    ierr = VecCopy(arc_ctx->x0,x00); CHKERRQ(ierr);

    if(step == 1) {

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
      ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
      ierr = arc_ctx->setAlphaBeta(0,1); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      double dlambda;
      ierr = arc_method.calculateInitDlambda(&dlambda); CHKERRQ(ierr);
      ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRQ(ierr);

    } else if(step == 2) {

      ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);
      ierr = arc_method.calculateDxAndDlambda(D); CHKERRQ(ierr);
      step_size = sqrt(arc_method.calculateLambdaInt());
      step_size0 = step_size;
      ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
      double dlambda = arc_ctx->dlambda;
      double dx_nrm;
      ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
        step,step_size,dlambda,dx_nrm,arc_ctx->dx2
      ); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
      ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRQ(ierr);

    } else {

      if(jj == 0) {
        step_size0 = step_size;
      }

      ierr = arc_method.calculateDxAndDlambda(D); CHKERRQ(ierr);
      step_size *= reduction;
      if(step_size > max_reudction*step_size0) {
        step_size = max_reudction*step_size0;
      } else if(step_size<min_reduction*step_size0) {
        step_size = min_reduction*step_size0;
      }
      ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
      double dlambda = reduction*arc_ctx->dlambda;
      double dx_nrm;
      ierr = VecScale(arc_ctx->dx,reduction); CHKERRQ(ierr);
      ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
        step,step_size,dlambda,dx_nrm,arc_ctx->dx2
      ); CHKERRQ(ierr);
      ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
      ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRQ(ierr);

    }

    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
    if(reason < 0) {

      ierr = VecCopy(D0,D); CHKERRQ(ierr);
      ierr = VecCopy(x00,arc_ctx->x0); CHKERRQ(ierr);

      double x0_nrm;
      ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dlambda);
      ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);


      reduction = 0.1;
      converged_state = false;

      continue;

    } else {

      if(step > 1 && converged_state) {

        reduction = pow((double)its_d/(double)(its+1),gamma);
        if(step_size >= max_reudction*step_size0 && reduction > 1) {
          reduction = 1;
        } else if(step_size <= min_reduction*step_size0 && reduction < 1) {
          reduction = 1;
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n",reduction); CHKERRQ(ierr);
      }

      //Save data on mesh
      ierr = m_field.set_global_ghost_vector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = m_field.set_other_global_ghost_vector(
        "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_REVERSE
      ); CHKERRQ(ierr);
      converged_state = true;

    }

    if(step % 1 == 0) {
      //Save restart file
      ostringstream sss;
      sss << "restart_" << step << ".h5m";
      rval = moab.write_file(sss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
      //Save data on mesh
      ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
      ostringstream o1;
      o1 << "out_" << step << ".h5m";
      rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
    }

    ierr = pre_post_method.potsProcessLoadPath(); CHKERRQ(ierr);

  }

  ierr = VecDestroy(&D0); CHKERRQ(ierr);
  ierr = VecDestroy(&x00); CHKERRQ(ierr);

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  delete mat_ctx;
  delete pc_ctx;
  delete arc_ctx;
  delete arc_method_ptr;

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
