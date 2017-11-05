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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include <ElasticMaterials.hpp>
#include <NeoHookean.hpp>

#include <SurfacePressureComplexForLazy.hpp>

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRG(ierr);
    if(flg != PETSC_TRUE) {
      order = 3;
    }

    // use this if your mesh is partitioned and you run code on parts,
    // you can solve very big problems
    PetscBool is_partitioned = PETSC_FALSE;
    ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRG(ierr);

    if(is_partitioned == PETSC_TRUE) {
      //Read mesh to MOAB
      const char *option;
      option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
      rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
    } else {
      const char *option;
      option = "";
      rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
    }

    //data stored on mesh for restart
    Tag th_step_size,th_step;
    double def_step_size = 1;
    rval = moab.tag_get_handle("_STEPSIZE",1,MB_TYPE_DOUBLE,th_step_size,MB_TAG_CREAT|MB_TAG_MESH,&def_step_size);
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERRG(rval);
    int def_step = 1;
    rval = moab.tag_get_handle("_STEP",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_MESH,&def_step);
    if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
    CHKERRG(rval);
    const void* tag_data_step_size[1];
    EntityHandle root = moab.get_root_set();
    rval = moab.tag_get_by_ptr(th_step_size,&root,1,tag_data_step_size); CHKERRG(rval);
    double& step_size = *(double *)tag_data_step_size[0];
    const void* tag_data_step[1];
    rval = moab.tag_get_by_ptr(th_step,&root,1,tag_data_step); CHKERRG(rval);
    int& step = *(int *)tag_data_step[0];
    //end of data stored for restart
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Start step %D and step_size = %6.4e\n",step,step_size); CHKERRG(ierr);

    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    //ref meshset ref level 0
    ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRG(ierr);
    std::vector<BitRefLevel> bit_levels;
    bit_levels.push_back(BitRefLevel().set(0));
    BitRefLevel problem_bit_level;

    if(step == 1) {

      problem_bit_level = bit_levels.back();

      //Fields
      ierr = m_field.add_field("SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);
      ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);

      ierr = m_field.add_field("LAMBDA",NOFIELD,NOBASE,1); CHKERRG(ierr);

      //Field for ArcLength
      ierr = m_field.add_field("X0_SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);

      //FE
      ierr = m_field.add_finite_element("ELASTIC"); CHKERRG(ierr);
      ierr = m_field.add_finite_element("ARC_LENGTH"); CHKERRG(ierr);

      //Define rows/cols and element data
      ierr = m_field.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRG(ierr); //this is for parmetis
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRG(ierr);

      //Define rows/cols and element data
      ierr = m_field.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRG(ierr);
      //elem data
      ierr = m_field.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRG(ierr);

      //define problems
      ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRG(ierr);

      //set finite elements for problems
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRG(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGTH"); CHKERRG(ierr);

      //set refinement level for problem
      ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRG(ierr);

      //add entitities (by tets) to the field
      ierr = m_field.add_ents_to_field_by_type(0,MBTET,"SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRG(ierr);

      // Setting up LAMBDA field and ARC_LENGTH interface
      {
        //Add dummy no-field vertex
        EntityHandle no_field_vertex;
        {
          const double coords[] = {0,0,0};
          rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERRG(rval);
          Range range_no_field_vertex;
          range_no_field_vertex.insert(no_field_vertex);
          ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRG(ierr);
          EntityHandle lambda_meshset = m_field.get_field_meshset("LAMBDA");
          rval = m_field.get_moab().add_entities(lambda_meshset,range_no_field_vertex); CHKERRG(rval);
        }
        //this entity will carray data for this finite element
        EntityHandle meshset_fe_arc_length;
        {
          rval = moab.create_meshset(MESHSET_SET,meshset_fe_arc_length); CHKERRG(rval);
          rval = moab.add_entities(meshset_fe_arc_length,&no_field_vertex,1); CHKERRG(rval);
          ierr = m_field.seed_ref_level_MESHSET(meshset_fe_arc_length,BitRefLevel().set()); CHKERRG(ierr);
        }
        //finally add created meshset to the ARC_LENGTH finite element
        ierr = m_field.add_ents_to_finite_element_by_MESHSET(meshset_fe_arc_length,"ARC_LENGTH",false); CHKERRG(ierr);
      }

      //set app. order
      ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRG(ierr);
      //
      ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRG(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRG(ierr);

      //add Neumman finite elements to add static boundary conditions
      ierr = m_field.add_finite_element("NEUAMNN_FE"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","MESH_NODE_POSITIONS"); CHKERRG(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRG(ierr);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
        Range tris;
        rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRG(rval);
        ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI,"NEUAMNN_FE"); CHKERRG(ierr);
      }
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
        Range tris;
        rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRG(rval);
        ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI,"NEUAMNN_FE"); CHKERRG(ierr);
      }
      //add nodal force element
      ierr = MetaNodalForces::addElement(m_field,"SPATIAL_POSITION"); CHKERRG(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRG(ierr);
    }

    PetscBool linear;
    ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-is_linear",&linear,&linear); CHKERRG(ierr);

    //NeoHookean<adouble> neo_hooke_adouble;
    //NeoHookean<double> neo_hooke_double;
    //NonlinearElasticElement elastic(m_field,2);
    //ierr = elastic.setBlocks(&neo_hooke_double,&neo_hooke_adouble); CHKERRG(ierr);
    NonlinearElasticElement elastic(m_field,2);
    ElasticMaterials elastic_materials(m_field);
    ierr = elastic_materials.setBlocks(elastic.setOfBlocks); CHKERRG(ierr);
    ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRG(ierr);
    ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRG(ierr);

    //post_processing
    PostProcVolumeOnRefinedMesh post_proc(m_field);
    ierr = post_proc.generateReferenceElementMesh(); CHKERRG(ierr);
    ierr = post_proc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRG(ierr);
    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRG(ierr);
    ierr = post_proc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRG(ierr);
    std::map<int,NonlinearElasticElement::BlockData>::iterator sit = elastic.setOfBlocks.begin();
    for(;sit!=elastic.setOfBlocks.end();sit++) {
      post_proc.getOpPtrVector().push_back(
        new PostProcStress(
          post_proc.postProcMesh,
          post_proc.mapGaussPts,
          "SPATIAL_POSITION",
          sit->second,
          post_proc.commonData)
        );
      }

      //build field
      ierr = m_field.build_fields(); CHKERRG(ierr);
      if(step==1) {
        //10 node tets
        Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
        ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material,0); CHKERRG(ierr);
        ierr = m_field.set_field(0,MBVERTEX,"SPATIAL_POSITION"); CHKERRG(ierr);
        ierr = m_field.set_field(0,MBEDGE,"SPATIAL_POSITION"); CHKERRG(ierr);
        ierr = m_field.field_axpy(1.,"MESH_NODE_POSITIONS","SPATIAL_POSITION"); CHKERRG(ierr);
        ierr = m_field.set_field(0,MBTRI,"SPATIAL_POSITION"); CHKERRG(ierr);
        ierr = m_field.set_field(0,MBTET,"SPATIAL_POSITION"); CHKERRG(ierr);
      }

      //build finite elements
      ierr = m_field.build_finite_elements(); CHKERRG(ierr);

      //build adjacencies
      ierr = m_field.build_adjacencies(problem_bit_level); CHKERRG(ierr);


      ProblemsManager *prb_mng_ptr;
      ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);
      //build database
      if(is_partitioned) {
        SETERRQ(PETSC_COMM_SELF,1,"Not implemented, problem with arc-length force multiplayer");
      } else {
        ierr = prb_mng_ptr->buildProblem("ELASTIC_MECHANICS",true); CHKERRG(ierr);
        ierr = prb_mng_ptr->partitionProblem("ELASTIC_MECHANICS"); CHKERRG(ierr);
        ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS"); CHKERRG(ierr);
      }
      ierr = prb_mng_ptr->partitionGhostDofs("ELASTIC_MECHANICS"); CHKERRG(ierr);

      //print bcs
      MeshsetsManager *mmanager_ptr;
      ierr = m_field.getInterface(mmanager_ptr); CHKERRG(ierr);
      ierr = mmanager_ptr->printDisplacementSet(); CHKERRG(ierr);
      ierr = mmanager_ptr->printForceSet(); CHKERRG(ierr);
      //print block sets with materials
      ierr = mmanager_ptr->printMaterialsSet(); CHKERRG(ierr);

      //create matrices
      Vec F;
      ierr = m_field.getInterface<VecManager>()->vecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRG(ierr);
      Vec D;
      ierr = VecDuplicate(F,&D); CHKERRG(ierr);
      Mat Aij;
      ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRG(ierr);

      ArcLengthCtx* arc_ctx = new ArcLengthCtx(m_field,"ELASTIC_MECHANICS");

      PetscInt M,N;
      ierr = MatGetSize(Aij,&M,&N); CHKERRG(ierr);
      PetscInt m,n;
      ierr = MatGetLocalSize(Aij,&m,&n); CHKERRG(ierr);
      ArcLengthMatShell* mat_ctx = new ArcLengthMatShell(Aij,arc_ctx,"ELASTIC_MECHANICS");
      Mat ShellAij;
      ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)mat_ctx,&ShellAij); CHKERRG(ierr);
      ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))ArcLengthMatMultShellOp); CHKERRG(ierr);

      ArcLengthSnesCtx snes_ctx(m_field,"ELASTIC_MECHANICS",arc_ctx);
      ///< do not very if element of given name exist when do loop over elements
      snes_ctx.bH = MF_ZERO;

      Range node_set;
      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"LoadPath",cit)) {
        EntityHandle meshset = cit->getMeshset();
        Range nodes;
        rval = moab.get_entities_by_type(meshset,MBVERTEX,nodes,true); MOAB_THROW(rval);
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
        ierr = fe_neumann.addForce(it->getMeshsetId()); CHKERRG(ierr);
      }
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
        ierr = fe_neumann.addPreassure(it->getMeshsetId()); CHKERRG(ierr);
      }
      DirichletSpatialPositionsBc my_dirichlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);
      ierr = m_field.get_problem("ELASTIC_MECHANICS",&my_dirichlet_bc.problemPtr); CHKERRG(ierr);
      ierr = my_dirichlet_bc.iNitalize(); CHKERRG(ierr);

      struct AssembleRhsVectors: public FEMethod {

        ArcLengthCtx *arcPtr;
        Range &nodeSet;

        AssembleRhsVectors(
          ArcLengthCtx *arc_ptr,
          Range &node_set
        ):
        arcPtr(arc_ptr),
        nodeSet(node_set) {}

        MoFEMErrorCode preProcess() {
          MoFEMFunctionBeginHot;

          //PetscAttachDebugger();
          switch(snes_ctx) {
            case CTX_SNESSETFUNCTION: {
              ierr = VecZeroEntries(snes_f); CHKERRG(ierr);
              ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
              ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
              ierr = VecZeroEntries(arcPtr->F_lambda); CHKERRG(ierr);
              ierr = VecGhostUpdateBegin(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
              ierr = VecGhostUpdateEnd(arcPtr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
            }
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }

          MoFEMFunctionReturnHot(0);
        }

        MoFEMErrorCode postProcess() {
          MoFEMFunctionBeginHot;
          switch(snes_ctx) {
            case CTX_SNESSETFUNCTION: {
              //snes_f
              ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
              ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
              ierr = VecAssemblyBegin(snes_f); CHKERRG(ierr);
              ierr = VecAssemblyEnd(snes_f); CHKERRG(ierr);
            }
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }
          MoFEMFunctionReturnHot(0);
        }

        MoFEMErrorCode potsProcessLoadPath() {
          MoFEMFunctionBeginHot;
          boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_rows = problemPtr->getNumeredDofsRows();
          Range::iterator nit = nodeSet.begin();
          for(;nit!=nodeSet.end();nit++) {
            NumeredDofEntityByEnt::iterator dit,hi_dit;
            dit = numered_dofs_rows->get<Ent_mi_tag>().lower_bound(*nit);
            hi_dit = numered_dofs_rows->get<Ent_mi_tag>().upper_bound(*nit);
            for(;dit!=hi_dit;dit++) {
              PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ","LAMBDA",0,arcPtr->getFieldData());
              PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e\n",dit->get()->getName().c_str(),dit->get()->getDofCoeffIdx(),dit->get()->getFieldData());
            }
          }
          MoFEMFunctionReturnHot(0);
        }

      };

      struct AddLambdaVectorToFinternal: public FEMethod {

        ArcLengthCtx *arcPtr;
        DirichletSpatialPositionsBc *bC;

        AddLambdaVectorToFinternal(
          ArcLengthCtx *arc_ptr,
          DirichletSpatialPositionsBc *bc
        ):
        arcPtr(arc_ptr),
        bC(bc) {}



        MoFEMErrorCode preProcess() {
          MoFEMFunctionBeginHot;
          MoFEMFunctionReturnHot(0);
        }
        MoFEMErrorCode operator()() {
          MoFEMFunctionBeginHot;
          MoFEMFunctionReturnHot(0);
        }
        MoFEMErrorCode postProcess() {
          MoFEMFunctionBeginHot;
          switch(snes_ctx) {
            case CTX_SNESSETFUNCTION: {
              //F_lambda
              ierr = VecGhostUpdateBegin(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
              ierr = VecGhostUpdateEnd(arcPtr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
              ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRG(ierr);
              ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRG(ierr);
              for(std::vector<int>::iterator vit = bC->dofsIndices.begin();vit!=bC->dofsIndices.end();vit++) {
                ierr = VecSetValue(arcPtr->F_lambda,*vit,0,INSERT_VALUES); CHKERRG(ierr);
              }
              ierr = VecAssemblyBegin(arcPtr->F_lambda); CHKERRG(ierr);
              ierr = VecAssemblyEnd(arcPtr->F_lambda); CHKERRG(ierr);
              ierr = VecDot(arcPtr->F_lambda,arcPtr->F_lambda,&arcPtr->F_lambda2); CHKERRG(ierr);
              PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arcPtr->F_lambda2);
              //add F_lambda
              ierr = VecAXPY(snes_f,arcPtr->getFieldData(),arcPtr->F_lambda); CHKERRG(ierr);
              PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arcPtr->getFieldData());
              double fnorm;
              ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRG(ierr);
              PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);
            }
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
          }
          MoFEMFunctionReturnHot(0);
        }

      };

      AssembleRhsVectors pre_post_method(arc_ctx,node_set);
      AddLambdaVectorToFinternal assemble_F_lambda(arc_ctx,&my_dirichlet_bc);

      SNES snes;
      ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRG(ierr);
      ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRG(ierr);
      ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRG(ierr);
      ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&snes_ctx); CHKERRG(ierr);
      ierr = SNESSetFromOptions(snes); CHKERRG(ierr);

      PetscReal my_tol;
      ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRG(ierr);
      if(flg == PETSC_TRUE) {
        PetscReal atol,rtol,stol;
        PetscInt maxit,maxf;
        ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRG(ierr);
        atol = my_tol;
        rtol = atol*1e2;
        ierr = SNESSetTolerances(snes,atol,rtol,stol,maxit,maxf); CHKERRG(ierr);
      }

      //
      /*ierr = SNESSetType(snes,SNESSHELL); CHKERRG(ierr);
      ierr = SNESShellSetContext(snes,&snes_ctx); CHKERRG(ierr);
      ierr = SNESShellSetSolve(snes,snes_apply_arc_length); CHKERRG(ierr);*/
      //

      KSP ksp;
      ierr = SNESGetKSP(snes,&ksp); CHKERRG(ierr);
      PC pc;
      ierr = KSPGetPC(ksp,&pc); CHKERRG(ierr);
      PCArcLengthCtx* pc_ctx = new PCArcLengthCtx(ShellAij,Aij,arc_ctx);
      ierr = PCSetType(pc,PCSHELL); CHKERRG(ierr);
      ierr = PCShellSetContext(pc,pc_ctx); CHKERRG(ierr);
      ierr = PCShellSetApply(pc,PCApplyArcLength); CHKERRG(ierr);
      ierr = PCShellSetSetUp(pc,PCSetupArcLength); CHKERRG(ierr);

      if(flg == PETSC_TRUE) {
        PetscReal rtol,atol,dtol;
        PetscInt maxits;
        ierr = KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits); CHKERRG(ierr);
        atol = my_tol*1e-2;
        rtol = atol*1e-2;
        ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRG(ierr);
      }

      SnesCtx::FEMethodsSequence& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
      snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
      snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_method);
      loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeRhs()));
      //surface forces and pressures
      loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("NEUAMNN_FE",&fe_neumann));

      //edge forces
      boost::ptr_map<std::string,EdgeForce> edge_forces;
      string fe_name_str ="FORCE_FE";
      edge_forces.insert(fe_name_str,new EdgeForce(m_field));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
        ierr = edge_forces.at(fe_name_str).addForce("SPATIAL_POSITION",arc_ctx->F_lambda,it->getMeshsetId());  CHKERRG(ierr);
      }
      for(
        boost::ptr_map<std::string,EdgeForce>::iterator eit = edge_forces.begin();
        eit!=edge_forces.end();eit++
      ) {
        loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr(eit->first,&eit->second->getLoopFe()));
      }

      //nodal forces
      boost::ptr_map<std::string,NodalForce> nodal_forces;
      // string fe_name_str ="FORCE_FE";
      nodal_forces.insert(fe_name_str,new NodalForce(m_field));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
        ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",arc_ctx->F_lambda,it->getMeshsetId());  CHKERRG(ierr);
      }
      for(
        boost::ptr_map<std::string,NodalForce>::iterator fit = nodal_forces.begin();
        fit!=nodal_forces.end();fit++) {
          loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr(fit->first,&fit->second->getLoopFe()));
        }

        //arc length
        loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("NONE",&assemble_F_lambda));
        loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ARC_LENGTH",&arc_method));
        snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_method);
        snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

        SnesCtx::FEMethodsSequence& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
        snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
        loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeLhs()));
        loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("NEUAMNN_FE",&fe_neumann));
        loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ARC_LENGTH",&arc_method));
        snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

        ierr = m_field.getInterface<VecManager>()->setLocalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRG(ierr);

        PetscScalar step_size_reduction;
        ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-my_sr",&step_size_reduction,&flg); CHKERRG(ierr);
        if(flg != PETSC_TRUE) {
          step_size_reduction = 1.;
        }

        PetscInt max_steps;
        ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_ms",&max_steps,&flg); CHKERRG(ierr);
        if(flg != PETSC_TRUE) {
          max_steps = 5;
        }

        int its_d;
        ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_its_d",&its_d,&flg); CHKERRG(ierr);
        if(flg != PETSC_TRUE) {
          its_d = 4;
        }
        PetscScalar max_reudction = 10,min_reduction = 0.1;
        ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-my_max_step_reduction",&max_reudction,&flg); CHKERRG(ierr);
        ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-my_min_step_reduction",&min_reduction,&flg); CHKERRG(ierr);

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
          ierr = m_field.getInterface<VecManager>()->setOtherGlobalGhostVector(
            "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",
            COL,arc_ctx->x0,INSERT_VALUES,SCATTER_FORWARD
          ); CHKERRG(ierr);
          double x0_nrm;
          ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRG(ierr);
          ierr = PetscPrintf(
            PETSC_COMM_WORLD,
            "\tRead x0_nrm = %6.4e dlambda = %6.4e\n",
            x0_nrm,arc_ctx->dLambda
          );
          ierr = arc_ctx->setAlphaBeta(1,0); CHKERRG(ierr);
        } else {
          ierr = arc_ctx->setS(step_size); CHKERRG(ierr);
          ierr = arc_ctx->setAlphaBeta(0,1); CHKERRG(ierr);
        }
        ierr = SnesRhs(snes,D,F,&snes_ctx); CHKERRG(ierr);

        Vec D0,x00;
        ierr = VecDuplicate(D,&D0); CHKERRG(ierr);
        ierr = VecDuplicate(arc_ctx->x0,&x00); CHKERRG(ierr);
        bool converged_state = false;

        for(int jj = 0;step<max_steps;step++,jj++) {

          ierr = VecCopy(D,D0); CHKERRG(ierr);
          ierr = VecCopy(arc_ctx->x0,x00); CHKERRG(ierr);

          if(step == 1) {

            ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D step_size = %6.4e\n",step,step_size); CHKERRG(ierr);
            ierr = arc_ctx->setS(step_size); CHKERRG(ierr);
            ierr = arc_ctx->setAlphaBeta(0,1); CHKERRG(ierr);
            ierr = VecCopy(D,arc_ctx->x0); CHKERRG(ierr);
            double dlambda;
            ierr = arc_method.calculateInitDlambda(&dlambda); CHKERRG(ierr);
            ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRG(ierr);

          } else if(step == 2) {

            ierr = arc_ctx->setAlphaBeta(1,0); CHKERRG(ierr);
            ierr = arc_method.calculateDxAndDlambda(D); CHKERRG(ierr);
            step_size = sqrt(arc_method.calculateLambdaInt());
            step_size0 = step_size;
            ierr = arc_ctx->setS(step_size); CHKERRG(ierr);
            double dlambda = arc_ctx->dLambda;
            double dx_nrm;
            ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRG(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,
              "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
              step,step_size,dlambda,dx_nrm,arc_ctx->dx2
            ); CHKERRG(ierr);
            ierr = VecCopy(D,arc_ctx->x0); CHKERRG(ierr);
            ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRG(ierr);
            ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRG(ierr);

          } else {

            if(jj == 0) {
              step_size0 = step_size;
            }

            ierr = arc_method.calculateDxAndDlambda(D); CHKERRG(ierr);
            step_size *= reduction;
            if(step_size > max_reudction*step_size0) {
              step_size = max_reudction*step_size0;
            } else if(step_size<min_reduction*step_size0) {
              step_size = min_reduction*step_size0;
            }
            ierr = arc_ctx->setS(step_size); CHKERRG(ierr);
            double dlambda = reduction*arc_ctx->dLambda;
            double dx_nrm;
            ierr = VecScale(arc_ctx->dx,reduction); CHKERRG(ierr);
            ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRG(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,
              "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
              step,step_size,dlambda,dx_nrm,arc_ctx->dx2
            ); CHKERRG(ierr);
            ierr = VecCopy(D,arc_ctx->x0); CHKERRG(ierr);
            ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRG(ierr);
            ierr = arc_method.setDlambdaToX(D,dlambda); CHKERRG(ierr);

          }

          ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRG(ierr);
          int its;
          ierr = SNESGetIterationNumber(snes,&its); CHKERRG(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRG(ierr);

          SNESConvergedReason reason;
          ierr = SNESGetConvergedReason(snes,&reason); CHKERRG(ierr);
          if(reason < 0) {

            ierr = VecCopy(D0,D); CHKERRG(ierr);
            ierr = VecCopy(x00,arc_ctx->x0); CHKERRG(ierr);

            double x0_nrm;
            ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRG(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dLambda);
            ierr = arc_ctx->setAlphaBeta(1,0); CHKERRG(ierr);


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
              ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n",reduction); CHKERRG(ierr);
            }

            //Save data on mesh
            ierr = m_field.getInterface<VecManager>()->setGlobalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRG(ierr);
            ierr = m_field.getInterface<VecManager>()->setOtherGlobalGhostVector(
              "ELASTIC_MECHANICS","SPATIAL_POSITION","X0_SPATIAL_POSITION",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_REVERSE
            ); CHKERRG(ierr);
            converged_state = true;

          }

          if(step % 1 == 0) {
            //Save restart file
            // #ifdef MOAB_HDF5_PARALLEL
            //   std::ostringstream sss;
            //   sss << "restart_" << step << ".h5m";
            //   rval = moab.write_file(sss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERRG(rval);
            // #else
            // #warning "No parallel HDF5, no writing restart file"
            // #endif
            //Save data on mesh
            ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRG(ierr);
            std::ostringstream o1;
            o1 << "out_" << step << ".h5m";
            ierr = post_proc.writeFile(o1.str().c_str()); CHKERRG(ierr);
          }

          ierr = pre_post_method.potsProcessLoadPath(); CHKERRG(ierr);

        }

        ierr = VecDestroy(&D0); CHKERRG(ierr);
        ierr = VecDestroy(&x00); CHKERRG(ierr);

        //detroy matrices
        ierr = VecDestroy(&F); CHKERRG(ierr);
        ierr = VecDestroy(&D); CHKERRG(ierr);
        ierr = MatDestroy(&Aij); CHKERRG(ierr);
        ierr = MatDestroy(&ShellAij); CHKERRG(ierr);
        ierr = SNESDestroy(&snes); CHKERRG(ierr);

        delete mat_ctx;
        delete pc_ctx;
        delete arc_ctx;
        delete arc_method_ptr;


      } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);

  return 0;

}
