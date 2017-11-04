/** \file arc_length_interface.cpp
  * \brief Example of arc-length with interface element

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

#include <InterfaceGapArcLengthControl.hpp>
#include <CohesiveInterfaceElement.hpp>
#include <Hooke.hpp>

using namespace boost::numeric;




static char help[] = "...\n\n";

#define DATAFILENAME "load_disp.txt"

namespace CohesiveElement {

  struct ArcLengthElement: public ArcLengthIntElemFEMethod {
    MoFEM::Interface& mField;
    Range postProcNodes;
    ArcLengthElement(MoFEM::Interface& m_field,ArcLengthCtx *arc_ptr):
    ArcLengthIntElemFEMethod(m_field.get_moab(),arc_ptr),
    mField(m_field) {

      for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"LoadPath",cit)) {
        EntityHandle meshset = cit->getMeshset();
        Range nodes;
        rval = mOab.get_entities_by_type(meshset,MBVERTEX,nodes,true); MOAB_THROW(rval);
        postProcNodes.merge(nodes);
      }

      PetscPrintf(PETSC_COMM_WORLD,"Nb. PostProcNodes %lu\n",postProcNodes.size());

    };

    MoFEMErrorCode postProcessLoadPath() {
      MoFEMFunctionBeginHot;
      FILE *datafile;
      PetscFOpen(PETSC_COMM_SELF,DATAFILENAME,"a+",&datafile);
      boost::shared_ptr<NumeredDofEntity_multiIndex> numered_dofs_rows = problemPtr->getNumeredDofsRows();
      NumeredDofEntityByFieldName::iterator lit;
      lit = numered_dofs_rows->get<FieldName_mi_tag>().find("LAMBDA");
      if(lit == numered_dofs_rows->get<FieldName_mi_tag>().end()) {
        fclose(datafile);
        MoFEMFunctionReturnHot(0);
      }
      Range::iterator nit = postProcNodes.begin();
      for(;nit!=postProcNodes.end();nit++) {
        NumeredDofEntityByEnt::iterator dit,hi_dit;
        dit = numered_dofs_rows->get<Ent_mi_tag>().lower_bound(*nit);
        hi_dit = numered_dofs_rows->get<Ent_mi_tag>().upper_bound(*nit);
        double coords[3];
        rval = mOab.get_coords(&*nit,1,coords);  MOAB_THROW(rval);
        for(;dit!=hi_dit;dit++) {
          PetscPrintf(
            PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ",
            lit->get()->getName().c_str(),lit->get()->getDofCoeffIdx(),lit->get()->getFieldData()
          );
          PetscPrintf(
            PETSC_COMM_WORLD,"%s [ %d ] %6.4e ",
            dit->get()->getName().c_str(),dit->get()->getDofCoeffIdx(),dit->get()->getFieldData()
          );
          PetscPrintf(
            PETSC_COMM_WORLD,"-> %3.4f %3.4f %3.4f\n",coords[0],coords[1],coords[2]
          );
          PetscFPrintf(
            PETSC_COMM_WORLD,datafile,"%6.4e %6.4e ",dit->get()->getFieldData(),lit->get()->getFieldData()
          );
        }
      }
      PetscFPrintf(PETSC_COMM_WORLD,datafile,"\n");
      fclose(datafile);
      MoFEMFunctionReturnHot(0);
    }

  };

  struct AssembleRhsVectors: public FEMethod {

    MoFEM::Interface& mField;
    Vec &bodyForce;
    ArcLengthCtx *arcPtr;

    AssembleRhsVectors(
      MoFEM::Interface& m_field,Vec &body_force,ArcLengthCtx *arc_ptr
    ):
    mField(m_field),
    bodyForce(body_force),
    arcPtr(arc_ptr) {
    }

    

    MoFEMErrorCode preProcess() {
      MoFEMFunctionBeginHot;

      switch(snes_ctx) {
        case CTX_SNESNONE: {}
        break;
        case CTX_SNESSETFUNCTION: {
          ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
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
        case CTX_SNESNONE: {}
        break;
        case CTX_SNESSETFUNCTION: {
          ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
          //add F_lambda
          ierr = VecAXPY(snes_f,arcPtr->getFieldData(),arcPtr->F_lambda); CHKERRQ(ierr);
          ierr = VecAXPY(snes_f,-1.,bodyForce); CHKERRQ(ierr);
          PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arcPtr->getFieldData());
          //snes_f norm
          double fnorm;
          ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);
          PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);
        }
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      MoFEMFunctionReturnHot(0);
    }

  };

}

using namespace CohesiveElement;

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    moab::Core mb_instance;
    moab::Interface& moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //Reade parameters from line command
    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    PetscScalar step_size_reduction;
    ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-my_sr",&step_size_reduction,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      step_size_reduction = 1.;
    }

    PetscInt max_steps;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_ms",&max_steps,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      max_steps = 5;
    }

    int its_d;
    ierr = PetscOptionsGetInt(PETSC_NULL,"","-my_its_d",&its_d,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      its_d = 6;
    }

    PetscInt order;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      order = 2;
    }

    //Check if new start or restart. If new start, delete previous load_disp.txt
    if (std::string(mesh_file_name).find("restart") == std::string::npos) {
      remove(DATAFILENAME);
    }

    //Read mesh to MOAB
    const char *option;
    option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

    //Data stored on mesh for restart
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
    EntityHandle root = 0;
    rval = moab.tag_get_by_ptr(th_step_size,&root,1,tag_data_step_size); CHKERRG(rval);
    double& step_size = *(double *)tag_data_step_size[0];
    const void* tag_data_step[1];
    rval = moab.tag_get_by_ptr(th_step,&root,1,tag_data_step); CHKERRG(rval);
    int& step = *(int *)tag_data_step[0];
    //end of data stored for restart
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Start step %D and step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);

    //Create MoFEM 2database
    MoFEM::Core core(moab);
    MoFEM::Interface& m_field = core;

    PrismInterface *interface_ptr;
    ierr = m_field.getInterface(interface_ptr); CHKERRQ(ierr);

    Tag th_my_ref_level;
    BitRefLevel def_bit_level = 0;
    rval = m_field.get_moab().tag_get_handle("_MY_REFINEMENT_LEVEL",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_my_ref_level,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level);
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    BitRefLevel *ptr_bit_level0;
    rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERRG(rval);
    BitRefLevel& bit_level0 = *ptr_bit_level0;
    BitRefLevel problem_bit_level = bit_level0;

    if(step == 1) {

      //ref meshset ref level 0
      ierr = m_field.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
      std::vector<BitRefLevel> bit_levels;
      bit_levels.push_back(BitRefLevel().set(0));

      int ll = 1;

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|INTERFACESET,cit)) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->getMeshsetId()); CHKERRQ(ierr);
        EntityHandle cubit_meshset = cit->getMeshset();
        {
          //get tet entities form back bit_level
          EntityHandle ref_level_meshset = 0;
          rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERRG(rval);
          ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
            bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset
          ); CHKERRQ(ierr);
          ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
            bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset
          ); CHKERRQ(ierr);
          Range ref_level_tets;
          rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERRG(rval);
          //get faces and test to split
          ierr = interface_ptr->getSides(cubit_meshset,bit_levels.back(),true,0); CHKERRQ(ierr);
          //set new bit level
          bit_levels.push_back(BitRefLevel().set(ll++));
          //split faces and
          ierr = interface_ptr->splitSides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
          //clean meshsets
          rval = moab.delete_entities(&ref_level_meshset,1); CHKERRG(rval);
        }
        // Update cubit meshsets
        for(_IT_CUBITMESHSETS_FOR_LOOP_(m_field,ciit)) {
          EntityHandle cubit_meshset = ciit->meshset;
          ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
          ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
          ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
          ierr = m_field.getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
        }
      }

      bit_level0 = bit_levels.back();
      problem_bit_level = bit_level0;

      /***/
      //Define problem

      //Fields
      ierr = m_field.add_field("DISPLACEMENT",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);
      ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);

      ierr = m_field.add_field("LAMBDA",NOFIELD,NOBASE,1); CHKERRQ(ierr);
      //Field for ArcLength
      ierr = m_field.add_field("X0_DISPLACEMENT",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRQ(ierr);

      //FE
      ierr = m_field.add_finite_element("ELASTIC"); CHKERRQ(ierr);
      //Define rows/cols and element data
      ierr = m_field.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("ELASTIC","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ELASTIC","LAMBDA"); CHKERRQ(ierr); //this is for paremtis
      ierr = m_field.modify_finite_element_add_field_data("ELASTIC","LAMBDA"); CHKERRQ(ierr);

      //FE Interface
      ierr = m_field.add_finite_element("INTERFACE"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("INTERFACE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      //FE ArcLength
      ierr = m_field.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);
      //Define rows/cols and element data
      ierr = m_field.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
      //elem data
      ierr = m_field.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

      //define problems
      ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

      //set finite elements for problem
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ARC_LENGTH"); CHKERRQ(ierr);

      //set refinement level for problem
      ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

      /***/
      //Declare problem

      //add entitities (by tets) to the field
      ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      //add finite elements entities
      ierr = m_field.add_ents_to_finite_element_by_bit_ref(problem_bit_level,BitRefLevel().set(),"ELASTIC",MBTET); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_finite_element_by_bit_ref(problem_bit_level,BitRefLevel().set(),"INTERFACE",MBPRISM); CHKERRQ(ierr);

      // Setting up LAMBDA field and ARC_LENGTH interface
      {
        //Add dummy no-field vertex
        EntityHandle no_field_vertex;
        {
          const double coords[] = {0,0,0};
          rval = m_field.get_moab().create_vertex(coords,no_field_vertex); CHKERRG(rval);
          Range range_no_field_vertex;
          range_no_field_vertex.insert(no_field_vertex);
          ierr = m_field.seed_ref_level(range_no_field_vertex,BitRefLevel().set()); CHKERRQ(ierr);
          EntityHandle lambda_meshset = m_field.get_field_meshset("LAMBDA");
          rval = m_field.get_moab().add_entities(lambda_meshset,range_no_field_vertex); CHKERRG(rval);
        }
        //this entity will carray data for this finite element
        EntityHandle meshset_fe_arc_length;
        {
          rval = moab.create_meshset(MESHSET_SET,meshset_fe_arc_length); CHKERRG(rval);
          rval = moab.add_entities(meshset_fe_arc_length,&no_field_vertex,1); CHKERRG(rval);
          ierr = m_field.seed_ref_level_MESHSET(meshset_fe_arc_length,BitRefLevel().set()); CHKERRQ(ierr);
        }
        //finally add created meshset to the ARC_LENGTH finite element
        ierr = m_field.add_ents_to_finite_element_by_MESHSET(meshset_fe_arc_length,"ARC_LENGTH",false); CHKERRQ(ierr);
      }

      //set app. order
      //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
      ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

      ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

      /*//reduce level of approximation for entities on interface
      Range prims;
      ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(problem_bit_level,BitRefLevel().set(),MBPRISM,prims); CHKERRQ(ierr);
      Range prims_faces;
      rval = m_field.get_moab().get_adjacencies(prims,2,false,prims_faces,moab::Interface::UNION); CHKERRG(rval);
      Range prims_faces_edges;
      rval = m_field.get_moab().get_adjacencies(prims_faces,1,false,prims_faces_edges,moab::Interface::UNION); CHKERRG(rval);
      ierr = m_field.set_field_order(prims_faces,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);
      ierr = m_field.set_field_order(prims_faces_edges,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);*/

      //Elements with boundary conditions
      ierr = MetaNeummanForces::addNeumannBCElements(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
      ierr = MetaNodalForces::addElement(m_field,"DISPLACEMENT");  CHKERRQ(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","PRESSURE_FE"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("FORCE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("FORCE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("FORCE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("PRESSURE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("PRESSURE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("PRESSURE_FE","LAMBDA"); CHKERRQ(ierr);
      ierr = m_field.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
      ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","BODY_FORCE"); CHKERRQ(ierr);
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
        Range tets;
        rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERRG(rval);
        ierr = m_field.add_ents_to_finite_element_by_type(tets,MBTET,"BODY_FORCE"); CHKERRQ(ierr);
      }

    }

    /****/
    //build database

    //build field
    ierr = m_field.build_fields(); CHKERRQ(ierr);
    Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

    //build finite elemnts
    ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

    //build adjacencies
    ierr = m_field.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

    /****/
    ProblemsManager *prb_mng_ptr;
    ierr = m_field.getInterface(prb_mng_ptr); CHKERRQ(ierr);
    //build problem
    ierr = prb_mng_ptr->buildProblem("ELASTIC_MECHANICS",true); CHKERRQ(ierr);
    //partition
    ierr = prb_mng_ptr->partitionProblem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("ELASTIC_MECHANICS",false,0,pcomm->size()); CHKERRQ(ierr);
    //what are ghost nodes, see Petsc Manual
    ierr = prb_mng_ptr->partitionGhostDofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

    //print bcs
    MeshsetsManager *mmanager_ptr;
    ierr = m_field.getInterface(mmanager_ptr); CHKERRQ(ierr);
    ierr = mmanager_ptr->printDisplacementSet(); CHKERRQ(ierr);
    ierr = mmanager_ptr->printForceSet(); CHKERRQ(ierr);
    //print block sets with materials
    ierr = mmanager_ptr->printMaterialsSet(); CHKERRQ(ierr);

    //create matrices
    Vec F,F_body_force,D;
    ierr = m_field.getInterface<VecManager>()->vecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
    ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
    ierr = VecDuplicate(F,&F_body_force); CHKERRQ(ierr);
    Mat Aij;
    ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

    //Assemble F and Aij
    double young_modulus=1;
    double poisson_ratio=0.0;

    boost::ptr_vector<CohesiveInterfaceElement::PhysicalEquation> interface_materials;

    //FIXME this in fact allow only for one type of interface,
    //problem is Young Modulus in interface mayoung_modulusterial
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      cout << std::endl << *it << std::endl;

      //Get block name
      string name = it->getName();

      if (name.compare(0,11,"MAT_ELASTIC") == 0) {
        Mat_Elastic mydata;
        ierr = it->getAttributeDataStructure(mydata); CHKERRQ(ierr);
        cout << mydata;
        young_modulus=mydata.data.Young;
        poisson_ratio=mydata.data.Poisson;
      } else if (name.compare(0,10,"MAT_INTERF") == 0) {
        Mat_Interf mydata;
        ierr = it->getAttributeDataStructure(mydata); CHKERRQ(ierr);
        cout << mydata;

        interface_materials.push_back(new CohesiveInterfaceElement::PhysicalEquation(m_field));
        interface_materials.back().h = 1;
        interface_materials.back().youngModulus = mydata.data.alpha;
        interface_materials.back().beta = mydata.data.beta;
        interface_materials.back().ft = mydata.data.ft;
        interface_materials.back().Gf = mydata.data.Gf;

        EntityHandle meshset = it->getMeshset();
        Range tris;
        rval = moab.get_entities_by_type(meshset,MBTRI,tris,true); CHKERRG(rval);
        Range ents3d;
        rval = moab.get_adjacencies(tris,3,false,ents3d,moab::Interface::UNION); CHKERRG(rval);
        interface_materials.back().pRisms = ents3d.subset_by_type(MBPRISM);

      }
    }

    { //FIXME
      boost::ptr_vector<CohesiveInterfaceElement::PhysicalEquation>::iterator pit = interface_materials.begin();
      for(; pit != interface_materials.end();pit++) {
        pit->youngModulus = young_modulus;
      }
    }


    ArcLengthCtx* arc_ctx = new ArcLengthCtx(m_field,"ELASTIC_MECHANICS");
    ArcLengthElement* my_arc_method_ptr = new ArcLengthElement(m_field,arc_ctx);
    ArcLengthElement& my_arc_method = *my_arc_method_ptr;
    ArcLengthSnesCtx snes_ctx(m_field,"ELASTIC_MECHANICS",arc_ctx);
    AssembleRhsVectors pre_post_proc_fe(m_field,F_body_force,arc_ctx);

    DirichletDisplacementBc my_dirichlet_bc(m_field,"DISPLACEMENT",Aij,D,F);
    ierr = m_field.get_problem("ELASTIC_MECHANICS",&my_dirichlet_bc.problemPtr); CHKERRQ(ierr);
    ierr = my_dirichlet_bc.iNitalize(); CHKERRQ(ierr);
    boost::shared_ptr<Hooke<adouble> > hooke_adouble_ptr(new Hooke<adouble>);
    boost::shared_ptr<Hooke<double> > hooke_double_ptr(new Hooke<double>);
    NonlinearElasticElement elastic(m_field,2);
    {
      int id  = 0;
      elastic.setOfBlocks[id].iD = id;
      elastic.setOfBlocks[id].E = young_modulus;
      elastic.setOfBlocks[id].PoissonRatio = poisson_ratio;
      elastic.setOfBlocks[id].materialDoublePtr = hooke_double_ptr;
      elastic.setOfBlocks[id].materialAdoublePtr = hooke_adouble_ptr;
      rval = m_field.get_moab().get_entities_by_type(
        m_field.get_finite_element_meshset("ELASTIC"),MBTET,
        elastic.setOfBlocks[id].tEts,true
      ); CHKERRG(rval);
      ierr = elastic.setOperators(
        "DISPLACEMENT","MESH_NODE_POSITIONS",false,true
      ); CHKERRQ(ierr);
    }
    CohesiveInterfaceElement cohesive_elements(m_field);
    ierr = cohesive_elements.addOps("DISPLACEMENT",interface_materials); CHKERRQ(ierr);

    PetscInt M,N;
    ierr = MatGetSize(Aij,&M,&N); CHKERRQ(ierr);
    PetscInt m,n;
    MatGetLocalSize(Aij,&m,&n);
    ArcLengthMatShell* mat_ctx = new ArcLengthMatShell(Aij,arc_ctx,"ELASTIC_MECHANICS");
    Mat ShellAij;
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)mat_ctx,&ShellAij); CHKERRQ(ierr);
    ierr = MatShellSetOperation(ShellAij,MATOP_MULT,(void(*)(void))ArcLengthMatMultShellOp); CHKERRQ(ierr);

    //body forces
    BodyForceConstantField body_forces_methods(m_field);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
      ierr = body_forces_methods.addBlock("DISPLACEMENT",F_body_force,it->getMeshsetId()); CHKERRQ(ierr);
    }
    ierr = VecZeroEntries(F_body_force); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F_body_force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_body_force,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","BODY_FORCE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F_body_force,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_body_force,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_body_force); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_body_force); CHKERRQ(ierr);

    //surface forces
    boost::ptr_map<std::string,NeummanForcesSurface> neumann_forces;
    string fe_name_str = "FORCE_FE";
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(m_field));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = neumann_forces.at(fe_name_str).addForce("DISPLACEMENT",arc_ctx->F_lambda,it->getMeshsetId());  CHKERRQ(ierr);
    }
    fe_name_str = "PRESSURE_FE";
    neumann_forces.insert(fe_name_str,new NeummanForcesSurface(m_field));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      ierr = neumann_forces.at(fe_name_str).addPreassure("DISPLACEMENT",arc_ctx->F_lambda,it->getMeshsetId()); CHKERRQ(ierr);
    }
    //add nodal forces
    boost::ptr_map<std::string,NodalForce> nodal_forces;
    fe_name_str ="FORCE_FE";
    nodal_forces.insert(fe_name_str,new NodalForce(m_field));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = nodal_forces.at(fe_name_str).addForce("DISPLACEMENT",arc_ctx->F_lambda,it->getMeshsetId());  CHKERRQ(ierr);
    }

    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,ShellAij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    KSP ksp;
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    PC pc;
    ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    PCArcLengthCtx* pc_ctx = new PCArcLengthCtx(ShellAij,Aij,arc_ctx);
    ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
    ierr = PCShellSetContext(pc,pc_ctx); CHKERRQ(ierr);
    ierr = PCShellSetApply(pc,PCApplyArcLength); CHKERRQ(ierr);
    ierr = PCShellSetSetUp(pc,PCSetupArcLength); CHKERRQ(ierr);

    //Rhs
    SnesCtx::FEMethodsSequence& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
    snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
    snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_proc_fe);
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("INTERFACE",&cohesive_elements.getFeRhs()));
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeRhs()));
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ARC_LENGTH",&my_arc_method));
    snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_proc_fe);
    snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

    //Mat
    SnesCtx::FEMethodsSequence& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
    snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("INTERFACE",&cohesive_elements.getFeLhs()));
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeLhs()));
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ARC_LENGTH",&my_arc_method));
    snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

    double gamma = 0.5,reduction = 1;
    //step = 1;
    if(step == 1) {
      step_size = step_size_reduction;
    } else {
      reduction = step_size_reduction;
      step++;
    }

    boost::ptr_map<std::string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
    ierr = VecZeroEntries(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(arc_ctx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(arc_ctx->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    for(;mit!=neumann_forces.end();mit++) {
      ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
    }
    ierr = VecGhostUpdateBegin(arc_ctx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(arc_ctx->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arc_ctx->F_lambda); CHKERRQ(ierr);
    for(std::vector<int>::iterator vit = my_dirichlet_bc.dofsIndices.begin();
    vit!=my_dirichlet_bc.dofsIndices.end();vit++) {
      ierr = VecSetValue(arc_ctx->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(arc_ctx->F_lambda); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(arc_ctx->F_lambda); CHKERRQ(ierr);
    //F_lambda2
    ierr = VecDot(arc_ctx->F_lambda,arc_ctx->F_lambda,&arc_ctx->F_lambda2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ctx->F_lambda2);

    if(step>1) {
      ierr = m_field.getInterface<VecManager>()->setLocalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field.getInterface<VecManager>()->setOtherGlobalGhostVector(
        "ELASTIC_MECHANICS","DISPLACEMENT","X0_DISPLACEMENT",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        double x0_nrm;
        ierr = VecNorm(arc_ctx->x0,NORM_2,&x0_nrm);  CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"\tRead x0_nrm = %6.4e dlambda = %6.4e\n",x0_nrm,arc_ctx->dLambda);
        ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);
      } else {
        ierr = arc_ctx->setS(0); CHKERRQ(ierr);
        ierr = arc_ctx->setAlphaBeta(0,1); CHKERRQ(ierr);
      }
      ierr = SnesRhs(snes,D,F,&snes_ctx); CHKERRQ(ierr);

      PostProcVolumeOnRefinedMesh post_proc(m_field);
      ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
      ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
      //add postpocessing for sresses
      post_proc.getOpPtrVector().push_back(
        new PostProcHookStress(
          m_field,
          post_proc.postProcMesh,
          post_proc.mapGaussPts,
          "DISPLACEMENT",
          post_proc.commonData
        )
      );

      bool converged_state  = false;
      for(;step<max_steps;step++) {

        if(step == 1) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Load Step %D step_size = %6.4e\n",step,step_size); CHKERRQ(ierr);
          ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
          ierr = arc_ctx->setAlphaBeta(0,1); CHKERRQ(ierr);
          ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
          double dlambda;
          ierr = my_arc_method.calculate_init_dlambda(&dlambda); CHKERRQ(ierr);
          ierr = my_arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
        } else if(step == 2) {
          ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);
          ierr = my_arc_method.calculate_dx_and_dlambda(D); CHKERRQ(ierr);
          ierr = my_arc_method.calculate_lambda_int(step_size); CHKERRQ(ierr);
          ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
          double dlambda = arc_ctx->dLambda;
          double dx_nrm;
          ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
            step,step_size,dlambda,dx_nrm,arc_ctx->dx2
          ); CHKERRQ(ierr);
          ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
          ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
          ierr = my_arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
        } else {
          ierr = my_arc_method.calculate_dx_and_dlambda(D); CHKERRQ(ierr);
          ierr = my_arc_method.calculate_lambda_int(step_size); CHKERRQ(ierr);
          //step_size0_1/step_size0 = step_stize1/step_size
          //step_size0_1 = step_size0*(step_stize1/step_size)
          step_size *= reduction;
          ierr = arc_ctx->setS(step_size); CHKERRQ(ierr);
          double dlambda = reduction*arc_ctx->dLambda;
          ierr = VecScale(arc_ctx->dx,reduction); CHKERRQ(ierr);
          double dx_nrm;
          ierr = VecNorm(arc_ctx->dx,NORM_2,&dx_nrm);  CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,
            "Load Step %D step_size = %6.4e dlambda0 = %6.4e dx_nrm = %6.4e dx2 = %6.4e\n",
            step,step_size,dlambda,dx_nrm,arc_ctx->dx2
          ); CHKERRQ(ierr);
          ierr = VecCopy(D,arc_ctx->x0); CHKERRQ(ierr);
          ierr = VecAXPY(D,1.,arc_ctx->dx); CHKERRQ(ierr);
          ierr = my_arc_method.set_dlambda_to_x(D,dlambda); CHKERRQ(ierr);
        }

        ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);

        //Distribute displacements on all processors
        ierr = m_field.getInterface<VecManager>()->setGlobalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = m_field.loop_finite_elements(
          "ELASTIC_MECHANICS","INTERFACE",cohesive_elements.getFeHistory(),0,m_field.get_comm_size()
        ); CHKERRQ(ierr);
        //Remove nodes of damaged prisms
        ierr = my_arc_method.remove_damaged_prisms_nodes(); CHKERRQ(ierr);

        int its;
        ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

        SNESConvergedReason reason;
        ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);

        if (reason < 0) {
          ierr = arc_ctx->setAlphaBeta(1,0); CHKERRQ(ierr);
          reduction =0.1;
          converged_state = false;
          continue;
        } else {
          if (step > 1 && converged_state) {
            reduction = pow((double)its_d/(double)(its+1),gamma);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"reduction step_size = %6.4e\n", reduction); CHKERRQ(ierr);
          }

          //Save data on mesh
          ierr = m_field.getInterface<VecManager>()->setGlobalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field.getInterface<VecManager>()->setOtherGlobalGhostVector(
            "ELASTIC_MECHANICS","DISPLACEMENT","X0_DISPLACEMENT",COL,arc_ctx->x0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            converged_state = true;
          }
          //
          if (reason > 0) {
            FILE *datafile;
            PetscFOpen(PETSC_COMM_SELF,DATAFILENAME,"a+",&datafile);
            PetscFPrintf(PETSC_COMM_WORLD,datafile,"%d %d ",reason,its);
            fclose(datafile);
            ierr = my_arc_method.postProcessLoadPath(); CHKERRQ(ierr);
          }

          if(step % 1 == 0) {

            // if(pcomm->rank()==0) {
            //   std::ostringstream sss;
            //   sss << "restart_" << step << ".h5m";
            //   rval = moab.write_file(sss.str().c_str()); CHKERRG(rval);
            // }

            ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",post_proc); CHKERRQ(ierr);
            std::ostringstream ss;
            ss << "out_values_" << step << ".h5m";
            ierr = post_proc.writeFile(ss.str().c_str()); CHKERRQ(ierr);

          }

        }

        //Save data on mesh
        ierr = m_field.getInterface<VecManager>()->setGlobalGhostVector("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

        //detroy matrices
        ierr = VecDestroy(&F); CHKERRQ(ierr);
        ierr = VecDestroy(&D); CHKERRQ(ierr);
        ierr = VecDestroy(&F_body_force); CHKERRQ(ierr);
        ierr = MatDestroy(&Aij); CHKERRQ(ierr);
        ierr = SNESDestroy(&snes); CHKERRQ(ierr);
        ierr = MatDestroy(&ShellAij); CHKERRQ(ierr);
        delete arc_ctx;
        delete mat_ctx;
        delete pc_ctx;
        delete my_arc_method_ptr;


      } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

  return 0;
}
