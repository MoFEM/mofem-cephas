using namespace MoFEM;

#include <adolc/adolc.h> 
#include <ConvectiveMassElement.hpp>
#include <TimeForceScale.hpp>
#include <PotsProcOnRefMesh.hpp>
#include <ConfigurationalFractureForDynamics.hpp>

struct MonitorRestart: public FEMethod {

  double *time;
  int *step;
  FieldInterface &mField;
  int pRT;

  MonitorRestart(FieldInterface &m_field,TS ts): 
    mField(m_field) {

    PetscErrorCode ierr;
    ErrorCode rval;
    double def_t_val = 0;

    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_time;
    rval = m_field.get_moab().tag_get_handle("_TsTime_",1,MB_TYPE_DOUBLE,th_time,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); CHKERR(rval);
      ierr = TSSetTime(ts,*time); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      rval = m_field.get_moab().tag_set_data(th_time,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); CHKERR(rval);
    }
    int zero = 0;
    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&zero); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    }

    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    //PetscErrorCode ierr;
    ErrorCode rval;
    (*time) = ts_t;
    if(pRT>0) {
      if((*step)%pRT==0) {
	ostringstream ss;
	ss << "restart_" << (*step) << ".h5m";
	rval = mField.get_moab().write_file(ss.str().c_str()); CHKERR_PETSC(rval);
      }
    }
    (*step)++;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct MonitorPostProc: public FEMethod {

  FieldInterface &mField;
  PostPocOnRefinedMesh postProc;

  bool iNit;

  int pRT;
  int *step;

  MonitorPostProc(FieldInterface &m_field): 
    FEMethod(),mField(m_field),postProc(m_field),iNit(false) { 
    
    ErrorCode rval;
    PetscErrorCode ierr;
    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_step;
    int zero = 0;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&zero); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    }


    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    if(!iNit) {
      ierr = postProc.generateRefereneElemenMesh(); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_VELOCITY"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);

      iNit = true;
    }

    if((*step)%pRT==0) {
      ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","MASS_ELEMENT",postProc); CHKERRQ(ierr);
      ostringstream sss;
      sss << "out_values_" << (*step) << ".h5m";
      rval = postProc.postProcMesh.write_file(sss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct MonitorUpdateFrezedNodes: public FEMethod {

  FieldInterface &mField;
  ConfigurationalFracturDynamics *confProb;
  FixBcAtEntities &bC;
  MyEshelby &feMaterial;
  ConvectiveMassElement &iNertia;
  TangentWithMeshSmoothingFrontConstrain &tangentConstrain;
  SnesConstrainSurfacGeometry &surfaceConstrains;
  SnesConstrainSurfacGeometry &crackSurfaceConstrains;


  MonitorUpdateFrezedNodes(FieldInterface &m_field,
    ConfigurationalFracturDynamics *conf_prob,
    FixBcAtEntities &bc,
    MyEshelby &fe_material,ConvectiveMassElement &inertia,
    TangentWithMeshSmoothingFrontConstrain &tangent_constrain,
    SnesConstrainSurfacGeometry &surface_constrains,
    SnesConstrainSurfacGeometry &crack_surface_constrains): 
    mField(m_field),confProb(conf_prob),
    bC(bc),feMaterial(fe_material),iNertia(inertia),
    tangentConstrain(tangent_constrain),
    surfaceConstrains(surface_constrains),
    crackSurfaceConstrains(crack_surface_constrains) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    Vec F_Material;
    ierr = mField.VecCreateGhost("COUPLED_DYNAMIC",ROW,&F_Material); CHKERRQ(ierr);
    ierr = VecZeroEntries(F_Material); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    surfaceConstrains.snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
    surfaceConstrains.snes_f = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","CandCT_SURFACE_ELEM",surfaceConstrains); CHKERRQ(ierr);
    crackSurfaceConstrains.snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
    crackSurfaceConstrains.snes_f = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","CandCT_CRACK_SURFACE_ELEM",crackSurfaceConstrains); CHKERRQ(ierr);
    feMaterial.snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
    feMaterial.snes_f = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","MATERIAL_COUPLED",feMaterial); CHKERRQ(ierr);
    tangentConstrain.meshFEPtr->snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
    tangentConstrain.meshFEPtr->snes_f = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","MESH_SMOOTHER",*tangentConstrain.meshFEPtr); CHKERRQ(ierr);
    tangentConstrain.snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
    tangentConstrain.snes_f = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","C_TANGENT_ELEM",tangentConstrain); CHKERRQ(ierr);
    iNertia.getLoopFeTRhs().ts_ctx = FEMethod::CTX_TSSETIFUNCTION;
    iNertia.getLoopFeTRhs().ts_F = F_Material;
    ierr = mField.loop_finite_elements("COUPLED_DYNAMIC","DYNAMIC_ESHELBY_TERM",iNertia.getLoopFeTRhs()); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_Material); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_Material); CHKERRQ(ierr);
    ierr = mField.set_other_global_ghost_vector("COUPLED_DYNAMIC","MESH_NODE_POSITIONS","MATERIAL_FORCE",
      ROW,F_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecDestroy(&F_Material); CHKERRQ(ierr);

    //caculate griffith forces
    ierr = confProb->project_force_vector(mField,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
    ierr = confProb->griffith_force_vector(mField,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
    ierr = confProb->griffith_g(mField,"COUPLED_DYNAMIC"); CHKERRQ(ierr);

    //fix nodes
    ierr = confProb->rate_of_work_of_fracture(mField); CHKERRQ(ierr);
    Range fix_nodes;
    ierr = confProb->fix_front_nodes(mField,fix_nodes); CHKERRQ(ierr);

    Range corners_edges,corners_nodes;
    ierr = mField.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
    Range corners_edges_nodes;
    ErrorCode rval;
    rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
    corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());
    corners_nodes.insert(fix_nodes.begin(),fix_nodes.end());

    bC.map_zero_rows.clear();
    bC.dofsIndices.clear();
    bC.dofsValues.clear();
    bC.eNts = corners_nodes;
    bC.problemPtr = problemPtr;
    bC.iNitalize();

    //update mesh
    ierr = confProb->set_coordinates_from_material_solution(mField,false); CHKERRQ(ierr);
    //reset material velocities, no inertaia component, thus not change RHS
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"DOT_MESH_NODE_POSITIONS",dof_ptr)) {
      dof_ptr->get_FieldData() = 0;
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct MonitorLoadPath: public FEMethod {

  FieldInterface &mField;
  TS tS;

  MonitorLoadPath(FieldInterface &m_field,TS ts): 
    FEMethod(),mField(m_field),tS(ts) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    double time;
    ierr = TSGetTime(tS,&time); CHKERRQ(ierr);

    //const MoFEMProblem *problemPtr;
    //ierr = mField.get_problem("COUPLED_DYNAMIC",&problemPtr); CHKERRQ(ierr);

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|UNKNOWNCUBITNAME,it)) {
      if(it->get_Cubit_name() != "LoadPath") continue;

      Range nodes;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
      for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
	double coords[3];
	rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problemPtr,*nit,dof)) {
	  if(dof->get_name()!="SPATIAL_POSITION") continue;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,
	      "load_path_disp ent %ld dim %d "
	      "coords ( %8.6f %8.6f %8.6f ) "
	      "val %6.4e Time %6.4e\n",
	      dof->get_ent(),dof->get_dof_rank(),
	      coords[0],coords[1],coords[2],
	      dof->get_FieldData()-coords[dof->get_dof_rank()],time); CHKERRQ(ierr);
	}
      }

    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct MyPrePostProcessDynamics: public FEMethod {
    
  FieldInterface& mField;
  string velocityField,spatialPositionField,meshPositionField;


  MyPrePostProcessDynamics(FieldInterface& _m_field): 
    mField(_m_field),
    velocityField("SPATIAL_VELOCITY"),
    spatialPositionField("SPATIAL_POSITION"),
    meshPositionField("MESH_NODE_POSITIONS") {}

  PetscErrorCode ierr;
    
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
      
    switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  snes_ctx = CTX_SNESSETFUNCTION;
	  snes_f = ts_F;
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  snes_ctx = CTX_SNESSETJACOBIAN;
	  snes_B = ts_B;
	  break;
	}
	default:
	break;
    }

    ierr = mField.set_other_local_ghost_vector(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = mField.set_other_local_ghost_vector(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = mField.set_other_local_ghost_vector(problemPtr,meshPositionField,"DOT_"+meshPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }
    
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: {
	  //snes_f
        ierr = VecGhostUpdateBegin(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
      }
      break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};


PetscErrorCode ConfigurationalFracturDynamics::coupled_dynamic_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;

  //Fields
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("GRIFFITH_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("GRIFFITH_FORCE_TANGENT",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("ELASTIC_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("MATERIAL_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("MESH_SMOOTHER",MF_ZERO); CHKERRQ(ierr);

  //fes definitions
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);

  //
  ierr = m_field.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("COUPLED_DYNAMIC",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","ELASTIC_COUPLED"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","NEUAMNN_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","MESH_SMOOTHER"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","BOTH_SIDE_OF_CRACK"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC",ss2.str()); CHKERRQ(ierr);
  }

  ierr = m_field.add_field("SPATIAL_VELOCITY",H1,3,MF_ZERO); CHKERRQ(ierr);


  PetscInt order;
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    Tag th_set_order;
    rval = m_field.get_moab().tag_get_handle("_SET_ORDER",th_set_order); CHKERR_PETSC(rval);
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    rval = m_field.get_moab().tag_get_data(th_set_order,&root_meshset,1,&order); CHKERR_PETSC(rval);
  }

  Range level_tets;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"SPATIAL_VELOCITY"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"MATERIAL_FORCE"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"GRIFFITH_FORCE"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"GRIFFITH_FORCE_TANGENT"); CHKERRQ(ierr);

  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  ierr = m_field.add_field("SPATIAL_VELOCITY",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = m_field.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MATERIAL_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MATERIAL_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MATERIAL_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MATERIAL_FORCE",1); CHKERRQ(ierr);

  ierr = m_field.add_field("GRIFFITH_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"GRIFFITH_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"GRIFFITH_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"GRIFFITH_FORCE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"GRIFFITH_FORCE",1); CHKERRQ(ierr);

  ierr = m_field.add_field("GRIFFITH_FORCE_TANGENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"GRIFFITH_FORCE_TANGENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"GRIFFITH_FORCE_TANGENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"GRIFFITH_FORCE_TANGENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"GRIFFITH_FORCE_TANGENT",1); CHKERRQ(ierr);

  ierr = m_field.add_field("DOT_SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_POSITION",1); CHKERRQ(ierr);

  ierr = m_field.add_field("DOT_SPATIAL_VELOCITY",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_VELOCITY",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  ierr = m_field.add_field("DOT_MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"DOT_MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  ierr = iNertia.setBlocks(); CHKERRQ(ierr);
  ierr = iNertia.addConvectiveMassElement("MASS_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION",
    "MESH_NODE_POSITIONS",true,*ptr_bit_level0); CHKERRQ(ierr);
  ierr = iNertia.addVelocityElement("VELOCITY_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION",
    "MESH_NODE_POSITIONS",true,*ptr_bit_level0); CHKERRQ(ierr);

  Range crack_front_edges,crack_front_nodes,crack_front_nodes_tets;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  rval = m_field.get_moab().get_adjacencies(crack_front_nodes,3,false,crack_front_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  ierr = iNertia.addEshelbyDynamicMaterialMomentum("DYNAMIC_ESHELBY_TERM","SPATIAL_VELOCITY","SPATIAL_POSITION",
    "MESH_NODE_POSITIONS",true,*ptr_bit_level0,&crack_front_nodes_tets); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","MASS_ELEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","VELOCITY_ELEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","DYNAMIC_ESHELBY_TERM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFracturDynamics::coupled_dynamic_partition_problems(FieldInterface& m_field) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //partition
  ierr = m_field.partition_problem("COUPLED_DYNAMIC"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("COUPLED_DYNAMIC"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("COUPLED_DYNAMIC"); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFracturDynamics::fix_front_nodes(FieldInterface& m_field,Range &fix_nodes,const double treshhold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  ErrorCode rval;
  Tag th_freez;
  const int def_order = 0;
  rval = m_field.get_moab().tag_get_handle("FROZEN_NODE",1,MB_TYPE_INTEGER,th_freez,MB_TAG_CREAT|MB_TAG_SPARSE,&def_order); CHKERR_PETSC(rval);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"freeze front nodes:\n");
  int nb_fix_nodes = 0;
  for(
    map<EntityHandle,double>::iterator mit = map_ent_g.begin();
    mit!=map_ent_g.end();mit++) {
    double g = mit->second;
    double rate = map_ent_rate_work[mit->first];
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "front node = %ld g/g_c = %4.3g rate = %4.3g\n",
      mit->first,mit->second/gc,rate); CHKERRQ(ierr);
    if( ((g/gc + treshhold) < 1) || (rate<0) || (g != g)) {
      fix_nodes.insert(mit->first);
      nb_fix_nodes++;
    }
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "nb. of fixed nodes %d\n\n",nb_fix_nodes); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFracturDynamics::rate_of_work_of_fracture(FieldInterface& m_field) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;
  map_ent_rate_work.clear();
  ublas::vector<double> griffith_force(3);
  ublas::vector<double> material_veloc(3);
  ublas::vector<double> coords(3);
  const MoFEMProblem *problem_ptr;
  ierr = m_field.get_problem("C_CRACKFRONT_MATRIX",&problem_ptr); CHKERRQ(ierr);
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problem_ptr,"LAMBDA_CRACKFRONT_AREA",diit)) {
    EntityHandle ent = diit->get_ent();
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(m_field,"GRIFFITH_FORCE",diit->get_ent(),diiit)) {
      griffith_force[diiit->get_dof_rank()] = diiit->get_FieldData();
    }
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(m_field,"DOT_MESH_NODE_POSITIONS",diit->get_ent(),diiit)) {
      material_veloc[diiit->get_dof_rank()] = diiit->get_FieldData();
    }
    rval = m_field.get_moab().get_coords(&ent,1,&coords[0]); CHKERR_PETSC(rval);
    double rate = inner_prod(griffith_force,material_veloc);
    map_ent_rate_work[ent] = rate;
    ostringstream ss;
    ss << "rate of fracture work at ent ";
    ss << setw(5) << ent;
    //coords
    ss << "\tcoords";
    ss << " " << setw(10) << setprecision(4) << coords[0];
    ss << " " << setw(10) << setprecision(4) << coords[1];
    ss << " " << setw(10) << setprecision(4) << coords[2];
    ss << " " << " rate of work at time step " << setprecision(4)  << rate;
    ss << " " << " griffith force " << setw(10) << setprecision(4) << griffith_force;
    ss << " " << " material veloc. " << setw(10) << setprecision(4) << material_veloc;
    ss << endl;
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFracturDynamics::solve_dynmaic_problem(FieldInterface& m_field,TS ts,double fraction_treshold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;
  //PetscBool flg;

  set_PhysicalEquationNumber(hooke);

  //create matrices
  Mat K;
  ierr = m_field.MatCreateMPIAIJWithArrays("COUPLED_DYNAMIC",&K); CHKERRQ(ierr);
  //create vectors
  Vec F;
  ierr = m_field.VecCreateGhost("COUPLED_DYNAMIC",ROW,&F); CHKERRQ(ierr);
  Vec D;
  ierr = m_field.VecCreateGhost("COUPLED_DYNAMIC",COL,&D); CHKERRQ(ierr);

  //create projection matrices
  ierr = front_projection_data(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);
  ierr = surface_projection_data(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);

  //spatial and material forces
  const double young_modulus = 1;
  const double poisson_ratio = 0.;
  MyNonLinearSpatialElastic fe_spatial(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_spatial.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_spatial.isCoupledProblem = true;
  MyEshelby fe_material(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_material.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_material.isCoupledProblem = true;
  //meshs moothing
  MyMeshSmoothing smoother(m_field);
  ierr = smoother.initCrackFrontData(m_field); CHKERRQ(ierr);
  set_qual_ver(0);
  //constrains
  SnesConstrainSurfacGeometry constrain_body_surface(m_field,"LAMBDA_SURFACE");
  constrain_body_surface.nonlinear = true;
  SnesConstrainSurfacGeometry constrain_crack_surface(m_field,"LAMBDA_CRACK_SURFACE");
  constrain_crack_surface.nonlinear = true;
  Snes_CTgc_CONSTANT_AREA ct_gc(m_field,projFrontCtx,"COUPLED_DYNAMIC","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA dct_gc(m_field,K,"LAMBDA_CRACKFRONT_AREA");
  TangentWithMeshSmoothingFrontConstrain tangent_constrain(m_field,&smoother,"LAMBDA_CRACK_TANGENT_CONSTRAIN");
  map<int,SnesConstrainSurfacGeometry*> other_body_surface_constrains;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    other_body_surface_constrains[msId] = new SnesConstrainSurfacGeometry(m_field,ss.str());
    other_body_surface_constrains[msId]->nonlinear = true;
  }
  //bothsieds constrains
  BothSurfaceConstrains  both_sides_constrains(m_field);
  //dirichlet constrains
  Range corners_nodes;
  {
    Range corners_edges;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
    ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
    Range corners_edges_nodes;
    rval = m_field.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
    corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());
  }
  FixBcAtEntities fix_material_pts(m_field,"MESH_NODE_POSITIONS",corners_nodes);
  fix_material_pts.fieldNames.push_back("LAMBDA_SURFACE");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_SURFACE");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_TANGENT_CONSTRAIN");
  fix_material_pts.fieldNames.push_back("LAMBDA_BOTH_SIDES");
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"SPATIAL_POSITION",K,D,F);
  my_dirichlet_bc.fixFields.push_back("SPATIAL_VELOCITY");
  my_dirichlet_bc.fixFields.push_back("MESH_NODE_POSITIONS");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_SURFACE");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_BOTH_SIDES");
  //boundary conditions
  Tag th_scale;
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_scale); CHKERR_PETSC(rval);
  double *force_scale;
  rval = m_field.get_moab().tag_get_by_ptr(th_scale,&root_meshset,1,(const void**)&force_scale); CHKERR_PETSC(rval);

  NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,K,F);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_forces = neumann_forces.getLoopSpatialFe();
  //fe_forces.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = fe_forces.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = fe_forces.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  fe_forces.methodsOp.push_back(new TimeForceScale());
  //portsproc
  MyPrePostProcessDynamics pre_post_method(m_field);

  //iNertia
  ierr = iNertia.setConvectiveMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",true,true); CHKERRQ(ierr);
  ierr = iNertia.setVelocityOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",true); CHKERRQ(ierr);

  Range crack_front_edges;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  crack_front_edges_nodes.merge(crack_front_edges);
  ierr = iNertia.setKinematicEshelbyOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",&crack_front_edges_nodes); CHKERRQ(ierr);

  /*//this is for debugging only
  struct StopAndViewMatrix: public FEMethod {
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  ierr = MatAssemblyBegin(ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	  MatView(ts_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	  std::string wait;
	  std::cin >> wait;
	  break;
	}
	default:
	break;
      }
      PetscFunctionReturn(0);
    }
  };
  StopAndViewMatrix stop_and_view;*/

  //TS
  TsCtx ts_ctx(m_field,"COUPLED_DYNAMIC");
  //rhs preProcess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&pre_post_method);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&fix_material_pts);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&ct_gc);

  //rhs loops
  TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->get_msId());  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("VELOCITY_ELEMENT",&iNertia.getLoopFeVelRhs()));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&iNertia.getLoopFeMassRhs()));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("DYNAMIC_ESHELBY_TERM",&iNertia.getLoopFeTRhs()));

  ts_ctx.get_postProcess_to_do_IFunction().push_back(&fix_material_pts);
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&pre_post_method);

  //lhs
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&fix_material_pts);
  TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&dct_gc));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("VELOCITY_ELEMENT",&iNertia.getLoopFeVelLhs()));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&iNertia.getLoopFeMassLhs()));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("DYNAMIC_ESHELBY_TERM",&iNertia.getLoopFeTLhs()));

  //lh postproc
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&fix_material_pts);
  //ts_ctx.get_postProcess_to_do_IJacobian().push_back(&stop_and_view);

  //monitor
  MonitorRestart monitor_restart(m_field,ts);
  MonitorUpdateFrezedNodes monitor_fix_matrial_nodes(
    m_field,this,fix_material_pts,fe_material,iNertia,tangent_constrain,
    constrain_body_surface,constrain_crack_surface);
  MonitorPostProc post_proc(m_field);
  MonitorLoadPath monitor_load_path(m_field,ts);

  TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&monitor_restart));
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&monitor_fix_matrial_nodes));
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&monitor_load_path));
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&post_proc));

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,K,K,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

  ierr = m_field.set_local_ghost_vector("COUPLED_DYNAMIC",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //fe_spatial.snes_ctx = FEMethod::CTX_SNESSETFUNCTION;
  //fe_spatial.snes_f = F;
  //fe_spatial.snes_B = K;
  //ierr = m_field.loop_finite_elements("COUPLED_DYNAMIC","ELASTIC_COUPLED",fe_spatial); CHKERRQ(ierr);

  /*const MoFEMProblem *problem_ptr;
  ierr = m_field.get_problem("COUPLED_DYNAMIC",&problem_ptr); CHKERRQ(ierr);
  const NumeredDofMoFEMEntity *dof_ptr;
  ierr = problem_ptr->get_col_dofs_by_petsc_gloabl_dof_idx(2947,&dof_ptr); CHKERRQ(ierr);
  cerr << *dof_ptr << endl;*/

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,D); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
  
  ierr = TSSolve(ts,D); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  ierr = delete_surface_projection_data(m_field); CHKERRQ(ierr);
  ierr = delete_surface_projection_data(m_field); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);



  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

