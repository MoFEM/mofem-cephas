using namespace MoFEM;

#include <adolc/adolc.h> 
#include <ConvectiveMassElement.hpp>

#include <ConfigurationalFractureForDynamics.hpp>

PetscErrorCode ConfigurationalFracturDynamics::coupled_dynamic_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;

  //Fields
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
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
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

  //
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

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
  bool cs = true;
  if(cs) {
    ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }
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

  ierr = iNertia.setBlocks(); CHKERRQ(ierr);
  ierr = iNertia.addConvectiveMassElement("MASS_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = iNertia.addVelocityElement("VELOCITY_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","MASS_ELEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_DYNAMIC","VELOCITY_ELEMENT"); CHKERRQ(ierr);

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

PetscErrorCode ConfigurationalFracturDynamics::solve_dynmaic_problem(FieldInterface& m_field,TS *ts,double fraction_treshold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;
  PetscBool flg;

  set_PhysicalEquationNumber(hooke);

  ierr = front_projection_data(m_field,"COUPLED_DYNAMIC"); CHKERRQ(ierr);

  //create matrices
  Mat K;
  ierr = m_field.MatCreateMPIAIJWithArrays("COUPLED_DYNAMIC",&K); CHKERRQ(ierr);
  //create vectors
  Vec F;
  ierr = m_field.VecCreateGhost("COUPLED_DYNAMIC",ROW,&F); CHKERRQ(ierr);
  Vec D;
  ierr = m_field.VecCreateGhost("COUPLED_DYNAMIC",COL,&D); CHKERRQ(ierr);

  Range corners_edges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edgesNodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edgesNodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edgesNodes.begin(),corners_edgesNodes.end());
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  Range crack_front_edges,crack_front_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  Range level_nodes;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  crack_front_nodes = intersect(crack_front_nodes,level_nodes);
  //corners_nodes.merge(subtract(level_nodes,crack_front_nodes));

  Range fix_nodes;
  ierr = fix_all_but_one(m_field,fix_nodes,fraction_treshold); CHKERRQ(ierr);
  corners_nodes.merge(fix_nodes);

  struct MyPrePostProcessFEMethod: public FEMethod {
    
    FieldInterface& mField;
    string velocityField,spatialPositionField;

    MyPrePostProcessFEMethod(FieldInterface& _m_field): 
      mField(_m_field),
      velocityField("SPATIAL_VELOCITY"),
      spatialPositionField("SPATIAL_POSITION") {}
  
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

      ierr = mField.set_other_local_VecCreateGhost(problemPtr,velocityField,"DOT_"+velocityField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = mField.set_other_local_VecCreateGhost(problemPtr,spatialPositionField,"DOT_"+spatialPositionField,COL,ts_u_t,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

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

  struct BothSurfaceConstrains: public FEMethod {

    FieldInterface& mField;
    BothSurfaceConstrains(FieldInterface& m_field): mField(m_field) {} 

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      switch(snes_ctx) {
	case CTX_SNESSETFUNCTION: { 
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	}
	break;
	case CTX_SNESSETJACOBIAN: 
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	break;
	default:
	break;
      }
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      //cerr << "AAAAAAAA\n";
      PetscErrorCode ierr;
      vector<int> lambda_dofs(9,-1);
      vector<double> lambda_vals(9,0);
      for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"LAMBDA_BOTH_SIDES",MBVERTEX,it)) {
	lambda_dofs[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] = it->get_petsc_gloabl_dof_idx();
	lambda_vals[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] = it->get_FieldData();
	//cerr << "l " 
	  //<< 3*it->side_number_ptr->side_number+it->get_EntDofIdx() << " " 
	  //<< lambda_dofs[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] << " " 
	  //<< lambda_vals[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] 
	  //<< endl;
      }
      vector<int> positions_dofs(18,-1);
      vector<double> positions_vals(18,0);
      for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"MESH_NODE_POSITIONS",MBVERTEX,it)) {
	int dd = 3*it->side_number_ptr->side_number+it->get_EntDofIdx();
	if(lambda_dofs[dd>8 ? dd -9 : dd] == -1) continue;
	positions_dofs[dd] = it->get_petsc_gloabl_dof_idx();
	positions_vals[dd] = it->get_FieldData();
	//cerr << "p " << dd << " " << positions_dofs[dd] << " " << positions_vals[dd] << endl;
      }
      //cerr << endl;
      const double alpha = 0;
      const double betha = 1e2;
      switch(snes_ctx) {
	case CTX_SNESSETFUNCTION:  
	  for(int ii = 0;ii<9;ii++) {
	    if(lambda_dofs[ii] == -1) continue;
	    double val1 = betha*(positions_vals[0+ii] - positions_vals[9+ii]) + alpha*lambda_vals[ii];
	    ierr = VecSetValue(snes_f,lambda_dofs[ii],val1,INSERT_VALUES); CHKERRQ(ierr);
	    double val2 = betha*lambda_vals[ii];
	    ierr = VecSetValue(snes_f,positions_dofs[0+ii],+val2,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = VecSetValue(snes_f,positions_dofs[9+ii],-val2,INSERT_VALUES); CHKERRQ(ierr);
	  }
	break;
	case CTX_SNESSETJACOBIAN: 
	  for(int ii = 0;ii<9;ii++) {
	    if(lambda_dofs[ii] == -1) continue;
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],positions_dofs[0+ii],+1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],positions_dofs[9+ii],-1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,positions_dofs[0+ii],lambda_dofs[ii],+1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,positions_dofs[9+ii],lambda_dofs[ii],-1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],lambda_dofs[ii],alpha,INSERT_VALUES); CHKERRQ(ierr);
	  }
	break;
	default:
	break;
      }

      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      switch(snes_ctx) {
	case CTX_SNESSETFUNCTION: { 
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	}
	break;
	case CTX_SNESSETJACOBIAN: 
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	break;
	default:
	break;
      }
      PetscFunctionReturn(0);
    }

  };

  //spatial and material forces
  const double young_modulus = 1;
  const double poisson_ratio = 0.;
  MyNonLinearSpatialElasticFEMthod fe_spatial(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_spatial.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_spatial.isCoupledProblem = true;
  MyEshelbyFEMethod fe_material(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_material.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_material.isCoupledProblem = true;
  //meshs moothing
  MyMeshSmoothingFEMethod smoother(m_field);
  ierr = smoother.initCrackFrontData(m_field); CHKERRQ(ierr);
  set_qual_ver(3);
  //constrains
  SnesConstrainSurfacGeometry constrain_body_surface(m_field,"LAMBDA_SURFACE");
  constrain_body_surface.nonlinear = true;
  SnesConstrainSurfacGeometry constrain_crack_surface(m_field,"LAMBDA_CRACK_SURFACE");
  constrain_crack_surface.nonlinear = true;
  Snes_CTgc_CONSTANT_AREA_FEMethod ct_gc(m_field,*projFrontCtx,"COUPLED_DYNAMIC","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod dct_gc(m_field,K,"LAMBDA_CRACKFRONT_AREA");
  TangentWithMeshSmoothingFrontConstrain_FEMethod tangent_constrain(m_field,&smoother,"LAMBDA_CRACK_TANGENT_CONSTRAIN");
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
  //portsproc
  MyPrePostProcessFEMethod pre_post_method(m_field);

  //TS
  TsCtx ts_ctx(m_field,"COUPLED_DYNAMIC");
  //rhs preProcess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&pre_post_method);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&fix_material_pts);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&ct_gc);

  //rhs loops
  TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->get_msId());  CHKERRQ(ierr);
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }

  ts_ctx.get_preProcess_to_do_IFunction().push_back(&fix_material_pts);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&pre_post_method);

  //lhs
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&fix_material_pts);
  TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&dct_gc));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Mat.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&fix_material_pts);
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);




  ierr = delete_surface_projection_data(m_field); CHKERRQ(ierr);

  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);



  PetscFunctionReturn(0);
}

