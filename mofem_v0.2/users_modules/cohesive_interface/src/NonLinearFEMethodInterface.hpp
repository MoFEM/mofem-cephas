/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

/** \brief Inteface element for damage with linear cohesive law
 *
 * Function, make a loop for all gauss points, and calculate gap ( separation
 * of interface ). We have three types of shape functions, for nodes, edges and
 * face of interface itself.  Values of shepe functions, for each gauss pt, are
 * stored in matrixes, nodeNTRI, _H1edgeN_, _H1edgeN_, _H1faceN_, for nodes,
 * edges and faces, respectively. 
 *
 * Since interface has to faces, top and bottom, which can have different
 * approximation order, shape functions for top and button surface are used
 * first to calculate displacements. Subtracting displacements from top and
 * button surface gap at each Gauss point is calculated. 
 *
 * Gap is stored in
 * vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap;
 * which is vector of vectors, size of vector is equal to number of Gauss
 * points, each element of that vector, is vector dimension 3, which represent
 * gap. 
 *
 * Gap at each Gauss pt, in local coordinates is in
 * vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap_loc;
 *
 * Having calculated gap, based on physical equation for interface, cohesive
 * forces are calculated. In particular,  RhsInt() function calculates vector of
 * internal forces . In LhsInt() { ... }, tangent element stiffness matrix is
 * calculated.
 *
 */
struct NonLinearInterfaceFEMethod: public ObosleteUsersModules::InterfaceFEMethod {

  enum interface_materials_context { CTX_INTLINEARSOFTENING, CTX_INTBILINEARSOFTENING, CTX_INTNONE };
  interface_materials_context intMatCtx;
  double h,beta,ft,Gf,E0,g0,kappa1;
  enum interface_context { CTX_KAPPAUPDATE = 1,  CTX_INTERFACENONE = 2 };
  interface_context ctxInt;

  Tag thDamagedPrism;

  NonLinearInterfaceFEMethod(
      FieldInterface& _mField,Mat _Aij,Vec _X,Vec _F,
      double _young_modulus, double _h,double _beta,double _ft,double _Gf, string _field_name,
      interface_materials_context _intMatCtx = CTX_INTLINEARSOFTENING): 
      InterfaceFEMethod(_mField,_Aij,_X,_F,_young_modulus,_field_name),
      intMatCtx(_intMatCtx),h(_h),beta(_beta),ft(_ft),Gf(_Gf),ctxInt(CTX_INTERFACENONE) {

    E0 = youngModulus/h;
    g0 = ft/E0;
    kappa1 = 2*Gf/ft;

    double def_damaged = 0;
    rval = mField.get_moab().tag_get_handle(
      "DAMAGED_PRISM",1,MB_TYPE_INTEGER,thDamagedPrism,MB_TAG_CREAT|MB_TAG_SPARSE,&def_damaged); CHKERR_THROW(rval);
    
  };

  ublas::vector<FieldData> g;

  /** \brief calculate gap norm
   *
   */
  virtual PetscErrorCode Calc_g() {
    PetscFunctionBegin;
    int g_dim = g_NTRI.size()/3;
    g.resize(g_dim);
    switch (intMatCtx) {
      case CTX_INTBILINEARSOFTENING:
      case CTX_INTLINEARSOFTENING: {
	for(int gg = 0;gg<g_dim;gg++) {
	  double g2 = pow(gap_loc[gg][0],2)+beta*(pow(gap_loc[gg][1],2)+pow(gap_loc[gg][1],2));
	  g[gg] = sqrt( g2 );
	}
      } break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  /** \brief calculate local and global material elastic stiffnes matrix for inteface
   *
   */
  virtual PetscErrorCode CalcDglob(const double _omega_) {
    PetscFunctionBegin;
    double E = (1-_omega_)*E0;
    ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
    Dloc(0,0) = E;
    Dloc(1,1) = E;
    Dloc(2,2) = E;
    Dglob = prod( Dloc, R );
    Dglob = prod( trans(R), Dglob );
    PetscFunctionReturn(0);
  }
  
  /** \brief calculate damage variable
   *
   */
  virtual PetscErrorCode Calc_omega(const double _kappa_,double& _omega_) {
    PetscFunctionBegin;
    switch (intMatCtx) {
      case CTX_INTLINEARSOFTENING: {
      _omega_ = 0;
      if(_kappa_>=kappa1) {
	_omega_ = 1;
	PetscFunctionReturn(0);
      } else if(_kappa_>0) {
	double a = (2.0*Gf*E0+ft*ft)*_kappa_;
	double b = (ft+E0*_kappa_)*Gf;
	_omega_ = 0.5*a/b;
      }
      } break;
      case CTX_INTBILINEARSOFTENING: {
	SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      } break;
      default:
	 SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  
  /** \brief calculate local and global material elastic stiffnes matrix for inteface
   *
   */
  virtual PetscErrorCode CalcTangetDglob(const double _omega_,const double _g_,const ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> >& _gap_loc_) {
    PetscFunctionBegin;
    switch (intMatCtx) {
      case CTX_INTLINEARSOFTENING: {
      double d_omega_ = 
	0.5*(2*Gf*E0+ft*ft)/((ft+(_g_-ft/E0)*E0)*Gf) - 0.5*((_g_-ft/E0)*(2*Gf*E0+ft*ft)*E0)/(pow(ft+(_g_-ft/E0)*E0,2)*Gf);
      ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
      Dloc(0,0) = (1-_omega_)*E0 - d_omega_*E0*_gap_loc_[0]*_gap_loc_[0]/_g_;
      Dloc(0,1) = -d_omega_*E0*_gap_loc_[0]*beta*_gap_loc_[1]/_g_;
      Dloc(0,2) = -d_omega_*E0*_gap_loc_[0]*beta*_gap_loc_[2]/_g_;
      //
      Dloc(1,0) = -d_omega_*E0*_gap_loc_[1]*_gap_loc_[0]/_g_;
      Dloc(1,1) = (1-_omega_)*E0 - d_omega_*E0*_gap_loc_[1]*beta*_gap_loc_[1]/_g_;
      Dloc(1,2) = -d_omega_*E0*_gap_loc_[1]*beta*_gap_loc_[2]/_g_;
      //
      Dloc(2,0) = -d_omega_*E0*_gap_loc_[2]*_gap_loc_[0]/_g_;
      Dloc(2,1) = -d_omega_*E0*_gap_loc_[2]*beta*_gap_loc_[1]/_g_;
      Dloc(2,2) = (1-_omega_)*E0 - d_omega_*E0*_gap_loc_[2]*beta*_gap_loc_[2]/_g_;
      Dglob = prod( Dloc, R );
      Dglob = prod( trans(R), Dglob );
      } break;
      case CTX_INTBILINEARSOFTENING: {
	SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      } break;
      default:
	 SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  /** \brief calculate internal force vector
   *
   */
  PetscErrorCode RhsInt() {
    PetscFunctionBegin;
    int g_dim = g_NTRI.size()/3;
    bool is_fully_damaged = true;
    for(int gg = 0;gg<g_dim;gg++) {
	double _kappa_ = fmax(g[gg]-g0,kappa[gg]);
	switch(ctxInt) {
	  case CTX_KAPPAUPDATE: {
	      double _omega_ = 0;
	      ierr = Calc_omega(_kappa_,_omega_); CHKERRQ(ierr);
	      if(_omega_ < 1.) is_fully_damaged = false;
	      kappa[gg] = _kappa_;
	    }
	    break;
	  default: {
	    double _omega_ = 0;
	    ierr = Calc_omega(_kappa_,_omega_); CHKERRQ(ierr);
	    //Dglob
	    ierr = CalcDglob(_omega_); CHKERRQ(ierr);
	    //Traction
	    ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > traction;
	    traction = prod(Dglob,gap[gg]);
	    if(traction.size()!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    double w = area3*G_TRI_W[gg];
	    for(int rr = 0;rr<row_mat;rr++) {
	      ublas::matrix<FieldData> &N = (rowNMatrices[rr])[gg];
	      ublas::vector<FieldData> f_int = prod(trans(N),w*traction);
	      if(RowGlob[rr].size()==0) continue;
	      if(RowGlob[rr].size()!=f_int.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(snes_f,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }}
	    break;
	}
    }

    if(ctxInt == CTX_KAPPAUPDATE) {

      if(is_fully_damaged) {

	EntityHandle ent = fePtr->get_ent();
	int set_prism_as_demaged = 1;
	rval = mField.get_moab().tag_set_data(thDamagedPrism,&ent,1,&set_prism_as_demaged); CHKERR_PETSC(rval);

      }

    }
    PetscFunctionReturn(0);
  }

  /** calculate elemement tangent matrix
   *
   */
  PetscErrorCode LhsInt() {
    PetscFunctionBegin;
    int g_dim = g_NTRI.size()/3;
    K.resize(row_mat,col_mat);
    for(int rr = 0;rr<row_mat;rr++) {
	if(RowGlob[rr].size()==0) continue;
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	    ///K matrices
	    if(gg == 0) {
	      K(rr,cc) = ublas::zero_matrix<FieldData>(row_Mat.size2(),col_Mat.size2());
	    }
	    double _kappa_ = fmax(g[gg]-g0,kappa[gg]);
	    double _omega_ = 0;
	    ierr = Calc_omega(_kappa_,_omega_); CHKERRQ(ierr);
	    //Dglob
	    if((_kappa_ <= kappa[gg])||(_kappa_>=kappa1)||(iter <= 1)) {
	      ierr = CalcDglob(_omega_); CHKERRQ(ierr);
	    } else {
	      ierr = CalcTangetDglob(_omega_,g[gg],gap_loc[gg]); CHKERRQ(ierr);
	    }
	    double w = area3*G_TRI_W[gg];
	    ublas::matrix<FieldData> NTD = prod( trans(row_Mat), w*Dglob );
	    K(rr,cc) += prod(NTD , col_Mat ); 
	  }
	}
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
    }
    PetscFunctionReturn(0);
  }

  Tag th_kappa;
  const void* tag_data_kappa[1];
  double* kappa;

  PetscErrorCode set_ctxInt(interface_context _ctx) {
    PetscFunctionBegin;
    ctxInt = _ctx;
    PetscFunctionReturn(0);
  }

  int iter;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

    //Rotation matrix
    ierr = CalcR(); CHKERRQ(ierr);
    //Calculate Matrices
    ierr = Matrices();    CHKERRQ(ierr);
    //Calcualte gap
    ierr = Calc_gap(); CHKERRQ(ierr);
    ierr = Calc_g(); CHKERRQ(ierr);

    //History
    int g_dim = g_NTRI.size()/3;
    EntityHandle fe_ent = fePtr->get_ent();
    vector<double> def_kappa(g_dim,0);
    rval = mField.get_moab().tag_get_handle("_KAPPA",g_dim,MB_TYPE_DOUBLE,th_kappa,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_EXCL,&def_kappa[0]);  
    if(rval==MB_ALREADY_ALLOCATED) {
      rval = MB_SUCCESS;
    } else {
      rval = mField.get_moab().tag_set_data(th_kappa,&fe_ent,1,&def_kappa[0]); 
    }
    CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_by_ptr(th_kappa,&fe_ent,1,tag_data_kappa); CHKERR_PETSC(rval);
    kappa = (double*)tag_data_kappa[0];

    switch(snes_ctx) {
      break;
      case CTX_SNESNONE: 
      case CTX_SNESSETFUNCTION: { 
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: { 
	ierr = SNESGetIterationNumber(snes,&iter); CHKERRQ(ierr);
	ierr = LhsInt(); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

struct ArcLengthIntElemFEMethod: public FEMethod {
  Interface& mOab;
  ErrorCode rval;
  PetscErrorCode ierr;

  ObosleteUsersModules::ArcLengthCtx* arcPtr;
  Vec GhostDiag,GhostLambdaInt;
  Range Faces3,Faces4;
  Range Edges3,Edges4;
  Range Nodes3,Nodes4;

  Tag thDamagedPrism;

  ArcLengthIntElemFEMethod(Interface& moab,
    ObosleteUsersModules::ArcLengthCtx *arcptr): 
    FEMethod(),mOab(moab),arcPtr(arcptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mOab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostDiag);
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostLambdaInt);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostLambdaInt);
    }
    Range prisms;
    rval = mOab.get_entities_by_type(0,MBPRISM,prisms,false); CHKERR_THROW(rval);
    for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
      EntityHandle f3,f4;
      rval = mOab.side_element(*pit,2,3,f3); CHKERR_THROW(rval);
      rval = mOab.side_element(*pit,2,4,f4); CHKERR_THROW(rval);
      Faces3.insert(f3);
      Faces4.insert(f4);
    }

    rval = mOab.get_adjacencies(Faces3,1,false,Edges3); CHKERR_THROW(rval);
    rval = mOab.get_adjacencies(Faces4,1,false,Edges4); CHKERR_THROW(rval);
    rval = mOab.get_connectivity(Faces3,Nodes3,true); CHKERR_THROW(rval);
    rval = mOab.get_connectivity(Faces4,Nodes4,true); CHKERR_THROW(rval);
    //Faces3.insert(Edges3.begin(),Edges3.end());
    Faces3.insert(Nodes3.begin(),Nodes3.end());
    //Faces4.insert(Edges4.begin(),Edges4.end());
    Faces4.insert(Nodes4.begin(),Nodes4.end());

    double def_damaged = 0;
    rval = mOab.tag_get_handle(
      "DAMAGED_PRISM",1,MB_TYPE_INTEGER,thDamagedPrism,MB_TAG_CREAT|MB_TAG_SPARSE,&def_damaged); CHKERR_THROW(rval);

  }
  ~ArcLengthIntElemFEMethod() {
    VecDestroy(&GhostDiag);
    VecDestroy(&GhostLambdaInt);
  }


  /** \brief remove nodes of prims which are fully damaged
    * 
    */
  PetscErrorCode remove_damaged_prisms_nodes() {
    PetscFunctionBegin;
    Range prisms;
    rval = mOab.get_entities_by_type(0,MBPRISM,prisms,false); CHKERR_PETSC(rval);
    vector<int> is_prism_damaged(prisms.size());
    rval = mOab.tag_get_data(thDamagedPrism,prisms,&*is_prism_damaged.begin()); CHKERR_PETSC(rval);
    Range::iterator pit = prisms.begin();
    vector<int>::iterator vit = is_prism_damaged.begin();
    for(;pit!=prisms.end();pit++,vit++) {
      if(*vit>0) {
	Range nodes;
	rval = mOab.get_connectivity(&*pit,1,nodes,true); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
	  Faces3.erase(*nit);
	  Faces4.erase(*nit);
	}
      }
    }
    PetscFunctionReturn(0);
  }

  double lambda_int;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = calculate_dx_and_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calculate_db(); CHKERRQ(ierr);
	ierr = calculate_lambda_int(lambda_int); CHKERRQ(ierr);
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_lambda_int(double &_lambda_int_) {
    PetscFunctionBegin;
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mOab,MYPCOMM_INDEX);
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = problemPtr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_dit = problemPtr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(problemPtr->get_nb_local_dofs_row());
    double *array;
    double *array_int_lambda;
    ierr = VecZeroEntries(GhostLambdaInt); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(arcPtr->dx,&array); CHKERRQ(ierr);
    ierr = VecGetArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    array_int_lambda[0] = 0;
    for(;dit!=hi_dit;dit++) {
      if(dit->get_ent_type() != MBVERTEX) continue;
      if(pcomm->rank() != dit->get_part()) continue;
      if(Nodes3.find(dit->get_ent())!=Nodes3.end()) {
	array_int_lambda[0] += array[dit->petsc_local_dof_idx];
      }
      if(Nodes4.find(dit->get_ent())!=Nodes4.end()) {
	array_int_lambda[0] -= array[dit->petsc_local_dof_idx];
      }
    }
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,
    //  "array_int_lambda[0] = %6.4ee\n",array_int_lambda[0]);
    ierr = VecRestoreArray(arcPtr->dx,&array); CHKERRQ(ierr);
    ierr = VecRestoreArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    _lambda_int_ = arcPtr->alpha*array_int_lambda[0] + arcPtr->dlambda*arcPtr->beta*sqrt(arcPtr->F_lambda2);
    /*PetscSynchronizedPrintf(PETSC_COMM_WORLD,
      "array_int_lambda[0] = %6.4e arcPtr->F_lambda2 = %6.4e\n",
      array_int_lambda[0],arcPtr->F_lambda2);
    PetscSynchronizedFlush(PETSC_COMM_WORLD);*/
    ierr = VecRestoreArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode calculate_db() {
    PetscFunctionBegin;
    ierr = VecZeroEntries(arcPtr->db); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(arcPtr->db,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(arcPtr->db,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = problemPtr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_dit = problemPtr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(
      problemPtr->get_nb_local_dofs_row()+problemPtr->get_nb_ghost_dofs_row());
    double *array;
    ierr = VecGetArray(arcPtr->db,&array); CHKERRQ(ierr);
    for(;dit!=hi_dit;dit++) {
      if(dit->get_ent_type() != MBVERTEX) {
	array[dit->petsc_local_dof_idx] = 0;
	continue;
      }
      if(Nodes3.find(dit->get_ent())!=Nodes3.end()) {
	array[dit->petsc_local_dof_idx] = +arcPtr->alpha;
      } else if(Nodes4.find(dit->get_ent())!=Nodes4.end()) {
	  array[dit->petsc_local_dof_idx] = -arcPtr->alpha;
      } else array[dit->petsc_local_dof_idx] = 0;
    }
    ierr = VecRestoreArray(arcPtr->db,&array); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: {
	//calculate residual for arc length row
	arcPtr->res_lambda = lambda_int - arcPtr->s;
	ierr = VecSetValue(snes_f,arcPtr->get_petsc_gloabl_dof_idx(),arcPtr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e lambda_int = %6.4e s = %6.4e\n",
	  arcPtr->res_lambda,lambda_int,arcPtr->s);
      }
      break; 
      case CTX_SNESSETJACOBIAN: {
	//calculate diagonal therm
	double diag = arcPtr->beta*sqrt(arcPtr->F_lambda2);
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(snes_B,arcPtr->get_petsc_gloabl_dof_idx(),arcPtr->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
      break;
    }	
    
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN: {
	ierr = VecAssemblyBegin(GhostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(GhostDiag,&diag); CHKERRQ(ierr);
	arcPtr->diag = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arcPtr->diag);
	ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_dx_and_dlambda(Vec &x) {
    PetscFunctionBegin;
    //dx
    ierr = VecCopy(x,arcPtr->dx); CHKERRQ(ierr);
    ierr = VecAXPY(arcPtr->dx,-1,arcPtr->x0); CHKERRQ(ierr);
    //if LAMBDA dof is on this partition
    if(arcPtr->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(arcPtr->dx,&array); CHKERRQ(ierr);
      arcPtr->dlambda = array[arcPtr->get_petsc_local_dof_idx()];
      array[arcPtr->get_petsc_local_dof_idx()] = 0;
      ierr = VecRestoreArray(arcPtr->dx,&array); CHKERRQ(ierr);
    }
    //brodcast dlambda
    int part = arcPtr->get_part();
    MPI_Bcast(&(arcPtr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
    //calculate dx2 (dot product)
    ierr = VecDot(arcPtr->dx,arcPtr->dx,&arcPtr->dx2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tdlambda = %6.4e dx2 = %6.4e\n",arcPtr->dlambda,arcPtr->dx2);
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_init_dlambda(double *dlambda) {

      PetscFunctionBegin;

      *dlambda = arcPtr->s/(arcPtr->beta*sqrt(arcPtr->F_lambda2));
      PetscPrintf(PETSC_COMM_WORLD,"\tInit dlambda = %6.4e s = %6.4e beta = %6.4e F_lambda2 = %6.4e\n",*dlambda,arcPtr->s,arcPtr->beta,arcPtr->F_lambda2);
      double a = *dlambda;
      if(a - a != 0) {
	ostringstream sss;
	sss << "s " << arcPtr->s << " " << arcPtr->beta << " " << arcPtr->F_lambda2;
	SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
      }

      PetscFunctionReturn(0);
  }

  PetscErrorCode set_dlambda_to_x(Vec &x,double dlambda) {
      PetscFunctionBegin;

      if(arcPtr->get_petsc_local_dof_idx()!=-1) {
	double *array;
	ierr = VecGetArray(x,&array); CHKERRQ(ierr);
	double lambda_old = array[arcPtr->get_petsc_local_dof_idx()];
	if(!(dlambda == dlambda)) {
	  ostringstream sss;
	  sss << "s " << arcPtr->s << " " << arcPtr->beta << " " << arcPtr->F_lambda2;
	  SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
	}
	array[arcPtr->get_petsc_local_dof_idx()] = lambda_old + dlambda;
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	  lambda_old, array[arcPtr->get_petsc_local_dof_idx()], dlambda);
	ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
  }

};
