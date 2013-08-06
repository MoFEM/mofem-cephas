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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ElasticFEMethodForInterface.hpp"
#include "ArcLeghtTools.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <petsc-private/snesimpl.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

struct ArcInterfaceElasticFEMethod: public InterfaceElasticFEMethod {

  ArcInterfaceElasticFEMethod(Interface& _moab): InterfaceElasticFEMethod(_moab) {};

  ArcLenghtCtx *arc_ptr;
  ArcInterfaceElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,
      ArcLenghtCtx *_arc_ptr): 
      InterfaceElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2,_SideSet3)
      ,arc_ptr(_arc_ptr) {};

  PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction2(3);
      traction2[0] = 0;
      traction2[1] = +1;
      traction2[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(arc_ptr->F_lambda,traction2,SideSet2); CHKERRQ(ierr);

      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction3(3);
      traction3[0] = 0;
      traction3[1] = -1;
      traction3[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(arc_ptr->F_lambda,traction3,SideSet3); CHKERRQ(ierr);

      PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
    // See FEAP - - A Finite Element Analysis Program
    D_lambda = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<3;rr++) {
	ublas::matrix_row<ublas::matrix<FieldData> > row_D_lambda(D_lambda,rr);
	for(int cc = 0;cc<3;cc++) {
	  row_D_lambda[cc] = 1;
	}
    }
    D_mu = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<6;rr++) {
	D_mu(rr,rr) = rr<3 ? 2 : 1;
    }
    D = lambda*D_lambda + mu*D_mu;


    switch(snes_ctx) {
      case ctx_SNESNone: 
      case ctx_SNESSetFunction: { 
	//F_lambda
	ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	Diagonal = PETSC_NULL;
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);
      DirihletBC.resize(0);

      switch(snes_ctx) {
	case ctx_SNESNone: {
	}
	break;
	case ctx_SNESSetJacobian: 
	case ctx_SNESSetFunction: { 
	  //Dirihlet Boundary Condition
	  ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	  if(Diagonal!=PETSC_NULL) {
	    if(DirihletBC.size()>0) {
	      DirihletBCDiagVal.resize(DirihletBC.size());
	      fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	      ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	    }
	  }
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      switch(snes_ctx) {
	case ctx_SNESNone: 
	case ctx_SNESSetFunction: {
	  //Assembly  F
	  ierr = Fint(F); CHKERRQ(ierr);
	  //Neumann Boundary Conditions
	  ierr = NeumannBC(); CHKERRQ(ierr);
	}
	break;
	case ctx_SNESSetJacobian: {
	  //Assembly  F
	  ierr = Lhs(); CHKERRQ(ierr);
	}
	break;
	default:
	  SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
  }


  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone: {
      }
      case ctx_SNESSetFunction: { 
	//F_lambda
	ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	//F_lambda2
	ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(Aij,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
	//Note MAT_FLUSH_ASSEMBLY
	ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

};

struct ArcInterfaceFEMethod: public InterfaceFEMethod {

  enum interface_materials_context { ctx_IntLinearSoftening, ctx_InTBILinearSoftening, ctx_IntNone };
  interface_materials_context int_mat_ctx;

  ArcInterfaceFEMethod(
      Interface& _moab,double _YoungModulus): 
      InterfaceFEMethod(_moab,_YoungModulus),int_mat_ctx(ctx_IntLinearSoftening) {};

  double h,beta,ft,Gf,E0,g0,kappa1;
  enum interface_context { ctx_KappaUpdate = 1,  ctx_InterfaceNone = 2 };
  interface_context ctx_int;


  Vec D;
  ArcInterfaceFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,Vec& _D,
      double _YoungModulus,double _h,double _beta,double _ft,double _Gf,
      Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,interface_materials_context _int_mat_ctx = ctx_IntLinearSoftening): 
      InterfaceFEMethod(_moab,_Aij,_F,_YoungModulus,_SideSet1,_SideSet2,_SideSet3),int_mat_ctx(_int_mat_ctx),
      h(_h),beta(_beta),ft(_ft),Gf(_Gf),D(_D) {
    
    E0 = YoungModulus/h;
    g0 = ft/E0;
    kappa1 = 2*Gf/ft;
    
  };

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 

    switch(snes_ctx) {
      case ctx_SNESNone: {}
      break;
      case ctx_SNESSetFunction: { }
      break;
      case ctx_SNESSetJacobian: { }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap;
  vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap_loc;
  ublas::vector<FieldData> g;

  PetscErrorCode Calc_gap() {
    PetscFunctionBegin;
    int g_dim = g_NTRI.size()/3;
    gap.resize(g_dim);
    gap_loc.resize(g_dim);
    g.resize(g_dim);
    for(int gg = 0;gg<g_dim;gg++) {
	gap[gg] = ublas::zero_vector<FieldData>(3);
	//nodes
	double *nodeNTRI = &g_NTRI[gg*3];
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,0));
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  (gap[gg])[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	}
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  assert(dit->side_number_ptr->side_number>=3);
	  assert(dit->side_number_ptr->side_number<=5);
	  (gap[gg])[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	}
	//edges
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,2));
	for(;dit!=hi_dit;dit++) {
	  int side_number = dit->side_number_ptr->side_number;	
	  assert(side_number >= 0);
	  assert(side_number <= 2);
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double *_H1edgeN_ = &*H1edgeN[side_number].begin();
	  double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	} 
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	} 
	//faces
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
	gap_loc[gg] = prod(R,gap[gg]);
	ierr = Calc_g(gap_loc[gg],g[gg]); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode CalcDglob(const double _omega_) {
    PetscFunctionBegin;
    double E = (1-_omega_)*E0;
    ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
    Dloc(0,0) = E;
    Dglob = prod( Dloc, R );
    Dglob = prod( trans(R), Dglob );
    PetscFunctionReturn(0);
  }

  PetscErrorCode Calc_g(const ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> >& gap_loc_at_GaussPt,double& g_at_GaussPt) {
    PetscFunctionBegin;
    switch (int_mat_ctx) {
      case ctx_InTBILinearSoftening:
      case ctx_IntLinearSoftening: {
	double g2 = pow(gap_loc_at_GaussPt[0],2)+beta*(pow(gap_loc_at_GaussPt[1],2)+pow(gap_loc_at_GaussPt[2],2));
	g_at_GaussPt = sqrt( g2 );
      } break;
      
      default:
	 SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }
  
  PetscErrorCode Calc_omega(const double _kappa_,double& _omega_) {
    PetscFunctionBegin;
    switch (int_mat_ctx) {
      case ctx_IntLinearSoftening: {
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
      case ctx_InTBILinearSoftening: {
	SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      } break;
      default:
	 SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode CalcTangetDglob(const double _omega_,const double _g_,const ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> >& _gap_loc_) {
    PetscFunctionBegin;
    switch (int_mat_ctx) {
      case ctx_IntLinearSoftening: {
      double d_omega_ = 
	0.5*(2*Gf*E0+ft*ft)/((ft+(_g_-ft/E0)*E0)*Gf) - 0.5*((_g_-ft/E0)*(2*Gf*E0+ft*ft)*E0)/(pow(ft+(_g_-ft/E0)*E0,2)*Gf);
      double Et = (1-_omega_)*E0 - d_omega_*E0*_g_;
      ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
      Dloc(0,0) = Et*_gap_loc_[0]/_g_;
      Dloc(0,1) = Et*beta*_gap_loc_[1]/_g_;
      Dloc(0,2) = Et*beta*_gap_loc_[2]/_g_;
      Dglob = prod( Dloc, R );
      Dglob = prod( trans(R), Dglob );
      } break;
      case ctx_InTBILinearSoftening: {
	SETERRQ(PETSC_COMM_SELF,1,"not implemented yet");
      } break;
      default:
	 SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode RhsInt() {
    PetscFunctionBegin;
    int g_dim = g_NTRI.size()/3;
    for(int rr = 0;rr<row_mat;rr++) {
      for(int gg = 0;gg<g_dim;gg++) {
	double _kappa_ = fmax(g[gg]-g0,kappa[gg]);
	switch(ctx_int) {
	  case ctx_KappaUpdate:
	    kappa[gg] = _kappa_;
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
	      double w = area3*G_TRI_W13[gg];
	      ublas::matrix<FieldData> &N = (rowNMatrices[rr])[gg];
	      ublas::vector<FieldData> f_int = prod(trans(N),w*traction);
	      if(RowGlob[rr].size()==0) continue;
	      if(RowGlob[rr].size()!=f_int.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	    break;
	}
      }
    }
    PetscFunctionReturn(0);
  }

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
	    double w = area3*G_TRI_W13[gg];
	    ublas::matrix<FieldData> NTD = prod( trans(row_Mat), w*Dglob );
	    K(rr,cc) += prod(NTD , col_Mat ); 
	  }
	}
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
    }
    PetscFunctionReturn(0);
  }

  Tag th_kappa;
  const void* tag_data_kappa[1];
  double* kappa;

  PetscErrorCode set_ctx_int(interface_context _ctx) {
    PetscFunctionBegin;
    ctx_int = _ctx;
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

    //History
    int g_dim = g_NTRI.size()/3;
    EntityHandle fe_ent = fe_ptr->get_ent();
    vector<double> def_kappa(g_dim,0);
    rval = moab.tag_get_handle("_KAPPA",g_dim,MB_TYPE_DOUBLE,th_kappa,MB_TAG_CREAT|MB_TAG_SPARSE,&def_kappa[0]);  
    if(rval==MB_ALREADY_ALLOCATED) {
      rval = MB_SUCCESS;
    } else {
      rval = moab.tag_set_data(th_kappa,&fe_ent,1,&def_kappa[0]); 
    }
    CHKERR_PETSC(rval);
    rval = moab.tag_get_by_ptr(th_kappa,&fe_ent,1,tag_data_kappa); CHKERR_PETSC(rval);
    kappa = (double*)tag_data_kappa[0];


    switch(snes_ctx) {
      case ctx_SNESNone: {
      }
      break;
      case ctx_SNESSetJacobian: 
      case ctx_SNESSetFunction: { 
	//Apply Dirihlet BC
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    switch(snes_ctx) {
      case ctx_SNESNone: {
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetFunction: { 
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: { 
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


  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    switch(snes_ctx) {
      case ctx_SNESNone: {
      }
      break;
      case ctx_SNESSetFunction: {}
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

};

struct ArcLenghtIntElemFEMethod: public moabField::FEMethod {
  Interface& moab;
  ErrorCode rval;
  PetscErrorCode ierr;

  Mat Aij;
  Vec F,D;
  ArcLenghtCtx* arc_ptr;
  Vec GhostDiag,GhostLambdaInt;
  Range Faces3,Faces4,PostProcNodes;
  Range Edges3,Edges4;
  Range Nodes3,Nodes4;

  ArcLenghtIntElemFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,Vec& _D,
    ArcLenghtCtx *_arc_ptr): FEMethod(),moab(_moab),Aij(_Aij),F(_F),D(_D),arc_ptr(_arc_ptr) {
    PetscInt ghosts[1] = { 0 };
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostDiag);
      VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostLambdaInt);
    } else {
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostDiag);
      VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostLambdaInt);
    }
    Range prisms;
    rval = moab.get_entities_by_type(0,MBPRISM,prisms,false); CHKERR_THROW(rval);
    for(Range::iterator pit = prisms.begin();pit!=prisms.end();pit++) {
      EntityHandle f3,f4;
      rval = moab.side_element(*pit,2,3,f3); CHKERR_THROW(rval);
      rval = moab.side_element(*pit,2,4,f4); CHKERR_THROW(rval);
      Faces3.insert(f3);
      Faces4.insert(f4);
    }

    rval = moab.get_adjacencies(Faces3,1,false,Edges3); CHKERR_THROW(rval);
    rval = moab.get_adjacencies(Faces4,1,false,Edges4); CHKERR_THROW(rval);
    rval = moab.get_connectivity(Faces3,Nodes3,true); CHKERR_THROW(rval);
    rval = moab.get_connectivity(Faces4,Nodes4,true); CHKERR_THROW(rval);
    Faces3.insert(Edges3.begin(),Edges3.end());
    Faces3.insert(Nodes3.begin(),Nodes3.end());
    Faces4.insert(Edges4.begin(),Edges4.end());
    Faces4.insert(Nodes4.begin(),Nodes4.end());
 
    Range all_nodes;
    rval = moab.get_entities_by_type(0,MBVERTEX,all_nodes,true); CHKERR_THROW(rval);
    for(Range::iterator nit = all_nodes.begin();nit!=all_nodes.end();nit++) {
      double coords[3];
      rval = moab.get_coords(&*nit,1,coords);  CHKERR_THROW(rval);
      if(fabs(coords[0]-5)<1e-6) {
	PostProcNodes.insert(*nit);
      }
    }
    PetscPrintf(PETSC_COMM_WORLD,"Nb. PostProcNodes %lu\n",PostProcNodes.size());

    
  }
  ~ArcLenghtIntElemFEMethod() {
    VecDestroy(&GhostDiag);
    VecDestroy(&GhostLambdaInt);
  }

  PetscErrorCode potsProcessLoadPath() {
    PetscFunctionBegin;
    NumeredDofMoFEMEntity_multiIndex &numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator lit;
    lit = numered_dofs_rows.get<FieldName_mi_tag>().find("LAMBDA");
    if(lit == numered_dofs_rows.get<FieldName_mi_tag>().end()) PetscFunctionReturn(0);
    Range::iterator nit = PostProcNodes.begin();
    for(;nit!=PostProcNodes.end();nit++) {
      NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
      dit = numered_dofs_rows.get<MoABEnt_mi_tag>().lower_bound(*nit);
      hi_dit = numered_dofs_rows.get<MoABEnt_mi_tag>().upper_bound(*nit);
      double coords[3];
      rval = moab.get_coords(&*nit,1,coords);  CHKERR_THROW(rval);
      for(;dit!=hi_dit;dit++) {
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e -> ",lit->get_name().c_str(),lit->get_dof_rank(),lit->get_FieldData());
	PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e ",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
      	PetscPrintf(PETSC_COMM_WORLD,"-> %3.4f %3.4f %3.4f\n",coords[0],coords[1],coords[2]);
      }
    }
    PetscFunctionReturn(0);
  }


  double lambda_int;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = calulate_dx_and_dlambda(D); CHKERRQ(ierr);
	ierr = calulate_db(); CHKERRQ(ierr);
	ierr = calulate_lambda_int(lambda_int); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_lambda_int(double &_lambda_int_) {
    PetscFunctionBegin;
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(problem_ptr->get_nb_local_dofs_row());
    double *array;
    double *array_int_lambda;
    ierr = VecZeroEntries(GhostLambdaInt); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    ierr = VecGetArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    array_int_lambda[0] = 0;
    for(;dit!=hi_dit;dit++) {
      if(Nodes3.find(dit->get_ent())!=Nodes3.end()) {
	array_int_lambda[0] += array[dit->petsc_local_dof_idx];
      }
      if(Nodes4.find(dit->get_ent())!=Nodes4.end()) {
	array_int_lambda[0] -= array[dit->petsc_local_dof_idx];
      }
    }
    ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    ierr = VecRestoreArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(GhostLambdaInt,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGetArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    _lambda_int_ = arc_ptr->alpha*array_int_lambda[0] + arc_ptr->dlambda*arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
    ierr = VecRestoreArray(GhostLambdaInt,&array_int_lambda); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_db() {
    PetscFunctionBegin;
    NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
    dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
    hi_dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(
      problem_ptr->get_nb_local_dofs_row()+problem_ptr->get_nb_ghost_dofs_row());
    double *array;
    ierr = VecGetArray(arc_ptr->db,&array); CHKERRQ(ierr);
    for(;dit!=hi_dit;dit++) {
      if(Nodes3.find(dit->get_ent())!=Nodes3.end()) {
	array[dit->petsc_local_dof_idx] = +arc_ptr->alpha;
      } else if(Nodes4.find(dit->get_ent())!=Nodes4.end()) {
	  array[dit->petsc_local_dof_idx] = -arc_ptr->alpha;
      } else array[dit->petsc_local_dof_idx] = 0;
    }
    ierr = VecRestoreArray(arc_ptr->db,&array); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("LAMBDA");
    hi_dit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("LAMBDA");
    //only one LAMBDA
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

    switch(snes_ctx) {
      case ctx_SNESSetFunction: {
	arc_ptr->res_lambda = lambda_int - arc_ptr->s;
	ierr = VecSetValue(F,dit->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,ADD_VALUES); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_SELF,"\tres_lambda = %6.4e\n",arc_ptr->res_lambda);
      }
      break; 
      case ctx_SNESSetJacobian: {
	double diag = arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
	ierr = VecSetValue(GhostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
	ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1,ADD_VALUES); CHKERRQ(ierr);
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
      case ctx_SNESSetFunction: { 
	//add F_lambda
	NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
	NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
	dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
	hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
	if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = VecAXPY(F,-dit->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",dit->get_FieldData());  
	//snes_f norm
	double fnorm;
	ierr = VecNormBegin(F,NORM_2,&fnorm); CHKERRQ(ierr);	
	ierr = VecNormEnd(F,NORM_2,&fnorm);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(GhostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(GhostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(GhostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(GhostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(GhostDiag,&diag); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
      }
      break;
      default:
      break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_dx_and_dlambda(Vec x) {
    PetscFunctionBegin;
    //dx
    ierr = VecCopy(x,arc_ptr->dx); CHKERRQ(ierr);
    ierr = VecAXPY(arc_ptr->dx,-1,arc_ptr->x0); CHKERRQ(ierr);
    //dlambda
    NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	  = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
    if(dofs_moabfield_no_const.size()==0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
    dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
    if(dit==dofs_moabfield_no_const.get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
    if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    if(dit->get_petsc_local_dof_idx()!=-1) {
      double *array;
      ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
      arc_ptr->dlambda = array[dit->get_petsc_local_dof_idx()];
      array[dit->get_petsc_local_dof_idx()] = 0;
      ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    }
    int part = dit->part;
    MPI_Bcast(&(arc_ptr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
    //dx2
    ierr = VecDot(arc_ptr->dx,arc_ptr->dx,&arc_ptr->dx2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"\tdlambda = %6.4e dx2 = %6.4e\n",arc_ptr->dlambda,arc_ptr->dx2);
    PetscFunctionReturn(0);
  }

  PetscErrorCode calculate_init_dlambda(double *dlambda) {

      PetscFunctionBegin;

      *dlambda = arc_ptr->s/(arc_ptr->beta*sqrt(arc_ptr->F_lambda2));
      PetscPrintf(PETSC_COMM_WORLD,"\tInit dlambda = %6.4e s = %6.4e beta = %6.4e F_lambda2 = %6.4e\n",*dlambda,arc_ptr->s,arc_ptr->beta,arc_ptr->F_lambda2);
      double a = *dlambda;
      if(a - a != 0) {
	ostringstream sss;
	sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
      }

      PetscFunctionReturn(0);
  }

  PetscErrorCode set_dlambda_to_x(Vec x,double dlambda) {
      PetscFunctionBegin;

      NumeredDofMoFEMEntity_multiIndex& dofs_moabfield_no_const 
	    = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
      NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
      dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().lower_bound("LAMBDA");
      if(dit == dofs_moabfield_no_const.get<FieldName_mi_tag>().end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      hi_dit = dofs_moabfield_no_const.get<FieldName_mi_tag>().upper_bound("LAMBDA");
      if(distance(dit,hi_dit)!=1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

      if(dit->get_petsc_local_dof_idx()!=-1) {
	    double *array;
	    ierr = VecGetArray(x,&array); CHKERRQ(ierr);
	    double lambda_old = array[dit->get_petsc_local_dof_idx()];
	    if(!(dlambda == dlambda)) {
	      ostringstream sss;
	      sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
	      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
	    }
	    array[dit->get_petsc_local_dof_idx()] = lambda_old + dlambda;
	    PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e, %6.4e (%6.4e)\n",
	      lambda_old, array[dit->get_petsc_local_dof_idx()], dlambda);
	    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
  }

};




}
