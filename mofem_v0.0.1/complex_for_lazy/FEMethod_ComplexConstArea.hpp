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

#ifndef __MOABFEMETHOD_CONSTAREA_HPP__
#define __MOABFEMETHOD_CONSTAREA_HPP__

#include <math.h>
#include <complex>

extern "C" {
void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);
}


namespace MoFEM {

/** 
  * 
  * dN/dX = (1/A) * [ Spin[dX/dksi]*dN/deta - Spin[dX/deta]*dN/dksi ]
  *
  */
struct C_CONSTANT_AREA_FEMethod: public FieldInterface::FEMethod {

  ErrorCode rval;
  PetscErrorCode ierr;
  FieldInterface& mField;
  BaseDirihletBC *dirihlet_bc_method_ptr;

  Interface& moab;

  Mat C,Q;
  string lambda_field_name;
  int verbose;

  C_CONSTANT_AREA_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _C,Mat _Q,string _lambda_field_name,int _verbose = 0):
    mField(_mField),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),moab(_mField.get_moab()),C(_C),Q(_Q),lambda_field_name(_lambda_field_name),verbose(_verbose) { 
    //calulate face shape functions direvatives
    diffNTRI.resize(3,2);
    ShapeDiffMBTRI(&*diffNTRI.data().begin());
    //shape functions Gauss integration weigths
    G_TRI_W = G_TRI_W1;
    //nodal material positions
    dofs_X.resize(9);
    //noal values of Lagrange multipliers
    lambda.resize(3);
    //dofs indices for Lagrnage multipliers
    lambda_dofs_row_indx.resize(3);
    //lambda_dofs_row_ents.resize(3);
    lambda_dofs_col_indx.resize(3);
    //dofs indices for rows and columns
    disp_dofs_row_idx.resize(1+6+4+1);
    disp_dofs_col_idx.resize(1+6+4+1);
    disp_dofs_row_idx[0].resize(9);
    disp_dofs_col_idx[0].resize(9);
    //face node coordinates
    coords.resize(9);
  }

  //elem data
  ublas::matrix<double> diffNTRI;
  const double *G_TRI_W;
  vector<DofIdx> DirihletBC;
  vector<vector<DofIdx> > disp_dofs_col_idx,disp_dofs_row_idx;
  ublas::vector<DofIdx> lambda_dofs_row_indx,lambda_dofs_col_indx;
  ublas::vector<double,ublas::bounded_array<double,9> > coords;
  ublas::vector<double,ublas::bounded_array<double,9> > dofs_X;
  ublas::vector<double,ublas::bounded_array<double,3> > lambda;

  /**
    * \brief get face data, indices and coors and nodal values
    *
    * \param is_that_C_otherwise_dC
    */
  PetscErrorCode getData(bool is_that_C_otherwise_dC,bool trans) {
    PetscFunctionBegin;
    try {
    EntityHandle face = fe_ptr->get_ent();
    try {
      fill(lambda_dofs_row_indx.begin(),lambda_dofs_row_indx.end(),-1);
      fill(lambda_dofs_col_indx.begin(),lambda_dofs_col_indx.end(),-1);
      //fill(lambda_dofs_row_ents.begin(),lambda_dofs_row_ents.end(),no_handle);
      fill(disp_dofs_row_idx[0].begin(),disp_dofs_row_idx[0].end(),-1);
      fill(disp_dofs_col_idx[0].begin(),disp_dofs_col_idx[0].end(),-1);
      fill(lambda.begin(),lambda.end(),0);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    const EntityHandle* conn_face; 
    int num_nodes; 
    rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
    if(num_nodes != 3) SETERRQ(PETSC_COMM_SELF,1,"face should have three nodes");
    ierr = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERRQ(ierr);
    for(int nn = 0;nn<num_nodes;nn++) {
      if(is_that_C_otherwise_dC) {
	try { 
	  // it is C
	  // get rows which are Lagrabge multipliers
	  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	  dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  if(distance(dit,hi_dit)>0) {
	    if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	    if(dit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	    lambda[nn] = dit->get_FieldData();
	    lambda_dofs_row_indx[nn] = dit->get_petsc_gloabl_dof_idx();
	    //lambda_dofs_row_ents[nn] = dit->get_ent();
	  }
	} catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	}
      }
      if((!is_that_C_otherwise_dC)||(trans)) {
	try {
	  // it is dC
	  // get rows which are material nodal positions
	  FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	  dit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = data_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  if(distance(dit,hi_dit)>0) {
	    if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	    lambda[nn] = dit->get_FieldData();
	    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator diit,hi_diit;
	    diit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	    hi_diit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	    if(distance(diit,hi_diit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	    for(;diit!=hi_diit;diit++) {
	      if(diit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	      assert(nn*3+diit->get_dof_rank()<9);
	      disp_dofs_row_idx[0][nn*3+diit->get_dof_rank()] = diit->get_petsc_gloabl_dof_idx();
	    }
	  }
	} catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	}

      }
      try {
	//get columns which are material nodal positions
	FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	hi_dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	for(;dit!=hi_dit;dit++) {
	  if(dit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	  dofs_X[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	  assert(nn*3+dit->get_dof_rank()<9);
	  disp_dofs_col_idx[0][nn*3+dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
	}
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      if(trans) {
	try {
	  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	  dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  if(distance(dit,hi_dit)>0) {
	    if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	    if(dit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	    lambda_dofs_col_indx[nn] = dit->get_petsc_gloabl_dof_idx();
	  }
	} catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
	}
      }
    }
    DirihletBC.clear();
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(
      this,disp_dofs_row_idx,disp_dofs_col_idx,DirihletBC); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief calulate direvatives
   *
   */
  PetscErrorCode calcDirevatives(double *diffNTRI,double *dofs_X,double *dofs_iX,
    double *C,double *iC,double *T,double *iT) {
    PetscFunctionBegin;
    //set complex material position vector
    __CLPK_doublecomplex x_dofs_X[9];
    for(int nn = 0;nn<3;nn++) {
      for(int dd = 0;dd<3;dd++) {
        x_dofs_X[nn*3+dd].r = dofs_X[nn*3+dd];
        if(dofs_iX != NULL) {
	  x_dofs_X[nn*3+dd].i = dofs_iX[nn*3+dd];
	} else {
	  x_dofs_X[nn*3+dd].i = 0;
	}
    }}
    //calulate complex face normal vector
    __CLPK_doublecomplex x_normal[3];
    ierr = ShapeFaceNormalMBTRI_complex(diffNTRI,x_dofs_X,x_normal); CHKERRQ(ierr);
    //calulare complex normal length
    double __complex__ x_nrm2 = csqrt(
      cpow((x_normal[0].r+I*x_normal[0].i),2)+
      cpow((x_normal[1].r+I*x_normal[1].i),2)+
      cpow((x_normal[2].r+I*x_normal[2].i),2));
    double diffX_xi[3],diffX_eta[3];
    bzero(diffX_xi,3*sizeof(double));
    bzero(diffX_eta,3*sizeof(double));
    double i_diffX_xi[3],i_diffX_eta[3];
    bzero(i_diffX_xi,3*sizeof(double));
    bzero(i_diffX_eta,3*sizeof(double));
    //calulate tangent vectors
    //those vectors are in plane of face
    for(int nn = 0; nn<3; nn++) {
      diffX_xi[0] += dofs_X[3*nn + 0]*diffNTRI[2*nn+0];
      diffX_xi[1] += dofs_X[3*nn + 1]*diffNTRI[2*nn+0];
      diffX_xi[2] += dofs_X[3*nn + 2]*diffNTRI[2*nn+0];
      diffX_eta[0] += dofs_X[3*nn + 0]*diffNTRI[2*nn+1];
      diffX_eta[1] += dofs_X[3*nn + 1]*diffNTRI[2*nn+1];
      diffX_eta[2] += dofs_X[3*nn + 2]*diffNTRI[2*nn+1];
      if( dofs_iX == NULL ) continue;
      i_diffX_xi[0] += dofs_iX[3*nn + 0]*diffNTRI[2*nn+0];
      i_diffX_xi[1] += dofs_iX[3*nn + 1]*diffNTRI[2*nn+0];
      i_diffX_xi[2] += dofs_iX[3*nn + 2]*diffNTRI[2*nn+0];
      i_diffX_eta[0] += dofs_iX[3*nn + 0]*diffNTRI[2*nn+1];
      i_diffX_eta[1] += dofs_iX[3*nn + 1]*diffNTRI[2*nn+1];
      i_diffX_eta[2] += dofs_iX[3*nn + 2]*diffNTRI[2*nn+1];
    }
    double SpinX_xi[9],SpinX_eta[9];
    ierr = Spin(SpinX_xi,diffX_xi); CHKERRQ(ierr);
    ierr = Spin(SpinX_eta,diffX_eta); CHKERRQ(ierr);
    double iSpinX_xi[9],iSpinX_eta[9];
    ierr = Spin(iSpinX_xi,i_diffX_xi); CHKERRQ(ierr);
    ierr = Spin(iSpinX_eta,i_diffX_eta); CHKERRQ(ierr);
    __CLPK_doublecomplex xSpinX_xi[9],xSpinX_eta[9];
    ierr = make_complex_matrix(SpinX_xi,iSpinX_xi,xSpinX_xi); CHKERRQ(ierr);
    ierr = make_complex_matrix(SpinX_eta,iSpinX_eta,xSpinX_eta); CHKERRQ(ierr);
   __CLPK_doublecomplex xNSpinX_xi[3],xNSpinX_eta[3];
    bzero(xNSpinX_xi,3*sizeof(__CLPK_doublecomplex));
    bzero(xNSpinX_eta,3*sizeof(__CLPK_doublecomplex));
    __CLPK_doublecomplex x_zero = { 0, 0 };
    // calulate dA/dX 
    __CLPK_doublecomplex x_scalar = { -creal(1./x_nrm2), -cimag(1./x_nrm2) }; // unit [ 1/m^2 ]
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_xi,3,x_normal,1,&x_zero,xNSpinX_xi,1);
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_eta,3,x_normal,1,&x_zero,xNSpinX_eta,1);
    for(int nn = 0;nn<3;nn++) {
      double A[3],iA[3];
      A[0] = xNSpinX_xi[0].r*diffNTRI[2*nn+1]-xNSpinX_eta[0].r*diffNTRI[2*nn+0]; // unit [ 1/m ]
      A[1] = xNSpinX_xi[1].r*diffNTRI[2*nn+1]-xNSpinX_eta[1].r*diffNTRI[2*nn+0];
      A[2] = xNSpinX_xi[2].r*diffNTRI[2*nn+1]-xNSpinX_eta[2].r*diffNTRI[2*nn+0];
      iA[0] = xNSpinX_xi[0].i*diffNTRI[2*nn+1]-xNSpinX_eta[0].i*diffNTRI[2*nn+0]; // unit [ 1/m ]
      iA[1] = xNSpinX_xi[1].i*diffNTRI[2*nn+1]-xNSpinX_eta[1].i*diffNTRI[2*nn+0];
      iA[2] = xNSpinX_xi[2].i*diffNTRI[2*nn+1]-xNSpinX_eta[2].i*diffNTRI[2*nn+0];
      if(C != NULL) {
        C[3*nn + 0] = A[0]; 
        C[3*nn + 1] = A[1];
        C[3*nn + 2] = A[2];
      }
      if(iC != NULL) {
	iC[3*nn + 0] = iA[0];
	iC[3*nn + 1] = iA[1];
	iC[3*nn + 2] = iA[2];
      }
      if((T != NULL)||(iT != NULL)) {
 	double SpinA[9];
	ierr = Spin(SpinA,A); CHKERRQ(ierr); // unit [1/m]
	double iSpinA[9];
	ierr = Spin(iSpinA,iA); CHKERRQ(ierr); // unit [1/m]
	__CLPK_doublecomplex xSpinA[9];
	ierr = make_complex_matrix(SpinA,iSpinA,xSpinA); CHKERRQ(ierr);
	__CLPK_doublecomplex xT[3];
	cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinA,3,x_normal,1,&x_zero,xT,1);
	/*cerr << "NODE : " << lambda_dofs_row_ents[nn] << endl;
	cerr << "A : " << A[0] << " " << A[1] << " " << A[2] << endl;
	cerr << "N : " << x_normal[0].r << " " << x_normal[1].r << " " << x_normal[2].r << endl;
	cerr << "T : " << T[0] << " " << T[1] << " " << T[2] << endl;*/
	if(T != NULL) {
	  T[3*nn + 0] = xT[0].r;
	  T[3*nn + 1] = xT[1].r;
	  T[3*nn + 2] = xT[2].r;
	}
	if(iT != NULL) {
	  iT[3*nn + 0] = xT[0].i;
	  iT[3*nn + 1] = xT[1].i;
	  iT[3*nn + 2] = xT[2].i;
	}
      }
    }
    if( C != NULL) cblas_dscal(9,0.25, C,1);
    if(iC != NULL) cblas_dscal(9,0.25,iC,1);
    if( T != NULL) {
      cblas_dscal(9,0.25, T,1);
    }
    if(iT != NULL) {
      cblas_dscal(9,0.25,iT,1);
    }
    PetscFunctionReturn(0);
  }

  map<DofIdx,Vec> mapV;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(Q != PETSC_NULL) {
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problem_ptr,lambda_field_name,dofs)) {
	ierr = MatGetVecs(C,&mapV[dofs->get_petsc_gloabl_dof_idx()],PETSC_NULL); CHKERRQ(ierr);
	//ierr = mField.VecCreateGhost(problem_ptr->get_name(),Col,&mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
	ierr = VecZeroEntries(mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    try {
      ierr = getData(true,false); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    try {
      ublas::vector<double,ublas::bounded_array<double,9> > ELEM_CONSTRAIN(9);
      ierr = calcDirevatives(
	&*diffNTRI.data().begin(),&*dofs_X.data().begin(),NULL,
	&*ELEM_CONSTRAIN.data().begin(),NULL,NULL,NULL); CHKERRQ(ierr);
      for(int nn = 0;nn<3;nn++) {
	if(lambda_dofs_row_indx[nn]==-1) continue;
	for(int NN = 0;NN<3;NN++) {
	  if(NN!=nn) continue;
	  if(Q == PETSC_NULL) {
	  ierr = MatSetValues(C,
	    1,&(lambda_dofs_row_indx.data()[nn]),
	    3,&(disp_dofs_col_idx[0].data()[3*NN]),
	    &ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(Q != PETSC_NULL) {
	    if(mapV.find(lambda_dofs_row_indx[nn])==mapV.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = VecSetValues(mapV[lambda_dofs_row_indx[nn]],
	      3,&(disp_dofs_col_idx[0].data()[3*NN]),
	      &ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	  }
	}	
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    } 
    PetscFunctionReturn(0);
  }
  
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    if(Q != PETSC_NULL) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
      Vec Qv;
      ierr = mField.VecCreateGhost(problem_ptr->get_name(),Col,&Qv); CHKERRQ(ierr);
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problem_ptr,lambda_field_name,dofs)) {
	if(mapV.find(dofs->get_petsc_gloabl_dof_idx())==mapV.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = VecAssemblyBegin(mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
	ierr = MatMult(Q,mapV[dofs->get_petsc_gloabl_dof_idx()],Qv); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(Qv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(Qv,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *array;
	ierr = VecGetArray(Qv,&array); CHKERRQ(ierr);
	if(dofs->get_part()==pcomm->rank()) {
	  vector<DofIdx> glob_idx;
	  vector<double> vals;
	  for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(problem_ptr,dofs->get_ent(),diit)) {
	    if(diit->get_name() != "MESH_NODE_POSITIONS") continue;
	    glob_idx.push_back(diit->get_petsc_gloabl_dof_idx());
	    vals.push_back(array[diit->get_petsc_local_dof_idx()]);
	  }
	  int row = dofs->get_petsc_gloabl_dof_idx();
	  ierr = MatSetValues(C,
		1,&row,glob_idx.size(),&(glob_idx[0]),
		&vals[0],INSERT_VALUES); CHKERRQ(ierr);
	}
	ierr = VecRestoreArray(Qv,&array); CHKERRQ(ierr);
	ierr = VecDestroy(&mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
      }
      ierr = VecDestroy(&Qv); CHKERRQ(ierr);
      mapV.clear();
    }
    PetscFunctionReturn(0);
  }
  
};

/** 
  * 
  * Calulate direcative of dA^2/dX^2 
  *
  */
struct dCTgc_CONSTANT_AREA_FEMethod: public C_CONSTANT_AREA_FEMethod {

  Mat dCT;
  double gc;  
  const double eps;

  dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _dCT,string _lambda_field_name,int _verbose = 0):
    C_CONSTANT_AREA_FEMethod(_mField,_dirihlet_bc_method_ptr,PETSC_NULL,PETSC_NULL,_lambda_field_name,_verbose),dCT(_dCT),eps(1e-10) {}

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = C_CONSTANT_AREA_FEMethod::preProcess(); CHKERRQ(ierr);
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    EntityHandle face = fe_ptr->get_ent();
    try {
	ierr = getData(false,false); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    try {
      double center[3]; 
      const EntityHandle* conn_face; 
      int num_nodes; 
      rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
      rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
      tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
      cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
      double r = cblas_dnrm2(3,center,1);
      for(int NN = 0;NN<3;NN++) {
	for(int dd = 0;dd<3;dd++) {
	  ublas::vector<double,ublas::bounded_array<double,9> > idofs_X(9,0);
	  idofs_X[NN*3+dd] = r*eps;
	  ublas::vector<double,ublas::bounded_array<double,9> > dELEM_CONSTRAIN(9);
	  ierr = calcDirevatives(&*diffNTRI.data().begin(),
	    &*dofs_X.data().begin(),
	    &*idofs_X.data().begin(),
	    NULL,&*dELEM_CONSTRAIN.data().begin(),NULL,NULL); CHKERRQ(ierr);
	  dELEM_CONSTRAIN /= r*eps;
	  /*cerr << idofs_X << endl;
	  cerr << dELEM_CONSTRAIN << endl;
	  cerr << lambda_dofs_row_indx << endl;*/
	  dELEM_CONSTRAIN *= gc;
	  //cerr << dELEM_CONSTRAIN << endl;
	  ierr = MatSetValues(dCT,
	    9,&(disp_dofs_row_idx[0].data()[0]),
	    1,&(disp_dofs_col_idx[0].data()[3*NN+dd]),
	    &dELEM_CONSTRAIN.data()[0],ADD_VALUES); CHKERRQ(ierr);
	}
      }
    } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    } 
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  
};

struct Snes_CTgc_CONSTANT_AREA_FEMethod: public FieldInterface::FEMethod {

  FieldInterface& mField;
  BaseDirihletBC *dirihlet_bc_method_ptr;
  matPROJ_ctx &proj_ctx;
  string problem;
  string lambda_field_name;
  double gc;
  int verbose;

  Snes_CTgc_CONSTANT_AREA_FEMethod(
    FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,matPROJ_ctx &_proj_all_ctx,
    string _problem,string _lambda_field_name,int _verbose = 0):
    mField(_mField),dirihlet_bc_method_ptr(_dirihlet_bc_method_ptr),proj_ctx(_proj_all_ctx),
    problem(_problem),lambda_field_name(_lambda_field_name),verbose(_verbose) {}

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is fracture energy ?)");
    }
    
    Vec D;
    ierr = mField.VecCreateGhost(problem,Col,&D); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost(problem,Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    Vec _D_;
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Col,&_D_); CHKERRQ(ierr);
    VecScatter scatter;
    string y_problem = "C_CRACKFRONT_MATRIX";
    ierr = mField.VecScatterCreate(D,problem,Col,_D_,y_problem,Col,&scatter); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter,D,_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter,D,_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("C_CRACKFRONT_MATRIX",Col,_D_,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecDestroy(&_D_); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);

    C_CONSTANT_AREA_FEMethod C_AREA_ELEM(mField,dirihlet_bc_method_ptr,proj_ctx.C,PETSC_NULL,lambda_field_name);

    ierr = MatSetOption(proj_ctx.C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
    ierr = MatSetOption(proj_ctx.C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

    ierr = MatZeroEntries(proj_ctx.C); CHKERRQ(ierr);
    ierr = mField.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_AREA_ELEM);  CHKERRQ(ierr);
    ierr = MatAssemblyBegin(proj_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(proj_ctx.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    /*{
      //Matrix View
      MatView(C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      std::string wait;
      std::cin >> wait;
    }*/

    Vec LambdaVec;
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Row,&LambdaVec); CHKERRQ(ierr);
    const MoFEMProblem *front_problem_ptr;
    ierr = mField.get_problem("C_CRACKFRONT_MATRIX",&front_problem_ptr); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(front_problem_ptr,"LAMBDA_CRACKFRONT_AREA",dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      ierr = VecSetValue(LambdaVec,dit->get_petsc_gloabl_dof_idx(),gc,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(LambdaVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(LambdaVec); CHKERRQ(ierr);
    //ierr = VecView(LambdaVec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    Vec _f_;
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",Col,&_f_); CHKERRQ(ierr);
    ierr = MatMultTranspose(proj_ctx.C,LambdaVec,_f_); CHKERRQ(ierr);
    //PetscReal _f_nrm2;
    //ierr = VecNorm(_f_, NORM_2,&_f_nrm2); CHKERRQ(ierr);
    //PetscPrintf(PETSC_COMM_WORLD,"\tfront f_nrm2 = %6.4e\n",_f_nrm2);

    ierr = mField.VecScatterCreate(snes_f,problem,Row,_f_,y_problem,Col,&scatter); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter,_f_,snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter,_f_,snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);
    ierr = VecDestroy(&_f_); CHKERRQ(ierr);

    ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

};

struct Snes_dCTgc_CONSTANT_AREA_FEMethod: public dCTgc_CONSTANT_AREA_FEMethod {

  matPROJ_ctx *proj_ctx;
  Mat K;

  Snes_dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,
    matPROJ_ctx &_proj_all_ctx,string _lambda_field_name,int _verbose = 0):
    dCTgc_CONSTANT_AREA_FEMethod(_mField,_dirihlet_bc_method_ptr,_proj_all_ctx.K,_lambda_field_name,_verbose),
    proj_ctx(&_proj_all_ctx),K(_proj_all_ctx.K) {}

  Snes_dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,Mat _K,string _lambda_field_name,int _verbose = 0):
    dCTgc_CONSTANT_AREA_FEMethod(_mField,_dirihlet_bc_method_ptr,_K,_lambda_field_name,_verbose),proj_ctx(NULL),K(_K) {}


  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = dCTgc_CONSTANT_AREA_FEMethod::preProcess(); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    ierr = dCTgc_CONSTANT_AREA_FEMethod::postProcess(); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(K,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct TangentFrontConstrain_FEMethod: public C_CONSTANT_AREA_FEMethod {

  const double eps;
  TangentFrontConstrain_FEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_bc_method_ptr,string _lambda_field_name,int _verbose = 0):
    C_CONSTANT_AREA_FEMethod(_mField,_dirihlet_bc_method_ptr,PETSC_NULL,PETSC_NULL,_lambda_field_name,_verbose),eps(1e-10) {}

  Tag th_front_tangent;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ierr = C_CONSTANT_AREA_FEMethod::preProcess(); CHKERRQ(ierr);
    //TAG  - only for one proc analysis
    /*double def[] = {0,0,0};
    rval = mField.get_moab().tag_get_handle("FRONT_TANGENT",3,MB_TYPE_DOUBLE,
      th_front_tangent,MB_TAG_CREAT|MB_TAG_SPARSE,&def); CHKERR_THROW(rval);*/
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	/*//resent tags - only for one proc analysis
	for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,lambda_field_name,dit)) {
	  EntityHandle ent = dit->get_ent();
	  rval = mField.get_moab().tag_set_data(th_front_tangent,&ent,1,def); CHKERR_PETSC(rval);
	}*/
      } break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      } break;
      default: {
      } break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    EntityHandle face = fe_ptr->get_ent();
    Range tet;
    rval = mField.get_moab().get_adjacencies(&face,1,3,false,tet); CHKERR_PETSC(rval);


    try {
      ierr = getData(true,true); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    try {
      //tag create
      Tag th_interface_side;
      rval = moab.tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
      int side;
      rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
      ublas::vector<double,ublas::bounded_array<double,9> > ELEM_CONSTRAIN1(9);
      ierr = calcDirevatives(
	&*diffNTRI.data().begin(),&*dofs_X.data().begin(),
	NULL,NULL,NULL,&*ELEM_CONSTRAIN1.data().begin(),NULL); CHKERRQ(ierr);
      //take in account face orientation in respect crack surface
      if(side == 1) {
	ELEM_CONSTRAIN1 *= -1;
      }
      //tangent
      if(snes_ctx == ctx_SNESSetJacobian) {
	double center[3]; 
	tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
	cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
	double r = cblas_dnrm2(3,center,1);
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    // ---> calulate tangent starts here
	    ublas::vector<double,ublas::bounded_array<double,9> > idofs_X(9,0);
	    idofs_X[nn*3+dd] = r*eps;
	    ublas::vector<double,ublas::bounded_array<double,9> > dELEM_CONSTRAIN1(9);
	    ierr = calcDirevatives(
	      &*diffNTRI.data().begin(),
	      &*dofs_X.data().begin(),
	      &*idofs_X.data().begin(),
	      NULL,NULL,NULL,&*dELEM_CONSTRAIN1.data().begin()); CHKERRQ(ierr);
	    if(side == 1) {
	      dELEM_CONSTRAIN1 /= -r*eps;
	    } else {
	      dELEM_CONSTRAIN1 /= +r*eps;
	    }
	    //dC
	    double g = 0;
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      for(int ddd = 0;ddd<3;ddd++) {
		g += dELEM_CONSTRAIN1[nnn*3+ddd]*(dofs_X[3*nnn+ddd]-coords[3*nnn+ddd]);
	      }
	    }
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      ierr = MatSetValues(*snes_B,
		1,&lambda_dofs_row_indx[nnn],1,&disp_dofs_col_idx[0][3*nn+dd],
		&g,ADD_VALUES); CHKERRQ(ierr);
	    }
	    //dCT_lambda
	    for(int nnn = 0;nnn<3;nnn++) {
	      for(int ddd = 0;ddd<3;ddd++) {
		dELEM_CONSTRAIN1[nnn*3+ddd] *= lambda[nnn];
	      }
	    }
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      ierr = MatSetValues(*snes_B,
		3,&disp_dofs_row_idx[0][3*nnn],
		1,&disp_dofs_col_idx[0][3*nn+dd],
		&dELEM_CONSTRAIN1[3*nnn],ADD_VALUES); CHKERRQ(ierr);
	    }
	    // ---> calulate tangent end here
	  }
	}
      }
      switch(snes_ctx) {
	case ctx_SNESSetFunction: { 
	  ublas::vector<double,ublas::bounded_array<double,3> > g(3);
	  for(int nn = 0;nn<3;nn++) {
	    g[nn] = cblas_ddot(3,&ELEM_CONSTRAIN1[3*nn],1,&dofs_X[3*nn],1);
	    g[nn] -= cblas_ddot(3,&ELEM_CONSTRAIN1[3*nn],1,&coords[3*nn],1);
	  }
	  //cerr << "g : " << g << endl;
	  ierr = VecSetValues(snes_f,3,&*lambda_dofs_row_indx.data().begin(),&*g.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  ublas::vector<double,ublas::bounded_array<double,9> > f(9);
	  for(int nn = 0;nn<3;nn++) {
	    for(int dd = 0;dd<3;dd++) {
	      f[nn*3+dd] = lambda[nn]*ELEM_CONSTRAIN1[3*nn+dd];
	    }
	  }
	  //cerr << "f : " << f << endl;
	  ierr = VecSetValues(snes_f,9,&disp_dofs_row_idx[0][0],&*f.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  /*//TAG - only for one proc analysis
	  for(int nn = 0;nn<3;nn++) { 
	    EntityHandle ent = lambda_dofs_row_ents[nn];
	    if(ent == no_handle) continue;
	    double *t;
	    rval = mField.get_moab().tag_get_by_ptr(th_front_tangent,&ent,1,(const void **)&t); CHKERR_PETSC(rval);
	    cblas_daxpy(3,+1,&ELEM_CONSTRAIN1[3*nn],1,t,1);
	  }*/
	} break;
	case ctx_SNESSetJacobian: {
	  for(int nn = 0;nn<3;nn++) {
	    int lambda_dof_idx = lambda_dofs_row_indx[nn];
	    ierr = MatSetValues(*snes_B,1,&lambda_dof_idx,3,&disp_dofs_col_idx[0][3*nn],&ELEM_CONSTRAIN1[3*nn],ADD_VALUES); CHKERRQ(ierr);
	    lambda_dof_idx = lambda_dofs_col_indx[nn];
	    ierr = MatSetValues(*snes_B,3,&disp_dofs_row_idx[0][3*nn],1,&lambda_dof_idx,&ELEM_CONSTRAIN1[3*nn],ADD_VALUES); CHKERRQ(ierr);
	  }
	} break;
	default:
	  break;
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    } 
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      } break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      } break;
      default: {
      } break;
    }
    PetscFunctionReturn(0);
  }
  
};

}

#endif // __MOABFEMETHOD_CONSTAREA_HPP__
