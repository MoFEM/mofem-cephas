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

#include "/usr/include/complex.h"

void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);
}

using namespace ObosleteUsersModules;

namespace MoFEM {

/** 
  * 
  * dN/dX = (1/A) * [ Spin[dX/dksi]*dN/deta - Spin[dX/deta]*dN/dksi ]
  *
  */
struct C_CONSTANT_AREA_FEMethod: public FEMethod {

  ErrorCode rval;
  PetscErrorCode ierr;
  FieldInterface& mField;

  Interface& moab;

  Mat C,Q;
  string lambda_field_name;
  int verbose;

  C_CONSTANT_AREA_FEMethod(FieldInterface& _mField,Mat _C,Mat _Q,string _lambda_field_name,int _verbose = 0):
    mField(_mField),moab(_mField.get_moab()),C(_C),Q(_Q),lambda_field_name(_lambda_field_name),verbose(_verbose) { 
    //calculate face shape functions direvatives
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
    disp_dofs_row_idx.resize(9);
    disp_dofs_col_idx.resize(9);
    local_disp_dofs_row_idx.resize(9);
    //face node coordinates
    coords.resize(9);
  }

  //elem data
  ublas::matrix<double> diffNTRI;
  const double *G_TRI_W;
  vector<DofIdx> DirichletBC;
  ublas::vector<DofIdx> disp_dofs_col_idx,disp_dofs_row_idx;
  ublas::vector<DofIdx> local_disp_dofs_row_idx;
  ublas::vector<DofIdx> lambda_dofs_row_indx,lambda_dofs_col_indx;
  ublas::vector<double,ublas::bounded_array<double,9> > coords;
  ublas::vector<double,ublas::bounded_array<double,9> > dofs_X;
  ublas::vector<double,ublas::bounded_array<double,3> > lambda;
  //vector<EntityHandle> lambda_dofs_row_ents;
  EntityHandle face;

  /**
    * \brief get face data, indices and coors and nodal values
    *
    * \param is_that_C_otherwise_dC
    */
  PetscErrorCode getData(bool is_that_C_otherwise_dC,bool trans) {
    PetscFunctionBegin;
    try {
    face = fePtr->get_ent();
    try {
      fill(lambda_dofs_row_indx.begin(),lambda_dofs_row_indx.end(),-1);
      fill(lambda_dofs_col_indx.begin(),lambda_dofs_col_indx.end(),-1);
      fill(disp_dofs_row_idx.begin(),disp_dofs_row_idx.end(),-1);
      fill(disp_dofs_col_idx.begin(),disp_dofs_col_idx.end(),-1);
      fill(local_disp_dofs_row_idx.begin(),local_disp_dofs_row_idx.end(),-1);
      fill(lambda.begin(),lambda.end(),0);
      //fill(lambda_dofs_row_ents.begin(),lambda_dofs_row_ents.end(),no_handle);
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
	  dit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
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
	  dit = dataPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = dataPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  if(distance(dit,hi_dit)>0) {
	    if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	    lambda[nn] = dit->get_FieldData();
	    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator diit,hi_diit;
	    diit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	    hi_diit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	    if(distance(diit,hi_diit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	    for(;diit!=hi_diit;diit++) {
	      if(diit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	      assert(nn*3+diit->get_dof_rank()<9);
	      disp_dofs_row_idx[nn*3+diit->get_dof_rank()] = diit->get_petsc_gloabl_dof_idx();
	      local_disp_dofs_row_idx[nn*3+diit->get_dof_rank()] = diit->get_petsc_local_dof_idx();
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
	dit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	hi_dit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
	if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
	for(;dit!=hi_dit;dit++) {
	  if(dit->get_petsc_local_dof_idx()<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
	  dofs_X[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
	  assert(nn*3+dit->get_dof_rank()<9);
	  disp_dofs_col_idx[nn*3+dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
	}
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      if(trans) {
	try {
	  FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	  dit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
	  hi_dit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
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
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  /**
   * \brief calculate direvatives
   *
   */
  PetscErrorCode calcDirevatives(double *diffNTRI,double *dofs_X,double *dofs_iX,
    double *C,double *iC,double *T,double *iT) {
    PetscFunctionBegin;
    double diffX_xi[3],diffX_eta[3];
    bzero(diffX_xi,3*sizeof(double));
    bzero(diffX_eta,3*sizeof(double));
    double i_diffX_xi[3],i_diffX_eta[3];
    bzero(i_diffX_xi,3*sizeof(double));
    bzero(i_diffX_eta,3*sizeof(double));
    Range adj_side_elems;
    BitRefLevel bit = problemPtr->get_BitRefLevel();
    ierr = mField.get_adjacencies(bit,&face,1,3,adj_side_elems); CHKERRQ(ierr);
    adj_side_elems = adj_side_elems.subset_by_type(MBTET);
    if(adj_side_elems.size()==0) {
      Range adj_tets_on_surface;
      BitRefLevel bit_tet_on_surface;
      bit_tet_on_surface.set(BITREFLEVEL_SIZE-2);
      ierr = mField.get_adjacencies(bit_tet_on_surface,&face,1,3,adj_tets_on_surface,Interface::INTERSECT,0); CHKERRQ(ierr);
      adj_side_elems.insert(*adj_tets_on_surface.begin());
    }
    if(adj_side_elems.size()!=1) {
      adj_side_elems.clear();
      ierr = mField.get_adjacencies(bit,&face,1,3,adj_side_elems,Interface::INTERSECT,5); CHKERRQ(ierr);
      Range::iterator it = adj_side_elems.begin();
      for(;it!=adj_side_elems.end();it++) {	
	Range nodes;
	rval = mField.get_moab().get_connectivity(&*it,1,nodes,true); CHKERR_PETSC(rval);
	PetscPrintf(PETSC_COMM_WORLD,"%lu %lu %lu %lu\n",nodes[0],nodes[1],nodes[2],nodes[3]);
      }
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {
	EntityHandle out_meshset;
	rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,adj_side_elems); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,&face,1); CHKERR_PETSC(rval);
	rval = mField.get_moab().write_file("debug_error.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      }
      SETERRQ1(PETSC_COMM_SELF,1,"expect 1 tet but is %u",adj_side_elems.size());
    }
    EntityHandle side_elem = *adj_side_elems.begin();
    int order[] = {0, 1, 2};
    if(side_elem != 0) {
      int side_number,sense,offset;
      rval = moab.side_number(side_elem,face,side_number,sense,offset); CHKERR_PETSC(rval);
      if(sense == -1) {
	order[0] = 1;
	order[1] = 0;
      }
    }
    //calculate tangent vectors
    //those vectors are in plane of face
    for(int nn = 0; nn<3; nn++) {
      diffX_xi[0] += dofs_X[3*order[nn] + 0]*diffNTRI[2*nn+0]; // unit [ m ]
      diffX_xi[1] += dofs_X[3*order[nn] + 1]*diffNTRI[2*nn+0];
      diffX_xi[2] += dofs_X[3*order[nn] + 2]*diffNTRI[2*nn+0];
      diffX_eta[0] += dofs_X[3*order[nn] + 0]*diffNTRI[2*nn+1];
      diffX_eta[1] += dofs_X[3*order[nn] + 1]*diffNTRI[2*nn+1];
      diffX_eta[2] += dofs_X[3*order[nn] + 2]*diffNTRI[2*nn+1];
      if( dofs_iX == NULL ) continue;
      i_diffX_xi[0] += dofs_iX[3*order[nn] + 0]*diffNTRI[2*nn+0];
      i_diffX_xi[1] += dofs_iX[3*order[nn] + 1]*diffNTRI[2*nn+0];
      i_diffX_xi[2] += dofs_iX[3*order[nn] + 2]*diffNTRI[2*nn+0];
      i_diffX_eta[0] += dofs_iX[3*order[nn] + 0]*diffNTRI[2*nn+1];
      i_diffX_eta[1] += dofs_iX[3*order[nn] + 1]*diffNTRI[2*nn+1];
      i_diffX_eta[2] += dofs_iX[3*order[nn] + 2]*diffNTRI[2*nn+1];
    }
    //spins
    double SpinX_xi[9],SpinX_eta[9];
    ierr = Spin(SpinX_xi,diffX_xi); CHKERRQ(ierr);
    ierr = Spin(SpinX_eta,diffX_eta); CHKERRQ(ierr);
    double iSpinX_xi[9],iSpinX_eta[9];
    ierr = Spin(iSpinX_xi,i_diffX_xi); CHKERRQ(ierr);
    ierr = Spin(iSpinX_eta,i_diffX_eta); CHKERRQ(ierr);
    __CLPK_doublecomplex xSpinX_xi[9],xSpinX_eta[9];
    ierr = make_complex_matrix(SpinX_xi,iSpinX_xi,xSpinX_xi); CHKERRQ(ierr); // unit [ m ]
    ierr = make_complex_matrix(SpinX_eta,iSpinX_eta,xSpinX_eta); CHKERRQ(ierr); // unit [ m ]
    //calculate complex face normal vector
    __CLPK_doublecomplex x_one = { 1, 0 };
    __CLPK_doublecomplex x_zero = { 0, 0 };
    __CLPK_doublecomplex x_normal[3];
    __CLPK_doublecomplex x_diffX_eta[3]; 
    x_diffX_eta[0].r = diffX_eta[0]; x_diffX_eta[0].i = i_diffX_eta[0];
    x_diffX_eta[1].r = diffX_eta[1]; x_diffX_eta[1].i = i_diffX_eta[1];
    x_diffX_eta[2].r = diffX_eta[2]; x_diffX_eta[2].i = i_diffX_eta[2];
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_one,xSpinX_xi,3,x_diffX_eta,1,&x_zero,x_normal,1); // unit [ m^2 ]
    /*//debug
    Tag th_normal1,th_normal2;
    double def_NORMAL[] = {0, 0, 0};
    rval = moab.tag_get_handle("NORMAL1",3,MB_TYPE_DOUBLE,th_normal1,MB_TAG_CREAT|MB_TAG_SPARSE,def_NORMAL); CHKERR_PETSC(rval);
    rval = moab.tag_get_handle("NORMAL2",3,MB_TYPE_DOUBLE,th_normal2,MB_TAG_CREAT|MB_TAG_SPARSE,def_NORMAL); CHKERR_PETSC(rval);
    double normal1[3] = { x_normal[0].r, x_normal[1].r, x_normal[2].r };
    double normal2[3];
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&*coords.data().begin(),normal2); CHKERRQ(ierr);
    if(side_elem!=0) {
      int side_number,sense,offset;
      rval = moab.side_number(side_elem,face,side_number,sense,offset); CHKERR_PETSC(rval);
      if(sense == -1) {
	cblas_dscal(3,-1,normal2,1);
      }
    }
    double normal1_nrm2 = cblas_dnrm2(3,normal1,1);
    cblas_dscal(3,1./normal1_nrm2,normal1,1);
    double normal2_nrm2 = cblas_dnrm2(3,normal2,1);
    cblas_dscal(3,1./normal2_nrm2,normal2,1);
    rval = moab.tag_set_data(th_normal1,&face,1,normal1); CHKERR_PETSC(rval);
    rval = moab.tag_set_data(th_normal2,&face,1,normal2); CHKERR_PETSC(rval);*/
    //calulare complex normal length
    double __complex__ x_nrm2 = csqrt(
      cpow((x_normal[0].r+I*x_normal[0].i),2+I*0)+
      cpow((x_normal[1].r+I*x_normal[1].i),2+I*0)+
      cpow((x_normal[2].r+I*x_normal[2].i),2+I*0));
    // calculate dA/dX 
    __CLPK_doublecomplex xNSpinX_xi[3],xNSpinX_eta[3];
    bzero(xNSpinX_xi,3*sizeof(__CLPK_doublecomplex));
    bzero(xNSpinX_eta,3*sizeof(__CLPK_doublecomplex));
    __CLPK_doublecomplex x_scalar = { -creal(1./x_nrm2), -cimag(1./x_nrm2) }; // unit [ 1/m^2 ]
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_xi,3,x_normal,1,&x_zero,xNSpinX_xi,1); // unit [ 1/m^2 * m = 1/m ]
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_eta,3,x_normal,1,&x_zero,xNSpinX_eta,1); // unit [ 1/m^2 * m = 1/m ]
    if( C!=NULL) bzero( C,9*sizeof(double));
    if(iC!=NULL) bzero(iC,9*sizeof(double));
    if( T!=NULL) bzero( T,9*sizeof(double));
    if(iT!=NULL) bzero(iT,9*sizeof(double));
    for(int nn = 0;nn<3;nn++) {
      double A[3],iA[3];
      A[0] = xNSpinX_xi[0].r*diffNTRI[2*order[nn]+1]-xNSpinX_eta[0].r*diffNTRI[2*order[nn]+0]; // unit [ 1/m ]
      A[1] = xNSpinX_xi[1].r*diffNTRI[2*order[nn]+1]-xNSpinX_eta[1].r*diffNTRI[2*order[nn]+0];
      A[2] = xNSpinX_xi[2].r*diffNTRI[2*order[nn]+1]-xNSpinX_eta[2].r*diffNTRI[2*order[nn]+0];
      iA[0] = xNSpinX_xi[0].i*diffNTRI[2*order[nn]+1]-xNSpinX_eta[0].i*diffNTRI[2*order[nn]+0]; // unit [ 1/m ]
      iA[1] = xNSpinX_xi[1].i*diffNTRI[2*order[nn]+1]-xNSpinX_eta[1].i*diffNTRI[2*order[nn]+0];
      iA[2] = xNSpinX_xi[2].i*diffNTRI[2*order[nn]+1]-xNSpinX_eta[2].i*diffNTRI[2*order[nn]+0];
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
	// make spin matrix to calculate cross product
	ierr = make_complex_matrix(SpinA,iSpinA,xSpinA); CHKERRQ(ierr); 
	__CLPK_doublecomplex xT[3];
	cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinA,3,x_normal,1,&x_zero,xT,1); // unit [1/m]
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
    if( T != NULL) cblas_dscal(9,0.25, T,1);
    if(iT != NULL) cblas_dscal(9,0.25,iT,1);
    PetscFunctionReturn(0);
  }

  map<DofIdx,Vec> mapV;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(Q != PETSC_NULL) {
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problemPtr,lambda_field_name,dofs)) {
	ierr = MatGetVecs(C,&mapV[dofs->get_petsc_gloabl_dof_idx()],PETSC_NULL); CHKERRQ(ierr);
	//ierr = mField.VecCreateGhost(problemPtr->get_name(),COL,&mapV[dofs->get_petsc_gloabl_dof_idx()]); CHKERRQ(ierr);
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
	    3,&(disp_dofs_col_idx.data()[3*NN]),
	    &ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(Q != PETSC_NULL) {
	    if(mapV.find(lambda_dofs_row_indx[nn])==mapV.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = VecSetValues(mapV[lambda_dofs_row_indx[nn]],
	      3,&(disp_dofs_col_idx.data()[3*NN]),
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
      ierr = mField.VecCreateGhost(problemPtr->get_name(),COL,&Qv); CHKERRQ(ierr);
      for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problemPtr,lambda_field_name,dofs)) {
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
	  for(_IT_NUMEREDDOFMOFEMENTITY_COL_BY_ENT_FOR_LOOP_(problemPtr,dofs->get_ent(),diit)) {
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

struct C_FRONT_TANGENT_FEMethod: public C_CONSTANT_AREA_FEMethod {

  C_FRONT_TANGENT_FEMethod(FieldInterface& _mField,Mat _C,Mat _Q,string _lambda_field_name,int _verbose = 0):
    C_CONSTANT_AREA_FEMethod(_mField,_C,_Q,_lambda_field_name,_verbose) {}

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
	NULL,NULL,&*ELEM_CONSTRAIN.data().begin(),NULL); CHKERRQ(ierr);
      //take in account face orientation in respect crack surface
      Tag th_interface_side;
      rval = moab.tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
      int side;
      rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
      if(side == 1) {
	ELEM_CONSTRAIN *= -1;
      }
      for(int nn = 0;nn<3;nn++) {
	if(lambda_dofs_row_indx[nn]==-1) continue;
	for(int NN = 0;NN<3;NN++) {
	  if(NN!=nn) continue;
	  if(Q == PETSC_NULL) {
	  ierr = MatSetValues(C,
	    1,&(lambda_dofs_row_indx.data()[nn]),
	    3,&(disp_dofs_col_idx.data()[3*NN]),
	    &ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(Q != PETSC_NULL) {
	    if(mapV.find(lambda_dofs_row_indx[nn])==mapV.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = VecSetValues(mapV[lambda_dofs_row_indx[nn]],
	      3,&(disp_dofs_col_idx.data()[3*NN]),
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

};

/** 
  * 
  * Calculate direcative of dA^2/dX^2 
  *
  */
struct dCTgc_CONSTANT_AREA_FEMethod: public C_CONSTANT_AREA_FEMethod {

  Mat dCT;
  double gc;  
  const double eps;

  dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,Mat _dCT,string _lambda_field_name,int _verbose = 0):
    C_CONSTANT_AREA_FEMethod(_mField,PETSC_NULL,PETSC_NULL,_lambda_field_name,_verbose),dCT(_dCT),eps(1e-10) {}

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
    EntityHandle face = fePtr->get_ent();
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
	    9,&(disp_dofs_row_idx.data()[0]),
	    1,&(disp_dofs_col_idx.data()[3*NN+dd]),
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

struct Snes_CTgc_CONSTANT_AREA_FEMethod: public FEMethod {

  FieldInterface& mField;
  matPROJ_ctx &proj_ctx;
  string problem;
  string lambda_field_name;
  double gc;
  int verbose;

  Snes_CTgc_CONSTANT_AREA_FEMethod(
    FieldInterface& _mField,matPROJ_ctx &_proj_all_ctx,
    string _problem,string _lambda_field_name,int _verbose = 0):
    mField(_mField),proj_ctx(_proj_all_ctx),
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
    ierr = mField.VecCreateGhost(problem,COL,&D); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost(problem,COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    Vec _D_;
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",COL,&_D_); CHKERRQ(ierr);
    VecScatter scatter;
    string y_problem = "C_CRACKFRONT_MATRIX";
    ierr = mField.VecScatterCreate(D,problem,COL,_D_,y_problem,COL,&scatter); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter,D,_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter,D,_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(_D_,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = mField.set_local_VecCreateGhost("C_CRACKFRONT_MATRIX",COL,_D_,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecDestroy(&_D_); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);

    C_CONSTANT_AREA_FEMethod C_AREA_ELEM(mField,proj_ctx.C,PETSC_NULL,lambda_field_name);

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
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,&LambdaVec); CHKERRQ(ierr);
    const MoFEMProblem *front_problemPtr;
    ierr = mField.get_problem("C_CRACKFRONT_MATRIX",&front_problemPtr); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(front_problemPtr,"LAMBDA_CRACKFRONT_AREA",dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      ierr = VecSetValue(LambdaVec,dit->get_petsc_gloabl_dof_idx(),gc,INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(LambdaVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(LambdaVec); CHKERRQ(ierr);
    //ierr = VecView(LambdaVec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    Vec _f_;
    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",COL,&_f_); CHKERRQ(ierr);
    ierr = MatMultTranspose(proj_ctx.C,LambdaVec,_f_); CHKERRQ(ierr);
    //PetscReal _f_nrm2;
    //ierr = VecNorm(_f_, NORM_2,&_f_nrm2); CHKERRQ(ierr);
    //PetscPrintf(PETSC_COMM_WORLD,"\tfront f_nrm2 = %6.4e\n",_f_nrm2);

    ierr = mField.VecScatterCreate(snes_f,problem,ROW,_f_,y_problem,COL,&scatter); CHKERRQ(ierr);
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

  Snes_dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,
    matPROJ_ctx &_proj_all_ctx,string _lambda_field_name,int _verbose = 0):
    dCTgc_CONSTANT_AREA_FEMethod(_mField,_proj_all_ctx.K,_lambda_field_name,_verbose),
    proj_ctx(&_proj_all_ctx),K(_proj_all_ctx.K) {}

  Snes_dCTgc_CONSTANT_AREA_FEMethod(FieldInterface& _mField,Mat _K,string _lambda_field_name,int _verbose = 0):
    dCTgc_CONSTANT_AREA_FEMethod(_mField,_K,_lambda_field_name,_verbose),proj_ctx(NULL),K(_K) {}


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

}

#endif // __MOABFEMETHOD_CONSTAREA_HPP__
