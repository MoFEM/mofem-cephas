
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

//#include "FieldInterface.hpp"
//#include "FieldCore.hpp"
//#include "complex_for_lazy.h"
#include <math.h>
//#include <complex>

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

  Interface& moab;

  Mat C,Q;
  Range surface;
  string lambda_field_name;
  int verbose;

  C_CONSTANT_AREA_FEMethod(FieldInterface& _mField,Range &_surface,Mat _C,Mat _Q,string _lambda_field_name,int _verbose = 0):
    mField(_mField),moab(_mField.get_moab()),C(_C),Q(_Q),surface(_surface),lambda_field_name(_lambda_field_name),verbose(_verbose) { 
    diffNTRI.resize(3,2);
    ShapeDiffMBTRI(&*diffNTRI.data().begin());
    g_NTRI3.resize(1,3);
    ShapeMBTRI(&*g_NTRI3.data().begin(),G_TRI_X1,G_TRI_Y1,1);
    G_TRI_W = G_TRI_W1;
  }

  ublas::matrix<double> diffNTRI;
  ublas::matrix<double> g_NTRI3;
  const double *G_TRI_W;
  ublas::vector<double> coords;
  ublas::vector<double> normal0;

  ublas::vector<DofIdx> ent_global_col_indices,ent_global_row_indices;
  ublas::vector<double,ublas::bounded_array<double,9> > dofs_X;
  ublas::vector<double,ublas::bounded_array<double,3> > Lambda;

  PetscErrorCode getData(SideNumber_multiIndex::nth_index<1>::type::iterator siit) {
    PetscFunctionBegin;
    try {
    EntityHandle face = siit->ent;
    dofs_X.resize(9);
    Lambda.resize(3);
    ent_global_row_indices.resize(3);
    ent_global_col_indices.resize(9);
    const EntityHandle* conn_face; 
    int num_nodes; 
    rval = moab.get_connectivity(face,conn_face,num_nodes,true); CHKERR_PETSC(rval);
    if(num_nodes != 3) SETERRQ(PETSC_COMM_SELF,1,"face should have three nodes");
    for(int nn = 0;nn<num_nodes;nn++) {
      FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator dit,hi_dit;
      dit = row_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
      hi_dit = row_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple(lambda_field_name,conn_face[nn]));
      int NN = nn;
      if(siit->sense!=1) {
	if(nn == 1) NN = 2;
	if(nn == 2) NN = 1;
      } 
      try {
      if(distance(dit,hi_dit)==0) {
  	ent_global_row_indices[NN] = -1;
      } else {
  	if(distance(dit,hi_dit)!=1) SETERRQ1(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for < %s > should be 1",lambda_field_name.c_str());
	Lambda[NN] = dit->get_FieldData();
  	int global_idx = dit->get_petsc_gloabl_dof_idx();
  	ent_global_row_indices[NN] = global_idx;
  	int local_idx = dit->get_petsc_local_dof_idx();
  	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
      }
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      dit = col_multiIndex->get<Composite_mi_tag3>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
      hi_dit = col_multiIndex->get<Composite_mi_tag3>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",conn_face[nn]));
      if(distance(dit,hi_dit)!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, number of dof on node for MESH_NODE_POSITIONS should be 3");
      try {
      for(;dit!=hi_dit;dit++) {
  	dofs_X[NN*3+dit->get_dof_rank()] = dit->get_FieldData();
  	int global_idx = dit->get_petsc_gloabl_dof_idx();
  	assert(NN*3+dit->get_dof_rank()<9);
	if(ent_global_row_indices[NN] == -1) {
	  ent_global_col_indices[NN*3+dit->get_dof_rank()] = -1;
	} else {
	  ent_global_col_indices[NN*3+dit->get_dof_rank()] = global_idx;
	}
  	int local_idx = dit->get_petsc_local_dof_idx();
  	if(local_idx<0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, negative index of local dofs on element");
      }
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
    }
    coords.resize(num_nodes*3);
    rval = moab.get_coords(conn_face,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
    normal0.resize(3);
    ierr = ShapeFaceNormalMBTRI(&*diffNTRI.data().begin(),&coords.data()[0],&normal0.data()[0]); CHKERRQ(ierr);
    normal0 *= siit->sense;
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    };
    PetscFunctionReturn(0);
  }

  PetscErrorCode calulate_CONROL_AREA(double *diffNTRI,double *dofs_iX,double *C,double *iC) {
    PetscFunctionBegin;
    __CLPK_doublecomplex x_dofs_X[9];
    for(int nn = 0;nn<3;nn++) {
      for(int dd = 0;dd<3;dd++) {
        x_dofs_X[nn*3+dd].r = dofs_X[nn*3+dd];
        if(dofs_iX!=NULL) {
	  x_dofs_X[nn*3+dd].i = dofs_iX[nn*3+dd];
	} else {
	  x_dofs_X[nn*3+dd].i = 0;
	}
    }}
    __CLPK_doublecomplex x_normal[3];
    ierr = ShapeFaceNormalMBTRI_complex(diffNTRI,x_dofs_X,x_normal); CHKERRQ(ierr);
    /*for(int nn = 0;nn<3;nn++) {
      x_normal[nn].r *= snes;
      x_normal[nn].i *= snes;
    }*/
    double complex x_nrm2 = csqrt(
      cpow((x_normal[0].r+I*x_normal[0].i),2)+
      cpow((x_normal[1].r+I*x_normal[1].i),2)+
      cpow((x_normal[2].r+I*x_normal[2].i),2));
    double diffX_xi[3],diffX_eta[3];
    bzero(diffX_xi,3*sizeof(double));
    bzero(diffX_eta,3*sizeof(double));
    double i_diffX_xi[3],i_diffX_eta[3];
    bzero(i_diffX_xi,3*sizeof(double));
    bzero(i_diffX_eta,3*sizeof(double));
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
    bzero(SpinX_xi,9*sizeof(double));
    ierr = Spin(SpinX_xi,diffX_xi); CHKERRQ(ierr);
    bzero(SpinX_eta,9*sizeof(double));
    ierr = Spin(SpinX_eta,diffX_eta); CHKERRQ(ierr);
    double iSpinX_xi[9],iSpinX_eta[9];
    bzero(iSpinX_xi,9*sizeof(double));
    ierr = Spin(iSpinX_xi,i_diffX_xi); CHKERRQ(ierr);
    bzero(iSpinX_eta,9*sizeof(double));
    ierr = Spin(iSpinX_eta,i_diffX_eta); CHKERRQ(ierr);
    __CLPK_doublecomplex xSpinX_xi[9],xSpinX_eta[9];
    ierr = make_complex_matrix(SpinX_xi,iSpinX_xi,xSpinX_xi); CHKERRQ(ierr);
    ierr = make_complex_matrix(SpinX_eta,iSpinX_eta,xSpinX_eta); CHKERRQ(ierr);
   __CLPK_doublecomplex xNSpinX_xi[3],xNSpinX_eta[3];
    bzero(xNSpinX_xi,3*sizeof(__CLPK_doublecomplex));
    bzero(xNSpinX_eta,3*sizeof(__CLPK_doublecomplex));
    __CLPK_doublecomplex x_zero = { 0, 0 };
    __CLPK_doublecomplex x_scalar = { -creal(1./x_nrm2), -cimag(1./x_nrm2) };
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_xi,3,x_normal,1,&x_zero,xNSpinX_xi,1);
    cblas_zgemv(CblasRowMajor,CblasNoTrans,3,3,&x_scalar,xSpinX_eta,3,x_normal,1,&x_zero,xNSpinX_eta,1);
    for(int nn = 0;nn<3;nn++) {
      if(C != NULL) {
        C[3*nn + 0] = xNSpinX_xi[0].r*diffNTRI[2*nn+1]-xNSpinX_eta[0].r*diffNTRI[2*nn+0];
        C[3*nn + 1] = xNSpinX_xi[1].r*diffNTRI[2*nn+1]-xNSpinX_eta[1].r*diffNTRI[2*nn+0];
        C[3*nn + 2] = xNSpinX_xi[2].r*diffNTRI[2*nn+1]-xNSpinX_eta[2].r*diffNTRI[2*nn+0];
      }
      if(iC == NULL) continue;
      iC[3*nn + 0] = xNSpinX_xi[0].i*diffNTRI[2*nn+1]-xNSpinX_eta[0].i*diffNTRI[2*nn+0];
      iC[3*nn + 1] = xNSpinX_xi[1].i*diffNTRI[2*nn+1]-xNSpinX_eta[1].i*diffNTRI[2*nn+0];
      iC[3*nn + 2] = xNSpinX_xi[2].i*diffNTRI[2*nn+1]-xNSpinX_eta[2].i*diffNTRI[2*nn+0];
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
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));
    ublas::vector<double> ELEM_CONSTRAIN(9);
    for(;siit!=hi_siit;siit++) {
      EntityHandle face = siit->ent;
      if(find(surface.begin(),surface.end(),face) == surface.end()) continue;
      try {
	ierr = getData(siit); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
	  ostringstream ss;
	  ss << "thorw in method: " << ex.what() << endl;
	  SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      };
      try {
	ierr = calulate_CONROL_AREA(&*diffNTRI.data().begin(),NULL,&*ELEM_CONSTRAIN.data().begin(),NULL); CHKERRQ(ierr);
	for(int nn = 0;nn<3;nn++) {
	  if(ent_global_row_indices[nn]==-1) continue;
	  for(int NN = 0;NN<3;NN++) {
	    if(NN!=nn) continue;
	    if(Q == PETSC_NULL) {
	      ierr = MatSetValues(C,
		1,&(ent_global_row_indices.data()[nn]),
		3,&(ent_global_col_indices.data()[3*NN]),
		&ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	    }
	    if(Q != PETSC_NULL) {
	      if(mapV.find(ent_global_row_indices[nn])==mapV.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(mapV[ent_global_row_indices[nn]],3,&(ent_global_col_indices.data()[3*NN]),
		&ELEM_CONSTRAIN.data()[3*NN],ADD_VALUES); CHKERRQ(ierr);
	    }
	  }	
	}
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }; 
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
/*struct dC_CONSTANT_AREA_FEMethod: public  dC_CONSTANT_AREA_FEMethod {

  Mat dC;
  dC_CONSTANT_AREA_FEMethod(Interface& _moab,Range &_surface,Mat _dC,string _lambda_field_name,int _verbose = 0);


};*/


}

#endif // __MOABFEMETHOD_CONSTAREA_HPP__
