/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
* --------------------------------------------------------------
*
* DESCRIPTION: FIXME
*
* This is not exactly procedure for linear elatic dynamics, since jacobian is
* evaluated at every time step and snes procedure is involved. However it is
* implemented like that, to test methodology for general nonlinear problem.
*
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

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <ForcesAndSurcesCore.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

PetscErrorCode DataOperator::opLhs(
    DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data,bool symm) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  //nodes
  ierr = doWork(
    -1,-1,MBVERTEX,MBVERTEX,
    row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
  for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
    if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
    ierr = doWork(
      -1,VV,MBVERTEX,MBTET,
      row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
  }
  if(!symm) {
    for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
      if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
      ierr = doWork(
	-1,EE,MBVERTEX,MBEDGE,
	row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
    }
    for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
	-1,FF,MBVERTEX,MBTRI,
	row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
    }
  }

  //edges
  for(unsigned int ee = 0;ee<row_data.dataOnEntities[MBEDGE].size();ee++) {
    if(row_data.dataOnEntities[MBEDGE][ee].getN().size1()==0) continue;
    ierr = doWork(
	ee,-1,MBEDGE,MBVERTEX,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
	ee,VV,MBEDGE,MBTET,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    unsigned int EE = 0;
    if(symm) EE = ee;
    for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
      if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
      ierr = doWork(
	ee,EE,MBEDGE,MBEDGE,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
    }
    for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
	ee,FF,MBEDGE,MBTRI,
	row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
    }
  }

  //faces
  for(unsigned int ff = 0;ff<row_data.dataOnEntities[MBTRI].size();ff++) {
    if(row_data.dataOnEntities[MBTRI][ff].getN().size1()==0) continue;
    ierr = doWork(
	ff,-1,MBTRI,MBVERTEX,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
	ff,VV,MBTRI,MBTET,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    unsigned int FF = 0;
    if(symm) FF = ff;
    for(;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
	ff,FF,MBTRI,MBTRI,
	row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
    }
    if(!symm) {
      unsigned int EE = 0;
      for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
	if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
	ierr = doWork(
	  ff,EE,MBTRI,MBEDGE,
	  row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
      }
    }
  }

  //volumes
  for(unsigned int vv = 0;vv<row_data.dataOnEntities[MBTET].size();vv++) {
    if(row_data.dataOnEntities[MBTET][vv].getN().size1()==0) continue;
    unsigned int VV = 0;
    if(symm) VV = vv;
    for(;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      ierr = doWork(
	vv,VV,MBTET,MBTET,
	row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTET][VV]); CHKERRQ(ierr);
    }
    if(!symm) {
      //vertex
      ierr = doWork(
	VV,-1,MBTET,MBVERTEX,
	row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBVERTEX][0]); CHKERRQ(ierr);
      //edges
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
	if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
	ierr = doWork(
	  VV,EE,MBTET,MBEDGE,
	  row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBEDGE][EE]); CHKERRQ(ierr);
      }
      //faces
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
	if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
	ierr = doWork(
	  VV,FF,MBTET,MBTRI,
	  row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTRI][FF]); CHKERRQ(ierr);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DataOperator::opRhs(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  for(unsigned int nn = 0;nn<data.dataOnEntities[MBVERTEX].size();nn++) {
    ierr = doWork(-1,MBVERTEX,data.dataOnEntities[MBVERTEX][nn]); CHKERRQ(ierr);
  }
  for(unsigned int ee = 0;ee<data.dataOnEntities[MBEDGE].size();ee++) {
    //if(data.dataOnEntities[MBEDGE][ee].getN().size1()==0) continue;
    ierr = doWork(ee,MBEDGE,data.dataOnEntities[MBEDGE][ee]); CHKERRQ(ierr);
  }
  for(unsigned int ff = 0;ff<data.dataOnEntities[MBTRI].size();ff++) {
    //if(data.dataOnEntities[MBTRI][ff].getN().size1()==0) continue;
    ierr = doWork(ff,MBTRI,data.dataOnEntities[MBTRI][ff]); CHKERRQ(ierr);
  }
  for(unsigned int vv = 0;vv<data.dataOnEntities[MBTET].size();vv++) {
    //if(data.dataOnEntities[MBTET][vv].getN().size1()==0) continue;
    ierr = doWork(-1,MBTET,data.dataOnEntities[MBTET][vv]); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacH1::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  try {

    if(data.getDiffN().size2()==0) PetscFunctionReturn(0);

    diffNinvJac.resize(data.getDiffN().size1(),data.getDiffN().size2());
    unsigned int nb_gauss_pts = data.getN().size1();
    unsigned int nb_dofs = data.getN().size2();
    if(type!=MBVERTEX) {
      if(nb_dofs != data.getDiffN().size2()/3) {
        SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
  	"data inconsistency nb_dofs != data.diffN.size2()/3 ( %u != %u/3 )",
  	nb_dofs,data.getDiffN().size2());
      }
    }
  
    //cerr << type << endl;
    //cerr << data.getDiffN() << endl;
    //cerr << endl;

    switch (type) {
  
      case MBVERTEX: {
        ierr = ShapeDiffMBTETinvJ(
  	&*data.getDiffN().data().begin(),&*invJac.data().begin(),&*diffNinvJac.data().begin()); CHKERRQ(ierr);
      }
      break;
      case MBEDGE:
      case MBTRI:
      case MBTET: {
        for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
	  for(unsigned int dd = 0;dd<nb_dofs;dd++) {
	    cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	      &*invJac.data().begin(),3,&data.getDiffN()(gg,3*dd),1,0.,&diffNinvJac(gg,3*dd),1); 
	  }
        }
      }
      break;
      default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  
    }
  
    data.getDiffN().data().swap(diffNinvJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacHdiv::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

    diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2());

    unsigned int nb_gauss_pts = data.getDiffHdivN().size1();
    unsigned int nb_dofs = data.getDiffHdivN().size2()/9;
    
    unsigned int gg = 0;
    for(;gg<nb_gauss_pts;gg++) {
      unsigned int dd = 0;
      for(;dd<nb_dofs;dd++) {
	const double *DiffHdivN = &((data.getDiffHdivN(gg))(dd,0));
	for(int kk = 0;kk<3;kk++) {
	  cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,
	    &*invJac.data().begin(),3,&DiffHdivN[kk],3,
	    0.,&diffHdiv_invJac(gg,9*dd+kk),3); 
	}
      }
    }

    data.getDiffHdivN().data().swap(diffHdiv_invJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetPiolaTransform::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

  const double c = 1./6.;

  unsigned int nb_gauss_pts = data.getHdivN().size1();
  unsigned int nb_dofs = data.getHdivN().size2()/3;
  unsigned int gg = 0;
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2());
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2());
  for(;gg<nb_gauss_pts;gg++) {
    unsigned int dd = 0;
    for(;dd<nb_dofs;dd++) {
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,c/vOlume,
	&*Jac.data().begin(),3,&data.getHdivN()(gg,3*dd),1,0.,&piolaN(gg,3*dd),1);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,c/vOlume,
	  &*Jac.data().begin(),3,&data.getDiffHdivN()(gg,9*dd+3*kk),1,0.,&piolaDiffN(gg,9*dd+3*kk),1);
      }
    }
  }
  data.getHdivN().data().swap(piolaN.data());
  data.getDiffHdivN().data().swap(piolaDiffN.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetHoInvJacH1::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  try {

  if(data.getDiffN().size2()==0) PetscFunctionReturn(0);

  unsigned int nb_gauss_pts = data.getN().size1();
  unsigned int nb_dofs = data.getN().size2();
  //note Vetex diffN row has size of number of gass dof
  diffNinvJac.resize(nb_gauss_pts,3*nb_dofs);
  unsigned int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    double *inv_H = &invHoJac(gg,0);
    for(unsigned dd = 0;dd<nb_dofs;dd++) {
      double *diff_N;
      if(type == MBVERTEX) {
	diff_N = &data.getDiffN()(dd,0);
      } else {
	diff_N = &data.getDiffN()(gg,3*dd);
      }
      double *diff_N_inv_Jac = &diffNinvJac(gg,3*dd);
      cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_H,3,diff_N,1,0.,diff_N_inv_Jac,1); 
    }
  }

  if(type == MBVERTEX) {
    data.getDiffN().resize(diffNinvJac.size1(),diffNinvJac.size2());
  }
  data.getDiffN().data().swap(diffNinvJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetHoInvJacHdiv::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data)  {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

  diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2());

  unsigned int nb_gauss_pts = data.getDiffHdivN().size1();
  unsigned int nb_dofs = data.getDiffHdivN().size2()/9;

  unsigned int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    double *inv_h = &invHoJac(gg,0);
    for(unsigned dd = 0;dd<nb_dofs;dd++) {
      const double *diff_hdiv = &(data.getDiffHdivN(gg)(dd,0));
      double *diff_hdiv_inv_jac = &diffHdiv_invJac(gg,9*dd);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_h,3,&diff_hdiv[kk],3,0.,&diff_hdiv_inv_jac[kk],3); 
      }
    }
  }

  data.getDiffHdivN().data().swap(diffHdiv_invJac.data());

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetHoPiolaTransform::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try{

  unsigned int nb_gauss_pts = data.getHdivN().size1();
  unsigned int nb_dofs = data.getHdivN().size2()/3;
  unsigned int gg = 0;
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2());
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2());

  for(;gg<nb_gauss_pts;gg++) {
    unsigned int dd = 0;
    for(;dd<nb_dofs;dd++) {
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1./detHoJac[gg],
	&hoJac(gg,0),3,&data.getHdivN()(gg,3*dd),1,0.,&piolaN(gg,3*dd),1);
      int kk = 0;
      for(;kk<3;kk++) {
	cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,1./detHoJac[gg],
	  &hoJac(gg,0),3,&data.getDiffHdivN()(gg,9*dd+3*kk),1,0.,&piolaDiffN(gg,9*dd+3*kk),1);
      }
    }
  }

  data.getHdivN().data().swap(piolaN.data());
  data.getDiffHdivN().data().swap(piolaDiffN.data());


  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpGetData::doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  try {

  if(data.getFieldData().size() == 0) {
    PetscFunctionReturn(0);
  }

  unsigned int nb_dofs = data.getFieldData().size();
  if(nb_dofs % rank != 0) {
    SETERRQ4(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
      "data inconsistency, type %d, side %d, nb_dofs %d, rank %d",
      type,side,nb_dofs,rank);
  }
  if(nb_dofs/rank > data.getN().size2()) {
    SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,
      "data inconsistency nb_dofs >= data.N.size2() %u >= %u",nb_dofs,data.getN().size2());
  }
  data_at_GaussPt.resize(data.getN().size1(),rank);
  dataGrad_at_GaussPt.resize(data.getN().size1(),rank*dim);
  if(type == MBVERTEX) {
    bzero(&*data_at_GaussPt.data().begin(),data.getN().size1()*rank*sizeof(FieldData));
    bzero(&*dataGrad_at_GaussPt.data().begin(),data.getN().size1()*rank*dim*sizeof(FieldData));
    for(int rr = 0;rr<rank;rr++) {
      for(unsigned int dd = 0;dd<dim;dd++) {
	dataGrad_at_GaussPt(0,dim*rr+dd) = cblas_ddot(nb_dofs/rank,&data.getDiffN()(0,dd),dim,&data.getFieldData()[rr],rank);
      }
    }
  }
  for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
    for(int rr = 0;rr<rank;rr++) {
      data_at_GaussPt(gg,rr) = cblas_ddot(nb_dofs/rank,&data.getN()(gg,0),1,&data.getFieldData()[rr],rank);
      for(unsigned int dd = 0;dd<dim;dd++) {
	if(type == MBVERTEX) {
	  if(gg == 0) continue;
	  dataGrad_at_GaussPt(gg,dim*rr+dd) += dataGrad_at_GaussPt(0,dim*rr+dd);
	} else {
	  dataGrad_at_GaussPt(gg,dim*rr+dd) += cblas_ddot(nb_dofs/rank,&data.getDiffN()(gg,dd),dim,&data.getFieldData()[rr],rank);
	}
      }
    }
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  try {

  if(fePtr->get_ent_type() != MBTET) PetscFunctionReturn(0);

  ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

  //H1

  if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
    ierr = getTetsOrder(dataH1,H1); CHKERRQ(ierr);
    ierr = getFaceNodes(dataH1); CHKERRQ(ierr);
  }

  //Hdiv
  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
    ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
    ierr = getTetsOrder(dataHdiv,HDIV); CHKERRQ(ierr);
    ierr = getFaceNodes(dataHdiv); CHKERRQ(ierr);
  }

  //L2
  if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
    ierr = getTetsOrder(dataL2,L2); CHKERRQ(ierr);
  }

  int order = 1;
  for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
    order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
  }
  for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
    order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
  }
  order = max(order,dataL2.dataOnEntities[MBTET][0].getOrder());

  int nb_gauss_pts;
  int rule = getRule(order);
  if(rule >= 0) {
    nb_gauss_pts = gm_rule_size(rule,3);
    gaussPts.resize(4,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_3D_TET(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);
  } else {
    ierr = setGaussPts(order); CHKERRQ(ierr);
    nb_gauss_pts = gaussPts.size2();
  }

  ierr = shapeTETFunctions_H1(dataH1,
    &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);

  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = shapeTETFunctions_Hdiv(dataHdiv,
      &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
  }

  if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
    ierr = shapeTETFunctions_L2(dataL2,
      &gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
  }

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);
  vOlume = Shape_intVolumeMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.data().begin()); 
  Jac.resize(3,3);
  invJac.resize(3,3);
  ierr = ShapeJacMBTET(&*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),&*coords.begin(),&*Jac.data().begin()); CHKERRQ(ierr);
  noalias(invJac) = Jac;
  ierr = Shape_invJac(&*invJac.data().begin()); CHKERRQ(ierr);

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(4,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  try {
    ierr = opSetInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
    if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
      ierr = opSetInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
    }
    if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
      ierr = opPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
      ierr = opSetInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  if(mField.check_field(meshPositionsFieldName)) {
    BitFieldId id = mField.get_field_structure(meshPositionsFieldName)->get_id();
    if((fePtr->get_BitFieldId_data()&id).none()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
    }
    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTetsOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    if(dataH1.dataOnEntities[MBVERTEX][0].getFieldData().size()!=12) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no MESH_NODE_POSITIONS in element data");
    }
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTetsFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTetsFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opHOatGaussPoints.opRhs(dataH1); CHKERRQ(ierr);
      hoGaussPtsInvJac.resize(hoGaussPtsJac.size1(),hoGaussPtsJac.size2());
      ublas::noalias(hoGaussPtsInvJac) = hoGaussPtsJac;
      ublas::matrix<double> jac(3,3);
      hoGaussPtsDetJac.resize(nb_gauss_pts);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
	cblas_dcopy(9,&hoGaussPtsJac(gg,0),1,&jac(0,0),1);
	hoGaussPtsDetJac[gg] = Shape_detJac(&jac(0,0));
	ierr = Shape_invJac(&hoGaussPtsInvJac(gg,0)); CHKERRQ(ierr);
      }
      ierr = opSetHoInvJacH1.opRhs(dataH1); CHKERRQ(ierr);
      if((dataH1.spacesOnEntities[MBTET]).test(L2)) {
	ierr = opSetHoInvJacH1.opRhs(dataL2); CHKERRQ(ierr);
      }
      if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
	ierr = opSetHoPiolaTransform.opRhs(dataHdiv); CHKERRQ(ierr);
	ierr = opSetHoInvJacHdiv.opRhs(dataHdiv); CHKERRQ(ierr);
      }
    } catch (exception& ex) {
      ostringstream ss;
      ss << "problem with indices in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
  } else {
    ublas::matrix<double> diffN(nb_gauss_pts,12);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int nn = 0;nn<4;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  diffN(gg,nn*3+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
	}
      }
    }
    dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2());
    dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());
  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId data_id = mField.get_field_structure(oit->row_field_name)->get_id();
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();
    
    DataForcesAndSurcesCore *op_data = NULL;
    switch(row_space) {
      case H1:
	op_data = &dataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	op_data = &dataHdiv;
	break;
      case L2:
	op_data = &dataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpNN.begin();
    oit != vecUserOpNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();

    DataForcesAndSurcesCore *row_op_data = NULL;
    switch(row_space) {
      case H1:
	row_op_data = &dataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	row_op_data = &dataHdiv;
	break;
      case L2:
	row_op_data = &dataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    FieldSpace col_space = mField.get_field_structure(oit->col_field_name)->get_space();
    DataForcesAndSurcesCore *col_op_data = NULL;
    switch(col_space) {
      case H1:
	col_op_data = &derivedDataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	col_op_data = &derivedDataHdiv;
	break;
      case L2:
	col_op_data = &derivedDataL2;
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(col_space) {
      case H1:
      ierr = getColNodesIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HDIV:
      ierr = getTrisColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case L2:
      ierr = getTetsColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTetsOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTetsFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opLhs(*row_op_data,*col_op_data,oit->symm); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TetElementForcesAndSourcesCore::UserDataOperator::getDivergenceMatrixOperato_Hdiv(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data,
      int gg,ublas::vector<FieldData> &div) {
  PetscFunctionBegin;

  try {

  int nb_dofs = data.getFieldData().size();
  if((unsigned int)nb_dofs != data.getDiffHdivN().size2()/9) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  }

  if(nb_dofs == 0) PetscFunctionReturn(0);

  int dd = 0;
  for(;dd<nb_dofs;dd++) {
    div[dd] = 
      (data.getDiffHdivN(dd,gg))(0,0)+
      (data.getDiffHdivN(dd,gg))(1,1)+
      (data.getDiffHdivN(dd,gg))(2,2);
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetNormals::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;
  
  if(data.getFieldData().size()==0)  PetscFunctionReturn(0);
  
  switch (type) {
    case MBVERTEX: {
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	for(int nn = 0;nn<3;nn++) {
	  tAngent1_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,0),2,&data.getFieldData()[nn],3);
	  tAngent2_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,1),2,&data.getFieldData()[nn],3);
	}
      }
    } 
    break;
    case MBEDGE:     
    case MBTRI: {
      /*cerr << side << " " << type << endl;
      cerr << data.getN() << endl;
      cerr << data.getDiffN() << endl;
      cerr << data.getFieldData() << endl;
      cerr << "t1 " << tAngent1_at_GaussPt << endl;
      cerr << "t2 " << tAngent2_at_GaussPt << endl;*/
      if(2*data.getN().size2() != data.getDiffN().size2()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      unsigned int nb_dofs = data.getFieldData().size();
      if(nb_dofs%3!=0) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      if(nb_dofs > 3*data.getN().size2()) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      }
      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
	for(int dd = 0;dd<3;dd++) {
	  tAngent1_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
	  tAngent2_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
	}
      }
      //cerr << "t1 " << tAngent1_at_GaussPt << endl;
      //cerr << "t2 " << tAngent2_at_GaussPt << endl;
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetPiolaTransoformOnTriangle::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(type != MBTRI) PetscFunctionReturn(0);
	
  double l0 = cblas_dnrm2(3,&normal[0],1);
  int nb_gauss_pts = data.getHdivN().size1();
  int nb_dofs = data.getHdivN().size2()/3;
  int gg = 0;
  for(;gg<nb_gauss_pts;gg++) {
    
    int dd = 0;
    for(;dd<nb_dofs;dd++) {
      double val = data.getHdivN()(gg,3*dd);
      if(nOrmals_at_GaussPt.size1()==(unsigned int)nb_gauss_pts) {
	double l = cblas_dnrm2(3,&nOrmals_at_GaussPt(gg,0),1);
	cblas_dcopy(3,&nOrmals_at_GaussPt(gg,0),1,&data.getHdivN()(gg,3*dd),1);
	cblas_dscal(3,val/pow(l,2),&data.getHdivN()(gg,3*dd),1);
      } else {
	cblas_dcopy(3,&normal[0],1,&data.getHdivN()(gg,3*dd),1);
	cblas_dscal(3,val/pow(l0,2),&data.getHdivN()(gg,3*dd),1);
      }
    }    

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetNormals::calculateNormals() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  sPin.resize(3,3);
  bzero(&*sPin.data().begin(),9*sizeof(FieldData));
  nOrmals_at_GaussPt.resize(tAngent1_at_GaussPt.size1(),3);
  for(unsigned int gg = 0;gg<tAngent1_at_GaussPt.size1();gg++) {
    ierr = Spin(&*sPin.data().begin(),&tAngent1_at_GaussPt(gg,0)); CHKERRQ(ierr);
    cblas_dgemv(
      CblasRowMajor,CblasNoTrans,3,3,1.,
      &*sPin.data().begin(),3,&tAngent2_at_GaussPt(gg,0),1,0.,
      &nOrmals_at_GaussPt(gg,0),1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode TriElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBTRI) PetscFunctionReturn(0);

  ierr = getSpacesOnEntities(dataH1); CHKERRQ(ierr);

  //H1

  if((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    ierr = getEdgesSense(dataH1); CHKERRQ(ierr);
    ierr = getTrisSense(dataH1); CHKERRQ(ierr);
    ierr = getEdgesOrder(dataH1,H1); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,H1); CHKERRQ(ierr);
  }

  //Hdiv
  if((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    ierr = getTrisSense(dataHdiv); CHKERRQ(ierr);
    ierr = getTrisOrder(dataHdiv,HDIV); CHKERRQ(ierr);
  }

  int order = 1;
  for(unsigned int ee = 0;ee<dataH1.dataOnEntities[MBEDGE].size();ee++) {
    order = max(order,dataH1.dataOnEntities[MBEDGE][ee].getOrder());
  }
  for(unsigned int ff = 0;ff<dataHdiv.dataOnEntities[MBTRI].size();ff++) {
    order = max(order,dataHdiv.dataOnEntities[MBTRI][ff].getOrder());
  }

  int nb_gauss_pts;
  int rule = getRule(order);
  if(rule >= 0) {
    nb_gauss_pts = gm_rule_size(rule,2);
    gaussPts.resize(3,nb_gauss_pts);
    ierr = Grundmann_Moeller_integration_points_2D_TRI(
      rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0)); CHKERRQ(ierr);
  } else {
    ierr = setGaussPts(order); CHKERRQ(ierr);
    nb_gauss_pts = gaussPts.size2();
  }

  ierr = shapeTRIFunctions_H1(dataH1,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr);

  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = shapeTRIFunctions_Hdiv(dataHdiv,&gaussPts(0,0),&gaussPts(1,0),nb_gauss_pts); CHKERRQ(ierr); CHKERRQ(ierr);
  }

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

  normal.resize(3);
  ierr = ShapeFaceNormalMBTRI(
    &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
    &*coords.data().begin(),&*normal.data().begin()); CHKERRQ(ierr);
  aRea = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) = cblas_ddot(3,&dataH1.dataOnEntities[MBVERTEX][0].getN()(gg,0),1,&coords[dd],3);
    }
  }

  if(mField.check_field(meshPositionsFieldName)) {
    nOrmals_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent1_at_GaussPt.resize(nb_gauss_pts,3);
    tAngent2_at_GaussPt.resize(nb_gauss_pts,3);
    ierr = getEdgesOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisOrder(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldData(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    ierr = getTrisFieldDofs(dataH1,meshPositionsFieldName); CHKERRQ(ierr);
    try {
      ierr = opHONormals.opRhs(dataH1); CHKERRQ(ierr);
      ierr = opHONormals.calculateNormals(); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
  } else {
    ublas::matrix<double> diffN(nb_gauss_pts,6);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<2;dd++) {
	  diffN(gg,nn*2+dd) = dataH1.dataOnEntities[MBVERTEX][0].getDiffN()(nn,dd);
	}
      }
    }
    dataH1.dataOnEntities[MBVERTEX][0].getDiffN().resize(diffN.size1(),diffN.size2());
    dataH1.dataOnEntities[MBVERTEX][0].getDiffN().data().swap(diffN.data());
  }

  if(dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    ierr = opSetPiolaTransoformOnTriangle.opRhs(dataHdiv); CHKERRQ(ierr);
  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId data_id = mField.get_field_structure(oit->row_field_name)->get_id();
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&data_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();
    
    DataForcesAndSurcesCore *op_data = NULL;
    switch(row_space) {
      case H1:
	op_data = &dataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	op_data = &dataHdiv;
	break;
      case L2:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisRowIndices(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opRhs(*op_data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    FieldSpace row_space = mField.get_field_structure(oit->row_field_name)->get_space();

    DataForcesAndSurcesCore *row_op_data = NULL;
    switch(row_space) {
      case H1:
	row_op_data = &dataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	row_op_data = &dataHdiv;
	break;
      case L2:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(row_space) {
      case H1:
      ierr = getRowNodesIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisRowIndices(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*row_op_data,oit->row_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    FieldSpace col_space = mField.get_field_structure(oit->col_field_name)->get_space();
    DataForcesAndSurcesCore *col_op_data = NULL;
    switch(col_space) {
      case H1:
	col_op_data = &derivedDataH1;
	break;
      case HCURL:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      case HDIV:
	col_op_data = &derivedDataHdiv;
	break;
      case L2:
	SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented yet");
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
      break;
    }

    switch(col_space) {
      case H1:
      ierr = getColNodesIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getNodesFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HCURL:
      ierr = getEdgesColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getEdgesFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      case HDIV:
      case L2:
      ierr = getTrisColIndices(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisOrder(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldData(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      ierr = getTrisFieldDofs(*col_op_data,oit->col_field_name); CHKERRQ(ierr);
      default:
      break;
    }

    try {
      ierr = oit->opLhs(*row_op_data,*col_op_data,oit->symm); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode EdgeElementForcesAndSurcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBEDGE) PetscFunctionReturn(0);

  PetscErrorCode ierr;

  //PetscAttachDebugger();

  ierr = getEdgesOrder(data,H1); CHKERRQ(ierr);

  int order = data.dataOnEntities[MBEDGE][0].getOrder();
  int rule = getRule(order);
  int nb_gauss_pts = gm_rule_size(rule,1);
  gaussPts.resize(2,nb_gauss_pts);

  ierr = Grundmann_Moeller_integration_points_1D_EDGE(rule,&gaussPts(0,0),&gaussPts(1,0)); CHKERRQ(ierr);
  ierr = shapeEDGEFunctions_H1(data,&gaussPts(0,0),nb_gauss_pts); CHKERRQ(ierr);

  EntityHandle ent = fePtr->get_ent();
  int num_nodes;
  const EntityHandle* conn;
  rval = mField.get_moab().get_connectivity(ent,conn,num_nodes,true); CHKERR_PETSC(rval);
  coords.resize(num_nodes*3);
  rval = mField.get_moab().get_coords(conn,num_nodes,&*coords.data().begin()); CHKERR_PETSC(rval);

  dIrection.resize(3);
  cblas_dcopy(3,&coords[3],1,&*dIrection.data().begin(),1);
  cblas_daxpy(3,-1.,&coords[0],1,&*dIrection.data().begin(),1);
  lEngth = cblas_dnrm2(3,&*dIrection.data().begin(),1);

  coordsAtGaussPts.resize(nb_gauss_pts,3);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    for(int dd = 0;dd<3;dd++) {
      coordsAtGaussPts(gg,dd) 
	= N_MBEDGE0(gaussPts(0,gg))*coords[dd] + N_MBEDGE1(gaussPts(0,gg))*coords[3+dd]; 
    }
  }
  //cerr << coordsAtGaussPts << endl;

  DataForcesAndSurcesCore *col_data = &derivedData;

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    //row indices
    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getEdgesRowIndices(data,oit->row_field_name); CHKERRQ(ierr);
    //col data
    ierr = getEdgesOrder(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opRhs(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_col()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    //row indices
    ierr = getEdgesOrder(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getEdgesRowIndices(data,oit->row_field_name); CHKERRQ(ierr);
    //col indices
    ierr = getEdgesOrder(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getColNodesIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    //col data
    ierr = getEdgesColIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getEdgesFieldDofs(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opLhs(data,*col_data,true); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode VertexElementForcesAndSourcesCore::operator()() {
  PetscFunctionBegin;

  if(fePtr->get_ent_type() != MBVERTEX) PetscFunctionReturn(0);

  EntityHandle ent = fePtr->get_ent();
  coords.resize(3);
  rval = mField.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERR_PETSC(rval);

  DataForcesAndSurcesCore *col_data = &derivedData;

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpN.begin();
    oit != vecUserOpN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->row_field_name.c_str());
    }

    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opRhs(data); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  for(
    boost::ptr_vector<UserDataOperator>::iterator oit = vecUserOpSymmNN.begin();
    oit != vecUserOpSymmNN.end(); oit++) {

    oit->setPtrFE(this);
    BitFieldId row_id = mField.get_field_structure(oit->row_field_name)->get_id();
    BitFieldId col_id = mField.get_field_structure(oit->col_field_name)->get_id();

    if((oit->getMoFEMFEPtr()->get_BitFieldId_row()&row_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no row field < %s > on finite elemeny",oit->row_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_col()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }
    if((oit->getMoFEMFEPtr()->get_BitFieldId_data()&col_id).none()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no data field < %s > on finite elemeny",oit->col_field_name.c_str());
    }

    ierr = getRowNodesIndices(data,oit->row_field_name); CHKERRQ(ierr);
    ierr = getColNodesIndices(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldData(*col_data,oit->col_field_name); CHKERRQ(ierr);
    ierr = getNodesFieldDofs(*col_data,oit->col_field_name); CHKERRQ(ierr);

    try {
      ierr = oit->opLhs(data,*col_data,true); CHKERRQ(ierr);
    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

  }

  PetscFunctionReturn(0);
}



}
