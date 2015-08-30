/** file DataOperators.cpp

  \brief implementation of Data Operators for Forces and Sources

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
  for(unsigned int nn = 0;nn!=row_data.dataOnEntities[MBVERTEX].size();nn++) {
    unsigned int NN = 0;
    if(symm) NN = nn;
    for(;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
      ierr = doWork(
        nn,NN,MBVERTEX,MBVERTEX,
        row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBVERTEX][0]
      ); CHKERRQ(ierr);
    }
    if(!symm) {
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          nn,EE,MBVERTEX,MBEDGE,
          row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
        if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
        ierr = doWork(
          nn,FF,MBVERTEX,MBTRI,
          row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTRI][FF]
        ); CHKERRQ(ierr);
      }
    }
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        nn,VV,MBVERTEX,MBTET,
        row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        row_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        row_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        nn,MM,MBVERTEX,MBENTITYSET,
        row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBENTITYSET][MM]
      ); CHKERRQ(ierr);
    }
  }

  //edges
  for(unsigned int ee = 0;ee<row_data.dataOnEntities[MBEDGE].size();ee++) {
    if(row_data.dataOnEntities[MBEDGE][ee].getN().size1()==0) continue;
    for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
      ierr = doWork(
        ee,NN,MBEDGE,MBVERTEX,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBVERTEX][0]
      ); CHKERRQ(ierr);
    }
    unsigned int EE = 0;
    if(symm) EE = ee;
    for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
      if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
      ierr = doWork(
        ee,EE,MBEDGE,MBEDGE,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBEDGE][EE]
      ); CHKERRQ(ierr);
    }
    for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
        ee,FF,MBEDGE,MBTRI,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTRI][FF]
      ); CHKERRQ(ierr);
    }
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        ee,VV,MBEDGE,MBTET,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        col_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        col_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        ee,MM,MBEDGE,MBENTITYSET,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBENTITYSET][MM]
      ); CHKERRQ(ierr);
    }
  }

  //faces
  for(unsigned int ff = 0;ff<row_data.dataOnEntities[MBTRI].size();ff++) {
    if(row_data.dataOnEntities[MBTRI][ff].getN().size1()==0) continue;
    for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
      ierr = doWork(
        ff,NN,MBTRI,MBVERTEX,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBVERTEX][0]
      ); CHKERRQ(ierr);
    }
    if(!symm) {
      unsigned int EE = 0;
      for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          ff,EE,MBTRI,MBEDGE,
          row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
    }
    unsigned int FF = 0;
    if(symm) FF = ff;
    for(;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
        ff,FF,MBTRI,MBTRI,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTRI][FF]
      ); CHKERRQ(ierr);
    }
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        ff,VV,MBTRI,MBTET,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        col_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        col_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        ff,MM,MBTRI,MBENTITYSET,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBENTITYSET][MM]
      ); CHKERRQ(ierr);
    }
  }

  //volumes
  for(unsigned int vv = 0;vv<row_data.dataOnEntities[MBTET].size();vv++) {
    if(row_data.dataOnEntities[MBTET][vv].getN().size1()==0) continue;
    if(!symm) {
      //vertex
      for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
        ierr = doWork(
          vv,NN,MBTET,MBVERTEX,
          row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBVERTEX][0]
        ); CHKERRQ(ierr);
      }
      //edges
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          vv,EE,MBTET,MBEDGE,
          row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
      //faces
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
        if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
        ierr = doWork(
          vv,FF,MBTET,MBTRI,
          row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTRI][FF]
        ); CHKERRQ(ierr);
      }
    }
    unsigned int VV = 0;
    if(symm) VV = vv;
    for(;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        vv,VV,MBTET,MBTET,
        row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        col_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        col_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        vv,MM,MBTET,MBENTITYSET,
        row_data.dataOnEntities[MBTET][vv],col_data.dataOnEntities[MBENTITYSET][MM]
      ); CHKERRQ(ierr);
    }
  }

  //meshsets
  for(unsigned int mm = 0;mm<row_data.dataOnEntities[MBENTITYSET].size();mm++) {
    if(
      row_data.dataOnEntities[MBENTITYSET][mm].getIndices().empty()&&
      row_data.dataOnEntities[MBENTITYSET][mm].getFieldData().empty()
    ) continue;
    if(!symm) {
      //vertex
      for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
        ierr = doWork(
          mm,NN,MBENTITYSET,MBVERTEX,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBVERTEX][0]
        ); CHKERRQ(ierr);
      }
      //edges
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          mm,EE,MBENTITYSET,MBEDGE,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
      //faces
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
        if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
        ierr = doWork(
          mm,FF,MBENTITYSET,MBTRI,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBTRI][FF]
        ); CHKERRQ(ierr);
      }
      //volume
      for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
        ierr = doWork(
          mm,VV,MBENTITYSET,MBTET,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBTET][VV]
        ); CHKERRQ(ierr);
      }
    }
    unsigned int MM = 0;
    if(symm) MM = mm;
    for(;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        row_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        row_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        mm,MM,MBENTITYSET,MBENTITYSET,
        row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBENTITYSET][MM]
      ); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DataOperator::opRhs(DataForcesAndSurcesCore &data) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  for(unsigned int nn = 0;nn<data.dataOnEntities[MBVERTEX].size();nn++) {
    ierr = doWork(nn,MBVERTEX,data.dataOnEntities[MBVERTEX][nn]); CHKERRQ(ierr);
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
    ierr = doWork(vv,MBTET,data.dataOnEntities[MBTET][vv]); CHKERRQ(ierr);
  }
  for(unsigned int mm = 0;mm<data.dataOnEntities[MBENTITYSET].size();mm++) {
    if(
        data.dataOnEntities[MBENTITYSET][mm].getIndices().empty()&&
        data.dataOnEntities[MBENTITYSET][mm].getFieldData().empty()
    ) continue;
    ierr = doWork(mm,MBENTITYSET,data.dataOnEntities[MBENTITYSET][mm]); CHKERRQ(ierr);
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

    diffNinvJac.resize(data.getDiffN().size1(),data.getDiffN().size2(),false);
    unsigned int nb_gauss_pts = data.getN().size1();
    unsigned int nb_dofs = data.getN().size2();
    if(type!=MBVERTEX) {
      if(nb_dofs != data.getDiffN().size2()/3) {
        SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
          "data inconsistency nb_dofs != data.diffN.size2()/3 ( %u != %u/3 )",
          nb_dofs,data.getDiffN().size2()
        );
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
                &*invJac.data().begin(),3,&data.getDiffN()(gg,3*dd),1,0.,&diffNinvJac(gg,3*dd),1
              );
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

    diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2(),false);

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
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2(),false);
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2(),false);
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
    // Note for Vetex diffN row has size of number of dof
    diffNinvJac.resize(nb_gauss_pts,3*nb_dofs,false);

    unsigned int gg = 0;
    for(;gg<nb_gauss_pts;gg++) {
      double *inv_h = &invHoJac(gg,0);
      for(unsigned dd = 0;dd<nb_dofs;dd++) {
        double *diff_n;
        if(type == MBVERTEX) {
          diff_n = &data.getDiffN()(dd,0);
        } else {
          diff_n = &data.getDiffN()(gg,3*dd);
        }
        double *diff_n_inv_jac = &diffNinvJac(gg,3*dd);
        cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_h,3,diff_n,1,0.,diff_n_inv_jac,1);
      }
    }

    if(type == MBVERTEX) {
      data.getDiffN().resize(diffNinvJac.size1(),diffNinvJac.size2(),false);
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

  diffHdiv_invJac.resize(data.getDiffHdivN().size1(),data.getDiffHdivN().size2(),false);

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
  piolaN.resize(nb_gauss_pts,data.getHdivN().size2(),false);
  piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN().size2(),false);

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


PetscErrorCode OpGetDataAndGradient::doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  try {

    if(data.getFieldData().size() == 0) {
      PetscFunctionReturn(0);
    }

    unsigned int nb_dofs = data.getFieldData().size();
    if(nb_dofs % rank != 0) {
      SETERRQ4(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
        "data inconsistency, type %d, side %d, nb_dofs %d, rank %d",
        type,side,nb_dofs,rank
      );
    }
    if(nb_dofs/rank > data.getN().size2()) {
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,
        "data inconsistency nb_dofs >= data.N.size2() %u >= %u",nb_dofs,data.getN().size2()
      );
    }

    data_at_GaussPt.resize(data.getN().size1(),rank,false);
    dataGrad_at_GaussPt.resize(data.getN().size1(),rank*dim,false);

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
      double *data_ptr,*n_ptr,*diff_n_ptr;
      n_ptr = &data.getN()(gg,0);
      data_ptr = &data.getFieldData()[0];
      if(type != MBVERTEX) {
        diff_n_ptr  = &data.getDiffN()(gg,0);
      } else {
        diff_n_ptr  = &data.getDiffN()(0,0);
      }

      for(int rr = 0;rr<rank;rr++,data_ptr++) {
        data_at_GaussPt(gg,rr) += cblas_ddot(nb_dofs/rank,n_ptr,1,data_ptr,rank);
        double *diff_n_ptr2 = diff_n_ptr;

        for(unsigned int dd = 0;dd<dim;dd++,diff_n_ptr2++) {
          if(type == MBVERTEX) {
            if(gg == 0) continue;
            dataGrad_at_GaussPt(gg,dim*rr+dd) += dataGrad_at_GaussPt(0,dim*rr+dd);
          } else {
            dataGrad_at_GaussPt(gg,dim*rr+dd) += cblas_ddot(nb_dofs/rank,diff_n_ptr2,dim,data_ptr,rank);
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


PetscErrorCode OpGetCoordsAndNormalsOnFace::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  unsigned int nb_dofs = data.getFieldData().size();
  if(nb_dofs==0)  PetscFunctionReturn(0);

  int nb_gauss_pts = data.getN().size1();
  cOords_at_GaussPt.resize(nb_gauss_pts,3,false);
  nOrmals_at_GaussPt.resize(nb_gauss_pts,3,false);
  tAngent1_at_GaussPt.resize(nb_gauss_pts,3,false);
  tAngent2_at_GaussPt.resize(nb_gauss_pts,3,false);

  switch (type) {
    case MBVERTEX: {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int nn = 0;nn<3;nn++) {
          cOords_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getN(gg)[0],1,&data.getFieldData()[nn],3);
          tAngent1_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,0),2,&data.getFieldData()[nn],3);
          tAngent2_at_GaussPt(gg,nn) = cblas_ddot(3,&data.getDiffN()(0,1),2,&data.getFieldData()[nn],3);
        }
      }
    }
    break;
    case MBEDGE:
    case MBTRI: {
      if(2*data.getN().size2() != data.getDiffN().size2()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      if(nb_dofs%3!=0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      if(nb_dofs > 3*data.getN().size2()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int dd = 0;dd<3;dd++) {
          cOords_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
          tAngent1_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
          tAngent2_at_GaussPt(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
        }
      }
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetCoordsAndNormalsOnFace::calculateNormals() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  sPin.resize(3,3);
  sPin.clear();
  nOrmals_at_GaussPt.resize(tAngent1_at_GaussPt.size1(),3,false);
  for(unsigned int gg = 0;gg<tAngent1_at_GaussPt.size1();gg++) {
    ierr = Spin(&*sPin.data().begin(),&tAngent1_at_GaussPt(gg,0)); CHKERRQ(ierr);
    cblas_dgemv(
      CblasRowMajor,CblasNoTrans,3,3,1.,
      &*sPin.data().begin(),3,&tAngent2_at_GaussPt(gg,0),1,0.,
      &nOrmals_at_GaussPt(gg,0),1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetCoordsAndNormalsOnPrism::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  if(data.getFieldData().size()==0)  PetscFunctionReturn(0);
  const int valid_edges3[] = { 1,1,1,0,0,0,0,0,0 };
  const int valid_faces3[] = { 0,0,0,1,0,0,0,0,0 };
  const int valid_edges4[] = { 0,0,0,0,0,0,1,1,1 };
  const int valid_faces4[] = { 0,0,0,0,1,0,0,0,0 };

  try {

    switch (type) {
      case MBVERTEX: {
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          for(int dd = 0;dd<3;dd++) {
            cOords_at_GaussPtF3(gg,dd) = cblas_ddot(3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
            tAngent1_at_GaussPtF3(gg,dd) = cblas_ddot(3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
            tAngent2_at_GaussPtF3(gg,dd) = cblas_ddot(3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
            cOords_at_GaussPtF4(gg,dd) = cblas_ddot(3,&data.getN(gg)[0],1,&data.getFieldData()[9+dd],3);
            tAngent1_at_GaussPtF4(gg,dd) = cblas_ddot(3,&data.getDiffN()(gg,6+0),2,&data.getFieldData()[9+dd],3);
            tAngent2_at_GaussPtF4(gg,dd) = cblas_ddot(3,&data.getDiffN()(gg,6+1),2,&data.getFieldData()[9+dd],3);
          }
        }
      }
      break;
      case MBEDGE:
      case MBTRI: {
        if(2*data.getN().size2() != data.getDiffN().size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }
        unsigned int nb_dofs = data.getFieldData().size();
        if(nb_dofs%3!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }
        if(nb_dofs > 3*data.getN().size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
        }
        for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
          for(int dd = 0;dd<3;dd++) {
            if ((type == MBTRI && valid_faces3[side]) || (type == MBEDGE && valid_edges3[side]))  {
              cOords_at_GaussPtF3(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
              tAngent1_at_GaussPtF3(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
              tAngent2_at_GaussPtF3(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
            } else if((type == MBTRI && valid_faces4[side]) || (type == MBEDGE && valid_edges4[side])) {
              cOords_at_GaussPtF4(gg,dd) += cblas_ddot(nb_dofs/3,&data.getN(gg)[0],1,&data.getFieldData()[dd],3);
              tAngent1_at_GaussPtF4(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,0),2,&data.getFieldData()[dd],3);
              tAngent2_at_GaussPtF4(gg,dd) += cblas_ddot(nb_dofs/3,&data.getDiffN()(gg,1),2,&data.getFieldData()[dd],3);
            }
          }
        }
      }
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetCoordsAndNormalsOnPrism::calculateNormals() {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  try {
    sPin.resize(3,3);
    sPin.clear();
    nOrmals_at_GaussPtF3.resize(tAngent1_at_GaussPtF3.size1(),3,false);
    for(unsigned int gg = 0;gg<tAngent1_at_GaussPtF3.size1();gg++) {
      ierr = Spin(&*sPin.data().begin(),&tAngent1_at_GaussPtF3(gg,0)); CHKERRQ(ierr);
      cblas_dgemv(
        CblasRowMajor,CblasNoTrans,3,3,1.,
        &*sPin.data().begin(),3,&tAngent2_at_GaussPtF3(gg,0),1,0.,
        &nOrmals_at_GaussPtF3(gg,0),1
      );
    }
    sPin.clear();
    nOrmals_at_GaussPtF4.resize(tAngent1_at_GaussPtF4.size1(),3,false);
    for(unsigned int gg = 0;gg<tAngent1_at_GaussPtF4.size1();gg++) {
      ierr = Spin(&*sPin.data().begin(),&tAngent1_at_GaussPtF4(gg,0)); CHKERRQ(ierr);
      cblas_dgemv(
        CblasRowMajor,CblasNoTrans,3,3,1.,
        &*sPin.data().begin(),3,&tAngent2_at_GaussPtF4(gg,0),1,0.,
        &nOrmals_at_GaussPtF4(gg,0),1
      );
    }
  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetPiolaTransoformOnTriangle::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
) {
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

PetscErrorCode OpGetHoTangentOnEdge::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if(nb_dofs == 0)  PetscFunctionReturn(0);

  try {

    int nb_gauss_pts = data.getN().size1();
    tAngent.resize(nb_gauss_pts,3,false);

    int nb_approx_fun = data.getN().size2();
    double *diff = &*data.getDiffN().data().begin();
    double *dofs[] = { &data.getFieldData()[0], &data.getFieldData()[1], &data.getFieldData()[2] };

    tAngent.resize(nb_gauss_pts,3,false);

    switch(type) {
      case MBVERTEX:
      for(int dd = 0;dd!=3;dd++) {
        for(int gg = 0;gg!=nb_gauss_pts;gg++) {
          tAngent(gg,dd) = cblas_ddot(2,diff,1,dofs[dd],3);
        }
      }
      break;
      case MBEDGE:
        if(nb_dofs%3) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"Approximated field should be rank 3, i.e. vector in 3d space");
        }
        for(int dd = 0;dd!=3;dd++) {
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            tAngent(gg,dd) += cblas_ddot(nb_dofs/3,&diff[gg*nb_approx_fun],1,dofs[dd],3);
          }
        }
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"This operator can calculate tangent vector only on edge");
    }


  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

}
