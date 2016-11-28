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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>

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
      for(unsigned int QQ = 0;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
        if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
        ierr = doWork(
          nn,QQ,MBVERTEX,MBQUAD,
          row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBQUAD][QQ]
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
    for(unsigned int PP = 0;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
      if(col_data.dataOnEntities[MBPRISM][PP].getN().size1()==0) continue;
      ierr = doWork(
        nn,PP,MBVERTEX,MBPRISM,
        row_data.dataOnEntities[MBVERTEX][0],col_data.dataOnEntities[MBPRISM][PP]
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
    //tris
    for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
      if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
      ierr = doWork(
        ee,FF,MBEDGE,MBTRI,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTRI][FF]
      ); CHKERRQ(ierr);
    }
    //quad
    for(unsigned int QQ = 0;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
      if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
      ierr = doWork(
        ee,QQ,MBEDGE,MBQUAD,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBQUAD][QQ]
      ); CHKERRQ(ierr);
    }
    //tet
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        ee,VV,MBEDGE,MBTET,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    //prism
    for(unsigned int PP = 0;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
      if(col_data.dataOnEntities[MBPRISM][PP].getN().size1()==0) continue;
      ierr = doWork(
        ee,PP,MBEDGE,MBPRISM,
        row_data.dataOnEntities[MBEDGE][ee],col_data.dataOnEntities[MBPRISM][PP]
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
    for(unsigned int QQ = 0;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
      if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
      ierr = doWork(
        ff,QQ,MBTRI,MBQUAD,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBQUAD][QQ]
      ); CHKERRQ(ierr);
    }
    for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
      if(col_data.dataOnEntities[MBTET][VV].getN().size1()==0) continue;
      ierr = doWork(
        ff,VV,MBTRI,MBTET,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBTET][VV]
      ); CHKERRQ(ierr);
    }
    for(unsigned int PP = 0;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
      if(col_data.dataOnEntities[MBPRISM][PP].getN().size1()==0) continue;
      ierr = doWork(
        ff,PP,MBTRI,MBPRISM,
        row_data.dataOnEntities[MBTRI][ff],col_data.dataOnEntities[MBPRISM][PP]
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

  //quads
  for(unsigned int qq = 0;qq<row_data.dataOnEntities[MBQUAD].size();qq++) {
    if(row_data.dataOnEntities[MBQUAD][qq].getN().size1()==0) continue;
    for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
      ierr = doWork(
        qq,NN,MBQUAD,MBVERTEX,
        row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBVERTEX][0]
      ); CHKERRQ(ierr);
    }
    if(!symm) {
      unsigned int EE = 0;
      for(;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          qq,EE,MBQUAD,MBEDGE,
          row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
      unsigned int FF = 0;
      for(;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
        if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
        ierr = doWork(
          qq,FF,MBQUAD,MBTRI,
          row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBTRI][FF]
        ); CHKERRQ(ierr);
      }
    }
    unsigned int QQ = 0;
    if(symm) QQ = qq;
    for(;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
      if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
      ierr = doWork(
        qq,QQ,MBQUAD,MBQUAD,
        row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBQUAD][QQ]
      ); CHKERRQ(ierr);
    }
    for(unsigned int PP = 0;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
      if(col_data.dataOnEntities[MBPRISM][PP].getN().size1()==0) continue;
      ierr = doWork(
        qq,PP,MBTRI,MBPRISM,
        row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBPRISM][PP]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        col_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        col_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        qq,MM,MBQUAD,MBENTITYSET,
        row_data.dataOnEntities[MBQUAD][qq],col_data.dataOnEntities[MBENTITYSET][MM]
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

  for(unsigned int pp = 0;pp<row_data.dataOnEntities[MBPRISM].size();pp++) {
    if(row_data.dataOnEntities[MBPRISM][pp].getN().size1()==0) continue;
    if(!symm) {
      //vertex
      for(unsigned int NN = 0;NN!=col_data.dataOnEntities[MBVERTEX].size();NN++) {
        ierr = doWork(
          pp,NN,MBPRISM,MBVERTEX,
          row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBVERTEX][0]
        ); CHKERRQ(ierr);
      }
      //edges
      for(unsigned int EE = 0;EE<col_data.dataOnEntities[MBEDGE].size();EE++) {
        if(col_data.dataOnEntities[MBEDGE][EE].getN().size1()==0) continue;
        ierr = doWork(
          pp,EE,MBPRISM,MBEDGE,
          row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBEDGE][EE]
        ); CHKERRQ(ierr);
      }
      //faces
      for(unsigned int FF = 0;FF<col_data.dataOnEntities[MBTRI].size();FF++) {
        if(col_data.dataOnEntities[MBTRI][FF].getN().size1()==0) continue;
        ierr = doWork(
          pp,FF,MBPRISM,MBTRI,
          row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBTRI][FF]
        ); CHKERRQ(ierr);
      }
      //quads
      for(unsigned int QQ = 0;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
        if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
        ierr = doWork(
          pp,QQ,MBPRISM,MBQUAD,
          row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBQUAD][QQ]
        ); CHKERRQ(ierr);
      }
    }
    unsigned int PP = 0;
    if(symm) PP = pp;
    for(;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
      if(col_data.dataOnEntities[MBPRISM][PP].getN().size1()==0) continue;
      ierr = doWork(
        pp,PP,MBPRISM,MBPRISM,
        row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBPRISM][PP]
      ); CHKERRQ(ierr);
    }
    for(unsigned int MM = 0;MM<col_data.dataOnEntities[MBENTITYSET].size();MM++) {
      if(
        col_data.dataOnEntities[MBENTITYSET][MM].getIndices().empty()&&
        col_data.dataOnEntities[MBENTITYSET][MM].getFieldData().empty()
      ) continue;
      ierr = doWork(
        pp,MM,MBPRISM,MBENTITYSET,
        row_data.dataOnEntities[MBPRISM][pp],col_data.dataOnEntities[MBENTITYSET][MM]
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
      //quad
      for(unsigned int QQ = 0;QQ<col_data.dataOnEntities[MBQUAD].size();QQ++) {
        if(col_data.dataOnEntities[MBQUAD][QQ].getN().size1()==0) continue;
        ierr = doWork(
          mm,QQ,MBENTITYSET,MBQUAD,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBQUAD][QQ]
        ); CHKERRQ(ierr);
      }
      //volume
      for(unsigned int VV = 0;VV<col_data.dataOnEntities[MBTET].size();VV++) {
        ierr = doWork(
          mm,VV,MBENTITYSET,MBTET,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBTET][VV]
        ); CHKERRQ(ierr);
      }
      for(unsigned int PP = 0;PP<col_data.dataOnEntities[MBPRISM].size();PP++) {
        ierr = doWork(
          mm,PP,MBENTITYSET,MBPRISM,
          row_data.dataOnEntities[MBENTITYSET][mm],col_data.dataOnEntities[MBPRISM][PP]
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

template<>
PetscErrorCode invertTensor3by3<3,double,ublas::row_major,ublas::unbounded_array<double> >(
  MatrixDouble &jac_data,
  VectorDouble &det_data,
  MatrixDouble &inv_jac_data
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  FTensor::Tensor2<double*,3,3> A = getTensor2FormData<3,3>(jac_data);
  int nb_gauss_pts = jac_data.size2();
  det_data.resize(nb_gauss_pts,false);
  inv_jac_data.resize(3,nb_gauss_pts,false);
  FTensor::Tensor0<double*> det = getTensor0FormData(det_data);
  FTensor::Tensor2<double*,3,3> I = getTensor2FormData<3,3>(inv_jac_data);
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    ierr = determinantTensor3by3(A,det); CHKERRQ(ierr);
    ierr = invertTensor3by3(A,det,I); CHKERRQ(ierr);
    ++A;
    ++det;
    ++I;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DataOperator::opRhs(
  DataForcesAndSurcesCore &data,
  const bool do_vertices,
  const bool do_edges,
  const bool do_quads,
  const bool do_tris,
  const bool do_tets,
  const bool do_prisms,
  const bool error_if_no_base
) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if(do_vertices) {
    for(unsigned int nn = 0;nn<data.dataOnEntities[MBVERTEX].size();nn++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBVERTEX][nn].getFieldData().size()&&
        (
          data.dataOnEntities[MBVERTEX][nn].getBase()==NOBASE||
          data.dataOnEntities[MBVERTEX][nn].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Vertex and side %d",nn);
      }
      ierr = doWork(nn,MBVERTEX,data.dataOnEntities[MBVERTEX][nn]); CHKERRQ(ierr);
    }
  }
  if(do_edges) {
    for(unsigned int ee = 0;ee<data.dataOnEntities[MBEDGE].size();ee++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBEDGE][ee].getFieldData().size()&&
        (
          data.dataOnEntities[MBEDGE][ee].getBase()==NOBASE||
          data.dataOnEntities[MBEDGE][ee].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Edge and side %d",ee);
      }
      ierr = doWork(ee,MBEDGE,data.dataOnEntities[MBEDGE][ee]); CHKERRQ(ierr);
    }
  }
  if(do_tris) {
    for(unsigned int ff = 0;ff<data.dataOnEntities[MBTRI].size();ff++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBTRI][ff].getFieldData().size()&&
        (
          data.dataOnEntities[MBTRI][ff].getBase()==NOBASE||
          data.dataOnEntities[MBTRI][ff].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Triangle and side %d",ff);
      }
      ierr = doWork(ff,MBTRI,data.dataOnEntities[MBTRI][ff]); CHKERRQ(ierr);
    }
  }
  if(do_quads) {
    for(unsigned int qq = 0;qq<data.dataOnEntities[MBQUAD].size();qq++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBQUAD][qq].getFieldData().size()&&
        (
          data.dataOnEntities[MBQUAD][qq].getBase()==NOBASE||
          data.dataOnEntities[MBQUAD][qq].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Quad and side %d",qq);
      }
      ierr = doWork(qq,MBQUAD,data.dataOnEntities[MBQUAD][qq]); CHKERRQ(ierr);
    }
  }
  if(do_tets) {
    for(unsigned int vv = 0;vv<data.dataOnEntities[MBTET].size();vv++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBTET][vv].getFieldData().size()&&
        (
          data.dataOnEntities[MBTET][vv].getBase()==NOBASE&&
          data.dataOnEntities[MBTET][vv].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Tet and side %d",vv);
      }
      ierr = doWork(vv,MBTET,data.dataOnEntities[MBTET][vv]); CHKERRQ(ierr);
    }
  }
  if(do_prisms) {
    for(unsigned int pp = 0;pp<data.dataOnEntities[MBPRISM].size();pp++) {
      if(
        error_if_no_base&&
        data.dataOnEntities[MBPRISM][pp].getFieldData().size()&&
        (
          data.dataOnEntities[MBPRISM][pp].getBase()==NOBASE||
          data.dataOnEntities[MBPRISM][pp].getBase()==LASTBASE
        )
      ) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No base on Prism and side %d",pp);
      }
      ierr = doWork(pp,MBPRISM,data.dataOnEntities[MBPRISM][pp]); CHKERRQ(ierr);
    }
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
    DataForcesAndSurcesCore::EntData &data
  ) {
  PetscFunctionBegin;
  // PetscErrorCode ierr;

  try {

    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];
      const unsigned int nb_base_functions = data.getN(base).size2();
      if(!nb_base_functions) continue;
      const unsigned int nb_gauss_pts = data.getN(base).size1();
      if(!nb_gauss_pts) continue;

      if(type!=MBVERTEX) {
        if(nb_base_functions != data.getDiffN(base).size2()/3) {
          SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
            "Data inconsistency nb_base_functions != data.diffN.size2()/3 ( %u != %u/3 )",
            nb_base_functions,data.getDiffN(base).size2()
          );
        }
      } else {
        if(
          data.getDiffN(base).size1()!=4||
          data.getDiffN(base).size2()!=3
        ) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data inconsistency");
        }
      }

      diffNinvJac.resize(
        data.getDiffN(base).size1(),data.getDiffN(base).size2(),false
      );

      double *t_diff_n_ptr = &*data.getDiffN(base).data().begin();
      FTensor::Tensor1<double*,3> t_diff_n(
        t_diff_n_ptr,&t_diff_n_ptr[1],&t_diff_n_ptr[2],3
      );
      double *t_inv_n_ptr = &*diffNinvJac.data().begin();
      FTensor::Tensor1<double*,3> t_inv_diff_n(
        t_inv_n_ptr,&t_inv_n_ptr[1],&t_inv_n_ptr[2],3
      );

      switch (type) {

        case MBVERTEX: {
          for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
            t_inv_diff_n(i) = t_diff_n(j)*tInvJac(j,i);
            ++t_diff_n;
            ++t_inv_diff_n;
          }
        }
        break;
        case MBEDGE:
        case MBTRI:
        case MBTET: {
          for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
            for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
              t_inv_diff_n(i) = t_diff_n(j)*tInvJac(j,i);
              ++t_diff_n;
              ++t_inv_diff_n;
            }
          }

        }
        break;
        default:
        SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }

      data.getDiffN(base).data().swap(diffNinvJac.data());


    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetInvJacHdivAndHcurl::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBEDGE && type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  if(
    (int)HDIV0_1!=(int)HCURL0_1 ||
    (int)HDIV0_2!=(int)HCURL0_2 ||
    (int)HDIV1_0!=(int)HCURL1_0 ||
    (int)HDIV1_2!=(int)HCURL1_2 ||
    (int)HDIV2_0!=(int)HCURL2_0 ||
    (int)HDIV2_1!=(int)HCURL2_1
  ) {
    SETERRQ(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "Data inconsistency between Hcurl and Hdiv struture of base functions"
    );
  }

  try {

    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];

      const unsigned int nb_gauss_pts = data.getDiffHdivN(base).size1();
      const unsigned int nb_base_functions = data.getDiffHdivN(base).size2()/9;
      if(!nb_base_functions) continue;

      diffHdivInvJac.resize(
        nb_gauss_pts,data.getDiffHdivN(base).size2(),false
      );

      FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHdivN<3,3>(base);
      double *t_inv_diff_n_ptr = &*diffHdivInvJac.data().begin();
      FTensor::Tensor2<double*,3,3> t_inv_diff_n(
        t_inv_diff_n_ptr,
        &t_inv_diff_n_ptr[HDIV0_1],
        &t_inv_diff_n_ptr[HDIV0_2],
        &t_inv_diff_n_ptr[HDIV1_0],
        &t_inv_diff_n_ptr[HDIV1_1],
        &t_inv_diff_n_ptr[HDIV1_2],
        &t_inv_diff_n_ptr[HDIV2_0],
        &t_inv_diff_n_ptr[HDIV2_1],
        &t_inv_diff_n_ptr[HDIV2_2],9
      );

      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          t_inv_diff_n(k,i) = t_diff_n(k,j)*tInvJac(j,i);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffHdivN(base).data().swap(diffHdivInvJac.data());

    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetContravariantPiolaTransform::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
)  {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {


    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];

      const unsigned int nb_base_functions = data.getHdivN(base).size2()/3;
      if(!nb_base_functions) continue;

      const double c = 1./6.;
      const unsigned int nb_gauss_pts = data.getHdivN(base).size1();
      piolaN.resize(nb_gauss_pts,data.getHdivN(base).size2(),false);
      piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN(base).size2(),false);

      FTensor::Tensor1<double*,3> t_n = data.getFTensor1HdivN<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<double*,3> t_transformed_n(
        t_transformed_n_ptr, //HDIV0
        &t_transformed_n_ptr[HDIV1],
        &t_transformed_n_ptr[HDIV2],3
      );
      FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHdivN<3,3>(base);
      double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
      FTensor::Tensor2<double*,3,3> t_transformed_diff_n(
        t_transformed_diff_n_ptr,
        &t_transformed_diff_n_ptr[HDIV0_1],
        &t_transformed_diff_n_ptr[HDIV0_2],
        &t_transformed_diff_n_ptr[HDIV1_0],
        &t_transformed_diff_n_ptr[HDIV1_1],
        &t_transformed_diff_n_ptr[HDIV1_2],
        &t_transformed_diff_n_ptr[HDIV2_0],
        &t_transformed_diff_n_ptr[HDIV2_1],
        &t_transformed_diff_n_ptr[HDIV2_2],9
      );

      double const a = c/vOlume;
      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          t_transformed_n(i) = a*tJac(i,k)*t_n(k);
          t_transformed_diff_n(i,k) = a*tJac(i,j)*t_diff_n(j,k);
          ++t_n;
          ++t_transformed_n;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
      }
      data.getHdivN(base).data().swap(piolaN.data());
      data.getDiffHdivN(base).data().swap(piolaDiffN.data());

    }

    // data.getBase() = base;

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetCovariantPiolaTransform::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
)  {
  PetscFunctionBegin;

  if(type != MBEDGE && type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];

      const unsigned int nb_base_functions = data.getHcurlN(base).size2()/3;
      if(!nb_base_functions) continue;

      const unsigned int nb_gauss_pts = data.getHcurlN(base).size1();
      piolaN.resize(nb_gauss_pts,data.getHcurlN(base).size2(),false);
      piolaDiffN.resize(nb_gauss_pts,data.getDiffHcurlN(base).size2(),false);

      FTensor::Tensor1<double*,3> t_n = data.getFTensor1HcurlN<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<double*,3> t_transformed_n(
        t_transformed_n_ptr, //HDIV0
        &t_transformed_n_ptr[HCURL1],
        &t_transformed_n_ptr[HCURL2],3
      );
      FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHdivN<3,3>(base);
      double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
      FTensor::Tensor2<double*,3,3> t_transformed_diff_n(
        t_transformed_diff_n_ptr,
        &t_transformed_diff_n_ptr[HCURL0_1],
        &t_transformed_diff_n_ptr[HCURL0_2],
        &t_transformed_diff_n_ptr[HCURL1_0],
        &t_transformed_diff_n_ptr[HCURL1_1],
        &t_transformed_diff_n_ptr[HCURL1_2],
        &t_transformed_diff_n_ptr[HCURL2_0],
        &t_transformed_diff_n_ptr[HCURL2_1],
        &t_transformed_diff_n_ptr[HCURL2_2],9
      );

      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          t_transformed_n(i) = tInvJac(k,i)*t_n(k);
          t_transformed_diff_n(i,k) = tInvJac(j,i)*t_diff_n(j,k);
          ++t_n;
          ++t_transformed_n;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
      }
      data.getHcurlN(base).data().swap(piolaN.data());
      data.getDiffHcurlN(base).data().swap(piolaDiffN.data());

    }

    // data.getBase() = base;

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetHoInvJacH1::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  )  {
    PetscFunctionBegin;

    try {

      for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

        FieldApproximationBase base = ApproximationBaseArray[b];
        if(data.getDiffN(base).size2()==0) continue;

        unsigned int nb_gauss_pts = data.getN(base).size1();
        if(nb_gauss_pts==0) continue;
        unsigned int nb_base_functions = data.getN(base).size2();
        if(nb_base_functions==0) continue;

        // Note for Vetex diffN row has size of number of dof
        diffNinvJac.resize(nb_gauss_pts,3*nb_base_functions,false);
        double *t_inv_n_ptr = &*diffNinvJac.data().begin();
        FTensor::Tensor1<double*,3> t_inv_diff_n(
          t_inv_n_ptr,&t_inv_n_ptr[1],&t_inv_n_ptr[2],3
        );
        if(invHoJac.size2()!=9) {
          SETERRQ1(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "It looks that ho inverse of jacobian is not calculated %d != 9",
            invHoJac.size2()
          );
        }
        if(invHoJac.size1()!=nb_gauss_pts) {
          SETERRQ2(
            PETSC_COMM_SELF,
            MOFEM_DATA_INCONSISTENCY,
            "It looks that ho inverse of jacobian is not calculated %d != %d",
            invHoJac.size1(),
            nb_gauss_pts
          );
        }
        double *t_inv_jac_ptr = &*invHoJac.data().begin();
        FTensor::Tensor2<double*,3,3> t_inv_jac(
          t_inv_jac_ptr,&t_inv_jac_ptr[1],&t_inv_jac_ptr[2],
          &t_inv_jac_ptr[3],&t_inv_jac_ptr[4],&t_inv_jac_ptr[5],
          &t_inv_jac_ptr[6],&t_inv_jac_ptr[7],&t_inv_jac_ptr[8],9
        );

        switch (type) {
          case MBVERTEX: {
            // std::cerr << data.getDiffN(base) << std::endl;
            double *t_diff_n_ptr = &*data.getDiffN(base).data().begin();
            for(int gg = 0;gg!=nb_gauss_pts;gg++) {
              FTensor::Tensor1<double*,3> t_diff_n(
                t_diff_n_ptr,&t_diff_n_ptr[1],&t_diff_n_ptr[2],3
              );
              for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
                t_inv_diff_n(i) = t_diff_n(j)*t_inv_jac(j,i);
                ++t_diff_n;
                ++t_inv_diff_n;
              }
              ++t_inv_jac;
            }
          }
          break;
          case MBEDGE:
          case MBTRI:
          case MBTET: {
            FTensor::Tensor1<double*,3> t_diff_n = data.getFTensor1DiffN<3>(base);
            for(unsigned int gg = 0;gg<nb_gauss_pts;gg++) {
              for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
                t_inv_diff_n(i) = t_diff_n(j)*t_inv_jac(j,i);
                ++t_diff_n;
                ++t_inv_diff_n;
              }
              ++t_inv_jac;
            }
          }
          break;
          default:
          SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
        }

        // unsigned int gg = 0;
        // for(;gg<nb_gauss_pts;gg++) {
        //   double *inv_h = &invHoJac(gg,0);
        //   for(unsigned dd = 0;dd<nb_base_functions;dd++) {
        //     double *diff_n;
        //     if(type == MBVERTEX) {
        //       diff_n = &data.getDiffN(base)(dd,0);
        //     } else {
        //       diff_n = &data.getDiffN(base)(gg,3*dd);
        //     }
        //     double *diff_n_inv_jac = &diffNinvJac(gg,3*dd);
        //     cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_h,3,diff_n,1,0.,diff_n_inv_jac,1);
        //   }
        // }
        if(type == MBVERTEX) {
          data.getDiffN(base).resize(diffNinvJac.size1(),diffNinvJac.size2(),false);
        }
        data.getDiffN(base).data().swap(diffNinvJac.data());

      }

    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode OpSetHoInvJacHdivAndHcurl::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
  ) {
    PetscFunctionBegin;

    if(type != MBEDGE && type != MBTRI && type != MBTET) PetscFunctionReturn(0);
    // if(data.getSpace() == HDIV && type == MBEDGE) PetscFunctionReturn(0);

    try {

      for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

        FieldApproximationBase base = ApproximationBaseArray[b];

        diffHdivInvJac.resize(data.getDiffHdivN(base).size1(),data.getDiffHdivN(base).size2(),false);

        unsigned int nb_gauss_pts = data.getDiffHdivN(base).size1();
        unsigned int nb_base_functions = data.getDiffHdivN(base).size2()/9;
        if(nb_base_functions == 0) continue;

        FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHdivN<3,3>(base);
        double *t_inv_diff_n_ptr = &*diffHdivInvJac.data().begin();
        FTensor::Tensor2<double*,3,3> t_inv_diff_n(
          t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HDIV0_1], &t_inv_diff_n_ptr[HDIV0_2],
          &t_inv_diff_n_ptr[HDIV1_0], &t_inv_diff_n_ptr[HDIV1_1], &t_inv_diff_n_ptr[HDIV1_2],
          &t_inv_diff_n_ptr[HDIV2_0], &t_inv_diff_n_ptr[HDIV2_1], &t_inv_diff_n_ptr[HDIV2_2],9
        );
        double *t_inv_jac_ptr = invHoJac.data().begin();
        // cerr << invHoJac << endl;
        FTensor::Tensor2<double*,3,3> t_inv_jac(
          t_inv_jac_ptr,&t_inv_jac_ptr[1],&t_inv_jac_ptr[2],
          &t_inv_jac_ptr[3],&t_inv_jac_ptr[4],&t_inv_jac_ptr[5],
          &t_inv_jac_ptr[6],&t_inv_jac_ptr[7],&t_inv_jac_ptr[8],9
        );

        for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
          for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
            t_inv_diff_n(i,j) = t_inv_jac(k,j)*t_diff_n(i,k);
            ++t_diff_n;
            ++t_inv_diff_n;
          }
          ++t_inv_jac;
        }

        // unsigned int gg = 0;
        // for(;gg<nb_gauss_pts;gg++) {
        //   double *inv_h = &invHoJac(gg,0);
        //   for(unsigned dd = 0;dd<nb_base_functions;dd++) {
        //     const double *diff_hdiv = &(data.getDiffHdivN(base,gg)(dd,0));
        //     double *diff_hdiv_inv_jac = &diffHdivInvJac(gg,9*dd);
        //     int kk = 0;
        //     for(;kk<3;kk++) {
        //       cblas_dgemv(CblasRowMajor,CblasTrans,3,3,1.,inv_h,3,&diff_hdiv[kk],3,0.,&diff_hdiv_inv_jac[kk],3);
        //     }
        //   }
        // }

        data.getDiffHdivN(base).data().swap(diffHdivInvJac.data());

      }

    } catch (std::exception& ex) {
      std::ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }


PetscErrorCode OpSetHoContravariantPiolaTransform::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];

      unsigned int nb_gauss_pts = data.getHdivN(base).size1();
      unsigned int nb_base_functions = data.getHdivN(base).size2()/3;
      piolaN.resize(nb_gauss_pts,data.getHdivN(base).size2(),false);
      piolaDiffN.resize(nb_gauss_pts,data.getDiffHdivN(base).size2(),false);

      FTensor::Tensor1<double*,3> t_n = data.getFTensor1HdivN<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<double*,3> t_transformed_n(
        t_transformed_n_ptr, //HDIV0
        &t_transformed_n_ptr[HDIV1],
        &t_transformed_n_ptr[HDIV2],3
      );
      FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHdivN<3,3>(base);
      double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
      FTensor::Tensor2<double*,3,3> t_transformed_diff_n(
        t_transformed_diff_n_ptr,
        &t_transformed_diff_n_ptr[HDIV0_1],
        &t_transformed_diff_n_ptr[HDIV0_2],
        &t_transformed_diff_n_ptr[HDIV1_0],
        &t_transformed_diff_n_ptr[HDIV1_1],
        &t_transformed_diff_n_ptr[HDIV1_2],
        &t_transformed_diff_n_ptr[HDIV2_0],
        &t_transformed_diff_n_ptr[HDIV2_1],
        &t_transformed_diff_n_ptr[HDIV2_2],9
      );

      FTensor::Tensor0<double*> t_det(&*detHoJac.data().begin());
      double *t_jac_ptr = hoJac.data().begin();
      FTensor::Tensor2<double*,3,3> t_jac(
        t_jac_ptr,&t_jac_ptr[1],&t_jac_ptr[2],
        &t_jac_ptr[3],&t_jac_ptr[4],&t_jac_ptr[5],
        &t_jac_ptr[6],&t_jac_ptr[7],&t_jac_ptr[8],9
      );

      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          const double a = 1./t_det;
          t_transformed_n(i) = a*t_jac(i,k)*t_n(k);
          t_transformed_diff_n(i,k) = a*t_jac(i,j)*t_diff_n(j,k);
          ++t_n;
          ++t_transformed_n;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
        ++t_det;
        ++t_jac;
      }

      data.getHdivN(base).data().swap(piolaN.data());
      data.getDiffHdivN(base).data().swap(piolaDiffN.data());

    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetHoCovariantPiolaTransform::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBEDGE && type != MBTRI && type != MBTET) PetscFunctionReturn(0);

  try {

    for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

      FieldApproximationBase base = ApproximationBaseArray[b];

      unsigned int nb_gauss_pts = data.getHcurlN(base).size1();
      unsigned int nb_base_functions = data.getHcurlN(base).size2()/3;
      piolaN.resize(nb_gauss_pts,data.getHcurlN(base).size2(),false);
      piolaDiffN.resize(nb_gauss_pts,data.getDiffHcurlN(base).size2(),false);

      FTensor::Tensor1<double*,3> t_n = data.getFTensor1HcurlN<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<double*,3> t_transformed_n(
        t_transformed_n_ptr, //HDIV0
        &t_transformed_n_ptr[HCURL1],
        &t_transformed_n_ptr[HCURL2],3
      );
      FTensor::Tensor2<double*,3,3> t_diff_n = data.getFTensor2DiffHcurlN<3,3>(base);
      double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
      FTensor::Tensor2<double*,3,3> t_transformed_diff_n(
        t_transformed_diff_n_ptr,
        &t_transformed_diff_n_ptr[HCURL0_1],
        &t_transformed_diff_n_ptr[HCURL0_2],
        &t_transformed_diff_n_ptr[HCURL1_0],
        &t_transformed_diff_n_ptr[HCURL1_1],
        &t_transformed_diff_n_ptr[HCURL1_2],
        &t_transformed_diff_n_ptr[HCURL2_0],
        &t_transformed_diff_n_ptr[HCURL2_1],
        &t_transformed_diff_n_ptr[HCURL2_2],9
      );

      double *t_inv_jac_ptr = hoInvJac.data().begin();
      FTensor::Tensor2<double*,3,3> t_inv_jac(
        t_inv_jac_ptr,&t_inv_jac_ptr[1],&t_inv_jac_ptr[2],
        &t_inv_jac_ptr[3],&t_inv_jac_ptr[4],&t_inv_jac_ptr[5],
        &t_inv_jac_ptr[6],&t_inv_jac_ptr[7],&t_inv_jac_ptr[8],9
      );

      for(unsigned int gg = 0;gg!=nb_gauss_pts;gg++) {
        for(unsigned int bb = 0;bb!=nb_base_functions;bb++) {
          t_transformed_n(i) = t_inv_jac(k,i)*t_n(k);
          t_transformed_diff_n(i,k) = t_inv_jac(j,i)*t_diff_n(j,k);
          ++t_n;
          ++t_transformed_n;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
        ++t_inv_jac;
      }

      data.getHcurlN(base).data().swap(piolaN.data());
      data.getDiffHcurlN(base).data().swap(piolaDiffN.data());

    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
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
    if(nb_dofs % rAnk != 0) {
      SETERRQ4(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
        "data inconsistency, type %d, side %d, nb_dofs %d, rAnk %d",
        type,side,nb_dofs,rAnk
      );
    }
    if(nb_dofs/rAnk > data.getN().size2()) {
      std::cerr << side << " " << type << " " << ApproximationBaseNames[data.getBase()] << std::endl;
      std::cerr << data.getN() << std::endl;
      std::cerr << data.getN(NOBASE) << std::endl;
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
        "data inconsistency nb_dofs >= data.N.size2() %u >= %u",nb_dofs,data.getN().size2()
      );
    }

    if(type == MBVERTEX) {
      data_at_GaussPt.resize(data.getN().size1(),rAnk,false);
      dataGrad_at_GaussPt.resize(data.getN().size1(),rAnk*dIm,false);
      data_at_GaussPt.clear();
      dataGrad_at_GaussPt.clear();
      // bzero(&*data_at_GaussPt.data().begin(),data.getN().size1()*rAnk*sizeof(FieldData));
      // bzero(&*dataGrad_at_GaussPt.data().begin(),data.getN().size1()*rAnk*dIm*sizeof(FieldData));
      for(int rr = 0;rr<rAnk;rr++) {
        for(unsigned int dd = 0;dd<dIm;dd++) {
          dataGrad_at_GaussPt(0,dIm*rr+dd) = cblas_ddot(nb_dofs/rAnk,&data.getDiffN()(0,dd),dIm,&data.getFieldData()[rr],rAnk);
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

      for(int rr = 0;rr<rAnk;rr++,data_ptr++) {
        data_at_GaussPt(gg,rr) += cblas_ddot(nb_dofs/rAnk,n_ptr,1,data_ptr,rAnk);
        double *diff_n_ptr2 = diff_n_ptr;

        for(unsigned int dd = 0;dd<dIm;dd++,diff_n_ptr2++) {
          if(type == MBVERTEX) {
            if(gg == 0) continue;
            dataGrad_at_GaussPt(gg,dIm*rr+dd) += dataGrad_at_GaussPt(0,dIm*rr+dd);
          } else {
            dataGrad_at_GaussPt(gg,dIm*rr+dd) += cblas_ddot(nb_dofs/rAnk,diff_n_ptr2,dIm,data_ptr,rAnk);
          }
        }
      }
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
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

  // std::cerr << type << " " << side << " " << ApproximationBaseNames[data.getBase()] << std::endl;

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
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(nb_dofs%3!=0) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
      }
      if(nb_dofs > 3*data.getN().size2()) {
        SETERRQ1(
          PETSC_COMM_SELF,
          MOFEM_DATA_INCONSISTENCY,
          "data inconsistency for base %s",
          ApproximationBaseNames[data.getBase()]
        );
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

    if(type == MBEDGE) {
      if(!valid_edges3[side]||valid_edges4[side]) PetscFunctionReturn(0);
    } else if(type == MBTRI) {
      if(!valid_faces3[side]||valid_faces4[side]) PetscFunctionReturn(0);
    }

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
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        unsigned int nb_dofs = data.getFieldData().size();
        if(nb_dofs%3!=0) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
        }
        if(nb_dofs > 3*data.getN().size2()) {
          SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, side %d type %d",side,type);
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

  } catch (std::exception& ex) {
    std::ostringstream ss;
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
  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OpSetContravariantPiolaTransoformOnTriangle::doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBTRI) PetscFunctionReturn(0);

  // FieldApproximationBase base = data.getBase();

  for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];
    // data.getBase() = ApproximationBaseArray[b];

    double l0 = cblas_dnrm2(3,&nOrmal[0],1);
    int nb_gauss_pts = data.getHdivN(base).size1();
    int nb_dofs = data.getHdivN(base).size2()/3;
    int gg = 0;
    for(;gg<nb_gauss_pts;gg++) {
      int dd = 0;
      for(;dd<nb_dofs;dd++) {
        double val = data.getHdivN(base)(gg,3*dd);
        if(normalsAtGaussPt.size1()==(unsigned int)nb_gauss_pts) {
          double l = cblas_dnrm2(3,&normalsAtGaussPt(gg,0),1);
          cblas_dcopy(3,&normalsAtGaussPt(gg,0),1,&data.getHdivN(base)(gg,3*dd),1);
          cblas_dscal(3,val/pow(l,2),&data.getHdivN(base)(gg,3*dd),1);
        } else {
          cblas_dcopy(3,&nOrmal[0],1,&data.getHdivN(base)(gg,3*dd),1);
          cblas_dscal(3,val/pow(l0,2),&data.getHdivN(base)(gg,3*dd),1);
        }
      }
    }
  }

  // data.getBase() = base;

  PetscFunctionReturn(0);
}

PetscErrorCode OpSetCovariantPiolaTransoformOnTriangle::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if(type != MBEDGE && type != MBTRI) PetscFunctionReturn(0);

  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  FTensor::Index<'k',2> k;

  double zero = 0;
  FTensor::Tensor2<const double*,3,3> t_m(
    &tAngent0[0],&tAngent1[0],&nOrmal[0],
    &tAngent0[1],&tAngent1[1],&nOrmal[1],
    &tAngent0[2],&tAngent1[2],&nOrmal[2],3
  );
  double det;
  FTensor::Tensor2<double,3,3> t_inv_m;
  ierr = determinantTensor3by3(t_m,det); CHKERRQ(ierr);
  ierr = invertTensor3by3(t_m,det,t_inv_m); CHKERRQ(ierr);

  for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    int nb_dofs = data.getHcurlN(base).size2()/3;
    int nb_gauss_pts = data.getHcurlN(base).size1();

    MatrixDouble piola_n(data.getHcurlN(base).size1(),data.getHcurlN(base).size2());
    MatrixDouble diff_piola_n(data.getDiffHcurlN(base).size1(),data.getDiffHcurlN(base).size2());

    if(nb_dofs>0 && nb_gauss_pts>0) {

      FieldApproximationBase base = ApproximationBaseArray[b];
      FTensor::Tensor1<double*,3> t_h_curl(
        &data.getHcurlN(base)(0,HCURL0),
        &data.getHcurlN(base)(0,HCURL1),
        &data.getHcurlN(base)(0,HCURL2),3
      );
      FTensor::Tensor2<double*,3,2> t_diff_h_curl(
        &data.getDiffHcurlN(base)(0,HCURL0_0),&data.getDiffHcurlN(base)(0,HCURL0_1),
        &data.getDiffHcurlN(base)(0,HCURL1_0),&data.getDiffHcurlN(base)(0,HCURL1_1),
        &data.getDiffHcurlN(base)(0,HCURL2_0),&data.getDiffHcurlN(base)(0,HCURL2_1),6
      );
      FTensor::Tensor1<double*,3> t_transformed_h_curl(
        &piola_n(0,HCURL0),
        &piola_n(0,HCURL1),
        &piola_n(0,HCURL2),3
      );
      FTensor::Tensor2<double*,3,2> t_transformed_diff_h_curl(
        &diff_piola_n(0,HCURL0_0),&diff_piola_n(0,HCURL0_1),
        &diff_piola_n(0,HCURL1_0),&diff_piola_n(0,HCURL1_1),
        &diff_piola_n(0,HCURL2_0),&diff_piola_n(0,HCURL2_1),6
      );

      int cc = 0;
      if(normalsAtGaussPt.size1()==(unsigned int)nb_gauss_pts) {
        // HO geometry is set, so jacobian is different at each gauss point
        FTensor::Tensor2<const double*,3,3> t_m_at_pts(
          &tangent0AtGaussPt(0,0),&tangent1AtGaussPt(0,0),&normalsAtGaussPt(0,0),
          &tangent0AtGaussPt(0,1),&tangent1AtGaussPt(0,1),&normalsAtGaussPt(0,1),
          &tangent0AtGaussPt(0,2),&tangent1AtGaussPt(0,2),&normalsAtGaussPt(0,2),3
        );
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = determinantTensor3by3(t_m_at_pts,det); CHKERRQ(ierr);
          ierr = invertTensor3by3(t_m_at_pts,det,t_inv_m); CHKERRQ(ierr);
          for(int ll = 0;ll!=nb_dofs;ll++) {
            t_transformed_h_curl(i) = t_inv_m(j,i)*t_h_curl(j);
            t_transformed_diff_h_curl(i,k) = t_inv_m(j,i)*t_diff_h_curl(j,k);
            ++t_h_curl;
            ++t_transformed_h_curl;
            ++t_diff_h_curl;
            ++t_transformed_diff_h_curl;
            ++cc;
          }
          ++t_m_at_pts;
        }
      } else {
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          for(int ll = 0;ll!=nb_dofs;ll++) {
            t_transformed_h_curl(i) = t_inv_m(j,i)*t_h_curl(j);
            t_transformed_diff_h_curl(i,k) = t_inv_m(j,i)*t_diff_h_curl(j,k);
            ++t_h_curl;
            ++t_transformed_h_curl;
            ++t_diff_h_curl;
            ++t_transformed_diff_h_curl;
            ++cc;
          }
        }
      }
      if(cc!=nb_gauss_pts*nb_dofs) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"Data inconsistency");
      }
      data.getHcurlN(base).data().swap(piola_n.data());
      data.getDiffHcurlN(base).data().swap(diff_piola_n.data());

    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode OpGetHoTangentOnEdge::doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if(nb_dofs == 0) PetscFunctionReturn(0);

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
          SETERRQ(
            PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"Approximated field should be rank 3, i.e. vector in 3d space"
          );
        }
        for(int dd = 0;dd!=3;dd++) {
          for(int gg = 0;gg!=nb_gauss_pts;gg++) {
            tAngent(gg,dd) += cblas_ddot(nb_dofs/3,&diff[gg*nb_approx_fun],1,dofs[dd],3);
          }
        }
      break;
      default:
      SETERRQ(
        PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"This operator can calculate tangent vector only on edge"
      );
    }


  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode OpSetCovariantPiolaTransoformOnEdge::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(type != MBEDGE) PetscFunctionReturn(0);

  FTensor::Index<'i',3> i;
  FTensor::Tensor1<const double*,3> t_m(
    &tAngent[0],&tAngent[1],&tAngent[2]
  );
  const double l0 = t_m(i)*t_m(i);
  std::vector<double> l1;
  {
    int nb_gauss_pts = tangentAtGaussPt.size1();
    if(nb_gauss_pts) {
      l1.resize(nb_gauss_pts);
      FTensor::Tensor1<const double*,3> t_m_at_pts(
        &tangentAtGaussPt(0,0),&tangentAtGaussPt(0,1),&tangentAtGaussPt(0,2),3
      );
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        l1[gg] = t_m_at_pts(i)*t_m_at_pts(i);
        ++t_m_at_pts;
      }
    }
  }

  for(int b = AINSWORTH_COLE_BASE; b!=USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];
    int nb_gauss_pts = data.getHcurlN(base).size1();
    int nb_dofs = data.getHcurlN(base).size2()/3;
    if(nb_gauss_pts>0 && nb_dofs>0) {
      FTensor::Tensor1<double*,3> t_h_curl(
        &data.getHcurlN(base)(0,HCURL0),
        &data.getHcurlN(base)(0,HCURL1),
        &data.getHcurlN(base)(0,HCURL2),3
      );
      int cc = 0;
      if(tangentAtGaussPt.size1()==(unsigned int)nb_gauss_pts) {
        FTensor::Tensor1<const double*,3> t_m_at_pts(
          &tangentAtGaussPt(0,0),&tangentAtGaussPt(0,1),&tangentAtGaussPt(0,2),3
        );
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          const double l0 = l1[gg];
          for(int ll = 0;ll!=nb_dofs;ll++) {
            const double val = t_h_curl(0);
            const double a = val/l0;
            t_h_curl(i) = t_m_at_pts(i)*a;
            ++t_h_curl;
            ++cc;
          }
          ++t_m_at_pts;
        }
      } else {
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          for(int ll = 0;ll!=nb_dofs;ll++) {
            const double val = t_h_curl(0);
            const double a = val/l0;
            t_h_curl(i) = t_m(i)*a;
            ++t_h_curl;
            ++cc;
          }
        }
      }
      if(cc!=nb_gauss_pts*nb_dofs) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"Data inconsistency");
      }
    }

  }

  PetscFunctionReturn(0);
}


}
