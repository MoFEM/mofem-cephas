/** \file HexPolynomialBase.cpp
\brief Implementation of hierarchical bases on tetrahedral

A l2, h1, h-div and h-curl spaces are implemented.

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
#include <config.h>
#include <definitions.h>
#include <version.h>

#include <Common.hpp>
#include <UnknownInterface.hpp>
#include <base_functions.h>
#include <fem_tools.h>
#include <h1_hdiv_hcurl_l2.h>
using namespace MoFEM;

#include <AdjacencyMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <DataStructures.hpp>
#include <DofsMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <FTensor.hpp>
#include <FieldMultiIndices.hpp>
#include <LoopMethods.hpp>
#include <ProblemsMultiIndices.hpp>
#include <TagMultiIndices.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <HexPolynomialBase.hpp>

#include <Hcurl.hpp>
#include <Hdiv.hpp>

MoFEMErrorCode
HexPolynomialBase::query_interface(const MOFEMuuid &uuid,
                                   MoFEM::UnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_TET_BASE_FUNCTION) {
    *iface = const_cast<HexPolynomialBase *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = BaseFunction::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

HexPolynomialBase::HexPolynomialBase() {}
HexPolynomialBase::~HexPolynomialBase() {}

MoFEMErrorCode HexPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }

   PetscErrorCode (*base_polynomials)(int p, double alpha, double x,
                                            double t, double *diff_x,
                                            double *diff_t, double *L,
                                            double *diffL, const int dim) =
      cTx->basePolynomialsType1;

  int nb_gauss_pts = pts.size2();

  int sense[12], order[12];
  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    // edges
    if (data.dataOnEntities[MBEDGE].size() != 12) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *h1_edge_n[12], *diff_h1_edge_n[12];
    for (int ee = 0; ee != 12; ++ee) {
      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                        false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                            3 * nb_dofs, false);
      h1_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_h1_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR H1_EdgeShapeFunctions_MBHEX(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        h1_edge_n, diff_h1_edge_n, nb_gauss_pts, base_polynomials);
  } else {
    for (int ee = 0; ee != 12; ++ee) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(0, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(0, 0, false);
    }
  }

  if (data.spacesOnEntities[MBQUAD].test(H1)) {
    // faces
    if (data.dataOnEntities[MBQUAD].size() != 6) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *h1_face_n[6], *diff_h1_face_n[6];
    for (int ff = 0; ff != 6; ++ff) {
      if (data.dataOnEntities[MBQUAD][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      int nb_dofs = NBFACEQUAD_H1(data.dataOnEntities[MBQUAD][ff].getDataOrder());
      order[ff] = data.dataOnEntities[MBQUAD][ff].getDataOrder();
      data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                       false);
      data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(nb_gauss_pts,
                                                           3 * nb_dofs, false);
      h1_face_n[ff] =
          &*data.dataOnEntities[MBQUAD][ff].getN(base).data().begin();
      diff_h1_face_n[ff] =
          &*data.dataOnEntities[MBQUAD][ff].getDiffN(base).data().begin();
    }
    if (data.facesNodes.size1() != 8) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    if (data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    // CHKERR H1_FaceShapeFunctions_MBHEX(
    //     &*data.facesNodes.data().begin(), order,
    //     &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
    //     &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
    //     h1_face_n, diff_h1_face_n, nb_gauss_pts, base_polynomials);

  } else {
    for (int ff = 0; ff != 6; ++ff) {
      data.dataOnEntities[MBQUAD][ff].getN(base).resize(0, false);
      data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(0, 0, false);
    }
  }

  // if (data.spacesOnEntities[MBHEX].test(H1)) {
  //   // volume
  //   int order = data.dataOnEntities[MBHEX][0].getDataOrder();
  //   int nb_vol_dofs = NBVOLUMETET_H1(order);
  //   data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, nb_vol_dofs,
  //                                                   false);
  //   data.dataOnEntities[MBHEX][0].getDiffN(base).resize(nb_gauss_pts,
  //                                                       3 * nb_vol_dofs, false);
  //   CHKERR H1_VolumeShapeFunctions_MBHEX(
  //       data.dataOnEntities[MBHEX][0].getDataOrder(),
  //       &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
  //       &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
  //       &*data.dataOnEntities[MBHEX][0].getN(base).data().begin(),
  //       &*data.dataOnEntities[MBHEX][0].getDiffN(base).data().begin(),
  //       nb_gauss_pts, base_polynomials);
  // } else {
  //   data.dataOnEntities[MBHEX][0].getN(base).resize(0, 0, false);
  //   data.dataOnEntities[MBHEX][0].getDiffN(base).resize(0, 0, false);
  // }

  MoFEMFunctionReturn(0);
}



MoFEMErrorCode
HexPolynomialBase::getValue(MatrixDouble &pts,
                            boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;

  MoFEM::UnknownInterface *iface;
  CHKERR ctx_ptr->query_interface(IDD_TET_BASE_FUNCTION, &iface);
  cTx = reinterpret_cast<EntPolynomialBaseCtx *>(iface);

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts) {
    MoFEMFunctionReturnHot(0);
  }

  if (pts.size1() < 3) {
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");
  }

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSourcesCore &data = cTx->dAta;
  if (cTx->copyNodeBase == LASTBASE) {
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 8, false);
    CHKERR ShapeMBHEX(
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &pts(0, 0), &pts(1, 0), &pts(2, 0), nb_gauss_pts);
  } else {
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
        data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
  }
  if (data.dataOnEntities[MBVERTEX][0].getN(base).size1() !=
      (unsigned int)nb_gauss_pts) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[base]);
  }
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(8, 3, false);
  CHKERR ShapeDiffMBHEX(
     &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &pts(0, 0), &pts(1, 0), &pts(2, 0), nb_gauss_pts);

  switch (cTx->sPace) {
  case H1:
    CHKERR getValueH1(pts);
    break;
  // case HDIV:
  //   CHKERR getValueHdiv(pts);
  //   break;
  // case HCURL:
  //   CHKERR getValueHcurl(pts);
  //   break;
  // case L2:
  //   CHKERR getValueL2(pts);
  //   break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown space");
  }

  MoFEMFunctionReturn(0);
}
