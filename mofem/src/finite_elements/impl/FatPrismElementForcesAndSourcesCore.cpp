/** \file FatPrismElementForcesAndSourcesCore.cpp

\brief Implementation of fat prism element

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

#include <base_functions.h>
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
#include <FEMultiIndices.hpp>
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
#include <ForcesAndSourcesCore.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FatPrismElementForcesAndSourcesCore.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <FlatPrismPolynomialBase.hpp>
#include <FatPrismPolynomialBase.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
  #include <quad.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

FatPrismElementForcesAndSourcesCore::FatPrismElementForcesAndSourcesCore(
    Interface &m_field)
    : VolumeElementForcesAndSourcesCore(m_field, MBPRISM),
      dataH1TrianglesOnly(MBPRISM), dataH1TroughThickness(MBPRISM),
      opHOCoordsAndNormals(hoCoordsAtGaussPtsF3, nOrmals_at_GaussPtF3,
                           tAngent1_at_GaussPtF3, tAngent2_at_GaussPtF3,
                           hoCoordsAtGaussPtsF4, nOrmals_at_GaussPtF4,
                           tAngent1_at_GaussPtF4, tAngent2_at_GaussPtF4) {}

MoFEMErrorCode FatPrismElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBPRISM)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  int num_nodes;
  const EntityHandle *conn;
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  {
    coords.resize(num_nodes * 3, false);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());
    double diff_n[6];
    CHKERR ShapeDiffMBTRI(diff_n);
    normal.resize(6, false);
    CHKERR ShapeFaceNormalMBTRI(diff_n, &coords[0], &normal[0]);
    CHKERR ShapeFaceNormalMBTRI(diff_n, &coords[9], &normal[3]);
    aRea[0] = cblas_dnrm2(3, &normal[0], 1) * 0.5;
    aRea[1] = cblas_dnrm2(3, &normal[3], 1) * 0.5;
  }

  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getSpacesAndBaseOnEntities(dataH1TrianglesOnly);
  CHKERR getSpacesAndBaseOnEntities(dataH1TroughThickness);

  // H1
  if ((dataH1.spacesOnEntities[MBEDGE]).test(H1)) {
    CHKERR getEntitySense<MBEDGE>(dataH1);
    CHKERR getEntitySense<MBTRI>(dataH1);
    CHKERR getEntitySense<MBQUAD>(dataH1);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
    CHKERR getEntityDataOrder<MBTRI>(dataH1, H1);
    CHKERR getEntityDataOrder<MBQUAD>(dataH1, H1);
    CHKERR getEntityDataOrder<MBPRISM>(dataH1, H1);
    // Triangles only
    CHKERR getEntitySense<MBEDGE>(dataH1TrianglesOnly);
    CHKERR getEntitySense<MBTRI>(dataH1TrianglesOnly);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1TrianglesOnly, H1);
    CHKERR getEntityDataOrder<MBTRI>(dataH1TrianglesOnly, H1);
    // Through thickness
    CHKERR getEntitySense<MBEDGE>(dataH1TroughThickness);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1TroughThickness, H1);
  }
  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");
  }
  // Hcurl
  if ((dataH1.spacesOnEntities[MBEDGE]).test(HCURL)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");
  }
  // L2
  if ((dataH1.spacesOnEntities[MBPRISM]).test(L2)) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");
  }

  // get approx. on triangles, i.e. faces 3 and 4
  int nb_gauss_pts_on_faces;
  try {
    int order_triangles_only = 1;
    int valid_edges[] = {1, 1, 1, 0, 0, 0, 1, 1, 1};
    for (unsigned int ee = 0; ee < 9; ee++) {
      if (!valid_edges[ee])
        continue;
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getDataOrder());
    }
    for (unsigned int ff = 3; ff <= 4; ff++) {
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getDataOrder());
    }
    for (unsigned int qq = 0; qq < 3; qq++) {
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getDataOrder());
    }
    order_triangles_only = std::max(
        order_triangles_only,
        dataH1TroughThickness.dataOnEntities[MBPRISM][0].getDataOrder());
    // integration pts on the triangles surfaces
    int rule = getRuleTrianglesOnly(order_triangles_only);
    if (rule >= 0) {
      if (rule < QUAD_2D_TABLE_SIZE) {
        if (QUAD_2D_TABLE[rule]->dim != 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
        }
        if (QUAD_2D_TABLE[rule]->order < rule) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "wrong order %d != %d", QUAD_2D_TABLE[rule]->order, rule);
        }
        nb_gauss_pts_on_faces = QUAD_2D_TABLE[rule]->npoints;
        gaussPtsTrianglesOnly.resize(3, nb_gauss_pts_on_faces, false);
        cblas_dcopy(nb_gauss_pts_on_faces, &QUAD_2D_TABLE[rule]->points[1], 3,
                    &gaussPtsTrianglesOnly(0, 0), 1);
        cblas_dcopy(nb_gauss_pts_on_faces, &QUAD_2D_TABLE[rule]->points[2], 3,
                    &gaussPtsTrianglesOnly(1, 0), 1);
        cblas_dcopy(nb_gauss_pts_on_faces, QUAD_2D_TABLE[rule]->weights, 1,
                    &gaussPtsTrianglesOnly(2, 0), 1);
        dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
            nb_gauss_pts_on_faces, 3, false);
        double *shape_ptr = &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0]
                                  .getN(NOBASE)
                                  .data()
                                  .begin();
        cblas_dcopy(3 * nb_gauss_pts_on_faces, QUAD_2D_TABLE[rule]->points, 1,
                    shape_ptr, 1);
      } else {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
        nb_gauss_pts_on_faces = 0;
      }
    } else {
      CHKERR setGaussPtsTrianglesOnly(order_triangles_only);
      nb_gauss_pts_on_faces = gaussPtsTrianglesOnly.size2();
      if (nb_gauss_pts_on_faces == 0)
        MoFEMFunctionReturnHot(0);
      dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(
          nb_gauss_pts_on_faces, 3, false);
      if (nb_gauss_pts_on_faces) {
        CHKERR ShapeMBTRI(&*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0]
                                .getN(NOBASE)
                                .data()
                                .begin(),
                          &gaussPtsTrianglesOnly(0, 0),
                          &gaussPtsTrianglesOnly(1, 0), nb_gauss_pts_on_faces);
      }
    }
  }
  CATCH_ERRORS;

  // approx. trough prism thickness
  int nb_gauss_pts_through_thickness;
  try {
    int order_thickness = 1;
    for (unsigned int ee = 3; ee <= 5; ee++) {
      order_thickness = std::max(
          order_thickness,
          dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getDataOrder());
    }
    for (unsigned int qq = 0; qq < 3; qq++) {
      order_thickness = std::max(
          order_thickness,
          dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getDataOrder());
    }
    order_thickness = std::max(
        order_thickness,
        dataH1TroughThickness.dataOnEntities[MBPRISM][0].getDataOrder());
    // integration points
    int rule = getRuleThroughThickness(order_thickness);
    if (rule >= 0) {
      if (rule < QUAD_1D_TABLE_SIZE) {
        if (QUAD_1D_TABLE[rule]->dim != 1) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
        }
        if (QUAD_1D_TABLE[rule]->order < rule) {
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "wrong order %d != %d", QUAD_1D_TABLE[rule]->order, rule);
        }
        nb_gauss_pts_through_thickness = QUAD_1D_TABLE[rule]->npoints;
        gaussPtsThroughThickness.resize(2, nb_gauss_pts_through_thickness,
                                        false);
        cblas_dcopy(nb_gauss_pts_through_thickness,
                    &QUAD_1D_TABLE[rule]->points[1], 2,
                    &gaussPtsThroughThickness(0, 0), 1);
        cblas_dcopy(nb_gauss_pts_through_thickness,
                    QUAD_1D_TABLE[rule]->weights, 1,
                    &gaussPtsThroughThickness(1, 0), 1);
      } else {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rule > quadrature order %d < %d", rule, QUAD_1D_TABLE_SIZE);
        nb_gauss_pts_through_thickness = 0;
      }
    } else {
      CHKERR setGaussPtsThroughThickness(order_thickness);
      nb_gauss_pts_through_thickness = gaussPtsThroughThickness.size2();
    }
    if (nb_gauss_pts_through_thickness == 0)
      MoFEMFunctionReturnHot(0);
  }
  CATCH_ERRORS;

  // Generate integration pts.
  int nb_gauss_pts = nb_gauss_pts_on_faces * nb_gauss_pts_through_thickness;
  gaussPts.resize(4, nb_gauss_pts, false);
  {
    int gg = 0;
    for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {
      for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
        gaussPts(0, gg) = gaussPtsTrianglesOnly(0, ggf);
        gaussPts(1, gg) = gaussPtsTrianglesOnly(1, ggf);
        gaussPts(2, gg) = gaussPtsThroughThickness(0, ggt);
        gaussPts(3, gg) =
            gaussPtsTrianglesOnly(2, ggf) * gaussPtsThroughThickness(1, ggt);
      }
    }
  }

  {
    coordsAtGaussPtsTrianglesOnly.resize(nb_gauss_pts_on_faces, 6, false);
    for (int gg = 0; gg < nb_gauss_pts_on_faces; gg++) {
      for (int dd = 0; dd < 3; dd++) {
        coordsAtGaussPtsTrianglesOnly(gg, dd) =
            cblas_ddot(3,
                       &dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN(
                           NOBASE)(gg, 0),
                       1, &coords[dd], 3);
        coordsAtGaussPtsTrianglesOnly(gg, 3 + dd) =
            cblas_ddot(3,
                       &dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN(
                           NOBASE)(gg, 0),
                       1, &coords[9 + dd], 3);
      }
    }
    // linear for xi,eta and zeta
    dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(
        1, 6, false);
    CHKERR ShapeDiffMBTRI(&*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0]
                                .getDiffN(NOBASE)
                                .data()
                                .begin());

    // Calculate "nobase" base functions on prism, this is cartesian product
    // of base functions on triangles with base functions through thickness
    // FIXME: This could be effectively implemented with tensors
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nb_gauss_pts, 6,
                                                           false);
    dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(nb_gauss_pts, 18,
                                                               false);
    for (int dd = 0; dd != 6; dd++) {
      int gg = 0;
      for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {
        int ddd = dd > 2 ? dd - 3 : dd;
        double tri_n = dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getN(
            NOBASE)(ggf, ddd);
        double dksi_tri_n =
            dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(
                0, 2 * ddd + 0);
        double deta_tri_n =
            dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(
                0, 2 * ddd + 1);
        for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
          double zeta = gaussPtsThroughThickness(0, ggt);
          double dzeta, edge_shape;
          if (dd < 3) {
            dzeta = diffN_MBEDGE0;
            edge_shape = N_MBEDGE0(zeta);
          } else {
            dzeta = diffN_MBEDGE1;
            edge_shape = N_MBEDGE1(zeta);
          }
          dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, dd) =
              tri_n * edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(gg, 3 * dd + 0) =
              dksi_tri_n * edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(gg, 3 * dd + 1) =
              deta_tri_n * edge_shape;
          dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE)(gg, 3 * dd + 2) =
              tri_n * dzeta;
        }
      }
    }

    /// Calculate coordinates at integration points
    coordsAtGaussPts.resize(nb_gauss_pts, 3, false);
    for (int gg = 0; gg < nb_gauss_pts; gg++) {
      for (int dd = 0; dd < 3; dd++) {
        coordsAtGaussPts(gg, dd) = cblas_ddot(
            6, &dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[dd], 3);
      }
    }
  }

  // Calculate base functions on prism
  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
    if (dataH1.bAse.test(b)) {
      switch (static_cast<FieldApproximationBase>(b)) {
      case AINSWORTH_LEGENDRE_BASE:
      case AINSWORTH_LOBATTO_BASE:
        if (dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
          CHKERR FatPrismPolynomialBase().getValue(
              gaussPts,
              boost::shared_ptr<BaseFunctionCtx>(new FatPrismPolynomialBaseCtx(
                  dataH1, dataH1TrianglesOnly, dataH1TroughThickness,
                  gaussPtsTrianglesOnly, gaussPtsThroughThickness,
                  mField.get_moab(), numeredEntFiniteElementPtr.get(), H1,
                  static_cast<FieldApproximationBase>(b), NOBASE)));
        }
        if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        if (dataH1.spacesOnEntities[MBTET].test(L2)) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Not yet implemented");
        }
        break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Not yet implemented");
      }
    }
  }

  // Calculate ho-geometry tangents and normals

  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {
    hoCoordsAtGaussPtsF3.resize(nb_gauss_pts_on_faces, 3, false);
    nOrmals_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
    tAngent1_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
    tAngent2_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
    hoCoordsAtGaussPtsF4.resize(nb_gauss_pts_on_faces, 3, false);
    nOrmals_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
    tAngent1_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
    tAngent2_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
    CHKERR getNodesFieldData(dataH1TrianglesOnly, meshPositionsFieldName);
    CHKERR getEntityFieldData(dataH1TrianglesOnly, meshPositionsFieldName, MBEDGE);
    CHKERR getEntityFieldData(dataH1TrianglesOnly, meshPositionsFieldName, MBEDGE);
    CHKERR opHOCoordsAndNormals.opRhs(dataH1TrianglesOnly);
    CHKERR opHOCoordsAndNormals.calculateNormals();
  } else {
    hoCoordsAtGaussPtsF3.resize(0, 0, false);
    nOrmals_at_GaussPtF3.resize(0, 0, false);
    tAngent1_at_GaussPtF3.resize(0, 0, false);
    tAngent2_at_GaussPtF3.resize(0, 0, false);
    hoCoordsAtGaussPtsF4.resize(0, 0, false);
    nOrmals_at_GaussPtF4.resize(0, 0, false);
    tAngent1_at_GaussPtF4.resize(0, 0, false);
    tAngent2_at_GaussPtF4.resize(0, 0, false);
  }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFatPrism::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX) {

    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();

    const int nb_gauss_pts = data.getN(NOBASE).size1();
    auto t_diff_n = data.getFTensor1DiffN<3>(NOBASE);
    invJac.resize(9, nb_gauss_pts, false);
    invJac.clear();
    auto t_inv_jac = getFTensor2FromMat<3, 3>(invJac);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Tensor2<double, 3, 3> t_jac;

    double &vol =
        const_cast<VolumeElementForcesAndSourcesCore *>(getVolumeFE())->vOlume;
    vol = 0;

    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      FTensor::Tensor1<double *, 3> t_coords(coords_ptr, &coords_ptr[1],
                                             &coords_ptr[2], 3);
      t_jac(i, j) = 0;
      for (int bb = 0; bb != 6; bb++) {
        t_jac(i, j) += t_coords(i) * t_diff_n(j);
        ++t_diff_n;
        ++t_coords;
      }

      double det;
      CHKERR determinantTensor3by3(t_jac, det);
      CHKERR invertTensor3by3(t_jac, det, t_inv_jac);
      ++t_inv_jac;

      vol += 0.5 * det * getGaussPts()(3, gg);
    }
  }

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets = false;
  doPrisms = false;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacH1ForFatPrism::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    if (data.getN(base).size2() == 0)
      continue;

    const int nb_gauss_pts = data.getN(base).size1();
    auto t_diff_n = data.getFTensor1DiffN<3>(base);
    diffNinvJac.resize(data.getDiffN(base).size1(), data.getDiffN(base).size2(),
                       false);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;

    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_inv_diff_n(
        &diffNinvJac(0, 0), &diffNinvJac(0, 1), &diffNinvJac(0, 2));
    auto t_inv_jac = getFTensor2FromMat<3, 3>(invJac);

    const int nb_dofs = data.getN(base).size2();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int bb = 0; bb != nb_dofs; bb++) {
        t_inv_diff_n(i) = t_diff_n(j) * t_inv_jac(j, i);
        ++t_inv_diff_n;
        ++t_diff_n;
      }
      ++t_inv_jac;
    }

    data.getDiffN(base).data().swap(diffNinvJac.data());
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
