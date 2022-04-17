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

  auto get_fe_coordinates = [&]() {
    MoFEMFunctionBegin;
    EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
    int num_nodes;
    const EntityHandle *conn;
    CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
    coords.resize(num_nodes * 3, false);
    CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                        &*coords.data().begin());
    MoFEMFunctionReturn(0);
  };

  auto calculate_area_of_triangles = [&] {
    MoFEMFunctionBegin;
    normal.resize(6, false);
    CHKERR Tools::getTriNormal(&coords[0], &normal[0]);
    CHKERR Tools::getTriNormal(&coords[9], &normal[3]);
    aRea[0] = cblas_dnrm2(3, &normal[0], 1) * 0.5;
    aRea[1] = cblas_dnrm2(3, &normal[3], 1) * 0.5;
    MoFEMFunctionReturn(0);
  };

  CHKERR get_fe_coordinates();
  CHKERR calculate_area_of_triangles();

  CHKERR getSpacesAndBaseOnEntities(dataH1);
  CHKERR getSpacesAndBaseOnEntities(dataH1TrianglesOnly);
  CHKERR getSpacesAndBaseOnEntities(dataH1TroughThickness);

  auto get_h1_base_data = [&](auto &dataH1) {
    MoFEMFunctionBegin;
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
    MoFEMFunctionReturn(0);
  };

  // H1
  if ((dataH1.spacesOnEntities[MBEDGE]).test(H1))
    CHKERR get_h1_base_data(dataH1);

  // Hdiv
  if ((dataH1.spacesOnEntities[MBTRI]).test(HDIV))
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");

  // Hcurl
  if ((dataH1.spacesOnEntities[MBEDGE]).test(HCURL))
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");

  // L2
  if ((dataH1.spacesOnEntities[MBPRISM]).test(L2))
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented yet");

  // get approx. on triangles, i.e. faces 3 and 4
  auto set_gauss_points_on_triangle = [&](int &nb_gauss_pts_on_faces) {
    MoFEMFunctionBegin;
    int order_triangles_only = 1;
    int valid_edges[] = {1, 1, 1, 0, 0, 0, 1, 1, 1};
    for (unsigned int ee = 0; ee < 9; ee++) {
      if (!valid_edges[ee])
        continue;
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TrianglesOnly.dataOnEntities[MBEDGE][ee].getOrder());
    }
    for (unsigned int ff = 3; ff <= 4; ff++) {
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TrianglesOnly.dataOnEntities[MBTRI][ff].getOrder());
    }
    for (unsigned int qq = 0; qq < 3; qq++) {
      order_triangles_only = std::max(
          order_triangles_only,
          dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getOrder());
    }
    order_triangles_only = std::max(
        order_triangles_only,
        dataH1TroughThickness.dataOnEntities[MBPRISM][0].getOrder());

    // integration pts on the triangles surfaces
    nb_gauss_pts_on_faces = 0;
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
        dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(
            1, 6, false);
        std::copy(Tools::diffShapeFunMBTRI.begin(),
                  Tools::diffShapeFunMBTRI.end(),
                  &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0]
                        .getDiffN(NOBASE)
                        .data()
                        .begin());
      } else
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);

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
        dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(
            1, 6, false);
        std::copy(Tools::diffShapeFunMBTRI.begin(),
                  Tools::diffShapeFunMBTRI.end(),
                  &*dataH1TrianglesOnly.dataOnEntities[MBVERTEX][0]
                        .getDiffN(NOBASE)
                        .data()
                        .begin());
      }
    }
    MoFEMFunctionReturn(0);
  };

  // approx. trough prism thickness
  auto set_gauss_points_through_thickness =
      [&](int &nb_gauss_pts_through_thickness) {
        MoFEMFunctionBegin;
        nb_gauss_pts_through_thickness = 0;
        int order_thickness = 1;
        for (unsigned int ee = 3; ee <= 5; ee++) {
          order_thickness = std::max(
              order_thickness,
              dataH1TroughThickness.dataOnEntities[MBEDGE][ee].getOrder());
        }
        for (unsigned int qq = 0; qq < 3; qq++) {
          order_thickness = std::max(
              order_thickness,
              dataH1TroughThickness.dataOnEntities[MBQUAD][qq].getOrder());
        }
        order_thickness = std::max(
            order_thickness,
            dataH1TroughThickness.dataOnEntities[MBPRISM][0].getOrder());
        // integration points
        int rule = getRuleThroughThickness(order_thickness);
        if (rule >= 0) {
          if (rule < QUAD_1D_TABLE_SIZE) {
            if (QUAD_1D_TABLE[rule]->dim != 1) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                      "wrong dimension");
            }
            if (QUAD_1D_TABLE[rule]->order < rule) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                       "wrong order %d != %d", QUAD_1D_TABLE[rule]->order,
                       rule);
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
                     "rule > quadrature order %d < %d", rule,
                     QUAD_1D_TABLE_SIZE);
            nb_gauss_pts_through_thickness = 0;
          }
        } else {
          CHKERR setGaussPtsThroughThickness(order_thickness);
          nb_gauss_pts_through_thickness = gaussPtsThroughThickness.size2();
        }
        MoFEMFunctionReturn(0);
      };

  // Generate integration pts.
  auto set_gauss_points_in_volume = [&](int nb_gauss_pts_on_faces,
                                        int nb_gauss_pts_through_thickness,
                                        int &nb_gauss_pts) {
    MoFEMFunctionBegin;
    nb_gauss_pts = nb_gauss_pts_on_faces * nb_gauss_pts_through_thickness;
    gaussPts.resize(4, nb_gauss_pts, false);
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
    MoFEMFunctionReturn(0);
  };

  int nb_gauss_pts, nb_gauss_pts_through_thickness, nb_gauss_pts_on_faces;
  CHKERR set_gauss_points_on_triangle(nb_gauss_pts_on_faces);
  if (!nb_gauss_pts_on_faces)
    MoFEMFunctionReturnHot(0);
  CHKERR set_gauss_points_through_thickness(nb_gauss_pts_through_thickness);
  if (!nb_gauss_pts_through_thickness)
    MoFEMFunctionReturnHot(0);
  CHKERR set_gauss_points_in_volume(
      nb_gauss_pts_on_faces, nb_gauss_pts_through_thickness, nb_gauss_pts);

  auto calc_coordinates_at_triangles = [&]() {
    MoFEMFunctionBegin;
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
    MoFEMFunctionReturn(0);
  };

  auto calc_vertex_base_on_prism = [&]() {
    MoFEMFunctionBegin;
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
    MoFEMFunctionReturn(0);
  };

  auto calc_base_on_prism = [&]() {
    MoFEMFunctionBegin;
    // Calculate base functions on prism
    for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
      if (dataH1.bAse.test(b)) {
        switch (static_cast<FieldApproximationBase>(b)) {
        case AINSWORTH_LEGENDRE_BASE:
        case AINSWORTH_LOBATTO_BASE:
          if (dataH1.spacesOnEntities[MBVERTEX].test(H1)) {
            CHKERR FatPrismPolynomialBase().getValue(
                gaussPts,
                boost::shared_ptr<BaseFunctionCtx>(
                    new FatPrismPolynomialBaseCtx(
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
    MoFEMFunctionReturn(0);
  };

  auto calc_coordinate_on_prism = [&]() {
    MoFEMFunctionBegin;
    /// Calculate coordinates at integration points
    coordsAtGaussPts.resize(nb_gauss_pts, 3, false);
    for (int gg = 0; gg < nb_gauss_pts; gg++) {
      for (int dd = 0; dd < 3; dd++) {
        coordsAtGaussPts(gg, dd) = cblas_ddot(
            6, &dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE)(gg, 0), 1,
            &coords[dd], 3);
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto calculate_volume = [&]() {
    auto get_t_w = [&] {
      return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(&gaussPts(2, 0));
    };

    auto get_t_coords = [&]() {
      return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
          &coords[0], &coords[1], &coords[2]);
    };

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;
    FTensor::Tensor2<double, 3, 3> t_jac;

    const size_t nb_gauss_pts = gaussPts.size2();
    auto t_diff_n =
        dataH1.dataOnEntities[MBVERTEX][0].getFTensor1DiffN<3>(NOBASE);

    double vol = 0;
    auto t_w = get_t_w();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {

      auto t_coords = get_t_coords();
      t_jac(i, j) = 0;
      for (size_t n = 0; n != 6; ++n) {
        t_jac(i, j) += t_coords(i) * t_diff_n(j);
        ++t_diff_n;
        ++t_coords;
      }

      double det;
      CHKERR determinantTensor3by3(t_jac, det);
      vol += det * t_w / 2;

      ++t_w;
    }

    return vol;
  };

  auto calc_ho_triangle_face_normals = [&]() {
    MoFEMFunctionBegin;

    auto check_field = [&]() {
      auto field_it =
          fieldsPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName);
      if (field_it != fieldsPtr->get<FieldName_mi_tag>().end())
        if ((numeredEntFiniteElementPtr->getBitFieldIdData() &
             (*field_it)->getId())
                .any())
          return true;
      return false;
    };

    // Check if field meshPositionsFieldName exist
    if (check_field()) {
      hoCoordsAtGaussPtsF3.resize(nb_gauss_pts_on_faces, 3, false);
      nOrmals_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
      tAngent1_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
      tAngent2_at_GaussPtF3.resize(nb_gauss_pts_on_faces, 3, false);
      hoCoordsAtGaussPtsF4.resize(nb_gauss_pts_on_faces, 3, false);
      nOrmals_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
      tAngent1_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
      tAngent2_at_GaussPtF4.resize(nb_gauss_pts_on_faces, 3, false);
      const auto bit_number =
          mField.get_field_bit_number(meshPositionsFieldName);
      CHKERR getNodesFieldData(dataH1TrianglesOnly, bit_number);
      CHKERR getEntityFieldData(dataH1TrianglesOnly, bit_number, MBEDGE);
      CHKERR getEntityFieldData(dataH1TrianglesOnly, bit_number, MBEDGE);
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
    MoFEMFunctionReturn(0);
  };

  CHKERR calc_coordinates_at_triangles();
  CHKERR calc_vertex_base_on_prism();
  CHKERR calc_base_on_prism();
  CHKERR calc_coordinate_on_prism();
  vOlume = calculate_volume();
  CHKERR calc_ho_triangle_face_normals();

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FatPrismElementForcesAndSourcesCore::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<FatPrismElementForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
