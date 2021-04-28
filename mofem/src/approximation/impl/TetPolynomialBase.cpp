/** \file TetPolynomialBase.cpp
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

using namespace MoFEM;

MoFEMErrorCode
TetPolynomialBase::query_interface(const MOFEMuuid &uuid,
                                   BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_TET_BASE_FUNCTION) {
    *iface = const_cast<TetPolynomialBase *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = BaseFunction::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

TetPolynomialBase::~TetPolynomialBase() {}
TetPolynomialBase::TetPolynomialBase() {}

MoFEMErrorCode TetPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    CHKERR getValueH1AinsworthBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
    CHKERR getValueH1BernsteinBezierBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueH1AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  int sense[6], order[6];
  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    // edges
    if (data.dataOnEntities[MBEDGE].size() != 6) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *h1_edge_n[6], *diff_h1_egde_n[6];
    for (int ee = 0; ee != 6; ++ee) {
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
      diff_h1_egde_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR H1_EdgeShapeFunctions_MBTET(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        h1_edge_n, diff_h1_egde_n, nb_gauss_pts, base_polynomials);
  } else {
    for (int ee = 0; ee != 6; ++ee) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(0, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(0, 0, false);
    }
  }

  if (data.spacesOnEntities[MBTRI].test(H1)) {
    // faces
    if (data.dataOnEntities[MBTRI].size() != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *h1_face_n[4], *diff_h1_face_n[4];
    for (int ff = 0; ff != 4; ++ff) {
      if (data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      int nb_dofs = NBFACETRI_H1(data.dataOnEntities[MBTRI][ff].getDataOrder());
      order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                       false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,
                                                           3 * nb_dofs, false);
      h1_face_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      diff_h1_face_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if (data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    if (data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    CHKERR H1_FaceShapeFunctions_MBTET(
        &*data.facesNodes.data().begin(), order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        h1_face_n, diff_h1_face_n, nb_gauss_pts, base_polynomials);

  } else {
    for (int ff = 0; ff != 4; ++ff) {
      data.dataOnEntities[MBTRI][ff].getN(base).resize(0, false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(0, 0, false);
    }
  }

  if (data.spacesOnEntities[MBTET].test(H1)) {
    // volume
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_H1(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, nb_vol_dofs,
                                                    false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        3 * nb_vol_dofs, false);
    CHKERR H1_VolumeShapeFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getDataOrder(),
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
        nb_gauss_pts, base_polynomials);
  } else {
    data.dataOnEntities[MBTET][0].getN(base).resize(0, 0, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(0, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueH1BernsteinBezierBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const std::string field_name = cTx->fieldName;
  const int nb_gauss_pts = pts.size2();

  if (data.dataOnEntities[MBVERTEX][0].getN(NOBASE).size1() !=
      (unsigned int)nb_gauss_pts)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[NOBASE]);
  auto &lambda = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);

  auto get_alpha = [field_name](auto &data) -> MatrixInt & {
    auto &ptr = data.getBBAlphaIndicesSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixInt());
    return *ptr;
  };

  auto get_base = [field_name](auto &data) -> MatrixDouble & {
    auto &ptr = data.getBBNSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixDouble());
    return *ptr;
  };

  auto get_diff_base = [field_name](auto &data) -> MatrixDouble & {
    auto &ptr = data.getBBDiffNSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixDouble());
    return *ptr;
  };

  auto get_alpha_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixInt> & {
    return data.getBBAlphaIndicesSharedPtr(field_name);
  };

  auto get_base_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBNSharedPtr(field_name);
  };

  auto get_diff_base_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBDiffNSharedPtr(field_name);
  };

  auto get_alpha_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixInt> & {
    return data.getBBAlphaIndicesByOrderSharedPtr(o);
  };

  auto get_base_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBNByOrderSharedPtr(o);
  };

  auto get_diff_base_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBDiffNByOrderSharedPtr(o);
  };

  auto &vert_ent_data = data.dataOnEntities[MBVERTEX][0];
  auto &vertex_alpha = get_alpha(vert_ent_data);
  vertex_alpha.resize(4, 4, false);
  vertex_alpha.clear();
  for (int n = 0; n != 4; ++n)
    vertex_alpha(n, n) = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[n];

  auto &vert_get_n = get_base(vert_ent_data);
  auto &vert_get_diff_n = get_diff_base(vert_ent_data);
  vert_get_n.resize(nb_gauss_pts, 4, false);
  vert_get_diff_n.resize(nb_gauss_pts, 12, false);
  CHKERR BernsteinBezier::baseFunctionsTet(
      1, lambda.size1(), vertex_alpha.size1(), &vertex_alpha(0, 0),
      &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &vert_get_n(0, 0),
      &vert_get_diff_n(0, 0));
  for (int n = 0; n != 4; ++n) {
    const double f = boost::math::factorial<double>(
        data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[n]);
    for (int g = 0; g != nb_gauss_pts; ++g) {
      vert_get_n(g, n) *= f;
      for (int d = 0; d != 3; ++d)
        vert_get_diff_n(g, 3 * n + d) *= f;
    }
  }

  // edges
  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    if (data.dataOnEntities[MBEDGE].size() != 6)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size of ent data");

    constexpr int edges_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                       {0, 3}, {1, 3}, {2, 3}};
    for (int ee = 0; ee != 6; ++ee) {
      auto &ent_data = data.dataOnEntities[MBEDGE][ee];
      const int sense = ent_data.getSense();
      if (sense == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Sense of the edge unknown");
      const int order = ent_data.getDataOrder();
      const int nb_dofs = NBEDGE_H1(order);

      if (nb_dofs) {
        if (get_alpha_by_order_ptr(ent_data, order)) {
          get_alpha_by_name_ptr(ent_data, field_name) =
              get_alpha_by_order_ptr(ent_data, order);
          get_base_by_name_ptr(ent_data, field_name) =
              get_base_by_order_ptr(ent_data, order);
          get_diff_base_by_name_ptr(ent_data, field_name) =
              get_diff_base_by_order_ptr(ent_data, order);
        } else {
          auto &get_n = get_base(ent_data);
          auto &get_diff_n = get_diff_base(ent_data);
          get_n.resize(nb_gauss_pts, nb_dofs, false);
          get_diff_n.resize(nb_gauss_pts, 3 * nb_dofs, false);

          auto &edge_alpha = get_alpha(data.dataOnEntities[MBEDGE][ee]);
          edge_alpha.resize(nb_dofs, 4, false);
          CHKERR BernsteinBezier::generateIndicesEdgeTet(ee, order,
                                                         &edge_alpha(0, 0));
          if (sense == -1) {
            for (int i = 0; i != edge_alpha.size1(); ++i) {
              int a = edge_alpha(i, edges_nodes[ee][0]);
              edge_alpha(i, edges_nodes[ee][0]) =
                  edge_alpha(i, edges_nodes[ee][1]);
              edge_alpha(i, edges_nodes[ee][1]) = a;
            }
          }
          CHKERR BernsteinBezier::baseFunctionsTet(
              order, lambda.size1(), edge_alpha.size1(), &edge_alpha(0, 0),
              &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &get_n(0, 0),
              &get_diff_n(0, 0));

          get_alpha_by_order_ptr(ent_data, order) =
              get_alpha_by_name_ptr(ent_data, field_name);
          get_base_by_order_ptr(ent_data, order) =
              get_base_by_name_ptr(ent_data, field_name);
          get_diff_base_by_order_ptr(ent_data, order) =
              get_diff_base_by_name_ptr(ent_data, field_name);
        }
      }
    }
  } else {
    for (int ee = 0; ee != 6; ++ee) {
      auto &ent_data = data.dataOnEntities[MBEDGE][ee];
      ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
      auto &get_n = get_base(ent_data);
      auto &get_diff_n = get_diff_base(ent_data);
      get_n.resize(nb_gauss_pts, 0, false);
      get_diff_n.resize(nb_gauss_pts, 0, false);
    }
  }

  // face
  if (data.spacesOnEntities[MBTRI].test(H1)) {
    if (data.dataOnEntities[MBTRI].size() != 4)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size of ent data");
    if (data.facesNodes.size1() != 4)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    if (data.facesNodes.size2() != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

    for (int ff = 0; ff != 4; ++ff) {
      auto &ent_data = data.dataOnEntities[MBTRI][ff];
      const int order = ent_data.getDataOrder();
      const int nb_dofs = NBFACETRI_H1(order);

      if (nb_dofs) {
        if (get_alpha_by_order_ptr(ent_data, order)) {
          get_alpha_by_name_ptr(ent_data, field_name) =
              get_alpha_by_order_ptr(ent_data, order);
          get_base_by_name_ptr(ent_data, field_name) =
              get_base_by_order_ptr(ent_data, order);
          get_diff_base_by_name_ptr(ent_data, field_name) =
              get_diff_base_by_order_ptr(ent_data, order);
        } else {

          auto &get_n = get_base(ent_data);
          auto &get_diff_n = get_diff_base(ent_data);
          get_n.resize(nb_gauss_pts, nb_dofs, false);
          get_diff_n.resize(nb_gauss_pts, 3 * nb_dofs, false);

          auto &face_alpha = get_alpha(ent_data);
          face_alpha.resize(nb_dofs, 4, false);

          CHKERR BernsteinBezier::generateIndicesTriTet(ff, order,
                                                        &face_alpha(0, 0));
          senseFaceAlpha.resize(face_alpha.size1(), face_alpha.size2(), false);
          senseFaceAlpha.clear();
          constexpr int tri_nodes[4][3] = {
              {0, 1, 3}, {1, 2, 3}, {0, 2, 3}, {0, 1, 2}};
          for (int d = 0; d != nb_dofs; ++d)
            for (int n = 0; n != 3; ++n)
              senseFaceAlpha(d, data.facesNodes(ff, n)) =
                  face_alpha(d, tri_nodes[ff][n]);
          face_alpha.swap(senseFaceAlpha);
          CHKERR BernsteinBezier::baseFunctionsTet(
              order, lambda.size1(), face_alpha.size1(), &face_alpha(0, 0),
              &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &get_n(0, 0),
              &get_diff_n(0, 0));

          get_alpha_by_order_ptr(ent_data, order) =
              get_alpha_by_name_ptr(ent_data, field_name);
          get_base_by_order_ptr(ent_data, order) =
              get_base_by_name_ptr(ent_data, field_name);
          get_diff_base_by_order_ptr(ent_data, order) =
              get_diff_base_by_name_ptr(ent_data, field_name);
        }
      }
    }
  } else {
    for (int ff = 0; ff != 4; ++ff) {
      auto &ent_data = data.dataOnEntities[MBTRI][ff];
      ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
      auto &get_n = get_base(ent_data);
      auto &get_diff_n = get_diff_base(ent_data);
      get_n.resize(nb_gauss_pts, 0, false);
      get_diff_n.resize(nb_gauss_pts, 0, false);
    }
  }

  if (data.spacesOnEntities[MBTET].test(H1)) {
    if (data.dataOnEntities[MBTET].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size ent of ent data");

    auto &ent_data = data.dataOnEntities[MBTET][0];
    const int order = ent_data.getDataOrder();
    const int nb_dofs = NBVOLUMETET_H1(order);
    if (get_alpha_by_order_ptr(ent_data, order)) {
      get_alpha_by_name_ptr(ent_data, field_name) =
          get_alpha_by_order_ptr(ent_data, order);
      get_base_by_name_ptr(ent_data, field_name) =
          get_base_by_order_ptr(ent_data, order);
      get_diff_base_by_name_ptr(ent_data, field_name) =
          get_diff_base_by_order_ptr(ent_data, order);
    } else {

      auto &get_n = get_base(ent_data);
      auto &get_diff_n = get_diff_base(ent_data);
      get_n.resize(nb_gauss_pts, nb_dofs, false);
      get_diff_n.resize(nb_gauss_pts, 3 * nb_dofs, false);
      if (nb_dofs) {
        auto &tet_alpha = get_alpha(ent_data);
        tet_alpha.resize(nb_dofs, 4, false);

        CHKERR BernsteinBezier::generateIndicesTetTet(order, &tet_alpha(0, 0));
        CHKERR BernsteinBezier::baseFunctionsTet(
            order, lambda.size1(), tet_alpha.size1(), &tet_alpha(0, 0),
            &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &get_n(0, 0),
            &get_diff_n(0, 0));

        get_alpha_by_order_ptr(ent_data, order) =
            get_alpha_by_name_ptr(ent_data, field_name);
        get_base_by_order_ptr(ent_data, order) =
            get_base_by_name_ptr(ent_data, field_name);
        get_diff_base_by_order_ptr(ent_data, order) =
            get_diff_base_by_name_ptr(ent_data, field_name);
      }
    }
  } else {
    auto &ent_data = data.dataOnEntities[MBTET][0];
    ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
    auto &get_n = get_base(ent_data);
    auto &get_diff_n = get_diff_base(ent_data);
    get_n.resize(nb_gauss_pts, 0, false);
    get_diff_n.resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    CHKERR getValueL2AinsworthBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
    CHKERR getValueL2BernsteinBezierBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueL2AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  data.dataOnEntities[MBTET][0].getN(base).resize(
      nb_gauss_pts,
      NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getDataOrder()), false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(
      nb_gauss_pts,
      3 * NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getDataOrder()), false);

  CHKERR L2_Ainsworth_ShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
      nb_gauss_pts, base_polynomials);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueL2BernsteinBezierBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const std::string field_name = cTx->fieldName;
  const int nb_gauss_pts = pts.size2();

  if (data.dataOnEntities[MBVERTEX][0].getN(NOBASE).size1() !=
      (unsigned int)nb_gauss_pts)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[NOBASE]);
  auto &lambda = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);

  auto get_alpha = [field_name](auto &data) -> MatrixInt & {
    auto &ptr = data.getBBAlphaIndicesSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixInt());
    return *ptr;
  };

  auto get_base = [field_name](auto &data) -> MatrixDouble & {
    auto &ptr = data.getBBNSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixDouble());
    return *ptr;
  };

  auto get_diff_base = [field_name](auto &data) -> MatrixDouble & {
    auto &ptr = data.getBBDiffNSharedPtr(field_name);
    if (!ptr)
      ptr.reset(new MatrixDouble());
    return *ptr;
  };

  auto get_alpha_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixInt> & {
    return data.getBBAlphaIndicesSharedPtr(field_name);
  };

  auto get_base_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBNSharedPtr(field_name);
  };

  auto get_diff_base_by_name_ptr =
      [](auto &data,
         const std::string &field_name) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBDiffNSharedPtr(field_name);
  };

  auto get_alpha_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixInt> & {
    return data.getBBAlphaIndicesByOrderSharedPtr(o);
  };

  auto get_base_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBNByOrderSharedPtr(o);
  };

  auto get_diff_base_by_order_ptr =
      [](auto &data, const size_t o) -> boost::shared_ptr<MatrixDouble> & {
    return data.getBBDiffNByOrderSharedPtr(o);
  };

  if (data.spacesOnEntities[MBTET].test(L2)) {
    if (data.dataOnEntities[MBTET].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size ent of ent data");

    auto &ent_data = data.dataOnEntities[MBTET][0];
    const int order = ent_data.getDataOrder();
    const int nb_dofs = NBVOLUMETET_L2(order);

    if (get_alpha_by_order_ptr(ent_data, order)) {
      get_alpha_by_name_ptr(ent_data, field_name) =
          get_alpha_by_order_ptr(ent_data, order);
      get_base_by_name_ptr(ent_data, field_name) =
          get_base_by_order_ptr(ent_data, order);
      get_diff_base_by_name_ptr(ent_data, field_name) =
          get_diff_base_by_order_ptr(ent_data, order);
    } else {

      auto &get_n = get_base(ent_data);
      auto &get_diff_n = get_diff_base(ent_data);
      get_n.resize(nb_gauss_pts, nb_dofs, false);
      get_diff_n.resize(nb_gauss_pts, 3 * nb_dofs, false);

      if (nb_dofs) {

        if (order == 0) {

          if (nb_dofs != 1)
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Inconsistent number of DOFs");

          auto &tri_alpha = get_alpha(ent_data);
          tri_alpha.clear();
          get_n(0, 0) = 1;
          get_diff_n.clear();

        } else {

          if (nb_dofs != 4 + 6 * NBEDGE_H1(order) + 4 * NBFACETRI_H1(order) +
                             NBVOLUMETET_H1(order) ||
              nb_dofs != NBVOLUMETET_L2(order))
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Inconsistent number of DOFs");


          auto &tet_alpha = get_alpha(ent_data);
          tet_alpha.resize(nb_dofs, 4, false);

          CHKERR BernsteinBezier::generateIndicesVertexTet(order,
                                                           &tet_alpha(0, 0));
          if(order > 1) {
            std::array<int, 6> edge_n{order, order, order, order, order, order};
            std::array<int *, 6> tet_edge_ptr{
                &tet_alpha(4, 0),
                &tet_alpha(4 + 1 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 2 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 3 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 4 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 5 * NBEDGE_H1(order), 0)};
            CHKERR BernsteinBezier::generateIndicesEdgeTri(
                edge_n.data(), tet_edge_ptr.data());
            if (order > 2) {
              std::array<int, 6> face_n{order, order, order, order};
              std::array<int *, 6> tet_face_ptr{
                  &tet_alpha(4 + 6 * NBEDGE_H1(order), 0),
                  &tet_alpha(4 + 6 * NBEDGE_H1(order) + 1 * NBFACETRI_H1(order),
                             0),
                  &tet_alpha(4 + 6 * NBEDGE_H1(order) + 2 * NBFACETRI_H1(order),
                             0),
                  &tet_alpha(4 + 6 * NBEDGE_H1(order) + 3 * NBFACETRI_H1(order),
                             0),
              };
              CHKERR BernsteinBezier::generateIndicesTriTet(
                  face_n.data(), tet_face_ptr.data());
              if(order > 3)
                CHKERR BernsteinBezier::generateIndicesTetTet(
                    order,
                    &tet_alpha(
                        4 + 6 * NBEDGE_H1(order) + 4 * NBFACETRI_H1(order), 0));
            }
          }

          CHKERR BernsteinBezier::baseFunctionsTet(
              order, lambda.size1(), tet_alpha.size1(), &tet_alpha(0, 0),
              &lambda(0, 0), Tools::diffShapeFunMBTET.data(), &get_n(0, 0),
              &get_diff_n(0, 0));

          get_alpha_by_order_ptr(ent_data, order) =
              get_alpha_by_name_ptr(ent_data, field_name);
          get_base_by_order_ptr(ent_data, order) =
              get_base_by_name_ptr(ent_data, field_name);
          get_diff_base_by_order_ptr(ent_data, order) =
              get_diff_base_by_name_ptr(ent_data, field_name);
        }
      }
    }
  } else {
    auto &ent_data = data.dataOnEntities[MBTET][0];
    ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
    auto &get_n = get_base(ent_data);
    auto &get_diff_n = get_diff_base(ent_data);
    get_n.resize(nb_gauss_pts, 0, false);
    get_diff_n.resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdivAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  // face shape functions

  double *phi_f_e[4][3];
  double *phi_f[4];
  double *diff_phi_f_e[4][3];
  double *diff_phi_f[4];

  N_face_edge.resize(4, 3, false);
  N_face_bubble.resize(4, false);
  diffN_face_edge.resize(4, 3, false);
  diffN_face_bubble.resize(4, false);

  int faces_order[4];
  for (int ff = 0; ff != 4; ++ff) {
    if (data.dataOnEntities[MBTRI][ff].getSense() == 0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    faces_order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
    // three edges on face
    for (int ee = 0; ee < 3; ee++) {
      N_face_edge(ff, ee).resize(
          nb_gauss_pts, 3 * NBFACETRI_AINSWORTH_EDGE_HDIV(faces_order[ff]),
          false);
      diffN_face_edge(ff, ee).resize(
          nb_gauss_pts, 9 * NBFACETRI_AINSWORTH_EDGE_HDIV(faces_order[ff]),
          false);
      phi_f_e[ff][ee] = &*N_face_edge(ff, ee).data().begin();
      diff_phi_f_e[ff][ee] = &*diffN_face_edge(ff, ee).data().begin();
    }
    N_face_bubble[ff].resize(nb_gauss_pts,
                             3 * NBFACETRI_AINSWORTH_FACE_HDIV(faces_order[ff]),
                             false);
    diffN_face_bubble[ff].resize(
        nb_gauss_pts, 9 * NBFACETRI_AINSWORTH_FACE_HDIV(faces_order[ff]),
        false);
    phi_f[ff] = &*(N_face_bubble[ff].data().begin());
    diff_phi_f[ff] = &*(diffN_face_bubble[ff].data().begin());
  }

  CHKERR Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET(
      &data.facesNodes(0, 0), faces_order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(), phi_f_e,
      diff_phi_f_e, nb_gauss_pts, base_polynomials);

  CHKERR Hdiv_Ainsworth_FaceBubbleShapeFunctions(
      &data.facesNodes(0, 0), faces_order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(), phi_f,
      diff_phi_f, nb_gauss_pts, base_polynomials);

  // volume shape functions

  double *phi_v_e[6];
  double *phi_v_f[4];
  double *phi_v;
  double *diff_phi_v_e[6];
  double *diff_phi_v_f[4];
  double *diff_phi_v;

  int volume_order = data.dataOnEntities[MBTET][0].getDataOrder();

  N_volume_edge.resize(6, false);
  diffN_volume_edge.resize(6, false);
  for (int ee = 0; ee != 6; ++ee) {
    N_volume_edge[ee].resize(
        nb_gauss_pts, 3 * NBVOLUMETET_AINSWORTH_EDGE_HDIV(volume_order), false);
    diffN_volume_edge[ee].resize(
        nb_gauss_pts, 9 * NBVOLUMETET_AINSWORTH_EDGE_HDIV(volume_order), false);
    phi_v_e[ee] = &*(N_volume_edge[ee].data().begin());
    diff_phi_v_e[ee] = &*(diffN_volume_edge[ee].data().begin());
  }
  CHKERR Hdiv_Ainsworth_EdgeBasedVolumeShapeFunctions_MBTET(
      volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), phi_v_e,
      diff_phi_v_e, nb_gauss_pts, base_polynomials);

  N_volume_face.resize(4, false);
  diffN_volume_face.resize(4, false);
  for (int ff = 0; ff != 4; ++ff) {
    N_volume_face[ff].resize(
        nb_gauss_pts, 3 * NBVOLUMETET_AINSWORTH_FACE_HDIV(volume_order), false);
    diffN_volume_face[ff].resize(
        nb_gauss_pts, 9 * NBVOLUMETET_AINSWORTH_FACE_HDIV(volume_order), false);
    phi_v_f[ff] = &*(N_volume_face[ff].data().begin());
    diff_phi_v_f[ff] = &*(diffN_volume_face[ff].data().begin());
  }
  CHKERR Hdiv_Ainsworth_FaceBasedVolumeShapeFunctions_MBTET(
      volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), phi_v_f,
      diff_phi_v_f, nb_gauss_pts, base_polynomials);

  N_volume_bubble.resize(
      nb_gauss_pts, 3 * NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order), false);
  diffN_volume_bubble.resize(
      nb_gauss_pts, 9 * NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order), false);
  phi_v = &*(N_volume_bubble.data().begin());
  diff_phi_v = &*(diffN_volume_bubble.data().begin());
  CHKERR Hdiv_Ainsworth_VolumeBubbleShapeFunctions_MBTET(
      volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), phi_v, diff_phi_v,
      nb_gauss_pts, base_polynomials);

  // Set shape functions into data structure Shape functions hast to be put
  // in arrays in order which guarantee hierarchical series of degrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  // faces
  if (data.dataOnEntities[MBTRI].size() != 4) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  for (int ff = 0; ff != 4; ff++) {
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_AINSWORTH_HDIV(faces_order[ff]), false);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_AINSWORTH_HDIV(faces_order[ff]), false);
    if (NBFACETRI_AINSWORTH_HDIV(faces_order[ff]) == 0)
      continue;
    // face
    double *base_ptr =
        &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
    FTensor::Tensor1<double *, 3> t_base(base_ptr, &base_ptr[HVEC1],
                                         &base_ptr[HVEC2], 3);
    double *diff_base_ptr =
        &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    FTensor::Tensor2<double *, 3, 3> t_diff_base(
        &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
        &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
        &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
        &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
        &diff_base_ptr[HVEC2_2], 9);
    // face-face
    boost::shared_ptr<FTensor::Tensor1<double *, 3>> t_base_f;
    boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>> t_diff_base_f;
    if (NBFACETRI_AINSWORTH_FACE_HDIV(faces_order[ff]) > 0) {
      base_ptr = phi_f[ff];
      t_base_f = boost::shared_ptr<FTensor::Tensor1<double *, 3>>(
          new FTensor::Tensor1<double *, 3>(base_ptr, &base_ptr[HVEC1],
                                            &base_ptr[HVEC2], 3));
      diff_base_ptr = diff_phi_f[ff];
      t_diff_base_f = boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>>(
          new FTensor::Tensor2<double *, 3, 3>(
              &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
              &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
              &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
              &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
              &diff_base_ptr[HVEC2_2], 9));
    }
    // edge-face
    base_ptr = phi_f_e[ff][0];
    FTensor::Tensor1<double *, 3> t_base_f_e0(base_ptr, &base_ptr[HVEC1],
                                              &base_ptr[HVEC2], 3);
    diff_base_ptr = diff_phi_f_e[ff][0];
    FTensor::Tensor2<double *, 3, 3> t_diff_base_f_e0(
        &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
        &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
        &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
        &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
        &diff_base_ptr[HVEC2_2], 9);
    base_ptr = phi_f_e[ff][1];
    FTensor::Tensor1<double *, 3> t_base_f_e1(base_ptr, &base_ptr[HVEC1],
                                              &base_ptr[HVEC2], 3);
    diff_base_ptr = diff_phi_f_e[ff][1];
    FTensor::Tensor2<double *, 3, 3> t_diff_base_f_e1(
        &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
        &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
        &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
        &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
        &diff_base_ptr[HVEC2_2], 9);
    base_ptr = phi_f_e[ff][2];
    FTensor::Tensor1<double *, 3> t_base_f_e2(base_ptr, &base_ptr[HVEC1],
                                              &base_ptr[HVEC2], 3);
    diff_base_ptr = diff_phi_f_e[ff][2];
    FTensor::Tensor2<double *, 3, 3> t_diff_base_f_e2(
        &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
        &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
        &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
        &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
        &diff_base_ptr[HVEC2_2], 9);
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int oo = 0; oo != faces_order[ff]; oo++) {
        for (int dd = NBFACETRI_AINSWORTH_EDGE_HDIV(oo);
             dd != NBFACETRI_AINSWORTH_EDGE_HDIV(oo + 1); dd++) {
          t_base(i) = t_base_f_e0(i);
          ++t_base;
          ++t_base_f_e0;
          t_diff_base(i, j) = t_diff_base_f_e0(i, j);
          ++t_diff_base;
          ++t_diff_base_f_e0;
          t_base(i) = t_base_f_e1(i);
          ++t_base;
          ++t_base_f_e1;
          t_diff_base(i, j) = t_diff_base_f_e1(i, j);
          ++t_diff_base;
          ++t_diff_base_f_e1;
          t_base(i) = t_base_f_e2(i);
          ++t_base;
          ++t_base_f_e2;
          t_diff_base(i, j) = t_diff_base_f_e2(i, j);
          ++t_diff_base;
          ++t_diff_base_f_e2;
        }
        for (int dd = NBFACETRI_AINSWORTH_FACE_HDIV(oo);
             dd != NBFACETRI_AINSWORTH_FACE_HDIV(oo + 1); dd++) {
          t_base(i) = (*t_base_f)(i);
          ++t_base;
          ++(*t_base_f);
          t_diff_base(i, j) = (*t_diff_base_f)(i, j);
          ++t_diff_base;
          ++(*t_diff_base_f);
        }
      }
    }
  }

  // volume
  data.dataOnEntities[MBTET][0].getN(base).resize(
      nb_gauss_pts, 3 * NBVOLUMETET_AINSWORTH_HDIV(volume_order), false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(
      nb_gauss_pts, 9 * NBVOLUMETET_AINSWORTH_HDIV(volume_order), false);
  if (NBVOLUMETET_AINSWORTH_HDIV(volume_order) > 0) {
    double *base_ptr =
        &*data.dataOnEntities[MBTET][0].getN(base).data().begin();
    FTensor::Tensor1<double *, 3> t_base(base_ptr, &base_ptr[HVEC1],
                                         &base_ptr[HVEC2], 3);
    double *diff_base_ptr =
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();
    FTensor::Tensor2<double *, 3, 3> t_diff_base(
        &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
        &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
        &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
        &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
        &diff_base_ptr[HVEC2_2], 9);
    // edges
    std::vector<FTensor::Tensor1<double *, 3>> t_base_v_e;
    t_base_v_e.reserve(6);
    std::vector<FTensor::Tensor2<double *, 3, 3>> t_diff_base_v_e;
    t_diff_base_v_e.reserve(6);
    for (int ee = 0; ee != 6; ee++) {
      base_ptr = phi_v_e[ee];
      diff_base_ptr = diff_phi_v_e[ee];
      t_base_v_e.push_back(FTensor::Tensor1<double *, 3>(
          base_ptr, &base_ptr[HVEC1], &base_ptr[HVEC2], 3));
      t_diff_base_v_e.push_back(FTensor::Tensor2<double *, 3, 3>(
          &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
          &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
          &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
          &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
          &diff_base_ptr[HVEC2_2], 9));
    }
    // faces
    std::vector<FTensor::Tensor1<double *, 3>> t_base_v_f;
    t_base_v_f.reserve(4);
    std::vector<FTensor::Tensor2<double *, 3, 3>> t_diff_base_v_f;
    t_diff_base_v_f.reserve(4);
    if (NBVOLUMETET_AINSWORTH_FACE_HDIV(volume_order) > 0) {
      for (int ff = 0; ff != 4; ff++) {
        base_ptr = phi_v_f[ff];
        diff_base_ptr = diff_phi_v_f[ff];
        t_base_v_f.push_back(FTensor::Tensor1<double *, 3>(
            base_ptr, &base_ptr[HVEC1], &base_ptr[HVEC2], 3));
        t_diff_base_v_f.push_back(FTensor::Tensor2<double *, 3, 3>(
            &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
            &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
            &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
            &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
            &diff_base_ptr[HVEC2_2], 9));
      }
    }
    boost::shared_ptr<FTensor::Tensor1<double *, 3>> t_base_v;
    boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>> t_diff_base_v;
    if (NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order) > 0) {
      base_ptr = phi_v;
      t_base_v = boost::shared_ptr<FTensor::Tensor1<double *, 3>>(
          new FTensor::Tensor1<double *, 3>(base_ptr, &base_ptr[HVEC1],
                                            &base_ptr[HVEC2], 3));
      diff_base_ptr = diff_phi_v;
      t_diff_base_v = boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>>(
          new FTensor::Tensor2<double *, 3, 3>(
              &diff_base_ptr[HVEC0_0], &diff_base_ptr[HVEC0_1],
              &diff_base_ptr[HVEC0_2], &diff_base_ptr[HVEC1_0],
              &diff_base_ptr[HVEC1_1], &diff_base_ptr[HVEC1_2],
              &diff_base_ptr[HVEC2_0], &diff_base_ptr[HVEC2_1],
              &diff_base_ptr[HVEC2_2], 9));
    }
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int oo = 0; oo < volume_order; oo++) {
        for (int dd = NBVOLUMETET_AINSWORTH_EDGE_HDIV(oo);
             dd < NBVOLUMETET_AINSWORTH_EDGE_HDIV(oo + 1); dd++) {
          for (int ee = 0; ee < 6; ee++) {
            t_base(i) = t_base_v_e[ee](i);
            ++t_base;
            ++t_base_v_e[ee];
            t_diff_base(i, j) = t_diff_base_v_e[ee](i, j);
            ++t_diff_base;
            ++t_diff_base_v_e[ee];
          }
        }
        for (int dd = NBVOLUMETET_AINSWORTH_FACE_HDIV(oo);
             dd < NBVOLUMETET_AINSWORTH_FACE_HDIV(oo + 1); dd++) {
          for (int ff = 0; ff < 4; ff++) {
            t_base(i) = t_base_v_f[ff](i);
            ++t_base;
            ++t_base_v_f[ff];
            t_diff_base(i, j) = t_diff_base_v_f[ff](i, j);
            ++t_diff_base;
            ++t_diff_base_v_f[ff];
          }
        }
        for (int dd = NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo);
             dd < NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo + 1); dd++) {
          t_base(i) = (*t_base_v)(i);
          ++t_base;
          ++(*t_base_v);
          t_diff_base(i, j) = (*t_diff_base_v)(i, j);
          ++t_diff_base;
          ++(*t_diff_base_v);
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdivDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }
  int nb_gauss_pts = pts.size2();

  int volume_order = data.dataOnEntities[MBTET][0].getDataOrder();

  int p_f[4];
  double *phi_f[4];
  double *diff_phi_f[4];

  // Calculate base function on tet faces
  for (int ff = 0; ff != 4; ff++) {
    int face_order = data.dataOnEntities[MBTRI][ff].getDataOrder();
    int order = volume_order > face_order ? volume_order : face_order;
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
    p_f[ff] = order;
    phi_f[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
    diff_phi_f[ff] =
        &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    if (NBFACETRI_DEMKOWICZ_HDIV(order) == 0)
      continue;
    CHKERR Hdiv_Demkowicz_Face_MBTET_ON_FACE(
        &data.facesNodes(ff, 0), order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        phi_f[ff], diff_phi_f[ff], nb_gauss_pts, 4);
  }

  // Calculate base functions in tet interior
  if (NBVOLUMETET_DEMKOWICZ_HDIV(volume_order) > 0) {
    data.dataOnEntities[MBTET][0].getN(base).resize(
        nb_gauss_pts, 3 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);
    double *phi_v = &*data.dataOnEntities[MBTET][0].getN(base).data().begin();
    double *diff_phi_v =
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();
    CHKERR Hdiv_Demkowicz_Interior_MBTET(
        volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
        &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), p_f, phi_f,
        diff_phi_f, phi_v, diff_phi_v, nb_gauss_pts);
  }

  // Set size of face base correctly
  for (int ff = 0; ff != 4; ff++) {
    int face_order = data.dataOnEntities[MBTRI][ff].getDataOrder();
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    return getValueHdivAinsworthBase(pts);
  case DEMKOWICZ_JACOBI_BASE:
    return getValueHdivDemkowiczBase(pts);
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueHcurlAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  // edges
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {
    int sense[6], order[6];
    if (data.dataOnEntities[MBEDGE].size() != 6) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *hcurl_edge_n[6], *diff_hcurl_edge_n[6];
    for (int ee = 0; ee != 6; ee++) {
      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_AINSWORTH_HCURL(
          data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,
                                                        3 * nb_dofs, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                            9 * nb_dofs, false);
      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR Hcurl_Ainsworth_EdgeBaseFunctions_MBTET(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts, base_polynomials);
  } else {
    for (int ee = 0; ee != 6; ee++) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  // triangles
  if (data.spacesOnEntities[MBTRI].test(HCURL)) {
    int order[4];
    // faces
    if (data.dataOnEntities[MBTRI].size() != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    double *hcurl_base_n[4], *diff_hcurl_base_n[4];
    for (int ff = 0; ff != 4; ff++) {
      if (data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "data inconsistency");
      }
      order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      int nb_dofs = NBFACETRI_AINSWORTH_HCURL(order[ff]);
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,
                                                       3 * nb_dofs, false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,
                                                           9 * nb_dofs, false);
      hcurl_base_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      diff_hcurl_base_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if (data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    if (data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    CHKERR Hcurl_Ainsworth_FaceFunctions_MBTET(
        &*data.facesNodes.data().begin(), order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_base_n, diff_hcurl_base_n, nb_gauss_pts, base_polynomials);
  } else {
    for (int ff = 0; ff != 4; ff++) {
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                           false);
    }
  }

  if (data.spacesOnEntities[MBTET].test(HCURL)) {

    // volume
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_AINSWORTH_HCURL(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,
                                                    3 * nb_vol_dofs, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        9 * nb_vol_dofs, false);
    CHKERR Hcurl_Ainsworth_VolumeFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getDataOrder(),
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
        nb_gauss_pts, base_polynomials);

  } else {
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueHcurlDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }

  int nb_gauss_pts = pts.size2();

  // edges
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {
    int sense[6], order[6];
    if (data.dataOnEntities[MBEDGE].size() != 6) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong size of data structure, expected space for six edges "
               "but is %d",
               data.dataOnEntities[MBEDGE].size());
    }
    double *hcurl_edge_n[6], *diff_hcurl_edge_n[6];
    for (int ee = 0; ee != 6; ee++) {
      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "orintation of edges is not set");
      }
      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_DEMKOWICZ_HCURL(
          data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,
                                                        3 * nb_dofs, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                            9 * nb_dofs, false);
      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR Hcurl_Demkowicz_EdgeBaseFunctions_MBTET(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts);
  } else {
    // No DOFs on edges, resize base function matrices, indicating that no
    // dofs on them.
    for (int ee = 0; ee != 6; ee++) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  // triangles
  if (data.spacesOnEntities[MBTRI].test(HCURL)) {
    int order[4];
    // faces
    if (data.dataOnEntities[MBTRI].size() != 4) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "data structure for storing face h-curl base have wrong size "
               "should be four but is %d",
               data.dataOnEntities[MBTRI].size());
    }
    double *hcurl_base_n[4], *diff_hcurl_base_n[4];
    for (int ff = 0; ff != 4; ff++) {
      if (data.dataOnEntities[MBTRI][ff].getSense() == 0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "orintation of face is not set");
      }
      order[ff] = data.dataOnEntities[MBTRI][ff].getDataOrder();
      int nb_dofs = NBFACETRI_DEMKOWICZ_HCURL(order[ff]);
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts,
                                                       3 * nb_dofs, false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,
                                                           9 * nb_dofs, false);
      hcurl_base_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
      diff_hcurl_base_n[ff] =
          &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    }
    if (data.facesNodes.size1() != 4) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, should be four faces");
    }
    if (data.facesNodes.size2() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency, should be three nodes on face");
    }
    CHKERR Hcurl_Demkowicz_FaceBaseFunctions_MBTET(
        &*data.facesNodes.data().begin(), order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_base_n, diff_hcurl_base_n, nb_gauss_pts);
  } else {
    // No DOFs on faces, resize base function matrices, indicating that no
    // dofs on them.
    for (int ff = 0; ff != 4; ff++) {
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                           false);
    }
  }

  if (data.spacesOnEntities[MBTET].test(HCURL)) {
    // volume
    int order = data.dataOnEntities[MBTET][0].getDataOrder();
    int nb_vol_dofs = NBVOLUMETET_DEMKOWICZ_HCURL(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,
                                                    3 * nb_vol_dofs, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        9 * nb_vol_dofs, false);
    CHKERR Hcurl_Demkowicz_VolumeBaseFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getDataOrder(),
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin(),
        nb_gauss_pts);
  } else {
    // No DOFs on faces, resize base function matrices, indicating that no
    // dofs on them.
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHcurl(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    CHKERR getValueHcurlAinsworthBase(pts);
    break;
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueHcurlDemkowiczBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValue(MatrixDouble &pts,
                            boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;

  BaseFunctionUnknownInterface *iface;
  CHKERR ctx_ptr->query_interface(IDD_TET_BASE_FUNCTION, &iface);
  cTx = reinterpret_cast<EntPolynomialBaseCtx *>(iface);

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 3)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  if (cTx->bAse != AINSWORTH_BERNSTEIN_BEZIER_BASE) {
    const FieldApproximationBase base = cTx->bAse;
    DataForcesAndSourcesCore &data = cTx->dAta;
    if (cTx->copyNodeBase == LASTBASE) {
      data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 4,
                                                         false);
      CHKERR Tools::shapeFunMBTET(
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
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(4, 3, false);
    std::copy(Tools::diffShapeFunMBTET.begin(), Tools::diffShapeFunMBTET.end(),
              data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin());
  }

  switch (cTx->sPace) {
  case H1:
    CHKERR getValueH1(pts);
    break;
  case HDIV:
    CHKERR getValueHdiv(pts);
    break;
  case HCURL:
    CHKERR getValueHcurl(pts);
    break;
  case L2:
    CHKERR getValueL2(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown space");
  }

  MoFEMFunctionReturn(0);
}
