/** \file TriPolynomialBase.cpp
\brief Implementation of bases on triangle

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

TriPolynomialBase::TriPolynomialBase() {}
TriPolynomialBase::~TriPolynomialBase() {}

MoFEMErrorCode
TriPolynomialBase::query_interface(boost::typeindex::type_index type_index,
                                   UnknownInterface **iface) const {
  MoFEMFunctionBegin;
  *iface = const_cast<TriPolynomialBase *>(this);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueH1(MatrixDouble &pts) {
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

MoFEMErrorCode
TriPolynomialBase::getValueH1BernsteinBezierBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;
  EntitiesFieldData &data = cTx->dAta;
  const std::string &field_name = cTx->fieldName;
  int nb_gauss_pts = pts.size2();

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

  auto &vert_get_n = get_base(data.dataOnEntities[MBVERTEX][0]);
  auto &vert_get_diff_n = get_diff_base(data.dataOnEntities[MBVERTEX][0]);
  vert_get_n.resize(nb_gauss_pts, 3, false);
  vert_get_diff_n.resize(nb_gauss_pts, 6, false);

  if (data.dataOnEntities[MBVERTEX][0].getN(NOBASE).size1() !=
      (unsigned int)nb_gauss_pts)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[NOBASE]);
  auto &lambda = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);

  auto &vertex_alpha = get_alpha(data.dataOnEntities[MBVERTEX][0]);
  vertex_alpha.resize(3, 3, false);
  vertex_alpha.clear();
  for (int n = 0; n != 3; ++n)
    vertex_alpha(n, n) = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[n];

  CHKERR BernsteinBezier::baseFunctionsTri(
      1, lambda.size1(), vertex_alpha.size1(), &vertex_alpha(0, 0),
      &lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &vert_get_n(0, 0),
      &vert_get_diff_n(0, 0));
  for (int n = 0; n != 3; ++n) {
    const int f = boost::math::factorial<double>(
        data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[n]);
    for (int g = 0; g != nb_gauss_pts; ++g) {
      vert_get_n(g, n) *= f;
      for (int d = 0; d != 2; ++d)
        vert_get_diff_n(g, 2 * n + d) *= f;
    }
  }

  // edges
  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    if (data.dataOnEntities[MBEDGE].size() != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size of ent data");

    constexpr int edges_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
    for (int ee = 0; ee != 3; ++ee) {
      auto &ent_data = data.dataOnEntities[MBEDGE][ee];

      if (ent_data.getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Sense of the edge unknown");

      const int sense = ent_data.getSense();
      const int order = ent_data.getDataOrder();
      const int nb_dofs = NBEDGE_H1(ent_data.getDataOrder());

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
          get_diff_n.resize(nb_gauss_pts, 2 * nb_dofs, false);

          auto &edge_alpha = get_alpha(data.dataOnEntities[MBEDGE][ee]);
          edge_alpha.resize(nb_dofs, 3, false);
          CHKERR BernsteinBezier::generateIndicesEdgeTri(ee, order,
                                                         &edge_alpha(0, 0));
          if (sense == -1) {
            for (int i = 0; i != edge_alpha.size1(); ++i) {
              int a = edge_alpha(i, edges_nodes[ee][0]);
              edge_alpha(i, edges_nodes[ee][0]) =
                  edge_alpha(i, edges_nodes[ee][1]);
              edge_alpha(i, edges_nodes[ee][1]) = a;
            }
          }
          CHKERR BernsteinBezier::baseFunctionsTri(
              order, lambda.size1(), edge_alpha.size1(), &edge_alpha(0, 0),
              &lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &get_n(0, 0),
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
    for (int ee = 0; ee != 3; ++ee) {
      auto &ent_data = data.dataOnEntities[MBEDGE][ee];
      ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
      auto &get_n = get_base(ent_data);
      auto &get_diff_n = get_diff_base(ent_data);
      get_n.resize(nb_gauss_pts, 0, false);
      get_diff_n.resize(nb_gauss_pts, 0, false);
    }
  }

  if (data.spacesOnEntities[MBTRI].test(H1)) {
    if (data.dataOnEntities[MBTRI].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size ent of ent data");

    auto &ent_data = data.dataOnEntities[MBTRI][0];
    const int order = ent_data.getDataOrder();
    const int nb_dofs = NBFACETRI_H1(order);
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
      get_diff_n.resize(nb_gauss_pts, 2 * nb_dofs, false);
      if (nb_dofs) {
        auto &tri_alpha = get_alpha(ent_data);
        tri_alpha.resize(nb_dofs, 3, false);

        CHKERR BernsteinBezier::generateIndicesTriTri(order, &tri_alpha(0, 0));
        CHKERR BernsteinBezier::baseFunctionsTri(
            order, lambda.size1(), tri_alpha.size1(), &tri_alpha(0, 0),
            &lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &get_n(0, 0),
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
    auto &ent_data = data.dataOnEntities[MBTRI][0];
    ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
    auto &get_n = get_base(ent_data);
    auto &get_diff_n = get_diff_base(ent_data);
    get_n.resize(nb_gauss_pts, 0, false);
    get_diff_n.resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueH1AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (cTx->basePolynomialsType0 == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Polynomial type not set");
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    // edges
    if (data.dataOnEntities[MBEDGE].size() != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

    int sense[3], order[3];
    double *H1edgeN[3], *diffH1edgeN[3];
    for (int ee = 0; ee < 3; ee++) {
      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "sense of the edge unknown");

      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                        false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                            2 * nb_dofs, false);
      H1edgeN[ee] = &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diffH1edgeN[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR H1_EdgeShapeFunctions_MBTRI(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        H1edgeN, diffH1edgeN, nb_gauss_pts, base_polynomials);
  }

  if (data.spacesOnEntities[MBTRI].test(H1)) {
    // face
    if (data.dataOnEntities[MBTRI].size() != 1) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    int nb_dofs = NBFACETRI_H1(data.dataOnEntities[MBTRI][0].getDataOrder());
    data.dataOnEntities[MBTRI][0].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                    false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(nb_gauss_pts,
                                                        2 * nb_dofs, false);
    const int face_nodes[] = {0, 1, 2};
    CHKERR H1_FaceShapeFunctions_MBTRI(
        face_nodes, data.dataOnEntities[MBTRI][0].getDataOrder(),
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getDiffN(base).data().begin(),
        nb_gauss_pts, base_polynomials);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueL2(MatrixDouble &pts) {
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

MoFEMErrorCode TriPolynomialBase::getValueL2AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (cTx->basePolynomialsType0 == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Polynomial type not set");
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  data.dataOnEntities[MBTRI][0].getN(base).resize(
      nb_gauss_pts, NBFACETRI_L2(data.dataOnEntities[MBTRI][0].getDataOrder()),
      false);
  data.dataOnEntities[MBTRI][0].getDiffN(base).resize(
      nb_gauss_pts,
      2 * NBFACETRI_L2(data.dataOnEntities[MBTRI][0].getDataOrder()), false);

  CHKERR L2_Ainsworth_ShapeFunctions_MBTRI(
      data.dataOnEntities[MBTRI][0].getDataOrder(),
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
      &*data.dataOnEntities[MBTRI][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBTRI][0].getDiffN(base).data().begin(),
      nb_gauss_pts, base_polynomials);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TriPolynomialBase::getValueL2BernsteinBezierBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const std::string &field_name = cTx->fieldName;
  int nb_gauss_pts = pts.size2();

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

  if (data.spacesOnEntities[MBTRI].test(L2)) {
    if (data.dataOnEntities[MBTRI].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size ent of ent data");

    if (data.dataOnEntities[MBVERTEX][0].getN(NOBASE).size1() !=
        (unsigned int)nb_gauss_pts)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Base functions or nodes has wrong number of integration points "
               "for base %s",
               ApproximationBaseNames[NOBASE]);
    auto &lambda = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);

    auto &ent_data = data.dataOnEntities[MBTRI][0];
    const int order = ent_data.getDataOrder();
    const int nb_dofs = NBFACETRI_L2(order);

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
      get_diff_n.resize(nb_gauss_pts, 2 * nb_dofs, false);
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

          if (nb_dofs != 3 + 3 * NBEDGE_H1(order) + NBFACETRI_H1(order))
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "Inconsistent number of DOFs");

          auto &tri_alpha = get_alpha(ent_data);
          tri_alpha.resize(nb_dofs, 3, false);

          CHKERR BernsteinBezier::generateIndicesVertexTri(order,
                                                           &tri_alpha(0, 0));

          if (order > 1) {
            std::array<int, 3> face_n{order, order, order};
            std::array<int *, 3> face_edge_ptr{
                &tri_alpha(3, 0), &tri_alpha(3 + NBEDGE_H1(order), 0),
                &tri_alpha(3 + 2 * NBEDGE_H1(order), 0)};
            CHKERR BernsteinBezier::generateIndicesEdgeTri(
                face_n.data(), face_edge_ptr.data());
            if (order > 2)
              CHKERR BernsteinBezier::generateIndicesTriTri(
                  order, &tri_alpha(3 + 3 * NBEDGE_H1(order), 0));
          }
          CHKERR BernsteinBezier::baseFunctionsTri(
              order, lambda.size1(), tri_alpha.size1(), &tri_alpha(0, 0),
              &lambda(0, 0), Tools::diffShapeFunMBTRI.data(), &get_n(0, 0),
              &get_diff_n(0, 0));
        }

        get_alpha_by_order_ptr(ent_data, order) =
            get_alpha_by_name_ptr(ent_data, field_name);
        get_base_by_order_ptr(ent_data, order) =
            get_base_by_name_ptr(ent_data, field_name);
        get_diff_base_by_order_ptr(ent_data, order) =
            get_diff_base_by_name_ptr(ent_data, field_name);
      }
    }
  } else {
    auto &ent_data = data.dataOnEntities[MBTRI][0];
    ent_data.getBBAlphaIndicesSharedPtr(field_name).reset();
    auto &get_n = get_base(ent_data);
    auto &get_diff_n = get_diff_base(ent_data);
    get_n.resize(nb_gauss_pts, 0, false);
    get_diff_n.resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueHdivAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (cTx->basePolynomialsType0 == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Polynomial type not set");
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  double *PHI_f_e[3];
  double *PHI_f;

  N_face_edge.resize(1, 3, false);
  N_face_bubble.resize(1, false);
  int face_order = data.dataOnEntities[MBTRI][0].getDataOrder();
  // three edges on face
  for (int ee = 0; ee < 3; ee++) {
    N_face_edge(0, ee).resize(
        nb_gauss_pts, 3 * NBFACETRI_AINSWORTH_EDGE_HDIV(face_order), false);
    PHI_f_e[ee] = &((N_face_edge(0, ee))(0, 0));
  }
  N_face_bubble[0].resize(nb_gauss_pts,
                          3 * NBFACETRI_AINSWORTH_FACE_HDIV(face_order), false);
  PHI_f = &*(N_face_bubble[0].data().begin());

  int face_nodes[3] = {0, 1, 2};
  CHKERR Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET_ON_FACE(
      face_nodes, face_order,
      &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0), NULL, PHI_f_e, NULL,
      nb_gauss_pts, 3, base_polynomials);
  CHKERR Hdiv_Ainsworth_FaceBubbleShapeFunctions_ON_FACE(
      face_nodes, face_order,
      &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0), NULL, PHI_f, NULL,
      nb_gauss_pts, 3, base_polynomials);

  // set shape functions into data structure
  if (data.dataOnEntities[MBTRI].size() != 1) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  data.dataOnEntities[MBTRI][0].getN(base).resize(
      nb_gauss_pts, 3 * NBFACETRI_AINSWORTH_HDIV(face_order), false);
  int col = 0;
  for (int oo = 0; oo < face_order; oo++) {
    for (int ee = 0; ee < 3; ee++) {
      for (int dd = 3 * NBFACETRI_AINSWORTH_EDGE_HDIV(oo);
           dd < 3 * NBFACETRI_AINSWORTH_EDGE_HDIV(oo + 1); dd++, col++) {
        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          data.dataOnEntities[MBTRI][0].getN(base)(gg, col) =
              N_face_edge(0, ee)(gg, dd);
        }
      }
    }
    for (int dd = 3 * NBFACETRI_AINSWORTH_FACE_HDIV(oo);
         dd < 3 * NBFACETRI_AINSWORTH_FACE_HDIV(oo + 1); dd++, col++) {
      for (int gg = 0; gg < nb_gauss_pts; gg++) {
        data.dataOnEntities[MBTRI][0].getN(base)(gg, col) =
            N_face_bubble[0](gg, dd);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueHdivDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  // set shape functions into data structure
  if (data.dataOnEntities[MBTRI].size() != 1) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  const FieldApproximationBase base = cTx->bAse;
  if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }
  int order = data.dataOnEntities[MBTRI][0].getDataOrder();
  int nb_gauss_pts = pts.size2();
  data.dataOnEntities[MBTRI][0].getN(base).resize(
      nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
  double *phi_f = &*data.dataOnEntities[MBTRI][0].getN(base).data().begin();
  if (NBFACETRI_DEMKOWICZ_HDIV(order) == 0)
    MoFEMFunctionReturnHot(0);
  int face_nodes[3] = {0, 1, 2};  
  CHKERR Hdiv_Demkowicz_Face_MBTET_ON_FACE(
      face_nodes, order,
      &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
      &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(), phi_f,
      NULL, nb_gauss_pts, 3);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    return getValueHdivAinsworthBase(pts);
  case DEMKOWICZ_JACOBI_BASE:
    return getValueHdivDemkowiczBase(pts);
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TriPolynomialBase::getValueHcurlAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (data.dataOnEntities[MBTRI].size() != 1)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  // Calculation H-curl on triangle faces
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {
    if (data.dataOnEntities[MBEDGE].size() != 3) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    int sense[3], order[3];
    double *hcurl_edge_n[3];
    double *diff_hcurl_edge_n[3];
    for (int ee = 0; ee < 3; ee++) {
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
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(
          nb_gauss_pts, 2 * 3 * nb_dofs, false);
      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_FACE(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts, base_polynomials);
  } else {
    for (int ee = 0; ee < 3; ee++) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  if (data.spacesOnEntities[MBTRI].test(HCURL)) {

    // cerr << data.dataOnEntities[MBVERTEX][0].getN(base) << endl;
    // cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << endl;
    //
    // face
    if (data.dataOnEntities[MBTRI].size() != 1) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    int order = data.dataOnEntities[MBTRI][0].getDataOrder();
    int nb_dofs = NBFACETRI_AINSWORTH_HCURL(order);
    data.dataOnEntities[MBTRI][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                    false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(nb_gauss_pts,
                                                        3 * 2 * nb_dofs, false);
    // cerr << data.dataOnEntities[MBVERTEX][0].getDiffN(base) << endl;
    int face_nodes[] = {0, 1, 2};
    CHKERR Hcurl_Ainsworth_FaceFunctions_MBTET_ON_FACE(
        face_nodes, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getDiffN(base).data().begin(),
        nb_gauss_pts, base_polynomials);
    // cerr << data.dataOnEntities[MBTRI][0].getN(base) << endl;
  } else {
    data.dataOnEntities[MBTRI][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TriPolynomialBase::getValueHcurlDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();

  // Calculation H-curl on triangle faces
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {

    if (data.dataOnEntities[MBEDGE].size() != 3)
      SETERRQ1(
          PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "wrong number of data structures on edges, should be three but is %d",
          data.dataOnEntities[MBEDGE].size());

    int sense[3], order[3];
    double *hcurl_edge_n[3];
    double *diff_hcurl_edge_n[3];

    for (int ee = 0; ee != 3; ++ee) {

      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "orientation (sense) of edge is not set");

      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
      int nb_dofs = NBEDGE_DEMKOWICZ_HCURL(
          data.dataOnEntities[MBEDGE][ee].getDataOrder());
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,
                                                        3 * nb_dofs, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(
          nb_gauss_pts, 2 * 3 * nb_dofs, false);

      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }

    CHKERR Hcurl_Demkowicz_EdgeBaseFunctions_MBTRI(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts);

  } else {

    // No DOFs on faces, resize base function matrices, indicating that no
    // dofs on them.
    for (int ee = 0; ee != 3; ++ee) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  if (data.spacesOnEntities[MBTRI].test(HCURL)) {

    // face
    if (data.dataOnEntities[MBTRI].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No data struture to keep base functions on face");

    int order = data.dataOnEntities[MBTRI][0].getDataOrder();
    int nb_dofs = NBFACETRI_DEMKOWICZ_HCURL(order);
    data.dataOnEntities[MBTRI][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                    false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(nb_gauss_pts,
                                                        3 * 2 * nb_dofs, false);

    int face_nodes[] = {0, 1, 2};
    CHKERR Hcurl_Demkowicz_FaceBaseFunctions_MBTRI(
        face_nodes, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBTRI][0].getDiffN(base).data().begin(),
        nb_gauss_pts);

  } else {

    // No DOFs on faces, resize base function matrices, indicating that no
    // dofs on them.
    data.dataOnEntities[MBTRI][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBTRI][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TriPolynomialBase::getValueHcurl(MatrixDouble &pts) {
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
TriPolynomialBase::getValue(MatrixDouble &pts,
                            boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;

  cTx = ctx_ptr->getInterface<EntPolynomialBaseCtx>();

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts) {
    MoFEMFunctionReturnHot(0);
  }

  if (pts.size1() < 2)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  const FieldApproximationBase base = cTx->bAse;
  EntitiesFieldData &data = cTx->dAta;

  if (base != AINSWORTH_BERNSTEIN_BEZIER_BASE) {
    if (cTx->copyNodeBase == LASTBASE) {
      data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 3,
                                                         false);
      CHKERR ShapeMBTRI(
          &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
          &pts(0, 0), &pts(1, 0), nb_gauss_pts);
      data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(3, 2, false);
      CHKERR ShapeDiffMBTRI(
          &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin());
    } else {
      data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
          data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
      data.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(base) =
          data.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(cTx->copyNodeBase);
    }
    if (data.dataOnEntities[MBVERTEX][0].getN(base).size1() !=
        (unsigned int)nb_gauss_pts) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Base functions or nodes has wrong number of integration points "
               "for base %s",
               ApproximationBaseNames[base]);
    }
  }
  auto set_node_derivative_for_all_gauss_pts = [&]() {
    MoFEMFunctionBegin;
    // In linear geometry derivatives are constant,
    // this in expense of efficiency makes implementation
    // consistent between vertices and other types of entities
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts, 6,
                                                           false);
    for (int gg = 0; gg != nb_gauss_pts; ++gg)
      std::copy(Tools::diffShapeFunMBTRI.begin(),
                Tools::diffShapeFunMBTRI.end(),
                &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg, 0));
    MoFEMFunctionReturn(0);
  };

  CHKERR set_node_derivative_for_all_gauss_pts();

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
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  MoFEMFunctionReturn(0);
}
