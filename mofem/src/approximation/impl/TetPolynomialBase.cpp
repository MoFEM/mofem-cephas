/** \file TetPolynomialBase.cpp
\brief Implementation of hierarchical bases on tetrahedral

A l2, h1, h-div and h-curl spaces are implemented.

*/

using namespace MoFEM;

struct TetBaseCache {

  struct BaseCacheItem {
    int order;
    int nb_gauss_pts;
    mutable MatrixDouble N;
    mutable MatrixDouble diffN;
  };

  struct HDivBaseCacheItem {

    int order;
    int nb_gauss_pts;

    // Number of permeations for tetrahedron
    // That is P(3, 4) = 24

    int n0;
    int n1;
    int n2;

    mutable MatrixDouble N;
    mutable MatrixDouble diffN;
  };

  using BaseCacheMI = boost::multi_index_container<
      BaseCacheItem,
      boost::multi_index::indexed_by<

          boost::multi_index::hashed_unique<

              composite_key<

                  BaseCacheItem,
                  member<BaseCacheItem, int, &BaseCacheItem::order>,
                  member<BaseCacheItem, int, &BaseCacheItem::nb_gauss_pts>>>>

      >;

  using HDivBaseFaceCacheMI = boost::multi_index_container<
      HDivBaseCacheItem,
      boost::multi_index::indexed_by<

          boost::multi_index::hashed_unique<

              composite_key<

                  HDivBaseCacheItem,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::order>,
                  member<HDivBaseCacheItem, int,
                         &HDivBaseCacheItem::nb_gauss_pts>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n0>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n1>,
                  member<HDivBaseCacheItem, int, &HDivBaseCacheItem::n2>>>>

      >;

  static std::map<const void *, BaseCacheMI> hdivBaseInteriorDemkowicz;
  static std::map<const void *, BaseCacheMI> hdivBrokenBaseInteriorDemkowicz;
  static std::map<const void *, HDivBaseFaceCacheMI> hDivBaseFaceDemkowicz;

  static std::map<const void *, BaseCacheMI> hdivBrokenBaseInteriorAinsworth;
};

std::map<const void *, TetBaseCache::BaseCacheMI>
    TetBaseCache::hdivBaseInteriorDemkowicz;
std::map<const void *, TetBaseCache::BaseCacheMI>
    TetBaseCache::hdivBrokenBaseInteriorDemkowicz;
std::map<const void *, TetBaseCache::HDivBaseFaceCacheMI>
    TetBaseCache::hDivBaseFaceDemkowicz;

std::map<const void *, TetBaseCache::BaseCacheMI>
    TetBaseCache::hdivBrokenBaseInteriorAinsworth;

MoFEMErrorCode
TetPolynomialBase::query_interface(boost::typeindex::type_index type_index,
                                   UnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = const_cast<TetPolynomialBase *>(this);
  MoFEMFunctionReturnHot(0);
}

TetPolynomialBase::TetPolynomialBase(const void *ptr) : vPtr(ptr) {}

TetPolynomialBase::~TetPolynomialBase() {
  if (vPtr) {
    if (TetBaseCache::hdivBaseInteriorDemkowicz.find(vPtr) !=
        TetBaseCache::hdivBaseInteriorDemkowicz.end())
      TetBaseCache::hdivBaseInteriorDemkowicz.erase(vPtr);
    if (TetBaseCache::hDivBaseFaceDemkowicz.find(vPtr) !=
        TetBaseCache::hDivBaseFaceDemkowicz.end())
      TetBaseCache::hDivBaseFaceDemkowicz.erase(vPtr);
  }
}

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

  EntitiesFieldData &data = cTx->dAta;
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
      order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
      int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getOrder());
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
      int nb_dofs = NBFACETRI_H1(data.dataOnEntities[MBTRI][ff].getOrder());
      order[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
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
    int order = data.dataOnEntities[MBTET][0].getOrder();
    int nb_vol_dofs = NBVOLUMETET_H1(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, nb_vol_dofs,
                                                    false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        3 * nb_vol_dofs, false);
    CHKERR H1_VolumeShapeFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getOrder(),
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

  EntitiesFieldData &data = cTx->dAta;
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
      const int order = ent_data.getOrder();
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
      const int order = ent_data.getOrder();
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
    const int order = ent_data.getOrder();
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
  case DEMKOWICZ_JACOBI_BASE:
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

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  data.dataOnEntities[MBTET][0].getN(base).resize(
      nb_gauss_pts, NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getOrder()),
      false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(
      nb_gauss_pts,
      3 * NBVOLUMETET_L2(data.dataOnEntities[MBTET][0].getOrder()), false);

  CHKERR L2_Ainsworth_ShapeFunctions_MBTET(
      data.dataOnEntities[MBTET][0].getOrder(),
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

  EntitiesFieldData &data = cTx->dAta;
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
    const int order = ent_data.getOrder();
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
          if (order > 1) {
            std::array<int, 6> edge_n{order, order, order, order, order, order};
            std::array<int *, 6> tet_edge_ptr{
                &tet_alpha(4, 0),
                &tet_alpha(4 + 1 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 2 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 3 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 4 * NBEDGE_H1(order), 0),
                &tet_alpha(4 + 5 * NBEDGE_H1(order), 0)};
            CHKERR BernsteinBezier::generateIndicesEdgeTet(edge_n.data(),
                                                           tet_edge_ptr.data());
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
              if (order > 3)
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

MoFEMErrorCode TetPolynomialBase::getValueHdivAinsworthBaseImpl(

    MatrixDouble &pts,

    MatrixDouble &shape_functions, MatrixDouble &diff_shape_functions,

    int volume_order, std::array<int, 4> &faces_order,
    std::array<int, 3 * 4> &faces_nodes

) {
  MoFEMFunctionBegin;

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

  for (int ff = 0; ff != 4; ++ff) {
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
      faces_nodes.data(), faces_order.data(), &*shape_functions.data().begin(),
      &*diff_shape_functions.data().begin(), phi_f_e, diff_phi_f_e,
      nb_gauss_pts, base_polynomials);

  CHKERR Hdiv_Ainsworth_FaceBubbleShapeFunctions(
      faces_nodes.data(), faces_order.data(), &*shape_functions.data().begin(),
      &*diff_shape_functions.data().begin(), phi_f, diff_phi_f, nb_gauss_pts,
      base_polynomials);

  // volume shape functions

  double *phi_v_e[6];
  double *phi_v_f[4];
  double *phi_v;
  double *diff_phi_v_e[6];
  double *diff_phi_v_f[4];
  double *diff_phi_v;

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
  if (NBVOLUMETET_AINSWORTH_EDGE_HDIV(volume_order))
    CHKERR Hdiv_Ainsworth_EdgeBasedVolumeShapeFunctions_MBTET(
        volume_order, &*shape_functions.data().begin(),
        &*diff_shape_functions.data().begin(), phi_v_e, diff_phi_v_e,
        nb_gauss_pts, base_polynomials);

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
  if (NBVOLUMETET_AINSWORTH_FACE_HDIV(volume_order))
    CHKERR Hdiv_Ainsworth_FaceBasedVolumeShapeFunctions_MBTET(
        volume_order, &*shape_functions.data().begin(),
        &*diff_shape_functions.data().begin(), phi_v_f, diff_phi_v_f,
        nb_gauss_pts, base_polynomials);

  N_volume_bubble.resize(
      nb_gauss_pts, 3 * NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order), false);
  diffN_volume_bubble.resize(
      nb_gauss_pts, 9 * NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order), false);
  phi_v = &*(N_volume_bubble.data().begin());
  diff_phi_v = &*(diffN_volume_bubble.data().begin());
  if (NBVOLUMETET_AINSWORTH_VOLUME_HDIV(volume_order))
    CHKERR Hdiv_Ainsworth_VolumeBubbleShapeFunctions_MBTET(
        volume_order, &*shape_functions.data().begin(),
        &*diff_shape_functions.data().begin(), phi_v, diff_phi_v, nb_gauss_pts,
        base_polynomials);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdivAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  std::array<int, 4> faces_order;
  std::array<int, 4 * 3> faces_nodes;

  FieldApproximationBase base = cTx->bAse;
  EntitiesFieldData &data = cTx->dAta;

  int volume_order = data.dataOnEntities[MBTET][0].getOrder();
  std::copy(data.facesNodes.data().begin(), data.facesNodes.data().end(),
            faces_nodes.begin());
  for (int ff = 0; ff != 4; ff++) {
    faces_order[ff] = cTx->dAta.dataOnEntities[MBTRI][ff].getOrder();
  }

  CHKERR getValueHdivAinsworthBaseImpl(
      pts, data.dataOnEntities[MBVERTEX][0].getN(base),
      data.dataOnEntities[MBVERTEX][0].getDiffN(base), volume_order,
      faces_order, faces_nodes);

  // Set shape functions into data structure Shape functions hast to be put
  // in arrays in order which guarantee hierarchical series of degrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  int nb_gauss_pts = pts.size2();

  // faces
  if (data.dataOnEntities[MBTRI].size() != 4) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }

  // face-face
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_f_f[] = {
      getFTensor1FromPtr<3>(&*(N_face_bubble[0].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[1].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[2].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[3].data().begin()))};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_f_f[] = {
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[0].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[1].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[2].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[3].data().begin()))};
  // face-edge
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_f_e[] = {
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 2).data().begin()))};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_f_e[] = {
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 2).data().begin()))};

  for (int ff = 0; ff != 4; ff++) {
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_AINSWORTH_HDIV(faces_order[ff]), false);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_AINSWORTH_HDIV(faces_order[ff]), false);
    if (NBFACETRI_AINSWORTH_HDIV(faces_order[ff]) == 0)
      continue;
    double *base_ptr =
        &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
    double *diff_base_ptr =
        &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
    auto t_base = getFTensor1FromPtr<3>(base_ptr);
    auto t_diff_base = getFTensor2HVecFromPtr<3, 3>(diff_base_ptr);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int oo = 0; oo != faces_order[ff]; oo++) {

        // face-edge
        for (int dd = NBFACETRI_AINSWORTH_EDGE_HDIV(oo);
             dd != NBFACETRI_AINSWORTH_EDGE_HDIV(oo + 1); dd++) {
          for (int ee = 0; ee != 3; ++ee) {
            t_base(i) = t_base_f_e[ff * 3 + ee](i);
            ++t_base;
            ++t_base_f_e[ff * 3 + ee];
          }
          for (int ee = 0; ee != 3; ++ee) {
            t_diff_base(i, j) = t_diff_base_f_e[ff * 3 + ee](i, j);
            ++t_diff_base;
            ++t_diff_base_f_e[ff * 3 + ee];
          }
        }

        // face-face
        for (int dd = NBFACETRI_AINSWORTH_FACE_HDIV(oo);
             dd != NBFACETRI_AINSWORTH_FACE_HDIV(oo + 1); dd++) {
          t_base(i) = t_base_f_f[ff](i);
          ++t_base;
          ++t_base_f_f[ff];
          t_diff_base(i, j) = t_diff_base_f_f[ff](i, j);
          ++t_diff_base;
          ++t_diff_base_f_f[ff];
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
    double *diff_base_ptr =
        &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();
    auto t_base = getFTensor1FromPtr<3>(base_ptr);
    auto t_diff_base = getFTensor2HVecFromPtr<3, 3>(diff_base_ptr);

    // volume-edge
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v_e[] = {
        getFTensor1FromPtr<3>(&*N_volume_edge[0].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_edge[1].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_edge[2].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_edge[3].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_edge[4].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_edge[5].data().begin())};
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v_e[] = {
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[0].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[1].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[2].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[3].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[4].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[5].data().begin())};

    // volume-faces
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v_f[] = {
        getFTensor1FromPtr<3>(&*N_volume_face[0].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_face[1].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_face[2].data().begin()),
        getFTensor1FromPtr<3>(&*N_volume_face[3].data().begin())};
    FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v_f[] = {
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[0].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[1].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[2].data().begin()),
        getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[3].data().begin())};

    // volume-bubble
    base_ptr = &*(N_volume_bubble.data().begin());
    diff_base_ptr = &*(diffN_volume_bubble.data().begin());
    auto t_base_v = getFTensor1FromPtr<3>(base_ptr);
    auto t_diff_base_v = getFTensor2HVecFromPtr<3, 3>(diff_base_ptr);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      for (int oo = 0; oo < volume_order; oo++) {

        // volume-edge
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

        // volume-face
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

        // volume-bubble
        for (int dd = NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo);
             dd < NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo + 1); dd++) {
          t_base(i) = t_base_v(i);
          ++t_base;
          ++t_base_v;
          t_diff_base(i, j) = t_diff_base_v(i, j);
          ++t_diff_base;
          ++t_diff_base_v;
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueHdivAinsworthBrokenBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  // Set shape functions into data structure Shape functions has to be put
  // in arrays in order which guarantee hierarchical series of degrees of
  // freedom, i.e. in other words dofs form sub-entities has to be group
  // by order.

  FieldApproximationBase base = cTx->bAse;
  EntitiesFieldData &data = cTx->dAta;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  int volume_order = data.dataOnEntities[MBTET][0].getOrder();
  int nb_gauss_pts = pts.size2();
  int nb_dofs_face = NBFACETRI_AINSWORTH_HDIV(volume_order);
  int nb_dofs_volume = NBVOLUMETET_AINSWORTH_HDIV(volume_order);
  int nb_dofs = 4 * nb_dofs_face + nb_dofs_volume;
  data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                  false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts, 9 * nb_dofs,
                                                      false);
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);


  auto get_interior_cache = [this]() -> TetBaseCache::BaseCacheMI * {
    if (vPtr) {
      auto it = TetBaseCache::hdivBrokenBaseInteriorAinsworth.find(vPtr);
      if (it != TetBaseCache::hdivBrokenBaseInteriorAinsworth.end()) {
        return &it->second;
      }
    }
    return nullptr;
  };

  auto interior_cache_ptr = get_interior_cache();

  if (interior_cache_ptr) {
    auto it =
        interior_cache_ptr->find(boost::make_tuple(volume_order, nb_gauss_pts));
    if (it != interior_cache_ptr->end()) {
      noalias(data.dataOnEntities[MBTET][0].getN(base)) = it->N;
      noalias(data.dataOnEntities[MBTET][0].getDiffN(base)) = it->diffN;
      MoFEMFunctionBeginHot(0);
    } 
  }

  std::array<int, 4 * 3> faces_nodes{
      0, 1, 3, // face 0

      1, 2, 3, // face 1

      0, 2, 3, // face 2

      0, 1, 2 // face 3

  };
  std::array<int, 4> faces_order{volume_order, volume_order, volume_order,
                                 volume_order};
  CHKERR getValueHdivAinsworthBaseImpl(
      pts, data.dataOnEntities[MBVERTEX][0].getN(base),
      data.dataOnEntities[MBVERTEX][0].getDiffN(base), volume_order,
      faces_order, faces_nodes);

  auto *base_ptr = &*data.dataOnEntities[MBTET][0].getN(base).data().begin();
  auto *diff_base_ptr =
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();
  auto t_base = getFTensor1FromPtr<3>(base_ptr);
  auto t_diff_base = getFTensor2HVecFromPtr<3, 3>(diff_base_ptr);

  // face-edge
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_f_e[] = {
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(0, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(1, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(2, 2).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 0).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 1).data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_edge(3, 2).data().begin()))};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_f_e[] = {
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(0, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(1, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(2, 2).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 0).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 1).data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_edge(3, 2).data().begin()))};

  // face-face
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_f_f[] = {
      getFTensor1FromPtr<3>(&*(N_face_bubble[0].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[1].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[2].data().begin())),
      getFTensor1FromPtr<3>(&*(N_face_bubble[3].data().begin()))};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_f_f[] = {
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[0].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[1].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[2].data().begin())),
      getFTensor2HVecFromPtr<3, 3>(&*(diffN_face_bubble[3].data().begin()))};

  // volume-edge
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v_e[] = {
      getFTensor1FromPtr<3>(&*N_volume_edge[0].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_edge[1].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_edge[2].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_edge[3].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_edge[4].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_edge[5].data().begin())};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v_e[] = {
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[0].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[1].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[2].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[3].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[4].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_edge[5].data().begin())};

  // volume-faces
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v_f[] = {
      getFTensor1FromPtr<3>(&*N_volume_face[0].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_face[1].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_face[2].data().begin()),
      getFTensor1FromPtr<3>(&*N_volume_face[3].data().begin())};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v_f[] = {
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[0].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[1].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[2].data().begin()),
      getFTensor2HVecFromPtr<3, 3>(&*diffN_volume_face[3].data().begin())};

  // volume-bubble
  auto *base_vol_ptr = &*(N_volume_bubble.data().begin());
  auto *diff_base_vol_ptr = &*(diffN_volume_bubble.data().begin());
  auto t_base_v = getFTensor1FromPtr<3>(base_vol_ptr);
  auto t_diff_base_v = getFTensor2HVecFromPtr<3, 3>(diff_base_vol_ptr);

  // int count_dofs = 0;

  for (int gg = 0; gg != nb_gauss_pts; gg++) {
    for (int oo = 0; oo < volume_order; oo++) {

      // faces-edge (((P) > 0) ? (P) : 0)
      for (int dd = NBFACETRI_AINSWORTH_EDGE_HDIV(oo);
           dd != NBFACETRI_AINSWORTH_EDGE_HDIV(oo + 1); dd++) {
        for (auto ff = 0; ff != 4; ++ff) {
          for (int ee = 0; ee != 3; ++ee) {
            t_base(i) = t_base_f_e[ff * 3 + ee](i);
            ++t_base;
            ++t_base_f_e[ff * 3 + ee];
            // ++count_dofs;
          }
          for (int ee = 0; ee != 3; ++ee) {
            t_diff_base(i, j) = t_diff_base_f_e[ff * 3 + ee](i, j);
            ++t_diff_base;
            ++t_diff_base_f_e[ff * 3 + ee];
          }
        }
      }

      // face-face (P - 1) * (P - 2) / 2
      for (int dd = NBFACETRI_AINSWORTH_FACE_HDIV(oo);
           dd != NBFACETRI_AINSWORTH_FACE_HDIV(oo + 1); dd++) {
        for (auto ff = 0; ff != 4; ++ff) {
          t_base(i) = t_base_f_f[ff](i);
          ++t_base;
          ++t_base_f_f[ff];
          t_diff_base(i, j) = t_diff_base_f_f[ff](i, j);
          ++t_diff_base;
          ++t_diff_base_f_f[ff];
          // ++count_dofs;
        }
      }

      // volume-edge (P - 1)
      for (int dd = NBVOLUMETET_AINSWORTH_EDGE_HDIV(oo);
           dd != NBVOLUMETET_AINSWORTH_EDGE_HDIV(oo + 1); dd++) {
        for (int ee = 0; ee < 6; ee++) {
          t_base(i) = t_base_v_e[ee](i);
          ++t_base;
          ++t_base_v_e[ee];
          t_diff_base(i, j) = t_diff_base_v_e[ee](i, j);
          ++t_diff_base;
          ++t_diff_base_v_e[ee];
          // ++count_dofs;
        }
      }
      // volume-face (P - 1) * (P - 2)
      for (int dd = NBVOLUMETET_AINSWORTH_FACE_HDIV(oo);
           dd < NBVOLUMETET_AINSWORTH_FACE_HDIV(oo + 1); dd++) {
        for (int ff = 0; ff < 4; ff++) {
          t_base(i) = t_base_v_f[ff](i);
          ++t_base;
          ++t_base_v_f[ff];
          t_diff_base(i, j) = t_diff_base_v_f[ff](i, j);
          ++t_diff_base;
          ++t_diff_base_v_f[ff];
          // ++count_dofs;
        }
      }
      // volume-bubble (P - 3) * (P - 2) * (P - 1) / 2
      for (int dd = NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo);
           dd < NBVOLUMETET_AINSWORTH_VOLUME_HDIV(oo + 1); dd++) {
        t_base(i) = t_base_v(i);
        ++t_base;
        ++t_base_v;
        t_diff_base(i, j) = t_diff_base_v(i, j);
        ++t_diff_base;
        ++t_diff_base_v;
        // ++count_dofs;
      }
    }
  }

// #ifdef NDEBUG
//   if (nb_dofs != count_dofs / nb_gauss_pts)
//     SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//              "Number of dofs %d is different than expected %d", count_dofs,
//              nb_dofs);
// #endif // NDEBUG

  if (interior_cache_ptr) {
    auto p = interior_cache_ptr->emplace(
        TetBaseCache::BaseCacheItem{volume_order, nb_gauss_pts});
    p.first->N = data.dataOnEntities[MBTET][0].getN(base);
    p.first->diffN = data.dataOnEntities[MBTET][0].getDiffN(base);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdivDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }
  int nb_gauss_pts = pts.size2();

  int volume_order = data.dataOnEntities[MBTET][0].getOrder();

  int p_f[4];
  double *phi_f[4];
  double *diff_phi_f[4];

  auto get_face_cache_ptr = [this]() -> TetBaseCache::HDivBaseFaceCacheMI * {
    if (vPtr) {
      auto it = TetBaseCache::hDivBaseFaceDemkowicz.find(vPtr);
      if (it != TetBaseCache::hDivBaseFaceDemkowicz.end()) {
        return &it->second;
      }
    }
    return nullptr;
  };

  auto face_cache_ptr = get_face_cache_ptr();

  // Calculate base function on tet faces
  for (int ff = 0; ff != 4; ff++) {
    int face_order = data.dataOnEntities[MBTRI][ff].getOrder();
    int order = volume_order > face_order ? volume_order : face_order;
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
    if (NBFACETRI_DEMKOWICZ_HDIV(order) == 0)
      continue;

    if (face_cache_ptr) {
      auto it = face_cache_ptr->find(boost::make_tuple(

          face_order, nb_gauss_pts,

          data.facesNodes(ff, 0), data.facesNodes(ff, 1), data.facesNodes(ff, 2)

              ));
      if (it != face_cache_ptr->end()) {
        noalias(data.dataOnEntities[MBTRI][ff].getN(base)) = it->N;
        noalias(data.dataOnEntities[MBTRI][ff].getDiffN(base)) = it->diffN;
        continue;
      }
    }

    p_f[ff] = order;
    phi_f[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
    diff_phi_f[ff] =
        &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();

    CHKERR Hdiv_Demkowicz_Face_MBTET_ON_FACE(
        &data.facesNodes(ff, 0), order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        phi_f[ff], diff_phi_f[ff], nb_gauss_pts, 4);
    if (face_cache_ptr) {
      auto p = face_cache_ptr->emplace(TetBaseCache::HDivBaseCacheItem{
          face_order, nb_gauss_pts, data.facesNodes(ff, 0),
          data.facesNodes(ff, 1), data.facesNodes(ff, 2)});
      p.first->N = data.dataOnEntities[MBTRI][ff].getN(base);
      p.first->diffN = data.dataOnEntities[MBTRI][ff].getDiffN(base);
    }
  }

  auto get_interior_cache = [this]() -> TetBaseCache::BaseCacheMI * {
    if (vPtr) {
      auto it = TetBaseCache::hdivBaseInteriorDemkowicz.find(vPtr);
      if (it != TetBaseCache::hdivBaseInteriorDemkowicz.end()) {
        return &it->second;
      }
    }
    return nullptr;
  };

  auto interior_cache_ptr = get_interior_cache();

  // Calculate base functions in tet interior
  if (NBVOLUMETET_DEMKOWICZ_HDIV(volume_order) > 0) {
    data.dataOnEntities[MBTET][0].getN(base).resize(
        nb_gauss_pts, 3 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);

    for (int v = 0; v != 1; ++v) {
      if (interior_cache_ptr) {
        auto it = interior_cache_ptr->find(
            boost::make_tuple(volume_order, nb_gauss_pts));
        if (it != interior_cache_ptr->end()) {
          noalias(data.dataOnEntities[MBTET][0].getN(base)) = it->N;
          noalias(data.dataOnEntities[MBTET][0].getDiffN(base)) = it->diffN;
          continue;
        }
      }

      double *phi_v = &*data.dataOnEntities[MBTET][0].getN(base).data().begin();
      double *diff_phi_v =
          &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();

      CHKERR Hdiv_Demkowicz_Interior_MBTET(
          volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
          &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), p_f, phi_f,
          diff_phi_f, phi_v, diff_phi_v, nb_gauss_pts);
      if (interior_cache_ptr) {
        auto p = interior_cache_ptr->emplace(
            TetBaseCache::BaseCacheItem{volume_order, nb_gauss_pts});
        p.first->N = data.dataOnEntities[MBTET][0].getN(base);
        p.first->diffN = data.dataOnEntities[MBTET][0].getDiffN(base);
      }
    }
  }

  // Set size of face base correctly
  for (int ff = 0; ff != 4; ff++) {
    int face_order = data.dataOnEntities[MBTRI][ff].getOrder();
    data.dataOnEntities[MBTRI][ff].getN(base).resize(
        nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
        nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueHdivDemkowiczBrokenBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (base != DEMKOWICZ_JACOBI_BASE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "This should be used only with DEMKOWICZ_JACOBI_BASE "
             "but base is %s",
             ApproximationBaseNames[base]);
  }
  int nb_gauss_pts = pts.size2();

  int volume_order = data.dataOnEntities[MBTET][0].getOrder();
  int nb_dofs_face = NBFACETRI_DEMKOWICZ_HDIV(volume_order);
  int nb_dofs_volume = NBVOLUMETET_DEMKOWICZ_HDIV(volume_order);
  int nb_dofs = 4 * nb_dofs_face + nb_dofs_volume;
  data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                  false);
  data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts, 9 * nb_dofs,
                                                      false);
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);

  auto get_interior_cache = [this]() -> TetBaseCache::BaseCacheMI * {
    if (vPtr) {
      auto it = TetBaseCache::hdivBrokenBaseInteriorDemkowicz.find(vPtr);
      if (it != TetBaseCache::hdivBrokenBaseInteriorDemkowicz.end()) {
        return &it->second;
      }
    }
    return nullptr;
  };

  auto interior_cache_ptr = get_interior_cache();

  if (interior_cache_ptr) {
    auto it =
        interior_cache_ptr->find(boost::make_tuple(volume_order, nb_gauss_pts));
    if (it != interior_cache_ptr->end()) {
      noalias(data.dataOnEntities[MBTET][0].getN(base)) = it->N;
      noalias(data.dataOnEntities[MBTET][0].getDiffN(base)) = it->diffN;
      MoFEMFunctionBeginHot(0);
    }
  }

  std::array<MatrixDouble, 4> face_base_fun{
      MatrixDouble(nb_gauss_pts, 3 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 3 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 3 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 3 * nb_dofs_face)};
  std::array<MatrixDouble, 4> face_diff_base{
      MatrixDouble(nb_gauss_pts, 9 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 9 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 9 * nb_dofs_face),
      MatrixDouble(nb_gauss_pts, 9 * nb_dofs_face)};

  int face_node[4][3] = {{0, 1, 3}, {1, 2, 3}, {0, 2, 3}, {0, 1, 2}};

  std::array<int, 4> p_f{volume_order, volume_order, volume_order,
                         volume_order};
  std::array<double *, 4> phi_f{
      &*face_base_fun[0].data().begin(), &*face_base_fun[1].data().begin(),
      &*face_base_fun[2].data().begin(), &*face_base_fun[3].data().begin()};
  std::array<double *, 4> diff_phi_f{
      &*face_diff_base[0].data().begin(), &*face_diff_base[1].data().begin(),
      &*face_diff_base[2].data().begin(), &*face_diff_base[3].data().begin()};

  // Calculate base function on tet faces
  for (int ff = 0; ff != 4; ff++) {
    CHKERR Hdiv_Demkowicz_Face_MBTET_ON_FACE(
        face_node[ff], p_f[ff],
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        phi_f[ff], diff_phi_f[ff], nb_gauss_pts, 4);
  }

  MatrixDouble vol_bases(nb_gauss_pts, 3 * nb_dofs_volume);
  MatrixDouble vol_diff_bases(nb_gauss_pts, 9 * nb_dofs_volume);
  auto *phi_v = &*vol_bases.data().begin();
  auto *diff_phi_v = &*vol_diff_bases.data().begin();
  CHKERR Hdiv_Demkowicz_Interior_MBTET(
      volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
      &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), p_f.data(),
      phi_f.data(), diff_phi_f.data(), phi_v, diff_phi_v, nb_gauss_pts);

  // faces
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v_f[] = {
      getFTensor1FromPtr<3>(phi_f[0]), getFTensor1FromPtr<3>(phi_f[1]),
      getFTensor1FromPtr<3>(phi_f[2]), getFTensor1FromPtr<3>(phi_f[3])};
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v_f[] = {
      getFTensor2HVecFromPtr<3, 3>(diff_phi_f[0]),
      getFTensor2HVecFromPtr<3, 3>(diff_phi_f[1]),
      getFTensor2HVecFromPtr<3, 3>(diff_phi_f[2]),
      getFTensor2HVecFromPtr<3, 3>(diff_phi_f[3])};

  // volumes
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_base_v =
      getFTensor1FromPtr<3>(&*vol_bases.data().begin());
  FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3> t_diff_base_v =
      getFTensor2HVecFromPtr<3, 3>(&*vol_diff_bases.data().begin());

  auto t_base = getFTensor1FromPtr<3>(
      &*data.dataOnEntities[MBTET][0].getN(base).data().begin());
  auto t_diff_base = getFTensor2HVecFromPtr<3, 3>(
      &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin());

  i_FTIndex<3> i;
  j_FTIndex<3> j;

  for (auto gg = 0; gg != nb_gauss_pts; ++gg) {
    for (int oo = 0; oo < volume_order; oo++) {
      // face
      for (auto dd = NBFACETRI_DEMKOWICZ_HDIV(oo);
           dd != NBFACETRI_DEMKOWICZ_HDIV(oo + 1); ++dd) {
        for (auto ff = 0; ff != 4; ++ff) {
          t_base(i) = t_base_v_f[ff](i);
          ++t_base;
          ++t_base_v_f[ff];
          t_diff_base(i, j) = t_diff_base_v_f[ff](i, j);
          ++t_diff_base;
          ++t_diff_base_v_f[ff];
        }
      }
      // volume
      for (auto dd = NBVOLUMETET_DEMKOWICZ_HDIV(oo);
           dd != NBVOLUMETET_DEMKOWICZ_HDIV(oo + 1); ++dd) {
        t_base(i) = t_base_v(i);
        ++t_base;
        ++t_base_v;
        t_diff_base(i, j) = t_diff_base_v(i, j);
        ++t_diff_base;
        ++t_diff_base_v;
      }
    }
  }

  if (interior_cache_ptr) {
    auto p = interior_cache_ptr->emplace(
        TetBaseCache::BaseCacheItem{volume_order, nb_gauss_pts});
    p.first->N = data.dataOnEntities[MBTET][0].getN(base);
    p.first->diffN = data.dataOnEntities[MBTET][0].getDiffN(base);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TetPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  switch (cTx->spaceContinuity) {
  case CONTINUOUS:
    switch (cTx->bAse) {
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
      return getValueHdivAinsworthBase(pts);
    case DEMKOWICZ_JACOBI_BASE:
      return getValueHdivDemkowiczBase(pts);
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
    }
    break;
  case DISCONTINUOUS:
    switch (cTx->bAse) {
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBATTO_BASE:
      return getValueHdivAinsworthBrokenBase(pts);
    case DEMKOWICZ_JACOBI_BASE:
      return getValueHdivDemkowiczBrokenBase(pts);
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
    }
    break;

  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown continuity");
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
TetPolynomialBase::getValueHcurlAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
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
      order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
      int nb_dofs =
          NBEDGE_AINSWORTH_HCURL(data.dataOnEntities[MBEDGE][ee].getOrder());
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
      order[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
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
    int order = data.dataOnEntities[MBTET][0].getOrder();
    int nb_vol_dofs = NBVOLUMETET_AINSWORTH_HCURL(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,
                                                    3 * nb_vol_dofs, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        9 * nb_vol_dofs, false);
    CHKERR Hcurl_Ainsworth_VolumeFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getOrder(),
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

  EntitiesFieldData &data = cTx->dAta;
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
      order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
      int nb_dofs =
          NBEDGE_DEMKOWICZ_HCURL(data.dataOnEntities[MBEDGE][ee].getOrder());
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
      order[ff] = data.dataOnEntities[MBTRI][ff].getOrder();
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
    int order = data.dataOnEntities[MBTET][0].getOrder();
    int nb_vol_dofs = NBVOLUMETET_DEMKOWICZ_HCURL(order);
    data.dataOnEntities[MBTET][0].getN(base).resize(nb_gauss_pts,
                                                    3 * nb_vol_dofs, false);
    data.dataOnEntities[MBTET][0].getDiffN(base).resize(nb_gauss_pts,
                                                        9 * nb_vol_dofs, false);
    CHKERR Hcurl_Demkowicz_VolumeBaseFunctions_MBTET(
        data.dataOnEntities[MBTET][0].getOrder(),
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

  cTx = ctx_ptr->getInterface<EntPolynomialBaseCtx>();

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 3)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  if (cTx->bAse != AINSWORTH_BERNSTEIN_BEZIER_BASE) {
    const FieldApproximationBase base = cTx->bAse;
    EntitiesFieldData &data = cTx->dAta;
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

  switch (cTx->spaceContinuity) {
  case CONTINUOUS:

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
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown space %s",
               FieldSpaceNames[cTx->sPace]);
    }
    break;

  case DISCONTINUOUS:

    switch (cTx->sPace) {
    case HDIV:
      CHKERR getValueHdiv(pts);
      break;
    default:
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown space %s",
               FieldSpaceNames[cTx->sPace]);
    }
    break;

  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Unknown continuity");
  }

  MoFEMFunctionReturn(0);
}

template <typename T>
auto tetCacheSwitch(const void *ptr, T &cache, std::string cache_name) {
  auto it = cache.find(ptr);
  if (it != cache.end()) {
    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::noisy, "TetPolynomialBase")
        << "Cache off " << cache_name << ": " << it->second.size();
    cache.erase(it);
    return false;
  } else {
    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_TAG_AND_LOG("WORLD", Sev::noisy, "TetPolynomialBase")
        << "Cache on " << cache_name;
    cache[ptr];
    return true;
  }
}

bool TetPolynomialBase::switchCacheHDivBaseFaceDemkowicz(const void *ptr) {
  return tetCacheSwitch(ptr, TetBaseCache::hDivBaseFaceDemkowicz,
                        "hDivBaseFaceDemkowicz");
}

bool TetPolynomialBase::switchCacheHdivBaseInteriorDemkowicz(const void *ptr) {
  return tetCacheSwitch(ptr, TetBaseCache::hdivBaseInteriorDemkowicz,
                        "hdivBaseInteriorDemkowicz");
}

bool TetPolynomialBase::switchCacheHdivBrokenBaseInteriorDemkowicz(
    const void *ptr) {
  return tetCacheSwitch(ptr, TetBaseCache::hdivBrokenBaseInteriorDemkowicz,
                        "hdivBrokenBaseInteriorDemkowicz");
}

void TetPolynomialBase::switchCacheHDivBaseDemkowiczOn(std::vector<void *> v) {
  for (auto fe_ptr : v) {
    if (!TetPolynomialBase::switchCacheHDivBaseFaceDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHDivBaseFaceDemkowicz(fe_ptr);
    }
    if (!TetPolynomialBase::switchCacheHdivBaseInteriorDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBaseInteriorDemkowicz(fe_ptr);
    }
    if (!TetPolynomialBase::switchCacheHdivBrokenBaseInteriorDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBrokenBaseInteriorDemkowicz(fe_ptr);
    }
  }
}

void TetPolynomialBase::switchCacheHDivBaseDemkowiczOff(std::vector<void *> v) {
  for (auto fe_ptr : v) {
    if (TetPolynomialBase::switchCacheHDivBaseFaceDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHDivBaseFaceDemkowicz(fe_ptr);
    }
    if (TetPolynomialBase::switchCacheHdivBaseInteriorDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBaseInteriorDemkowicz(fe_ptr);
    }
    if (TetPolynomialBase::switchCacheHdivBrokenBaseInteriorDemkowicz(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBrokenBaseInteriorDemkowicz(fe_ptr);
    }
  }
}

bool TetPolynomialBase::switchCacheHdivBrokenBaseInteriorAinsworth(
    const void *ptr) {
  return tetCacheSwitch(ptr, TetBaseCache::hdivBrokenBaseInteriorAinsworth,
                        "hdivBrokenBaseInteriorAinsworth");
}

void TetPolynomialBase::switchCacheHDivBaseAinsworthOn(std::vector<void *> v) {
  for (auto fe_ptr : v) {
    if (!TetPolynomialBase::switchCacheHdivBrokenBaseInteriorAinsworth(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBrokenBaseInteriorAinsworth(fe_ptr);
    }
  }
}

void TetPolynomialBase::switchCacheHDivBaseAinsworthOff(std::vector<void *> v) {
  for (auto fe_ptr : v) {
    if (TetPolynomialBase::switchCacheHdivBrokenBaseInteriorAinsworth(fe_ptr)) {
      TetPolynomialBase::switchCacheHdivBrokenBaseInteriorAinsworth(fe_ptr);
    }
  }
}

void TetPolynomialBase::switchCacheHDivBaseOn(std::vector<void *> v) {
  switchCacheHDivBaseDemkowiczOn(v);
  switchCacheHDivBaseAinsworthOn(v);
}

void TetPolynomialBase::switchCacheHDivBaseOff(std::vector<void *> v) {
  switchCacheHDivBaseDemkowiczOff(v);
  switchCacheHDivBaseAinsworthOff(v);
}
