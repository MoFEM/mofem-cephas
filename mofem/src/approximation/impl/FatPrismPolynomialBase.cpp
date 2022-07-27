/** \file FatPrismPolynomialBase.cpp
\brief Implementation of Ainsworth-Cole H1 base on edge

\todo Prism element can be integrated exploiting tonsorial product. Current
implementation do not take that opportunity. That can be viewed as a bug. 

*/



using namespace MoFEM;

MoFEMErrorCode FatPrismPolynomialBaseCtx::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<FatPrismPolynomialBaseCtx *>(this);
  return 0;
}

FatPrismPolynomialBaseCtx::FatPrismPolynomialBaseCtx(
    EntitiesFieldData &data,
    EntitiesFieldData &data_triangles_only,
    EntitiesFieldData &data_trough_thickness,
    MatrixDouble &gauss_pts_triangles_only,
    MatrixDouble &gauss_pts_through_thickness, moab::Interface &moab,
    const NumeredEntFiniteElement *fe_ptr, const FieldSpace space,
    const FieldApproximationBase base,
    const FieldApproximationBase copy_node_base)
    : EntPolynomialBaseCtx(data, space, base, copy_node_base),
      dataTrianglesOnly(data_triangles_only),
      dataTroughThickness(data_trough_thickness),
      gaussPtsTrianglesOnly(gauss_pts_triangles_only),
      gaussPtsThroughThickness(gauss_pts_through_thickness), mOab(moab),
      fePtr(fe_ptr) {

  ierr = setBase();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}
FatPrismPolynomialBaseCtx::~FatPrismPolynomialBaseCtx() {}

MoFEMErrorCode
FatPrismPolynomialBase::query_interface(boost::typeindex::type_index type_index,
                                        UnknownInterface **iface) const {
  *iface = const_cast<FatPrismPolynomialBase *>(this);
  return 0;
}

FatPrismPolynomialBase::~FatPrismPolynomialBase() {}
FatPrismPolynomialBase::FatPrismPolynomialBase() {}

MoFEMErrorCode
FatPrismPolynomialBase::getValue(MatrixDouble &pts,
                                 boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {

  MoFEMFunctionBegin;

  cTx = ctx_ptr->getInterface<FatPrismPolynomialBaseCtx>();

  if (!cTx->fePtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer to element should be given "
            "when EntPolynomialBaseCtx is constructed "
            "(use different constructor)");

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 3)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  const FieldApproximationBase base = cTx->bAse;
  EntitiesFieldData &data = cTx->dAta;

  if (cTx->copyNodeBase == LASTBASE)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "It is assumed that base for vertices is calculated");
  else
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
        data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 6, false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts, 12,
                                                         false);
  if (data.dataOnEntities[MBVERTEX][0].getN(base).size1() !=
      (unsigned int)nb_gauss_pts)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[base]);

  if (cTx->gaussPtsTrianglesOnly.size2() *
          cTx->gaussPtsThroughThickness.size2() !=
      pts.size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  switch (cTx->sPace) {
  case H1:
    CHKERR getValueH1TrianglesOnly();
    CHKERR getValueH1ThroughThickness();
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

MoFEMErrorCode FatPrismPolynomialBase::getValueH1TrianglesOnly() {
  MoFEMFunctionBegin;
  const FieldApproximationBase base = cTx->bAse;
  CHKERR FlatPrismPolynomialBase().getValue(
      cTx->gaussPtsTrianglesOnly,
      boost::shared_ptr<BaseFunctionCtx>(new FlatPrismPolynomialBaseCtx(
          cTx->dataTrianglesOnly, cTx->mOab, cTx->fePtr, H1, base, NOBASE)));
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FatPrismPolynomialBase::getValueH1ThroughThickness() {
  MoFEMFunctionBegin;

  // EntitiesFieldData& data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts_through_thickness = cTx->gaussPtsThroughThickness.size2();
  // Calculate Legendre approx. on edges
  for (unsigned int ee = 3; ee <= 5; ee++) {
    auto &ent_data = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee];
    int sense = ent_data.getSense();
    int order = ent_data.getOrder();
    int nb_dofs = NBEDGE_H1(order);
    ent_data.getN(base).resize(nb_gauss_pts_through_thickness, nb_dofs, false);
    ent_data.getDiffN(base).resize(nb_gauss_pts_through_thickness, nb_dofs,
                                   false);

    if (nb_dofs > 0) {
      double diff_s = 1.; // s = s(xi), ds/dxi = 2., because change of basis
      for (int gg = 0; gg < nb_gauss_pts_through_thickness; gg++) {
        double s =
            2 * cTx->gaussPtsThroughThickness(0, gg) - 1; // makes form -1..1
        if (sense == -1) {
          s *= -2;
          diff_s *= -2;
        }
        // calculate Legendre polynomials at integration points on edges
        // thorough thickness
        CHKERR base_polynomials(order - 2, s, &diff_s,
                                &ent_data.getN(base)(gg, 0),
                                &ent_data.getDiffN(base)(gg, 0), 1);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FatPrismPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  EntitiesFieldData &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();
  int nb_gauss_pts_on_faces = cTx->gaussPtsTrianglesOnly.size2();
  int nb_gauss_pts_through_thickness = cTx->gaussPtsThroughThickness.size2();

  // nodes
  // linear for xi, eta and zeta
  auto &vert_dat = data.dataOnEntities[MBVERTEX][0];
  vert_dat.getN(base).resize(nb_gauss_pts, 6, false);
  vert_dat.getDiffN(base).resize(nb_gauss_pts, 18);
  noalias(vert_dat.getN(base)) = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);
  noalias(vert_dat.getDiffN(base)) =
      data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE);

  auto &vert_n = vert_dat.getN(base);
  auto &vert_diff_n  = vert_dat.getDiffN(base);

  constexpr int prism_edge_map[9][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 4},
                                        {2, 5}, {3, 4}, {4, 5}, {5, 3}};

  auto edge_through_thickness = [&](const int ee) {
    MoFEMFunctionBegin;
    // through thickness ho approximation
    // linear xi,eta, ho terms for zeta

    auto &thickness_ent = cTx->dataTroughThickness.dataOnEntities[MBEDGE][ee];
    auto &prism_ent = data.dataOnEntities[MBEDGE][ee];

    int order = thickness_ent.getOrder();
    int nb_dofs = NBEDGE_H1(order);
    if ((unsigned int)nb_dofs != thickness_ent.getN(base).size2())
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "nb_dofs != nb_dofs %d != %d", nb_dofs,
               thickness_ent.getN(base).size2());
    prism_ent.getN(base).resize(nb_gauss_pts, nb_dofs, false);
    prism_ent.getDiffN(base).resize(nb_gauss_pts, 3 * nb_dofs, false);
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);


    int gg = 0;
    for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {

      for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
        double extrude = vert_n(gg, prism_edge_map[ee][0]) +
                         vert_n(gg, prism_edge_map[ee][1]);

        double diff_extrude[3];
        for (auto d : {0, 1, 2})
          diff_extrude[d] = vert_diff_n(gg, 3 * prism_edge_map[ee][0] + d) +
                            vert_diff_n(gg, 3 * prism_edge_map[ee][1] + d);

        double n0 = vert_n(gg, 0) + vert_n(gg, 1) + vert_n(gg, 2);
        double n1 = vert_n(gg, 3) + vert_n(gg, 4) + vert_n(gg, 5);

        // Calculate base through thickness directly from shape fuctions. I
        // leave that bit of the code, could be useful for integration on the
        // skeleton.
        
        // double l[order + 1], diff_l[3 * (order + 1)];
        // double ksi = n1 - n0;
        // double diff_ksi[3];
        // for (auto d : {0, 1, 2})
        //   diff_ksi[d] =
        //       vert_diff_n(gg, 3 * 3 + d) + vert_diff_n(gg, 3 * 4 + d) +
        //       vert_diff_n(gg, 3 * 5 + d) - vert_diff_n(gg, 3 * 0 + d) -
        //       vert_diff_n(gg, 3 * 1 + d) - vert_diff_n(gg, 3 * 2 + d);
        // if(sense == -1) {
        //   ksi *= -1;
        //   for (auto d : {0, 1, 2})
        //     diff_ksi[d] *= -1;
        // }
        // CHKERR base_polynomials(order - 2, ksi, diff_ksi, l, diff_l, 3);

        double *l = &(thickness_ent.getN(base)(ggt, 0));
        double *diff_l_1d = &(thickness_ent.getDiffN(base)(ggt, 0));

        double bubble = n0 * n1;
        double diff_bubble[3];
        for (auto d : {0, 1, 2})
          diff_bubble[d] =
              n0 * (vert_diff_n(gg, 3 * 3 + d) + vert_diff_n(gg, 3 * 4 + d) +
                    vert_diff_n(gg, 3 * 5 + d)) +

              n1 * (vert_diff_n(gg, 3 * 0 + d) + vert_diff_n(gg, 3 * 1 + d) +
                    vert_diff_n(gg, 3 * 2 + d));

        for (int dd = 0; dd != nb_dofs; dd++) {

          prism_ent.getN(base)(gg, dd) = extrude * bubble * l[dd];
          for (auto d : {0, 1, 2})
            prism_ent.getDiffN(base)(gg, 3 * dd + d) =
                diff_extrude[d] * bubble * l[dd] +
                extrude * diff_bubble[d] * l[dd];

          prism_ent.getDiffN(base)(gg, 3 * dd + 2) +=
              extrude * bubble * 2 * diff_l_1d[dd];

        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto edge_on_the_triangle = [&](const int ee) {
    MoFEMFunctionBegin;
    // on triangles ho approximation
    // ho terms on edges, linear zeta
    int nb_dofs =
        cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getN(base).size2();
    data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                      false);
    data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                          3 * nb_dofs, false);
    for (int dd = 0; dd < nb_dofs; dd++) {
      int gg = 0;
      for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {
        double tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getN(
            base)(ggf, dd);
        double dksi_tri_n =
            cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN(base)(
                ggf, 2 * dd + 0);
        double deta_tri_n =
            cTx->dataTrianglesOnly.dataOnEntities[MBEDGE][ee].getDiffN(base)(
                ggf, 2 * dd + 1);
        for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
          double zeta = cTx->gaussPtsThroughThickness(0, ggt);
          double dzeta, edge_shape;
          if (ee < 3) {
            dzeta = diffN_MBEDGE0;
            edge_shape = N_MBEDGE0(zeta);
          } else {
            dzeta = diffN_MBEDGE1;
            edge_shape = N_MBEDGE1(zeta);
          }
          data.dataOnEntities[MBEDGE][ee].getN(base)(gg, dd) =
              tri_n * edge_shape;
          data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg, 3 * dd + 0) =
              dksi_tri_n * edge_shape;
          data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg, 3 * dd + 1) =
              deta_tri_n * edge_shape;
          data.dataOnEntities[MBEDGE][ee].getDiffN(base)(gg, 3 * dd + 2) =
              tri_n * dzeta;
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto trinagle_through_thickness = [&](const int ff) {
    MoFEMFunctionBegin;
    int nb_dofs;
    nb_dofs =
        cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getN(base).size2();
    data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                     false);
    data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,
                                                         3 * nb_dofs, false);
    for (int dd = 0; dd < nb_dofs; dd++) {
      int gg = 0;
      for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {
        double tri_n = cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getN(
            base)(ggf, dd);
        double dksi_tri_n[2];
        for (int kk = 0; kk < 2; kk++) {
          dksi_tri_n[kk] =
              cTx->dataTrianglesOnly.dataOnEntities[MBTRI][ff].getDiffN(base)(
                  ggf, 2 * dd + kk);
        }
        for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
          double zeta = cTx->gaussPtsThroughThickness(0, ggt);
          double dzeta, edge_shape;
          if (ff == 3) {
            dzeta = diffN_MBEDGE0;
            edge_shape = N_MBEDGE0(zeta);
          } else {
            dzeta = diffN_MBEDGE1;
            edge_shape = N_MBEDGE1(zeta);
          }
          data.dataOnEntities[MBTRI][ff].getN(base)(gg, dd) =
              tri_n * edge_shape;
          for (auto kk : {0, 1}) 
            data.dataOnEntities[MBTRI][ff].getDiffN(base)(gg, 3 * dd + kk) =
                dksi_tri_n[kk] * edge_shape;
          data.dataOnEntities[MBTRI][ff].getDiffN(base)(gg, 3 * dd + 2) =
              tri_n * dzeta;
        }
      }
    }
    MoFEMFunctionReturn(0);
  };

  auto quads_base = [&]() {
    MoFEMFunctionBegin;
    int quads_nodes[3 * 4];
    int quad_order[3] = {0, 0, 0};
    double *quad_n[3], *diff_quad_n[3];
    SideNumber_multiIndex &side_table =
        const_cast<SideNumber_multiIndex &>(cTx->fePtr->getSideNumberTable());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit;
    siit = side_table.get<1>().lower_bound(boost::make_tuple(MBQUAD, 0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit;
    hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBQUAD, 3));
    if (std::distance(siit, hi_siit) != 3)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected three quads on prisms");
    EntityHandle ent = cTx->fePtr->getEnt();
    int num_nodes;
    const EntityHandle *conn;
    CHKERR cTx->mOab.get_connectivity(ent, conn, num_nodes, true);
    for (; siit != hi_siit; ++siit) {
      int side = siit->get()->side_number;
      if(side < 0 && side > 2)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Side in range 0,1,2 expected");
      int num_nodes_quad;
      const EntityHandle *conn_quad;
      EntityHandle quad = siit->get()->ent;
      CHKERR cTx->mOab.get_connectivity(quad, conn_quad, num_nodes_quad, true);
      for (int nn = 0; nn < num_nodes_quad; nn++) {
        quads_nodes[4 * side + nn] =
            std::distance(conn, std::find(conn, conn + 6, conn_quad[nn]));
      }
      int order = data.dataOnEntities[MBQUAD][side].getOrder();
      quad_order[side] = order;
      data.dataOnEntities[MBQUAD][side].getN(base).resize(
          nb_gauss_pts, NBFACEQUAD_H1(order), false);
      data.dataOnEntities[MBQUAD][side].getDiffN(base).resize(
          nb_gauss_pts, 3 * NBFACEQUAD_H1(order), false);
      if (data.dataOnEntities[MBQUAD][side].getN(base).size2() > 0) {
        quad_n[side] =
            &*data.dataOnEntities[MBQUAD][side].getN(base).data().begin();
        diff_quad_n[side] =
            &*data.dataOnEntities[MBQUAD][side].getDiffN(base).data().begin();
      } else {
        quad_n[side] = NULL;
        diff_quad_n[side] = NULL;
      }
    }
    if (quad_order[0] > 0 || quad_order[1] > 0 || quad_order[2] > 0) {
      double *vertex_n =
          &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin();
      double *diff_vertex_n =
          &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin();
      CHKERR H1_QuadShapeFunctions_MBPRISM(quads_nodes, quad_order, vertex_n,
                                           diff_vertex_n, quad_n, diff_quad_n,
                                           nb_gauss_pts, base_polynomials);
    }
    MoFEMFunctionReturn(0);
  };

  auto prim_base = [&]() {
    MoFEMFunctionBegin;
    int order = data.dataOnEntities[MBPRISM][0].getOrder();
    double *vertex_n = &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0);
    double *diff_vertex_n =
        &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0);
    data.dataOnEntities[MBPRISM][0].getN(base).resize(
        nb_gauss_pts, NBVOLUMEPRISM_H1(order), false);
    data.dataOnEntities[MBPRISM][0].getDiffN(base).resize(
        nb_gauss_pts, 3 * NBVOLUMEPRISM_H1(order), false);
    if (NBVOLUMEPRISM_H1(order) > 0) {
      CHKERR H1_VolumeShapeFunctions_MBPRISM(
          order, vertex_n, diff_vertex_n,
          &data.dataOnEntities[MBPRISM][0].getN(base)(0, 0),
          &data.dataOnEntities[MBPRISM][0].getDiffN(base)(0, 0), nb_gauss_pts,
          base_polynomials);
    }
    MoFEMFunctionReturn(0);
  };

  // edges on triangles
  int ee = 0;
  for (; ee <= 2; ++ee)
    CHKERR edge_on_the_triangle(ee);
  for (; ee <= 5; ++ee)
    CHKERR edge_through_thickness(ee);
  for (; ee <= 8; ++ee)
    CHKERR edge_on_the_triangle(ee);

  // triangles
  // ho on triangles, linear zeta
  for (int ff = 3; ff <= 4; ++ff)
    CHKERR trinagle_through_thickness(ff);

  // quads
  // higher order edges and zeta
  CHKERR quads_base();

  // prism
  CHKERR prim_base();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FatPrismPolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FatPrismPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
}

MoFEMErrorCode FatPrismPolynomialBase::getValueHcurl(MatrixDouble &pts) {
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
}
