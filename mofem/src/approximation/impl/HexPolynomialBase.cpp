/** \file HexPolynomialBase.cpp
\brief Implementation of hierarchical bases on tetrahedral

A l2, h1, h-div and h-curl spaces are implemented.

*/



using namespace MoFEM;

MoFEMErrorCode
HexPolynomialBase::query_interface(boost::typeindex::type_index type_index,
                                   UnknownInterface **iface) const {
  MoFEMFunctionBegin;
  *iface = const_cast<HexPolynomialBase *>(this);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueH1DemkowiczBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueH1DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  auto &data = cTx->dAta;
  const auto base = cTx->bAse;
  const auto copy_base = cTx->copyNodeBase;
  int nb_gauss_pts = pts.size2();

  auto &copy_base_fun = data.dataOnEntities[MBVERTEX][0].getN(copy_base);
  auto &copy_diff_base_fun =
      data.dataOnEntities[MBVERTEX][0].getDiffN(copy_base);
  if (copy_base_fun.size1() != nb_gauss_pts)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Inconsistent number of integration pts");

  auto add_base_on_verts = [&] {
    MoFEMFunctionBeginHot;
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(
        nb_gauss_pts, copy_base_fun.size2(), false);
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(
        nb_gauss_pts, copy_diff_base_fun.size2(), false);
    noalias(data.dataOnEntities[MBVERTEX][0].getN(base)) = copy_base_fun;
    noalias(data.dataOnEntities[MBVERTEX][0].getDiffN(base)) =
        copy_diff_base_fun;
    MoFEMFunctionReturnHot(0);
  };

  // Edges
  auto add_base_on_edges = [&] {
    MoFEMFunctionBeginHot;
    std::array<int, 12> sense;
    std::array<int, 12> order;
    if (data.spacesOnEntities[MBEDGE].test(H1)) {
      if (data.dataOnEntities[MBEDGE].size() != 12)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected 12 data on entities");

      std::array<double *, 12> h1_edge_n;
      std::array<double *, 12> diff_h1_egde_n;
      bool nb_dofs_sum = false;
      for (int ee = 0; ee != 12; ++ee) {
        if (data.dataOnEntities[MBEDGE][ee].getSense() == 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Sense of entity not set");

        sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
        order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();

        int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getOrder());
        data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                          false);
        data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(
            nb_gauss_pts, 3 * nb_dofs, false);
        h1_edge_n[ee] =
            &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
        diff_h1_egde_n[ee] =
            &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();

        nb_dofs_sum |= (nb_dofs > 0);
      }
      if (nb_dofs_sum) {
        CHKERR DemkowiczHexAndQuad::H1_EdgeShapeFunctions_ONHEX(
            sense.data(), order.data(), &*copy_base_fun.data().begin(),
            &*copy_diff_base_fun.data().begin(), h1_edge_n.data(),
            diff_h1_egde_n.data(), nb_gauss_pts);
      }
    } else {
      for (int ee = 0; ee != 12; ++ee) {
        data.dataOnEntities[MBEDGE][ee].getN(base).resize(0, 0, false);
        data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(0, 0, false);
      }
    }
    MoFEMFunctionReturnHot(0);
  };

  // Face
  auto add_base_on_quads = [&]() {
    MoFEMFunctionBeginHot;
    std::array<int, 6> order;
    if (data.spacesOnEntities[MBQUAD].test(H1)) {
      // faces
      if (data.dataOnEntities[MBQUAD].size() != 6)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected six faces on hex");

      std::array<double *, 6> h1_face_n;
      std::array<double *, 6> diff_h1_face_n;
      bool nb_dofs_sum = false;
      for (int ff = 0; ff != 6; ++ff) {

        order[ff] = data.dataOnEntities[MBQUAD][ff].getOrder();
        const int nb_dofs = NBFACEQUAD_H1(order[ff]);

        data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                          false);
        data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(
            nb_gauss_pts, 3 * nb_dofs, false);

        h1_face_n[ff] =
            &*data.dataOnEntities[MBQUAD][ff].getN(base).data().begin();
        diff_h1_face_n[ff] =
            &*data.dataOnEntities[MBQUAD][ff].getDiffN(base).data().begin();

        nb_dofs_sum |= (nb_dofs > 0);
      }
      if (data.facesNodes.size1() != 6)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected six face nodes");
      if (data.facesNodes.size2() != 4)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected four nodes on face");

      if (nb_dofs_sum) {
        CHKERR DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONHEX(
            &*data.facesNodes.data().begin(),
            &*data.facesNodesOrder.data().begin(), order.data(),
            &*copy_base_fun.data().begin(), &*copy_diff_base_fun.data().begin(),
            h1_face_n.data(), diff_h1_face_n.data(), nb_gauss_pts);
      }

    } else {
      for (int ff = 0; ff != 6; ++ff) {
        data.dataOnEntities[MBQUAD][ff].getN(base).resize(0, false);
        data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(0, 0, false);
      }
    }

    MoFEMFunctionReturnHot(0);
  };

  // Face
  auto add_base_on_volume = [&]() {
    MoFEMFunctionBeginHot;

    if (data.spacesOnEntities[MBHEX].test(H1)) {
      // volume
      int order = data.dataOnEntities[MBHEX][0].getOrder();
      int nb_vol_dofs = NBVOLUMEHEX_H1(order);
      data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, nb_vol_dofs,
                                                      false);
      data.dataOnEntities[MBHEX][0].getDiffN(base).resize(
          nb_gauss_pts, 3 * nb_vol_dofs, false);

      if (nb_vol_dofs) {
        const std::array<int, 3> p{order, order, order};
        CHKERR DemkowiczHexAndQuad::H1_InteriorShapeFunctions_ONHEX(
            p.data(), &*copy_base_fun.data().begin(),
            &*copy_diff_base_fun.data().begin(),
            &*data.dataOnEntities[MBHEX][0].getN(base).data().begin(),
            &*data.dataOnEntities[MBHEX][0].getDiffN(base).data().begin(),
            nb_gauss_pts);
      }
    } else {
      data.dataOnEntities[MBHEX][0].getN(base).resize(0, 0, false);
      data.dataOnEntities[MBHEX][0].getDiffN(base).resize(0, 0, false);
    }

    MoFEMFunctionReturnHot(0);
  };

  CHKERR add_base_on_verts();
  CHKERR add_base_on_edges();
  CHKERR add_base_on_quads();
  CHKERR add_base_on_volume();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueL2DemkowiczBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueL2DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  auto &data = cTx->dAta;
  const auto base = cTx->bAse;
  const auto copy_base = cTx->copyNodeBase;

  int nb_gauss_pts = pts.size2();

  auto &copy_base_fun = data.dataOnEntities[MBVERTEX][0].getN(copy_base);
  auto &copy_diff_base_fun =
      data.dataOnEntities[MBVERTEX][0].getDiffN(copy_base);

  if (nb_gauss_pts != copy_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  if (nb_gauss_pts != copy_diff_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  auto &base_fun = data.dataOnEntities[MBHEX][0].getN(base);
  auto &diff_base_fun = data.dataOnEntities[MBHEX][0].getDiffN(base);
  const int vol_order = data.dataOnEntities[MBHEX][0].getOrder();

  const int nb_dofs = NBVOLUMEHEX_L2(vol_order);
  base_fun.resize(nb_gauss_pts, nb_dofs, false);
  diff_base_fun.resize(nb_gauss_pts, 3 * nb_dofs, false);

  if (nb_dofs) {
    const std::array<int, 3> p{vol_order, vol_order, vol_order};
    CHKERR DemkowiczHexAndQuad::L2_InteriorShapeFunctions_ONHEX(
        p.data(), &*copy_base_fun.data().begin(),
        &*copy_diff_base_fun.data().begin(), &*base_fun.data().begin(),
        &*diff_base_fun.data().begin(), nb_gauss_pts);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueHdivDemkowiczBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueHdivDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  auto &data = cTx->dAta;
  const auto base = cTx->bAse;
  const auto copy_base = cTx->copyNodeBase;
  const int nb_gauss_pts = pts.size2();

  auto &copy_base_fun = data.dataOnEntities[MBVERTEX][0].getN(copy_base);
  auto &copy_diff_base_fun =
      data.dataOnEntities[MBVERTEX][0].getDiffN(copy_base);

  if (nb_gauss_pts != copy_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  if (nb_gauss_pts != copy_diff_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  // Quad
  if (data.spacesOnEntities[MBQUAD].test(HDIV)) {

    if (data.dataOnEntities[MBQUAD].size() != 6)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected six data structures on Hex");

    std::array<int, 6> order;
    std::array<double *, 6> hdiv_face_n;
    std::array<double *, 6> diff_hdiv_face_n;

    bool sum_nb_dofs = false;
    for (int ff = 0; ff != 6; ff++) {
      if (data.dataOnEntities[MBQUAD][ff].getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Sense pn quad <%d> not set", ff);

      order[ff] = data.dataOnEntities[MBQUAD][ff].getOrder();
      if (data.facesNodes.size1() != 6)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected six faces");
      if (data.facesNodes.size2() != 4)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected four  nodes on face");

      const int nb_dofs = NBFACEQUAD_DEMKOWICZ_HDIV(order[ff]);
      auto &face_n = data.dataOnEntities[MBQUAD][ff].getN(base);
      auto &diff_face_n = data.dataOnEntities[MBQUAD][ff].getDiffN(base);
      face_n.resize(nb_gauss_pts, 3 * nb_dofs, false);
      diff_face_n.resize(nb_gauss_pts, 9 * nb_dofs, false);

      hdiv_face_n[ff] = &*face_n.data().begin();
      diff_hdiv_face_n[ff] = &*diff_face_n.data().begin();

      sum_nb_dofs |= (nb_dofs > 0);
    }

    if (sum_nb_dofs) {
      DemkowiczHexAndQuad::Hdiv_FaceShapeFunctions_ONHEX(
          &*data.facesNodes.data().begin(),
          &*data.facesNodesOrder.data().begin(), order.data(),
          &*copy_base_fun.data().begin(), &*copy_diff_base_fun.data().begin(),
          hdiv_face_n.data(), diff_hdiv_face_n.data(), nb_gauss_pts);

    } else {
      for (int ff = 0; ff != 6; ff++) {
        data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, 0,
                                                          false);
        data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                              false);
      }
    }

  } else {
    for (int ff = 0; ff != 6; ff++) {
      data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  // Hex
  if (data.spacesOnEntities[MBHEX].test(HDIV)) {

    const int order = data.dataOnEntities[MBHEX][0].getOrder();
    const int nb_dofs_family =
        NBVOLUMEHEX_DEMKOWICZ_FAMILY_HDIV(order, order, order);

    volFamily.resize(3, 3 * nb_dofs_family * nb_gauss_pts);
    diffVolFamily.resize(3, 9 * nb_dofs_family * nb_gauss_pts);
    if (nb_dofs_family) {

      std::array<double *, 3> family_ptr = {&volFamily(0, 0), &volFamily(1, 0),
                                            &volFamily(2, 0)};
      std::array<double *, 3> diff_family_ptr = {
          &diffVolFamily(0, 0), &diffVolFamily(1, 0), &diffVolFamily(2, 0)};

      std::array<int, 3> p{order, order, order};
      DemkowiczHexAndQuad::Hdiv_InteriorShapeFunctions_ONHEX(
          p.data(), &*copy_base_fun.data().begin(),
          &*copy_diff_base_fun.data().begin(), family_ptr.data(),
          diff_family_ptr.data(), nb_gauss_pts);

      const int nb_vol_dofs = NBVOLUMEHEX_DEMKOWICZ_HDIV(order);
      auto &face_n = data.dataOnEntities[MBHEX][0].getN(base);
      auto &diff_face_n = data.dataOnEntities[MBHEX][0].getDiffN(base);
      face_n.resize(nb_gauss_pts, 3 * nb_vol_dofs, false);
      diff_face_n.resize(nb_gauss_pts, 9 * nb_vol_dofs, false);

      auto ptr_f0 = &(volFamily(0, 0));
      auto ptr_f1 = &(volFamily(1, 0));
      auto ptr_f2 = &(volFamily(2, 0));
      double *ptr = &face_n(0, 0);
      for (int n = 0; n != volFamily.size2() / 3; ++n) {
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f0;
          ++ptr;
          ++ptr_f0;
        }
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f1;
          ++ptr;
          ++ptr_f1;
        }
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f2;
          ++ptr;
          ++ptr_f2;
        }
      }

      auto diff_ptr_f0 = &(diffVolFamily(0, 0));
      auto diff_ptr_f1 = &(diffVolFamily(1, 0));
      auto diff_ptr_f2 = &(diffVolFamily(2, 0));
      double *diff_ptr = &diff_face_n(0, 0);
      for (int n = 0; n != diffVolFamily.size2() / 9; ++n) {
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f0;
          ++diff_ptr;
          ++diff_ptr_f0;
        }
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f1;
          ++diff_ptr;
          ++diff_ptr_f1;
        }
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f2;
          ++diff_ptr;
          ++diff_ptr_f2;
        }
      }
    } else {
      data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBHEX][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                          false);
    }

  } else {
    data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBHEX][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode HexPolynomialBase::getValueHcurl(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueHcurlDemkowiczBase(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
HexPolynomialBase::getValueHcurlDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  auto &data = cTx->dAta;
  const auto base = cTx->bAse;
  const auto copy_base = cTx->copyNodeBase;
  const int nb_gauss_pts = pts.size2();

  auto &copy_base_fun = data.dataOnEntities[MBVERTEX][0].getN(copy_base);
  auto &copy_diff_base_fun =
      data.dataOnEntities[MBVERTEX][0].getDiffN(copy_base);

  if (nb_gauss_pts != copy_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  if (nb_gauss_pts != copy_diff_base_fun.size1())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong number of gauss pts");

  // edges
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {
    std::array<int, 12> sense;
    std::array<int, 12> order;
    if (data.dataOnEntities[MBEDGE].size() != 12)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected 12 edges data structures on Hex");

    std::array<double *, 12> hcurl_edge_n;
    std::array<double *, 12> diff_hcurl_edge_n;
    bool sum_nb_dofs = false;

    for (int ee = 0; ee != 12; ++ee) {
      if (data.dataOnEntities[MBEDGE][ee].getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Sense on edge <%d> on Hex not set", ee);

      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getOrder();
      const int nb_dofs = NBEDGE_DEMKOWICZ_HCURL(order[ee]);
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts,
                                                        3 * nb_dofs, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts,
                                                            9 * nb_dofs, false);
      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();

      sum_nb_dofs |= (nb_dofs > 0);
    }

    if (sum_nb_dofs) {
      CHKERR DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONHEX(
          sense.data(), order.data(), &*copy_base_fun.data().begin(),
          &*copy_diff_base_fun.data().begin(), hcurl_edge_n.data(),
          diff_hcurl_edge_n.data(), nb_gauss_pts);
    }

  } else {
    for (int ee = 0; ee != 12; ee++) {
      data.dataOnEntities[MBEDGE][ee].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBEDGE][ee].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  // Quad
  if (data.spacesOnEntities[MBQUAD].test(HCURL)) {

    if (data.dataOnEntities[MBQUAD].size() != 6)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Expected six data structures on Hex");

    std::array<int, 6> order;
    double *face_family_ptr[6][2];
    double *diff_face_family_ptr[6][2];

    bool sum_nb_dofs = false;
    for (int ff = 0; ff != 6; ff++) {
      if (data.dataOnEntities[MBQUAD][ff].getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Sense pn quad <%d> not set", ff);

      order[ff] = data.dataOnEntities[MBQUAD][ff].getOrder();
      if (data.facesNodes.size1() != 6)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected six faces");
      if (data.facesNodes.size2() != 4)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Expected four  nodes on face");

      const int nb_family_dofs =
          NBFACEQUAD_DEMKOWICZ_FAMILY_HCURL(order[ff], order[ff]);
      faceFamily[ff].resize(2, 3 * nb_family_dofs * nb_gauss_pts, false);
      diffFaceFamily[ff].resize(2, 9 * nb_family_dofs * nb_gauss_pts, false);

      if (nb_family_dofs) {
        face_family_ptr[ff][0] = &(faceFamily[ff](0, 0));
        face_family_ptr[ff][1] = &(faceFamily[ff](1, 0));
        diff_face_family_ptr[ff][0] = &(diffFaceFamily[ff](0, 0));
        diff_face_family_ptr[ff][1] = &(diffFaceFamily[ff](1, 0));
      }

      sum_nb_dofs |= (nb_family_dofs > 0);
    }

    if (sum_nb_dofs) {
      DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONHEX(
          &*data.facesNodes.data().begin(),
          &*data.facesNodesOrder.data().begin(), order.data(),
          &*copy_base_fun.data().begin(), &*copy_diff_base_fun.data().begin(),
          face_family_ptr, diff_face_family_ptr, nb_gauss_pts);

      for (int ff = 0; ff != 6; ++ff) {

        // put family back

        int nb_dofs = NBFACEQUAD_DEMKOWICZ_HCURL(order[ff]);
        if (nb_dofs) {
          auto &face_n = data.dataOnEntities[MBQUAD][ff].getN(base);
          auto &diff_face_n = data.dataOnEntities[MBQUAD][ff].getDiffN(base);
          face_n.resize(nb_gauss_pts, 3 * nb_dofs, false);
          diff_face_n.resize(nb_gauss_pts, 9 * nb_dofs, false);

          auto ptr_f0 = &(faceFamily[ff](0, 0));
          auto ptr_f1 = &(faceFamily[ff](1, 0));
          double *ptr = &face_n(0, 0);
          for (int n = 0; n != faceFamily[ff].size2() / 3; ++n) {
            for (int j = 0; j != 3; ++j) {
              *ptr = *ptr_f0;
              ++ptr;
              ++ptr_f0;
            }
            for (int j = 0; j != 3; ++j) {
              *ptr = *ptr_f1;
              ++ptr;
              ++ptr_f1;
            }
          }

          auto diff_ptr_f0 = &(diffFaceFamily[ff](0, 0));
          auto diff_ptr_f1 = &(diffFaceFamily[ff](1, 0));
          double *diff_ptr = &diff_face_n(0, 0);
          for (int n = 0; n != diffFaceFamily[ff].size2() / 9; ++n) {
            for (int j = 0; j != 9; ++j) {
              *diff_ptr = *diff_ptr_f0;
              ++diff_ptr;
              ++diff_ptr_f0;
            }
            for (int j = 0; j != 9; ++j) {
              *diff_ptr = *diff_ptr_f1;
              ++diff_ptr;
              ++diff_ptr_f1;
            }
          }
        }
      }
    } else {
      for (int ff = 0; ff != 6; ff++) {
        data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, 0,
                                                          false);
        data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                              false);
      }
    }

  } else {
    for (int ff = 0; ff != 6; ff++) {
      data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(nb_gauss_pts, 0,
                                                            false);
    }
  }

  // Hex
  if (data.spacesOnEntities[MBHEX].test(HCURL)) {

    const int order = data.dataOnEntities[MBHEX][0].getOrder();
    const int nb_dofs = NBVOLUMEHEX_DEMKOWICZ_FAMILY_HCURL(order, order, order);

    volFamily.resize(3, 3 * nb_dofs * nb_gauss_pts);
    diffVolFamily.resize(3, 9 * nb_dofs * nb_gauss_pts);
    if (nb_dofs) {

      std::array<double *, 3> family_ptr = {&volFamily(0, 0), &volFamily(1, 0),
                                            &volFamily(2, 0)};
      std::array<double *, 3> diff_family_ptr = {
          &diffVolFamily(0, 0), &diffVolFamily(1, 0), &diffVolFamily(2, 0)};

      std::array<int, 3> p{order, order, order};
      DemkowiczHexAndQuad::Hcurl_InteriorShapeFunctions_ONHEX(
          p.data(), &*copy_base_fun.data().begin(),
          &*copy_diff_base_fun.data().begin(), family_ptr.data(),
          diff_family_ptr.data(), nb_gauss_pts);

      const int nb_vol_dofs = NBVOLUMEHEX_DEMKOWICZ_HCURL(order);
      auto &face_n = data.dataOnEntities[MBHEX][0].getN(base);
      auto &diff_face_n = data.dataOnEntities[MBHEX][0].getDiffN(base);
      face_n.resize(nb_gauss_pts, 3 * nb_vol_dofs, false);
      diff_face_n.resize(nb_gauss_pts, 9 * nb_vol_dofs, false);

      auto ptr_f0 = &(volFamily(0, 0));
      auto ptr_f1 = &(volFamily(1, 0));
      auto ptr_f2 = &(volFamily(2, 0));
      double *ptr = &face_n(0, 0);
      for (int n = 0; n != volFamily.size2() / 3; ++n) {
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f0;
          ++ptr;
          ++ptr_f0;
        }
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f1;
          ++ptr;
          ++ptr_f1;
        }
        for (int j = 0; j != 3; ++j) {
          *ptr = *ptr_f2;
          ++ptr;
          ++ptr_f2;
        }
      }

      auto diff_ptr_f0 = &(diffVolFamily(0, 0));
      auto diff_ptr_f1 = &(diffVolFamily(1, 0));
      auto diff_ptr_f2 = &(diffVolFamily(2, 0));
      double *diff_ptr = &diff_face_n(0, 0);
      for (int n = 0; n != diffVolFamily.size2() / 9; ++n) {
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f0;
          ++diff_ptr;
          ++diff_ptr_f0;
        }
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f1;
          ++diff_ptr;
          ++diff_ptr_f1;
        }
        for (int j = 0; j != 9; ++j) {
          *diff_ptr = *diff_ptr_f2;
          ++diff_ptr;
          ++diff_ptr_f2;
        }
      }
    } else {
      data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, 0, false);
      data.dataOnEntities[MBHEX][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                          false);
    }

  } else {
    data.dataOnEntities[MBHEX][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBHEX][0].getDiffN(base).resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
HexPolynomialBase::getValue(MatrixDouble &pts,
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
