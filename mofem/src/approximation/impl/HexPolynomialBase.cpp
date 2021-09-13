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
    data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts,
                                                       copy_base_fun.size2());
    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(
        nb_gauss_pts, copy_diff_base_fun.size2());
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
        order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();

        int nb_dofs = NBEDGE_H1(data.dataOnEntities[MBEDGE][ee].getDataOrder());
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

        order[ff] = data.dataOnEntities[MBQUAD][ff].getDataOrder();
        const int nb_dofs = NBFACEQUAD_H1(order[ff]);

        data.dataOnEntities[MBQUAD][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                          false);
        data.dataOnEntities[MBQUAD][ff].getDiffN(base).resize(
            nb_gauss_pts, 3 * nb_dofs, false);

        h1_face_n[ff] =
            &*data.dataOnEntities[MBQUAD][ff].getN(base).data().begin();
        diff_h1_face_n[ff] =
            &*data.dataOnEntities[MBQUAD][ff].getDiffN(base).data().begin();

        nb_dofs_sum |= (nb_dofs_sum > 0);
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
      int order = data.dataOnEntities[MBHEX][0].getDataOrder();
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
  const int vol_order = data.dataOnEntities[MBHEX][0].getDataOrder();

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

  // DataForcesAndSourcesCore &data = cTx->dAta;
  // const FieldApproximationBase base = cTx->bAse;
  // if (base != DEMKOWICZ_JACOBI_BASE) {
  //   SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
  //            "This should be used only with DEMKOWICZ_JACOBI_BASE "
  //            "but base is %s",
  //            ApproximationBaseNames[base]);
  // }
  // int nb_gauss_pts = pts.size2();

  // int volume_order = data.dataOnEntities[MBTET][0].getDataOrder();

  // int p_f[4];
  // double *phi_f[4];
  // double *diff_phi_f[4];

  // // Calculate base function on tet faces
  // for (int ff = 0; ff != 4; ff++) {
  //   int face_order = data.dataOnEntities[MBTRI][ff].getDataOrder();
  //   int order = volume_order > face_order ? volume_order : face_order;
  //   data.dataOnEntities[MBTRI][ff].getN(base).resize(
  //       nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
  //   data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
  //       nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(order), false);
  //   p_f[ff] = order;
  //   phi_f[ff] = &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin();
  //   diff_phi_f[ff] =
  //       &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin();
  //   if (NBFACETRI_DEMKOWICZ_HDIV(order) == 0)
  //     continue;
  //   CHKERR Hdiv_Demkowicz_Face_MBTET_ON_FACE(
  //       &data.facesNodes(ff, 0), order,
  //       &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
  //       &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
  //       phi_f[ff], diff_phi_f[ff], nb_gauss_pts, 4);
  // }

  // // Calculate base functions in tet interior
  // if (NBVOLUMETET_DEMKOWICZ_HDIV(volume_order) > 0) {
  //   data.dataOnEntities[MBTET][0].getN(base).resize(
  //       nb_gauss_pts, 3 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);
  //   data.dataOnEntities[MBTET][0].getDiffN(base).resize(
  //       nb_gauss_pts, 9 * NBVOLUMETET_DEMKOWICZ_HDIV(volume_order), false);
  //   double *phi_v =
  //   &*data.dataOnEntities[MBTET][0].getN(base).data().begin(); double
  //   *diff_phi_v =
  //       &*data.dataOnEntities[MBTET][0].getDiffN(base).data().begin();
  //   CHKERR Hdiv_Demkowicz_Interior_MBTET(
  //       volume_order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
  //       &data.dataOnEntities[MBVERTEX][0].getDiffN(base)(0, 0), p_f, phi_f,
  //       diff_phi_f, phi_v, diff_phi_v, nb_gauss_pts);
  // }

  // // Set size of face base correctly
  // for (int ff = 0; ff != 4; ff++) {
  //   int face_order = data.dataOnEntities[MBTRI][ff].getDataOrder();
  //   data.dataOnEntities[MBTRI][ff].getN(base).resize(
  //       nb_gauss_pts, 3 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
  //   data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(
  //       nb_gauss_pts, 9 * NBFACETRI_DEMKOWICZ_HDIV(face_order), true);
  // }

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
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Sense on edge <%d> on Hex not set", ee);

      sense[ee] = data.dataOnEntities[MBEDGE][ee].getSense();
      order[ee] = data.dataOnEntities[MBEDGE][ee].getDataOrder();
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
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Sense pn quad <%d> not set", ff);

      order[ff] = data.dataOnEntities[MBQUAD][ff].getDataOrder();
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

    const int order = data.dataOnEntities[MBHEX][0].getDataOrder();
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
