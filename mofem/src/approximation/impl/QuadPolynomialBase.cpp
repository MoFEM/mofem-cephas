/** \file QuadPolynomialBase.cpp
\brief Implementation of bases on a quad face

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


MoFEMErrorCode QuadPolynomialBase::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {
  MoFEMFunctionBegin;
  *iface = NULL;
  if (uuid == IDD_QUAD_BASE_FUNCTION) {
    *iface = const_cast<QuadPolynomialBase *>(this);
    MoFEMFunctionReturnHot(0);
  } else
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");

  CHKERR BaseFunction::query_interface(uuid, iface);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    CHKERR getValueH1AinsworthBase(pts);
    break;
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueH1DemkowiczBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueH1AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  if (cTx->basePolynomialsType0 == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Polynomial type not set");
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();
  auto &vert_dat = data.dataOnEntities[MBVERTEX][0];

  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    // edges
    if (data.dataOnEntities[MBEDGE].size() != 4)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "should be four edges on quad");

    int sense[4], order[4];
    double *H1edgeN[4], *diffH1edgeN[4];
    for (int ee = 0; ee != 4; ++ee) {
      auto &ent_dat = data.dataOnEntities[MBEDGE][ee];
      if (ent_dat.getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "sense not set");

      sense[ee] = ent_dat.getSense();
      order[ee] = ent_dat.getDataOrder();
      int nb_dofs = NBEDGE_H1(ent_dat.getDataOrder());
      ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
      ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);
      H1edgeN[ee] = &*ent_dat.getN(base).data().begin();
      diffH1edgeN[ee] = &*ent_dat.getDiffN(base).data().begin();
    }
    CHKERR H1_EdgeShapeFunctions_MBQUAD(
        sense, order, &*vert_dat.getN(base).data().begin(),
        &*vert_dat.getDiffN(base).data().begin(), H1edgeN, diffH1edgeN,
        nb_gauss_pts, base_polynomials);
  }

  if (data.spacesOnEntities[MBQUAD].test(H1)) {

    // face
    if (data.dataOnEntities[MBQUAD].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "should be one quad to store bubble base on quad");

    auto &ent_dat = data.dataOnEntities[MBQUAD][0];
    int nb_dofs = NBFACEQUAD_H1(ent_dat.getDataOrder());
    ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
    ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);
    int face_nodes[] = {0, 1, 2, 3};
    CHKERR H1_QuadShapeFunctions_MBQUAD(
        face_nodes, ent_dat.getDataOrder(),
        &*vert_dat.getN(base).data().begin(),
        &*vert_dat.getDiffN(base).data().begin(),
        &*ent_dat.getN(base).data().begin(),
        &*ent_dat.getDiffN(base).data().begin(), nb_gauss_pts,
        base_polynomials);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueH1DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();
  auto &vert_dat = data.dataOnEntities[MBVERTEX][0];

  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    // edges
    if (data.dataOnEntities[MBEDGE].size() != 4)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "should be four edges on quad");

    int sense[4], order[4];
    double *H1edgeN[4], *diffH1edgeN[4];
    for (int ee = 0; ee != 4; ++ee) {
      auto &ent_dat = data.dataOnEntities[MBEDGE][ee];
      if (ent_dat.getSense() == 0)
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "sense not set");

      sense[ee] = ent_dat.getSense();
      order[ee] = ent_dat.getDataOrder();
      int nb_dofs = NBEDGE_H1(ent_dat.getDataOrder());
      ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
      ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);
      H1edgeN[ee] = &*ent_dat.getN(base).data().begin();
      diffH1edgeN[ee] = &*ent_dat.getDiffN(base).data().begin();
    }
    CHKERR DemkowiczHexAndQuad::H1_EdgeShapeFunctions_ONQUAD(
        sense, order, &*vert_dat.getN(base).data().begin(),
        &*vert_dat.getDiffN(base).data().begin(), H1edgeN, diffH1edgeN,
        nb_gauss_pts);
  }

  if (data.spacesOnEntities[MBQUAD].test(H1)) {

    // face
    if (data.dataOnEntities[MBQUAD].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "should be one quad to store bubble base on quad");

    auto &ent_dat = data.dataOnEntities[MBQUAD][0];
    int nb_dofs = NBFACEQUAD_H1(ent_dat.getDataOrder());
    int p = ent_dat.getDataOrder();
    int order[2] = {p, p};
    ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
    ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);

    int face_nodes[] = {0, 1, 2, 3};
    CHKERR DemkowiczHexAndQuad::H1_FaceShapeFunctions_ONQUAD(
        face_nodes, order, &*vert_dat.getN(base).data().begin(),
        &*vert_dat.getDiffN(base).data().begin(),
        &*ent_dat.getN(base).data().begin(),
        &*ent_dat.getDiffN(base).data().begin(), nb_gauss_pts);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
    CHKERR getValueL2DemkowiczBase(pts);
    break;
  case AINSWORTH_LOBATTO_BASE:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Ainsworth not implemented on quad for L2 base");
    break;
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueL2DemkowiczBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueL2DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();
  auto &vert_dat = data.dataOnEntities[MBVERTEX][0];

  auto &ent_dat = data.dataOnEntities[MBQUAD][0];
  int p = ent_dat.getDataOrder();
  int order[2] = {p, p};
  int nb_dofs = NBFACEQUAD_L2(p);
  ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
  ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);

  CHKERR DemkowiczHexAndQuad::L2_FaceShapeFunctions_ONQUAD(
      order, &*vert_dat.getN(base).data().begin(),
      &*vert_dat.getDiffN(base).data().begin(),
      &*ent_dat.getN(base).data().begin(),
      &*ent_dat.getDiffN(base).data().begin(), nb_gauss_pts);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueHcurl(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Ainsworth not implemented on quad for Hcurl base");
    break;
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueHcurlDemkowiczBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
QuadPolynomialBase::getValueHcurlDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();

  // Calculation H-curl on quad edges
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {

    if (data.dataOnEntities[MBEDGE].size() != 4)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong number of edges on quad, should be 4 but is %d",
               data.dataOnEntities[MBEDGE].size());

    int sense[4], order[4];
    double *hcurl_edge_n[4];
    double *diff_hcurl_edge_n[4];

    for (int ee = 0; ee != 4; ++ee) {

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
          nb_gauss_pts, 3 * 2 * nb_dofs, false);

      hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getN(base).data().begin();
      diff_hcurl_edge_n[ee] =
          &*data.dataOnEntities[MBEDGE][ee].getDiffN(base).data().begin();
    }
    CHKERR DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONQUAD(
        sense, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts);
  }

  if (data.spacesOnEntities[MBQUAD].test(HCURL)) {

    // face
    if (data.dataOnEntities[MBQUAD].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No data struture to keep base functions on face");

    int p = data.dataOnEntities[MBQUAD][0].getDataOrder();

    MatrixDouble face_family(
        2, 3 * NBFACEQUAD_DEMKOWICZ_FAMILY_QUAD_HCURL(p, p) * nb_gauss_pts);
    MatrixDouble diff_face_family(
        2, 3 * 2 * NBFACEQUAD_DEMKOWICZ_FAMILY_QUAD_HCURL(p, p) * nb_gauss_pts);

    int order[2] = {p, p};
    double *face_family_ptr[] = {&face_family(0, 0), &face_family(1, 0)};
    double *diff_face_family_ptr[] = {&diff_face_family(0, 0),
                                      &diff_face_family(1, 0)};
    int face_nodes[] = {0, 1, 2, 3};
    CHKERR DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONQUAD(
        face_nodes, order,
        &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        face_family_ptr, diff_face_family_ptr, nb_gauss_pts);

    // put family back

    int nb_dofs = NBFACEQUAD_DEMKOWICZ_HCURL(p);
    auto &face_n = data.dataOnEntities[MBQUAD][0].getN(base);
    auto &diff_face_n = data.dataOnEntities[MBQUAD][0].getDiffN(base);
    face_n.resize(nb_gauss_pts, 3 * nb_dofs, false);
    diff_face_n.resize(nb_gauss_pts, 3 * 2 * nb_dofs, false);
    
    double *ptr_f0 = &face_family(0, 0);
    double *ptr_f1 = &face_family(1, 0);
    double *ptr = &face_n(0, 0);
    for (int n = 0; n != face_family.size2() / 3; ++n) {
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

    double *diff_ptr_f0 = &diff_face_family(0, 0);
    double *diff_ptr_f1 = &diff_face_family(1, 0);
    double *diff_ptr = &diff_face_n(0, 0);
    for (int n = 0; n != diff_face_family.size2() / 6; ++n) {
      for (int j = 0; j != 6; ++j) {
        *diff_ptr = *diff_ptr_f0;
        ++diff_ptr;
        ++diff_ptr_f0;
      }
      for (int j = 0; j != 6; ++j) {
        *diff_ptr = *diff_ptr_f1;
        ++diff_ptr;
        ++diff_ptr_f1;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode QuadPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
  case AINSWORTH_LOBATTO_BASE:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Ainsworth not implemented on quad for Hdiv base");
    break;
  case DEMKOWICZ_JACOBI_BASE:
    CHKERR getValueHdivDemkowiczBase(pts);
    break;
  case AINSWORTH_BERNSTEIN_BEZIER_BASE:
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
QuadPolynomialBase::getValueHdivDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
QuadPolynomialBase::getValue(MatrixDouble &pts,
                             boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;

  BaseFunctionUnknownInterface *iface;
  CHKERR ctx_ptr->query_interface(IDD_QUAD_BASE_FUNCTION, &iface);
  cTx = reinterpret_cast<EntPolynomialBaseCtx *>(iface);

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 2)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSourcesCore &data = cTx->dAta;
  if (cTx->copyNodeBase != NOBASE)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Shape base has to be on NOBASE", ApproximationBaseNames[base]);

  data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
      data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);
  data.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(base) =
      data.dataOnEntities[MBVERTEX][0].getDiffNSharedPtr(cTx->copyNodeBase);

  auto &base_shape = data.dataOnEntities[MBVERTEX][0].getN(base);
  auto &diff_base = data.dataOnEntities[MBVERTEX][0].getDiffN(base);

  if (base_shape.size1() != pts.size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of base functions integration points not equal number of "
            "set integration point");
  if (base_shape.size2() != 4)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of shape functions should be four");
  if (diff_base.size1() != pts.size2())
    SETERRQ2(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Number of diff base functions integration points not equal number of "
        "set integration point %d != %d",
        diff_base.size1(), pts.size2());
  if (diff_base.size2() != 8)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of shape functions should be four");

  switch (cTx->sPace) {
  case H1:
    CHKERR getValueH1(pts);
    break;
  case L2:
    CHKERR getValueL2(pts);
    break;
  case HCURL:
    CHKERR getValueHcurl(pts);
    break;
  case HDIV:
    CHKERR getValueHdiv(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  MoFEMFunctionReturn(0);
}
