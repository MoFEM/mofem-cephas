/** \file EdgePolynomialBase.cpp
\brief Implementation of Ainsworth-Cole H1 base on edge
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

MoFEMErrorCode EdgePolynomialBase::query_interface(
    const MOFEMuuid &uuid, BaseFunctionUnknownInterface **iface) const {

  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_EDGE_BASE_FUNCTION) {
    *iface = const_cast<EdgePolynomialBase *>(this);
    MoFEMFunctionReturnHot(0);
  } else {
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "wrong interference");
  }
  ierr = BaseFunction::query_interface(uuid, iface);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

EdgePolynomialBase::~EdgePolynomialBase() {}
EdgePolynomialBase::EdgePolynomialBase() {}

MoFEMErrorCode
EdgePolynomialBase::getValue(MatrixDouble &pts,
                             boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBegin;

  BaseFunctionUnknownInterface *iface;
  CHKERR ctx_ptr->query_interface(IDD_EDGE_BASE_FUNCTION, &iface);
  cTx = reinterpret_cast<EntPolynomialBaseCtx *>(iface);

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 1)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSourcesCore &data = cTx->dAta;

  if (base != AINSWORTH_BERNSTEIN_BEZIER_BASE) {
    if (cTx->copyNodeBase == LASTBASE) {
      data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 2,
                                                         false);
      CHKERR ShapeMBEDGE(
          &*data.dataOnEntities[MBVERTEX][0].getN(base).data().begin(),
          &pts(0, 0), nb_gauss_pts);
    } else
      data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
          data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);

    if (data.dataOnEntities[MBVERTEX][0].getN(base).size1() !=
        (unsigned int)nb_gauss_pts)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Base functions or nodes has wrong number of integration points "
               "for base %s",
               ApproximationBaseNames[base]);

    data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(2, 1, false);
    std::copy(Tools::diffShapeFunMBEDGE.begin(),
              Tools::diffShapeFunMBEDGE.end(),
              &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin());
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
    CHKERR getValueL2DemkowiczBase(pts);
    // CHKERR getValueL2(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBegin;




  switch (cTx->bAse) {
  case AINSWORTH_LEGENDRE_BASE:
    CHKERR getValueH1DemkowiczBase(pts);
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
EdgePolynomialBase::getValueH1BernsteinBezierBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const std::string &field_name = cTx->fieldName;
  const int nb_gauss_pts = pts.size2();

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

  auto &get_n = get_base(data.dataOnEntities[MBVERTEX][0]);
  auto &get_diff_n = get_diff_base(data.dataOnEntities[MBVERTEX][0]);
  get_n.resize(nb_gauss_pts, 2, false);
  get_diff_n.resize(nb_gauss_pts, 2, false);
  get_n.clear();
  get_diff_n.clear();

  if (data.dataOnEntities[MBVERTEX][0].getN(NOBASE).size1() !=
      (unsigned int)nb_gauss_pts)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[NOBASE]);
  auto &lambda = data.dataOnEntities[MBVERTEX][0].getN(NOBASE);

  auto &vertex_alpha = get_alpha(data.dataOnEntities[MBVERTEX][0]);
  vertex_alpha.resize(2, 2, false);
  vertex_alpha(0, 0) = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[0];
  vertex_alpha(0, 1) = 0;
  vertex_alpha(1, 0) = 0;
  vertex_alpha(1, 1) = data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[1];

  CHKERR BernsteinBezier::baseFunctionsEdge(
      1, nb_gauss_pts, vertex_alpha.size1(), &vertex_alpha(0, 0), &lambda(0, 0),
      Tools::diffShapeFunMBEDGE.data(), &get_n(0, 0), &get_diff_n(0, 0));
  std::array<double, 2> f = {
      boost::math::factorial<double>(
          data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[0]),
      boost::math::factorial<double>(
          data.dataOnEntities[MBVERTEX][0].getBBNodeOrder()[1])};

  for (int g = 0; g != nb_gauss_pts; ++g)
    for (int n = 0; n != 2; ++n) {
      get_n(g, n) *= f[n];
      get_diff_n(g, n) *= f[n];
    }

  if (data.spacesOnEntities[MBEDGE].test(H1)) {
    if (data.dataOnEntities[MBEDGE].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Wrong size ent of ent data");

    int order = data.dataOnEntities[MBEDGE][0].getDataOrder();
    const int nb_dofs = NBEDGE_H1(order);

    auto &get_n = get_base(data.dataOnEntities[MBEDGE][0]);
    auto &get_diff_n = get_diff_base(data.dataOnEntities[MBEDGE][0]);
    get_n.resize(nb_gauss_pts, nb_dofs, false);
    get_diff_n.resize(nb_gauss_pts, nb_dofs, false);

    if (nb_dofs) {
      auto &edge_alpha = get_alpha(data.dataOnEntities[MBEDGE][0]);
      edge_alpha.resize(nb_dofs, 2);
      CHKERR BernsteinBezier::generateIndicesEdgeEdge(order, &edge_alpha(0, 0));
      CHKERR BernsteinBezier::baseFunctionsEdge(
          order, nb_gauss_pts, edge_alpha.size1(), &edge_alpha(0, 0),
          &lambda(0, 0), Tools::diffShapeFunMBEDGE.data(), &get_n(0, 0),
          &get_diff_n(0, 0));
    }
  } else {
    get_n.resize(nb_gauss_pts, 0, false);
    get_diff_n.resize(nb_gauss_pts, 0, false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueH1AinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  const int side_number = 0;
  int sense = data.dataOnEntities[MBEDGE][side_number].getSense();
  int order = data.dataOnEntities[MBEDGE][side_number].getDataOrder();
  data.dataOnEntities[MBEDGE][side_number].getN(base).resize(
      nb_gauss_pts, NBEDGE_H1(order), false);
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).resize(
      nb_gauss_pts, NBEDGE_H1(order), false);

  data.dataOnEntities[MBEDGE][side_number].getN(base).clear();
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).clear();

  L.resize(NBEDGE_H1(order), false);
  diffL.resize(NBEDGE_H1(order), false);

  if (data.dataOnEntities[MBEDGE][side_number].getDataOrder() > 1) {

    double diff_s = 2.; // s = s(xi), ds/dxi = 2., because change of basis
    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      double s = 2 * pts(0, gg) - 1; // makes form -1..1
      if (!sense) {
        s *= -1;
        diff_s *= -1;
      }

      // calculate Legendre polynomials at integration points
      CHKERR base_polynomials(NBEDGE_H1(order) - 1, s, &diff_s,
                              &*L.data().begin(), &*diffL.data().begin(), 1);

      for (unsigned int pp = 0;
           pp < data.dataOnEntities[MBEDGE][side_number].getN(base).size2();
           pp++) {

        // Calculate edge shape functions N0*N1*L(p), where N0 and N1 are nodal
        // shape functions
        double v = data.dataOnEntities[MBVERTEX][0].getN(base)(gg, 0) *
                   data.dataOnEntities[MBVERTEX][0].getN(base)(gg, 1);
        data.dataOnEntities[MBEDGE][side_number].getN(base)(gg, pp) = v * L(pp);

        // Calculate derivative edge shape functions
        // dN/dksi = dN0/dxi*N1*L + N0*dN1/ksi*L + N0*N1*dL/dxi
        data.dataOnEntities[MBEDGE][side_number].getDiffN(base)(gg, pp) =
            ((+1.) * data.dataOnEntities[MBVERTEX][0].getN(base)(gg, 1) +
             data.dataOnEntities[MBVERTEX][0].getN(base)(gg, 0) * (-1.)) *
                L(pp) +
            v * diffL(pp);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueH1DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();

  const int side_number = 0;
  int order = data.dataOnEntities[MBEDGE][side_number].getDataOrder();

  data.dataOnEntities[MBEDGE][side_number].getN(base).resize(
                                 nb_gauss_pts, NBEDGE_H1(order), false);
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).resize(
                                 nb_gauss_pts, NBEDGE_H1(order), false);

  data.dataOnEntities[MBEDGE][side_number].getN(base).clear();
  data.dataOnEntities[MBEDGE][side_number].getDiffN(base).clear();
  
  MatrixDouble Nshape(nb_gauss_pts, 2);
  Nshape.clear();
  if (data.dataOnEntities[MBEDGE][side_number].getDataOrder() > 1) {
    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      double mu0 = 1.0 - pts(0, gg);  // pts ranges over [0, 1]
      double mu1 = pts(0, gg);

      Nshape(gg, 0) = mu0;
      Nshape(gg, 1) = mu1;
    }

    CHKERR DemkowiczHexAndQuad::H1_BubbleShapeFunctions_ONSEGMENT(
        order, &*Nshape.data().begin(),
        &*data.dataOnEntities[MBEDGE][side_number].getN(base).data().begin(),
        &*data.dataOnEntities[MBEDGE][side_number]
              .getDiffN(base)
              .data()
              .begin(),
        nb_gauss_pts);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueL2DemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();

  const int side_number = 0;
  int order = data.dataOnEntities[MBEDGE][side_number].getDataOrder();

  data.dataOnEntities[MBEDGE][side_number].getN(base).resize(
      nb_gauss_pts, NBEDGE_EXACT_L2(order), false);

  data.dataOnEntities[MBEDGE][side_number].getN(base).clear();

  MatrixDouble Nshape(nb_gauss_pts, 2);
  Nshape.clear();
  if (data.dataOnEntities[MBEDGE][side_number].getDataOrder() > 1) {
    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      double mu0 = 1.0 - pts(0, gg); // pts ranges over [0, 1]
      double mu1 = pts(0, gg);

      Nshape(gg, 0) = mu0;
      Nshape(gg, 1) = mu1;
    }

    CHKERR DemkowiczHexAndQuad::L2_ShapeFunctions_ONSEGMENT(
        order, &*Nshape.data().begin(),
        &*data.dataOnEntities[MBEDGE][side_number].getN(base).data().begin(),
        nb_gauss_pts);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "Make no sense, unless problem is 1d (1d not implemented yet)");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueHdiv(MatrixDouble &pts) {
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "Make no sense, unless problem is 2d (2d not implemented yet)");
}

MoFEMErrorCode
EdgePolynomialBase::getValueHcurlAinsworthBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {

    if (data.dataOnEntities[MBEDGE].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

    int sense = data.dataOnEntities[MBEDGE][0].getSense();
    int order = data.dataOnEntities[MBEDGE][0].getDataOrder();
    int nb_dofs =
        NBEDGE_AINSWORTH_HCURL(data.dataOnEntities[MBEDGE][0].getDataOrder());
    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                     false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                         false);

    CHKERR Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_EDGE(
        sense, order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBEDGE][0].getN(base).data().begin(), nullptr,
        nb_gauss_pts, base_polynomials);

  } else {
    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                         false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
EdgePolynomialBase::getValueHcurlDemkowiczBase(MatrixDouble &pts) {
  MoFEMFunctionBegin;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;

  int nb_gauss_pts = pts.size2();
  if (data.spacesOnEntities[MBEDGE].test(HCURL)) {

    if (data.dataOnEntities[MBEDGE].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "No data structure to store base functions");

    int sense = data.dataOnEntities[MBEDGE][0].getSense();
    int order = data.dataOnEntities[MBEDGE][0].getDataOrder();
    int nb_dofs =
        NBEDGE_DEMKOWICZ_HCURL(data.dataOnEntities[MBEDGE][0].getDataOrder());
    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts, 3 * nb_dofs,
                                                     false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                         false);
    CHKERR Hcurl_Demkowicz_EdgeBaseFunctions_MBEDGE(
        sense, order, &data.dataOnEntities[MBVERTEX][0].getN(base)(0, 0),
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(base).data().begin(),
        &*data.dataOnEntities[MBEDGE][0].getN(base).data().begin(), NULL,
        nb_gauss_pts);
  } else {

    data.dataOnEntities[MBEDGE][0].getN(base).resize(nb_gauss_pts, 0, false);
    data.dataOnEntities[MBEDGE][0].getDiffN(base).resize(nb_gauss_pts, 0,
                                                         false);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EdgePolynomialBase::getValueHcurl(MatrixDouble &pts) {
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