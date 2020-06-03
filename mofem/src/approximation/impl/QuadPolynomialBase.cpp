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

QuadPolynomialBase::QuadPolynomialBase() {}
QuadPolynomialBase::~QuadPolynomialBase() {}

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
    CHKERR getValueH1AinsworthBase(pts);
    break;
  case AINSWORTH_LOBATTO_BASE:
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
    int pp[2] = {order[0], order[1]};
    CHKERR H1_EdgeShapeFunctions_ONQUAD(
        sense, pp, &*vert_dat.getN(base).data().begin(), H1edgeN, diffH1edgeN,
        nb_gauss_pts);
  }

  if (data.spacesOnEntities[MBQUAD].test(H1)) {

    // face
    if (data.dataOnEntities[MBQUAD].size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "should be one quad to store bubble base on quad");

    auto &ent_dat = data.dataOnEntities[MBQUAD][0];
    int nb_dofs = NBFACEQUAD_FULL_H1(ent_dat.getDataOrder());
    int p = ent_dat.getDataOrder();
    int order[2] = {p, p};
    ent_dat.getN(base).resize(nb_gauss_pts, nb_dofs, false);
    ent_dat.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);

    CHKERR H1_FaceShapeFunctions_ONQUAD(order,
        &*vert_dat.getN(base).data().begin(),
        &*ent_dat.getN(base).data().begin(),
        &*ent_dat.getDiffN(base).data().begin(), nb_gauss_pts);
  }

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
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Number of diff base functions integration points not equal number of "
        "set integration point");
  if (diff_base.size2() != 8)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Number of shape functions should be four");

  switch (cTx->sPace) {
  case H1:
    CHKERR getValueH1(pts);
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  MoFEMFunctionReturn(0);
}
