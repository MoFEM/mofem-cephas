/** \file FlatPrismPolynomialBase.cpp
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

MoFEMErrorCode FlatPrismPolynomialBaseCtx::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<FlatPrismPolynomialBaseCtx *>(this);
  return 0;
}

FlatPrismPolynomialBaseCtx::FlatPrismPolynomialBaseCtx(
    DataForcesAndSourcesCore &data, moab::Interface &moab,
    const NumeredEntFiniteElement *fe_ptr, const FieldSpace space,
    const FieldApproximationBase base,
    const FieldApproximationBase copy_node_base)
    : EntPolynomialBaseCtx(data, space, base, copy_node_base), mOab(moab),
      fePtr(fe_ptr) {
  CHKERR setBase();
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}
FlatPrismPolynomialBaseCtx::~FlatPrismPolynomialBaseCtx() {}

MoFEMErrorCode FlatPrismPolynomialBase::query_interface(
    boost::typeindex::type_index type_index, UnknownInterface **iface) const {
  *iface = const_cast<FlatPrismPolynomialBase *>(this);
  return 0;
}

FlatPrismPolynomialBase::~FlatPrismPolynomialBase() {}
FlatPrismPolynomialBase::FlatPrismPolynomialBase() {}

MoFEMErrorCode
FlatPrismPolynomialBase::getValue(MatrixDouble &pts,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr) {
  MoFEMFunctionBeginHot;

  cTx = ctx_ptr->getInterface<FlatPrismPolynomialBaseCtx>();

  if (!cTx->fePtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer to element should be given "
            "when EntPolynomialBaseCtx is constructed "
            "(use different constructor)");

  int nb_gauss_pts = pts.size2();
  if (!nb_gauss_pts)
    MoFEMFunctionReturnHot(0);

  if (pts.size1() < 1)
    SETERRQ(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "Wrong dimension of pts, should be at least 3 rows with coordinates");

  const FieldApproximationBase base = cTx->bAse;
  DataForcesAndSourcesCore &data = cTx->dAta;

  if (cTx->copyNodeBase == LASTBASE)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  else
    data.dataOnEntities[MBVERTEX][0].getNSharedPtr(base) =
        data.dataOnEntities[MBVERTEX][0].getNSharedPtr(cTx->copyNodeBase);

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 6, false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts, 12,
                                                         false);
  if (data.dataOnEntities[MBVERTEX][0].getN(base).size1() !=
      (unsigned int)nb_gauss_pts) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Base functions or nodes has wrong number of integration points "
             "for base %s",
             ApproximationBaseNames[base]);
  }

  data.dataOnEntities[MBVERTEX][0].getN(base).resize(nb_gauss_pts, 6, false);
  data.dataOnEntities[MBVERTEX][0].getDiffN(base).resize(nb_gauss_pts, 12,
                                                         false);
  N.resize(nb_gauss_pts, 3, false);
  diffN.resize(3, 2, false);
  CHKERR ShapeMBTRI(&*N.data().begin(), &pts(0, 0), &pts(1, 0), nb_gauss_pts);
  std::copy(Tools::diffShapeFunMBTRI.begin(), Tools::diffShapeFunMBTRI.end(),
            &*diffN.data().begin());

  // This is needed to have proper order of nodes on faces
  CHKERR cTx->mOab.get_connectivity(cTx->fePtr->getEnt(), connPrism, numNodes,
                                    true);
  SideNumber_multiIndex &side_table =
      const_cast<SideNumber_multiIndex &>(cTx->fePtr->getSideNumberTable());
  SideNumber_multiIndex::nth_index<1>::type::iterator siit3 =
      side_table.get<1>().find(boost::make_tuple(MBTRI, 3));
  SideNumber_multiIndex::nth_index<1>::type::iterator siit4 =
      side_table.get<1>().find(boost::make_tuple(MBTRI, 4));
  if (siit3 == side_table.get<1>().end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  if (siit4 == side_table.get<1>().end())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  CHKERR cTx->mOab.get_connectivity(siit3->get()->ent, connFace3, numNodes,
                                    true);
  CHKERR cTx->mOab.get_connectivity(siit4->get()->ent, connFace4, numNodes,
                                    true);

  for (int nn = 0; nn < 3; nn++) {
    faceNodes[0][nn] = std::distance(
        connPrism, std::find(connPrism, connPrism + 3, connFace3[nn]));
    faceNodes[1][nn] = std::distance(
        connPrism + 3, std::find(connPrism + 3, connPrism + 6, connFace4[nn]));
    for (int gg = 0; gg < nb_gauss_pts; gg++) {
      double val = N(gg, nn);
      double val_x = diffN(nn, 0);
      double val_y = diffN(nn, 1);
      data.dataOnEntities[MBVERTEX][0].getN(base)(gg, nn) = val;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg, 2 * nn + 0) = val_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg, 2 * nn + 1) = val_y;
      data.dataOnEntities[MBVERTEX][0].getN(base)(gg, 3 + nn) = val;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg, 6 + 2 * nn + 0) =
          val_x;
      data.dataOnEntities[MBVERTEX][0].getDiffN(base)(gg, 6 + 2 * nn + 1) =
          val_y;
    }
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
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FlatPrismPolynomialBase::getValueH1(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;

  DataForcesAndSourcesCore &data = cTx->dAta;
  const FieldApproximationBase base = cTx->bAse;
  PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s, double *L,
                                     double *diffL, const int dim) =
      cTx->basePolynomialsType0;

  int nb_gauss_pts = pts.size2();

  // edges
  if (data.dataOnEntities[MBEDGE].size() != 9)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  int sense[9], order[9];
  double *H1edgeN[9], *diffH1edgeN[9];

  auto set_edge_base_data = [&](const int ee) {
    MoFEMFunctionBegin;
    auto &ent_data = data.dataOnEntities[MBEDGE][ee];
    if (ent_data.getSense() == 0)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    sense[ee] = ent_data.getSense();
    order[ee] = ent_data.getDataOrder();
    int nb_dofs = NBEDGE_H1(ent_data.getDataOrder());
    ent_data.getN(base).resize(nb_gauss_pts, nb_dofs, false);
    ent_data.getDiffN(base).resize(nb_gauss_pts, 2 * nb_dofs, false);
    H1edgeN[ee] = &*ent_data.getN(base).data().begin();
    diffH1edgeN[ee] = &*ent_data.getDiffN(base).data().begin();
    MoFEMFunctionReturn(0);
  };

  if ((data.spacesOnEntities[MBEDGE]).test(H1)) {

    for (int ee = 0; ee != 3; ++ee)
      CHKERR set_edge_base_data(ee);
    for (int ee = 6; ee != 9; ++ee)
      CHKERR set_edge_base_data(ee);

    // shape functions on face 3
    CHKERR H1_EdgeShapeFunctions_MBTRI(
        &sense[0], &order[0], &*N.data().begin(), &*diffN.data().begin(),
        &H1edgeN[0], &diffH1edgeN[0], nb_gauss_pts, base_polynomials);
    // shape functions on face 4
    CHKERR H1_EdgeShapeFunctions_MBTRI(
        &sense[6], &order[6], &*N.data().begin(), &*diffN.data().begin(),
        &H1edgeN[6], &diffH1edgeN[6], nb_gauss_pts, base_polynomials);
  }

  // face
  if (data.dataOnEntities[MBTRI].size() != 5)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");

  if ((data.spacesOnEntities[MBTRI]).test(H1)) {
    for (int ff = 3; ff <= 4; ff++) {
      int nb_dofs = NBFACETRI_H1(data.dataOnEntities[MBTRI][ff].getDataOrder());
      data.dataOnEntities[MBTRI][ff].getN(base).resize(nb_gauss_pts, nb_dofs,
                                                       false);
      data.dataOnEntities[MBTRI][ff].getDiffN(base).resize(nb_gauss_pts,
                                                           2 * nb_dofs, false);
      CHKERR H1_FaceShapeFunctions_MBTRI(
          faceNodes[ff - 3], data.dataOnEntities[MBTRI][ff].getDataOrder(),
          &*N.data().begin(), &*diffN.data().begin(),
          &*data.dataOnEntities[MBTRI][ff].getN(base).data().begin(),
          &*data.dataOnEntities[MBTRI][ff].getDiffN(base).data().begin(),
          nb_gauss_pts, base_polynomials);
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FlatPrismPolynomialBase::getValueL2(MatrixDouble &pts) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FlatPrismPolynomialBase::getValueHdiv(MatrixDouble &pts) {
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
}

MoFEMErrorCode FlatPrismPolynomialBase::getValueHcurl(MatrixDouble &pts) {
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not yet implemented");
}
