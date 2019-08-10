/** \file FaceElementForcesAndSourcesCore.cpp

\brief Implementation of face element

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

#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <lapack_wrap.h>
// #include <gm_rule.h>
#include <quad.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

FaceElementForcesAndSourcesCoreBase::FaceElementForcesAndSourcesCoreBase(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field),
      meshPositionsFieldName("MESH_NODE_POSITIONS"),
      opHOCoordsAndNormals(hoCoordsAtGaussPts, normalsAtGaussPts,
                           tangentOneAtGaussPts, tangentTwoAtGaussPts),
      opContravariantTransform(nOrmal, normalsAtGaussPts),
      opCovariantTransform(nOrmal, normalsAtGaussPts, tangentOne,
                           tangentOneAtGaussPts, tangentTwo,
                           tangentTwoAtGaussPts) {
  getElementPolynomialBase() =
      boost::shared_ptr<BaseFunction>(new TriPolynomialBase());
}

MoFEMErrorCode
FaceElementForcesAndSourcesCoreBase::UserDataOperator::loopSideVolumes(
    const string &fe_name, VolumeElementForcesAndSourcesCoreOnSide &method) {
  MoFEMFunctionBegin;

  const EntityHandle ent = getNumeredEntFiniteElementPtr()->getEnt();
  const Problem *problem_ptr = getFEMethod()->problemPtr;
  Range adjacent_volumes;
  CHKERR getFaceFE()->mField.getInterface<BitRefManager>()->getAdjacenciesAny(
      ent, 3, adjacent_volumes);
  typedef NumeredEntFiniteElement_multiIndex::index<
      Composite_Name_And_Ent_mi_tag>::type FEByComposite;
  FEByComposite &numered_fe =
      problem_ptr->numeredFiniteElements->get<Composite_Name_And_Ent_mi_tag>();

  method.feName = fe_name;

  CHKERR method.setSideFEPtr(getFaceFE());
  CHKERR method.copyBasicMethod(*getFEMethod());
  CHKERR method.copyKsp(*getFEMethod());
  CHKERR method.copySnes(*getFEMethod());
  CHKERR method.copyTs(*getFEMethod());

  CHKERR method.preProcess();

  int nn = 0;
  method.loopSize = adjacent_volumes.size();
  for (Range::iterator vit = adjacent_volumes.begin();
       vit != adjacent_volumes.end(); vit++) {
    FEByComposite::iterator miit =
        numered_fe.find(boost::make_tuple(fe_name, *vit));
    if (miit != numered_fe.end()) {
      method.nInTheLoop = nn++;
      method.numeredEntFiniteElementPtr = *miit;
      method.dataFieldEntsPtr = (*miit)->sPtr->data_field_ents_view;
      method.rowFieldEntsPtr = (*miit)->sPtr->row_field_ents_view;
      method.colFieldEntsPtr = (*miit)->sPtr->col_field_ents_view;
      method.dataPtr = (*miit)->sPtr->data_dofs;
      method.rowPtr = (*miit)->rows_dofs;
      method.colPtr = (*miit)->cols_dofs;
      CHKERR method();
    }
  }

  CHKERR method.postProcess();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCoreBase::calculateAreaAndNormal() {
  MoFEMFunctionBegin;
  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
  coords.resize(num_nodes * 3, false);
  CHKERR mField.get_moab().get_coords(conn, num_nodes, &*coords.data().begin());
  double diff_n[6];
  CHKERR ShapeDiffMBTRI(diff_n);

  nOrmal.resize(3, false);
  tangentOne.resize(3, false);
  tangentTwo.resize(3, false);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_coords(
      &coords[0], &coords[1], &coords[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_normal(
      &nOrmal[0], &nOrmal[1], &nOrmal[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_t1(
      &tangentOne[0], &tangentOne[1], &tangentOne[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_t2(
      &tangentTwo[0], &tangentTwo[1], &tangentTwo[2]);
  FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff(&diff_n[0],
                                                            &diff_n[1]);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;
  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  t_t1(i) = 0;
  t_t2(i) = 0;
  for (int nn = 0; nn != 3; ++nn) {
    t_t1(i) += t_coords(i) * t_diff(N0);
    t_t2(i) += t_coords(i) * t_diff(N1);
    ++t_coords;
    ++t_diff;
  }
  t_normal(j) = FTensor::levi_civita(i, j, k) * t_t1(k) * t_t2(i);
  aRea = sqrt(t_normal(i) * t_normal(i)) / 2.;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCoreBase::setIntegrationPts() {
  MoFEMFunctionBegin;
  // Set integration points
  int order_data = getMaxDataOrder();
  int order_row = getMaxRowOrder();
  int order_col = getMaxColOrder();
  int rule = getRule(order_row, order_col, order_data);

  if (rule >= 0) {
    if (rule < QUAD_2D_TABLE_SIZE) {
      if (QUAD_2D_TABLE[rule]->dim != 2) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
      }
      if (QUAD_2D_TABLE[rule]->order < rule) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", QUAD_2D_TABLE[rule]->order, rule);
      }
      nbGaussPts = QUAD_2D_TABLE[rule]->npoints;
      gaussPts.resize(3, nbGaussPts, false);
      cblas_dcopy(nbGaussPts, &QUAD_2D_TABLE[rule]->points[1], 3,
                  &gaussPts(0, 0), 1);
      cblas_dcopy(nbGaussPts, &QUAD_2D_TABLE[rule]->points[2], 3,
                  &gaussPts(1, 0), 1);
      cblas_dcopy(nbGaussPts, QUAD_2D_TABLE[rule]->weights, 1, &gaussPts(2, 0),
                  1);
      dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                             false);
      double *shape_ptr =
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
      cblas_dcopy(3 * nbGaussPts, QUAD_2D_TABLE[rule]->points, 1, shape_ptr, 1);
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_2D_TABLE_SIZE);
      nbGaussPts = 0;
    }
  } else {
    // If rule is negative, set user defined integration points
    CHKERR setGaussPts(order_row, order_col, order_data);
    nbGaussPts = gaussPts.size2();
    dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(nbGaussPts, 3,
                                                           false);
    if (nbGaussPts) {
      CHKERR ShapeMBTRI(
          &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
          &gaussPts(0, 0), &gaussPts(1, 0), nbGaussPts);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCoreBase::getSpaceBaseAndOrderOnElement() {
  MoFEMFunctionBegin;
  // Get spaces order/base and sense of entities.

  CHKERR getSpacesAndBaseOnEntities(dataH1);

  // H1
  if (dataH1.spacesOnEntities[MBEDGE].test(H1)) {
    CHKERR getEntitySense<MBEDGE>(dataH1);
    CHKERR getEntityDataOrder<MBEDGE>(dataH1, H1);
  }
  if (dataH1.spacesOnEntities[MBTRI].test(H1)) {
    CHKERR getEntitySense<MBTRI>(dataH1);
    CHKERR getEntityDataOrder<MBTRI>(dataH1, H1);
  }

  // Hcurl
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR getEntitySense<MBEDGE>(dataHcurl);
    CHKERR getEntityDataOrder<MBEDGE>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBEDGE].set(HCURL);
  }
  if (dataH1.spacesOnEntities[MBTRI].test(HCURL)) {
    CHKERR getEntitySense<MBTRI>(dataHcurl);
    CHKERR getEntityDataOrder<MBTRI>(dataHcurl, HCURL);
    dataHcurl.spacesOnEntities[MBTRI].set(HCURL);
  }

  // Hdiv
  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    CHKERR getEntitySense<MBTRI>(dataHdiv);
    CHKERR getEntityDataOrder<MBTRI>(dataHdiv, HDIV);
    dataHdiv.spacesOnEntities[MBTRI].set(HDIV);
  }

  // L2
  if (dataH1.spacesOnEntities[MBTRI].test(L2)) {
    CHKERR getEntitySense<MBTRI>(dataL2);
    CHKERR getEntityDataOrder<MBTRI>(dataL2, L2);
    dataL2.spacesOnEntities[MBTRI].set(L2);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCoreBase::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBeginHot;

  double *shape_functions =
      &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
  coordsAtGaussPts.resize(nbGaussPts, 3, false);
  for (int gg = 0; gg != nbGaussPts; ++gg)
    for (int dd = 0; dd != 3; ++dd)
      coordsAtGaussPts(gg, dd) =
          cblas_ddot(3, &shape_functions[3 * gg], 1, &coords[dd], 3);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCoreBase::calculateHoNormal() {
  MoFEMFunctionBegin;
  // Check if field for high-order geometry is set and if it is set calculate
  // higher-order normals and face tangent vectors.
  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {

    const Field *field_struture =
        mField.get_field_structure(meshPositionsFieldName);
    BitFieldId id = field_struture->getId();

    if ((numeredEntFiniteElementPtr->getBitFieldIdData() & id).none()) 
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
              "no MESH_NODE_POSITIONS in element data");

    // Calculate normal for high-order geometry

    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData(dataH1, meshPositionsFieldName, MBEDGE);
    CHKERR getEntityFieldData(dataH1, meshPositionsFieldName, MBEDGE);
    CHKERR opHOCoordsAndNormals.opRhs(dataH1);
    CHKERR opHOCoordsAndNormals.calculateNormals();

  } else {
    hoCoordsAtGaussPts.resize(0, 0, false);
    normalsAtGaussPts.resize(0, 0, false);
    tangentOneAtGaussPts.resize(0, 0, false);
    tangentTwoAtGaussPts.resize(0, 0, false);
  }
  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
