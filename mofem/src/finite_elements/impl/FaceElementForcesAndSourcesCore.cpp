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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <base_functions.h>
#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ForcesAndSourcesCore.hpp>
#include <VolumeElementForcesAndSourcesCore.hpp>
#include <FaceElementForcesAndSourcesCore.hpp>

#include <BaseFunction.hpp>
#include <EntPolynomialBaseCtx.hpp>
#include <TriPolynomialBase.hpp>

#include <BitRefManager.hpp>

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

FaceElementForcesAndSourcesCore::FaceElementForcesAndSourcesCore(
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
FaceElementForcesAndSourcesCore::UserDataOperator::loopSideVolumes(
    const string &fe_name, VolumeElementForcesAndSourcesCoreOnSide &method) {
  MoFEMFunctionBegin;

  const EntityHandle ent = getNumeredEntFiniteElementPtr()->getEnt();
  const Problem *problem_ptr = getFEMethod()->problemPtr;
  Range adjacent_volumes;
  CHKERR getFaceFE()->mField.getInterface<BitRefManager>()->getAdjacenciesAny(
      ent, 3, adjacent_volumes);
  typedef NumeredEntFiniteElement_multiIndex::index<
      Composite_Name_And_Ent_mi_tag>::type FEByComposite;
  FEByComposite &numered_fe = (const_cast<NumeredEntFiniteElement_multiIndex &>(
                                   problem_ptr->numeredFiniteElements))
                                  .get<Composite_Name_And_Ent_mi_tag>();

  method.feName = fe_name;

  CHKERR method.setFaceFEPtr(getFaceFE());
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
      // cerr << **miit << endl;
      // cerr << &(**miit) << endl;
      // cerr << (*miit)->getEnt() << endl;
      method.nInTheLoop = nn++;
      method.numeredEntFiniteElementPtr = *miit;
      method.dataPtr = (*miit)->sPtr->data_dofs;
      method.rowPtr = (*miit)->rows_dofs;
      method.colPtr = (*miit)->cols_dofs;
      CHKERR method();
    }
  }

  CHKERR method.postProcess();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FaceElementForcesAndSourcesCore::calculateAreaAndNormal() {
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

MoFEMErrorCode FaceElementForcesAndSourcesCore::setIntegrationPts() {
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
FaceElementForcesAndSourcesCore::getSpaceBaseAndOrderOnElement() {
  MoFEMFunctionBegin;
  // Get spaces order/base and sense of entities.
  
  
  DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];
  DataForcesAndSourcesCore &data_div = *dataOnElement[HDIV];
  DataForcesAndSourcesCore &data_l2 = *dataOnElement[L2]; 

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
    CHKERR getEntitySense<MBEDGE>(data_curl);
    CHKERR getEntityDataOrder<MBEDGE>(data_curl, HCURL);
    data_curl.spacesOnEntities[MBEDGE].set(HCURL);
  }
  if (dataH1.spacesOnEntities[MBTRI].test(HCURL)) {
    CHKERR getEntitySense<MBTRI>(data_curl);
    CHKERR getEntityDataOrder<MBTRI>(data_curl, HCURL);
    data_curl.spacesOnEntities[MBTRI].set(HCURL);
  }

  // Hdiv
  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    CHKERR getEntitySense<MBTRI>(data_div);
    CHKERR getEntityDataOrder<MBTRI>(data_div, HDIV);
    data_div.spacesOnEntities[MBTRI].set(HDIV);
  }

  // L2
  if (dataH1.spacesOnEntities[MBTRI].test(L2)) {
    CHKERR getEntitySense<MBTRI>(data_l2);
    CHKERR getEntityDataOrder<MBTRI>(data_l2, L2);
    data_l2.spacesOnEntities[MBTRI].set(L2);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
FaceElementForcesAndSourcesCore::calculateCoordinatesAtGaussPts() {
  MoFEMFunctionBeginHot;
  

  double *shape_functions =
      &*dataH1.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin();
  coordsAtGaussPts.resize(nbGaussPts, 3, false);
  for (int gg = 0; gg != nbGaussPts; ++gg) {
    for (int dd = 0; dd != 3; ++dd) {
      coordsAtGaussPts(gg, dd) =
          cblas_ddot(3, &shape_functions[3 * gg], 1, &coords[dd], 3);
    }
  }
  MoFEMFunctionReturnHot(0);
}


MoFEMErrorCode FaceElementForcesAndSourcesCore::calculateHoNormal() {
  MoFEMFunctionBegin;
  // Check if field for high-order geometry is set and if it is set calculate
  // higher-order normals and face tangent vectors.
  if (dataPtr->get<FieldName_mi_tag>().find(meshPositionsFieldName) !=
      dataPtr->get<FieldName_mi_tag>().end()) {

    const Field *field_struture =
        mField.get_field_structure(meshPositionsFieldName);
    BitFieldId id = field_struture->getId();

    if ((numeredEntFiniteElementPtr->getBitFieldIdData() & id).none()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
              "no MESH_NODE_POSITIONS in element data");
    }

    // Calculate normal for high-order geometry
    

    CHKERR getEntityDataOrderSpaceAndBase<MBEDGE>(dataH1, meshPositionsFieldName);
    CHKERR getEntityDataOrderSpaceAndBase<MBTRI>(dataH1, meshPositionsFieldName);
    CHKERR getNodesFieldData(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData<MBEDGE>(dataH1, meshPositionsFieldName);
    CHKERR getEntityFieldData<MBTRI>(dataH1, meshPositionsFieldName);
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

MoFEMErrorCode FaceElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBTRI)
    MoFEMFunctionReturnHot(0);
  CHKERR createDataOnElement();

  // Calculate normal and tangent vectors for face geometry given by 3 nodes.
  CHKERR calculateAreaAndNormal();
  CHKERR getSpaceBaseAndOrderOnElement();

  CHKERR setIntegrationPts();
  if (nbGaussPts == 0)
    MoFEMFunctionReturnHot(0);

  
  DataForcesAndSourcesCore &data_curl = *dataOnElement[HCURL];
  DataForcesAndSourcesCore &data_div = *dataOnElement[HDIV];

  dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(3, 2, false);
  CHKERR ShapeDiffMBTRI(
      &*dataH1.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin());

  /// Use the some node base
  CHKERR calculateCoordinatesAtGaussPts();
  CHKERR calculateBaseFunctionsOnElement();
  CHKERR calculateHoNormal();

  // Apply Piola transform to HDiv and HCurl spaces, uses previously calculated
  // faces normal and tangent vectors.
  if (dataH1.spacesOnEntities[MBTRI].test(HDIV)) {
    CHKERR opContravariantTransform.opRhs(data_div);
  }
  if (dataH1.spacesOnEntities[MBEDGE].test(HCURL)) {
    CHKERR opCovariantTransform.opRhs(data_curl);
  }

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFace::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBeginHot;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");
  }

  try {

    if (type == MBVERTEX) {
      VectorDouble &coords = getCoords();
      double *coords_ptr = &*coords.data().begin();
      double diff_n[6];
      ierr = ShapeDiffMBTRI(diff_n);
      CHKERRG(ierr);
      double j00, j01, j10, j11;
      for (int gg = 0; gg < 1; gg++) {
        // this is triangle, derivative of nodal shape functions is constant.
        // So only need to do one node.
        j00 = cblas_ddot(3, &coords_ptr[0], 3, &diff_n[0], 2);
        j01 = cblas_ddot(3, &coords_ptr[0], 3, &diff_n[1], 2);
        j10 = cblas_ddot(3, &coords_ptr[1], 3, &diff_n[0], 2);
        j11 = cblas_ddot(3, &coords_ptr[1], 3, &diff_n[1], 2);
      }
      double det = j00 * j11 - j01 * j10;
      invJac.resize(2, 2, false);
      invJac(0, 0) = j11 / det;
      invJac(0, 1) = -j01 / det;
      invJac(1, 0) = -j10 / det;
      invJac(1, 1) = j00 / det;
    }
  } catch (std::exception &ex) {
    std::ostringstream ss;
    ss << "throw in method: " << ex.what() << " at line " << __LINE__
       << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
  }

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets = false;
  doPrisms = false;

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
OpSetInvJacH1ForFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI &&
      getNumeredEntFiniteElementPtr()->getEntType() != MBQUAD) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");
  }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    unsigned int nb_dofs = data.getN(base).size2();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    unsigned int nb_gauss_pts = data.getN(base).size1();
    diffNinvJac.resize(nb_gauss_pts, 2 * nb_dofs, false);

    if (type != MBVERTEX) {
      if (nb_dofs != data.getDiffN(base).size2() / 2) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "data inconsistency nb_dofs != data.diffN.size2()/2 ( %u != "
                 "%u/2 )",
                 nb_dofs, data.getDiffN(base).size2());
      }
    }

    FTensor::Tensor2<double, 2, 2> t_inv_jac;
    t_inv_jac(0, 0) = invJac(0, 0);
    t_inv_jac(0, 1) = invJac(0, 1);
    t_inv_jac(1, 0) = invJac(1, 0);
    t_inv_jac(1, 1) = invJac(1, 1);

    switch (type) {
    case MBVERTEX:
    case MBEDGE:
    case MBTRI: {
      FTensor::Index<'i', 2> i;
      FTensor::Index<'j', 2> j;
      FTensor::Index<'k', 2> k;
      FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n(
          &diffNinvJac(0, 0), &diffNinvJac(0, 1));
      FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n_ref(
          &data.getDiffN(base)(0, 0), &data.getDiffN(base)(0, 1));
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (unsigned int dd = 0; dd != nb_dofs; ++dd) {
          t_diff_n(i) = t_inv_jac(i, k) * t_diff_n_ref(k);
          ++t_diff_n;
          ++t_diff_n_ref;
        }
      }
      data.getDiffN(base).data().swap(diffNinvJac.data());
    } break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHcurlFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI)
    MoFEMFunctionReturnHot(0);

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");
  }

  FTensor::Tensor2<double *, 2, 2> t_inv_jac = FTensor::Tensor2<double *, 2, 2>(
      &invJac(0, 0), &invJac(0, 1), &invJac(1, 0), &invJac(1, 1));

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = ApproximationBaseArray[b];

    const unsigned int nb_base_functions = data.getVectorDiffN(base).size2() / 6;
    if (!nb_base_functions)
      continue;
    const unsigned int nb_gauss_pts = data.getVectorDiffN(base).size1();

    diffHcurlInvJac.resize(nb_gauss_pts, data.getVectorDiffN(base).size2(),
                           false);

    auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
    double *t_inv_diff_n_ptr = &*diffHcurlInvJac.data().begin();
    FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_inv_diff_n(
        t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],
        &t_inv_diff_n_ptr[HVEC1_0], &t_inv_diff_n_ptr[HVEC1_1],
        &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1]);

    for (unsigned int gg = 0; gg != nb_gauss_pts; gg++) {
      for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
        t_inv_diff_n(i, j) = t_diff_n(i, k) * t_inv_jac(k, j);
        ++t_diff_n;
        ++t_inv_diff_n;
      }
    }

    data.getVectorDiffN(base).data().swap(diffHcurlInvJac.data());
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
