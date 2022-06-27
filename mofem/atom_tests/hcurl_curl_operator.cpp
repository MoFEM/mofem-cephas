/**
 * \file hcurl_curl_operator.cpp
 * \brief Testich curl-curl operator by applying Stokes theorem
 * \example hcurl_curl_operator.cpp
 *
 * Using PipelineManager interface calculate the curl of base functions, and
 * integral of the vector tangent vector with normal on the boundary. Since the
 * h-curl space is used, volume integral and boundary integral should give the
 * same result, as a result, as we are applying Stokes theorem on h-curl space.
 *
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

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

struct OpVolCurl : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  FTensor::Tensor1<double, 3> &tCurl;
  OpVolCurl(FTensor::Tensor1<double, 3> &t_curl)
      : VolumeElementForcesAndSourcesCore::UserDataOperator(
            "HCURL", UserDataOperator::OPROW),
        tCurl(t_curl) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

struct OpFacesRot : public FaceElementForcesAndSourcesCore::UserDataOperator {

  FTensor::Tensor1<double, 3> &tCurl;
  OpFacesRot(FTensor::Tensor1<double, 3> &t_curl)
      : FaceElementForcesAndSourcesCore::UserDataOperator(
            "HCURL", UserDataOperator::OPROW),
        tCurl(t_curl) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    enum bases { AINSWORTH, DEMKOWICZ, LASTOP };

    const char *list[] = {"ainsworth", "demkowicz"};

    PetscBool flg;
    PetscInt choise_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list, LASTOP,
                                &choise_value, &flg);
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    PetscBool ho_geometry = PETSC_FALSE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-ho_geometry", &ho_geometry,
                               PETSC_NULL);

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Create MoAB
    moab::Core mb_instance;              ///< database
    moab::Interface &moab = mb_instance; ///< interface

    // Create MoFEM
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    Simple *simple_interface = m_field.getInterface<Simple>();
    PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();
    CHKERR simple_interface->getOptions();
    CHKERR simple_interface->loadFile("");

    // fields
    switch (choise_value) {
    case AINSWORTH:
      CHKERR simple_interface->addDomainField("HCURL", HCURL,
                                              AINSWORTH_LEGENDRE_BASE, 1);
      CHKERR simple_interface->addBoundaryField("HCURL", HCURL,
                                                AINSWORTH_LEGENDRE_BASE, 1);
      break;
    case DEMKOWICZ:
      CHKERR simple_interface->addDomainField("HCURL", HCURL,
                                              DEMKOWICZ_JACOBI_BASE, 1);
      CHKERR simple_interface->addBoundaryField("HCURL", HCURL,
                                                DEMKOWICZ_JACOBI_BASE, 1);
      break;
    }

    if (ho_geometry == PETSC_TRUE)
      CHKERR simple_interface->addDataField("MESH_NODE_POSITIONS", H1,
                                            AINSWORTH_LEGENDRE_BASE, 3);

    constexpr int order = 3;
    CHKERR simple_interface->setFieldOrder("HCURL", order);
    // Range ents;
    // CHKERR moab.get_entities_by_dimension(0, 2, ents, true);
    // CHKERR simple_interface->setFieldOrder("HCURL", 1, &ents);

    if (ho_geometry == PETSC_TRUE)
      CHKERR simple_interface->setFieldOrder("MESH_NODE_POSITIONS", 2);
    CHKERR simple_interface->setUp();

    auto integration_rule = [](int, int, int p_data) { return 2 * p_data; };
    CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
    CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(integration_rule);

    FTensor::Tensor1<double, 3> t_curl_vol;
    FTensor::Tensor1<double, 3> t_curl_skin;

    auto material_grad_mat = boost::make_shared<MatrixDouble>();
    auto material_det_vec = boost::make_shared<VectorDouble>();
    auto material_inv_grad_mat = boost::make_shared<MatrixDouble>();
    auto jac_ptr = boost::make_shared<MatrixDouble>();
    auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
    auto det_ptr = boost::make_shared<VectorDouble>();

    boost::dynamic_pointer_cast<VolumeElementForcesAndSourcesCore>(
        pipeline_mng->getDomainRhsFE())
        ->meshPositionsFieldName = "none";
    boost::dynamic_pointer_cast<PipelineManager::FaceEle>(
        pipeline_mng->getBoundaryRhsFE())
        ->meshPositionsFieldName = "none";

    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpCalculateHOJac<3>(jac_ptr));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpInvertMatrix<3>(jac_ptr, det_ptr, inv_jac_ptr));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpSetHOWeights(det_ptr));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpSetHOCovariantPiolaTransform(HCURL, inv_jac_ptr));
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpSetHOInvJacVectorBase(HCURL, inv_jac_ptr));

    if (ho_geometry) {
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateVectorFieldGradient<3, 3>("MESH_NODE_POSITIONS",
                                                   material_grad_mat));
      pipeline_mng->getOpDomainRhsPipeline().push_back(new OpInvertMatrix<3>(
          material_grad_mat, material_det_vec, material_inv_grad_mat));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetHOWeights(material_det_vec));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetHOCovariantPiolaTransform(HCURL, material_inv_grad_mat));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetHOInvJacVectorBase(HCURL, material_inv_grad_mat));
    } 
    pipeline_mng->getOpDomainRhsPipeline().push_back(new OpVolCurl(t_curl_vol));

    if (m_field.check_field("MESH_NODE_POSITIONS"))
      pipeline_mng->getOpBoundaryRhsPipeline().push_back(
          new OpGetHONormalsOnFace("MESH_NODE_POSITIONS"));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpHOSetCovariantPiolaTransformOnFace3D(HCURL));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpFacesRot(t_curl_skin));

    FTensor::Index<'i', 3> i;
    t_curl_vol(i) = 0;
    t_curl_skin(i) = 0;

    // project geometry form 10 node tets on higher order approx. functions
    if (ho_geometry == PETSC_TRUE) {
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
    }

    // Run pipelines on mesh
    CHKERR pipeline_mng->loopFiniteElements();

    std::cout.precision(12);

    t_curl_vol(i) -= t_curl_skin(i);
    double nrm2 = sqrt(t_curl_vol(i) * t_curl_vol(i));

    constexpr double eps = 1e-8;
    if (fabs(nrm2) > eps)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Curl operator not passed test\n");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}

MoFEMErrorCode OpVolCurl::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (data.getFieldData().size() == 0)
    MoFEMFunctionReturnHot(0);

  const unsigned int nb_gauss_pts = data.getDiffN().size1();
  const unsigned int nb_dofs = data.getFieldData().size();

  MatrixDouble curl_mat;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;

  auto t_curl_base = data.getFTensor2DiffN<3, 3>();

  unsigned int gg = 0;
  for (; gg < nb_gauss_pts; gg++) {
    double w = getGaussPts()(3, gg) * getVolume();
    FTensor::Tensor1<double, 3> t_curl;
    for (unsigned int dd = 0; dd != nb_dofs; dd++) {
      t_curl(i) = levi_civita(j, i, k) * t_curl_base(j, k);
      tCurl(i) += w * t_curl(i);
      ++t_curl_base;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpFacesRot::doWork(int side, EntityType type,
                                  EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  int nb_dofs = data.getFieldData().size();
  if (nb_dofs == 0)
    MoFEMFunctionReturnHot(0);
  int nb_gauss_pts = data.getN().size1();

  auto t_curl_base = data.getFTensor1N<3>();

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  for (int gg = 0; gg < nb_gauss_pts; gg++) {
    for (int dd = 0; dd < nb_dofs; dd++) {
      double w = getGaussPts()(2, gg);
      const double n0 = getNormalsAtGaussPts(gg)[0];
      const double n1 = getNormalsAtGaussPts(gg)[1];
      const double n2 = getNormalsAtGaussPts(gg)[2];
      if (getFEType() == MBTRI) {
        w *= 0.5;
      }

      tCurl(0) += (n1 * t_curl_base(2) - n2 * t_curl_base(1)) * w;
      tCurl(1) += (n2 * t_curl_base(0) - n0 * t_curl_base(2)) * w;
      tCurl(2) += (n0 * t_curl_base(1) - n1 * t_curl_base(0)) * w;
      ++t_curl_base;
    }
  }

  MoFEMFunctionReturn(0);
}