/**
 * \file hdiv_divergence_operator.cpp
 * \example hdiv_divergence_operator.cpp
 *
 * Using PipelineManager interface calculate the divergence of base functions, and
 * integral of flux on the boundary. Since the h-div space is used, volume
 * integral and boundary integral should give the same result.
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

struct OpTetDivergence
    : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  double &dIv;
  OpTetDivergence(double &div)
      : VolumeElementForcesAndSourcesCore::UserDataOperator(
            "HDIV", UserDataOperator::OPROW),
        dIv(div) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

struct OpFacesFluxes
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  double &dIv;
  OpFacesFluxes(double &div)
      : FaceElementForcesAndSourcesCore::UserDataOperator(
            "HDIV", UserDataOperator::OPROW),
        dIv(div) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    enum bases { AINSWORTH, DEMKOWICZ, LASTOP };

    const char *list[] = {"ainsworth", "demkowicz"};

    PetscBool flg;
    PetscInt choice_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list, LASTOP,
                                &choice_value, &flg);
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
    switch (choice_value) {
    case AINSWORTH:
      CHKERR simple_interface->addDomainField("HDIV", HDIV,
                                             AINSWORTH_LEGENDRE_BASE, 1);
      CHKERR simple_interface->addBoundaryField("HDIV", HDIV,
                                               AINSWORTH_LEGENDRE_BASE, 1);
      break;
    case DEMKOWICZ:
      CHKERR simple_interface->addDomainField("HDIV", HDIV,
                                             DEMKOWICZ_JACOBI_BASE, 1);
      CHKERR simple_interface->addBoundaryField("HDIV", HDIV,
                                               DEMKOWICZ_JACOBI_BASE, 1);
      break;
    }

    if (ho_geometry == PETSC_TRUE)
      CHKERR simple_interface->addDataField("MESH_NODE_POSITIONS", H1,
                                           AINSWORTH_LEGENDRE_BASE, 3);

    constexpr int order = 5;
    CHKERR simple_interface->setFieldOrder("HDIV", order);
    if (ho_geometry == PETSC_TRUE)
      CHKERR simple_interface->setFieldOrder("MESH_NODE_POSITIONS", 2);
    CHKERR simple_interface->setUp();

    /// This has no real effect, folling line are only for atom test purpose
    pipeline_mng->getDomainLhsFE().reset();
    pipeline_mng->getDomainRhsFE().reset();
    pipeline_mng->getBoundaryLhsFE().reset();
    pipeline_mng->getBoundaryRhsFE().reset();
    pipeline_mng->getSkeletonLhsFE().reset();
    pipeline_mng->getSkeletonRhsFE().reset();

    auto integration_rule = [](int, int, int p_data) { return 2 * p_data; };
    CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
    CHKERR pipeline_mng->setBoundaryRhsIntegrationRule(integration_rule);

    double divergence_vol = 0;
    double divergence_skin = 0;
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpTetDivergence(divergence_vol));

    boost::dynamic_pointer_cast<PipelineManager::FaceEle>(
        pipeline_mng->getBoundaryRhsFE())
        ->meshPositionsFieldName = "none";

    if (m_field.check_field("MESH_NODE_POSITIONS"))
      pipeline_mng->getOpBoundaryRhsPipeline().push_back(
          new OpGetHONormalsOnFace("MESH_NODE_POSITIONS"));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpHOSetContravariantPiolaTransformOnFace3D(HDIV));
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpFacesFluxes(divergence_skin));

    // project geometry form 10 node tets on higher order approx. functions
    if (ho_geometry == PETSC_TRUE) {
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
    }
    CHKERR pipeline_mng->loopFiniteElements();

    std::cout.precision(12);

    std::cout << "divergence_vol " << divergence_vol << std::endl;
    std::cout << "divergence_skin " << divergence_skin << std::endl;

    constexpr double eps = 1e-8;
    if (fabs(divergence_skin - divergence_vol) > eps)
      SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "invalid surface flux or divergence or both\n", divergence_skin,
               divergence_vol);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}

MoFEMErrorCode
OpTetDivergence::doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBTRI && type != MBTET)
    MoFEMFunctionReturnHot(0);

  if (data.getFieldData().size() == 0)
    MoFEMFunctionReturnHot(0);

  int nb_gauss_pts = data.getDiffN().size1();
  int nb_dofs = data.getFieldData().size();

  VectorDouble div_vec;
  div_vec.resize(nb_dofs, 0);

  FTensor::Index<'i', 3> i;
  auto t_base_diff_hdiv = data.getFTensor2DiffN<3, 3>();

  for (size_t gg = 0; gg < nb_gauss_pts; gg++) {
    for (size_t dd = 0; dd != div_vec.size(); dd++) {
      double w = getGaussPts()(3, gg) * getVolume();
      dIv += t_base_diff_hdiv(i, i) * w;
      ++t_base_diff_hdiv;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpFacesFluxes::doWork(int side, EntityType type,
                                     DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;

  if (type != MBTRI)
    MoFEMFunctionReturnHot(0);

  int nb_gauss_pts = data.getN().size1();
  int nb_dofs = data.getFieldData().size();

  int gg = 0;
  for (; gg < nb_gauss_pts; gg++) {
    int dd = 0;
    for (; dd < nb_dofs; dd++) {
      double area;
      VectorDouble n;
      if (getNormalsAtGaussPts().size1() == (unsigned int)nb_gauss_pts) {
        n = getNormalsAtGaussPts(gg);
        area = norm_2(getNormalsAtGaussPts(gg)) * 0.5;
      } else {
        n = getNormal();
        area = getArea();
      }
      n /= norm_2(n);
      dIv += (n[0] * data.getVectorN<3>(gg)(dd, 0) +
              n[1] * data.getVectorN<3>(gg)(dd, 1) +
              n[2] * data.getVectorN<3>(gg)(dd, 2)) *
             getGaussPts()(2, gg) * area;
    }
  }

  MoFEMFunctionReturnHot(0);
}
