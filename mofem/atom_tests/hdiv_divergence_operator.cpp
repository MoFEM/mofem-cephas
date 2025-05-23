/**
 * \file hdiv_divergence_operator.cpp
 * \example hdiv_divergence_operator.cpp
 *
 * Using PipelineManager interface calculate the divergence of base functions,
 * and integral of flux on the boundary. Since the h-div space is used, volume
 * integral and boundary integral should give the same result.
 */



#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

struct OpVolDivergence
    : public VolumeElementForcesAndSourcesCore::UserDataOperator {

  double &dIv;
  OpVolDivergence(double &div)
      : VolumeElementForcesAndSourcesCore::UserDataOperator(
            "HDIV", UserDataOperator::OPROW),
        dIv(div) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

struct OpFacesFluxes
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  double &dIv;
  OpFacesFluxes(double &div)
      : FaceElementForcesAndSourcesCore::UserDataOperator(
            "HDIV", UserDataOperator::OPROW),
        dIv(div) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
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
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "base not set");
    }

    PetscBool ho_geometry = PETSC_FALSE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-ho_geometry", &ho_geometry,
                               PETSC_NULL);

    PetscInt ho_choice_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-ho_base", list, LASTOP,
                                &ho_choice_value, &flg);

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

    if (ho_geometry == PETSC_TRUE) {
      switch (ho_choice_value) {
      case AINSWORTH:
        CHKERR simple_interface->addDataField("MESH_NODE_POSITIONS", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      case DEMKOWICZ:
        CHKERR simple_interface->addDataField("MESH_NODE_POSITIONS", H1,
                                              DEMKOWICZ_JACOBI_BASE, 3);
      }
    }

    constexpr int order = 3;
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

    if (m_field.check_field("MESH_NODE_POSITIONS")) {
      CHKERR AddHOOps<3, 3, 3>::add(pipeline_mng->getOpDomainRhsPipeline(),
                                    {HDIV}, "MESH_NODE_POSITIONS");
      CHKERR AddHOOps<2, 3, 3>::add(pipeline_mng->getOpBoundaryRhsPipeline(),
                                    {HDIV}, "MESH_NODE_POSITIONS");
    } else {
      CHKERR AddHOOps<3, 3, 3>::add(pipeline_mng->getOpDomainRhsPipeline(),
                                    {HDIV});
      CHKERR AddHOOps<2, 3, 3>::add(pipeline_mng->getOpBoundaryRhsPipeline(),
                                    {HDIV});
    }

    // domain
    pipeline_mng->getOpDomainRhsPipeline().push_back(
        new OpVolDivergence(divergence_vol));
    // boundary
    pipeline_mng->getOpBoundaryRhsPipeline().push_back(
        new OpFacesFluxes(divergence_skin));

    // project geometry form 10 node tets on higher order approx. functions
    if (ho_geometry == PETSC_TRUE) {
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);
    }
    CHKERR pipeline_mng->loopFiniteElements();

    MOFEM_LOG("WORLD", Sev::inform)
        << "divergence_vol " << std::setprecision(12) << divergence_vol;
    MOFEM_LOG("WORLD", Sev::inform)
        << "divergence_skin " << std::setprecision(12) << divergence_skin;

    constexpr double eps = 1e-8;
    const double error = divergence_skin - divergence_vol;
    if (fabs(error) > eps)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
               "invalid surface flux or divergence or both error = %3.4e",
               error);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}

MoFEMErrorCode
OpVolDivergence::doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (CN::Dimension(type) < 2)
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
                                     EntitiesFieldData::EntData &data) {
  MoFEMFunctionBeginHot;

  if (CN::Dimension(type) != 2)
    MoFEMFunctionReturnHot(0);

  int nb_gauss_pts = data.getN().size1();
  int nb_dofs = data.getFieldData().size();

  for (int gg = 0; gg < nb_gauss_pts; gg++) {
    for (int dd = 0; dd < nb_dofs; dd++) {

      double w = getGaussPts()(2, gg);
      const double n0 = getNormalsAtGaussPts(gg)[0];
      const double n1 = getNormalsAtGaussPts(gg)[1];
      const double n2 = getNormalsAtGaussPts(gg)[2];
      if (getFEType() == MBTRI) {
        w *= 0.5;
      }

      dIv += (n0 * data.getVectorN<3>(gg)(dd, 0) +
              n1 * data.getVectorN<3>(gg)(dd, 1) +
              n2 * data.getVectorN<3>(gg)(dd, 2)) *
             w;
    }
  }

  MoFEMFunctionReturnHot(0);
}
