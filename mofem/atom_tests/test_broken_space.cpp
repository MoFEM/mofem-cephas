/**
 * \example test_broken_space.cpp
 *
 * Testing broken spaces
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using EleOnSide = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::FaceSideEle;

using EntData = EntitiesFieldData::EntData;
using DomainEleOp = DomainEle::UserDataOperator;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "AT"));
    LogManager::setLog("AT");
    MOFEM_LOG_TAG("AT", "atom_test");

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    auto *simple = m_field.getInterface<Simple>();
    CHKERR simple->getOptions();

    simple->getAddBoundaryFE() = true;

    CHKERR simple->loadFile("", "");

    // Declare elements
    enum bases {
      AINSWORTH,
      AINSWORTH_LOBATTO,
      DEMKOWICZ,
      BERNSTEIN,
      LASBASETOP
    };
    const char *list_bases[] = {"ainsworth", "ainsworth_labatto", "demkowicz",
                                "bernstein"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);

    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH_LOBATTO)
      base = AINSWORTH_LOBATTO_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    else if (choice_base_value == BERNSTEIN)
      base = AINSWORTH_BERNSTEIN_BEZIER_BASE;

    enum spaces { hdiv, hcurl, last_space };
    const char *list_spaces[] = {"hdiv", "hcurl"};
    PetscInt choice_space_value = hdiv;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                last_space, &choice_space_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSSIBLE_CASE, "space not set");
    FieldSpace space = HDIV;
    if (choice_space_value == hdiv)
      space = HDIV;
    else if (choice_space_value == hcurl)
      space = HCURL;

    int approx_order = 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &approx_order,
                              PETSC_NULL);

    CHKERR simple->addDomainBrokenField("BROKEN", space, base, 1);
    CHKERR simple->addSkeletonField("HYBRID", L2, base, SPACE_DIM);
    CHKERR simple->addBoundaryField("HYBRID", L2, base, SPACE_DIM);
    CHKERR simple->setFieldOrder("BROKEN", approx_order);
    CHKERR simple->setFieldOrder("HYBRID", approx_order);

    CHKERR simple->setUp();

    auto assemble_domain = [&]() {
      MoFEMFunctionBegin;
      auto *pip_mng = m_field.getInterface<PipelineManager>();

      using OpSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<3, 3>;
      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<
          SCHUR>::BiLinearForm<GAUSS>::OpMass<3, 3>;

      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          pip_mng->getOpDomainLhsPipeline(), {space});
      CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
          pip_mng->getOpDomainRhsPipeline(), {space});

      pip_mng->getOpDomainLhsPipeline().push_back(new OpMass(
          "BROKEN", "BROKEN", [](double, double, double) { return 1.; }));

      pip_mng->getOpDomainRhsPipeline().push_back(
          new OpSource("BROKEN", [](double, double, double) {
            return FTensor::Tensor1<double, 3>{1., 0., 0.};
          }));

      auto integration_rule = [](int, int, int p_data) {
        return 2 * p_data + 1;
      };
      CHKERR pip_mng->setDomainRhsIntegrationRule(integration_rule);
      CHKERR pip_mng->setDomainLhsIntegrationRule(integration_rule);

      MoFEMFunctionReturn(0);
    };

    auto assemble_skeleton = [&](auto &pip_lhs) {
      MoFEMFunctionBegin;
			MoFEMFunctionReturn(0);
    };
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
