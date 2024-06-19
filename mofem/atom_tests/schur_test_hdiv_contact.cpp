/**
 * @file schur_test_hdiv_contact.cpp
 * @example schur_test_hdiv_contact.cpp
 * @brief Test block matrix and Schur complement matrix for contact problem with
 * H(div) approximation of LM
 *
 * @copyright Copyright (c) 2024
 *
 */

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
};

//! [Define dimension]
constexpr int SPACE_DIM = 2;
constexpr int FIELD_DIM = SPACE_DIM;

constexpr AssemblyType A = AssemblyType::BLOCK_MAT; //< selected assembly type
constexpr IntegrationType I =
    IntegrationType::GAUSS; //< selected integration type

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;

SmartPetscObj<Mat> petsc_mat;
SmartPetscObj<Mat> block_mat;

template <>
typename MoFEM::OpBaseImpl<PETSC, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<PETSC, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<PETSC>>(
              petsc_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_MAT, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_MAT, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_SCHUR, DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_SCHUR, DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_SCHUR>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

template <>
typename MoFEM::OpBaseImpl<BLOCK_PRECONDITIONER_SCHUR,
                           DomainEleOp>::MatSetValuesHook
    MoFEM::OpBaseImpl<BLOCK_PRECONDITIONER_SCHUR,
                      DomainEleOp>::matSetValuesHook =
        [](ForcesAndSourcesCore::UserDataOperator *op_ptr,
           const EntitiesFieldData::EntData &row_data,
           const EntitiesFieldData::EntData &col_data, MatrixDouble &m) {
          return MatSetValues<AssemblyTypeSelector<BLOCK_PRECONDITIONER_SCHUR>>(
              block_mat, row_data, col_data, m, ADD_VALUES);
        };

constexpr bool debug = false;

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface &moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface &m_field = mofem_core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Add logging channel for example
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "Timeline"));
    LogManager::setLog("Timeline");
    MOFEM_LOG_TAG("Timeline", "Timeline");

    // Simple interface
    auto simple = m_field.getInterface<Simple>();

    // get options from command line
    CHKERR simple->getOptions();
    // load mesh file
    CHKERR simple->loadFile();

    CHKERR simple->addDomainField("U", H1, AINSWORTH_LEGENDRE_BASE, SPACE_DIM);
    CHKERR simple->addBoundaryField("U", H1, AINSWORTH_LEGENDRE_BASE,
                                    SPACE_DIM);
    CHKERR simple->addDomainField("SIGMA", HCURL, DEMKOWICZ_JACOBI_BASE,
                                  SPACE_DIM);
    CHKERR simple->addBoundaryField("SIGMA", HCURL, DEMKOWICZ_JACOBI_BASE,
                                    SPACE_DIM);

    int order = 2;       //< Order of displacements in the domain
    int sigma_order = 1; //< Order of Lagrange multiplier in side elements
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);
    CHKERR simple->setFieldOrder("U", order);

    sigma_order = order - 1;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-sigma_order", &sigma_order,
                              PETSC_NULL);

    MOFEM_LOG("Timeline", Sev::inform) << "Order " << order;
    MOFEM_LOG("Timeline", Sev::inform) << "Sigma order " << sigma_order;

    auto get_skin = [&]() {
      Range body_ents;
      CHKERR m_field.get_moab().get_entities_by_dimension(0, SPACE_DIM,
                                                          body_ents);
      Skinner skin(&m_field.get_moab());
      Range skin_ents;
      CHKERR skin.find_skin(0, body_ents, false, skin_ents);
      return skin_ents;
    };

    auto filter_blocks = [&](auto skin) {
      bool is_contact_block = false;
      Range contact_range;
      for (auto m : m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
               std::regex(

                   (boost::format("%s(.*)") % "CONTACT").str()

                       ))

      ) {
        is_contact_block = true; ///< blocs interation is collective, so that is
                                 ///< set irrespective if there are entities in
                                 ///< given rank or not in the block
        MOFEM_LOG("Timeline", Sev::inform)
            << "Find contact block set:  " << m->getName();
        auto meshset = m->getMeshset();
        Range contact_meshset_range;
        CHKERR m_field.get_moab().get_entities_by_dimension(
            meshset, SPACE_DIM - 1, contact_meshset_range, true);

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            contact_meshset_range);
        contact_range.merge(contact_meshset_range);
      }
      if (is_contact_block) {
        MOFEM_LOG("Timeline", Sev::inform)
            << "Nb entities in contact surface: " << contact_range.size();
        MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
        skin = intersect(skin, contact_range);
      }
      return skin;
    };

    auto filter_true_skin = [&](auto skin) {
      Range boundary_ents;
      ParallelComm *pcomm =
          ParallelComm::get_pcomm(&m_field.get_moab(), MYPCOMM_INDEX);
      CHKERR pcomm->filter_pstatus(skin, PSTATUS_SHARED | PSTATUS_MULTISHARED,
                                   PSTATUS_NOT, -1, &boundary_ents);
      return boundary_ents;
    };

    auto boundary_ents = filter_true_skin(filter_blocks(get_skin()));
    CHKERR simple->setFieldOrder("SIGMA", 0);
    CHKERR simple->setFieldOrder("SIGMA", sigma_order, &boundary_ents);

    // setup problem
    CHKERR simple->setUp();

    auto schur_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(schur_dm, simple->getDM(), "SCHUR");
    CHKERR DMMoFEMSetSquareProblem(schur_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(schur_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(schur_dm, "U");
    CHKERR DMMoFEMAddSubFieldCol(schur_dm, "U");
    CHKERR DMSetUp(schur_dm);

    auto block_dm = createDM(m_field.get_comm(), "DMMOFEM");
    CHKERR DMMoFEMCreateSubDM(block_dm, simple->getDM(), "BLOCK");
    CHKERR DMMoFEMSetSquareProblem(block_dm, PETSC_TRUE);
    CHKERR DMMoFEMAddElement(block_dm, simple->getDomainFEName());
    CHKERR DMMoFEMAddSubFieldRow(block_dm, "SIGMA");
    CHKERR DMMoFEMAddSubFieldCol(block_dm, "SIGMA");

    CHKERR DMSetUp(block_dm);

    petsc_mat = createDMMatrix(simple->getDM());
    auto S = createDMMatrix(schur_dm);

    auto shell_data = createBlockMat(
        simple->getDM(),

        createBlockMatStructure(
            simple->getDM(),

            {{simple->getDomainFEName(),

              {{"U", "U"}, {"SIGMA", "SIGMA"}, {"U", "SIGMA"}, {"SIGMA", "U"}}}}

            )

    );

    auto [mat, block_data_ptr] = shell_data;
    block_mat = mat;

    std::vector<std::string> fields{"SIGMA"};

    auto [nested_mat, nested_data_ptr] = createSchurNestedMatrix(

        getNestSchurData(

            {schur_dm, block_dm}, block_data_ptr,

            fields, {nullptr, nullptr, nullptr}, true

            )

    );

    MOFEM_LOG("Timeline", Sev::inform) << "Save block mesh to vtk";

    if (m_field.get_comm_rank() == 0) {
      CHKERR schurSaveBlockMesh(block_data_ptr, "block_mesh_hdiv_contact.vtk");
    }

    petsc_mat.reset();
    block_mat.reset();
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}