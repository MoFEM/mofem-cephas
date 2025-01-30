/**
 * @file elastic.cpp
 * @brief elastic example
 * @version 0.13.2
 * @date 2022-09-19
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

constexpr int BASE_DIM = 1; //< Dimension of the base functions

//! [Define dimension]
constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh
//! [Define dimension]
constexpr AssemblyType A = (SCHUR_ASSEMBLE)
                               ? AssemblyType::BLOCK_SCHUR
                               : AssemblyType::PETSC; //< selected assembly type

constexpr IntegrationType I =
    IntegrationType::GAUSS; //< selected integration type

//! [Define entities]
using EntData = EntitiesFieldData::EntData;
using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
//! [Define entities]

//! [OpK]
using OpK = FormsIntegrators<DomainEleOp>::Assembly<A>::BiLinearForm<
    I>::OpGradSymTensorGrad<BASE_DIM, SPACE_DIM, SPACE_DIM, 0>;
//! [OpK]
//! [OpInternalForce]
using OpInternalForce = FormsIntegrators<DomainEleOp>::Assembly<A>::LinearForm<
    I>::OpGradTimesSymTensor<BASE_DIM, SPACE_DIM, SPACE_DIM>;
//! [OpInternalForce]
struct DomainBCs {};
struct BoundaryBCs {};

using DomainRhsBCs = NaturalBC<DomainEleOp>::Assembly<A>::LinearForm<I>;
using OpDomainRhsBCs = DomainRhsBCs::OpFlux<DomainBCs, 1, SPACE_DIM>;
using BoundaryRhsBCs = NaturalBC<BoundaryEleOp>::Assembly<A>::LinearForm<I>;
using OpBoundaryRhsBCs = BoundaryRhsBCs::OpFlux<BoundaryBCs, 1, SPACE_DIM>;
using BoundaryLhsBCs = NaturalBC<BoundaryEleOp>::Assembly<A>::BiLinearForm<I>;
using OpBoundaryLhsBCs = BoundaryLhsBCs::OpFlux<BoundaryBCs, 1, SPACE_DIM>;
using MeshsetSideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;

template <int DIM> struct PostProcEleByDim;

template <> struct PostProcEleByDim<2> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBaseCont<DomainEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBaseCont<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;
  using MeshsetSideEle = PipelineManager::ElementsAndOpsByDim<2>::FaceSideEle;
};

template <> struct PostProcEleByDim<3> {
  using PostProcEleDomain = PostProcBrokenMeshInMoabBaseCont<DomainEle>;
  using PostProcEleBdy = PostProcBrokenMeshInMoabBaseCont<BoundaryEle>;
  using SideEle = PipelineManager::ElementsAndOpsByDim<3>::FaceSideEle;
};

using PostProcEleDomain = PostProcEleByDim<SPACE_DIM>::PostProcEleDomain;
using SideEle = PostProcEleByDim<SPACE_DIM>::SideEle;
using PostProcEleBdy = PostProcEleByDim<SPACE_DIM>::PostProcEleBdy;
using MeshsetSideEleOp = MeshsetSideEle::UserDataOperator;

#include <ElasticSpring.hpp>
#include <FluidLevel.hpp>
#include <CalculateTraction.hpp>
#include <NaturalDomainBC.hpp>
#include <NaturalBoundaryBC.hpp>
#include <RigidBodyTieConstraint.hpp>

constexpr double young_modulus = 1;
constexpr double poisson_ratio = 0.3;
constexpr double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
constexpr double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));

PetscBool is_plane_strain = PETSC_FALSE;

struct Example {

  Example(MoFEM::Interface &m_field) : mField(m_field) {}

  MoFEMErrorCode runProblem();

private:
  MoFEM::Interface &mField;

  PetscBool doEvalField;
  std::array<double, SPACE_DIM> fieldEvalCoords;
  boost::shared_ptr<FieldEvaluatorInterface::SetPtsData> fieldEvalData;
  boost::shared_ptr<MatrixDouble> vectorFieldPtr;

  MoFEMErrorCode readMesh();
  MoFEMErrorCode setupProblem();
  MoFEMErrorCode boundaryCondition();
  MoFEMErrorCode assembleSystem();
  MoFEMErrorCode solveSystem();
  MoFEMErrorCode outputResults();
  MoFEMErrorCode checkResults();

  MoFEMErrorCode addMatBlockOps(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      std::string field_name, std::string block_name,
      boost::shared_ptr<MatrixDouble> mat_D_Ptr, Sev sev);

  struct TieBlock {
    Range tieFaces; // tie faces
    FTensor::Tensor1<double, 3> tieCoord;
    FTensor::Tensor1<double, 3> tieDirection;
  };

  std::vector<TieBlock> tieBlocks; //< Store infomation about tie blocks
};

MoFEMErrorCode Example::addMatBlockOps(
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string field_name, std::string block_name,
    boost::shared_ptr<MatrixDouble> mat_D_Ptr, Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatBlocks : public DomainEleOp {
    OpMatBlocks(std::string field_name, boost::shared_ptr<MatrixDouble> m,
                double bulk_modulus_K, double shear_modulus_G,
                MoFEM::Interface &m_field, Sev sev,
                std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(field_name, DomainEleOp::OPROW), matDPtr(m),
          bulkModulusKDefault(bulk_modulus_K),
          shearModulusGDefault(shear_modulus_G) {
      std::fill(&(doEntities[MBEDGE]), &(doEntities[MBMAXTYPE]), false);
      CHK_THROW_MESSAGE(extractBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          CHKERR getMatDPtr(matDPtr, b.bulkModulusK, b.shearModulusG);
          MoFEMFunctionReturnHot(0);
        }
      }

      CHKERR getMatDPtr(matDPtr, bulkModulusKDefault, shearModulusGDefault);
      MoFEMFunctionReturn(0);
    }

  private:
    boost::shared_ptr<MatrixDouble> matDPtr;

    struct BlockData {
      double bulkModulusK;
      double shearModulusG;
      Range blockEnts;
    };

    double bulkModulusKDefault;
    double shearModulusGDefault;
    std::vector<BlockData> blockData;

    MoFEMErrorCode
    extractBlockData(MoFEM::Interface &m_field,
                     std::vector<const CubitMeshSets *> meshset_vec_ptr,
                     Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() < 2) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        double young_modulus = block_data[0];
        double poisson_ratio = block_data[1];
        double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
        double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));

        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock")
            << "E = " << young_modulus << " nu = " << poisson_ratio;

        blockData.push_back(
            {bulk_modulus_K, shear_modulus_G, get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode getMatDPtr(boost::shared_ptr<MatrixDouble> mat_D_ptr,
                              double bulk_modulus_K, double shear_modulus_G) {
      MoFEMFunctionBegin;
      //! [Calculate elasticity tensor]
      auto set_material_stiffness = [&]() {
        FTensor::Index<'i', SPACE_DIM> i;
        FTensor::Index<'j', SPACE_DIM> j;
        FTensor::Index<'k', SPACE_DIM> k;
        FTensor::Index<'l', SPACE_DIM> l;
        constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
        double A = 1.;
        if (SPACE_DIM == 2 && !is_plane_strain) {
          A = 2 * shear_modulus_G /
              (bulk_modulus_K + (4. / 3.) * shear_modulus_G);
        }
        auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM, 0>(*mat_D_ptr);
        t_D(i, j, k, l) =
            2 * shear_modulus_G * ((t_kd(i, k) ^ t_kd(j, l)) / 4.) +
            A * (bulk_modulus_K - (2. / 3.) * shear_modulus_G) * t_kd(i, j) *
                t_kd(k, l);
      };
      //! [Calculate elasticity tensor]
      constexpr auto size_symm = (SPACE_DIM * (SPACE_DIM + 1)) / 2;
      mat_D_ptr->resize(size_symm * size_symm, 1);
      set_material_stiffness();
      MoFEMFunctionReturn(0);
    }
  };

  pipeline.push_back(new OpMatBlocks(
      field_name, mat_D_Ptr, bulk_modulus_K, shear_modulus_G, mField, sev,

      // Get blockset using regular expression
      mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

//! [Run problem]
MoFEMErrorCode Example::runProblem() {
  MoFEMFunctionBegin;
  CHKERR readMesh();
  CHKERR setupProblem();
  CHKERR boundaryCondition();
  CHKERR assembleSystem();
  CHKERR solveSystem();
  CHKERR outputResults();
  CHKERR checkResults();
  MoFEMFunctionReturn(0);
}
//! [Run problem]

//! [Read mesh]
MoFEMErrorCode Example::readMesh() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  CHKERR simple->getOptions();
  CHKERR simple->loadFile();

  if (A == PETSC) {
    Range tie_ents;
    for (auto m :
         mField.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

             (boost::format("%s(.*)") % "TIE_MATRIX").str()

                 ))

    ) {
      auto meshset = m->getMeshset();
      Range tie_meshset_range;
      CHKERR mField.get_moab().get_entities_by_dimension(
          meshset, SPACE_DIM - 1, tie_meshset_range, true);
      std::vector<double> attributes;
      CHKERR m->getAttributes(attributes);
      if (attributes.size() != 6) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
                 "Wrong number of head parameters %d", attributes.size());
      }
      tieBlocks.push_back({tie_meshset_range,
                           FTensor::Tensor1<double, 3>(
                               attributes[0], attributes[1], attributes[2]),
                           FTensor::Tensor1<double, 3>(
                               attributes[3], attributes[4], attributes[5])});
    }
  }
  MoFEMFunctionReturn(0);
}
//! [Read mesh]

//! [Set up problem]
MoFEMErrorCode Example::setupProblem() {
  MoFEMFunctionBegin;
  Simple *simple = mField.getInterface<Simple>();

  enum bases { AINSWORTH, DEMKOWICZ, LASBASETOPT };
  const char *list_bases[LASBASETOPT] = {"ainsworth", "demkowicz"};
  PetscInt choice_base_value = AINSWORTH;
  CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                              LASBASETOPT, &choice_base_value, PETSC_NULL);

  FieldApproximationBase base;
  switch (choice_base_value) {
  case AINSWORTH:
    base = AINSWORTH_LEGENDRE_BASE;
    MOFEM_LOG("WORLD", Sev::inform)
        << "Set AINSWORTH_LEGENDRE_BASE for displacements";
    break;
  case DEMKOWICZ:
    base = DEMKOWICZ_JACOBI_BASE;
    MOFEM_LOG("WORLD", Sev::inform)
        << "Set DEMKOWICZ_JACOBI_BASE for displacements";
    break;
  default:
    base = LASTBASE;
    break;
  }

  // Add field
  CHKERR simple->addDomainField("U", H1, base, SPACE_DIM);
  CHKERR simple->addBoundaryField("U", H1, base, SPACE_DIM);
  int order = 3;
  CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);

  CHKERR simple->addDataField("GEOMETRY", H1, base, SPACE_DIM);

  CHKERR simple->setFieldOrder("U", order);
  CHKERR simple->setFieldOrder("GEOMETRY", 2);

  //Range rigid_body_ents;
  //EntityHandle rigid_body_meshset;
  Range rigid_body_surface_ents;
  auto add_tie_lagrange_multiplier = [&]() {
    MoFEMFunctionBegin;
    CHKERR simple->addBoundaryField("LAMBDA", H1, base, SPACE_DIM);

    CHKERR simple->addMeshsetField("RIGID_BODY_LAMBDA", NOFIELD, NOBASE, 3);
    CHKERR simple->addMeshsetField("RIGID_BODY_THETA", NOFIELD, NOBASE, 3);
    CHKERR simple->addMeshsetField("LAMBDA", H1, base, SPACE_DIM);

    
    for (auto &t : tieBlocks) {
      // add tie nodes to rigid body meshset
      Range tie_nodes;
      CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 0, false, tie_nodes,
                                               moab::Interface::UNION);
      // add tie edges to rigid body meshset
      Range tie_edges;
      CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 1, false, tie_edges,
                                               moab::Interface::UNION);
      // add to range of rigid body entities
      rigid_body_surface_ents.merge(tie_nodes);
      rigid_body_surface_ents.merge(tie_edges);
      rigid_body_surface_ents.merge(t.tieFaces);
    }
    std::cout << "rigid_body_surface_ents = "
              << rigid_body_surface_ents << std::endl;
    simple->getMeshsetFiniteElementEntities().push_back(
        rigid_body_surface_ents);

    // // Create global fields for rigid body motion
    // CHKERR simple->addDataField("RIGID_BODY_LAMBDA", NOFIELD, NOBASE, 3);
    // CHKERR simple->addDataField("RIGID_BODY_THETA", NOFIELD, NOBASE, 3);

    // CHKERR mField.add_finite_element("RIGID_BODY");

    // auto add_rigid_body_fe_meshset = [&]() {
    //   MoFEMFunctionBegin;

    //   // Define row/cols and element data
    //   CHKERR mField.modify_finite_element_add_field_row("RIGID_BODY",
    //                                                     "RIGID_BODY_LAMBDA");
    //   CHKERR mField.modify_finite_element_add_field_col("RIGID_BODY",
    //                                                     "RIGID_BODY_LAMBDA");

    //   CHKERR mField.modify_finite_element_add_field_row("RIGID_BODY",
    //                                                     "RIGID_BODY_THETA");
    //   CHKERR mField.modify_finite_element_add_field_col("RIGID_BODY",
    //                                                     "RIGID_BODY_THETA");

    //   CHKERR mField.modify_finite_element_add_field_row("RIGID_BODY", "U");
    //   CHKERR mField.modify_finite_element_add_field_col("RIGID_BODY", "U");

    //   CHKERR mField.modify_finite_element_add_field_row("RIGID_BODY", "U");
    //   CHKERR mField.modify_finite_element_add_field_col("RIGID_BODY", "U");

    //   CHKERR mField.modify_finite_element_add_field_row("RIGID_BODY",
    //                                                      "LAMBDA");
    //   CHKERR mField.modify_finite_element_add_field_col("RIGID_BODY",
    //                                                      "LAMBDA");

    //   CHKERR mField.modify_finite_element_add_field_data("RIGID_BODY", "RIGID_BODY_LAMBDA");
    //   CHKERR mField.modify_finite_element_add_field_data("RIGID_BODY", "RIGID_BODY_THETA");
    //   CHKERR mField.modify_finite_element_add_field_data("RIGID_BODY", "U");
    //   CHKERR mField.modify_finite_element_add_field_data("RIGID_BODY", "LAMBDA");

    //   auto translation_meshset = mField.get_field_meshset("RIGID_BODY_LAMBDA");
    //   auto rotation_meshset = mField.get_field_meshset("RIGID_BODY_THETA");

    //   CHKERR mField.get_moab().get_entities_by_handle(translation_meshset,
    //                                                   rigid_body_ents, true);
    //   CHKERR mField.get_moab().get_entities_by_handle(rotation_meshset,
    //                                                   rigid_body_ents, true);

    //   {
    //     CHKERR mField.get_moab().create_meshset(MESHSET_SET,
    //                                             rigid_body_meshset);
    //     CHKERR mField.get_moab().add_entities(rigid_body_meshset, rigid_body_ents);

    //     // This implementation assumes only one rigid body
    //     // add tie faces to rigid body meshset
    //     for (auto &t : tieBlocks) {
    //       CHKERR mField.get_moab().add_entities(rigid_body_meshset, t.tieFaces);
    //       // add tie nodes to rigid body meshset
    //       Range tie_nodes;
    //       CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 0, false, tie_nodes,
    //                                               moab::Interface::UNION);
    //       CHKERR mField.get_moab().add_entities(rigid_body_meshset, tie_nodes);
    //       // add tie edges to rigid body meshset
    //       Range tie_edges;
    //       CHKERR mField.get_moab().get_adjacencies(t.tieFaces, 1, false, tie_edges,
    //                                               moab::Interface::UNION);
    //       CHKERR mField.get_moab().add_entities(rigid_body_meshset, tie_edges);
    //       // add to range of rigid body entities
    //       rigid_body_ents.merge(tie_nodes);
    //       rigid_body_ents.merge(tie_edges);
    //       rigid_body_ents.merge(t.tieFaces);

    //     }

    //     CHKERR mField.getInterface<BitRefManager>()->setBitLevelToMeshset(
    //         rigid_body_meshset, BitRefLevel().set());
    //   }
    //   //CHKERR mField.add_ents_to_finite_element_by_MESHSET(rigid_body_meshset,
    //   //                                                    "RIGID_BODY", false);

    //   MoFEMFunctionReturn(0);
    // };


    // CHKERR add_rigid_body_fe_meshset();

    // simple->getOtherFiniteElements().push_back("RIGID_BODY");

    CHKERR simple->setFieldOrder("LAMBDA", 0);
    // CHKERR simple->setFieldOrder("LAMBDA_ROT", 0);
    for (auto &t : tieBlocks) {
      CHKERR simple->setFieldOrder("LAMBDA", order, &t.tieFaces);
      // CHKERR simple->setFieldOrder("LAMBDA_ROT", order, &t.tieFaces);
    }
    MoFEMFunctionReturn(0);
  };

  if (A == PETSC) {
    CHKERR add_tie_lagrange_multiplier();
  }
  CHKERR simple->defineFiniteElements();

  // auto add_no_field = [&](auto fe_name, auto field) {
  //   MoFEMFunctionBegin;
  //   CHKERR mField.modify_finite_element_add_field_row(fe_name, field);
  //   CHKERR mField.modify_finite_element_add_field_col(fe_name, field);
  //   CHKERR mField.modify_finite_element_add_field_data(fe_name, field);
  //   MoFEMFunctionReturn(0);
  // };
  // CHKERR add_no_field("bFE", "RIGID_BODY_LAMBDA");
  // CHKERR add_no_field("bFE", "RIGID_BODY_THETA");
  //   auto add_no_field = [&](auto fe_name, auto field) {
  //   MoFEMFunctionBegin;
  //   CHKERR mField.modify_finite_element_add_field_row(fe_name, field);
  //   CHKERR mField.modify_finite_element_add_field_col(fe_name, field);
  //   CHKERR mField.modify_finite_element_add_field_data(fe_name, field);
  //   MoFEMFunctionReturn(0);
  // };
  // CHKERR add_no_field("mFE", "LAMBDA");

  CHKERR simple->defineProblem(PETSC_TRUE);
  CHKERR simple->buildFields();
  // add vertex to fe here from meshsets
  // Range rigid_body_ents;
  // auto translation_meshset = mField.get_field_meshset("RIGID_BODY_LAMBDA");
  // auto rotation_meshset = mField.get_field_meshset("RIGID_BODY_THETA");
  // CHKERR mField.get_moab().get_entities_by_handle(translation_meshset,
  //                                                 rigid_body_ents, true);
  // CHKERR mField.get_moab().get_entities_by_handle(rotation_meshset,
  //                                                 rigid_body_ents, true);

  // std::cout << "rigid_body_ents = " << rigid_body_ents << std::endl;

  // CHKERR mField.getInterface<BitRefManager>()->setBitLevelToMeshset(
  //     rigid_body_meshset, BitRefLevel().set());

  // CHKERR mField.add_ents_to_finite_element_by_MESHSET(rigid_body_meshset,
  //                                                     "RIGID_BODY", false);


  CHKERR simple->buildFiniteElements();
  CHKERR simple->buildProblem();

  // std::cout<< "Rigid body surface ents = " << rigid_body_surface_ents << std::endl;
  // // Remove prescirbed dofs from rigid body assuming just x translation
  // CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
  //     simple->getProblemName(), "RIGID_BODY_LAMBDA", rigid_body_surface_ents, 0, 0);
  // CHKERR mField.getInterface<ProblemsManager>()->removeDofsOnEntities(
  //     simple->getProblemName(), "RIGID_BODY_LAMBDA", rigid_body_ents, 0, 0);

  // auto set_translation = [this](boost::shared_ptr<FieldEntity>
  // field_entity_ptr) {
  //   for(auto &t : tieBlocks) {
  //     // Scale with time?
  //     field_entity_ptr->getEntFieldData()[0] = t.tieDirection(0);
  //   }

  //   return 0;
  // };
  // CHKERR
  // mField.getInterface<FieldBlas>()->fieldLambdaOnEntities(set_translation,
  //                                                                "RIGID_BODY_LAMBDA");

  // CHKERR simple->setUp();

  auto project_ho_geometry = [&]() {
    Projection10NodeCoordsOnField ent_method(mField, "GEOMETRY");
    return mField.loop_dofs("GEOMETRY", ent_method);
  };
  CHKERR project_ho_geometry();

  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-plane_strain", &is_plane_strain,
                             PETSC_NULL);

  int coords_dim = SPACE_DIM;
  CHKERR PetscOptionsGetRealArray(NULL, NULL, "-field_eval_coords",
                                  fieldEvalCoords.data(), &coords_dim,
                                  &doEvalField);

  if (doEvalField) {
    vectorFieldPtr = boost::make_shared<MatrixDouble>();
    fieldEvalData =
        mField.getInterface<FieldEvaluatorInterface>()->getData<DomainEle>();

    if constexpr (SPACE_DIM == 3) {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->buildTree3D(
          fieldEvalData, simple->getDomainFEName());
    } else {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->buildTree2D(
          fieldEvalData, simple->getDomainFEName());
    }

    fieldEvalData->setEvalPoints(fieldEvalCoords.data(), 1);
    auto no_rule = [](int, int, int) { return -1; };
    auto field_eval_fe_ptr = fieldEvalData->feMethodPtr.lock();
    field_eval_fe_ptr->getRuleHook = no_rule;

    field_eval_fe_ptr->getOpPtrVector().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", vectorFieldPtr));
    }

    MoFEMFunctionReturn(0);
  }
//! [Set up problem]

//! [Boundary condition]
MoFEMErrorCode Example::boundaryCondition() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto bc_mng = mField.getInterface<BcManager>();

  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_X",
                                           "U", 0, 0);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Y",
                                           "U", 1, 1);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(), "REMOVE_Z",
                                           "U", 2, 2);
  CHKERR bc_mng->removeBlockDOFsOnEntities(simple->getProblemName(),
                                           "REMOVE_ALL", "U", 0, 3);
  CHKERR bc_mng->pushMarkDOFsOnEntities<DisplacementCubitBcData>(
      simple->getProblemName(), "U");

  // adding MPCs
  CHKERR bc_mng->addBlockDOFsToMPCs(simple->getProblemName(), "U");

  MoFEMFunctionReturn(0);
}
//! [Boundary condition]

//! [Push operators to pipeline]
MoFEMErrorCode Example::assembleSystem() {
  MoFEMFunctionBegin;
  auto pip = mField.getInterface<PipelineManager>();
  auto simple = mField.getInterface<Simple>();

  //! [Integration rule]
  auto integration_rule = [](int, int, int approx_order) {
    return 2 * approx_order + 1;
  };

  auto integration_rule_bc = [](int, int, int approx_order) {
    return 2 * approx_order + 1;
  };

  CHKERR pip->setDomainRhsIntegrationRule(integration_rule);
  CHKERR pip->setDomainLhsIntegrationRule(integration_rule);
  CHKERR pip->setBoundaryRhsIntegrationRule(integration_rule_bc);
  CHKERR pip->setBoundaryLhsIntegrationRule(integration_rule_bc);
  //! [Integration rule]

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpDomainLhsPipeline(), {H1}, "GEOMETRY");
  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpDomainRhsPipeline(), {H1}, "GEOMETRY");
  CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpBoundaryRhsPipeline(), {NOSPACE}, "GEOMETRY");
  CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpBoundaryLhsPipeline(), {NOSPACE}, "GEOMETRY");

  //! [Push domain stiffness matrix to pipeline]
  auto mat_D_ptr = boost::make_shared<MatrixDouble>();

  // Assemble domain stiffness matrix
  CHKERR addMatBlockOps(pip->getOpDomainLhsPipeline(), "U", "MAT_ELASTIC",
                        mat_D_ptr, Sev::verbose);
  pip->getOpDomainLhsPipeline().push_back(new OpK("U", "U", mat_D_ptr));
  //! [Push domain stiffness matrix to pipeline]

  //! [Push Internal forces]
  auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
  auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
  auto mat_stress_ptr = boost::make_shared<MatrixDouble>();
  pip->getOpDomainRhsPipeline().push_back(
      new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("U",
                                                               mat_grad_ptr));
  CHKERR addMatBlockOps(pip->getOpDomainRhsPipeline(), "U", "MAT_ELASTIC",
                        mat_D_ptr, Sev::inform);
  pip->getOpDomainRhsPipeline().push_back(
      new OpSymmetrizeTensor<SPACE_DIM>(mat_grad_ptr, mat_strain_ptr));
  pip->getOpDomainRhsPipeline().push_back(
      new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
          mat_strain_ptr, mat_stress_ptr, mat_D_ptr));

  pip->getOpDomainRhsPipeline().push_back(
      new OpInternalForce("U", mat_stress_ptr,
                          [](double, double, double) constexpr { return -1; }));
  //! [Push Internal forces]



  //! [Constraint matrix]
  if (A == PETSC) {
    auto add_constrain_lhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      auto u_ptr = boost::make_shared<MatrixDouble>();
      auto lambda_ptr = boost::make_shared<MatrixDouble>();
      auto translation_ptr = boost::make_shared<VectorDouble>();
      auto theta_ptr = boost::make_shared<VectorDouble>();

      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
      pip.push_back(new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA",
                                                       translation_ptr));
      pip.push_back(
          new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

      for (auto &t : tieBlocks) {
        pip.push_back(new OpTieTermConstrainRigidBodyLhs_dU<SPACE_DIM>(
            "LAMBDA", "U", u_ptr, t.tieCoord, t.tieDirection,
            boost::make_shared<Range>(t.tieFaces), translation_ptr, theta_ptr));
        pip.push_back(
            new OpTieTermConstrainRigidBodyLhs_dTranslation<SPACE_DIM>(
                "LAMBDA", "RIGID_BODY_LAMBDA", u_ptr, t.tieCoord,
                t.tieDirection, boost::make_shared<Range>(t.tieFaces)));
        pip.push_back(new OpTieTermConstrainRigidBodyLhs_dRotation<SPACE_DIM>(
            "LAMBDA", "RIGID_BODY_THETA", u_ptr, t.tieCoord, t.tieDirection,
            boost::make_shared<Range>(t.tieFaces)));
      }
      MoFEMFunctionReturn(0);
    };

    auto add_constrain_rhs = [&](auto &pip) {
      MoFEMFunctionBegin;
      auto u_ptr = boost::make_shared<MatrixDouble>();
      auto lambda_ptr = boost::make_shared<MatrixDouble>();
      auto translation_ptr = boost::make_shared<VectorDouble>();
      auto theta_ptr = boost::make_shared<VectorDouble>();

      pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
      pip.push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));

      pip.push_back(new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA",
                                                       translation_ptr));
      pip.push_back(
          new OpCalculateNoFieldVectorValues("RIGID_BODY_THETA", theta_ptr));

      for (auto &t : tieBlocks) {
        pip.push_back(new OpTieTermConstrainRigidBodyRhs<SPACE_DIM>(
            "LAMBDA", u_ptr, t.tieCoord, t.tieDirection,
            boost::make_shared<Range>(t.tieFaces), translation_ptr, theta_ptr));
      }
      MoFEMFunctionReturn(0);
    };


    CHKERR add_constrain_lhs(pip->getOpBoundaryLhsPipeline());
    CHKERR add_constrain_rhs(pip->getOpBoundaryRhsPipeline());


    // Add Rigid body element to pipeline
    auto add_rigid_body_rhs = [&](auto &pip) {
      MoFEMFunctionBegin;

      // add ops
      auto lambda_ptr = boost::make_shared<MatrixDouble>();
      auto translation_ptr = boost::make_shared<VectorDouble>();
      auto u_ptr = boost::make_shared<MatrixDouble>();
      auto int_translation_ptr =
          boost::make_shared<VectorDouble>(VectorDouble(3));
      auto int_rotation_ptr =
          boost::make_shared<MatrixDouble>(MatrixDouble(3, 3));

      Range rigid_body_ents;
      auto v_rigid_body_ents = simple->getMeshsetFiniteElementEntities();

      // loop over all rigid body entities
      for(auto ents : v_rigid_body_ents) {
        rigid_body_ents.merge(ents);
      }



      auto rigid_body_ents_ptr = boost::make_shared<Range>(rigid_body_ents);

      for (auto &t : tieBlocks) {
        auto op_loop_side = new OpLoopSide<BoundaryEle>(
            mField, "bFE", SPACE_DIM - 1, rigid_body_ents_ptr);
        CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
            op_loop_side->getOpPtrVector(), {});
        op_loop_side->getOpPtrVector().push_back(
            new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
        op_loop_side->getOpPtrVector().push_back(
            new OpTieTermConstrainRigidBodyGlobalTranslationIntegralRhs<
                SPACE_DIM>("LAMBDA", u_ptr, t.tieCoord, t.tieDirection,
                           boost::make_shared<Range>(t.tieFaces), lambda_ptr,
                           int_translation_ptr));
        op_loop_side->getOpPtrVector().push_back(
            new OpTieTermConstrainRigidBodyGlobalRotationIntegralRhs<SPACE_DIM>(
                "LAMBDA", t.tieCoord, lambda_ptr, int_rotation_ptr));
        pip.push_back(op_loop_side);

        pip.push_back(
            new OpTieTermConstrainRigidBodyGlobalTranslationRhs<SPACE_DIM>(
                "RIGID_BODY_LAMBDA", t.tieCoord, t.tieDirection,
                boost::make_shared<Range>(t.tieFaces), lambda_ptr,
                int_translation_ptr));
        pip.push_back(
            new OpTieTermConstrainRigidBodyGlobalRotationRhs<SPACE_DIM>(
                "RIGID_BODY_THETA", int_rotation_ptr));
      }
      // MoFEMFunctionReturn(0);
      //};

      // CHKERR add_rigid_body_ops(fe_ptr);
      MoFEMFunctionReturn(0);
    };

    // auto add_rigid_body_lhs = [&](auto &pip) {
    //   MoFEMFunctionBegin;
    //   auto lambda_ptr = boost::make_shared<MatrixDouble>();
    //   auto translation_ptr = boost::make_shared<VectorDouble>();

    //   for (auto &t : tieBlocks) {
    //     pip.push_back(new OpTieTermSetBc<SPACE_DIM>(
    //         "RIGID_BODY_LAMBDA", "RIGID_BODY_LAMBDA", boost::make_shared<Range>(t.tieFaces)));
    //   }
    //   MoFEMFunctionReturn(0);
    // };

    // CHKERR add_rigid_body_lhs(pip->getOpMeshsetLhsPipeline());

    CHKERR add_rigid_body_rhs(pip->getOpMeshsetRhsPipeline());
  }
  //! [Constraint matrix]

  //! [Push Body forces]
  CHKERR DomainRhsBCs::AddFluxToPipeline<OpDomainRhsBCs>::add(
      pip->getOpDomainRhsPipeline(), mField, "U", Sev::inform);
  //! [Push Body forces]

  //! [Push natural boundary conditions]
  // Add force boundary condition
  CHKERR BoundaryRhsBCs::AddFluxToPipeline<OpBoundaryRhsBCs>::add(
      pip->getOpBoundaryRhsPipeline(), mField, "U", -1, Sev::inform);
  // Add case for mix type of BCs
  CHKERR BoundaryLhsBCs::AddFluxToPipeline<OpBoundaryLhsBCs>::add(
      pip->getOpBoundaryLhsPipeline(), mField, "U", Sev::verbose);
  //! [Push natural boundary conditions]
  MoFEMFunctionReturn(0);
}
//! [Push operators to pipeline]

struct SetUpSchur {
  static boost::shared_ptr<SetUpSchur>
  createSetUpSchur(MoFEM::Interface &m_field);
  virtual MoFEMErrorCode setUp(SmartPetscObj<KSP> solver) = 0;

protected:
  SetUpSchur() = default;
};

//! [Solve]
MoFEMErrorCode Example::solveSystem() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();
  auto solver = pip->createKSP();
  CHKERR KSPSetFromOptions(solver);

  auto dm = simple->getDM();
  auto D = createDMVector(dm);
  auto F = vectorDuplicate(D);

  auto set_essential_bc = [&]() {
    MoFEMFunctionBegin;
    // This is low level pushing finite elements (pipelines) to solver
    auto ksp_ctx_ptr = getDMKspCtx(dm);

    auto pre_proc_rhs = boost::make_shared<FEMethod>();
    auto post_proc_rhs = boost::make_shared<FEMethod>();
    auto post_proc_lhs = boost::make_shared<FEMethod>();

    auto get_pre_proc_hook = [&]() {
      return EssentialPreProc<DisplacementCubitBcData>(mField, pre_proc_rhs,
                                                       {});
    };


    auto get_pre_proc_tie_hook = [&](auto &tie_direction) {
      return SetTieBcPreProc(mField, pre_proc_rhs, tie_direction);
    };

    pre_proc_rhs->preProcessHook = get_pre_proc_hook();

    auto get_post_proc_hook_rhs = [this, post_proc_rhs]() {
      MoFEMFunctionBeginHot;

      CHKERR EssentialPostProcRhs<DisplacementCubitBcData>(mField,
                                                           post_proc_rhs, 1.)();
      CHKERR EssentialPostProcRhs<MPCsType>(mField, post_proc_rhs, 1.)();
      MoFEMFunctionReturnHot(0);
    };

    auto get_post_proc_hook_lhs = [this, post_proc_lhs]() {
      MoFEMFunctionBeginHot;

      CHKERR EssentialPostProcLhs<DisplacementCubitBcData>(mField,
                                                           post_proc_lhs, 1.)();
      CHKERR EssentialPostProcLhs<MPCsType>(mField, post_proc_lhs)();
      MoFEMFunctionReturnHot(0);
    };

    post_proc_rhs->postProcessHook = get_post_proc_hook_rhs;
    post_proc_lhs->postProcessHook = get_post_proc_hook_lhs;

    ksp_ctx_ptr->getPreProcComputeRhs().push_front(pre_proc_rhs);
    ksp_ctx_ptr->getPostProcComputeRhs().push_back(post_proc_rhs);
    ksp_ctx_ptr->getPostProcSetOperators().push_back(post_proc_lhs);
    MoFEMFunctionReturn(0);
  };

  auto set_global_bc = [&]() {
    MoFEMFunctionBegin; 
    auto ksp_ctx_ptr = getDMKspCtx(dm);
    auto fe_pre_proc_rhs = boost::make_shared<FEMethod>();
    auto fe_post_proc_rhs = boost::make_shared<FEMethod>();
    auto fe_post_proc_lhs = boost::make_shared<FEMethod>();
    

    auto set_pre_rhs = [this, fe_pre_proc_rhs, &dm]() {
      MoFEMFunctionBeginHot;
      auto is_mng = mField.getInterface<ISManager>();
      auto simple = mField.getInterface<Simple>();
      SmartPetscObj<IS> is;
      CHKERR is_mng->isCreateProblemFieldAndRankLocal(
          simple->getProblemName(), ROW, "RIGID_BODY_LAMBDA", 0,
          MAX_DOFS_ON_ENTITY, is, nullptr);
      const int *i_ptr;
      CHKERR ISGetIndices(is, &i_ptr);
      int size;
      CHKERR ISGetLocalSize(is, &size);
      std::cout << "size = " << size << std::endl;
      // if (size != SPACE_DIM)
      //   SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Expected
      //   SPACE_DIM");
      double *x_ptr;

      if (fe_pre_proc_rhs->x == NULL) {
        fprintf(stderr, "Error: fe_pre_proc_rhs is NULL\n");
      }
      CHKERR VecGetArray(fe_pre_proc_rhs->x, &x_ptr);
      for (auto i = 0; i != size; ++i) {
        if (i % SPACE_DIM == 0) { // set x
          x_ptr[i_ptr[i]] = 0.01; // times (fe_pre_proc_rhs->ts_t);
        }
      }
      CHKERR VecRestoreArray(fe_pre_proc_rhs->x, &x_ptr);
      CHKERR ISRestoreIndices(is, &i_ptr);
      CHKERR DMoFEMMeshToGlobalVector(dm, fe_pre_proc_rhs->x, INSERT_VALUES,
                                      SCATTER_REVERSE);
      MoFEMFunctionReturnHot(0);
    };

    auto get_pre_proc_tie_hook = [&](auto &tie_direction) {
      return SetTieBcPreProc(mField, fe_pre_proc_rhs, tie_direction);
    };

    auto set_post_rhs = [this, fe_post_proc_rhs]() {
      MoFEMFunctionBeginHot;
      auto is_mng = mField.getInterface<ISManager>();
      auto simple = mField.getInterface<Simple>();
      SmartPetscObj<IS> is;
      CHKERR is_mng->isCreateProblemFieldAndRankLocal(
          simple->getProblemName(), ROW, "RIGID_BODY_LAMBDA", 0,
          MAX_DOFS_ON_ENTITY, is, nullptr);
          
      const int *i_ptr;
      CHKERR ISGetIndices(is, &i_ptr);
      int size;
      CHKERR ISGetLocalSize(is, &size);
      std::cout << "size = " << size << std::endl;
      double *f_ptr;
      CHKERR VecGetArray(fe_post_proc_rhs->f, &f_ptr);
      for (auto i = 0; i != size; ++i) {
        if (i % SPACE_DIM == 0) {
          //f_ptr[i_ptr[i]] = x_ptr[i_ptr[i]] - f_ptr[i_ptr[i]]
          f_ptr[i_ptr[i]] = 0.1;
        }
      }
      CHKERR VecRestoreArray(fe_post_proc_rhs->f, &f_ptr);
      CHKERR ISRestoreIndices(is, &i_ptr);
      MoFEMFunctionReturnHot(0);
    };



    auto set_post_lhs = [this, fe_post_proc_lhs]() {
      MoFEMFunctionBeginHot;
      auto is_mng = mField.getInterface<ISManager>();
      auto simple = mField.getInterface<Simple>();
      SmartPetscObj<IS> is;
      CHKERR is_mng->isCreateProblemFieldAndRankLocal(
          simple->getProblemName(), ROW, "RIGID_BODY_LAMBDA", 0, SPACE_DIM, is,
          nullptr);
      *(fe_post_proc_lhs->matAssembleSwitch) = PETSC_TRUE;
      CHKERR MatAssemblyBegin(fe_post_proc_lhs->B, MAT_FINAL_ASSEMBLY);
      CHKERR MatAssemblyEnd(fe_post_proc_lhs->B, MAT_FINAL_ASSEMBLY);
      const int *i_ptr;
      CHKERR ISGetIndices(is, &i_ptr);
      const int idx = 0;
      // CHKERR MatZeroRowsColumns(fe_post_proc_lhs->B, 1, &i_ptr[idx], 1, PETSC_NULL,
      //                           PETSC_NULL);
      CHKERR MatZeroRows(fe_post_proc_lhs->B, 1, &i_ptr[idx], 1, PETSC_NULL,
                         PETSC_NULL);
      CHKERR ISRestoreIndices(is, &i_ptr);

      MoFEMFunctionReturnHot(0);
    };

    // for (auto &t : tieBlocks) {
    //   fe_post_proc_rhs->preProcessHook = get_pre_proc_tie_hook(t.tieDirection);
    // }
    fe_post_proc_rhs->postProcessHook = set_post_rhs;
    //ksp_ctx_ptr->getPreProcComputeRhs().push_back(fe_pre_proc_rhs);
    ksp_ctx_ptr->getPostProcComputeRhs().push_back(fe_post_proc_rhs);
    fe_post_proc_lhs->postProcessHook = set_post_lhs;
    ksp_ctx_ptr->getPostProcSetOperators().push_back(fe_post_proc_lhs);
    MoFEMFunctionReturn(0);
  };

  auto setup_and_solve = [&]() {
    MoFEMFunctionBegin;
    BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
    MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp";
    CHKERR KSPSetUp(solver);
    MOFEM_LOG("TIMER", Sev::inform) << "KSPSetUp <= Done";
    MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve";
    CHKERR KSPSolve(solver, F, D);
    MOFEM_LOG("TIMER", Sev::inform) << "KSPSolve <= Done";
    MoFEMFunctionReturn(0);
  };

  MOFEM_LOG_CHANNEL("TIMER");
  MOFEM_LOG_TAG("TIMER", "timer");

  CHKERR set_essential_bc();
  CHKERR set_global_bc();

  if (A == AssemblyType::BLOCK_SCHUR || A == AssemblyType::SCHUR) {
    auto schur_ptr = SetUpSchur::createSetUpSchur(mField);
    CHKERR schur_ptr->setUp(solver);
    CHKERR setup_and_solve();
  } else {
    CHKERR setup_and_solve();
  }

  CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);

  if (doEvalField) {
    if constexpr (SPACE_DIM == 3) {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->evalFEAtThePoint3D(
          fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
          simple->getDomainFEName(), fieldEvalData, mField.get_comm_rank(),
          mField.get_comm_rank(), nullptr, MF_EXIST, QUIET);
    } else {
      CHKERR mField.getInterface<FieldEvaluatorInterface>()->evalFEAtThePoint2D(
          fieldEvalCoords.data(), 1e-12, simple->getProblemName(),
          simple->getDomainFEName(), fieldEvalData, mField.get_comm_rank(),
          mField.get_comm_rank(), nullptr, MF_EXIST, QUIET);
    }

    if (vectorFieldPtr->size1()) {
      auto t_disp = getFTensor1FromMat<SPACE_DIM>(*vectorFieldPtr);
      if constexpr (SPACE_DIM == 2)
        MOFEM_LOG("FieldEvaluator", Sev::inform)
            << "U_X: " << t_disp(0) << " U_Y: " << t_disp(1);
      else
        MOFEM_LOG("FieldEvaluator", Sev::inform)
            << "U_X: " << t_disp(0) << " U_Y: " << t_disp(1)
            << " U_Z: " << t_disp(2);
    }

    MOFEM_LOG_SYNCHRONISE(mField.get_comm());
  }

  MoFEMFunctionReturn(0);
}
//! [Solve]

//! [Postprocess results]
MoFEMErrorCode Example::outputResults() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();
  auto det_ptr = boost::make_shared<VectorDouble>();
  auto jac_ptr = boost::make_shared<MatrixDouble>();
  auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
  //! [Postprocess clean]
  pip->getDomainRhsFE().reset();
  pip->getDomainLhsFE().reset();
  pip->getBoundaryRhsFE().reset();
  pip->getBoundaryLhsFE().reset();
  //! [Postprocess clean]

  //! [Postprocess initialise]
  auto post_proc_mesh = boost::make_shared<moab::Core>();
  auto post_proc_begin = boost::make_shared<PostProcBrokenMeshInMoabBaseBegin>(
      mField, post_proc_mesh);
  auto post_proc_end = boost::make_shared<PostProcBrokenMeshInMoabBaseEnd>(
      mField, post_proc_mesh);
  //! [Postprocess initialise]

  auto calculate_stress_ops = [&](auto &pip) {
    AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(pip, {H1}, "GEOMETRY");
    auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
    auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
    auto mat_stress_ptr = boost::make_shared<MatrixDouble>();
    pip.push_back(new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>(
        "U", mat_grad_ptr));
    auto mat_D_ptr = boost::make_shared<MatrixDouble>();
    CHKERR addMatBlockOps(pip, "U", "MAT_ELASTIC", mat_D_ptr, Sev::verbose);
    pip.push_back(
        new OpSymmetrizeTensor<SPACE_DIM>(mat_grad_ptr, mat_strain_ptr));
    pip.push_back(new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
        mat_strain_ptr, mat_stress_ptr, mat_D_ptr));
    auto u_ptr = boost::make_shared<MatrixDouble>();
    pip.push_back(new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));
    auto x_ptr = boost::make_shared<MatrixDouble>();
    pip.push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("GEOMETRY", x_ptr));
    return boost::make_tuple(u_ptr, x_ptr, mat_strain_ptr, mat_stress_ptr);
  };

  auto post_proc_domain = [&](auto post_proc_mesh) {
    auto post_proc_fe =
        boost::make_shared<PostProcEleDomain>(mField, post_proc_mesh);
    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    auto [u_ptr, x_ptr, mat_strain_ptr, mat_stress_ptr] =
        calculate_stress_ops(post_proc_fe->getOpPtrVector());

    post_proc_fe->getOpPtrVector().push_back(

        new OpPPMap(

            post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

            {},

            {{"U", u_ptr}, {"GEOMETRY", x_ptr}},

            {},

            {{"STRAIN", mat_strain_ptr}, {"STRESS", mat_stress_ptr}}

            )

    );
    return post_proc_fe;
  };

  auto post_proc_boundary = [&](auto post_proc_mesh) {
    auto post_proc_fe =
        boost::make_shared<PostProcEleBdy>(mField, post_proc_mesh);
    AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
        post_proc_fe->getOpPtrVector(), {}, "GEOMETRY");
    auto op_loop_side =
        new OpLoopSide<SideEle>(mField, simple->getDomainFEName(), SPACE_DIM);
    // push ops to side element, through op_loop_side operator
    auto [u_ptr, x_ptr, mat_strain_ptr, mat_stress_ptr] =
        calculate_stress_ops(op_loop_side->getOpPtrVector());
    post_proc_fe->getOpPtrVector().push_back(op_loop_side);
    auto mat_traction_ptr = boost::make_shared<MatrixDouble>();
    post_proc_fe->getOpPtrVector().push_back(
        new ElasticExample::OpCalculateTraction(mat_stress_ptr,
                                                mat_traction_ptr));

    using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

    post_proc_fe->getOpPtrVector().push_back(

        new OpPPMap(

            post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

            {},

            {{"U", u_ptr}, {"GEOMETRY", x_ptr}, {"T", mat_traction_ptr}},

            {},

            {{"STRAIN", mat_strain_ptr}, {"STRESS", mat_stress_ptr}}

            )

    );
    return post_proc_fe;
  };

  PetscBool post_proc_skin_only = PETSC_FALSE;
  if (SPACE_DIM == 3) {
    post_proc_skin_only = PETSC_TRUE;
    CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-post_proc_skin_only",
                               &post_proc_skin_only, PETSC_NULL);
  }
  if (post_proc_skin_only == PETSC_FALSE) {
    pip->getDomainRhsFE() = post_proc_domain(post_proc_mesh);
  } else {
    pip->getBoundaryRhsFE() = post_proc_boundary(post_proc_mesh);
  }

  CHKERR DMoFEMPreProcessFiniteElements(simple->getDM(),
                                        post_proc_begin->getFEMethod());
  CHKERR pip->loopFiniteElements();
  CHKERR DMoFEMPostProcessFiniteElements(simple->getDM(),
                                         post_proc_end->getFEMethod());

  CHKERR post_proc_end->writeFile("out_elastic.h5m");
  MoFEMFunctionReturn(0);
}
//! [Postprocess results]

//! [Check]
MoFEMErrorCode Example::checkResults() {
  MOFEM_LOG_CHANNEL("WORLD");
  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();
  MoFEMFunctionBegin;
  pip->getDomainRhsFE().reset();
  pip->getDomainLhsFE().reset();
  pip->getBoundaryRhsFE().reset();
  pip->getBoundaryLhsFE().reset();

  auto integration_rule = [](int, int, int p_data) { return 2 * p_data + 1; };
  CHKERR pip->setDomainRhsIntegrationRule(integration_rule);
  CHKERR pip->setBoundaryRhsIntegrationRule(integration_rule);

  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpDomainRhsPipeline(), {H1}, "GEOMETRY");
  CHKERR AddHOOps<SPACE_DIM - 1, SPACE_DIM, SPACE_DIM>::add(
      pip->getOpBoundaryRhsPipeline(), {}, "GEOMETRY");

  auto mat_grad_ptr = boost::make_shared<MatrixDouble>();
  auto mat_strain_ptr = boost::make_shared<MatrixDouble>();
  auto mat_stress_ptr = boost::make_shared<MatrixDouble>();

  pip->getOpDomainRhsPipeline().push_back(
      new OpCalculateVectorFieldGradient<SPACE_DIM, SPACE_DIM>("U",
                                                               mat_grad_ptr));
  pip->getOpDomainRhsPipeline().push_back(
      new OpSymmetrizeTensor<SPACE_DIM>(mat_grad_ptr, mat_strain_ptr));

  auto mat_D_ptr = boost::make_shared<MatrixDouble>();
  CHKERR addMatBlockOps(pip->getOpDomainRhsPipeline(), "U", "MAT_ELASTIC",
                        mat_D_ptr, Sev::verbose);
  pip->getOpDomainRhsPipeline().push_back(
      new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
          mat_strain_ptr, mat_stress_ptr, mat_D_ptr));

  pip->getOpDomainRhsPipeline().push_back(
      new OpInternalForce("U", mat_stress_ptr));

  pip->getOpBoundaryRhsPipeline().push_back(
      new OpTensorTimesSymmetricTensor<SPACE_DIM, SPACE_DIM>(
          mat_strain_ptr, mat_stress_ptr, mat_D_ptr));
  CHKERR DomainRhsBCs::AddFluxToPipeline<OpDomainRhsBCs>::add(
      pip->getOpDomainRhsPipeline(), mField, "U", Sev::verbose);
  CHKERR BoundaryRhsBCs::AddFluxToPipeline<OpBoundaryRhsBCs>::add(
      pip->getOpBoundaryRhsPipeline(), mField, "U", -1, Sev::verbose);

  if (A == AssemblyType::PETSC) {

    auto lambda_ptr = boost::make_shared<MatrixDouble>();
    //auto lambda_rotation_ptr = boost::make_shared<VectorDouble>();
    auto u_ptr = boost::make_shared<MatrixDouble>();
    auto translation_ptr = boost::make_shared<VectorDouble>();

    pip->getOpBoundaryRhsPipeline().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_ptr));

    pip->getOpBoundaryRhsPipeline().push_back(
        new OpCalculateVectorFieldValues<SPACE_DIM>("LAMBDA", lambda_ptr));
    pip->getOpBoundaryRhsPipeline().push_back(
        new OpCalculateNoFieldVectorValues("RIGID_BODY_LAMBDA", translation_ptr));    
    //pip->getOpBoundaryRhsPipeline().push_back(
        //new OpCalculateScalarFieldValues("LAMBDA_ROT", lambda_rotation_ptr));

    for (auto &t : tieBlocks) {
      // pip->getOpBoundaryRhsPipeline().push_back(new OpTieSetDofs("RIGID_BODY_LAMBDA", translation_ptr,
      //                                t.tieDirection));
      pip->getOpBoundaryRhsPipeline().push_back(
          new OpTieTermConstrainRigidBodyRhs_du<SPACE_DIM>(
              "U", "LAMBDA", u_ptr, lambda_ptr, t.tieCoord, t.tieDirection,
              boost::make_shared<Range>(t.tieFaces)));
      // pip->getOpBoundaryRhsPipeline().push_back(
      //     new OpTieTermConstrainRotationNormalRhs_du<SPACE_DIM>(
      //         "U", "LAMBDA_ROT", u_ptr, lambda_rotation_ptr, t.tieCoord, t.tieDirection,
      //         boost::make_shared<Range>(t.tieFaces), "y"));
    }
  }
  auto dm = simple->getDM();
  auto res = createDMVector(dm);
  CHKERR VecSetDM(res, PETSC_NULL);

  pip->getDomainRhsFE()->f = res;
  pip->getBoundaryRhsFE()->f = res;

  CHKERR VecZeroEntries(res);

  CHKERR mField.getInterface<FieldBlas>()->fieldScale(-1, "U");
  CHKERR pip->loopFiniteElements();
  CHKERR mField.getInterface<FieldBlas>()->fieldScale(-1, "U");

  CHKERR VecGhostUpdateBegin(res, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecGhostUpdateEnd(res, ADD_VALUES, SCATTER_REVERSE);
  CHKERR VecAssemblyBegin(res);
  CHKERR VecAssemblyEnd(res);

  auto zero_residual_at_constrains = [&]() {
    MoFEMFunctionBegin;
    auto fe_post_proc_ptr = boost::make_shared<FEMethod>();
    auto get_post_proc_hook_rhs = [&]() {
      MoFEMFunctionBegin;
      CHKERR EssentialPreProcReaction<DisplacementCubitBcData>(
          mField, fe_post_proc_ptr, res)();
      CHKERR EssentialPostProcRhs<DisplacementCubitBcData>(
          mField, fe_post_proc_ptr, 0, res)();
      CHKERR EssentialPostProcRhs<MPCsType>(mField, fe_post_proc_ptr, 0, res)();
      MoFEMFunctionReturn(0);
    };
    fe_post_proc_ptr->postProcessHook = get_post_proc_hook_rhs;
    CHKERR DMoFEMPostProcessFiniteElements(dm, fe_post_proc_ptr.get());
    MoFEMFunctionReturn(0);
  };

  CHKERR zero_residual_at_constrains();

  double nrm2;
  CHKERR VecNorm(res, NORM_2, &nrm2);
  MOFEM_LOG_CHANNEL("WORLD");
  MOFEM_LOG_C("WORLD", Sev::inform, "residual = %3.4e\n", nrm2);

  PetscBool test = PETSC_FALSE;
  CHKERR PetscOptionsGetBool(PETSC_NULL, "", "-test", &test, PETSC_NULL);
  if (test == PETSC_TRUE) {

    auto post_proc_residual = [&](auto dm, auto f_res, auto out_name) {
      MoFEMFunctionBegin;
      auto post_proc_fe =
          boost::make_shared<PostProcBrokenMeshInMoab<DomainEle>>(mField);
      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;
      auto u_vec = boost::make_shared<MatrixDouble>();
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateVectorFieldValues<SPACE_DIM>("U", u_vec, f_res));
      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {},

              {{"RES", u_vec}},

              {}, {})

      );

      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      post_proc_fe);
      post_proc_fe->writeFile(out_name);
      MoFEMFunctionReturn(0);
    };

    CHKERR post_proc_residual(simple->getDM(), res, "res.h5m");

    constexpr double eps = 1e-8;
    if (nrm2 > eps)
      SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
              "Residual is not zero");
  }

  MoFEMFunctionReturn(0);
}
//! [Check]

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  // Initialisation of MoFEM/PETSc and MOAB data structures
  const char param_file[] = "param_file.petsc";
  MoFEM::Core::Initialize(&argc, &argv, param_file, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmWorld(), "TIMER"));

  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSync(), "FieldEvaluator"));
  LogManager::setLog("FieldEvaluator");
  MOFEM_LOG_TAG("FieldEvaluator", "field_eval");

  try {

    //! [Register MoFEM discrete manager in PETSc]
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);
    DMType dm_name_mg = "DMMOFEM_MG";
    CHKERR DMRegister_MGViaApproxOrders(dm_name_mg);
    //! [Register MoFEM discrete manager in PETSc

    //! [Create MoAB]
    moab::Core mb_instance;              ///< mesh database
    moab::Interface &moab = mb_instance; ///< mesh database interface
    //! [Create MoAB]

    //! [Create MoFEM]
    MoFEM::Core core(moab);           ///< finite element database
    MoFEM::Interface &m_field = core; ///< finite element database interface
    //! [Create MoFEM]

    //! [Example]
    Example ex(m_field);
    CHKERR ex.runProblem();
    //! [Example]
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}

struct SetUpSchurImpl : public SetUpSchur {

  SetUpSchurImpl(MoFEM::Interface &m_field) : SetUpSchur(), mField(m_field) {
    if (S) {
      CHK_THROW_MESSAGE(
          MOFEM_DATA_INCONSISTENCY,
          "Is expected that schur matrix is not allocated. This is "
          "possible only is PC is set up twice");
    }
  }
  virtual ~SetUpSchurImpl() = default;

  MoFEMErrorCode setUp(SmartPetscObj<KSP> solver);

private:
  MoFEMErrorCode setEntities();
  MoFEMErrorCode createSubDM();
  MoFEMErrorCode setOperator();
  MoFEMErrorCode setPC(PC pc);
  MoFEMErrorCode setDiagonalPC(PC pc);

  SmartPetscObj<Mat> S;

  MoFEM::Interface &mField;

  SmartPetscObj<DM> schurDM;
  SmartPetscObj<DM> blockDM;
  Range volEnts;
  Range subEnts;
};

MoFEMErrorCode SetUpSchurImpl::setUp(SmartPetscObj<KSP> solver) {
  MoFEMFunctionBegin;
  auto pip = mField.getInterface<PipelineManager>();
  PC pc;
  CHKERR KSPGetPC(solver, &pc);
  PetscBool is_pcfs = PETSC_FALSE;
  PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT, &is_pcfs);
  if (is_pcfs) {
    if (S) {
      CHK_THROW_MESSAGE(
          MOFEM_DATA_INCONSISTENCY,
          "Is expected that schur matrix is not allocated. This is "
          "possible only is PC is set up twice");
    }
    CHKERR setEntities();
    CHKERR createSubDM();

    // Add data to DM storage
    S = createDMMatrix(schurDM);
    CHKERR MatSetDM(S, PETSC_NULL);
    CHKERR MatSetBlockSize(S, SPACE_DIM);
    CHKERR MatSetOption(S, MAT_SYMMETRIC, PETSC_TRUE);

    CHKERR setOperator();
    CHKERR setPC(pc);

    if constexpr (A == AssemblyType::BLOCK_SCHUR) {
      // Set DM to use shell block matrix
      DM solver_dm;
      CHKERR KSPGetDM(solver, &solver_dm);
      CHKERR DMSetMatType(solver_dm, MATSHELL);
    }
    CHKERR KSPSetUp(solver);
    CHKERR setDiagonalPC(pc);

  } else {
    pip->getOpBoundaryLhsPipeline().push_front(createOpSchurAssembleBegin());
    pip->getOpBoundaryLhsPipeline().push_back(createOpSchurAssembleEnd({}, {}));
    pip->getOpDomainLhsPipeline().push_front(createOpSchurAssembleBegin());
    pip->getOpDomainLhsPipeline().push_back(createOpSchurAssembleEnd({}, {}));
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setEntities() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  CHKERR mField.get_moab().get_entities_by_dimension(simple->getMeshset(),
                                                     SPACE_DIM, volEnts);
  CHKERR mField.get_moab().get_entities_by_handle(simple->getMeshset(),
                                                  subEnts);
  subEnts = subtract(subEnts, volEnts);
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode SetUpSchurImpl::createSubDM() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();

  auto create_dm = [&](const char *name, auto &ents, auto dm_type) {
    auto dm = createDM(mField.get_comm(), dm_type);
    auto create_dm_imp = [&]() {
      MoFEMFunctionBegin;
      CHKERR DMMoFEMCreateSubDM(dm, simple->getDM(), name);
      CHKERR DMMoFEMSetSquareProblem(dm, PETSC_TRUE);
      CHKERR DMMoFEMAddElement(dm, simple->getDomainFEName());
      auto sub_ents_ptr = boost::make_shared<Range>(ents);
      CHKERR DMMoFEMAddSubFieldRow(dm, "U", sub_ents_ptr);
      CHKERR DMMoFEMAddSubFieldCol(dm, "U", sub_ents_ptr);
      CHKERR DMSetUp(dm);
      MoFEMFunctionReturn(0);
    };
    CHK_THROW_MESSAGE(create_dm_imp(),
                      "Error in creating schurDM. It is possible that schurDM is "
                      "already created");
    return dm;
  };

  schurDM = create_dm("SCHUR", subEnts, "DMMOFEM_MG");
  blockDM = create_dm("BLOCK", volEnts, "DMMOFEM");

  if constexpr (A == AssemblyType::BLOCK_SCHUR) {

    auto get_nested_mat_data = [&]() -> boost::shared_ptr<NestSchurData> {
      auto block_mat_data =
          createBlockMatStructure(simple->getDM(),

                                  {{

                                      simple->getDomainFEName(),

                                      {{"U", "U"}

                                      }}}

          );

      return createSchurNestedMatrixStruture(

          {schurDM, blockDM}, block_mat_data,

          {"U"}, {boost::make_shared<Range>(volEnts)}

      );
    };

    auto nested_mat_data = get_nested_mat_data();

    CHKERR DMMoFEMSetNestSchurData(simple->getDM(), nested_mat_data);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setOperator() {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  auto pip = mField.getInterface<PipelineManager>();

  // Boundary
  auto dm_is = getDMSubData(schurDM)->getSmartRowIs();
  auto ao_up = createAOMappingIS(dm_is, PETSC_NULL);
  pip->getOpBoundaryLhsPipeline().push_front(createOpSchurAssembleBegin());
  pip->getOpBoundaryLhsPipeline().push_back(createOpSchurAssembleEnd(
      {"U"}, {boost::make_shared<Range>(volEnts)}, ao_up, S, true, true));
  // Domain
  pip->getOpDomainLhsPipeline().push_front(createOpSchurAssembleBegin());
  pip->getOpDomainLhsPipeline().push_back(createOpSchurAssembleEnd(
      {"U"}, {boost::make_shared<Range>(volEnts)}, ao_up, S, true, true));

  auto pre_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();
  auto post_proc_schur_lhs_ptr = boost::make_shared<FEMethod>();

  pre_proc_schur_lhs_ptr->preProcessHook = [this]() {
    MoFEMFunctionBegin;
    if (S) {
      CHKERR MatZeroEntries(S);
    }
    MOFEM_LOG("TIMER", Sev::inform) << "Lhs Assemble Begin";
    MoFEMFunctionReturn(0);
  };

  post_proc_schur_lhs_ptr->postProcessHook = [this, post_proc_schur_lhs_ptr,
                                              ao_up]() {
    MoFEMFunctionBegin;
    CHKERR MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
    CHKERR MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
    CHKERR EssentialPostProcLhs<DisplacementCubitBcData>(
        mField, post_proc_schur_lhs_ptr, 1, S, ao_up)();
    MOFEM_LOG("TIMER", Sev::inform) << "Lhs Assemble End";
    MoFEMFunctionReturn(0);
  };

  auto ksp_ctx_ptr = getDMKspCtx(simple->getDM());
  ksp_ctx_ptr->getPreProcSetOperators().push_front(pre_proc_schur_lhs_ptr);
  ksp_ctx_ptr->getPostProcSetOperators().push_back(post_proc_schur_lhs_ptr);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setPC(PC pc) {
  MoFEMFunctionBegin;
  auto simple = mField.getInterface<Simple>();
  SmartPetscObj<IS> vol_is;
  mField.getInterface<ISManager>()->isCreateProblemFieldAndRank(
      simple->getProblemName(), ROW, "U", 0, SPACE_DIM, vol_is, &volEnts);
  CHKERR PCFieldSplitSetIS(pc, NULL, vol_is);
  CHKERR PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, S);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SetUpSchurImpl::setDiagonalPC(PC pc) {
  MoFEMFunctionBegin;

  auto get_pc = [](auto ksp) {
    PC pc_raw;
    CHKERR KSPGetPC(ksp, &pc_raw);
    return SmartPetscObj<PC>(pc_raw, true); // bump reference
  };

  if constexpr (A == AssemblyType::BLOCK_SCHUR) {
    auto simple = mField.getInterface<Simple>();
    auto A = createDMBlockMat(simple->getDM());
    auto P = createDMNestSchurMat(simple->getDM());
    CHKERR PCSetOperators(pc, A, P);

    KSP *subksp;
    CHKERR PCFieldSplitSchurGetSubKSP(pc, PETSC_NULL, &subksp);
    CHKERR setSchurA00MatSolvePC(get_pc(subksp[0]));

    auto set_pc_p_mg = [](auto dm, auto pc, auto S) {
      MoFEMFunctionBegin;
      CHKERR PCSetDM(pc, dm);
      PetscBool same = PETSC_FALSE;
      PetscObjectTypeCompare((PetscObject)pc, PCMG, &same);
      if (same) {
        CHKERR PCMGSetUpViaApproxOrders(
            pc, createPCMGSetUpViaApproxOrdersCtx(dm, S, true));
        CHKERR PCSetFromOptions(pc);
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR set_pc_p_mg(schurDM, get_pc(subksp[1]), S);

    CHKERR PetscFree(subksp);
  }
  MoFEMFunctionReturn(0);
}

boost::shared_ptr<SetUpSchur>
SetUpSchur::createSetUpSchur(MoFEM::Interface &m_field) {
  return boost::shared_ptr<SetUpSchur>(new SetUpSchurImpl(m_field));
}
