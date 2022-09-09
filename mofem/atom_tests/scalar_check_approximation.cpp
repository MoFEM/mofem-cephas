/**
 * \file scalar_check_approximation.cpp
 * \example scalar_check_approximation.cpp
 *
 * Approximate polynomial in 2D and check derivatives
 *
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

static int approx_order = 4;

template <int DIM> struct ApproxFunctionsImpl {};

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using DomainEleOp = DomainEle::UserDataOperator;
  using BoundaryEle = PipelineManager::EdgeEle;
  using EleOnSide = FaceElementForcesAndSourcesCoreOnSide;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
  using BoundaryEle = PipelineManager::FaceEle;
  using EleOnSide = VolumeElementForcesAndSourcesCoreOnSide;
};

constexpr int SPACE_DIM =
    EXECUTABLE_DIMENSION; //< Space dimension of problem, mesh

using EntData = EntitiesFieldData::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = ElementsAndOps<SPACE_DIM>::DomainEleOp;
using BoundaryEle = ElementsAndOps<SPACE_DIM>::BoundaryEle;
using EleOnSide = ElementsAndOps<SPACE_DIM>::EleOnSide;

template <> struct ApproxFunctionsImpl<2> {
  static double fUn(const double x, const double y, double z) {
    double r = 1;
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        int j = o - i;
        if (j >= 0)
          r += pow(x, i) * pow(y, j);
      }
    }
    return r;
  }

  static FTensor::Tensor1<double, 2> diffFun(const double x, const double y,
                                             double z) {
    FTensor::Tensor1<double, 2> r{0., 0.};
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        int j = o - i;
        if (j >= 0) {
          r(0) += i > 0 ? i * pow(x, i - 1) * pow(y, j) : 0;
          r(1) += j > 0 ? j * pow(x, i) * pow(y, j - 1) : 0;
        }
      }
    }
    return r;
  }
};

template <> struct ApproxFunctionsImpl<3> {
  static double fUn(const double x, const double y, double z) {
    double r = 1;
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        for (int j = 0; j <= o - i; j++) {
          int k = o - i - j;
          if (k >= 0) {
            r += pow(x, i) * pow(y, j) * pow(z, k);
          }
        }
      }
    }
    return r;
  }

  static FTensor::Tensor1<double, 3> diffFun(const double x, const double y,
                                             double z) {
    FTensor::Tensor1<double, 3> r{0., 0., 0.};
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        for (int j = 0; j <= o - i; j++) {
          int k = o - i - j;
          if (k >= 0) {
            r(0) += i > 0 ? i * pow(x, i - 1) * pow(y, j) * pow(z, k) : 0;
            r(1) += j > 0 ? j * pow(x, i) * pow(y, j - 1) * pow(z, k) : 0;
            r(2) += k > 0 ? k * pow(x, i) * pow(y, j) * pow(z, k - 1) : 0;
          }
        }
      }
    }
    return r;
  }
};

using ApproxFunctions = ApproxFunctionsImpl<SPACE_DIM>;

struct OpValsDiffVals : public DomainEleOp {
  VectorDouble &vAls;
  MatrixDouble &diffVals;
  const bool checkGradients;
  OpValsDiffVals(VectorDouble &vals, MatrixDouble &diff_vals, bool check_grads)
      : DomainEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals),
        checkGradients(check_grads) {}

  FTensor::Index<'i', SPACE_DIM> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_gauss_pts = getGaussPts().size2();
    if (type == MBVERTEX) {
      vAls.resize(nb_gauss_pts, false);
      diffVals.resize(SPACE_DIM, nb_gauss_pts, false);
      vAls.clear();
      diffVals.clear();
    }

    const int nb_dofs = data.getIndices().size();
    if (nb_dofs) {

      MOFEM_LOG("AT", Sev::noisy)
          << "Type  " << moab::CN::EntityTypeName(type) << " side " << side;
      MOFEM_LOG("AT", Sev::noisy) << data.getN();
      MOFEM_LOG("AT", Sev::noisy) << data.getDiffN();

      auto t_vals = getFTensor0FromVec(vAls);
      auto t_base_fun = data.getFTensor0N();
      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        auto t_data = data.getFTensor0FieldData();
        for (int bb = 0; bb != nb_dofs; bb++) {
          t_vals += t_base_fun * t_data;
          ++t_base_fun;
          ++t_data;
        }
        const double v = t_vals;
        if (!std::isnormal(v))
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Not a number");
        ++t_vals;
      }

      if (checkGradients) {
        auto t_diff_vals = getFTensor1FromMat<SPACE_DIM>(diffVals);
        auto t_diff_base_fun = data.getFTensor1DiffN<SPACE_DIM>();
        for (int gg = 0; gg != nb_gauss_pts; gg++) {
          auto t_data = data.getFTensor0FieldData();
          for (int bb = 0; bb != nb_dofs; bb++) {
            t_diff_vals(i) += t_diff_base_fun(i) * t_data;
            ++t_diff_base_fun;
            ++t_data;
          }
          for (int d = 0; d != SPACE_DIM; ++d)
            if (!std::isnormal(t_diff_vals(d)))
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Not a number");
          ++t_diff_vals;
        }
      }
    }
    MoFEMFunctionReturn(0);
  }
};

struct OpCheckValsDiffVals : public DomainEleOp {
  VectorDouble &vAls;
  MatrixDouble &diffVals;
  boost::shared_ptr<VectorDouble> ptrVals;
  boost::shared_ptr<MatrixDouble> ptrDiffVals;
  const bool checkGradients;

  OpCheckValsDiffVals(VectorDouble &vals, MatrixDouble &diff_vals,
                      boost::shared_ptr<VectorDouble> &ptr_vals,
                      boost::shared_ptr<MatrixDouble> &ptr_diff_vals,
                      bool check_grads)
      : DomainEleOp("FIELD1", OPROW), vAls(vals), diffVals(diff_vals),
        ptrVals(ptr_vals), ptrDiffVals(ptr_diff_vals),
        checkGradients(check_grads) {
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
  }

  FTensor::Index<'i', SPACE_DIM> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    const double eps = 1e-6;
    const int nb_gauss_pts = data.getN().size1();

    auto t_vals = getFTensor0FromVec(vAls);

    auto t_ptr_vals = getFTensor0FromVec(*ptrVals);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {

      double err_val;

      // Check user data operators
      err_val = std::abs(t_vals - t_ptr_vals);
      MOFEM_LOG("AT", Sev::noisy) << "Val op error " << err_val;

      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value from operator %4.3e", err_val);

      const double x = getCoordsAtGaussPts()(gg, 0);
      const double y = getCoordsAtGaussPts()(gg, 1);
      const double z = getCoordsAtGaussPts()(gg, 2);

      // Check approximation
      const double delta_val = t_vals - ApproxFunctions::fUn(x, y, z);

      err_val = std::fabs(delta_val * delta_val);
      MOFEM_LOG("AT", Sev::verbose) << err_val << " : " << t_vals;
      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value %4.3e",
                 err_val);

      ++t_vals;
      ++t_ptr_vals;
    }

    if (checkGradients) {

      auto t_diff_vals = getFTensor1FromMat<SPACE_DIM>(diffVals);
      auto t_ptr_diff_vals = getFTensor1FromMat<SPACE_DIM>(*ptrDiffVals);

      for (int gg = 0; gg != nb_gauss_pts; gg++) {

        FTensor::Tensor1<double, SPACE_DIM> t_delta_diff_val;
        double err_diff_val;

        t_delta_diff_val(i) = t_diff_vals(i) - t_ptr_diff_vals(i);
        err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        MOFEM_LOG("AT", Sev::noisy) << "Diff val op error " << err_diff_val;

        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivatives from operator %4.3e", err_diff_val);

        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);
        const double z = getCoordsAtGaussPts()(gg, 2);

        // Check approximation
        auto t_diff_anal = ApproxFunctions::diffFun(x, y, z);
        t_delta_diff_val(i) = t_diff_vals(i) - t_diff_anal(i);

        err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        if (SPACE_DIM == 3)
          MOFEM_LOG("AT", Sev::noisy)
              << "Diff val " << err_diff_val << " : "
              << sqrt(t_diff_vals(i) * t_diff_vals(i)) << " :  "
              << t_diff_vals(0) << " (" << t_diff_anal(0) << ") "
              << t_diff_vals(1) << " (" << t_diff_anal(1) << ")  "
              << t_diff_vals(2) << " (" << t_diff_anal(2) << ")";
        else
          MOFEM_LOG("AT", Sev::noisy)
              << "Diff val " << err_diff_val << " : "
              << sqrt(t_diff_vals(i) * t_diff_vals(i)) << " :  "
              << t_diff_vals(0) << " (" << t_diff_anal(0) << ") "
              << t_diff_vals(1) << " (" << t_diff_anal(1) << ")";

        MOFEM_LOG("AT", Sev::verbose)
            << getCoords()(3 * 1 + 0) - getCoords()(3 * 0 + 0);
        MOFEM_LOG("AT", Sev::verbose)
            << getCoords()(3 * 1 + 1) - getCoords()(3 * 0 + 1);
        MOFEM_LOG("AT", Sev::verbose)
            << getCoords()(3 * 1 + 2) - getCoords()(3 * 0 + 2);

        MOFEM_LOG("AT", Sev::verbose) << "Diff val error " << err_diff_val;
        if (err_diff_val > eps)
          SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivative of value %4.3e %4.3e", err_diff_val,
                   t_diff_anal.l2());

        ++t_diff_vals;
        ++t_ptr_diff_vals;
      }
    }

    MoFEMFunctionReturn(0);
  }
};

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

    Simple *simple = m_field.getInterface<Simple>();
    PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();
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
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH_LOBATTO)
      base = AINSWORTH_LOBATTO_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    else if (choice_base_value == BERNSTEIN)
      base = AINSWORTH_BERNSTEIN_BEZIER_BASE;

    enum spaces { H1SPACE, L2SPACE, LASBASETSPACE };
    const char *list_spaces[] = {"h1", "l2"};
    PetscInt choice_space_value = H1SPACE;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                LASBASETSPACE, &choice_space_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "space not set");
    FieldSpace space = H1;
    if (choice_space_value == H1SPACE)
      space = H1;
    else if (choice_space_value == L2SPACE)
      space = L2;

    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &approx_order,
                              PETSC_NULL);

    CHKERR simple->addDomainField("FIELD1", space, base, 1);
    CHKERR simple->setFieldOrder("FIELD1", approx_order);

    Range edges, faces;
    CHKERR moab.get_entities_by_dimension(0, 1, edges);
    CHKERR moab.get_entities_by_dimension(0, 2, faces);

    if (choice_base_value != BERNSTEIN) {
      Range rise_order;

      int i = 0;
      for (auto e : edges) {
        if (!(i % 2)) {
          rise_order.insert(e);
        }
        ++i;
      }

      for (auto f : faces) {
        if (!(i % 3)) {
          rise_order.insert(f);
        }
        ++i;
      }

      Range rise_twice;
      for (auto e : rise_order) {
        if (!(i % 2)) {
          rise_twice.insert(e);
        }
        ++i;
      }

      CHKERR simple->setFieldOrder("FIELD1", approx_order + 1, &rise_order);

      CHKERR simple->setFieldOrder("FIELD1", approx_order + 2, &rise_twice);
    }

    CHKERR simple->defineFiniteElements();

    auto volume_adj = [](moab::Interface &moab, const Field &field,
                         const EntFiniteElement &fe,
                         std::vector<EntityHandle> &adjacency) {
      MoFEMFunctionBegin;
      EntityHandle fe_ent = fe.getEnt();
      switch (field.getSpace()) {
      case H1:
        CHKERR moab.get_connectivity(&fe_ent, 1, adjacency, true);
      case HCURL:
        CHKERR moab.get_adjacencies(&fe_ent, 1, 1, false, adjacency,
                                    moab::Interface::UNION);
      case HDIV:
      case L2:
        CHKERR moab.get_adjacencies(&fe_ent, 1, 2, false, adjacency,
                                    moab::Interface::UNION);
        adjacency.push_back(fe_ent);
        // build side table
        for (auto ent : adjacency)
          fe.getSideNumberPtr(ent);
        break;
      case NOFIELD: {
        CHKERR moab.get_entities_by_handle(field.getMeshset(), adjacency,
                                           false);
        for (auto ent : adjacency) {
          const_cast<SideNumber_multiIndex &>(fe.getSideNumberTable())
              .insert(
                  boost::shared_ptr<SideNumber>(new SideNumber(ent, -1, 0, 0)));
        }
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "this field is not implemented for TRI finite element");
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR m_field.modify_finite_element_adjacency_table(
        simple->getDomainFEName(), MBTET, volume_adj);
    CHKERR m_field.modify_finite_element_adjacency_table(
        simple->getDomainFEName(), MBHEX, volume_adj);

    CHKERR simple->defineProblem(PETSC_TRUE);
    CHKERR simple->buildFields();
    CHKERR simple->buildFiniteElements();
    CHKERR simple->buildProblem();

    auto dm = simple->getDM();

    VectorDouble vals;
    auto jac_ptr = boost::make_shared<MatrixDouble>();
    auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
    auto det_ptr = boost::make_shared<VectorDouble>();
    MatrixDouble diff_vals;

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;

      using OpSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;

      if (SPACE_DIM == 2) {
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpSetHOWeightsOnFace());
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpSetHOWeightsOnFace());
      }

      if (SPACE_DIM == 3) {
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, nullptr));
        pipeline_mng->getOpDomainLhsPipeline().push_back(
            new OpSetHOWeights(det_ptr));
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, nullptr));
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpSetHOWeights(det_ptr));
      }

      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<GAUSS>::OpMass<1, 1>;
      pipeline_mng->getOpDomainLhsPipeline().push_back(new OpMass(
          "FIELD1", "FIELD1", [](double, double, double) { return 1.; }));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSource("FIELD1", ApproxFunctions::fUn));

      auto integration_rule = [](int, int, int p_data) {
        return 2 * p_data + 1;
      };
      CHKERR pipeline_mng->setDomainRhsIntegrationRule(integration_rule);
      CHKERR pipeline_mng->setDomainLhsIntegrationRule(integration_rule);

      MoFEMFunctionReturn(0);
    };

    auto solve_problem = [&] {
      MoFEMFunctionBegin;
      auto solver = pipeline_mng->createKSP();
      CHKERR KSPSetFromOptions(solver);
      CHKERR KSPSetUp(solver);

      auto dm = simple->getDM();
      auto D = smartCreateDMVector(dm);
      auto F = smartVectorDuplicate(D);

      CHKERR KSPSolve(solver, F, D);
      CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR DMoFEMMeshToLocalVector(dm, D, INSERT_VALUES, SCATTER_REVERSE);

      MoFEMFunctionReturn(0);
    };

    auto check_solution = [&] {
      MoFEMFunctionBegin;

      auto ptr_values = boost::make_shared<VectorDouble>();
      auto ptr_diff_vals = boost::make_shared<MatrixDouble>();

      pipeline_mng->getDomainLhsFE().reset();
      pipeline_mng->getOpDomainRhsPipeline().clear();

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetHOInvJacToScalarBases<SPACE_DIM>(space, inv_jac_ptr));

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpValsDiffVals(vals, diff_vals, true));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldValues("FIELD1", ptr_values));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>("FIELD1",
                                                        ptr_diff_vals));
      pipeline_mng->getOpDomainRhsPipeline().push_back(new OpCheckValsDiffVals(
          vals, diff_vals, ptr_values, ptr_diff_vals, true));

      CHKERR pipeline_mng->loopFiniteElements();

      MoFEMFunctionReturn(0);
    };

    auto post_proc = [&] {
      MoFEMFunctionBegin;

      auto post_proc_fe =
          boost::make_shared<PostProcBrokenMeshInMoab<DomainEle>>(m_field);

      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      post_proc_fe->getOpPtrVector().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      post_proc_fe->getOpPtrVector().push_back(
          new OpSetHOInvJacToScalarBases<SPACE_DIM>(space, inv_jac_ptr));

      auto ptr_values = boost::make_shared<VectorDouble>();
      auto ptr_diff_vals = boost::make_shared<MatrixDouble>();

      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("FIELD1", ptr_values));
      post_proc_fe->getOpPtrVector().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>("FIELD1",
                                                        ptr_diff_vals));

      using OpPPMap = OpPostProcMapInMoab<SPACE_DIM, SPACE_DIM>;

      post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              post_proc_fe->getPostProcMesh(), post_proc_fe->getMapGaussPts(),

              {{"FIELD1", ptr_values}},

              {{"FIELD1_GRAD", ptr_diff_vals}},

              {}, {})

      );

      CHKERR DMoFEMLoopFiniteElements(dm, simple->getDomainFEName(),
                                      post_proc_fe);

      post_proc_fe->writeFile("out_post_proc.h5m");

      auto bdy_post_proc_fe =
          make_post_proc_fe_in_moab<PostProcBrokenMeshInMoab<BoundaryEle>>(
              m_field);

      auto op_loop_side = new OpLoopSide<EleOnSide>(
          m_field, simple->getDomainFEName(), SPACE_DIM);

      // push operators to side element
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateHOJac<SPACE_DIM>(jac_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpInvertMatrix<SPACE_DIM>(jac_ptr, det_ptr, inv_jac_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpSetHOInvJacToScalarBases<SPACE_DIM>(space, inv_jac_ptr));
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("FIELD1", ptr_values));
      op_loop_side->getOpPtrVector().push_back(
          new OpCalculateScalarFieldGradient<SPACE_DIM>("FIELD1",
                                                        ptr_diff_vals));
      // push op to boundary element
      bdy_post_proc_fe->getOpPtrVector().push_back(op_loop_side);

      bdy_post_proc_fe->getOpPtrVector().push_back(

          new OpPPMap(

              bdy_post_proc_fe->getPostProcMesh(),
              bdy_post_proc_fe->getMapGaussPts(),

              {{"FIELD1", ptr_values}},

              {{"FIELD1_GRAD", ptr_diff_vals}},

              {}, {})

      );

      CHKERR DMoFEMLoopFiniteElements(dm, simple->getBoundaryFEName(),
                                      bdy_post_proc_fe);

      bdy_post_proc_fe->writeFile("out_post_proc_bdy.h5m");

      MoFEMFunctionReturn(0);
    };

    CHKERR assemble_matrices_and_vectors();
    CHKERR solve_problem();
    CHKERR check_solution();
    CHKERR post_proc();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
