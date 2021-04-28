/**
 * \file scalar_check_approximation_2d.cpp
 * \example scalar_check_approximation_2d.cpp
 *
 * Approximate polynomial in 2D and check direvatives
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

static constexpr int approx_order = 5;
template <int DIM> struct ApproxFunctionsImpl {};

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle2D;
  using DomainEleOp = DomainEle::UserDataOperator;
};

template <> struct ElementsAndOps<3> {
  using DomainEle = VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
};

constexpr int SPACE_DIM = 2; //< Space dimension of problem, mesh

using EntData = DataForcesAndSourcesCore::EntData;
using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = ElementsAndOps<SPACE_DIM>::DomainEleOp;

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
                        DataForcesAndSourcesCore::EntData &data) {
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
      auto t_vals = getFTensor0FromVec(vAls);
      auto t_base_fun = data.getFTensor0N();
      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        auto t_data = data.getFTensor0FieldData();
        for (int bb = 0; bb != nb_dofs; bb++) {
          t_vals += t_base_fun * t_data;
          ++t_base_fun;
          ++t_data;
        }
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
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const double eps = 1e-6;
    const int nb_gauss_pts = data.getN().size1();

    auto t_vals = getFTensor0FromVec(vAls);

    auto t_ptr_vals = getFTensor0FromVec(*ptrVals);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double x = getCoordsAtGaussPts()(gg, 0);
      const double y = getCoordsAtGaussPts()(gg, 1);
      const double z = getCoordsAtGaussPts()(gg, 2);

      // Check approximation
      const double delta_val = t_vals - ApproxFunctions::fUn(x, y, z);

      double err_val = std::fabs(delta_val * delta_val);
      MOFEM_LOG("AT", Sev::verbose) << err_val << " : " << t_vals;
      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID, "Wrong value %4.3e",
                 err_val);

      // Check H1 user data operators
      err_val = std::abs(t_vals - t_ptr_vals);

      if (err_val > eps)
        SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value from operator %4.3e", err_val);

      ++t_vals;
      ++t_ptr_vals;
    }

    if (checkGradients) {

      auto t_diff_vals = getFTensor1FromMat<SPACE_DIM>(diffVals);
      auto t_ptr_diff_vals = getFTensor1FromMat<SPACE_DIM>(*ptrDiffVals);

      for (int gg = 0; gg != nb_gauss_pts; gg++) {
        const double x = getCoordsAtGaussPts()(gg, 0);
        const double y = getCoordsAtGaussPts()(gg, 1);
        const double z = getCoordsAtGaussPts()(gg, 2);

        // Check approximation
        FTensor::Tensor1<double, SPACE_DIM> t_delta_diff_val;
        auto t_diff_anal = ApproxFunctions::diffFun(x, y, z);
        t_delta_diff_val(i) = t_diff_vals(i) - t_diff_anal(i);

        double err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        MOFEM_LOG("AT", Sev::verbose)
            << err_diff_val << " : " << sqrt(t_diff_vals(i) * t_diff_vals(i))
            << " :  " << t_diff_vals(0) / t_diff_anal(0) << " "
            << t_diff_vals(1) / t_diff_anal(1) << "  "
            << t_diff_vals(2) / t_diff_anal(2);
        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong derivative of value %4.3e", err_diff_val);

        t_delta_diff_val(i) = t_diff_vals(i) - t_ptr_diff_vals(i);
        err_diff_val = sqrt(t_delta_diff_val(i) * t_delta_diff_val(i));
        if (err_diff_val > eps)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong direvatives from operator %4.3e", err_diff_val);

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
    core_log->add_sink(LogManager::createSink(LogManager::getStrmSelf(), "AT"));
    LogManager::setLog("AT");
    MOFEM_LOG_TAG("AT", "atom_test");

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    Simple *simple_interface = m_field.getInterface<Simple>();
    PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();
    CHKERR simple_interface->getOptions();
    CHKERR simple_interface->loadFile("", "");

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

    CHKERR simple_interface->addDomainField("FIELD1", space, base, 1);
    CHKERR simple_interface->setFieldOrder("FIELD1", approx_order);
    CHKERR simple_interface->setUp();
    auto dm = simple_interface->getDM();

    VectorDouble vals;
    MatrixDouble jac, inv_jac, diff_vals;

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpMakeHighOrderGeometryWeightsOnFace());

      using OpSource = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSource("FIELD1", ApproxFunctions::fUn));
      // pipeline_mng->getOpDomainRhsPipeline().push_back(new OpAssembleVec());

      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpCalculateInvJacForFace(inv_jac));
      pipeline_mng->getOpDomainLhsPipeline().push_back(
          new OpMakeHighOrderGeometryWeightsOnFace());
      using OpMass = FormsIntegrators<DomainEleOp>::Assembly<
          PETSC>::BiLinearForm<GAUSS>::OpMass<1, 1>;
      pipeline_mng->getOpDomainLhsPipeline().push_back(new OpMass(
          "FIELD1", "FIELD1", [](double, double, double) { return 1.; }));

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

      auto dm = simple_interface->getDM();
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

      boost::shared_ptr<VectorDouble> ptr_values(new VectorDouble());
      boost::shared_ptr<MatrixDouble> ptr_diff_vals(new MatrixDouble());

      pipeline_mng->getDomainLhsFE().reset();
      pipeline_mng->getOpDomainRhsPipeline().clear();

      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateInvJacForFace(inv_jac));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpSetInvJacH1ForFace(inv_jac));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpValsDiffVals(vals, diff_vals, choice_space_value == H1SPACE));
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCalculateScalarFieldValues("FIELD1", ptr_values));
      if (choice_space_value == H1SPACE) {
        pipeline_mng->getOpDomainRhsPipeline().push_back(
            new OpCalculateScalarFieldGradient<SPACE_DIM>("FIELD1",
                                                          ptr_diff_vals));
      }
      pipeline_mng->getOpDomainRhsPipeline().push_back(
          new OpCheckValsDiffVals(vals, diff_vals, ptr_values, ptr_diff_vals,
                                  choice_space_value == H1SPACE));

      CHKERR pipeline_mng->loopFiniteElements();

      MoFEMFunctionReturn(0);
    };

    CHKERR assemble_matrices_and_vectors();
    CHKERR solve_problem();
    CHKERR check_solution();
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
