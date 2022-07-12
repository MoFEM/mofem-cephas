/** \file prism_polynomial_approximation.cpp
  \example prism_polynomial_approximation.cpp
  \brief Checking approximation functions on prism
*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

static constexpr int approx_order = 6;

struct ApproxFunction {
  static inline double fun(double x, double y, double z) {
    double r = 1;
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        for (int j = 0; j <= (o - i); ++j) {
          int k = o - i - j;
          if (k >= 0) {
            r += pow(x, i) * pow(y, j) * pow(z, k);
          }
        }
      }
    }
    return r;
  }

  static inline VectorDouble3 diff_fun(double x, double y, double z) {
    VectorDouble3 r(3);
    r.clear();
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        for (int j = 0; j <= (o - i); ++j) {
          int k = o - i - j;
          if (k >= 0) {
            r[0] += i > 0 ? i * pow(x, i - 1) * pow(y, j) * pow(z, k) : 0;
            r[1] += j > 0 ? j * pow(x, i) * pow(y, j - 1) * pow(z, k) : 0;
            r[2] += k > 0 ? k * pow(x, i) * pow(y, j) * pow(z, k - 1) : 0;
          }
        }
      }
    }
    return r;
  }
};

struct PrismOpCheck
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  PrismOpCheck(boost::shared_ptr<VectorDouble> &field_vals,
               boost::shared_ptr<MatrixDouble> &diff_field_vals);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<VectorDouble> fieldVals;
  boost::shared_ptr<MatrixDouble> diffFieldVals;
};

struct PrismOpRhs
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  PrismOpRhs(SmartPetscObj<Vec> &f);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  SmartPetscObj<Vec> F;
};

struct PrismOpLhs
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  PrismOpLhs(SmartPetscObj<Mat> &a);
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data);

private:
  SmartPetscObj<Mat> A;
};

struct PrismFE : public FatPrismElementForcesAndSourcesCore {

  using FatPrismElementForcesAndSourcesCore::
      FatPrismElementForcesAndSourcesCore;
  int getRuleTrianglesOnly(int order);
  int getRuleThroughThickness(int order);
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    auto moab_comm_wrap =
        boost::make_shared<WrapMPIComm>(PETSC_COMM_WORLD, false);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, moab_comm_wrap->get_comm());

    std::array<double, 18> one_prism_coords = {0, 0, 0, 1, 0, 0, 0, 1, 0,
                                               0, 0, 1, 1, 0, 1, 0, 1, 1};
    std::array<EntityHandle, 6> one_prism_nodes;
    for (int n = 0; n != 6; ++n)
      CHKERR moab.create_vertex(&one_prism_coords[3 * n], one_prism_nodes[n]);
    EntityHandle one_prism;
    CHKERR moab.create_element(MBPRISM, one_prism_nodes.data(), 6, one_prism);
    Range one_prism_range;
    one_prism_range.insert(one_prism);
    Range one_prism_adj_ents;
    for (int d = 1; d != 3; ++d)
      CHKERR moab.get_adjacencies(one_prism_range, d, true, one_prism_adj_ents,
                                  moab::Interface::UNION);

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    PrismsFromSurfaceInterface *prisms_from_surface_interface;
    CHKERR m_field.getInterface(prisms_from_surface_interface);
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setEntitiesBitRefLevel(
        one_prism_range, bit_level0);
    CHKERR prisms_from_surface_interface->seedPrismsEntities(one_prism_range,
                                                             bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_ents_to_field_by_type(0, MBPRISM, "FIELD1");

    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD1", approx_order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD1", approx_order);
    CHKERR m_field.set_field_order(0, MBQUAD, "FIELD1", approx_order);
    CHKERR m_field.set_field_order(0, MBPRISM, "FIELD1", approx_order);
    CHKERR m_field.build_fields();

    // FE
    CHKERR m_field.add_finite_element("PRISM");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("PRISM", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("PRISM", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("PRISM", "FIELD1");
    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBPRISM, "PRISM");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "PRISM");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    // Create matrices
    SmartPetscObj<Mat> A;
    CHKERR m_field.getInterface<MatrixManager>()
        ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>("TEST_PROBLEM", A);
    SmartPetscObj<Vec> F;
    CHKERR m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",
                                                              ROW, F);
    SmartPetscObj<Vec> D;
    CHKERR m_field.getInterface<VecManager>()->vecCreateGhost("TEST_PROBLEM",
                                                              COL, D);

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;
      PrismFE fe(m_field);
      MatrixDouble inv_jac;
      fe.getOpPtrVector().push_back(
          new OpMultiplyDeterminantOfJacobianAndWeightsForFatPrisms());
      fe.getOpPtrVector().push_back(
          new MoFEM::OpCalculateInvJacForFatPrism(inv_jac));
      fe.getOpPtrVector().push_back(
          new MoFEM::OpSetInvJacH1ForFatPrism(inv_jac));
      fe.getOpPtrVector().push_back(new PrismOpRhs(F));
      fe.getOpPtrVector().push_back(new PrismOpLhs(A));
      CHKERR VecZeroEntries(F);
      CHKERR MatZeroEntries(A);
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "PRISM", fe);
      CHKERR VecAssemblyBegin(F);
      CHKERR VecAssemblyEnd(F);
      CHKERR MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      CHKERR MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
      MoFEMFunctionReturn(0);
    };

    auto solve_problem = [&] {
      MoFEMFunctionBegin;
      auto solver = createKSP(PETSC_COMM_WORLD);
      CHKERR KSPSetOperators(solver, A, A);
      CHKERR KSPSetFromOptions(solver);
      CHKERR KSPSetUp(solver);
      CHKERR KSPSolve(solver, F, D);
      CHKERR VecGhostUpdateBegin(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR VecGhostUpdateEnd(D, INSERT_VALUES, SCATTER_FORWARD);
      CHKERR m_field.getInterface<VecManager>()->setLocalGhostVector(
          "TEST_PROBLEM", COL, D, INSERT_VALUES, SCATTER_REVERSE);
      MoFEMFunctionReturn(0);
    };

    auto check_solution = [&] {
      MoFEMFunctionBegin;
      PrismFE fe(m_field);
      boost::shared_ptr<VectorDouble> field_vals_ptr(new VectorDouble());
      boost::shared_ptr<MatrixDouble> diff_field_vals_ptr(new MatrixDouble());
      MatrixDouble inv_jac;
      fe.getOpPtrVector().push_back(
          new OpCalculateInvJacForFatPrism(inv_jac));
      fe.getOpPtrVector().push_back(
          new OpSetInvJacH1ForFatPrism(inv_jac));
      fe.getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("FIELD1", field_vals_ptr));
      fe.getOpPtrVector().push_back(
          new OpCalculateScalarFieldGradient<3>("FIELD1", diff_field_vals_ptr));
      fe.getOpPtrVector().push_back(
          new PrismOpCheck(field_vals_ptr, diff_field_vals_ptr));
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "PRISM", fe);
      MoFEMFunctionReturn(0);
    };

    CHKERR assemble_matrices_and_vectors();
    CHKERR solve_problem();
    CHKERR check_solution();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}

PrismOpCheck::PrismOpCheck(boost::shared_ptr<VectorDouble> &field_vals,
                           boost::shared_ptr<MatrixDouble> &diff_field_vals)
    : FatPrismElementForcesAndSourcesCore::UserDataOperator(
          "FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      fieldVals(field_vals), diffFieldVals(diff_field_vals) {}

MoFEMErrorCode PrismOpCheck::doWork(int side, EntityType type,
                                    EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  if (type == MBVERTEX) {
    const int nb_gauss_pts = data.getN().size2();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      double f = ApproxFunction::fun(t_coords(0), t_coords(1), t_coords(2));
      VectorDouble3 diff_f =
          ApproxFunction::diff_fun(t_coords(0), t_coords(1), t_coords(2));

      std::cout << f - (*fieldVals)[gg] << " : ";
      for (auto d : {0, 1, 2}) 
        std::cout << diff_f[d] - (*diffFieldVals)(d, gg) << " ";
      std::cout << std::endl;
          
      constexpr double eps = 1e-6;
      if (std::abs(f - (*fieldVals)[gg]) > eps ||
          !std::isnormal((*fieldVals)[gg]))
        SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value %6.4e != %6.4e (%6.4e)", f, (*fieldVals)[gg],
                 f - (*fieldVals)[gg]);
      
      for (auto d : {0, 1, 2})
        if (std::abs(diff_f[d] - (*diffFieldVals)(d, gg)) > eps ||
            !std::isnormal((*diffFieldVals)(d, gg)))
          SETERRQ3(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong diff value %6.4e != %6.4e (%6.4e)", diff_f[d],
                   (*diffFieldVals)(d, gg),
                   diff_f[d] - (*diffFieldVals)(d, gg));

      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
}

PrismOpRhs::PrismOpRhs(SmartPetscObj<Vec> &f)
    : FatPrismElementForcesAndSourcesCore::UserDataOperator(
          "FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      F(f) {}

MoFEMErrorCode PrismOpRhs::doWork(int side, EntityType type,
                                  EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_dofs = data.getN().size2();
  if (nb_dofs) {
    const int nb_gauss_pts = data.getN().size1();
    VectorDouble nf(nb_dofs);
    nf.clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    auto t_w = getFTensor0IntegrationWeight();
    double vol = getVolume();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      double f = ApproxFunction::fun(t_coords(0), t_coords(1), t_coords(2));
      double v = t_w * vol * f;
      double *val = &*nf.begin();
      for (int bb = 0; bb != nb_dofs; ++bb) {
        *val += v * t_base;
        ++t_base;
        ++val;
      }
      ++t_coords;
      ++t_w;
    }
    CHKERR VecSetValues(F, data, &*nf.data().begin(), ADD_VALUES);
  }
  MoFEMFunctionReturn(0);
}

PrismOpLhs::PrismOpLhs(SmartPetscObj<Mat> &a)
    : FatPrismElementForcesAndSourcesCore::UserDataOperator(
          "FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROWCOL),
      A(a) {
  // FIXME: Can be symmetric, is not for simplicity
  sYmm = false;
}

MoFEMErrorCode PrismOpLhs::doWork(int row_side, int col_side,
                                  EntityType row_type, EntityType col_type,
                                  EntitiesFieldData::EntData &row_data,
                                  EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;
  const int row_nb_dofs = row_data.getN().size2();
  const int col_nb_dofs = col_data.getN().size2();
  if (row_nb_dofs && col_nb_dofs) {
    const int nb_gauss_pts = row_data.getN().size1();
    MatrixDouble m(row_nb_dofs, col_nb_dofs);
    m.clear();
    auto t_w = getFTensor0IntegrationWeight();
    double vol = getVolume();
    double *row_base_ptr = &*row_data.getN().data().begin();
    double *col_base_ptr = &*col_data.getN().data().begin();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {

      double v = t_w * vol;
      cblas_dger(CblasRowMajor, row_nb_dofs, col_nb_dofs, v, row_base_ptr, 1,
                 col_base_ptr, 1, &*m.data().begin(), col_nb_dofs);

      row_base_ptr += row_nb_dofs;
      col_base_ptr += col_nb_dofs;
      ++t_w;
    }
    CHKERR MatSetValues(A, row_data, col_data, &*m.data().begin(), ADD_VALUES);
  }
  MoFEMFunctionReturn(0);
}

int PrismFE::getRuleTrianglesOnly(int order) { return 2 * (order + 1); };
int PrismFE::getRuleThroughThickness(int order) { return 2 * (order + 1); };
