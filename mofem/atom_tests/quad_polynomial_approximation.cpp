/** \file quad_polynomial_approximation.cpp
  \example quad_polynomial_approximation.cpp
  \brief Checking approximation functions for quad

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

using Ele = FaceElementForcesAndSourcesCore;
using OpEle = FaceElementForcesAndSourcesCore::UserDataOperator;
using EntData = DataForcesAndSourcesCore::EntData;

static char help[] = "...\n\n";
static int debug = 1;

static constexpr int approx_order = 5;
struct ApproxFunction {
  static inline double fun(double x, double y) {
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

  static inline VectorDouble3 diff_fun(double x, double y) {
    VectorDouble3 r(2);
    r.clear();
    for (int o = 1; o <= approx_order; ++o) {
      for (int i = 0; i <= o; ++i) {
        int j = o - i;
        if (j >= 0) {
          r[0] += i > 0 ? i * pow(x, i - 1) * pow(y, j) : 0;
          r[1] += j > 0 ? j * pow(x, i) * pow(y, j - 1) : 0;
        }
      }
    }
    return r;
  }
};

struct QuadOpCheck : public OpEle {

  QuadOpCheck(boost::shared_ptr<VectorDouble> &field_vals,
              boost::shared_ptr<MatrixDouble> &diff_field_vals);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<VectorDouble> fieldVals;
  boost::shared_ptr<MatrixDouble> diffFieldVals;
};

struct QuadOpRhs : public OpEle {

  QuadOpRhs(SmartPetscObj<Vec> &f);
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  SmartPetscObj<Vec> F;
};

struct QuadOpLhs : public OpEle {

  QuadOpLhs(SmartPetscObj<Mat> &a);
  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        DataForcesAndSourcesCore::EntData &row_data,
                        DataForcesAndSourcesCore::EntData &col_data);

private:
  SmartPetscObj<Mat> A;
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    std::array<double, 12> one_quad_coords = {0, 0, 0,

                                              1, 0, 0,

                                              1, 1, 0,

                                              0, 1, 0};
    std::array<EntityHandle, 4> one_quad_nodes;
    for (int n = 0; n != 4; ++n)
      CHKERR moab.create_vertex(&one_quad_coords[3 * n], one_quad_nodes[n]);
    EntityHandle one_quad;
    CHKERR moab.create_element(MBQUAD, one_quad_nodes.data(), 4, one_quad);
    Range one_quad_range;
    one_quad_range.insert(one_quad);
    Range one_quad_adj_ents;
    CHKERR moab.get_adjacencies(one_quad_range, 1, true, one_quad_adj_ents,
                                moab::Interface::UNION);

    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    BitRefLevel bit_level0 = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_ents_to_field_by_type(0, MBQUAD, "FIELD1");

    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD1", 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD1", approx_order);
    CHKERR m_field.set_field_order(0, MBQUAD, "FIELD1", approx_order);
    CHKERR m_field.build_fields();

    // FE
    CHKERR m_field.add_finite_element("QUAD");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("QUAD", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("QUAD", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("QUAD", "FIELD1");
    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBQUAD, "QUAD");

    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // //build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "QUAD");
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

    auto rule = [&](int, int, int p) { return 2 * (p + 1); };

    auto assemble_matrices_and_vectors = [&]() {
      MoFEMFunctionBegin;
      Ele fe(m_field);
      fe.getRuleHook = rule;
      MatrixDouble inv_jac;
      fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace(inv_jac));
      fe.getOpPtrVector().push_back(new OpSetInvJacH1ForFace(inv_jac));
      fe.getOpPtrVector().push_back(new QuadOpRhs(F));
      fe.getOpPtrVector().push_back(new QuadOpLhs(A));
      CHKERR VecZeroEntries(F);
      CHKERR MatZeroEntries(A);
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "QUAD", fe);
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
      Ele fe(m_field);
      fe.getRuleHook = rule;
      boost::shared_ptr<VectorDouble> field_vals_ptr(new VectorDouble());
      boost::shared_ptr<MatrixDouble> diff_field_vals_ptr(new MatrixDouble());
      MatrixDouble inv_jac;
      fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace(inv_jac));
      fe.getOpPtrVector().push_back(new OpSetInvJacH1ForFace(inv_jac));
      fe.getOpPtrVector().push_back(
          new OpCalculateScalarFieldValues("FIELD1", field_vals_ptr));
      fe.getOpPtrVector().push_back(
          new OpCalculateScalarFieldGradient<2>("FIELD1", diff_field_vals_ptr));
      fe.getOpPtrVector().push_back(
          new QuadOpCheck(field_vals_ptr, diff_field_vals_ptr));
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "QUAD", fe);
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

QuadOpCheck::QuadOpCheck(boost::shared_ptr<VectorDouble> &field_vals,
                         boost::shared_ptr<MatrixDouble> &diff_field_vals)
    : OpEle("FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      fieldVals(field_vals), diffFieldVals(diff_field_vals) {}

MoFEMErrorCode QuadOpCheck::doWork(int side, EntityType type,
                                   DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type == MBVERTEX) {
    const int nb_gauss_pts = data.getN().size1();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      double f = ApproxFunction::fun(t_coords(0), t_coords(1));
      constexpr double eps = 1e-6;
      if (std::abs(f - (*fieldVals)[gg]) > eps)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Wrong value %6.4e != %6.4e", f, (*fieldVals)[gg]);
      VectorDouble3 diff_f = ApproxFunction::diff_fun(t_coords(0), t_coords(1));
      for (auto d : {0, 1}) {
        if (std::abs(diff_f[d] - (*diffFieldVals)(d, gg)) > eps)
          SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Wrong directive value (%d) %6.4e != %6.4e", diff_f[d],
                   (*diffFieldVals)(d, gg));
      }

      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
}

QuadOpRhs::QuadOpRhs(SmartPetscObj<Vec> &f)
    : OpEle("FIELD1", "FIELD1", ForcesAndSourcesCore::UserDataOperator::OPROW),
      F(f) {}

MoFEMErrorCode QuadOpRhs::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_dofs = data.getN().size2();
  if (nb_dofs) {
    const int nb_gauss_pts = data.getN().size1();
    VectorDouble nf(nb_dofs);
    nf.clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    auto t_w = getFTensor0IntegrationWeight();
    double vol = getArea();
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      double f = ApproxFunction::fun(t_coords(0), t_coords(1));
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

QuadOpLhs::QuadOpLhs(SmartPetscObj<Mat> &a)
    : OpEle("FIELD1", "FIELD1",
            ForcesAndSourcesCore::UserDataOperator::OPROWCOL),
      A(a) {
  // FIXME: Can be symmetric, is not for simplicity
  sYmm = false;
}

MoFEMErrorCode QuadOpLhs::doWork(int row_side, int col_side,
                                 EntityType row_type, EntityType col_type,
                                 DataForcesAndSourcesCore::EntData &row_data,
                                 DataForcesAndSourcesCore::EntData &col_data) {
  MoFEMFunctionBegin;
  const int row_nb_dofs = row_data.getN().size2();
  const int col_nb_dofs = col_data.getN().size2();
  if (row_nb_dofs && col_nb_dofs) {
    const int nb_gauss_pts = row_data.getN().size1();
    MatrixDouble m(row_nb_dofs, col_nb_dofs);
    m.clear();
    auto t_w = getFTensor0IntegrationWeight();
    double vol = getArea();
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
