/** \file SurfaceSlidingConstrains.hpp
 * \brief Implementing surface sliding constrains
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

#ifndef __SURFACE_SLIDING_CONSTRAINS_HPP__
#define __SURFACE_SLIDING_CONSTRAINS_HPP__

#ifndef WITH_ADOL_C
#error "MoFEM need to be compiled with ADOL-C"
#endif

/** \brief Shape preserving constrains, i.e. nodes sliding on body surface.

  Derivation and implementation of constrains preserving body surface,
  i.e. body shape and volume.

  The idea starts form observation that body shape can be globally characterized
  by constant calculated as volume over its area
  \f[
  \frac{V}{A} = C
  \f]
  Above equation expressed in integral form is
  \f[
  \int_\Omega \textrm{d}V = C \int_\Gamma \textrm{d}S
  \f]
  where notting that,
  \f[
  \frac{1}{3}
  \int_\Omega \textrm{div}[\mathbf{X}] \textrm{d}\Omega
  =
  C \int_\Gamma  \textrm{d}S
  \f]
  and applying Gauss theorem we get
  \f[
  \int_\Gamma
  \mathbf{X}\cdot \frac{\mathbf{N}}{\|\mathbf{N}\|}
  \textrm{d}\Gamma
  =
  3C \int_\Gamma  \textrm{d}S.
  \f]
  Drooping integrals on both sides, and linearizing equation, we get
  \f[
  \frac{\mathbf{N}}{\|\mathbf{N}\|} \cdot \delta \mathbf{X}
  =
  3C - \frac{\mathbf{N}}{\|\mathbf{N}\|}\cdot \mathbf{X}
  \f]
  where \f$\delta \mathbf{X}\f$ is displacement sub-inctrement. Above equation
  is a constrain if satisfied in body shape and volume is conserved. Final form
  of constrain equation is \f[ \mathcal{r} =
  \frac{\mathbf{N}}{\|\mathbf{N}\|}\cdot \mathbf{X}
  -
  \frac{\mathbf{N_0}}{\|\mathbf{N_0}\|}\cdot \mathbf{X}_0 =
  \frac{\mathbf{N}}{\|\mathbf{N}\|}\cdot (\mathbf{X}-\mathbf{X}_0)
  \f]

  In the spirit of finite element method the constrain equation is multiplied
  by shape functions and enforce using Lagrange multiplier method
  \f[
  \int_\Gamma \mathbf{N}^\mathsf{T}_\lambda
   \left(
     \frac{\mathbf{N}}{\|\mathbf{N}\|}\mathbf{N}_\mathbf{X}\cdot
     (\overline{\mathbf{X}}-\overline{\mathbf{X}}_0)
  \right)
   \|\mathbf{N}\|
  \textrm{d}\Gamma
   =
  \mathbf{0}.
  \f]
  Above equation is nonlinear, applying to it Taylor expansion, we can get form
  which can be used with Newton interactive method \f[ \begin{split}
   &\int_\Gamma \mathbf{N}^\mathsf{T}_\lambda
    \left\{
    \mathbf{N}\mathbf{N}_\mathbf{X}
    +
    \left(\mathbf{X}-\mathbf{X}_0\right) \cdot
    \left(
    \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\xi}\right]\cdot\mathbf{B}_\eta
    -
    \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\eta}\right]\cdot\mathbf{B}_\xi
    \right)
    \right\}
    \textrm{d}\Gamma
    \cdot
    \delta
    \overline{\mathbf{X}}\\
    =
    &\int_\Gamma \mathbf{N}^\mathsf{T}_\lambda
    \mathbf{N}\cdot(\mathbf{X}-\mathbf{X}_0)
    \textrm{d}\Gamma
  \end{split}.
  \f]
  Equation expressing forces on shape as result of constrains, as result
  Lagrange multiplier method have following form \f[ \begin{split}
  &\int_\Gamma
  \mathbf{N}^\mathsf{T}_\mathbf{X} \cdot \mathbf{N}
  \mathbf{N}_\lambda
  \textrm{d}\Gamma
  \cdot
  \delta\overline{\lambda}\\
  +
  &\int_\Gamma
  \lambda
  \mathbf{N}^\mathsf{T}_\mathbf{X}
  \left(
  \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\xi}\right]\cdot\mathbf{B}_\eta
  -
  \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\eta}\right]\cdot\mathbf{B}_\xi
  \right)
  \textrm{d}\Gamma
  \delta\overline{\mathbf{X}}\\
  =
  &\int_\Gamma
  \lambda
  \mathbf{N}^\mathsf{T}_\mathbf{X} \cdot \mathbf{N}
  \textrm{d}\Gamma
  \end{split}
  \f]

  Above equations are assembled into global system of equations as following
  \f[
  \left[
    \begin{array}{cc}
        \mathbf{K} + \mathbf{B} & \mathbf{C}^\mathsf{T} \\
        \mathbf{C} + \mathbf{A} & 0
    \end{array}
  \right]
  \left\{
    \begin{array}{c}
      \delta \overline{\mathbf{X}} \\
      \delta \overline{\lambda}
    \end{array}
  \right\}=
  \left[
    \begin{array}{c}
      \mathbf{f} - \mathbf{C}^\mathsf{T}\overline{\lambda} \\
      \overline{\mathbf{r}}
    \end{array}
  \right]
  \f]
  where
  \f[
  \mathbf{C}=
  \int_\Gamma
  \mathbf{N}_\lambda^\mathsf{T}
   \mathbf{N} \cdot
   \mathbf{N}_\mathbf{X}
  \textrm{d}\Gamma,
  \f]
  \f[
  \mathbf{B}=
  \int_\Gamma
  \lambda
  (\mathbf{X}-\mathbf{X}_0)\cdot
  \left(
    \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\xi}\right]\cdot\mathbf{B}_\eta
    -
    \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\eta}\right]\cdot\mathbf{B}_\xi
  \right)
  \textrm{d}\Gamma
  \f]
  and
  \f[
  \mathbf{A}=
  \int_\Gamma
  \mathbf{N}^\mathsf{T}_\lambda
  \left(\mathbf{X}-\mathbf{X}_0\right) \cdot
  \left(
  \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\xi}\right]\cdot\mathbf{B}_\eta
  -
  \textrm{Spin}\left[\frac{\partial\mathbf{X}}{\partial\eta}\right]\cdot\mathbf{B}_\xi
  \right).
  \f]

*/
struct SurfaceSlidingConstrains {

  MoFEM::Interface &mField;

  struct MyTriangleFE : public MoFEM::FaceElementForcesAndSourcesCore {

    Mat B;
    Vec F;

    MyTriangleFE(MoFEM::Interface &m_field)
        : MoFEM::FaceElementForcesAndSourcesCore(m_field), B(PETSC_NULL),
          F(PETSC_NULL) {}
    int getRule(int order) { return 2 * order; };

    MoFEMErrorCode preProcess() {
      MoFEMFunctionBeginHot;

      ierr = MoFEM::FaceElementForcesAndSourcesCore::preProcess();
      CHKERRG(ierr);

      if (B != PETSC_NULL) {
        snes_B = B;
      }

      if (F != PETSC_NULL) {
        snes_f = F;
      }

      switch (ts_ctx) {
      case CTX_TSSETIFUNCTION: {
        if (!F) {
          snes_ctx = CTX_SNESSETFUNCTION;
          snes_f = ts_F;
        }
        break;
      }
      case CTX_TSSETIJACOBIAN: {
        if (!B) {
          snes_ctx = CTX_SNESSETJACOBIAN;
          snes_B = ts_B;
        }
        break;
      }
      default:
        break;
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  boost::shared_ptr<MyTriangleFE> feRhsPtr, feLhsPtr;

  MyTriangleFE &feRhs;
  MyTriangleFE &getLoopFeRhs() { return feRhs; }
  MyTriangleFE &feLhs;
  MyTriangleFE &getLoopFeLhs() { return feLhs; }

  /** \brief Class implemented by user to detect face orientation

   If mesh generated is with surface mesher, usually you don't have to do
   nothing, all elements on the surface have consistent orientation. In case of
   internal faces or if you do something with mesh connectivity which breaks
   orientation on the face, you have to implement method which will set
   orientation to face.

  */
  struct DriverElementOrientation {

    int elementOrientation;

    virtual MoFEMErrorCode
    getElementOrientation(MoFEM::Interface &m_field,
                          const FEMethod *fe_method_ptr) {
      MoFEMFunctionBeginHot;
      elementOrientation = 1;
      MoFEMFunctionReturnHot(0);
    }
  };

  DriverElementOrientation &crackFrontOrientation;

  SurfaceSlidingConstrains(MoFEM::Interface &m_field,
                           DriverElementOrientation &orientation)
      : mField(m_field), feRhsPtr(new MyTriangleFE(m_field)),
        feLhsPtr(new MyTriangleFE(m_field)), feRhs(*feRhsPtr), feLhs(*feLhsPtr),
        crackFrontOrientation(orientation) {}

  struct OpGetActiveDofsLambda
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
    boost::shared_ptr<VectorDouble> activeVariablesPtr;
    OpGetActiveDofsLambda(
        const std::string field_name,
        boost::shared_ptr<VectorDouble> &active_variables_ptr)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPCOL),
          activeVariablesPtr(active_variables_ptr) {}
    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBegin;
      if(type == MBVERTEX) {
        for (unsigned int dd = 0; dd != data.getFieldData().size(); ++dd)
          (*activeVariablesPtr)[dd] = data.getFieldData()[dd];
      }
      MoFEMFunctionReturn(0);
    }
  };

  struct OpGetActiveDofsPositions
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {
    boost::shared_ptr<VectorDouble> activeVariablesPtr;
    OpGetActiveDofsPositions(
        const std::string field_name,
        boost::shared_ptr<VectorDouble> &active_variables_ptr)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPCOL),
          activeVariablesPtr(active_variables_ptr) {}
    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBegin;
      if(type == MBVERTEX) {
        for (unsigned int dd = 0; dd != data.getFieldData().size(); ++dd)
          (*activeVariablesPtr)[3+dd] = data.getFieldData()[dd];
      }
      MoFEMFunctionReturn(0);
    }
  };

  struct OpJacobian
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    const int tAg;
    boost::shared_ptr<VectorDouble> activeVariablesPtr;
    boost::shared_ptr<VectorDouble> resultsPtr;
    boost::shared_ptr<MatrixDouble> jacobianPtr;
    DriverElementOrientation &oRientation;
    bool evaluateJacobian;

    OpJacobian(int tag, const std::string field_name,
               boost::shared_ptr<VectorDouble> &active_variables_ptr,
               boost::shared_ptr<VectorDouble> &results_ptr,
               boost::shared_ptr<MatrixDouble> &jacobian_ptr,
               DriverElementOrientation &orientation, bool evaluate_jacobian)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPCOL),
          tAg(tag), activeVariablesPtr(active_variables_ptr),
          resultsPtr(results_ptr), jacobianPtr(jacobian_ptr),
          oRientation(orientation), evaluateJacobian(evaluate_jacobian) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBegin;
      if (type != MBVERTEX)
        MoFEMFunctionReturnHot(0);

      VectorInt &indices = data.getIndices();

      trace_on(tAg);

      ublas::vector<adouble> lambda_dofs(3);
      for (int dd = 0; dd != 3; ++dd) {
        lambda_dofs[dd] <<= (*activeVariablesPtr)[dd];
      }
      ublas::vector<adouble> position_dofs(9);
      for (int dd = 0; dd != 9; ++dd) {
        position_dofs[dd] <<= (*activeVariablesPtr)[3+dd];
      }

      FTensor::Index<'i', 3> i;
      FTensor::Number<0> N0;
      FTensor::Number<1> N1;
      FTensor::Number<2> N2;

      CHKERR oRientation.getElementOrientation(getFaceFE()->mField,
                                                     getFEMethod());
      int eo = oRientation.elementOrientation;

      int nb_gauss_pts = data.getN().size1();
      int nb_base_functions = data.getN().size2();
      FTensor::Tensor0<double *> t_base1 = data.getFTensor0N();
      FTensor::Tensor0<double *> t_base2 = data.getFTensor0N();
      FTensor::Tensor1<double *, 2> t_diff_base = data.getFTensor1DiffN<2>();
      FTensor::Tensor1<adouble, 3> t_position;
      FTensor::Tensor1<adouble, 3> t_position_ksi;
      FTensor::Tensor1<adouble, 3> t_position_eta;
      adouble lambda;
      FTensor::Tensor1<adouble, 3> t_normal;
      FTensor::Tensor1<adouble, 3> t_delta;

      ublas::vector<adouble> c_vec(3);
      ublas::vector<adouble> f_vec(9);
      c_vec.clear();
      f_vec.clear();

      for (int gg = 0; gg != nb_gauss_pts; ++gg) {

        FTensor::Tensor1<adouble *, 3> t_position_dofs(
            &position_dofs[0], &position_dofs[1], &position_dofs[2], 3);
        FTensor::Tensor0<adouble *> t_lambda_dof(&lambda_dofs[0]);

        t_position(i) = 0;
        t_position_ksi(i) = 0;
        t_position_eta(i) = 0;
        lambda = 0;

        for (int bb = 0; bb != nb_base_functions; ++bb) {
          t_position(i) += t_base1 * t_position_dofs(i);
          t_position_ksi(i) += t_diff_base(N0) * t_position_dofs(i);
          t_position_eta(i) += t_diff_base(N1) * t_position_dofs(i);
          lambda += t_base1 * t_lambda_dof;
          ++t_base1;
          ++t_position_dofs;
          ++t_lambda_dof;
          ++t_diff_base;
        }

        t_delta(i) = t_position(i);
        for (int dd = 0; dd != 3; ++dd) {
          t_delta(dd) -= getCoordsAtGaussPts()(gg,dd);
        }

        t_normal(0) = t_position_ksi(1) * t_position_eta(2) -
                      t_position_ksi(2) * t_position_eta(1);
        t_normal(1) = t_position_ksi(2) * t_position_eta(0) -
                      t_position_ksi(0) * t_position_eta(2);
        t_normal(2) = t_position_ksi(0) * t_position_eta(1) -
                      t_position_ksi(1) * t_position_eta(0);

        double w = getGaussPts()(2, gg) * 0.5;
        adouble val;
        FTensor::Tensor0<adouble *> t_c(&c_vec[0]);
        FTensor::Tensor1<adouble *, 3> t_f(&f_vec[0], &f_vec[1], &f_vec[2], 3);

        for (int bb = 0; bb != nb_base_functions; ++bb) {
          if (indices[bb] != -1) {
            val = w * eo * t_base2;
            t_c += val * t_normal(i) * t_delta(i);
            val *= lambda;
            t_f(i) += val * t_normal(i);
          }
          ++t_c;
          ++t_f;
          ++t_base2;
        }
      }

      for (int rr = 0; rr != 3; ++rr) {
        c_vec[rr] >>= (*resultsPtr)[rr];
      }
      for (int rr = 0; rr != 9; ++rr) {
        f_vec(rr) >>= (*resultsPtr)[3 + rr];
      }

      trace_off();

      if (evaluateJacobian) {
        double *jac_ptr[3 + 9];
        for (int rr = 0; rr != 3 + 9; ++rr) {
          jac_ptr[rr] = &(*jacobianPtr)(rr, 0);
        }
        // play recorder for jacobians
        int r =
            ::jacobian(tAg, 3 + 9, 3 + 9, &(*activeVariablesPtr)[0], jac_ptr);
        if (r < 0) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
                  "ADOL-C function evaluation with error");
        }
      }

      MoFEMFunctionReturn(0);
    }
  };

  struct OpAssembleRhs
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    boost::shared_ptr<VectorDouble> resultsPtr;

    OpAssembleRhs(const std::string field_name,
                           boost::shared_ptr<VectorDouble> &results_ptr)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPROW),
          resultsPtr(results_ptr) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBegin;
      if (type != MBVERTEX)
        MoFEMFunctionReturnHot(0);
      VectorInt& indices = data.getIndices();
      int shift = 0;
      if (indices.empty()) {
        MoFEMFunctionReturnHot(0);
      } else if (indices.size() == 3) {
        shift = 0;
      } else if (indices.size() == 9) {
        shift = 3;
      } else {
        SETERRQ1(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                 "Data inconsistency nb of indices %d", indices.size());
      }
      CHKERR VecSetOption(getFEMethod()->snes_f, VEC_IGNORE_NEGATIVE_INDICES,
                          PETSC_TRUE);
      CHKERR VecSetValues(getFEMethod()->snes_f, indices.size(), &indices[0],
                          &(*resultsPtr)[shift], ADD_VALUES);
      MoFEMFunctionReturn(0);
    }
  };

  struct OpAssembleLhs
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    boost::shared_ptr<MatrixDouble> jacobianPtr;

    OpAssembleLhs(const std::string field_name_row,
                  const std::string field_name_col,
                  boost::shared_ptr<MatrixDouble> &jacobian_ptr)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name_row, field_name_col, UserDataOperator::OPROWCOL),
          jacobianPtr(jacobian_ptr) {
      sYmm = false;
    }

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type,
                          DataForcesAndSourcesCore::EntData &row_data,
                          DataForcesAndSourcesCore::EntData &col_data) {
      MoFEMFunctionBegin;
      if (row_type != MBVERTEX)
        MoFEMFunctionReturnHot(0);
      if (col_type != MBVERTEX)
        MoFEMFunctionReturnHot(0);
      VectorInt &row_indices = row_data.getIndices();
      VectorInt &col_indices = col_data.getIndices();
      if (row_indices.empty() || col_indices.empty())
        MoFEMFunctionReturnHot(0);
      int shift_row = 0;
      if (row_indices.size() == 3) {
        shift_row = 0;
      } else if (row_indices.size() == 9) {
        shift_row = 3;
      } else {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      int shift_col = 0;
      if (col_indices.size() == 3) {
        shift_col = 0;
      } else if (col_indices.size() == 9) {
        shift_col = 3;
      } else {
        SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      MatrixDouble jac(row_indices.size(), col_indices.size());
      for(int rr = 0;rr!=row_indices.size();++rr) {
        for(int cc = 0;cc!=col_indices.size();++cc) {
          jac(rr, cc) = (*jacobianPtr)(shift_row + rr, shift_col + cc);
        }
      }
      ierr = MatSetValues(getFEMethod()->snes_B, row_indices.size(),
                          &row_indices[0], col_indices.size(), &col_indices[0],
                          &jac(0, 0), ADD_VALUES);
      MoFEMFunctionReturn(0);
    }
  };

  MoFEMErrorCode setOperators(int tag,
                              const std::string lagrange_multipliers_field_name,
                              const std::string material_field_name) {
    MoFEMFunctionBegin;

    boost::shared_ptr<VectorDouble> active_variables_ptr(new VectorDouble(3+9));
    boost::shared_ptr<VectorDouble> results_ptr(new VectorDouble(3+9));
    boost::shared_ptr<MatrixDouble> jacobian_ptr(new MatrixDouble(3+9,3+9));

    feRhs.getOpPtrVector().clear();
    feRhs.getOpPtrVector().push_back(new OpGetActiveDofsLambda(
        lagrange_multipliers_field_name, active_variables_ptr));
    feRhs.getOpPtrVector().push_back(new OpGetActiveDofsPositions(
        material_field_name, active_variables_ptr));
    feRhs.getOpPtrVector().push_back(new OpJacobian(
        tag, lagrange_multipliers_field_name, active_variables_ptr, results_ptr,
        jacobian_ptr, crackFrontOrientation, false));
    feRhs.getOpPtrVector().push_back(
        new OpAssembleRhs(lagrange_multipliers_field_name, results_ptr));
    feRhs.getOpPtrVector().push_back(new OpAssembleRhs(
        material_field_name, results_ptr));

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().clear();
    feLhs.getOpPtrVector().push_back(new OpGetActiveDofsLambda(
        lagrange_multipliers_field_name, active_variables_ptr));
    feLhs.getOpPtrVector().push_back(new OpGetActiveDofsPositions(
        material_field_name, active_variables_ptr));
    feLhs.getOpPtrVector().push_back(new OpJacobian(
        tag, lagrange_multipliers_field_name, active_variables_ptr, results_ptr,
        jacobian_ptr, crackFrontOrientation, true));
    feLhs.getOpPtrVector().push_back(new OpAssembleLhs(
        lagrange_multipliers_field_name, material_field_name, jacobian_ptr));
    feLhs.getOpPtrVector().push_back(new OpAssembleLhs(
        material_field_name, lagrange_multipliers_field_name, jacobian_ptr));
    feLhs.getOpPtrVector().push_back(new OpAssembleLhs(
        material_field_name, material_field_name, jacobian_ptr));

    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode
  setOperatorsConstrainOnly(int tag,
                            const std::string lagrange_multipliers_field_name,
                            const std::string material_field_name) {
    MoFEMFunctionBeginHot;

    boost::shared_ptr<VectorDouble> active_variables_ptr(
        new VectorDouble(3 + 9));
    boost::shared_ptr<VectorDouble> results_ptr(new VectorDouble(3 + 9));
    boost::shared_ptr<MatrixDouble> jacobian_ptr(
        new MatrixDouble(3 + 9, 3 + 9));

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().clear();
    feLhs.getOpPtrVector().push_back(new OpGetActiveDofsLambda(
        lagrange_multipliers_field_name, active_variables_ptr));
    feLhs.getOpPtrVector().push_back(new OpGetActiveDofsPositions(
        material_field_name, active_variables_ptr));
    feLhs.getOpPtrVector().push_back(new OpJacobian(
        tag, lagrange_multipliers_field_name, active_variables_ptr, results_ptr,
        jacobian_ptr, crackFrontOrientation, true));
    feLhs.getOpPtrVector().push_back(new OpAssembleLhs(
        lagrange_multipliers_field_name, material_field_name, jacobian_ptr));

    MoFEMFunctionReturnHot(0);
  }
  
};

#endif // __SURFACE_SLIDING_CONSTRAINS_HPP__
