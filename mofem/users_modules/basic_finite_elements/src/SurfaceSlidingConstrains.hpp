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
   inetranl faces or if you do something with mesh connectivity which breaks
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

  struct AuxFunctions {

    bool useProjectionFromCrackFront;
    AuxFunctions() : useProjectionFromCrackFront(false) {}

    MatrixDouble N;
    MatrixDouble Bksi;
    MatrixDouble Beta;

    VectorDouble pOsition;
    VectorDouble dXdKsi;
    VectorDouble dXdEta;
    MatrixDouble spinKsi;
    MatrixDouble spinEta;
    VectorDouble nOrmal;
    MatrixDouble sPin;
    double aRea;
    double lAmbda;

    std::vector<bool> nodesWithoutLambda;

    static MoFEMErrorCode calcSpin(MatrixDouble &spin, VectorDouble &vec) {
      MoFEMFunctionBeginHot;
      spin.resize(3, 3, false);
      spin.clear();
      spin(0, 1) = -vec[2];
      spin(0, 2) = +vec[1];
      spin(1, 0) = +vec[2];
      spin(1, 2) = -vec[0];
      spin(2, 0) = -vec[1];
      spin(2, 1) = +vec[0];
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode matrixN(int gg, DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;
      try {

        int nb_dofs = data.getN().size2();
        N.resize(3, 3 * nb_dofs, false);
        N.clear();
        for (int ii = 0; ii < nb_dofs; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            N(jj, ii * 3 + jj) = data.getN(gg)[ii];
          }
        }

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode matrixB(int gg, DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;
      try {
        int nb_dofs = data.getN().size2();
        Bksi.resize(3, 3 * nb_dofs, false);
        Bksi.clear();
        for (int ii = 0; ii < nb_dofs; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Bksi(jj, ii * 3 + jj) = data.getDiffN(gg)(ii, 0);
          }
        }
        Beta.resize(3, 3 * nb_dofs);
        Beta.clear();
        for (int ii = 0; ii < nb_dofs; ii++) {
          for (int jj = 0; jj < 3; jj++) {
            Beta(jj, ii * 3 + jj) = data.getDiffN(gg)(ii, 1);
          }
        }
      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode calculateNormal() {
      MoFEMFunctionBeginHot;

      try {
        sPin.resize(3, 3, false);
        ierr = calcSpin(sPin, dXdKsi);
        CHKERRG(ierr);
        nOrmal.resize(3, false);
        noalias(nOrmal) = 0.5 * prod(sPin, dXdEta);
        aRea = norm_2(nOrmal);
      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  std::vector<AuxFunctions> cUrrent;

  /** \brief Operator calculate material positions and tangent vectors to
   * element surface
   */
  struct OpPositions
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpPositions(const std::string field_name, std::vector<AuxFunctions> &aux,
                DriverElementOrientation &orientation)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPCOL),
          aUx(aux), oRientation(orientation) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;

      try {

        int nb_dofs = data.getFieldData().size();
        if (nb_dofs == 0) {
          MoFEMFunctionReturnHot(0);
        }

        int nb_gauss_pts = data.getN().size1();
        if (type == MBVERTEX) {
          aUx.resize(nb_gauss_pts);
          ierr = oRientation.getElementOrientation(getFaceFE()->mField,
                                                   getFEMethod());
          CHKERRG(ierr);
          for (int gg = 0; gg < nb_gauss_pts; gg++) {
            aUx[gg].pOsition.resize(3, false);
            aUx[gg].dXdKsi.resize(3, false);
            aUx[gg].dXdEta.resize(3, false);
            aUx[gg].pOsition.clear();
            aUx[gg].dXdKsi.clear();
            aUx[gg].dXdEta.clear();
          }
        }

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          ierr = aUx[gg].matrixN(gg, data);
          CHKERRG(ierr);
          ierr = aUx[gg].matrixB(gg, data);
          CHKERRG(ierr);
          noalias(aUx[gg].pOsition) += prod(aUx[gg].N, data.getFieldData());
          noalias(aUx[gg].dXdKsi) += prod(aUx[gg].Bksi, data.getFieldData());
          noalias(aUx[gg].dXdEta) += prod(aUx[gg].Beta, data.getFieldData());
        }

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Operator calculate Lagrange multiplier values at integration points
   */
  struct OpLambda
      : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;

    OpLambda(const std::string field_name, vector<AuxFunctions> &aux)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPROW),
          aUx(aux) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;
      //

      try {

        int nb_gauss_pts = data.getN().size1();
        if (type == MBVERTEX) {
          if (aUx.size() != nb_gauss_pts) {
            SETERRQ2(getFaceFE()->mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
                     "Size of aUx should be equal to number of integration "
                     "point but is %d != %d",
                     aUx.size(), nb_gauss_pts);
          }
          for (int gg = 0; gg != nb_gauss_pts; gg++) {
            aUx[gg].lAmbda = 0;
          }
        }

        int nb_dofs = data.getFieldData().size();
        if (nb_dofs == 0) {
          MoFEMFunctionReturnHot(0);
        }

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          // cerr << data.getN(gg) << endl;
          // cerr << data.getFieldData() << endl;

          aUx[gg].lAmbda +=
              inner_prod(data.getN(gg, nb_dofs), data.getFieldData());

          if (aUx[gg].lAmbda != aUx[gg].lAmbda) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "NaN value");
          }
        }

        if (type == MBVERTEX) {
          aUx[0].nodesWithoutLambda.resize(nb_dofs);
          for (int dd = 0; dd < nb_dofs; dd++) {
            if (data.getIndices()[dd] == -1) {
              aUx[0].nodesWithoutLambda[dd] = true;
            } else {
              aUx[0].nodesWithoutLambda[dd] = false;
            }
          }
        }

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Operator calculate \f$\overline{\lambda}\mathbf{C}^\mathsf{T}\f$
   */
  struct OpF : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpF(const std::string field_name, std::vector<AuxFunctions> &aux,
        DriverElementOrientation &orientation)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPROW),
          aUx(aux), oRientation(orientation) {}

    VectorDouble c;
    VectorDouble nF;
    ublas::vector<int> rowIndices;

    MoFEMErrorCode doWork(int row_side, EntityType row_type,
                          DataForcesAndSourcesCore::EntData &row_data) {

      try {

        unsigned int nb_dofs = row_data.getIndices().size();
        if (nb_dofs == 0) {
          MoFEMFunctionReturnHot(0);
        }
        // if (nb_dofs != row_data.getIndices().size()) {
        //   SETERRQ4(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        //            "FE %s Field %s data inconsistency %d != %d",
        //            getNumeredEntFiniteElementPtr()->getName().c_str(),
        //            rowFieldName.c_str(), nb_dofs, row_data.getIndices().size());
        // }
        int nb_gauss_pts = row_data.getN().size1();

        if (row_type == MBVERTEX) {
          for (int gg = 0; gg < nb_gauss_pts; gg++) {
            aUx[gg].nOrmal.resize(3, false);
            aUx[gg].nOrmal.clear();
            aUx[gg].calculateNormal();
          }
        }

        int eo = oRientation.elementOrientation;

        c.resize(nb_dofs, false);
        nF.resize(nb_dofs, false);

        c.clear();
        nF.clear();

        for (int gg = 0; gg < nb_gauss_pts; gg++) {

          double val = getGaussPts()(2, gg);
          noalias(c) = prod(aUx[gg].nOrmal, aUx[gg].N);
          noalias(nF) += val * eo * c * aUx[gg].lAmbda;
        }

        rowIndices.resize(nb_dofs, false);
        noalias(rowIndices) = row_data.getIndices();

        if (row_type == MBVERTEX) {
          int nb_dofs = aUx[0].nodesWithoutLambda.size();
          for (int dd = 0; dd < nb_dofs; dd++) {
            if (aUx[0].nodesWithoutLambda[dd]) {
              for (int jj = 0; jj < 3; jj++) {
                rowIndices[dd * 3 + jj] = -1;
              }
            }
          }
        }

        int *indices_ptr = &rowIndices[0];

        ierr = VecSetOption(getFEMethod()->snes_f, VEC_IGNORE_NEGATIVE_INDICES,
                            PETSC_TRUE);
        CHKERRG(ierr);
        ierr = VecSetValues(getFEMethod()->snes_f, nb_dofs, indices_ptr, &nF[0],
                            ADD_VALUES);
        CHKERRG(ierr);

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Calculate residual \$ \mathbf{g} = \int_\Gamma \mathbf{N}_\lambda
   * \mathbf{N}\left( \mathbf{X}-\mathbf{X}_0 \right) \$
   */
  struct OpG : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpG(const std::string field_name, std::vector<AuxFunctions> &aux,
        DriverElementOrientation &orientation)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, UserDataOperator::OPROW),
          aUx(aux), oRientation(orientation) {}

    VectorDouble g, dElta;

    MoFEMErrorCode doWork(int row_side, EntityType row_type,
                          DataForcesAndSourcesCore::EntData &row_data) {

      try {

        int nb_dofs = row_data.getFieldData().size();
        if (nb_dofs == 0) {
          MoFEMFunctionReturnHot(0);
        }
        int nb_gauss_pts = row_data.getN().size1();

        if (row_type == MBVERTEX) {
          for (int gg = 0; gg < nb_gauss_pts; gg++) {
            aUx[gg].nOrmal.resize(3, false);
            aUx[gg].nOrmal.clear();
            aUx[gg].calculateNormal();
          }
        }

        int eo = oRientation.elementOrientation;

        dElta.resize(3);
        dElta.clear();
        g.resize(nb_dofs, false);
        g.clear();

        for (int gg = 0; gg < nb_gauss_pts; gg++) {

          noalias(dElta) = aUx[gg].pOsition;
          for (int dd = 0; dd < 3; dd++) {
            dElta[dd] -= getCoordsAtGaussPts()(gg, dd);
          }

          double val = getGaussPts()(2, gg);
          double r = inner_prod(aUx[gg].nOrmal, dElta);
          noalias(g) += val * eo * row_data.getN(gg) * r;
        }

        int *indices_ptr = &row_data.getIndices()[0];

        ierr = VecSetOption(getFEMethod()->snes_f, VEC_IGNORE_NEGATIVE_INDICES,
                            PETSC_TRUE);
        CHKERRG(ierr);

        ierr = VecSetValues(getFEMethod()->snes_f, nb_dofs, indices_ptr, &g[0],
                            ADD_VALUES);
        CHKERRG(ierr);

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Operator calculating matrix \b C
   */
  struct OpC : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;
    bool assembleTranspose;

    OpC(const std::string lambda_field_name,
        const std::string positions_field_name, std::vector<AuxFunctions> &aux,
        DriverElementOrientation &orientation, bool assemble_transpose)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              lambda_field_name, positions_field_name,
              UserDataOperator::OPROWCOL),
          aUx(aux), oRientation(orientation),
          assembleTranspose(assemble_transpose) {
      sYmm = false;
    }

    VectorDouble c;
    MatrixDouble C;
    MatrixDouble transC;
    ublas::vector<int> transRowIndices;

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type,
                          DataForcesAndSourcesCore::EntData &row_data,
                          DataForcesAndSourcesCore::EntData &col_data) {
      MoFEMFunctionBeginHot;

      if (col_type != MBVERTEX) {
        MoFEMFunctionReturnHot(0);
      }

      try {

        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if (!nb_row || !nb_col) {
          MoFEMFunctionReturnHot(0);
        }

        int nb_gauss_pts = row_data.getN().size1();

        if (row_type == MBVERTEX) {
          for (int gg = 0; gg < nb_gauss_pts; gg++) {
            aUx[gg].nOrmal.resize(3, false);
            aUx[gg].nOrmal.clear();
            aUx[gg].calculateNormal();
          }
        }

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          ierr = aUx[gg].matrixN(gg, col_data);
          CHKERRG(ierr);
        }

        int eo = oRientation.elementOrientation;

        c.resize(nb_col, false);
        C.resize(nb_row, nb_col, false);
        C.clear();
        for (int gg = 0; gg < nb_gauss_pts; gg++) {

          noalias(c) = prod(aUx[gg].nOrmal, aUx[gg].N);
          double val = getGaussPts()(2, gg);
          noalias(C) += val * eo * outer_prod(row_data.getN(gg), c);
        }

        int *row_indices_ptr = &row_data.getIndices()[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        for (unsigned int n1 = 0; n1 != C.size1(); n1++) {
          for (unsigned int n2 = 0; n2 != C.size1(); n2++) {
            if (C(n1, n2) != C(n1, n2)) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "NaN value");
            }
          }
        }

        ierr = MatSetValues(getFEMethod()->snes_B, nb_row, row_indices_ptr,
                            nb_col, col_indices_ptr, &C(0, 0), ADD_VALUES);
        CHKERRG(ierr);

        if (assembleTranspose) {

          transC.resize(nb_col, nb_row);
          noalias(transC) = trans(C);

          transRowIndices.resize(nb_col, false);
          noalias(transRowIndices) = col_data.getIndices();
          int *trans_row_indices_ptr = &transRowIndices[0];
          if (row_type == MBVERTEX) {
            int nb_dofs = aUx[0].nodesWithoutLambda.size();
            for (int dd = 0; dd < nb_dofs; dd++) {
              if (aUx[0].nodesWithoutLambda[dd]) {
                for (int jj = 0; jj < 3; jj++) {
                  transRowIndices[dd * 3 + jj] = -1;
                }
              }
            }
          }

          ierr =
              MatSetValues(getFEMethod()->snes_B, nb_col, trans_row_indices_ptr,
                           nb_row, row_indices_ptr, &transC(0, 0), ADD_VALUES);
          CHKERRG(ierr);
        }

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Operator calculating matrix \b B
   */
  struct OpB : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpB(const std::string field_name, std::vector<AuxFunctions> &aux,
        DriverElementOrientation &orientation)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              field_name, field_name, UserDataOperator::OPROWCOL),
          aUx(aux), oRientation(orientation) {
      sYmm = false;
    }

    MatrixDouble spindXdKsi, spindXdEta;
    MatrixDouble NdNormal;
    MatrixDouble dNormal;

    MatrixDouble B;
    ublas::vector<int> rowIndices;

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type,
                          DataForcesAndSourcesCore::EntData &row_data,
                          DataForcesAndSourcesCore::EntData &col_data) {
      MoFEMFunctionBeginHot;

      try {

        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if (!nb_row || !nb_col) {
          MoFEMFunctionReturnHot(0);
        }

        int nb_gauss_pts = row_data.getN().size1();
        int eo = oRientation.elementOrientation;

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          ierr = aUx[gg].matrixN(gg, row_data);
          CHKERRG(ierr);
          ierr = aUx[gg].matrixB(gg, col_data);
          CHKERRG(ierr);
        }

        spindXdKsi.resize(3, 3, false);
        spindXdEta.resize(3, 3, false);
        dNormal.resize(3, nb_col, false);
        NdNormal.resize(nb_row, nb_col, false);

        B.resize(nb_row, nb_col, false);
        B.clear();
        for (int gg = 0; gg < nb_gauss_pts; gg++) {

          ierr = AuxFunctions::calcSpin(spindXdKsi, aUx[gg].dXdKsi);
          CHKERRG(ierr);
          ierr = AuxFunctions::calcSpin(spindXdEta, aUx[gg].dXdEta);
          CHKERRG(ierr);

          noalias(dNormal) = eo * (prod(spindXdKsi, aUx[gg].Beta) -
                                   prod(spindXdEta, aUx[gg].Bksi));
          noalias(NdNormal) = prod(trans(aUx[gg].N), dNormal);

          double val = getGaussPts()(2, gg);
          noalias(B) += val * aUx[gg].lAmbda * NdNormal;
        }

        rowIndices.resize(nb_row, false);
        noalias(rowIndices) = row_data.getIndices();
        if (row_type == MBVERTEX) {
          int nb_dofs = aUx[0].nodesWithoutLambda.size();
          for (int dd = 0; dd < nb_dofs; dd++) {
            if (aUx[0].nodesWithoutLambda[dd]) {
              for (int jj = 0; jj < 3; jj++) {
                rowIndices[dd * 3 + jj] = -1;
              }
            }
          }
        }

        int *row_indices_ptr = &rowIndices[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        ierr = MatSetValues(getFEMethod()->snes_B, nb_row, row_indices_ptr,
                            nb_col, col_indices_ptr, &B(0, 0), ADD_VALUES);
        CHKERRG(ierr);

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Operator calculating matrix \b A
   */
  struct OpA : public MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator {

    std::vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpA(const std::string lagrange_multipliers_field_name,
        const std::string field_name, std::vector<AuxFunctions> &aux,
        DriverElementOrientation &orientation)
        : MoFEM::FaceElementForcesAndSourcesCore::UserDataOperator(
              lagrange_multipliers_field_name, field_name,
              UserDataOperator::OPROWCOL),
          aUx(aux), oRientation(orientation) {
      sYmm = false;
    }

    VectorDouble XdNormal;
    MatrixDouble spindXdKsi, spindXdEta;
    MatrixDouble dNormal, NXdNormal;
    VectorDouble dElta;

    MatrixDouble A;

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type,
                          DataForcesAndSourcesCore::EntData &row_data,
                          DataForcesAndSourcesCore::EntData &col_data) {
      MoFEMFunctionBeginHot;

      try {

        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if (!nb_row || !nb_col) {
          MoFEMFunctionReturnHot(0);
        }

        int nb_gauss_pts = row_data.getN().size1();
        int eo = oRientation.elementOrientation;

        for (int gg = 0; gg < nb_gauss_pts; gg++) {
          ierr = aUx[gg].matrixB(gg, col_data);
          CHKERRG(ierr);
        }

        XdNormal.resize(nb_col, false);
        spindXdKsi.resize(3, 3, false);
        spindXdEta.resize(3, 3, false);
        dNormal.resize(3, nb_col, false);
        NXdNormal.resize(nb_row, nb_col, false);
        dElta.resize(3);

        A.resize(nb_row, nb_col, false);
        A.clear();
        for (int gg = 0; gg < nb_gauss_pts; gg++) {

          noalias(dElta) = aUx[gg].pOsition;
          for (int dd = 0; dd < 3; dd++) {
            dElta[dd] -= getCoordsAtGaussPts()(gg, dd);
          }

          ierr = AuxFunctions::calcSpin(spindXdKsi, aUx[gg].dXdKsi);
          CHKERRG(ierr);
          ierr = AuxFunctions::calcSpin(spindXdEta, aUx[gg].dXdEta);
          CHKERRG(ierr);

          noalias(dNormal) = eo * (prod(spindXdKsi, aUx[gg].Beta) -
                                   prod(spindXdEta, aUx[gg].Bksi));
          noalias(XdNormal) = prod(trans(dElta), dNormal);
          noalias(NXdNormal) = outer_prod(row_data.getN(gg), XdNormal);

          double val = getGaussPts()(2, gg);
          noalias(A) += val * NXdNormal;
        }

        int *row_indices_ptr = &row_data.getIndices()[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        ierr = MatSetValues(getFEMethod()->snes_B, nb_row, row_indices_ptr,
                            nb_col, col_indices_ptr, &A(0, 0), ADD_VALUES);
        CHKERRG(ierr);

      } catch (const std::exception &ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF, 1, ss.str().c_str());
      }

      MoFEMFunctionReturnHot(0);
    }
  };

  /** \brief Driver function setting operators to calculate \b C matrix only
   */
  MoFEMErrorCode
  setOperatorsCOnly(const std::string lagrange_multipliers_field_name,
                    const std::string material_field_name) {
    MoFEMFunctionBeginHot;

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
        new OpPositions(material_field_name, cUrrent, crackFrontOrientation));
    feLhs.getOpPtrVector().push_back(
        new OpLambda(lagrange_multipliers_field_name, cUrrent));
    feLhs.getOpPtrVector().push_back(new OpC(lagrange_multipliers_field_name,
                                             material_field_name, cUrrent,
                                             crackFrontOrientation, false));

    MoFEMFunctionReturnHot(0);
  }

  /** \brief Driver function setting operators to calculate nonlinear problems
   * with sliding points on the surface
   */
  MoFEMErrorCode setOperatorsWithLinearGeometry(
      const std::string lagrange_multipliers_field_name,
      const std::string material_field_name, bool assemble_transpose,
      bool add_nonlinear_term) {
    MoFEMFunctionBeginHot;

    // Adding operators to calculate the right hand side
    feRhs.getOpPtrVector().push_back(
        new OpPositions(material_field_name, cUrrent, crackFrontOrientation));
    feRhs.getOpPtrVector().push_back(
        new OpLambda(lagrange_multipliers_field_name, cUrrent));
    feRhs.getOpPtrVector().push_back(
        new OpF(material_field_name, cUrrent, crackFrontOrientation));
    feRhs.getOpPtrVector().push_back(new OpG(lagrange_multipliers_field_name,
                                             cUrrent, crackFrontOrientation));

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
        new OpPositions(material_field_name, cUrrent, crackFrontOrientation));
    feLhs.getOpPtrVector().push_back(
        new OpLambda(lagrange_multipliers_field_name, cUrrent));
    feLhs.getOpPtrVector().push_back(
        new OpC(lagrange_multipliers_field_name, material_field_name, cUrrent,
                crackFrontOrientation, assemble_transpose));
    feLhs.getOpPtrVector().push_back(
        new OpB(material_field_name, cUrrent, crackFrontOrientation));
    feLhs.getOpPtrVector().push_back(new OpA(lagrange_multipliers_field_name,
                                             material_field_name, cUrrent,
                                             crackFrontOrientation));

    MoFEMFunctionReturnHot(0);
  }
};

#endif // __SURFACE_SLIDING_CONSTRAINS_HPP__
