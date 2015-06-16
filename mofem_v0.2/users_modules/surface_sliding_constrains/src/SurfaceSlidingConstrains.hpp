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
  where \f$\delta \mathbf{X}\f$ is displacement sub-inctrement. Above equation is a
  constrain if satisfied in body shape and volume is conserved. Final form of constrain equation
  is
  \f[
  \mathcal{r} =
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
  Above equation is nonlinear, applying to it Taylor expansion, we can get form which
  can be used with Newton interactive method
  \f[
  \begin{split}
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
  Equation expressing forces on shape as result of constrains, as result Lagrange multiplier
  method have following form
  \f[
  \begin{split}
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
  \mathbf{N}^\mathsf{T}_\mathbf{X}
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

  FieldInterface &mField;

  struct MyTriangleFE: public FaceElementForcesAndSourcesCore {

    Mat B;
    Vec F;

    MyTriangleFE(FieldInterface &m_field):
    FaceElementForcesAndSourcesCore(m_field),
    B(PETSC_NULL),
    F(PETSC_NULL)
    {}
    int getRule(int order) { return order; };

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      ierr = FaceElementForcesAndSourcesCore::preProcess(); CHKERRQ(ierr);

      if(B != PETSC_NULL) {
        snes_B = B;
      }

      if(F != PETSC_NULL) {
        snes_f = F;
      }

      switch (ts_ctx) {
        case CTX_TSSETIFUNCTION: {
          if(!F) {
            snes_ctx = CTX_SNESSETFUNCTION;
            snes_f = ts_F;
          }
          break;
        }
        case CTX_TSSETIJACOBIAN: {
          if(!B) {
            snes_ctx = CTX_SNESSETJACOBIAN;
            snes_B = ts_B;
          }
          break;
        }
        default:
        break;
      }
      PetscFunctionReturn(0);
    }

  };

  MyTriangleFE feRhs;
  MyTriangleFE& getLoopFeRhs() { return feRhs; }
  MyTriangleFE feLhs;
  MyTriangleFE& getLoopFeLhs() { return feLhs ; }

  /** \brief Class implemented by user to detect face orientation

   If mesh generated is with surface mesher, usually you don't have to do nothing, all elements
   on the surface have consistent orientation. In case of inetranl faces or if you do
   something with mesh connectivity which breaks orientation on the face, you have to
   implement method which will set orientation to face.

  */
  struct DriverElementOrientation {

    int elementOrientation;

    virtual PetscErrorCode getElementOrientation(FieldInterface &m_field,const FEMethod *fe_method_ptr) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  };
  DriverElementOrientation &crackFrontOrientation;

  SurfaceSlidingConstrains(FieldInterface &m_field,DriverElementOrientation &orientation):
  mField(m_field),
  feRhs(m_field),
  feLhs(m_field),
  crackFrontOrientation(orientation)
  {}

  struct AuxFunctions {

    bool useProjectionFromCrackFront;
    AuxFunctions():
    useProjectionFromCrackFront(false)
    {}

    ublas::matrix<double> N;
    ublas::matrix<double> Bksi;
    ublas::matrix<double> Beta;

    ublas::vector<double> pOsition;
    ublas::vector<double> dXdKsi;
    ublas::vector<double> dXdEta;
    ublas::matrix<double> spinKsi;
    ublas::matrix<double> spinEta;
    ublas::vector<double> nOrmal;
    ublas::matrix<double> sPin;
    double aRea;
    double lAmbda;

    vector<bool> nodesWithoutLambda;

    static PetscErrorCode calcSpin(
      ublas::matrix<double> &spin,ublas::vector<double> &vec
    ) {
      PetscFunctionBegin;
      spin.resize(3,3,false);
      spin.clear();
      spin(0,1) = -vec[2];
      spin(0,2) = +vec[1];
      spin(1,0) = +vec[2];
      spin(1,2) = -vec[0];
      spin(2,0) = -vec[1];
      spin(2,1) = +vec[0];
      PetscFunctionReturn(0);
    }


    PetscErrorCode matrixN(int gg,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {

        int nb_dofs = data.getN().size2();
        N.resize(3,3*nb_dofs,false);
        N.clear();
        for(int ii = 0;ii<nb_dofs;ii++) {
          for(int jj = 0;jj<3;jj++) {
            N(jj,ii*3+jj) = data.getN(gg)[ii];
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode matrixB(int gg,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      try {
        int nb_dofs = data.getN().size2();
        Bksi.resize(3,3*nb_dofs,false);
        Bksi.clear();
        for(int ii = 0;ii<nb_dofs;ii++) {
          for(int jj = 0;jj<3;jj++) {
            Bksi(jj,ii*3+jj) = data.getDiffN(gg)(ii,0);
          }
        }
        Beta.resize(3,3*nb_dofs);
        Beta.clear();
        for(int ii = 0;ii<nb_dofs;ii++) {
          for(int jj = 0;jj<3;jj++) {
            Beta(jj,ii*3+jj) = data.getDiffN(gg)(ii,1);
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode calculateNormal() {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      try {
        sPin.resize(3,3,false);
        ierr = calcSpin(sPin,dXdKsi); CHKERRQ(ierr);
        nOrmal.resize(3,false);
        noalias(nOrmal) = 0.5*prod(sPin,dXdEta);
        aRea = norm_2(nOrmal);
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }
  };

  vector<AuxFunctions> cUrrent;

  /** \brief Operator calculate material positions and tangent vectors to element surface
   */
  struct OpPositions: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpPositions(const string field_name,vector<AuxFunctions> &aux,DriverElementOrientation &orientation):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    aUx(aux),
    oRientation(orientation)
    {}

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {

        int nb_dofs = data.getFieldData().size();
        if(nb_dofs == 0) {
          PetscFunctionReturn(0);
        }
        int nb_gauss_pts = data.getN().size1();

        aUx.resize(nb_gauss_pts);

        if(type == MBVERTEX) {
          oRientation.getElementOrientation(getTriFE()->mField,getFEMethod());
        }

        if(type == MBVERTEX) {
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            aUx[gg].pOsition.resize(3,false);
            aUx[gg].pOsition.clear();
            aUx[gg].dXdKsi.resize(3,false);
            aUx[gg].dXdKsi.clear();
            aUx[gg].dXdEta.resize(3,false);
            aUx[gg].dXdEta.clear();
          }
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = aUx[gg].matrixN(gg,data); CHKERRQ(ierr);
          noalias(aUx[gg].pOsition) += prod(aUx[gg].N,data.getFieldData());
          ierr = aUx[gg].matrixB(gg,data); CHKERRQ(ierr);
          noalias(aUx[gg].dXdKsi) += prod(aUx[gg].Bksi,data.getFieldData());
          noalias(aUx[gg].dXdEta) += prod(aUx[gg].Beta,data.getFieldData());
        }


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Operator calculate Lagrange multiplier values at integration points
  */
  struct OpLambda: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;

    OpLambda(const string field_name,vector <AuxFunctions> &aux):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    aUx(aux) {
    }


    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      //PetscErrorCode ierr;

      try {

        int nb_dofs = data.getFieldData().size();
        if(nb_dofs == 0) {
          PetscFunctionReturn(0);
        }
        int nb_gauss_pts = data.getN().size1();
        aUx.resize(nb_gauss_pts);

        if(type == MBVERTEX) {
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            aUx[gg].lAmbda = 0;
          }
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          aUx[gg].lAmbda += inner_prod(data.getN(gg),data.getFieldData());

          if(aUx[gg].lAmbda!=aUx[gg].lAmbda) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"NaN value");
          }

        }

        if(type==MBVERTEX) {
          aUx[0].nodesWithoutLambda.resize(nb_dofs);
          for(int dd = 0;dd<nb_dofs;dd++) {
            if(data.getIndices()[dd]==-1) {
              aUx[0].nodesWithoutLambda[dd] = true;
            } else {
              aUx[0].nodesWithoutLambda[dd] = false;
            }
          }
        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Operator calculate \f$\overline{\lambda}\mathbf{C}^\mathsf{T}\f$
  */
  struct OpF: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpF(const string field_name,vector<AuxFunctions> &aux,DriverElementOrientation &orientation):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    aUx(aux),
    oRientation(orientation)
    {}

    ublas::vector<double> c;
    ublas::vector<double> nF;
    ublas::vector<int> rowIndices;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscErrorCode ierr;

      try {

        int nb_dofs = row_data.getFieldData().size();
        if(nb_dofs == 0) {
          PetscFunctionReturn(0);
        }
        int nb_gauss_pts = row_data.getN().size1();

        if(row_type == MBVERTEX) {
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            aUx[gg].nOrmal.resize(3,false);
            aUx[gg].nOrmal.clear();
            aUx[gg].calculateNormal();
          }
        }

        int eo = oRientation.elementOrientation;

        c.resize(nb_dofs,false);
        nF.resize(nb_dofs,false);
        nF.clear();

        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          double val = getGaussPts()(2,gg);
          noalias(c) = prod(aUx[gg].nOrmal,aUx[gg].N);
          noalias(nF) += val*eo*aUx[gg].lAmbda*c;

        }


        rowIndices.resize(nb_dofs,false);
        noalias(rowIndices) = row_data.getIndices();

        if(row_type == MBVERTEX) {
          int nb_dofs = aUx[0].nodesWithoutLambda.size();
          for(int dd = 0;dd<nb_dofs;dd++) {
            if(aUx[0].nodesWithoutLambda[dd]) {
              for(int jj = 0;jj<3;jj++) {
                rowIndices[dd*3+jj] = -1;
              }
            }
          }
        }

        int *indices_ptr = &rowIndices[0];

        ierr = VecSetOption(
          getFEMethod()->snes_f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE
        );  CHKERRQ(ierr);
        ierr = VecSetValues(
          getFEMethod()->snes_f,
          nb_dofs,
          indices_ptr,
          &nF[0],
          ADD_VALUES
        ); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpG: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpG(const string field_name,vector<AuxFunctions> &aux,DriverElementOrientation &orientation):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
    aUx(aux),
    oRientation(orientation)
    {}

    ublas::vector<double> g,dElta;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {

      PetscErrorCode ierr;

      try {

        int nb_dofs = row_data.getFieldData().size();
        if(nb_dofs == 0) {
          PetscFunctionReturn(0);
        }
        int nb_gauss_pts = row_data.getN().size1();

        int eo = oRientation.elementOrientation;

        dElta.resize(3);
        g.resize(nb_dofs,false);
        g.clear();

        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          noalias(dElta) =  aUx[gg].pOsition;
          for(int dd = 0;dd<3;dd++) {
            dElta[dd] -= getCoordsAtGaussPts()(gg,dd);
          }

          double val = getGaussPts()(2,gg);
          double r = inner_prod(aUx[gg].nOrmal,dElta);
          noalias(g) += val*eo*row_data.getN(gg)*r;

        }

        int *indices_ptr = &row_data.getIndices()[0];

        ierr = VecSetOption(
          getFEMethod()->snes_f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE
        );  CHKERRQ(ierr);

        ierr = VecSetValues(
          getFEMethod()->snes_f,
          nb_dofs,
          indices_ptr,
          &g[0],
          ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Operator calculating matrix \b C
  */
  struct OpC: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;
    bool assembleTranspose;

    OpC(
      const string lambda_field_name,
      const string positions_field_name,
      vector<AuxFunctions> &aux,
      DriverElementOrientation &orientation,
      bool assemble_transpose):
    FaceElementForcesAndSourcesCore::UserDataOperator(
      lambda_field_name,
      positions_field_name,
      UserDataOperator::OPROWCOL
    ),
    aUx(aux),
    oRientation(orientation),
    assembleTranspose(assemble_transpose) {
      sYmm = false;
    }


    ublas::vector<double> c;
    ublas::matrix<double> C;
    ublas::matrix<double> transC;
    ublas::vector<int> transRowIndices;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {

        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if(!nb_row || !nb_col) {
          PetscFunctionReturn(0);
        }

        int nb_gauss_pts = row_data.getN().size1();

        if(row_type == MBVERTEX) {
          for(int gg = 0;gg<nb_gauss_pts;gg++) {
            aUx[gg].nOrmal.resize(3,false);
            aUx[gg].nOrmal.clear();
            aUx[gg].calculateNormal();
          }
        }

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = aUx[gg].matrixN(gg,col_data); CHKERRQ(ierr);
        }

        int eo = oRientation.elementOrientation;

        c.resize(nb_col,false);
        C.resize(nb_row,nb_col,false);
        C.clear();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          noalias(c) = prod(aUx[gg].nOrmal,aUx[gg].N);
          double val = getGaussPts()(2,gg);
          noalias(C) += val*eo*outer_prod(row_data.getN(gg),c);

        }

        int *row_indices_ptr = &row_data.getIndices()[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        for(int n1 = 0; n1 != C.size1();n1++) {
          for(int n2 = 0; n2 != C.size1();n2++) {
            if(C(n1,n2)!=C(n1,n2)) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_INVALID_DATA,"NaN value");
            }
          }
        }

        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_row,row_indices_ptr,
          nb_col,col_indices_ptr,
          &C(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        if(assembleTranspose) {

          transC.resize(nb_col,nb_row);
          noalias(transC) = trans(C);

          transRowIndices.resize(nb_col,false);
          noalias(transRowIndices) = col_data.getIndices();
          int *trans_row_indices_ptr = &transRowIndices[0];
          if(row_type == MBVERTEX) {
            int nb_dofs = aUx[0].nodesWithoutLambda.size();
            for(int dd = 0;dd<nb_dofs;dd++) {
              if(aUx[0].nodesWithoutLambda[dd]) {
                for(int jj = 0;jj<3;jj++) {
                  transRowIndices[dd*3+jj] = -1;
                }
              }
            }
          }

          ierr = MatSetValues(
            getFEMethod()->snes_B,
            nb_col,trans_row_indices_ptr,
            nb_row,row_indices_ptr,
            &transC(0,0),ADD_VALUES
          ); CHKERRQ(ierr);

        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  /** \brief Operator calculating matrix \b B
  */
  struct OpB: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpB(const string field_name,vector<AuxFunctions> &aux,DriverElementOrientation &orientation):
    FaceElementForcesAndSourcesCore::UserDataOperator(
      field_name,field_name,UserDataOperator::OPROWCOL
    ),
    aUx(aux),
    oRientation(orientation)
    {
      sYmm = false;
    }

    ublas::matrix<double> spindXdKsi,spindXdEta;
    ublas::matrix<double> NdNormal;
    ublas::matrix<double> dNormal;

    ublas::matrix<double> B;
    ublas::vector<int> rowIndices;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {


        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if(!nb_row || !nb_col) {
          PetscFunctionReturn(0);
        }

        int nb_gauss_pts = row_data.getN().size1();
        int eo = oRientation.elementOrientation;

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = aUx[gg].matrixN(gg,row_data); CHKERRQ(ierr);
          ierr = aUx[gg].matrixB(gg,col_data); CHKERRQ(ierr);
        }

        spindXdKsi.resize(3,3,false);
        spindXdEta.resize(3,3,false);
        dNormal.resize(3,nb_col,false);
        NdNormal.resize(nb_row,nb_col,false);

        B.resize(nb_row,nb_col,false);
        B.clear();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          ierr = AuxFunctions::calcSpin(spindXdKsi,aUx[gg].dXdKsi); CHKERRQ(ierr);
          ierr = AuxFunctions::calcSpin(spindXdEta,aUx[gg].dXdEta); CHKERRQ(ierr);

          noalias(dNormal) = eo*(prod(spindXdKsi,aUx[gg].Beta)-prod(spindXdEta,aUx[gg].Bksi));
          noalias(NdNormal) = prod(trans(aUx[gg].N),dNormal);

          double val = getGaussPts()(2,gg);
          noalias(B) += val*aUx[gg].lAmbda*NdNormal;

        }


        rowIndices.resize(nb_row,false);
        noalias(rowIndices) = row_data.getIndices();
        if(row_type == MBVERTEX) {
          int nb_dofs = aUx[0].nodesWithoutLambda.size();
          for(int dd = 0;dd<nb_dofs;dd++) {
            if(aUx[0].nodesWithoutLambda[dd]) {
              for(int jj = 0;jj<3;jj++) {
                rowIndices[dd*3+jj] = -1;
              }
            }
          }
        }

        int *row_indices_ptr = &rowIndices[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_row,row_indices_ptr,
          nb_col,col_indices_ptr,
          &B(0,0),ADD_VALUES
        ); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Operator calculating matrix \b A
  */
  struct OpA: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    DriverElementOrientation &oRientation;

    OpA(
      const string lagrange_multipliers_field_name,
      const string field_name,
      vector<AuxFunctions> &aux,
      DriverElementOrientation &orientation
    ):
    FaceElementForcesAndSourcesCore::UserDataOperator(
      lagrange_multipliers_field_name,field_name,UserDataOperator::OPROWCOL
    ),
    aUx(aux),
    oRientation(orientation)
    {
      sYmm = false;
    }


    ublas::vector<double> XdNormal;
    ublas::matrix<double> spindXdKsi,spindXdEta;
    ublas::matrix<double> dNormal,NXdNormal;
    ublas::vector<double> dElta;

    ublas::matrix<double> A;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      try {

        int nb_row = row_data.getIndices().size();
        int nb_col = col_data.getIndices().size();
        if(!nb_row || !nb_col) {
          PetscFunctionReturn(0);
        }

        int nb_gauss_pts = row_data.getN().size1();
        int eo = oRientation.elementOrientation;

        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          ierr = aUx[gg].matrixB(gg,col_data); CHKERRQ(ierr);
        }

        XdNormal.resize(nb_col,false);
        spindXdKsi.resize(3,3,false);
        spindXdEta.resize(3,3,false);
        dNormal.resize(3,nb_col,false);
        NXdNormal.resize(nb_row,nb_col,false);
        dElta.resize(3);

        A.resize(nb_row,nb_col,false);
        A.clear();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          noalias(dElta) =  aUx[gg].pOsition;
          for(int dd = 0;dd<3;dd++) {
            dElta[dd] -= getCoordsAtGaussPts()(gg,dd);
          }

          ierr = AuxFunctions::calcSpin(spindXdKsi,aUx[gg].dXdKsi); CHKERRQ(ierr);
          ierr = AuxFunctions::calcSpin(spindXdEta,aUx[gg].dXdEta); CHKERRQ(ierr);

          noalias(dNormal) = eo*(prod(spindXdKsi,aUx[gg].Beta)-prod(spindXdEta,aUx[gg].Bksi));
          noalias(XdNormal) = prod(trans(dElta),dNormal);
          noalias(NXdNormal) = outer_prod(row_data.getN(gg),XdNormal);

          double val = getGaussPts()(2,gg);
          noalias(A) += val*NXdNormal;

        }

        int *row_indices_ptr = &row_data.getIndices()[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_row,row_indices_ptr,
          nb_col,col_indices_ptr,
          &A(0,0),ADD_VALUES
        ); CHKERRQ(ierr);


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /** \brief Driver function setting operators to calculate \b C matrix only
  */
  PetscErrorCode setOperatorsCOnly(
    const string lagrange_multipliers_field_name,
    const string material_field_name) {
    PetscFunctionBegin;

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent,crackFrontOrientation)
    );
    feLhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpC(
        lagrange_multipliers_field_name,material_field_name,cUrrent,crackFrontOrientation,false
      )
    );

    PetscFunctionReturn(0);
  }


  /** \brief Driver function setting operators to calculate nonlinear problems with sliding points on the surface
  */
  PetscErrorCode setOperatorsWithLinearGeometry(
    const string lagrange_multipliers_field_name,
    const string material_field_name,
    bool assemble_transpose,
    bool add_nonlinear_term
  ) {
    PetscFunctionBegin;

    // Adding operators to calculate the right hand side
    feRhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent,crackFrontOrientation)
    );
    feRhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feRhs.getOpPtrVector().push_back(
      new OpF(material_field_name,cUrrent,crackFrontOrientation)
    );
    feRhs.getOpPtrVector().push_back(
      new OpG(lagrange_multipliers_field_name,cUrrent,crackFrontOrientation)
    );

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent,crackFrontOrientation)
    );
    feLhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpC(
        lagrange_multipliers_field_name,material_field_name,cUrrent,crackFrontOrientation,assemble_transpose
      )
    );
    feLhs.getOpPtrVector().push_back(
      new OpB(material_field_name,cUrrent,crackFrontOrientation)
    );
    feLhs.getOpPtrVector().push_back(
      new OpA(lagrange_multipliers_field_name,material_field_name,cUrrent,crackFrontOrientation)
    );

    PetscFunctionReturn(0);
  }

};

#endif // __SURFACE_SLIDING_CONSTRAINS_HPP__
