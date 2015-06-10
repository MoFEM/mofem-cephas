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

/** \brief Surface sliding constrains.

  Displacements on the body tangential direction are friction less.
  Displacements in direction normal are restricted.

  Constrains are derived form equation:
  \f[
  \frac{V}{A} = C
  \f]

  \f[
  \int_\Omega \textrm{d}V = C \int_\Gamma \textrm{d}S
  \f]

  \f[
  \frac{1}{3}
  \int_\Omega \textrm{div}[\mathbf{X}] \textrm{d}\Omega
  =
  C \int_\Gamma  \textrm{d}S
  \f]

  \f[
  \int_\Gamma
  \mathbf{X}\cdot \frac{\mathbf{N}}{\|\mathbf{N}\|}
  \textrm{d}\Gamma
  =
  3C \int_\Gamma  \textrm{d}S
  \f]

  Drooping integrals on both sides, and linearizing equation, we get
  \f[
  \frac{\mathbf{n}}{\|\mathbf{n}\|} \cdot \delta \mathbf{X}
  =
  3c - \frac{\mathbf{n}}{\|\mathbf{n}\|}\cdot \mathbf{x}
  \f]
  where \f$\delta \mathbf{X}\f$ is displacement sub-inctrement.

  Above equation is a constrain which need to be enforced.

  \f[
  \mathcal{r} =
  \frac{\mathbf{N}}{\|\mathbf{N}\|}\cdot \mathbf{X}
  -
  \frac{\mathbf{N_0}}{\|\mathbf{N_0}\|}\cdot \mathbf{X}_0
  =
  \mathbf{C}\overline{\mathbf{X}}-\mathbf{C}_0\overline{\mathbf{X}_0}
  \f]

  \f[
  \int_\Gamma \mathbf{N}^\mathsf{T}_\lambda
  \left(
  \frac{\mathbf{N}}{\|\mathbf{N}\|}\cdot \mathbf{N}_\mathbf{X}
  \right)
  \textrm{d}\Gamma
  \cdot
  \overline{\mathbf{X}}
  =
  \int_\Gamma \mathbf{N}^\mathsf{T}_\lambda r \textrm{d}\Gamma
  \f]

  \f[
  \mathbf{C}\overline{\mathbf{X}} = \overline{\mathbf{r}}
  \f]

  What result in additional terms in global system of equations
  \f[
  \left[
    \begin{array}{cc}
        \mathbf{K} + \lambda\frac{\partial \mathbf{C}}{\partial\overline{\mathbf{X}}} & \mathbf{C}^\mathsf{T} \\
        \mathbf{C} & 0
    \end{array}
  \right]
  \left\{
    \begin{array}{c}
      \delta \overline{\mathbf{X}} \\
      \delta \lambda
    \end{array}
  \right\}=
  \left[
    \begin{array}{c}
      \mathbf{f} - \mathbf{C}^\mathsf{T}\delta\lambda \\
      \overline{\mathbf{r}}
    \end{array}
  \right]
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
    int getRule(int order) { return order-1; };

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

  SurfaceSlidingConstrains(FieldInterface &m_field):
  mField(m_field),
  feRhs(m_field),
  feLhs(m_field) {}

  struct AuxFunctions {

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

  vector<AuxFunctions> cUrrent,rEference;

  struct OpPositions: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;

    OpPositions(const string field_name,vector<AuxFunctions> &aux):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    aUx(aux) {
    }

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

        }

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpF: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;

    OpF(const string field_name,vector<AuxFunctions> &aux):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROWCOL),
    aUx(aux) {}

    ublas::vector<double> c;
    ublas::vector<double> nf;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscErrorCode ierr;

      int nb_dofs = row_data.getFieldData().size();
      if(nb_dofs == 0) {
        PetscFunctionReturn(0);
      }
      int nb_gauss_pts = row_data.getN().size1();

      c.resize(nb_dofs,false);
      nf.resize(nb_dofs,false);
      nf.clear();

      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        double val = getGaussPts()(2,gg)*getArea();
        noalias(c) = prod(aUx[gg].nOrmal,aUx[gg].N);

        noalias(nf) += val*aUx[gg].lAmbda*c;

      }

      int *indices_ptr = &row_data.getIndices()[0];

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        indices_ptr,
        &nf[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpG: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    double aLpha;

    OpG(const string field_name,vector<AuxFunctions> &aux,double alpha):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROWCOL),
    aUx(aux),
    aLpha(alpha) {}

    ublas::vector<double> g;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {

      PetscErrorCode ierr;

      int nb_dofs = row_data.getFieldData().size();
      if(nb_dofs == 0) {
        PetscFunctionReturn(0);
      }
      int nb_gauss_pts = row_data.getN().size1();

      g.resize(nb_dofs,false);
      g.clear();

      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        double val = getGaussPts()(2,gg);
        double r = inner_prod(aUx[gg].nOrmal,aUx[gg].pOsition);
        noalias(g) += val*aLpha*row_data.getN(gg)*r;

      }

      int *indices_ptr = &row_data.getIndices()[0];

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        indices_ptr,
        &g[0],
        ADD_VALUES
      ); CHKERRQ(ierr);


      PetscFunctionReturn(0);
    }

  };

  struct OpGFromMeshCoords: FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    double aLpha;

    OpGFromMeshCoords(const string field_name,vector<AuxFunctions> &aux,double alpha):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROWCOL),
    aUx(aux),
    aLpha(alpha) {}

    ublas::vector<double> g;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {

      PetscErrorCode ierr;

      int nb_dofs = row_data.getFieldData().size();
      if(nb_dofs == 0) {
        PetscFunctionReturn(0);
      }
      int nb_gauss_pts = row_data.getN().size1();

      g.resize(nb_dofs,false);
      g.clear();

      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        double val = getGaussPts()(2,gg);
        double r;
        for(int dd = 0;dd!=3;dd++) {
          r = 0.5*getNormal()[dd]*getCoordsAtGaussPts()(gg,dd);
        }
        noalias(g) += val*aLpha*row_data.getN(gg)*r;

      }

      int *indices_ptr = &row_data.getIndices()[0];

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        indices_ptr,
        &g[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpC: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;
    bool assembleTranspose;

    OpC(
      const string lambda_field_name,
      const string positions_field_name,
      vector<AuxFunctions> &aux,
      bool assemble_transpose):
    FaceElementForcesAndSourcesCore::UserDataOperator(
      lambda_field_name,
      positions_field_name,
      UserDataOperator::OPROWCOL
    ),
    aUx(aux),
    assembleTranspose(assemble_transpose) {
      sYmm = false;
    }


    ublas::vector<double> c;
    ublas::matrix<double> C;
    ublas::matrix<double> transC;

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


        c.resize(nb_col,false);
        C.resize(nb_row,nb_col,false);
        C.clear();
        for(int gg = 0;gg<nb_gauss_pts;gg++) {

          //ierr = aUx[gg].matrixN(gg,col_data); CHKERRQ(ierr);
          noalias(c) = prod(aUx[gg].nOrmal,aUx[gg].N);
          double val = getGaussPts()(2,gg);
          noalias(C) += val*outer_prod(row_data.getN(gg),c);

        }

        int *row_indices_ptr = &row_data.getIndices()[0];
        int *col_indices_ptr = &col_data.getIndices()[0];

        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_row,row_indices_ptr,
          nb_col,col_indices_ptr,
          &C(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

        if(assembleTranspose) {

          transC.resize(nb_col,nb_row);
          noalias(transC) = trans(C);

          ierr = MatSetValues(
            getFEMethod()->snes_B,
            nb_col,col_indices_ptr,
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

  struct OpDiffC: public FaceElementForcesAndSourcesCore::UserDataOperator {

    vector<AuxFunctions> &aUx;

    OpDiffC(const string field_name,vector<AuxFunctions> &aux):
    FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROWCOL),
    aUx(aux) {
      sYmm = false;
    }

    ublas::matrix<double> spindXdKsi,spindXdEta;
    ublas::matrix<double> NdNormal;
    ublas::matrix<double> dNormal;

    ublas::matrix<double> diffC;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      PetscErrorCode ierr;

      int nb_dofs = row_data.getFieldData().size();
      if(nb_dofs == 0) {
        PetscFunctionReturn(0);
      }
      int nb_gauss_pts = row_data.getN().size1();


      spindXdKsi.resize(3,3,false);
      spindXdEta.resize(3,3,false);
      dNormal.resize(nb_dofs,nb_dofs,false);
      NdNormal.resize(nb_dofs,false);

      diffC.resize(nb_dofs,nb_dofs,false);
      diffC.clear();
      for(int gg = 0;gg<nb_gauss_pts;gg++) {

        double val = getGaussPts()(2,gg);
        ierr = AuxFunctions::calcSpin(spindXdKsi,aUx[gg].dXdKsi); CHKERRQ(ierr);
        ierr = AuxFunctions::calcSpin(spindXdEta,aUx[gg].dXdEta); CHKERRQ(ierr);

        noalias(dNormal) = prod(spindXdKsi,aUx[gg].Bksi)-prod(spindXdEta,aUx[gg].Beta);
        noalias(NdNormal) = prod(aUx[gg].N,dNormal);
        noalias(diffC) += val*aUx[gg].lAmbda*NdNormal;

      }

      int nb_row = row_data.getIndices().size();
      int nb_col = col_data.getIndices().size();

      int *row_indices_ptr = &row_data.getIndices()[0];
      int *col_indices_ptr = &col_data.getIndices()[0];

      ierr = MatSetValues(
        getFEMethod()->snes_B,
        nb_row,row_indices_ptr,
        nb_col,col_indices_ptr,
        &diffC(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode setOperatorsCOnly(
    const string lagrange_multipliers_field_name,
    const string material_field_name) {
    PetscFunctionBegin;

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpC(lagrange_multipliers_field_name,material_field_name,cUrrent,false)
    );

    PetscFunctionReturn(0);
  }

  PetscErrorCode setOperatorsWithLinearGeometry(
    const string lagrange_multipliers_field_name,
    const string material_field_name,
    bool assemble_transpose,
    bool add_nonlinear_term
  ) {
    PetscFunctionBegin;

    // Adding operators to calculate the right hand side
    feRhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent)
    );
    feRhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feRhs.getOpPtrVector().push_back(
      new OpF(material_field_name,cUrrent)
    );
    feRhs.getOpPtrVector().push_back(
      new OpG(lagrange_multipliers_field_name,cUrrent,+1)
    );
    feRhs.getOpPtrVector().push_back(
      new OpGFromMeshCoords(lagrange_multipliers_field_name,cUrrent,-1)
    );

    // Adding operators to calculate the left hand side
    feLhs.getOpPtrVector().push_back(
      new OpPositions(material_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpLambda(lagrange_multipliers_field_name,cUrrent)
    );
    feLhs.getOpPtrVector().push_back(
      new OpC(lagrange_multipliers_field_name,material_field_name,cUrrent,assemble_transpose)
    );
    feLhs.getOpPtrVector().push_back(
      new OpDiffC(material_field_name,cUrrent)
    );

    PetscFunctionReturn(0);
  }

};

#endif // __SURFACE_SLIDING_CONSTRAINS_HPP__
