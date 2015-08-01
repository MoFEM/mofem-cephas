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

#ifndef __GRIFFITH_FORCE_ELEMENT_HPP__
#define __GRIFFITH_FORCE_ELEMENT_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

struct GriffithForceElement {

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

  GriffithForceElement(FieldInterface &m_field):
    mField(m_field),
    feRhs(m_field),
    feLhs(m_field) {}

  struct CommonData {

    ublas::vector<double> griffithForce;
    ublas::matrix<double> tangentGriffithForce;
    vector<double*> tangentGriffithForceRowPtr;

    ublas::vector<double> penaltyForce;
    ublas::matrix<double> tangentPenaltyForce;
    vector<double*> tangentPenaltyForceRowPtr;

  };
  CommonData commonData;

  struct BlockData {
    double gc,penalty;
    Range frontEdges;
    Range frontNodes;
  };
  map<int,BlockData> blockData;

  struct OpJacobian: public FaceElementForcesAndSourcesCore::UserDataOperator {

    int tAg;
    BlockData blockData;
    CommonData &commonData;
    bool isTapeRecorded;

    OpJacobian(int tag,BlockData &data,CommonData &common_data):
    FaceElementForcesAndSourcesCore::UserDataOperator("MESH_NODE_POSITIONS",UserDataOperator::OPROW),
    tAg(tag),
    blockData(data),
    commonData(common_data),
    isTapeRecorded(false)
    {}

    /**

    \f[
    A = \| \mathbf{N} \| =
      \left\|
        \textrm{Spin}\left[
          \frac{\partial\mathbf{N}\mathbf{X}}{\partial\xi}
        \right]
        \frac{\partial\mathbf{N}\mathbf{X}}{\partial\eta}
        \right\|
    \f]

    */
    template<typename TYPE>
    struct AuxFunctions {

      ublas::matrix<double> N,NTN;
      ublas::matrix<double> Bksi;
      ublas::matrix<double> Beta;

      ublas::vector<TYPE> referenceCoords;

      ublas::vector<TYPE> referenceXdKsi;
      ublas::vector<TYPE> referenceXdEta;
      ublas::matrix<TYPE> referenceSpinKsi;
      ublas::matrix<TYPE> referenceSpinEta;
      ublas::vector<TYPE> referenceNormal;
      TYPE referenceArea;

      ublas::vector<TYPE> currentCoords;
      ublas::vector<TYPE> currentXdKsi;
      ublas::vector<TYPE> currentXdEta;
      ublas::matrix<TYPE> currentSpinKsi;
      ublas::matrix<TYPE> currentSpinEta;
      ublas::vector<TYPE> currentNormal;
      TYPE currentArea;

      ublas::matrix<TYPE> A;
      ublas::vector<TYPE> griffithForce;
      ublas::vector<TYPE> penaltyForce;
      TYPE dElta;

      PetscErrorCode sPin(ublas::matrix<TYPE> &spin,ublas::vector<TYPE> &vec) {
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

          N.resize(3,9,false);
          N.clear();
          for(int ii = 0;ii<3;ii++) {
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
        Bksi.resize(3,9);
        Bksi.clear();
        for(int ii = 0;ii<3;ii++) {
          for(int jj = 0;jj<3;jj++) {
            Bksi(jj,ii*3+jj) = data.getDiffN(gg)(ii,0);
          }
        }
        Beta.resize(3,9);
        Beta.clear();
        for(int ii = 0;ii<3;ii++) {
          for(int jj = 0;jj<3;jj++) {
            Beta(jj,ii*3+jj) = data.getDiffN(gg)(ii,1);
          }
        }
        PetscFunctionReturn(0);
      }

      PetscErrorCode dIffX() {
        PetscFunctionBegin;
        currentXdKsi.resize(3,false);
        currentXdEta.resize(3,false);
        noalias(currentXdKsi) = prod(Bksi,currentCoords);
        noalias(currentXdEta) = prod(Beta,currentCoords);
        PetscFunctionReturn(0);
      }

      PetscErrorCode nOrmal() {
        PetscFunctionBegin;
        PetscErrorCode ierr;
        ierr = sPin(currentSpinKsi,currentXdKsi); CHKERRQ(ierr);
        ierr = sPin(currentSpinEta,currentXdEta); CHKERRQ(ierr);
        currentNormal.resize(3,false);
        noalias(currentNormal) = 0.5*prod(currentSpinKsi,currentXdEta);
        currentArea = 0;
        for(int dd = 0;dd!=3;dd++) {
          currentArea += currentNormal[dd]*currentNormal[dd];
        }
        currentArea = sqrt(currentArea);
        PetscFunctionReturn(0);
      }

      PetscErrorCode matrixA() {
        PetscFunctionBegin;
        A.resize(3,9,false);
        noalias(A) = 0.5*(prod(currentSpinKsi,Beta) - prod(currentSpinEta,Bksi));
        PetscFunctionReturn(0);
      }

      PetscErrorCode calculateGrifthForce(
        double gc,double beta
      ) {
        PetscFunctionBegin;
        for(int dd = 0;dd!=9;dd++) {
          for(int ii = 0;ii!=3;ii++) {
            griffithForce[dd] += gc*beta*A(ii,dd)*currentNormal[ii]/currentArea;
          }
        }
        PetscFunctionReturn(0);
      }

      PetscErrorCode calculateReferenceNormal() {
        PetscFunctionBegin;

        PetscErrorCode ierr;

        referenceXdKsi.resize(3,false);
        referenceXdEta.resize(3,false);
        noalias(referenceXdKsi) = prod(Bksi,referenceCoords);
        noalias(referenceXdEta) = prod(Beta,referenceCoords);

        ierr = sPin(referenceSpinKsi,referenceXdKsi); CHKERRQ(ierr);
        //ierr = sPin(referenceSpinEta,referenceXdEta); CHKERRQ(ierr);
        referenceNormal.resize(3,false);
        noalias(referenceNormal) = 0.5*prod(referenceSpinKsi,referenceXdEta);
        referenceArea = 0;
        for(int dd = 0;dd!=3;dd++) {
          referenceArea += referenceNormal[dd]*referenceNormal[dd];
        }
        referenceArea = sqrt(referenceArea);

        PetscFunctionReturn(0);
      }

      PetscErrorCode calculatePenalty(double beta) {
        PetscFunctionBegin;

        dElta = (currentArea-referenceArea)/currentArea;
        dElta = -fmin(0,dElta);
        //dElta *= dElta*dElta;
        dElta *= currentArea;

        NTN.resize(9,9,false);
        noalias(NTN) = beta*prod(trans(N),N);
        noalias(penaltyForce) += dElta*prod(NTN,currentCoords-referenceCoords);

        PetscFunctionReturn(0);
      }

    };

    AuxFunctions<adouble> auxFun;

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(isTapeRecorded) {
        PetscFunctionReturn(0);
      }
      if(type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      try {

        PetscErrorCode ierr;
        int nb_dofs = data.getFieldData().size();
        int nb_gauss_pts = data.getN().size1();

        auxFun.referenceCoords.resize(9,false);
        auxFun.currentCoords.resize(9,false);

        trace_on(tAg);

        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.referenceCoords[dd] <<= getCoords()[dd];
        }
        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.currentCoords[dd] <<= data.getFieldData()[dd];
        }

        auxFun.griffithForce.resize(9,false);
        auxFun.griffithForce.clear();

        for(int gg = 0;gg!=nb_gauss_pts;gg++) {

          double val = getGaussPts()(2,gg)*0.5;

          ierr = auxFun.matrixB(gg,data); CHKERRQ(ierr);
          ierr = auxFun.dIffX(); CHKERRQ(ierr);
          ierr = auxFun.nOrmal(); CHKERRQ(ierr);
          ierr = auxFun.matrixA(); CHKERRQ(ierr);
          ierr = auxFun.calculateGrifthForce(blockData.gc,val); CHKERRQ(ierr);

          /*cerr << "gg: " << gg << endl;
          cerr << auxFun.Bksi << endl;
          cerr << auxFun.Beta << endl;
          cerr << auxFun.currentXdEta << endl;
          cerr << auxFun.currentXdEta << endl;
          cerr << "area " << auxFun.currentArea << endl;
          cerr << "normal " << auxFun.currentNormal << endl;
          cerr << "A " << auxFun.A << endl;
          cerr << "Griffith Force " << auxFun.griffithForce << endl;*/

        }

        //cerr << "Griffith Force " << auxFun.griffithForce << endl;

        commonData.griffithForce.resize(nb_dofs,false);
        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.griffithForce[dd] >>= commonData.griffithForce[dd];
        }

        trace_off();

        isTapeRecorded = true;


      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct AuxOp {

    int tAg;
    BlockData &blockData;
    CommonData &commonData;
    AuxOp(int tag,BlockData &block_data,CommonData &common_data):
    tAg(tag),
    blockData(block_data),
    commonData(common_data) {
    };

    ublas::vector<int> rowIndices;
    ublas::vector<double> activeVariables;

    PetscErrorCode setIndices(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      int nb_dofs = data.getIndices().size();
      rowIndices.resize(nb_dofs,false);
      noalias(rowIndices) = data.getIndices();
      ublas::vector<const FEDofMoFEMEntity*>& dofs = data.getFieldDofs();
      ublas::vector<const FEDofMoFEMEntity*>::iterator dit,hi_dit;
      dit = dofs.begin();
      hi_dit = dofs.end();
      for(int ii = 0;dit!=hi_dit;dit++,ii++) {
        if(blockData.frontNodes.find((*dit)->get_ent())==blockData.frontNodes.end()) {
          rowIndices[ii] = -1;
        }
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode setVariables(
      FaceElementForcesAndSourcesCore::UserDataOperator *fe_ptr,
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    ) {
      PetscFunctionBegin;
      int nb_dofs = data.getIndices().size();
      activeVariables.resize(18);
      for(int dd = 0;dd!=nb_dofs;dd++) {
        activeVariables[dd] = fe_ptr->getCoords()[dd];
      }
      for(int dd = 0;dd!=nb_dofs;dd++) {
        activeVariables[9+dd] = data.getFieldData()[dd];
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpRhs: public FaceElementForcesAndSourcesCore::UserDataOperator,AuxOp {

    OpRhs(int tag,BlockData &block_data,CommonData &common_data):
    FaceElementForcesAndSourcesCore::UserDataOperator("MESH_NODE_POSITIONS",UserDataOperator::OPROW),
    AuxOp(tag,block_data,common_data) {
    }

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      if(type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      PetscErrorCode ierr;

      ierr = setIndices(side,type,data); CHKERRQ(ierr);
      ierr = setVariables(this,side,type,data); CHKERRQ(ierr);

      int nb_dofs = data.getIndices().size();
      commonData.griffithForce.resize(9,false);
      int r;
      //play recorder for values
      r = ::function(tAg,nb_dofs,18,&activeVariables[0],&commonData.griffithForce[0]);
      if(r<3) {
        SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
      }

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        &rowIndices[0],
        &commonData.griffithForce[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpLhs: public FaceElementForcesAndSourcesCore::UserDataOperator,AuxOp {

    OpLhs(int tag,BlockData &block_data,CommonData &common_data):
    FaceElementForcesAndSourcesCore::UserDataOperator(
      "MESH_NODE_POSITIONS",
      "MESH_NODE_POSITIONS",
      UserDataOperator::OPROWCOL
    ),
    AuxOp(tag,block_data,common_data) {
    }

    ublas::matrix<double> k;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      if(row_type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      PetscErrorCode ierr;

      ierr = setIndices(row_side,row_type,row_data);
      ierr = setVariables(this,row_side,row_type,row_data);

      commonData.tangentGriffithForce.resize(9,18,false);
      commonData.tangentGriffithForceRowPtr.resize(9);
      for(int dd = 0;dd<9;dd++) {
        commonData.tangentGriffithForceRowPtr[dd] = &commonData.tangentGriffithForce(dd,0);
      }

      int row_nb_dofs = row_data.getIndices().size();

      int r;
      //play recorder for jacobian
      r = jacobian(
        tAg,row_nb_dofs,18,
        &activeVariables[0],
        &commonData.tangentGriffithForceRowPtr[0]
      );
      if(r<3) { // function is locally analytic
        SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
      }
      k.resize(9,9,false);
      for(int dd1 = 0;dd1<9;dd1++) {
        for(int dd2 = 0;dd2<9;dd2++) {
          k(dd1,dd2) = commonData.tangentGriffithForce(dd1,9+dd2);
        }
      }

      int col_nb_dofs = row_data.getIndices().size();

      ierr = MatSetValues(
        getFEMethod()->snes_B,
        row_nb_dofs,&rowIndices[0],
        col_nb_dofs,&col_data.getIndices()[0],
        &k(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };


  struct OpJacobianPenalty: public OpJacobian {

    OpJacobianPenalty(int tag,BlockData &data,CommonData &common_data):
    OpJacobian(tag,data,common_data) {
    }

  PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(isTapeRecorded) {
        PetscFunctionReturn(0);
      }
      if(type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      try {

        PetscErrorCode ierr;
        int nb_dofs = data.getFieldData().size();
        int nb_gauss_pts = data.getN().size1();

        auxFun.referenceCoords.resize(9,false);
        auxFun.currentCoords.resize(9,false);

        trace_on(tAg);

        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.referenceCoords[dd] <<= getCoords()[dd];
        }
        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.currentCoords[dd] <<= data.getFieldData()[dd];
        }

        auxFun.penaltyForce.resize(9,false);
        auxFun.penaltyForce.clear();

        for(int gg = 0;gg!=nb_gauss_pts;gg++) {

          double val = blockData.penalty*getGaussPts()(2,gg)*0.5;
          ierr = auxFun.matrixB(gg,data); CHKERRQ(ierr);
          ierr = auxFun.matrixN(gg,data); CHKERRQ(ierr);
          ierr = auxFun.dIffX(); CHKERRQ(ierr);
          ierr = auxFun.nOrmal(); CHKERRQ(ierr);
          ierr = auxFun.matrixA(); CHKERRQ(ierr);
          ierr = auxFun.calculateReferenceNormal(); CHKERRQ(ierr);
          ierr = auxFun.calculatePenalty(val);

        }

        commonData.penaltyForce.resize(nb_dofs,false);
        for(int dd = 0;dd!=nb_dofs;dd++) {
          auxFun.penaltyForce[dd] >>= commonData.penaltyForce[dd];
        }

        trace_off();

        isTapeRecorded = true;

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpPenaltyRhs: public OpRhs {

    OpPenaltyRhs(int tag,BlockData &block_data,CommonData &common_data):
    OpRhs(tag,block_data,common_data) {
    }

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      if(type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      PetscErrorCode ierr;

      ierr = setIndices(side,type,data); CHKERRQ(ierr);
      ierr = setVariables(this,side,type,data); CHKERRQ(ierr);

      int nb_dofs = data.getIndices().size();
      commonData.penaltyForce.resize(9,false);
      int r;
      //play recorder for values
      r = ::function(tAg,nb_dofs,18,&activeVariables[0],&commonData.penaltyForce[0]);
      if(r<1) { // function is locally analytic
        SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
      }

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        &rowIndices[0],
        &commonData.penaltyForce[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpPenaltyLhs: public OpLhs {

    OpPenaltyLhs(int tag,BlockData &block_data,CommonData &common_data):
    OpLhs(tag,block_data,common_data) {

    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      if(row_type != MBVERTEX) {
        PetscFunctionReturn(0);
      }

      PetscErrorCode ierr;

      ierr = setIndices(row_side,row_type,row_data);
      ierr = setVariables(this,row_side,row_type,row_data);

      commonData.tangentPenaltyForce.resize(9,18,false);
      commonData.tangentPenaltyForceRowPtr.resize(9);
      for(int dd = 0;dd<9;dd++) {
        commonData.tangentPenaltyForceRowPtr[dd] = &commonData.tangentPenaltyForce(dd,0);
      }

      int row_nb_dofs = row_data.getIndices().size();

      int r;
      //play recorder for jacobian
      r = jacobian(
        tAg,row_nb_dofs,18,
        &activeVariables[0],
        &commonData.tangentPenaltyForceRowPtr[0]
      );
      if(r<1) { // function is locally analytic
        SETERRQ1(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"ADOL-C function evaluation with error r = %d",r);
      }
      k.resize(9,9,false);
      for(int dd1 = 0;dd1<9;dd1++) {
        for(int dd2 = 0;dd2<9;dd2++) {
          k(dd1,dd2) = commonData.tangentPenaltyForce(dd1,9+dd2);
        }
      }

      int col_nb_dofs = row_data.getIndices().size();

      ierr = MatSetValues(
        getFEMethod()->snes_B,
        row_nb_dofs,&rowIndices[0],
        col_nb_dofs,&col_data.getIndices()[0],
        &k(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

};

#endif //__GRIFFITH_FORCE_ELEMENT_HPP__
