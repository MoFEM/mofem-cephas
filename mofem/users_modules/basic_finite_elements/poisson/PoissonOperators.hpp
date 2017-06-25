/**
 * \file PoissonOperators.cpp
 * \example PoissonOperators.hpp
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

#ifndef __POISSONOPERATORS_HPP__
#define __POISSONOPERATORS_HPP__

namespace PoissonOperators {

  struct OpGradGrad: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    OpGradGrad():
    VolumeElementForcesAndSourcesCore::UserDataOperator("U","U",OPROWCOL,true) {
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      nbRows = row_data.getIndices().size();
      if(!nbRows) PetscFunctionReturn(0);
      nbCols = col_data.getIndices().size();
      if(!nbCols) PetscFunctionReturn(0);
      nbIntegrationPts = getGaussPts().size2();
      if(
        row_side==col_side&&
        row_type==col_type
      ) {
        isDiag = true;
      } else {
        isDiag = false;
      }
      ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    PetscErrorCode ierr;

    int nbRows;
    int nbCols;
    int nbIntegrationPts;
    bool isDiag;

    FTensor::Index<'i',3> i;
    MatrixDouble locMat;

    inline PetscErrorCode iNtegrte(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      locMat.resize(nbRows,nbCols,false);
      locMat.clear();
      double vol = getVolume();
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      FTensor::Tensor1<double*,3> t_row_grad = row_data.getFTensor1DiffN<3>();
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        const double alpha = t_w*vol;
        FTensor::Tensor0<double*> a(&*locMat.data().begin());
        for(int rr = 0;rr!=nbRows;rr++) {
          FTensor::Tensor1<double*,3> t_col_grad = col_data.getFTensor1DiffN<3>(gg,0);
          for(int cc = 0;cc!=nbCols;cc++) {
            a += alpha*(t_row_grad(i)*t_col_grad(i));
            ++t_col_grad;
            ++a;
          }
          ++t_row_grad;
        }
        ++t_w;
      }
      PetscFunctionReturn(0);
    }

    inline PetscErrorCode aSsemble(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      const int* row_indices = &*row_data.getIndices().data().begin();
      const int* col_indices = &*col_data.getIndices().data().begin();
      ierr = MatSetValues(
        getFEMethod()->ksp_B,
        nbRows,row_indices,
        nbCols,col_indices,
        &*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      if(!isDiag) {
        locMat = trans(locMat);
        ierr = MatSetValues(
          getFEMethod()->ksp_B,
          nbCols,col_indices,
          nbRows,row_indices,
          &*locMat.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  template<typename OPBASE>
  struct OpBaseRhs: public OPBASE {

    OpBaseRhs(const std::string field_name):
    OPBASE(field_name,OPBASE::OPROW) {
    }

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;
      nbRows = row_data.getIndices().size();
      if(!nbRows) PetscFunctionReturn(0);
      nbIntegrationPts = OPBASE::getGaussPts().size2();
      ierr = iNtegrte(row_data); CHKERRQ(ierr);
      ierr = aSsemble(row_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) = 0;
    virtual PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) = 0;

  protected:

    PetscErrorCode ierr;
    int nbRows;
    int nbIntegrationPts;

  };

  struct OpVF: public OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> FSource;

    OpVF(FSource f_source):
    OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator>("U"),
    fSource(f_source) {
    }

  private:

    PetscErrorCode ierr;
    FTensor::Number<0> NX;
    FTensor::Number<1> NY;
    FTensor::Number<2> NZ;
    FSource fSource;

    Vec F;
    VectorDouble locVec;

    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      // cerr << nbRows << " " << nbIntegrationPts << endl;
      locVec.resize(nbRows,false);
      locVec.clear();
      double vol = getVolume();
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      FTensor::Tensor0<double*> t_v = data.getFTensor0N();
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        double alpha = vol*t_w*fSource(t_coords(NX),t_coords(NY),t_coords(NZ));
        FTensor::Tensor0<double*> t_a(&*locVec.data().begin());
        for(int rr = 0;rr!=nbRows;rr++) {
          t_a -= alpha*t_v;
          ++t_a;
          ++t_v;
        }
        ++t_w;
        ++t_coords;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      const int* indices = &*data.getIndices().data().begin();
      const double* vals = &*locVec.data().begin();
      ierr = VecSetValues(
        getFEMethod()->ksp_f,nbRows,indices,vals,ADD_VALUES
      ); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  struct OpLUHdiv: public FaceElementForcesAndSourcesCore::UserDataOperator {

    OpLUHdiv(const bool assemble_transpose):
    FaceElementForcesAndSourcesCore::UserDataOperator("L","U",OPROWCOL,false),
    assembleTraspose(assemble_transpose) {
    }

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      nbRows = row_data.getIndices().size();
      if(!nbRows) PetscFunctionReturn(0);
      nbCols = col_data.getIndices().size();
      if(!nbCols) PetscFunctionReturn(0);
      nbIntegrationPts = getGaussPts().size2();
      ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
      ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  private:

    PetscErrorCode ierr;

    int nbRows;
    int nbCols;
    int nbIntegrationPts;
    const bool assembleTraspose;

    FTensor::Index<'i',3> i;
    MatrixDouble locMat;

    inline PetscErrorCode iNtegrte(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      locMat.resize(nbRows,nbCols,false);
      locMat.clear();
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      FTensor::Tensor1<double*,3> t_normal = getTensor1Normal();
      FTensor::Tensor1<double*,3> t_row = row_data.getFTensor1HdivN<3>();
      // cerr << nbRows << " :::  " << nbIntegrationPts << " ::: " << nbCols << endl;
      // cerr << row_data.getHdivN() << endl;
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        const double alpha = 0.5*t_w;
        FTensor::Tensor0<double*> c(&*locMat.data().begin());
        for(int rr = 0;rr!=nbRows;rr++) {
          FTensor::Tensor0<double*> t_col = col_data.getFTensor0N(gg,0);
          for(int cc = 0;cc!=nbCols;cc++) {
            c += alpha*t_normal(i)*t_row(i)*t_col;
            ++t_col;
            ++c;
          }
          ++t_row;
        }
        ++t_w;
      }
      PetscFunctionReturn(0);
    }

    inline PetscErrorCode aSsemble(
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;
      const int* row_indices = &*row_data.getIndices().data().begin();
      const int* col_indices = &*col_data.getIndices().data().begin();
      ierr = MatSetValues(
        getFEMethod()->ksp_B,
        nbRows,row_indices,
        nbCols,col_indices,
        &*locMat.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      // cerr << locMat << endl;
      if(assembleTraspose) {
        locMat = trans(locMat);
        ierr = MatSetValues(
          getFEMethod()->ksp_B,
          nbCols,col_indices,
          nbRows,row_indices,
          &*locMat.data().begin(),ADD_VALUES
        ); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpLgHdiv: public OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> FVal;

    OpLgHdiv(FVal f_value):
    OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator>("L"),
    fValue(f_value) {
    }

  private:

    FTensor::Number<0> NX;
    FTensor::Number<1> NY;
    FTensor::Number<2> NZ;
    FTensor::Index<'i',3> i;
    FVal fValue;

    VectorDouble locVec;

    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      locVec.resize(nbRows,false);
      locVec.clear();
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      FTensor::Tensor1<double*,3> t_l = data.getFTensor1HdivN<3>();
      FTensor::Tensor1<double*,3> t_normal = getTensor1Normal();
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        double alpha = 0.5*t_w*fValue(t_coords(NX),t_coords(NY),t_coords(NZ));
        FTensor::Tensor0<double*> t_a(&*locVec.data().begin());
        for(int rr = 0;rr!=nbRows;rr++) {
          t_a += alpha*t_l(i)*t_normal(i);
          ++t_a;
          ++t_l;
        }
        ++t_w;
        ++t_coords;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      const int* indices = &*data.getIndices().data().begin();
      const double* vals = &*locVec.data().begin();
      ierr = VecSetValues(
        getFEMethod()->ksp_f,nbRows,indices,&*locVec.data().begin(),ADD_VALUES
      ); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

  };

  struct OpErrorL2: public OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator> {

    typedef boost::function<double (const double,const double,const double)> FVal;

    OpErrorL2(
      FVal f_value,boost::shared_ptr<VectorDouble>& field_vals
    ):
    OpBaseRhs<VolumeElementForcesAndSourcesCore::UserDataOperator>("ERROR"),
    fValue(f_value),
    fieldVals(field_vals) {
    }

  private:

    FTensor::Number<0> NX;
    FTensor::Number<1> NY;
    FTensor::Number<2> NZ;
    FTensor::Index<'i',3> i;
    FVal fValue;

    boost::shared_ptr<VectorDouble>& fieldVals;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;
      nbRows = row_data.getFieldData().size();
      if(!nbRows) PetscFunctionReturn(0);
      nbIntegrationPts = getGaussPts().size2();
      ierr = iNtegrte(row_data); CHKERRQ(ierr);
      ierr = aSsemble(row_data); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      const double vol = getVolume();
      FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
      FTensor::Tensor0<double*> t_e = data.getFTensor0N();
      FTensor::Tensor0<double*> t_u = getTensor0FormData(*fieldVals);
      FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
      FTensor::Tensor0<double*> t_a = data.getFTensor0FieldData();
      for(int rr = 0;rr!=nbRows;rr++) {
        t_a = 0;
        ++t_a;
      }
      for(int gg = 0;gg!=nbIntegrationPts;gg++) {
        double alpha = vol*t_w;
        FTensor::Tensor0<double*> t_a = data.getFTensor0FieldData();
        for(int rr = 0;rr!=nbRows;rr++) {
          t_a += alpha*t_e*(t_u-fValue(t_coords(NX),t_coords(NY),t_coords(NZ)));
          ++t_a;
          ++t_e;
        }
        ++t_w;
        ++t_u;
        ++t_coords;
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      FTensor::Tensor0<double*> t_a = data.getFTensor0FieldData();
      for(int rr = 0;rr!=nbRows;rr++) {
        data.getFieldDofs()[rr]->getFieldData() = t_a;
        ++t_a;
      }
      PetscFunctionReturn(0);
    }


  };


  // struct OpLU: public FaceElementForcesAndSourcesCore::UserDataOperator {
  //
  //   OpLU(const bool assemble_transpose):
  //   FaceElementForcesAndSourcesCore::UserDataOperator("L","U",OPROWCOL,false),
  //   assembleTraspose(assemble_transpose) {
  //   }
  //
  //   PetscErrorCode doWork(
  //     int row_side,int col_side,
  //     EntityType row_type,EntityType col_type,
  //     DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
  //   ) {
  //     PetscFunctionBegin;
  //     nbRows = row_data.getIndices().size();
  //     if(!nbRows) PetscFunctionReturn(0);
  //     nbCols = col_data.getIndices().size();
  //     if(!nbCols) PetscFunctionReturn(0);
  //     nbIntegrationPts = getGaussPts().size2();
  //     ierr = iNtegrte(row_data,col_data); CHKERRQ(ierr);
  //     ierr = aSsemble(row_data,col_data); CHKERRQ(ierr);
  //     PetscFunctionReturn(0);
  //   }
  //
  // private:
  //
  //   PetscErrorCode ierr;
  //
  //   int nbRows;
  //   int nbCols;
  //   int nbIntegrationPts;
  //   const bool assembleTraspose;
  //
  //   MatrixDouble locMat;
  //
  //   inline PetscErrorCode iNtegrte(
  //     DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
  //   ) {
  //     PetscFunctionBegin;
  //     locMat.resize(nbRows,nbCols,false);
  //     locMat.clear();
  //     const double area = getArea();
  //     FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
  //     FTensor::Tensor0<double*> t_row = row_data.getFTensor0N();
  //     for(int gg = 0;gg!=nbIntegrationPts;gg++) {
  //       const double alpha = area*t_w;
  //       FTensor::Tensor0<double*> c(&*locMat.data().begin());
  //       for(int rr = 0;rr!=nbRows;rr++) {
  //         FTensor::Tensor0<double*> t_col = col_data.getFTensor0N(gg,0);
  //         for(int cc = 0;cc!=nbCols;cc++) {
  //           c += alpha*t_row*t_col;
  //           ++t_col;
  //           ++c;
  //         }
  //         ++t_row;
  //       }
  //       ++t_w;
  //     }
  //     PetscFunctionReturn(0);
  //   }
  //
  //   inline PetscErrorCode aSsemble(
  //     DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
  //   ) {
  //     PetscFunctionBegin;
  //     const int* row_indices = &*row_data.getIndices().data().begin();
  //     const int* col_indices = &*col_data.getIndices().data().begin();
  //     ierr = MatSetValues(
  //       getFEMethod()->ksp_B,
  //       nbRows,row_indices,
  //       nbCols,col_indices,
  //       &*locMat.data().begin(),ADD_VALUES
  //     ); CHKERRQ(ierr);
  //     // cerr << locMat << endl;
  //     if(assembleTraspose) {
  //       locMat = trans(locMat);
  //       ierr = MatSetValues(
  //         getFEMethod()->ksp_B,
  //         nbCols,col_indices,
  //         nbRows,row_indices,
  //         &*locMat.data().begin(),ADD_VALUES
  //       ); CHKERRQ(ierr);
  //     }
  //     PetscFunctionReturn(0);
  //   }
  //
  // };

  // struct OpLg: public OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator> {
  //
  //   typedef boost::function<double (const double,const double,const double)> FVal;
  //
  //   OpLg(FVal f_value):
  //   OpBaseRhs<FaceElementForcesAndSourcesCore::UserDataOperator>("L"),
  //   fValue(f_value) {
  //   }
  //
  // private:
  //
  //   FTensor::Number<0> NX;
  //   FTensor::Number<1> NY;
  //   FTensor::Number<2> NZ;
  //   FVal fValue;
  //
  //   VectorDouble locVec;
  //
  //   PetscErrorCode iNtegrte(DataForcesAndSurcesCore::EntData &data) {
  //     PetscFunctionBegin;
  //     locVec.resize(nbRows,false);
  //     locVec.clear();
  //     const double area = getArea();
  //     FTensor::Tensor0<double*> t_w = getFTensor0IntegrationWeight();
  //     FTensor::Tensor0<double*> t_l = data.getFTensor0N();
  //     FTensor::Tensor1<double*,3> t_coords = getTensor1CoordsAtGaussPts();
  //     for(int gg = 0;gg!=nbIntegrationPts;gg++) {
  //       double alpha = area*t_w*fValue(t_coords(NX),t_coords(NY),t_coords(NZ));
  //       FTensor::Tensor0<double*> t_a(&*locVec.data().begin());
  //       for(int rr = 0;rr!=nbRows;rr++) {
  //         t_a += alpha*t_l;
  //         ++t_a;
  //         ++t_l;
  //       }
  //       ++t_w;
  //       ++t_coords;
  //     }
  //     PetscFunctionReturn(0);
  //   }
  //
  //   PetscErrorCode aSsemble(DataForcesAndSurcesCore::EntData &data) {
  //     PetscFunctionBegin;
  //     const int* indices = &*data.getIndices().data().begin();
  //     const double* vals = &*locVec.data().begin();
  //     ierr = VecSetValues(
  //       getFEMethod()->ksp_f,nbRows,indices,&*locVec.data().begin(),ADD_VALUES
  //     ); CHKERRQ(ierr);
  //     PetscFunctionReturn(0);
  //   }
  //
  // };


}

#endif //__POISSONOPERATORS_HPP__
