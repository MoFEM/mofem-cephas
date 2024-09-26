#ifndef __MIXED_NONLINEAR_POISSON_HPP__
#define __MIXED_NONLINEAR_POISSON_HPP__

#include <stdlib.h>
#include <MoFEM.hpp>
using namespace MoFEM;
constexpr int SPACE_DIM = EXECUTABLE_DIMENSION;

using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle =
    PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;
using EntData = EntitiesFieldData::EntData;

using OpHdivU = FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
    GAUSS>::OpMixDivTimesScalar<SPACE_DIM>;

using OpBaseTimesScalarRhs = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpBaseTimesScalar<1>;
using OpDomainSource = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;
using OpHDivTimesScalarRhs = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpMixDivTimesU<3, 1, SPACE_DIM>;
using OpBoundaryRhsSource = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpNormalMixVecTimesScalar<SPACE_DIM>;

using AssemblyDomainEleOp =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::OpBase;
using AssemblyBoundaryEleOp =
    FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase;

using EntData = EntitiesFieldData::EntData;

namespace MixedNonlinearPoissonOps {

FTensor::Index<'i', 2> i;
const double nonlinearConstant = 1;

typedef boost::function<double(const double, const double, const double)>
    ScalarFunc;

struct OpHdivUHdivRhs : public AssemblyDomainEleOp {
public:
  OpHdivUHdivRhs(std::string field_name,
                 boost::shared_ptr<VectorDouble> ufield_vec,
                 boost::shared_ptr<MatrixDouble> qfield_vec)
      : AssemblyDomainEleOp(field_name, field_name, AssemblyDomainEleOp::OPROW),
        ufieldVec(ufield_vec), qfieldVec(qfield_vec) {}

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;

    auto &nf = AssemblyDomainEleOp::locF;
    const size_t nb_base_functions = data.getN().size2() / 3;

    // get element area
    const double area = getMeasure();

    // get number of integration points
    const int nb_integration_points = getGaussPts().size2();
    // get integration weights
    auto t_w = getFTensor0IntegrationWeight();
    // get solution (ufield value) at integration point
    auto u_field = getFTensor0FromVec(*ufieldVec);
    // get solution (ufield value) at integration point
    auto q_field = getFTensor1FromMat<3>(*qfieldVec);

    // get base function
    auto q_base = data.getFTensor1N<3>();

    // START THE LOOP OVER INTEGRATION POINTS TO CALCULATE LOCAL VECTOR
    for (int gg = 0; gg != nb_integration_points; gg++) {
      const double a = t_w * area;
      int rr = 0;
      // calculate the local vector
      for (; rr != AssemblyDomainEleOp::nbRows; rr++) {
        // nf[rr] += q_base * (1 / (1 + u_field * u_field)) * q_field(i) * a;
        nf[rr] += q_base(i) *
                  (1 / (1 + nonlinearConstant * u_field * u_field)) *
                  q_field(i) * a;
        // move to the next base function
        ++q_base;
      }
      for (; rr < nb_base_functions; ++rr)
        ++q_base;
      // move to the weight of the next integration point
      ++t_w;
      // move to the solution (field value) at the next integration point
      ++u_field;
      // move to the gradient of field value at the next integration point
      ++q_field;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> ufieldVec;
  boost::shared_ptr<MatrixDouble> qfieldVec;
};

struct OpDomainQULhs : public AssemblyDomainEleOp {
public:
  OpDomainQULhs(std::string row_field_name, std::string col_field_name,
                boost::shared_ptr<VectorDouble> ufield_vec,
                boost::shared_ptr<MatrixDouble> qfield_vec)
      : AssemblyDomainEleOp(row_field_name, col_field_name,
                            DomainEleOp::OPROWCOL),
        ufieldVec(ufield_vec), qfieldVec(qfield_vec) {}

  MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
    MoFEMFunctionBegin;

    auto &locLhs = AssemblyDomainEleOp::locMat;

    const size_t nb_base_functions = row_data.getN().size2() / 3;

    const int nb_row_dofs = row_data.getIndices().size();
    const int nb_col_dofs = col_data.getIndices().size();
    // get element area
    const double area = getMeasure();

    // get number of integration points
    const int nb_integration_points = getGaussPts().size2();
    // get integration weights
    auto t_w = getFTensor0IntegrationWeight();
    // get solution (ufield value) at integration points
    auto u_field = getFTensor0FromVec(*ufieldVec);
    // get solution (ufield value) at integration points
    auto q_field = getFTensor1FromMat<3>(*qfieldVec);

    // get base functions on row
    auto q_row_base = row_data.getFTensor1N<3>();

    // START THE LOOP OVER INTEGRATION POINTS TO CALCULATE LOCAL MATRIX
    for (int gg = 0; gg != nb_integration_points; gg++) {
      const double a = t_w * area;

      int rr = 0;
      for (; rr != nb_row_dofs; ++rr) {
        // get base functions on column
        auto u_col_base = col_data.getFTensor0N(gg, 0);

        for (int cc = 0; cc != nb_col_dofs; cc++) {
          locLhs(rr, cc) += q_row_base(i) *
                            ((-2 * nonlinearConstant * u_field) /
                             ((1 + nonlinearConstant * u_field * u_field) *
                              (1 + nonlinearConstant * u_field * u_field))) *
                            u_col_base * q_field(i) * a;
          // move to the next base functions on column
          ++u_col_base;
        }

        // move to the next base function on row
        ++q_row_base;
      }
      for (; rr < nb_base_functions; ++rr)
        ++q_row_base;

      // move to the weight of the next integration point
      ++t_w;
      // move to the solution (field value) at the next integration point
      ++u_field;
      // move to the gradient of field value at the next integration point
      ++q_field;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> ufieldVec;
  boost::shared_ptr<MatrixDouble> qfieldVec;
};

struct OpDomainQQLhs : public AssemblyDomainEleOp {
public:
  OpDomainQQLhs(std::string row_field_name, std::string col_field_name,
                boost::shared_ptr<VectorDouble> ufield_vec)
      : AssemblyDomainEleOp(row_field_name, col_field_name,
                            DomainEleOp::OPROWCOL),
        ufieldVec(ufield_vec) {}

  MoFEMErrorCode iNtegrate(EntData &row_data, EntData &col_data) {
    MoFEMFunctionBegin;

    auto &locLhs = AssemblyDomainEleOp::locMat;
    const size_t nb_base_functions = row_data.getN().size2() / 3;

    const int nb_row_dofs = row_data.getIndices().size();
    const int nb_col_dofs = col_data.getIndices().size();
    // get element area
    const double area = getMeasure();

    // get number of integration points
    const int nb_integration_points = getGaussPts().size2();
    // get integration weights
    auto t_w = getFTensor0IntegrationWeight();
    // get solution (ufield value) at integration points
    auto u_field = getFTensor0FromVec(*ufieldVec);

    // get base functions on row
    auto q_row_base = row_data.getFTensor1N<3>();

    // START THE LOOP OVER INTEGRATION POINTS TO CALCULATE LOCAL MATRIX
    for (int gg = 0; gg != nb_integration_points; gg++) {
      const double a = t_w * area;
      int rr = 0;
      for (; rr != nb_row_dofs; ++rr) {
        // get base functions on column
        auto q_col_base = col_data.getFTensor1N<3>(gg, 0);

        for (int cc = 0; cc != nb_col_dofs; cc++) {
          locLhs(rr, cc) += q_row_base(i) *
                            (1 / (1 + nonlinearConstant * u_field * u_field)) *
                            q_col_base(i) * a;
          // move to the next base functions on column
          ++q_col_base;
        }

        // move to the next base function on row
        ++q_row_base;
      }
      for (; rr < nb_base_functions; ++rr)
        ++q_row_base;

      // move to the weight of the next integration point
      ++t_w;
      // move to the solution (field value) at the next integration point
      ++u_field;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> ufieldVec;
};

//! [OpFluxError]
struct OpErrorIndicator : public DomainEleOp {
public:
  OpErrorIndicator(boost::shared_ptr<VectorDouble> ufield_vec,
                   boost::shared_ptr<MatrixDouble> qfield_vec,
                   boost::shared_ptr<MatrixDouble> ufield_grad,
                   SmartPetscObj<Vec> data_vec, const int index)
      : DomainEleOp("U", OPROW), ufieldVec(ufield_vec), qfieldVec(qfield_vec),
        ufieldGradMat(ufield_grad), dataVec(data_vec), iNdex(index) {};
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
    // check if entity is in the range

    const int nb_integration_pts = getGaussPts().size2();
    const double area = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_val = getFTensor0FromVec(*ufieldVec);
    auto t_val_grad = getFTensor1FromMat<2>(*ufieldGradMat);
    auto t_flux = getFTensor1FromMat<3>(*qfieldVec);
    FTensor::Tensor1<double, 2> t_diff;
    FTensor::Index<'i', 2> i;

    double error_ind = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      const double alpha = t_w * area;

      t_diff(i) =
          t_val_grad(i) * (1 + nonlinearConstant * t_val * t_val) - t_flux(i);
      error_ind += alpha * t_diff(i) * t_diff(i);

      ++t_w;
      ++t_val;
      ++t_val_grad;
      ++t_flux;
    }
    error_ind *= area;

    CHKERR VecSetValue(dataVec, iNdex, error_ind, ADD_VALUES);
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> ufieldVec;
  boost::shared_ptr<MatrixDouble> qfieldVec;
  boost::shared_ptr<MatrixDouble> ufieldGradMat;
  SmartPetscObj<Vec> dataVec;
  const int iNdex;
};

//! [OpFluxError]

}; // namespace MixedNonlinearPoissonOps

#endif //__NONLINEARPOISSON2D_HPP__