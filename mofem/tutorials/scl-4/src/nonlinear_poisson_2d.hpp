#ifndef __NONLINEARPOISSON2D_HPP__
#define __NONLINEARPOISSON2D_HPP__

#include <stdlib.h>
#include <MoFEM.hpp>
using namespace MoFEM;

using DomainEle = PipelineManager::FaceEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle = PipelineManager::EdgeEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;


using OpBoundaryLhs = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::BiLinearForm<GAUSS>::OpMass<1, 1>;
using OpBoundaryRhs = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpBaseTimesScalar<1>;
using OpBoundaryRhsSource = FormsIntegrators<BoundaryEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<1, 1>;


using AssemblyDomainEleOp =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::OpBase;
using AssemblyBoundaryEleOp =
    FormsIntegrators<BoundaryEleOp>::Assembly<PETSC>::OpBase;


using EntData = EntitiesFieldData::EntData;

namespace NonlinearPoissonOps {

FTensor::Index<'i', 2> i;

typedef boost::function<double(const double, const double, const double)>
    ScalarFunc;


struct OpDomainLhs : public AssemblyDomainEleOp {
public:
  OpDomainLhs(
      std::string row_field_name, std::string col_field_name,
      boost::shared_ptr<VectorDouble> field_vec,
      boost::shared_ptr<MatrixDouble> field_grad_mat)
      : AssemblyDomainEleOp(row_field_name, col_field_name, 
      DomainEleOp::OPROWCOL),
      fieldVec(field_vec), fieldGradMat(field_grad_mat) {}

  MoFEMErrorCode iNtegrate(EntData &row_data,
                        EntData &col_data) {
    MoFEMFunctionBegin;

    auto &locLhs = AssemblyDomainEleOp::locMat;

    const int nb_row_dofs = row_data.getIndices().size();
    const int nb_col_dofs = col_data.getIndices().size();
    // get element area
    const double area = getMeasure();

    // get number of integration points
    const int nb_integration_points = getGaussPts().size2();
    // get integration weights
    auto t_w = getFTensor0IntegrationWeight();
    // get solution (field value) at integration points
    auto t_field = getFTensor0FromVec(*fieldVec);
    // get gradient of field at integration points
    auto t_field_grad = getFTensor1FromMat<2>(*fieldGradMat);

    // get base functions on row
    auto t_row_base = row_data.getFTensor0N();
    // get derivatives of base functions on row
    auto t_row_diff_base = row_data.getFTensor1DiffN<2>();

    // START THE LOOP OVER INTEGRATION POINTS TO CALCULATE LOCAL MATRIX
    for (int gg = 0; gg != nb_integration_points; gg++) {
      const double a = t_w * area;

      for (int rr = 0; rr != nb_row_dofs; ++rr) {
        // get base functions on column
        auto t_col_base = col_data.getFTensor0N(gg, 0);
        // get derivatives of base functions on column
        auto t_col_diff_base = col_data.getFTensor1DiffN<2>(gg, 0);

        for (int cc = 0; cc != nb_col_dofs; cc++) {
          locLhs(rr, cc) += (((1 + t_field * t_field) * t_row_diff_base(i) *
                              t_col_diff_base(i)) +
                              (2.0 * t_field * t_field_grad(i) *
                              t_row_diff_base(i) * t_col_base)) *
                            a;

          // move to the next base functions on column
          ++t_col_base;
          // move to the derivatives of the next base function on column
          ++t_col_diff_base;
        }

        // move to the next base function on row
        ++t_row_base;
        // move to the derivatives of the next base function on row
        ++t_row_diff_base;
      }

      // move to the weight of the next integration point
      ++t_w;
      // move to the solution (field value) at the next integration point
      ++t_field;
      // move to the gradient of field value at the next integration point
      ++t_field_grad;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> fieldVec;
  boost::shared_ptr<MatrixDouble> fieldGradMat;
};

struct OpDomainRhs : public AssemblyDomainEleOp {
public:
  OpDomainRhs(
      std::string field_name, ScalarFunc source_term_function,
      boost::shared_ptr<VectorDouble> field_vec, 
      boost::shared_ptr<MatrixDouble> field_grad_mat)
      : AssemblyDomainEleOp(field_name, field_name, AssemblyDomainEleOp::OPROW),
        sourceTermFunc(source_term_function), fieldVec(field_vec),
        fieldGradMat(field_grad_mat) {}

  MoFEMErrorCode iNtegrate(EntData &data) {
    MoFEMFunctionBegin;

    auto &nf = AssemblyDomainEleOp::locF;


    // get element area
    const double area = getMeasure();

    // get number of integration points
    const int nb_integration_points = getGaussPts().size2();
    // get integration weights
    auto t_w = getFTensor0IntegrationWeight();
    // get coordinates of the integration point
    auto t_coords = getFTensor1CoordsAtGaussPts();
    // get solution (field value) at integration point
    auto t_field = getFTensor0FromVec(*fieldVec);
    // get gradient of field value of integration point
    auto t_field_grad = getFTensor1FromMat<2>(*fieldGradMat);

    // get base function
    auto t_base = data.getFTensor0N();
    // get derivatives of base function
    auto t_grad_base = data.getFTensor1DiffN<2>();

    // START THE LOOP OVER INTEGRATION POINTS TO CALCULATE LOCAL VECTOR
    for (int gg = 0; gg != nb_integration_points; gg++) {
      const double a = t_w * area;
      double body_source =
          sourceTermFunc(t_coords(0), t_coords(1), t_coords(2));

      // calculate the local vector
      for (int rr = 0; rr != AssemblyDomainEleOp::nbRows; rr++) {
        nf[rr] +=
            (-t_base * body_source +
              t_grad_base(i) * t_field_grad(i) * (1 + t_field * t_field)) *
            a;

        // move to the next base function
        ++t_base;
        // move to the derivatives of the next base function
        ++t_grad_base;
      }

      // move to the weight of the next integration point
      ++t_w;
      // move to the coordinates of the next integration point
      ++t_coords;
      // move to the solution (field value) at the next integration point
      ++t_field;
      // move to the gradient of field value at the next integration point
      ++t_field_grad;

    }

    MoFEMFunctionReturn(0);
  }

private:
  ScalarFunc sourceTermFunc;
  boost::shared_ptr<VectorDouble> fieldVec;
  boost::shared_ptr<MatrixDouble> fieldGradMat;

};

}; // namespace NonlinearPoissonOps

#endif //__NONLINEARPOISSON2D_HPP__