/** \file ElementsOnEntities.cpp

\brief Implementation of Elements on Entities for Forces and Sources

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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>
#include <DataOperators.hpp>
#include <ElementsOnEntities.hpp>
#include <UserDataOperators.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

template<>
FTensor::Tensor0<double*> getTensor0FormData<double,ublas::unbounded_array<double> >(
  boost::shared_ptr<ublas::vector<double,ublas::unbounded_array<double> > > data_ptr
) {
  return FTensor::Tensor0<double*>(&*data_ptr->data().begin());
}

template<>
FTensor::Tensor1<double*,3> getTensor1FormData<3,double,ublas::row_major,ublas::unbounded_array<double> >(
  boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > data_ptr
) {
  if(data_ptr->size1()!=3) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<double*,3>(
    &(*data_ptr)(0,0),&(*data_ptr)(1,0),&(*data_ptr)(2,0)
  );
}

template<>
FTensor::Tensor1<double*,2> getTensor1FormData<2,double,ublas::row_major,ublas::unbounded_array<double> >(
  boost::shared_ptr<MatrixDouble> data_ptr
) {
  if(data_ptr->size1()!=2) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<double*,2>(
    &(*data_ptr)(0,0),&(*data_ptr)(1,0)
  );
}

template<>
FTensor::Tensor2<double*,3,3> getTensor2FormData(
  boost::shared_ptr<MatrixDouble> data_ptr
) {
  if(data_ptr->size1()!=9) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  MatrixDouble &mat = *data_ptr;
  return FTensor::Tensor2<double*,3,3>(
    &mat(0,0),&mat(1,0),&mat(2,0),&mat(3,0),&mat(4,0),&mat(5,0),&mat(6,0),&mat(7,0),&mat(8,0)
  );
}

template<>
PetscErrorCode OpCalculateScalarFieldVaues_General<double,ublas::unbounded_array<double> >::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  if(!nb_dofs) {
    dataPtr->resize(0,false);
    PetscFunctionReturn(0);
  }
  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();
  VectorDouble &vec = *dataPtr;
  if(type == zeroType) {
    vec.resize(nb_gauss_pts,false);
    vec.clear();
  }
  FTensor::Tensor0<double*> base_function = data.getFTensor0N();
  FTensor::Tensor0<double*> values_at_gauss_pts = getTensor0FormData(dataPtr);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    Tensor0<double*> field_data = data.getFTensor0FieldData();
    for(int bb = 0;bb<nb_base_functions;bb++) {
       if(bb < nb_dofs) { // Number of dofs can be smaller than number of base functions
         values_at_gauss_pts += field_data*base_function;
         ++field_data;
       }
       ++base_function;
    }
    ++values_at_gauss_pts;
  }
  PetscFunctionReturn(0);
}

OpCalculateScalarFieldVaues::OpCalculateScalarFieldVaues(
  const std::string &field_name,
  boost::shared_ptr<VectorDouble> data_ptr,
  EntityType zero_type
):
OpCalculateScalarFieldVaues_General<double,ublas::unbounded_array<double> >(
  field_name,data_ptr,zero_type
) {
}


}
