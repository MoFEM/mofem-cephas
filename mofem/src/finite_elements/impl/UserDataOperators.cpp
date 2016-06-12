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
PetscErrorCode OpCalculateScalarFieldValues_General<double,ublas::unbounded_array<double> >::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  // cerr <<  data.getFieldData() << endl;

  if(!dataPtr) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Data pointer not allocated");
  }

  if(!nb_dofs && type == this->zeroType) {
    dataPtr->resize(0,false);
  }
  if(!nb_dofs) {
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
  FTensor::Tensor0<double*> values_at_gauss_pts = getTensor0FormData(vec);
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    FTensor::Tensor0<double*> field_data = data.getFTensor0FieldData();
    int bb = 0;
    for(;bb<nb_dofs;bb++) {
      values_at_gauss_pts += field_data*base_function;
      ++field_data;
      ++base_function;
    }
    // It is possible to have more base functions than dofs
    for(;bb!=nb_base_functions;bb++) ++base_function;
    ++values_at_gauss_pts;
  }
  PetscFunctionReturn(0);
}

OpCalculateScalarFieldValues::OpCalculateScalarFieldValues(
  const std::string &field_name,
  boost::shared_ptr<VectorDouble> data_ptr,
  EntityType zero_type
):
OpCalculateScalarFieldValues_General<double,ublas::unbounded_array<double> >(
  field_name,data_ptr,zero_type
) {
}

}
