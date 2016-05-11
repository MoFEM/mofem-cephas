/** \file UsersDataOperators.hpp

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

#ifndef __USER_DATA_OPERATORS_HPP__
#define __USER_DATA_OPERATORS_HPP__

namespace MoFEM {

template<int Tensor_Dim, class T, class L, class A>
FTensor::Tensor1<T*,Tensor_Dim> getTensor1FormData(
  boost::shared_ptr<ublas::matrix<T,L,A> > data_ptr
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  s << " and dim = " << Tensor_Dim;
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor1<T*,Tensor_Dim>();
}

template<int Tensor_Dim>
FTensor::Tensor1<double*,Tensor_Dim> getTensor1FormData(
  boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > data_ptr
) {
  return getTensor1FormData<
  Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double>
  >(data_ptr);
}

template<>
FTensor::Tensor1<double*,3> getTensor1FormData<3,double,ublas::row_major,ublas::unbounded_array<double> >(
  boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > data_ptr
) {
  PetscFunctionBegin;
  if(data_ptr->size1()!=3) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<double*,3>(
    &(*data_ptr)(0,0),&(*data_ptr)(1,0),&(*data_ptr)(2,0)
  );
}

/** \brief Calculate field values
*/
template<int Tensor_Dim, class T, class L, class A>
struct OpCalculateFieldValues_Tensor1_General: public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T,L,A> > dataPtr;
  EntityHandle zeroType;

  OpCalculateFieldValues_Tensor1_General(
    const std::string &field_name,
    boost::shared_ptr<ublas::matrix<T,L,A> > data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  ForcesAndSurcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
  dataPtr(data_ptr),
  zeroType(zero_type) {
  }

  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );

  // {
  //   PetscFunctionBegin;
  //   SETERRQ2(
  //     PETSC_COMM_SELF,
  //     MOFEM_NOT_IMPLEMENTED,
  //     "Not implemented for T = %s and dim = %d",
  //     typeid(T).name(),
  //     Tensor_Dim
  //   );
  //   PetscFunctionReturn(0);
  // }

};

template<int Tensor_Dim, class T, class L, class A>
PetscErrorCode OpCalculateFieldValues_Tensor1_General<Tensor_Dim,T,L,A>::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  SETERRQ2(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "Not implemented for T = %s and dim = %d",
    typeid(T).name(),
    Tensor_Dim
  );
  PetscFunctionReturn(0);
}


/** \brief Calculate field values
*/
template<int Tensor_Dim>
struct OpCalculateFieldValues_Tensor1_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >:
  public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  EntityHandle zeroType;

  OpCalculateFieldValues_Tensor1_General(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  ForcesAndSurcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
  dataPtr(data_ptr),
  zeroType(zero_type) {
  }

  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );

};

template<int Tensor_Dim>
PetscErrorCode OpCalculateFieldValues_Tensor1_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  const int nb_dofs = data.getFieldData().size();
  if(!nb_dofs) {
    dataPtr->resize(Tensor_Dim,0,false);
    PetscFunctionReturn(0);
  }

  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();

  ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > &mat = *dataPtr;
  if(type == zeroType) {
    mat.resize(Tensor_Dim,nb_gauss_pts,false);
    mat.clear();
  }

  FTensor::Tensor0<double*> base_function = data.getFTensor0N();
  FTensor::Tensor1<double*,Tensor_Dim> values_at_gauss_pts = getTensor1FormData<Tensor_Dim>(dataPtr);

  FTensor::Index<'I',Tensor_Dim> I;

  for(int gg = 0;gg<nb_gauss_pts;gg++) {

    Tensor1<double*,Tensor_Dim> field_data = data.getFTensor1FieldData<3>();
    for(int bb = 0;bb<nb_base_functions;bb++) {
      if(bb*Tensor_Dim < nb_dofs) { // Number of dofs can be smaller than number of 3 x base functions
        values_at_gauss_pts(I) += field_data(I)*base_function;
        ++field_data;
      }
      ++base_function;
    }
    ++values_at_gauss_pts;

  }

  PetscFunctionReturn(0);
}

template<int Tensor_Dim>
struct OpCalculateFieldValues_Tensor1:
public OpCalculateFieldValues_Tensor1_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> > {

  OpCalculateFieldValues_Tensor1(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  OpCalculateFieldValues_Tensor1_General<
  Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double>
  >(field_name,data_ptr,zero_type) {
  }

};


}

#endif // __USER_DATA_OPERATORS_HPP__
