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

/**
* \brief Get tensor rank 0 (scalar) form data vcetor
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<class T, class A>
FTensor::Tensor0<T*> getTensor0FormData(
  boost::shared_ptr<ublas::vector<T,A> > data_ptr
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor0<T*>();
}

template<>
FTensor::Tensor0<double*> getTensor0FormData<double,ublas::unbounded_array<double> >(
  boost::shared_ptr<ublas::vector<double,ublas::unbounded_array<double> > > data_ptr
);

/**
 * \brief Get tensor rank 1 (vector) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
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

/**
 * \brief Get tensor rank 1 (vector) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
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
);

template<>
FTensor::Tensor1<double*,2> getTensor1FormData<2,double,ublas::row_major,ublas::unbounded_array<double> >(
  boost::shared_ptr<MatrixDouble> data_ptr
);

/**
 * \brief Get tensor rank 2 (matrix) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
FTensor::Tensor2<T*,Tensor_Dim0,Tensor_Dim1> getTensor2FormData(
  boost::shared_ptr<ublas::matrix<T,L,A> > data_ptr
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  s << " and dim0 = " << Tensor_Dim0;
  s << " dim1 = " << Tensor_Dim1;
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor1<T*,Tensor_Dim>();
}

template<>
FTensor::Tensor2<double*,3,3> getTensor2FormData(
  boost::shared_ptr<MatrixDouble> data_ptr
);

/**
 * \brief Get tensor rank 2 (matrix) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getTensor1FormData(
  boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > data_ptr
) {
  return getTensor2FormData<
  Tensor_Dim0,Tensor_Dim1,double,ublas::row_major,ublas::unbounded_array<double>
  >(data_ptr);
}

// GET VALUES AT GAUSS PTS

// TENSOR0

/** \brief Calculate field values for tenor field rank 0, i.e. scalar field
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<class T, class A>
struct OpCalculateScalarFieldVaues_General: public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<ublas::vector<T,A> > dataPtr;
  EntityHandle zeroType;

  OpCalculateScalarFieldVaues_General(
    const std::string &field_name,
    boost::shared_ptr<ublas::vector<T,A> > data_ptr,
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

/**
* \brief Specialization of member function
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<class T, class A>
PetscErrorCode OpCalculateScalarFieldVaues_General<T,A>::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  SETERRQ1(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "Not implemented for T = %s",
    typeid(T).name()
  );
  PetscFunctionReturn(0);
}


/**
 * \brief Get value for scalar field
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
struct OpCalculateScalarFieldVaues:
public OpCalculateScalarFieldVaues_General<double,ublas::unbounded_array<double> > {

  OpCalculateScalarFieldVaues(
    const std::string &field_name,
    boost::shared_ptr<VectorDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  );

};

// TENSOR1

/** \brief Calculate field values for tenor field rank 1, i.e. vector field
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim, class T, class L, class A>
struct OpCalculateVectorFieldValues_General: public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T,L,A> > dataPtr;
  EntityHandle zeroType;

  OpCalculateVectorFieldValues_General(
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

};

template<int Tensor_Dim, class T, class L, class A>
PetscErrorCode OpCalculateVectorFieldValues_General<Tensor_Dim,T,L,A>::doWork(
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

/** \brief Calculate field values (template specialization) for tensor field rank 1, i.e. vector field
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim>
struct OpCalculateVectorFieldValues_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >:
  public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  EntityHandle zeroType;

  OpCalculateVectorFieldValues_General(
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

/**
 * \brief Member function specialization calculating values for tenor field rank 1
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim>
PetscErrorCode OpCalculateVectorFieldValues_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >::doWork(
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
  MatrixDouble &mat = *dataPtr;
  if(type == zeroType) {
    mat.resize(Tensor_Dim,nb_gauss_pts,false);
    mat.clear();
  }
  FTensor::Tensor0<double*> base_function = data.getFTensor0N();
  FTensor::Tensor1<double*,Tensor_Dim> values_at_gauss_pts = getTensor1FormData<Tensor_Dim>(dataPtr);
  FTensor::Index<'I',Tensor_Dim> I;
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    FTensor::Tensor1<double*,Tensor_Dim> field_data = data.getFTensor1FieldData<Tensor_Dim>();
    for(int bb = 0;bb!=nb_base_functions;bb++) {
      if(bb*Tensor_Dim < nb_dofs) { // Number of dofs can be smaller than number of 3 x base functions
        values_at_gauss_pts(I) = field_data(I)*base_function;
        ++field_data;
      }
      ++base_function;
    }
    ++values_at_gauss_pts;
  }
  PetscFunctionReturn(0);
}

/** \brief Get values at integration pts for tensor filed rank 1, i.e. vector field
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim>
struct OpCalculateVectorFieldValues:
public OpCalculateVectorFieldValues_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> > {

  OpCalculateVectorFieldValues(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  OpCalculateVectorFieldValues_General<
  Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double>
  >(field_name,data_ptr,zero_type) {
  }

};

/** \brief Calculate field values for tenor field rank 2.
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim0,int Tensor_Dim1, class T, class L, class A>
struct OpCalculateTensor2FieldValues_General: public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T,L,A> > dataPtr;
  EntityHandle zeroType;

  OpCalculateTensor2FieldValues_General(
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

};

template<int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
PetscErrorCode OpCalculateTensor2FieldValues_General<Tensor_Dim0,Tensor_Dim1,T,L,A>::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  SETERRQ3(
    PETSC_COMM_SELF,
    MOFEM_NOT_IMPLEMENTED,
    "Not implemented for T = %s, dim0 = %d and dim1 = %d",
    typeid(T).name(),
    Tensor_Dim0,
    Tensor_Dim1
  );
  PetscFunctionReturn(0);
}

template<int Tensor_Dim0,int Tensor_Dim1>
struct OpCalculateTensor2FieldValues_General<
Tensor_Dim0,Tensor_Dim1,double,ublas::row_major,ublas::unbounded_array<double>
>: public ForcesAndSurcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > dataPtr;
  EntityHandle zeroType;

  OpCalculateTensor2FieldValues_General(
    const std::string &field_name,
    boost::shared_ptr<ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > > data_ptr,
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

template<int Tensor_Dim0,int Tensor_Dim1>
PetscErrorCode OpCalculateTensor2FieldValues_General<
Tensor_Dim0,Tensor_Dim1, double, ublas::row_major, ublas::unbounded_array<double>
>::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"Not implemented yet");
  PetscFunctionReturn(0);
}

// GET GRADIENTS AT GAUSS POINTS

/**
 * \brief Evaluate field gradient values for scalar field, i.e. gradient is tensor rank 1 (vector)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim, class T, class L, class A>
struct OpCalculateScalarFieldGradient_General:
public OpCalculateVectorFieldValues_General<Tensor_Dim,T,L,A> {

  OpCalculateScalarFieldGradient_General(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  OpCalculateVectorFieldValues_General<Tensor_Dim,T,L,A>(
    field_name,data_ptr,zero_type
  ) {
  }

};

/** \brief Evaluate field gradient values for scalar field, i.e. gradient is tensor rank 1 (vector), specialization
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim>
struct OpCalculateScalarFieldGradient_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >:
public OpCalculateVectorFieldValues_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> > {

  OpCalculateScalarFieldGradient_General(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  OpCalculateVectorFieldValues_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >(
    field_name,data_ptr,zero_type
  ) {
  }

  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );

};

/**
 * \brief Member function specialization calculating scalar field gradients for tenor field rank 1
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim>
PetscErrorCode OpCalculateScalarFieldGradient_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> >::doWork(
  int side,EntityType type,DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  if(!nb_dofs) {
    this->dataPtr->resize(Tensor_Dim,0,false);
    PetscFunctionReturn(0);
  }
  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();
  ublas::matrix<double,ublas::row_major,ublas::unbounded_array<double> > &mat = *this->dataPtr;
  if(type == this->zeroType) {
    mat.resize(Tensor_Dim,nb_gauss_pts,false);
    mat.clear();
  }
  FTensor::Tensor1<double*,Tensor_Dim> diff_base_function = data.getFTensor1DiffN<Tensor_Dim>();
  FTensor::Tensor1<double*,Tensor_Dim> gradients_at_gauss_pts = getTensor1FormData<Tensor_Dim>(this->dataPtr);
  FTensor::Index<'I',Tensor_Dim> I;
  for(int gg = 0;gg<nb_gauss_pts;gg++) {
    FTensor::Tensor0<double*> field_data = data.getFTensor0FieldData();
    for(int bb = 0;bb<nb_base_functions;bb++) {
      if(bb*Tensor_Dim < nb_dofs) { // Number of dofs can be smaller than number of 3 x base functions
        gradients_at_gauss_pts(I) += field_data*diff_base_function(I);
        ++field_data;
      }
      ++diff_base_function;
    }
    ++gradients_at_gauss_pts;
  }
  PetscFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for scalar filed rank 0, i.e. vector field
* \ingroup mofem_forces_and_sources_user_data_operators
*/
template<int Tensor_Dim>
struct OpCalculateScalarFieldGradient:
public OpCalculateScalarFieldGradient_General<Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double> > {

  OpCalculateScalarFieldGradient(
    const std::string &field_name,
    boost::shared_ptr<MatrixDouble> data_ptr,
    EntityType zero_type = MBVERTEX
  ):
  OpCalculateScalarFieldGradient_General<
  Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double>
  >(field_name,data_ptr,zero_type) {
  }

};

}

#endif // __USER_DATA_OPERATORS_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_user_data_operators Users Operators
 *
 * \brief Classes and functions used to evaluate fields at integration pts, jacobians, etc..
 *
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
