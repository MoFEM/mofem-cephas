/** \file UserDataOperators.hpp
  * \brief User data Operators

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

// GET VALUES AT GAUSS PTS

// TENSOR0

/** \brief Calculate field values for tenor field rank 0, i.e. scalar field
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <class T, class A>
struct OpCalculateScalarFieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<ublas::vector<T, A>> dataPtr;
  const EntityHandle zeroType;

  OpCalculateScalarFieldValues_General(
      const std::string &field_name,
      boost::shared_ptr<ublas::vector<T, A>> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Specialization of member function
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <class T, class A>
MoFEMErrorCode OpCalculateScalarFieldValues_General<T, A>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented for T = %s",
           typeid(T).name() // boost::core::demangle(typeid(T).name()).c_str()
  );
  MoFEMFunctionReturn(0);
}

/**
 * \brief Get value at integration points for scalar field
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
struct OpCalculateScalarFieldValues
    : public OpCalculateScalarFieldValues_General<double, DoubleAllocator> {

  OpCalculateScalarFieldValues(const std::string &field_name,
                               boost::shared_ptr<VectorDouble> &data_ptr,
                               const EntityType zero_type = MBVERTEX)
      : OpCalculateScalarFieldValues_General<double, DoubleAllocator>(
            field_name, data_ptr, zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief calculate values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    VectorDouble &vec = *dataPtr;
    const int nb_gauss_pts = getGaussPts().size2();
    if (type == zeroType) {
      vec.resize(nb_gauss_pts, false);
      vec.clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (nb_dofs) {
      const int nb_gauss_pts = getGaussPts().size2();
      const int nb_base_functions = data.getN().size2();
      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts = getFTensor0FromVec(vec);
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor0FieldData();
        int bb = 0;
        for (; bb != nb_dofs; ++bb) {
          values_at_gauss_pts += field_data * base_function;
          ++field_data;
          ++base_function;
        }
        // It is possible to have more base functions than dofs
        for (; bb != nb_base_functions; ++bb)
          ++base_function;
        ++values_at_gauss_pts;
      }
    }
    MoFEMFunctionReturn(0);
  }
};

struct OpCalculateScalarValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroAtType;
  VectorDouble dotVector;

  OpCalculateScalarValuesDot(const std::string field_name,
                                  boost::shared_ptr<VectorDouble> &data_ptr,
                                  const EntityType zero_at_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroAtType(zero_at_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;

    const int nb_gauss_pts = getGaussPts().size2();
    VectorDouble &vec = *dataPtr;
    if (type == zeroAtType) {
      vec.resize(nb_gauss_pts, false);
      vec.clear();
    }

    auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    if (nb_dofs) {

      dotVector.resize(nb_dofs, false);
      const double *array;
      CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
      for (int i = 0; i != local_indices.size(); ++i)
        if (local_indices[i] != -1)
          dotVector[i] = array[local_indices[i]];
        else
          dotVector[i] = 0;
      CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

      const int nb_base_functions = data.getN().size2();
      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts = getFTensor0FromVec(vec);

      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = getFTensor0FromVec(dotVector);
        int bb = 0;
        for (; bb != nb_dofs; ++bb) {
          values_at_gauss_pts += field_data * base_function;
          ++field_data;
          ++base_function;
        }
        // Number of dofs can be smaller than number of Tensor_Dim x base
        // functions
        for (; bb != nb_base_functions; ++bb)
          ++base_function;
        ++values_at_gauss_pts;
      }
    }
    MoFEMFunctionReturn(0);
  }
};

// TENSOR1

/** \brief Calculate field values for tenor field rank 1, i.e. vector field
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim, class T, class L, class A>
struct OpCalculateVectorFieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;

  OpCalculateVectorFieldValues_General(
      const std::string &field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief calculate values of vector field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

template <int Tensor_Dim, class T, class L, class A>
MoFEMErrorCode
OpCalculateVectorFieldValues_General<Tensor_Dim, T, L, A>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  SETERRQ2(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
           "Not implemented for T = %s and dim = %d",
           typeid(T).name(), // boost::core::demangle(typeid(T).name()),
           Tensor_Dim);
  MoFEMFunctionReturnHot(0);
}

/** \brief Calculate field values (template specialization) for tensor field
 * rank 1, i.e. vector field \ingroup
 * mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateVectorFieldValues_General<Tensor_Dim, double,
                                            ublas::row_major, DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;

  OpCalculateVectorFieldValues_General(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Member function specialization calculating values for tenor field rank
 * 1 \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
MoFEMErrorCode OpCalculateVectorFieldValues_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  if (!nb_dofs && type == this->zeroType) {
    dataPtr->resize(Tensor_Dim, 0, false);
    MoFEMFunctionReturnHot(0);
  }
  if (!nb_dofs) {
    MoFEMFunctionReturnHot(0);
  }
  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();
  MatrixDouble &mat = *dataPtr;
  if (type == zeroType) {
    mat.resize(Tensor_Dim, nb_gauss_pts, false);
    mat.clear();
  }
  auto base_function = data.getFTensor0N();
  auto values_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);
  FTensor::Index<'I', Tensor_Dim> I;
  const int size = nb_dofs / Tensor_Dim;
  if (nb_dofs % Tensor_Dim) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    auto field_data = data.getFTensor1FieldData<Tensor_Dim>();
    int bb = 0;
    for (; bb != size; ++bb) {
      values_at_gauss_pts(I) += field_data(I) * base_function;
      ++field_data;
      ++base_function;
    }
    // Number of dofs can be smaller than number of Tensor_Dim x base functions
    for (; bb != nb_base_functions; ++bb)
      ++base_function;
    ++values_at_gauss_pts;
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get values at integration pts for tensor filed rank 1, i.e. vector
 * field \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateVectorFieldValues
    : public OpCalculateVectorFieldValues_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  OpCalculateVectorFieldValues(const std::string &field_name,
                               boost::shared_ptr<MatrixDouble> &data_ptr,
                               const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldValues_General<Tensor_Dim, double,
                                             ublas::row_major, DoubleAllocator>(
            field_name, data_ptr, zero_type) {}
};

/** \brief Get time direvatives of values at integration pts for tensor filed
 * rank 1, i.e. vector field \ingroup
 * mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateVectorFieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroAtType;
  VectorDouble dotVector;

  template <int Dim> inline auto getFTensorDotData() {
    static_assert(Dim || !Dim, "not implemented");
  }

  OpCalculateVectorFieldValuesDot(const std::string field_name,
                                  boost::shared_ptr<MatrixDouble> &data_ptr,
                                  const EntityType zero_at_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroAtType(zero_at_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    if (!nb_dofs && type == zeroAtType) {
      dataPtr->resize(Tensor_Dim, 0, false);
      MoFEMFunctionReturnHot(0);
    }
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    dotVector.resize(nb_dofs, false);
    const double *array;
    CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;
    CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

    const int nb_gauss_pts = data.getN().size1();
    const int nb_base_functions = data.getN().size2();
    MatrixDouble &mat = *dataPtr;
    if (type == zeroAtType) {
      mat.resize(Tensor_Dim, nb_gauss_pts, false);
      mat.clear();
    }
    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);
    FTensor::Index<'I', Tensor_Dim> I;
    const int size = nb_dofs / Tensor_Dim;
    if (nb_dofs % Tensor_Dim) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
    }
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = getFTensorDotData<3>();
      int bb = 0;
      for (; bb != size; ++bb) {
        values_at_gauss_pts(I) += field_data(I) * base_function;
        ++field_data;
        ++base_function;
      }
      // Number of dofs can be smaller than number of Tensor_Dim x base
      // functions
      for (; bb != nb_base_functions; ++bb)
        ++base_function;
      ++values_at_gauss_pts;
    }
    MoFEMFunctionReturn(0);
  }
};

template <>
template <>
inline auto OpCalculateVectorFieldValuesDot<3>::getFTensorDotData<3>() {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      &dotVector[0], &dotVector[1], &dotVector[2]);
}

/** \brief Calculate field values for tenor field rank 2.
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
struct OpCalculateTensor2FieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;

  OpCalculateTensor2FieldValues_General(
      const std::string &field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
MoFEMErrorCode
OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T, L, A>::
    doWork(int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  SETERRQ3(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
           "Not implemented for T = %s, dim0 = %d and dim1 = %d",
           typeid(T).name(), // boost::core::demangle(typeid(T).name()),
           Tensor_Dim0, Tensor_Dim1);
  MoFEMFunctionReturnHot(0);
}

template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, double,
                                             ublas::row_major, DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<double, ublas::row_major, DoubleAllocator>>
      dataPtr;
  const EntityHandle zeroType;

  OpCalculateTensor2FieldValues_General(
      const std::string &field_name,
      boost::shared_ptr<
          ublas::matrix<double, ublas::row_major, DoubleAllocator>> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

template <int Tensor_Dim0, int Tensor_Dim1>
MoFEMErrorCode OpCalculateTensor2FieldValues_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  const int nb_dofs = data.getFieldData().size();
  const int nb_gauss_pts = data.getN().size1();
  if (!nb_dofs && type == this->zeroType) {
    dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, 0, false);
    MoFEMFunctionReturnHot(0);
  }
  if (!nb_dofs) {
    MoFEMFunctionReturnHot(0);
  }
  const int nb_base_functions = data.getN().size2();
  MatrixDouble &mat = *dataPtr;
  if (type == zeroType) {
    mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
    mat.clear();
  }
  auto base_function = data.getFTensor0N();
  auto values_at_gauss_pts = getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
  FTensor::Index<'i', Tensor_Dim0> i;
  FTensor::Index<'j', Tensor_Dim1> j;
  const int size = nb_dofs / (Tensor_Dim0 * Tensor_Dim1);
  for (int gg = 0; gg != nb_gauss_pts; ++gg) {
    auto field_data = data.getFTensor2FieldData<Tensor_Dim0, Tensor_Dim1>();
    int bb = 0;
    for (; bb != size; ++bb) {
      values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
      ++field_data;
      ++base_function;
    }
    for (; bb != nb_base_functions; ++bb)
      ++base_function;
    ++values_at_gauss_pts;
  }
  MoFEMFunctionReturnHot(0);
}

/** \brief Get values at integration pts for tensor filed rank 2, i.e. matrix
 * field \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2FieldValues
    : public OpCalculateTensor2FieldValues_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  OpCalculateTensor2FieldValues(const std::string &field_name,
                                boost::shared_ptr<MatrixDouble> &data_ptr,
                                const EntityType zero_type = MBVERTEX)
      : OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, double,
                                              ublas::row_major,
                                              DoubleAllocator>(
            field_name, data_ptr, zero_type) {}
};

/** \brief Get time direvarive values at integration pts for tensor filed rank
 * 2, i.e. matrix field
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2FieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> &dataPtr; ///< Data computed into this matrix
  EntityType zeroAtType;  ///< Zero values at Gauss point at this type
  VectorDouble dotVector; ///< Keeps temoorary values of time directives

  template <int Dim0, int Dim1> auto getFTensorDotData() {
    static_assert(Dim0 || !Dim0 || Dim1 || !Dim1, "not implemented");
  }

  OpCalculateTensor2FieldValuesDot(const std::string field_name,
                                   boost::shared_ptr<MatrixDouble> &data_ptr,
                                   const EntityType zero_at_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroAtType(zero_at_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBeginHot;
    const auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    const int nb_gauss_pts = data.getN().size1();
    if (!nb_dofs && type == zeroAtType) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, 0, false);
      MoFEMFunctionReturnHot(0);
    }
    if (!nb_dofs) {
      MoFEMFunctionReturnHot(0);
    }

    dotVector.resize(nb_dofs, false);
    const double *array;
    CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;
    CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

    const int nb_base_functions = data.getN().size2();
    MatrixDouble &mat = *dataPtr;
    if (type == zeroAtType) {
      mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
      mat.clear();
    }
    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts =
        getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    const int size = nb_dofs / (Tensor_Dim0 * Tensor_Dim1);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = getFTensorDotData<Tensor_Dim0, Tensor_Dim1>();
      int bb = 0;
      for (; bb != size; ++bb) {
        values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
        ++field_data;
        ++base_function;
      }
      for (; bb != nb_base_functions; ++bb)
        ++base_function;
      ++values_at_gauss_pts;
    }
    MoFEMFunctionReturnHot(0);
  }
};

template <>
template <>
inline auto OpCalculateTensor2FieldValuesDot<3, 3>::getFTensorDotData<3, 3>() {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      &dotVector[0], &dotVector[1], &dotVector[2],

      &dotVector[3], &dotVector[4], &dotVector[5],

      &dotVector[6], &dotVector[7], &dotVector[8]);
}

template <int Tensor_Dim>
struct OpCalculateTensor2SymmetricFieldValues
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateTensor2SymmetricFieldValues(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBTRI, const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    MatrixDouble &mat = *dataPtr;
    const int nb_gauss_pts = getGaussPts().size2();
    if (type == this->zeroType && side == zeroSide) {
      mat.resize((Tensor_Dim * (Tensor_Dim + 1)) / 2, nb_gauss_pts, false);
      mat.clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs) {
      MoFEMFunctionReturnHot(0);
    }
    const int nb_base_functions = data.getN().size2();
    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts = getFTensor2SymmetricFromMat<Tensor_Dim>(mat);
    FTensor::Index<'i', Tensor_Dim> i;
    FTensor::Index<'j', Tensor_Dim> j;
    const int size = nb_dofs / ((Tensor_Dim * (Tensor_Dim + 1)) / 2);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = data.getFTensor2SymmetricFieldData<Tensor_Dim>();
      int bb = 0;
      for (; bb != size; ++bb) {
        values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
        ++field_data;
        ++base_function;
      }
      for (; bb != nb_base_functions; ++bb)
        ++base_function;
      ++values_at_gauss_pts;
    }

    MoFEMFunctionReturn(0);
  }
};

template <int Tensor_Dim>
struct OpCalculateTensor2SymmetricFieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
  VectorDouble dotVector;

  template <int Dim> inline auto getFTensorDotData() {
    static_assert(Dim || !Dim, "not implemented");
  }

  OpCalculateTensor2SymmetricFieldValuesDot(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBTRI, const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_gauss_pts = getGaussPts().size2();
    MatrixDouble &mat = *dataPtr;
    if (type == zeroType && side == zeroSide) {
      mat.resize((Tensor_Dim * (Tensor_Dim + 1)) / 2, nb_gauss_pts, false);
      mat.clear();
    }
    auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    dotVector.resize(nb_dofs, false);
    const double *array;
    CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;
    CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

    const int nb_base_functions = data.getN().size2();

    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts = getFTensor2SymmetricFromMat<Tensor_Dim>(mat);
    FTensor::Index<'i', Tensor_Dim> i;
    FTensor::Index<'j', Tensor_Dim> j;
    const int size = nb_dofs / ((Tensor_Dim * (Tensor_Dim + 1)) / 2);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = getFTensorDotData<Tensor_Dim>();
      int bb = 0;
      for (; bb != size; ++bb) {
        values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
        ++field_data;
        ++base_function;
      }
      for (; bb != nb_base_functions; ++bb)
        ++base_function;
      ++values_at_gauss_pts;
    }

    MoFEMFunctionReturn(0);
  }
};

template <>
template <>
inline auto
OpCalculateTensor2SymmetricFieldValuesDot<3>::getFTensorDotData<3>() {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>(
      &dotVector[0], &dotVector[1], &dotVector[2], &dotVector[3], &dotVector[4],
      &dotVector[5]);
}

// GET GRADIENTS AT GAUSS POINTS

/**
 * \brief Evaluate field gradient values for scalar field, i.e. gradient is
 * tensor rank 1 (vector) \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim, class T, class L, class A>
struct OpCalculateScalarFieldGradient_General
    : public OpCalculateVectorFieldValues_General<Tensor_Dim, T, L, A> {

  OpCalculateScalarFieldGradient_General(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldValues_General<Tensor_Dim, T, L, A>(
            field_name, data_ptr, zero_type) {}
};

/** \brief Evaluate field gradient values for scalar field, i.e. gradient is
 * tensor rank 1 (vector), specialization \ingroup
 * mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateScalarFieldGradient_General<Tensor_Dim, double,
                                              ublas::row_major, DoubleAllocator>
    : public OpCalculateVectorFieldValues_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  OpCalculateScalarFieldGradient_General(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldValues_General<Tensor_Dim, double,
                                             ublas::row_major, DoubleAllocator>(
            field_name, data_ptr, zero_type) {}

  /**
   * \brief calculate gradient values of scalar field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Member function specialization calculating scalar field gradients for
 * tenor field rank 1 \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
MoFEMErrorCode OpCalculateScalarFieldGradient_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  if (!this->dataPtr) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Data pointer not allocated");
  }
  if (!nb_dofs && type == this->zeroType) {
    this->dataPtr->resize(Tensor_Dim, 0, false);
    MoFEMFunctionReturnHot(0);
  }
  if (!nb_dofs) {
    MoFEMFunctionReturnHot(0);
  }
  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();
  ublas::matrix<double, ublas::row_major, DoubleAllocator> &mat =
      *this->dataPtr;
  if (type == this->zeroType) {
    mat.resize(Tensor_Dim, nb_gauss_pts, false);
    mat.clear();
  }
  auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim>();
  auto gradients_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);
  FTensor::Index<'I', Tensor_Dim> I;
  for (int gg = 0; gg < nb_gauss_pts; ++gg) {
    auto field_data = data.getFTensor0FieldData();
    int bb = 0;
    for (; bb < nb_dofs; ++bb) {
      gradients_at_gauss_pts(I) += field_data * diff_base_function(I);
      ++field_data;
      ++diff_base_function;
    }
    // Number of dofs can be smaller than number of base functions
    for (; bb != nb_base_functions; ++bb)
      ++diff_base_function;
    ++gradients_at_gauss_pts;
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for scalar filed rank 0, i.e.
 * vector field \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateScalarFieldGradient
    : public OpCalculateScalarFieldGradient_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  OpCalculateScalarFieldGradient(const std::string &field_name,
                                 boost::shared_ptr<MatrixDouble> &data_ptr,
                                 const EntityType zero_type = MBVERTEX)
      : OpCalculateScalarFieldGradient_General<
            Tensor_Dim, double, ublas::row_major, DoubleAllocator>(
            field_name, data_ptr, zero_type) {}
};

/**
 * \brief Evaluate field gradient values for vector field, i.e. gradient is
 * tensor rank 2 \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
struct OpCalculateVectorFieldGradient_General
    : public OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T,
                                                   L, A> {

  OpCalculateVectorFieldGradient_General(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T, L,
                                              A>(field_name, data_ptr,
                                                 zero_type) {}
};

template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradient_General<Tensor_Dim0, Tensor_Dim1, double,
                                              ublas::row_major, DoubleAllocator>
    : public OpCalculateTensor2FieldValues_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  OpCalculateVectorFieldGradient_General(
      const std::string &field_name, boost::shared_ptr<MatrixDouble> &data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, double,
                                              ublas::row_major,
                                              DoubleAllocator>(
            field_name, data_ptr, zero_type) {}

  /**
   * \brief calculate values of vector field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/**
 * \brief Member function specialization calculating vector field gradients for
 * tenor field rank 2 \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
MoFEMErrorCode OpCalculateVectorFieldGradient_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_dofs = data.getFieldData().size();
  if (!nb_dofs && type == this->zeroType) {
    this->dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, 0, false);
    MoFEMFunctionReturnHot(0);
  }
  if (!nb_dofs) {
    MoFEMFunctionReturnHot(0);
  }
  const int nb_gauss_pts = data.getN().size1();
  const int nb_base_functions = data.getN().size2();
  ublas::matrix<double, ublas::row_major, DoubleAllocator> &mat =
      *this->dataPtr;
  if (type == this->zeroType) {
    mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
    mat.clear();
  }
  auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim1>();
  auto gradients_at_gauss_pts =
      getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
  FTensor::Index<'I', Tensor_Dim0> I;
  FTensor::Index<'J', Tensor_Dim1> J;
  int size = nb_dofs / Tensor_Dim0;
  if (nb_dofs % Tensor_Dim0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
  }
  for (int gg = 0; gg < nb_gauss_pts; ++gg) {
    auto field_data = data.getFTensor1FieldData<Tensor_Dim0>();
    int bb = 0;
    for (; bb < size; ++bb) {
      gradients_at_gauss_pts(I, J) += field_data(I) * diff_base_function(J);
      ++field_data;
      ++diff_base_function;
    }
    // Number of dofs can be smaller than number of Tensor_Dim0 x base functions
    for (; bb != nb_base_functions; ++bb)
      ++diff_base_function;
    ++gradients_at_gauss_pts;
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for scalar filed rank 0, i.e.
 * vector field \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradient
    : public OpCalculateVectorFieldGradient_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  OpCalculateVectorFieldGradient(const std::string &field_name,
                                 boost::shared_ptr<MatrixDouble> &data_ptr,
                                 const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldGradient_General<Tensor_Dim0, Tensor_Dim1, double,
                                               ublas::row_major,
                                               DoubleAllocator>(
            field_name, data_ptr, zero_type) {}
};

/** \brief Get field gradients time derivative at integration pts for scalar
 * filed rank 0, i.e. vector field \ingroup
 * mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradientDot
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> &dataPtr; ///< Data computed into this matrix
  EntityType zeroAtType;  ///< Zero values at Gauss point at this type
  VectorDouble dotVector; ///< Keeps temoorary values of time directives

  template <int Dim> inline auto getFTensorDotData() {
    static_assert(Dim || !Dim, "not implemented");
  }

  OpCalculateVectorFieldGradientDot(const std::string field_name,
                                    boost::shared_ptr<MatrixDouble> &data_ptr,
                                    const EntityType zero_at_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroAtType(zero_at_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;

    const auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    if (!nb_dofs && type == zeroAtType) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, 0, false);
      MoFEMFunctionReturnHot(0);
    }
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    dotVector.resize(nb_dofs, false);
    const double *array;
    CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;
    CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

    const int nb_gauss_pts = data.getN().size1();
    const int nb_base_functions = data.getN().size2();
    ublas::matrix<double, ublas::row_major, DoubleAllocator> &mat = *dataPtr;
    if (type == zeroAtType) {
      mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
      mat.clear();
    }
    auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim1>();
    auto gradients_at_gauss_pts =
        getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
    FTensor::Index<'I', Tensor_Dim0> I;
    FTensor::Index<'J', Tensor_Dim1> J;
    int size = nb_dofs / Tensor_Dim0;
    if (nb_dofs % Tensor_Dim0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
    }
    for (int gg = 0; gg < nb_gauss_pts; ++gg) {
      auto field_data = getFTensorDotData<3>();
      int bb = 0;
      for (; bb < size; ++bb) {
        gradients_at_gauss_pts(I, J) += field_data(I) * diff_base_function(J);
        ++field_data;
        ++diff_base_function;
      }
      // Number of dofs can be smaller than number of Tensor_Dim0 x base
      // functions
      for (; bb != nb_base_functions; ++bb)
        ++diff_base_function;
      ++gradients_at_gauss_pts;
    }
    MoFEMFunctionReturn(0);
  }
};

template <>
template <>
inline auto OpCalculateVectorFieldGradientDot<3, 3>::getFTensorDotData<3>() {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      &dotVector[0], &dotVector[1], &dotVector[2]);
}

/** \brief Get vector field for H-div approximation
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, class T, class L, class A>
struct OpCalculateHdivVectorField_General
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHdivVectorField_General(
      const std::string &field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> &data_ptr,
      const EntityType zero_type = MBTRI, const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(0) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief calculate values of vector field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

template <int Tensor_Dim, class T, class L, class A>
MoFEMErrorCode OpCalculateHdivVectorField_General<Tensor_Dim, T, L, A>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  SETERRQ2(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
           "Not implemented for T = %s and dim = %d",
           typeid(T).name(), // boost::core::demangle(typeid(T).name()),
           Tensor_Dim);
  MoFEMFunctionReturnHot(0);
}

/** \brief Get vector field for H-div approximation
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateHdivVectorField_General<Tensor_Dim, double, ublas::row_major,
                                          DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHdivVectorField_General(const std::string &field_name,
                                     boost::shared_ptr<MatrixDouble> &data_ptr,
                                     const EntityType zero_type = MBTRI,
                                     const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief Calculate values of vector field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

template <int Tensor_Dim>
MoFEMErrorCode OpCalculateHdivVectorField_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const int nb_integration_points = getGaussPts().size2();
  if (type == zeroType && side == zeroSide) {
    dataPtr->resize(Tensor_Dim, nb_integration_points, false);
    dataPtr->clear();
  }
  const int nb_dofs = data.getFieldData().size();
  if (!nb_dofs)
    MoFEMFunctionReturnHot(0);
  const int nb_base_functions = data.getN().size2() / Tensor_Dim;
  FTensor::Index<'i', Tensor_Dim> i;
  auto t_n_hdiv = data.getFTensor1N<Tensor_Dim>();
  auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
  for (int gg = 0; gg != nb_integration_points; ++gg) {
    auto t_dof = data.getFTensor0FieldData();
    int bb = 0;
    for (; bb != nb_dofs; ++bb) {
      t_data(i) += t_n_hdiv(i) * t_dof;
      ++t_n_hdiv;
      ++t_dof;
    }
    for (; bb != nb_base_functions; ++bb)
      ++t_n_hdiv;
    ++t_data;
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get vector field for H-div approximation
 * \ingroup mofem_forces_and_sources_user_data_operators
 * \note Not tested
 * \FIXME Test this
 */
template <int Tensor_Dim>
struct OpCalculateHdivVectorField
    : public OpCalculateHdivVectorField_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  OpCalculateHdivVectorField(const std::string &field_name,
                             boost::shared_ptr<MatrixDouble> &data_ptr,
                             const EntityType zero_type = MBTRI,
                             const int zero_side = 0)
      : OpCalculateHdivVectorField_General<Tensor_Dim, double, ublas::row_major,
                                           DoubleAllocator>(
            field_name, data_ptr, zero_type, zero_side) {}
};

/**
 * @brief Calculate divergence of vector field
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim dimension of space
 */
template <int Tensor_Dim>
struct OpCalculateHdivVectorDivergence
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHdivVectorDivergence(const std::string &field_name,
                                  boost::shared_ptr<VectorDouble> &data_ptr,
                                  const EntityType zero_type = MBTRI,
                                  const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(nb_integration_points, false);
      dataPtr->clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const int nb_base_functions = data.getN().size2() / Tensor_Dim;
    FTensor::Index<'i', Tensor_Dim> i;
    auto t_n_diff_hdiv = data.getFTensor2DiffN<Tensor_Dim, Tensor_Dim>();
    auto t_data = getFTensor0FromVec(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {
        double div = t_n_diff_hdiv(i, i);
        t_data += t_dof * div;
        ++t_n_diff_hdiv;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_diff_hdiv;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Calculate curl   of vector field
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim dimension of space
 */
template <int Tensor_Dim>
struct OpCalculateHcurlVectorCurl
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHcurlVectorCurl(const std::string &field_name,
                             boost::shared_ptr<MatrixDouble> &data_ptr,
                             const EntityType zero_type = MBEDGE,
                             const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim, nb_integration_points, false);
      dataPtr->clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    MatrixDouble curl_mat;

    const int nb_base_functions = data.getN().size2() / Tensor_Dim;
    FTensor::Index<'i', Tensor_Dim> i;
    FTensor::Index<'j', Tensor_Dim> j;
    FTensor::Index<'k', Tensor_Dim> k;
    auto t_n_diff_hcurl = data.getFTensor2DiffN<Tensor_Dim, Tensor_Dim>();
    auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {

      // // get curl of base functions
      // CHKERR getCurlOfHCurlBaseFunctions(side, type, data, gg, curl_mat);
      // FTensor::Tensor1<double *, 3> t_base_curl(
      //     &curl_mat(0, HVEC0), &curl_mat(0, HVEC1), &curl_mat(0, HVEC2), 3);

      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {

        // FTensor::Tensor1<double, 3> t_test;
        // t_test(k) = levi_civita(j, i, k) * t_n_diff_hcurl(i, j);
        // t_base_curl(k) -= t_test(k);
        // std::cerr << "error: " << t_base_curl(0) << " " << t_base_curl(1) <<
        // " "
        //           << t_base_curl(2) << std::endl;
        // ++t_base_curl;

        t_data(k) += t_dof * (levi_civita(j, i, k) * t_n_diff_hcurl(i, j));
        ++t_n_diff_hcurl;
        ++t_dof;
      }

      for (; bb != nb_base_functions; ++bb)
        ++t_n_diff_hcurl;
      ++t_data;
    }

    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Calculate tenor field using vectorial base, i.e. Hdiv/Hcurl
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim0 rank of the filed
 * @tparam Tensor_Dim1 dimension of space
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateHVecTensorField
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHVecTensorField(const std::string &field_name,
                             boost::shared_ptr<MatrixDouble> &data_ptr,
                             const EntityType zero_type = MBTRI,
                             const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, nb_integration_points, false);
      dataPtr->clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const int nb_base_functions = data.getN().size2() / Tensor_Dim1;
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    auto t_n_hvec = data.getFTensor1N<Tensor_Dim1>();
    auto t_data = getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor1FieldData<Tensor_Dim0>();
      int bb = 0;
      for (; bb != nb_dofs / Tensor_Dim0; ++bb) {
        t_data(i, j) += t_dof(i) * t_n_hvec(j);
        ++t_n_hvec;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_hvec;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Calculate tenor field using vectorial base, i.e. Hdiv/Hcurl
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim0 rank of the filed
 * @tparam Tensor_Dim1 dimension of space
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateHTensorTensorField
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHTensorTensorField(const std::string &field_name,
                                boost::shared_ptr<MatrixDouble> &data_ptr,
                                const EntityType zero_type = MBTRI,
                                const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, nb_integration_points, false);
      dataPtr->clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const int nb_base_functions =
        data.getN().size2() / (Tensor_Dim0 * Tensor_Dim1);
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    auto t_n_hten = data.getFTensor2N<Tensor_Dim0, Tensor_Dim1>();
    auto t_data = getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_data(i, j) += t_dof * t_n_hten(i, j);
        ++t_n_hten;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_hten;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Calculate divergence of tonsorial field using vectorial base
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim0 rank of the field
 * @tparam Tensor_Dim1 dimension of space
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateHVecTensorDivergence
    : public ForcesAndSourcesCore::UserDataOperator {

  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;

  OpCalculateHVecTensorDivergence(const std::string &field_name,
                                  boost::shared_ptr<MatrixDouble> &data_ptr,
                                  const EntityType zero_type = MBTRI,
                                  const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == 0) {
      dataPtr->resize(Tensor_Dim0, nb_integration_points, false);
      dataPtr->clear();
    }
    const int nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const int nb_base_functions = data.getN().size2() / Tensor_Dim1;
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    auto t_n_diff_hvec = data.getFTensor2DiffN<Tensor_Dim1, Tensor_Dim1>();
    auto t_data = getFTensor1FromMat<Tensor_Dim0>(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor1FieldData<Tensor_Dim0>();
      int bb = 0;
      for (; bb != nb_dofs / Tensor_Dim0; ++bb) {
        double div = t_n_diff_hvec(j, j);
        t_data(i) += t_dof(i) * div;
        ++t_n_diff_hvec;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_diff_hvec;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }
};

// Face opeartors 

/** \brief Calculate jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

  \ingroup mofem_forces_and_sources_tri_element

*/
struct OpCalculateJacForFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {
  MatrixDouble &jac;

  OpCalculateJacForFace(MatrixDouble &jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(H1), jac(jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

  \ingroup mofem_forces_and_sources_tri_element

*/
struct OpCalculateInvJacForFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {
  MatrixDouble &invJac;

  OpCalculateInvJacForFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(H1),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Transform local reference derivatives of shape functions to global
derivatives

\ingroup mofem_forces_and_sources_tri_element

*/
struct OpSetInvJacH1ForFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetInvJacH1ForFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(H1),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  MatrixDouble &invJac;
  MatrixDouble diffNinvJac;
};

/**
 * \brief brief Transform local reference derivatives of shape function to
 global derivatives for face

 * \ingroup mofem_forces_and_sources_tri_element
 */
struct OpSetInvJacHcurlFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetInvJacHcurlFace(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(HCURL),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  MatrixDouble &invJac;
  MatrixDouble diffHcurlInvJac;
};

/**
 * @brief Make Hdiv space from Hcurl space in 2d
 * @ingroup mofem_forces_and_sources_tri_element
 */
struct OpMakeHdivFromHcurl
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpMakeHdivFromHcurl()
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(HCURL) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Apply contravariant (Piola) transfer to Hdiv space on face
 * 
 * \note Hdiv space is generated by Hcurl space in 2d.
 * 
 * Contravariant Piola transformation
 * \f[
 * \psi_i|_t = \frac{1}{\textrm{det}(J)}J_{ij}\hat{\psi}_j\\
 * \left.\frac{\partial \psi_i}{\partial \xi_j}\right|_t
 * =
 * \frac{1}{\textrm{det}(J)}J_{ik}\frac{\partial \hat{\psi}_k}{\partial \xi_j}
 * \f]
 * 
 * \ingroup mofem_forces_and_sources
 *
  */
struct OpSetContravariantPiolaTransformFace
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  FTensor::Tensor2<double *, 2, 2> tJac;
  FTensor::Index<'i', 2> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  OpSetContravariantPiolaTransformFace(MatrixDouble &jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(HCURL),
        tJac(

            &jac(0, 0), &jac(0, 1), &jac(1, 0), &jac(1, 1)

        ) {}

  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

// Edge

struct OpSetContrariantPiolaTransformOnEdge
    : public EdgeElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetContrariantPiolaTransformOnEdge()
      : EdgeElementForcesAndSourcesCoreBase::UserDataOperator(HCURL) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

} // namespace MoFEM

#endif // __USER_DATA_OPERATORS_HPP__

/**
 * \defgroup mofem_forces_and_sources_user_data_operators Users Operators
 *
 * \brief Classes and functions used to evaluate fields at integration pts,
 *jacobians, etc..
 *
 * \ingroup mofem_forces_and_sources
 **/
