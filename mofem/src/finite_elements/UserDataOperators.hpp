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

/** \name Get values at Gauss pts */

/**@{*/

/** \name Scalar values */

/**@{*/

/** \brief Scalar field values at integration points
 *
 */
template <class T, class A>
struct OpCalculateScalarFieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateScalarFieldValues_General(
      const std::string field_name,
      boost::shared_ptr<ublas::vector<T, A>> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  OpCalculateScalarFieldValues_General(
      const std::string field_name,
      boost::shared_ptr<ublas::vector<T, A>> data_ptr,
      SmartPetscObj<Vec> data_vec, const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(field_name, OPROW),
        dataPtr(data_ptr), zeroType(zero_type), dataVec(data_vec) {
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

protected:
  boost::shared_ptr<ublas::vector<T, A>> dataPtr;
  const EntityHandle zeroType;
  SmartPetscObj<Vec> dataVec;
  VectorDouble dotVector;
};

/**
 * \brief Specialization of member function
 *
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
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
struct OpCalculateScalarFieldValues
    : public OpCalculateScalarFieldValues_General<double, DoubleAllocator> {

  using OpCalculateScalarFieldValues_General<
      double, DoubleAllocator>::OpCalculateScalarFieldValues_General;

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
    const size_t nb_gauss_pts = getGaussPts().size2();
    if (type == zeroType) {
      vec.resize(nb_gauss_pts, false);
      vec.clear();
    }

    const size_t nb_dofs = data.getFieldData().size();

    if (nb_dofs) {

      if (dataVec.use_count()) {
        dotVector.resize(nb_dofs, false);
        const double *array;
        CHKERR VecGetArrayRead(dataVec, &array);
        const auto &local_indices = data.getLocalIndices();
        for (int i = 0; i != local_indices.size(); ++i)
          if (local_indices[i] != -1)
            dotVector[i] = array[local_indices[i]];
          else
            dotVector[i] = 0;
        CHKERR VecRestoreArrayRead(dataVec, &array);
        data.getFieldData().swap(dotVector);
      }

      const size_t nb_base_functions = data.getN().size2();
      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts = getFTensor0FromVec(vec);
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor0FieldData();
        size_t bb = 0;
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

      if (dataVec.use_count()) {
        data.getFieldData().swap(dotVector);
      }
    }

    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Get rate of scalar field at integration points
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
struct OpCalculateScalarFieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateScalarFieldValuesDot(const std::string field_name,
                                  boost::shared_ptr<VectorDouble> data_ptr,
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

    const size_t nb_gauss_pts = getGaussPts().size2();
    VectorDouble &vec = *dataPtr;
    if (type == zeroAtType) {
      vec.resize(nb_gauss_pts, false);
      vec.clear();
    }

    auto &local_indices = data.getLocalIndices();
    const size_t nb_dofs = local_indices.size();
    if (nb_dofs) {

      std::array<double, MAX_DOFS_ON_ENTITY> dot_dofs_vector;
      const double *array;
      CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
      for (size_t i = 0; i != local_indices.size(); ++i)
        if (local_indices[i] != -1)
          dot_dofs_vector[i] = array[local_indices[i]];
        else
          dot_dofs_vector[i] = 0;
      CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

      const size_t nb_base_functions = data.getN().size2();
      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts = getFTensor0FromVec(vec);

      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        size_t bb = 0;
        for (; bb != nb_dofs; ++bb) {
          values_at_gauss_pts += dot_dofs_vector[bb] * base_function;
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

private:
  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroAtType;
};

/**
 * \depreacted Name inconstent with other operators
 *
 */
using OpCalculateScalarValuesDot = OpCalculateScalarFieldValuesDot;

/**@}*/

/** \name Vector field values at integration points */

/**@{*/

/** \brief Calculate field values for tenor field rank 1, i.e. vector field
 *
 */
template <int Tensor_Dim, class T, class L, class A>
struct OpCalculateVectorFieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateVectorFieldValues_General(
      const std::string field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> data_ptr,
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

protected:
  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;
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
 * rank 1, i.e. vector field
 *
 */
template <int Tensor_Dim>
struct OpCalculateVectorFieldValues_General<Tensor_Dim, double,
                                            ublas::row_major, DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateVectorFieldValues_General(const std::string field_name,
                                       boost::shared_ptr<MatrixDouble> data_ptr,
                                       const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
};

/**
 * \brief Member function specialization calculating values for tenor field rank
 *
 */
template <int Tensor_Dim>
MoFEMErrorCode OpCalculateVectorFieldValues_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const size_t nb_gauss_pts = getGaussPts().size2();
  auto &mat = *dataPtr;
  if (type == zeroType) {
    mat.resize(Tensor_Dim, nb_gauss_pts, false);
    mat.clear();
  }

  const size_t nb_dofs = data.getFieldData().size();
  if (nb_dofs) {

    if (nb_gauss_pts) {
      const size_t nb_base_functions = data.getN().size2();
      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);
      FTensor::Index<'I', Tensor_Dim> I;
      const size_t size = nb_dofs / Tensor_Dim;
      if (nb_dofs % Tensor_Dim) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor1FieldData<Tensor_Dim>();
        size_t bb = 0;
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
    }
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get values at integration pts for tensor filed rank 1, i.e. vector
 * field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateVectorFieldValues
    : public OpCalculateVectorFieldValues_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  using OpCalculateVectorFieldValues_General<
      Tensor_Dim, double, ublas::row_major,
      DoubleAllocator>::OpCalculateVectorFieldValues_General;
};

/**@}*/

/** \name Vector field values at integration points */

/**@{*/

/** \brief Calculate divergence of vector field values for tenor field rank 1,
 * i.e. vector field
 *
 */
template <int Tensor_Dim, class T, class L, class A>
struct OpCalculateDivergenceVectorFieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateDivergenceVectorFieldValues_General(
      const std::string field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  /**
   * \brief calculate divergence values of vector field at integration points
   * @param  side side entity number
   * @param  type side entity type
   * @param  data entity data
   * @return      error code
   */
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;
};

template <int Tensor_Dim, class T, class L, class A>
MoFEMErrorCode
OpCalculateDivergenceVectorFieldValues_General<Tensor_Dim, T, L, A>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;
  SETERRQ2(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
           "Not implemented for T = %s and dim = %d",
           typeid(T).name(), // boost::core::demangle(typeid(T).name()),
           Tensor_Dim);
  MoFEMFunctionReturnHot(0);
}

/** \brief Calculate field values (template specialization) for tensor field
 * rank 1, i.e. vector field
 *
 */
template <int Tensor_Dim>
struct OpCalculateDivergenceVectorFieldValues_General<
    Tensor_Dim, double, ublas::row_major, DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateDivergenceVectorFieldValues_General(
      const std::string field_name, boost::shared_ptr<VectorDouble> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroType;
};

/**
 * \brief Member function specialization calculating values for tenor field rank
 *
 */
template <int Tensor_Dim>
MoFEMErrorCode OpCalculateDivergenceVectorFieldValues_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const size_t nb_gauss_pts = getGaussPts().size2();
  auto &vec = *dataPtr;
  if (type == zeroType) {
    vec.resize(nb_gauss_pts, false);
    vec.clear();
  }

  const size_t nb_dofs = data.getFieldData().size();
  if (nb_dofs) {

    if (nb_gauss_pts) {
      const size_t nb_base_functions = data.getN().size2();
      auto values_at_gauss_pts = getFTensor0FromVec(vec);
      FTensor::Index<'I', Tensor_Dim> I;
      const size_t size = nb_dofs / Tensor_Dim;
      if (nb_dofs % Tensor_Dim) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim>();
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor1FieldData<Tensor_Dim>();
        size_t bb = 0;
        for (; bb != size; ++bb) {
          values_at_gauss_pts += field_data(I) * diff_base_function(I);
          ++field_data;
          ++diff_base_function;
        }
        // Number of dofs can be smaller than number of Tensor_Dim x base
        // functions
        for (; bb != nb_base_functions; ++bb)
          ++diff_base_function;
        ++values_at_gauss_pts;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get values at integration pts for tensor filed rank 1, i.e. vector
 * field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateDivergenceVectorFieldValues
    : public OpCalculateDivergenceVectorFieldValues_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  using OpCalculateDivergenceVectorFieldValues_General<
      Tensor_Dim, double, ublas::row_major,
      DoubleAllocator>::OpCalculateDivergenceVectorFieldValues_General;
};

/** \brief Approximate field valuse for given petsc vector
 *
 * \note Look at PetscData to see what vectors could be extarcted with that user
 * data opetaor.
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim, PetscData::DataContext CTX>
struct OpCalculateVectorFieldValuesFromPetscVecImpl
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateVectorFieldValuesFromPetscVecImpl(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
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

    const double *array;

    auto get_array = [&](const auto ctx, auto vec) {
      MoFEMFunctionBegin;
      if ((getFEMethod()->data_ctx & ctx).none())
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Vector not set");
      CHKERR VecGetArrayRead(vec, &array);
      MoFEMFunctionReturn(0);
    };

    auto restore_array = [&](auto vec) {
      return VecRestoreArrayRead(vec, &array);
    };

    switch (CTX) {
    case PetscData::CTX_SET_X:
      CHKERR get_array(PetscData::CtxSetX, getFEMethod()->ts_u);
      break;
    case PetscData::CTX_SET_X_T:
      CHKERR get_array(PetscData::CtxSetX_T, getFEMethod()->ts_u_t);
      break;
    case PetscData::CTX_SET_X_TT:
      CHKERR get_array(PetscData::CtxSetX_TT, getFEMethod()->ts_u_tt);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "That case is not implemented");
    }

    dotVector.resize(local_indices.size());
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;

    switch (CTX) {
    case PetscData::CTX_SET_X:
      CHKERR restore_array(getFEMethod()->ts_u);
      break;
    case PetscData::CTX_SET_X_T:
      CHKERR restore_array(getFEMethod()->ts_u_t);
      break;
    case PetscData::CTX_SET_X_TT:
      CHKERR restore_array(getFEMethod()->ts_u_tt);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
              "That case is not implemented");
    }

    const size_t nb_gauss_pts = data.getN().size1();
    const size_t nb_base_functions = data.getN().size2();
    MatrixDouble &mat = *dataPtr;
    if (type == zeroAtType) {
      mat.resize(Tensor_Dim, nb_gauss_pts, false);
      mat.clear();
    }
    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);
    FTensor::Index<'I', Tensor_Dim> I;
    const size_t size = nb_dofs / Tensor_Dim;
    if (nb_dofs % Tensor_Dim) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Data inconsistency");
    }
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = getFTensor1FromArray<Tensor_Dim, Tensor_Dim>(dotVector);
      size_t bb = 0;
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

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroAtType;
  VectorDouble dotVector;
};

/** \brief Get time direvatives of values at integration pts for tensor filed
 * rank 1, i.e. vector field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
using OpCalculateVectorFieldValuesDot =
    OpCalculateVectorFieldValuesFromPetscVecImpl<Tensor_Dim,
                                                 PetscData::CTX_SET_X_T>;

/** \brief Get second time direvatives of values at integration pts for tensor
 * filed rank 1, i.e. vector field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
using OpCalculateVectorFieldValuesDotDot =
    OpCalculateVectorFieldValuesFromPetscVecImpl<Tensor_Dim,
                                                 PetscData::CTX_SET_X_TT>;

/**@}*/

/** \name Tensor field values at integration points */

/**@{*/

/** \brief Calculate field values for tenor field rank 2.
 *
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
struct OpCalculateTensor2FieldValues_General
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateTensor2FieldValues_General(
      const std::string field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;
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

  OpCalculateTensor2FieldValues_General(
      const std::string field_name,
      boost::shared_ptr<
          ublas::matrix<double, ublas::row_major, DoubleAllocator>>
          data_ptr,
      const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  OpCalculateTensor2FieldValues_General(
      const std::string field_name,
      boost::shared_ptr<
          ublas::matrix<double, ublas::row_major, DoubleAllocator>>
          data_ptr,
      SmartPetscObj<Vec> data_vec, const EntityType zero_type = MBVERTEX)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), dataVec(data_vec) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  boost::shared_ptr<ublas::matrix<double, ublas::row_major, DoubleAllocator>>
      dataPtr;
  const EntityHandle zeroType;
  SmartPetscObj<Vec> dataVec;
  VectorDouble dotVector;
};

template <int Tensor_Dim0, int Tensor_Dim1>
MoFEMErrorCode OpCalculateTensor2FieldValues_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  MatrixDouble &mat = *dataPtr;
  const size_t nb_gauss_pts = data.getN().size1();
  if (type == zeroType) {
    mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
    mat.clear();
  }

  const size_t nb_dofs = data.getFieldData().size();

  if (dataVec.use_count()) {
    dotVector.resize(nb_dofs, false);
    const double *array;
    CHKERR VecGetArrayRead(dataVec, &array);
    const auto &local_indices = data.getLocalIndices();
    for (int i = 0; i != local_indices.size(); ++i)
      if (local_indices[i] != -1)
        dotVector[i] = array[local_indices[i]];
      else
        dotVector[i] = 0;
    CHKERR VecRestoreArrayRead(dataVec, &array);
    data.getFieldData().swap(dotVector);
  }

  if (nb_dofs) {
    const size_t nb_base_functions = data.getN().size2();
    auto base_function = data.getFTensor0N();
    auto values_at_gauss_pts =
        getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    const size_t size = nb_dofs / (Tensor_Dim0 * Tensor_Dim1);
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      auto field_data = data.getFTensor2FieldData<Tensor_Dim0, Tensor_Dim1>();
      size_t bb = 0;
      for (; bb != size; ++bb) {
        values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
        ++field_data;
        ++base_function;
      }
      for (; bb != nb_base_functions; ++bb)
        ++base_function;
      ++values_at_gauss_pts;
    }

    if (dataVec.use_count()) {
      data.getFieldData().swap(dotVector);
    }
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get values at integration pts for tensor filed rank 2, i.e. matrix
 * field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2FieldValues
    : public OpCalculateTensor2FieldValues_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  using OpCalculateTensor2FieldValues_General<
      Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
      DoubleAllocator>::OpCalculateTensor2FieldValues_General;
};

/** \brief Get time direvarive values at integration pts for tensor filed rank
 * 2, i.e. matrix field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2FieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateTensor2FieldValuesDot(const std::string field_name,
                                   boost::shared_ptr<MatrixDouble> data_ptr,
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

    const size_t nb_gauss_pts = getGaussPts().size2();
    MatrixDouble &mat = *dataPtr;
    if (type == zeroAtType) {
      mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
      mat.clear();
    }
    const auto &local_indices = data.getLocalIndices();
    const size_t nb_dofs = local_indices.size();
    if (nb_dofs) {
      dotVector.resize(nb_dofs, false);
      const double *array;
      CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
      for (size_t i = 0; i != local_indices.size(); ++i)
        if (local_indices[i] != -1)
          dotVector[i] = array[local_indices[i]];
        else
          dotVector[i] = 0;
      CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

      const size_t nb_base_functions = data.getN().size2();

      auto base_function = data.getFTensor0N();
      auto values_at_gauss_pts =
          getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
      FTensor::Index<'i', Tensor_Dim0> i;
      FTensor::Index<'j', Tensor_Dim1> j;
      const size_t size = nb_dofs / (Tensor_Dim0 * Tensor_Dim1);
      for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
        auto field_data = getFTensorDotData<Tensor_Dim0, Tensor_Dim1>();
        size_t bb = 0;
        for (; bb != size; ++bb) {
          values_at_gauss_pts(i, j) += field_data(i, j) * base_function;
          ++field_data;
          ++base_function;
        }
        for (; bb != nb_base_functions; ++bb)
          ++base_function;
        ++values_at_gauss_pts;
      }
    }
    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<MatrixDouble> dataPtr; ///< Data computed into this matrix
  EntityType zeroAtType;  ///< Zero values at Gauss point at this type
  VectorDouble dotVector; ///< Keeps temoorary values of time directives

  template <int Dim0, int Dim1> auto getFTensorDotData() {
    static_assert(Dim0 || !Dim0 || Dim1 || !Dim1, "not implemented");
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

/**
 * @brief Calculate symmetric tensor field values at integration pts.
 *
 * @tparam Tensor_Dim

 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateTensor2SymmetricFieldValues
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateTensor2SymmetricFieldValues(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBEDGE, const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  OpCalculateTensor2SymmetricFieldValues(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      SmartPetscObj<Vec> data_vec, const EntityType zero_type = MBEDGE,
      const int zero_side = 0)
      : ForcesAndSourcesCore::UserDataOperator(
            field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
        dataPtr(data_ptr), zeroType(zero_type), zeroSide(zero_side),
        dataVec(data_vec) {
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
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    if (dataVec.use_count()) {
      dotVector.resize(nb_dofs, false);
      const double *array;
      CHKERR VecGetArrayRead(dataVec, &array);
      const auto &local_indices = data.getLocalIndices();
      for (int i = 0; i != local_indices.size(); ++i)
        if (local_indices[i] != -1)
          dotVector[i] = array[local_indices[i]];
        else
          dotVector[i] = 0;
      CHKERR VecRestoreArrayRead(dataVec, &array);
      data.getFieldData().swap(dotVector);
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

    if (dataVec.use_count()) {
      data.getFieldData().swap(dotVector);
    }

    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
  SmartPetscObj<Vec> dataVec;
  VectorDouble dotVector;
};

/**
 * @brief Calculate symmetric tensor field rates ant integratio pts.
 *
 * @tparam Tensor_Dim
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateTensor2SymmetricFieldValuesDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateTensor2SymmetricFieldValuesDot(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBEDGE, const int zero_side = 0)
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

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
  VectorDouble dotVector;

  template <int Dim> inline auto getFTensorDotData() {
    static_assert(Dim || !Dim, "not implemented");
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

template <>
template <>
inline auto
OpCalculateTensor2SymmetricFieldValuesDot<2>::getFTensorDotData<2>() {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>(
      &dotVector[0], &dotVector[1], &dotVector[2]);
}

/**@}*/

/** \name Gradients of scalar fields at integration points */

/**@{*/

/**
 * \brief Evaluate field gradient values for scalar field, i.e. gradient is
 * tensor rank 1 (vector)
 *
 */
template <int Tensor_Dim, class T, class L, class A>
struct OpCalculateScalarFieldGradient_General
    : public OpCalculateVectorFieldValues_General<Tensor_Dim, T, L, A> {

  OpCalculateScalarFieldGradient_General(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldValues_General<Tensor_Dim, T, L, A>(
            field_name, data_ptr, zero_type) {}
};

/** \brief Evaluate field gradient values for scalar field, i.e. gradient is
 * tensor rank 1 (vector), specialization
 *
 */
template <int Tensor_Dim>
struct OpCalculateScalarFieldGradient_General<Tensor_Dim, double,
                                              ublas::row_major, DoubleAllocator>
    : public OpCalculateVectorFieldValues_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {

  OpCalculateScalarFieldGradient_General(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateVectorFieldValues_General<Tensor_Dim, double,
                                             ublas::row_major, DoubleAllocator>(
            field_name, data_ptr, zero_type) {
    if (!this->dataPtr)
      THROW_MESSAGE("Data pointer not allocated");
  }

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
 * tenor field rank 1
 *
 */
template <int Tensor_Dim>
MoFEMErrorCode OpCalculateScalarFieldGradient_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  const size_t nb_gauss_pts = this->getGaussPts().size2();
  auto &mat = *this->dataPtr;
  if (type == this->zeroType) {
    mat.resize(Tensor_Dim, nb_gauss_pts, false);
    mat.clear();
  }

  const int nb_dofs = data.getFieldData().size();
  if (nb_dofs) {

    if (nb_gauss_pts) {
      const int nb_base_functions = data.getN().size2();
      auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim>();
      auto gradients_at_gauss_pts = getFTensor1FromMat<Tensor_Dim>(mat);

      FTensor::Index<'I', Tensor_Dim> I;
      for (int gg = 0; gg < nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor0FieldData();
        int bb = 0;
        for (; bb != nb_dofs; ++bb) {
          gradients_at_gauss_pts(I) += field_data * diff_base_function(I);
          ++field_data;
          ++diff_base_function;
        }
        // Number of dofs can be smaller than number of base functions
        for (; bb < nb_base_functions; ++bb)
          ++diff_base_function;
        ++gradients_at_gauss_pts;
      }
    }
  }

  MoFEMFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for scalar filed rank 0, i.e.
 * vector field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateScalarFieldGradient
    : public OpCalculateScalarFieldGradient_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {
  using OpCalculateScalarFieldGradient_General<
      Tensor_Dim, double, ublas::row_major,
      DoubleAllocator>::OpCalculateScalarFieldGradient_General;
};

/**}*/

/** \name Gradients of tensor fields at integration points */

/**@{*/

/**
 * \brief Evaluate field gradient values for vector field, i.e. gradient is
 * tensor rank 2
 *
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
struct OpCalculateVectorFieldGradient_General
    : public OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T,
                                                   L, A> {

  using OpCalculateTensor2FieldValues_General<
      Tensor_Dim0, Tensor_Dim1, T, L, A>::OpCalculateTensor2FieldValues_General;
};

template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradient_General<Tensor_Dim0, Tensor_Dim1, double,
                                              ublas::row_major, DoubleAllocator>
    : public OpCalculateTensor2FieldValues_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  using OpCalculateTensor2FieldValues_General<
      Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
      DoubleAllocator>::OpCalculateTensor2FieldValues_General;

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
 * tenor field rank 2
 *
 */
template <int Tensor_Dim0, int Tensor_Dim1>
MoFEMErrorCode OpCalculateVectorFieldGradient_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  if (!this->dataPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Data pointer not allocated");

  const size_t nb_gauss_pts = this->getGaussPts().size2();
  auto &mat = *this->dataPtr;
  if (type == this->zeroType) {
    mat.resize(Tensor_Dim0 * Tensor_Dim1, nb_gauss_pts, false);
    mat.clear();
  }

  if (nb_gauss_pts) {
    const size_t nb_dofs = data.getFieldData().size();

    if (nb_dofs) {

      if (this->dataVec.use_count()) {
        this->dotVector.resize(nb_dofs, false);
        const double *array;
        CHKERR VecGetArrayRead(this->dataVec, &array);
        const auto &local_indices = data.getLocalIndices();
        for (int i = 0; i != local_indices.size(); ++i)
          if (local_indices[i] != -1)
            this->dotVector[i] = array[local_indices[i]];
          else
            this->dotVector[i] = 0;
        CHKERR VecRestoreArrayRead(this->dataVec, &array);
        data.getFieldData().swap(this->dotVector);
      }

      const int nb_base_functions = data.getN().size2();
      auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim1>();
      auto gradients_at_gauss_pts =
          getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(mat);
      FTensor::Index<'I', Tensor_Dim0> I;
      FTensor::Index<'J', Tensor_Dim1> J;
      int size = nb_dofs / Tensor_Dim0;
      if (nb_dofs % Tensor_Dim0) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      for (int gg = 0; gg < nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor1FieldData<Tensor_Dim0>();
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

      if (this->dataVec.use_count()) {
        data.getFieldData().swap(this->dotVector);
      }
    }
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for scalar filed rank 0, i.e.
 * vector field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradient
    : public OpCalculateVectorFieldGradient_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  using OpCalculateVectorFieldGradient_General<
      Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
      DoubleAllocator>::OpCalculateVectorFieldGradient_General;
};

/** \brief Get field gradients time derivative at integration pts for scalar
 * filed rank 0, i.e. vector field
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateVectorFieldGradientDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateVectorFieldGradientDot(const std::string field_name,
                                    boost::shared_ptr<MatrixDouble> data_ptr,
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
      auto field_data = getFTensorDotData<Tensor_Dim0>();
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

private:
  boost::shared_ptr<MatrixDouble> dataPtr; ///< Data computed into this matrix
  EntityType zeroAtType;  ///< Zero values at Gauss point at this type
  VectorDouble dotVector; ///< Keeps temoorary values of time directives

  template <int Dim> inline auto getFTensorDotData() {
    static_assert(Dim || !Dim, "not implemented");
  }
};

template <>
template <>
inline auto OpCalculateVectorFieldGradientDot<3, 3>::getFTensorDotData<3>() {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(
      &dotVector[0], &dotVector[1], &dotVector[2]);
}

template <>
template <>
inline auto OpCalculateVectorFieldGradientDot<2, 2>::getFTensorDotData<2>() {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>(&dotVector[0],
                                                            &dotVector[1]);
}

/**
 * \brief Evaluate field gradient values for symmetric 2nd order tensor field,
 * i.e. gradient is tensor rank 3
 *
 */
template <int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
struct OpCalculateTensor2SymmetricFieldGradient_General
    : public OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T,
                                                   L, A> {

  OpCalculateTensor2SymmetricFieldGradient_General(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateTensor2FieldValues_General<Tensor_Dim0, Tensor_Dim1, T, L,
                                              A>(field_name, data_ptr,
                                                 zero_type) {}
};

template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2SymmetricFieldGradient_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator>
    : public OpCalculateTensor2FieldValues_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  OpCalculateTensor2SymmetricFieldGradient_General(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
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
 * \brief Member function specialization calculating tensor field gradients for
 * symmetric tensor field rank 2
 *
 */
template <int Tensor_Dim0, int Tensor_Dim1>
MoFEMErrorCode OpCalculateTensor2SymmetricFieldGradient_General<
    Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  if (!this->dataPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Data pointer not allocated");

  const size_t nb_gauss_pts = this->getGaussPts().size2();
  constexpr size_t msize = (Tensor_Dim0 * (Tensor_Dim0 + 1)) / 2;
  auto &mat = *this->dataPtr;
  if (type == this->zeroType) {
    mat.resize(msize * Tensor_Dim1, nb_gauss_pts, false);
    mat.clear();
  }

  if (nb_gauss_pts) {
    const size_t nb_dofs = data.getFieldData().size();

    if (nb_dofs) {

      const int nb_base_functions = data.getN().size2();
      auto diff_base_function = data.getFTensor1DiffN<Tensor_Dim1>();
      auto gradients_at_gauss_pts =
          getFTensor3DgFromMat<Tensor_Dim0, Tensor_Dim1>(mat);
      FTensor::Index<'I', Tensor_Dim0> I;
      FTensor::Index<'J', Tensor_Dim0> J;
      FTensor::Index<'K', Tensor_Dim1> K;
      int size = nb_dofs / msize;
      if (nb_dofs % msize) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Data inconsistency");
      }
      for (int gg = 0; gg < nb_gauss_pts; ++gg) {
        auto field_data = data.getFTensor2SymmetricFieldData<Tensor_Dim0>();
        int bb = 0;
        for (; bb < size; ++bb) {
          gradients_at_gauss_pts(I, J, K) +=
              field_data(I, J) * diff_base_function(K);
          ++field_data;
          ++diff_base_function;
        }
        // Number of dofs can be smaller than number of Tensor_Dim0 x base
        // functions
        for (; bb != nb_base_functions; ++bb)
          ++diff_base_function;
        ++gradients_at_gauss_pts;
      }
    }
  }
  MoFEMFunctionReturn(0);
}

/** \brief Get field gradients at integration pts for symmetric tensorial field
 * rank 2
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, int Tensor_Dim1>
struct OpCalculateTensor2SymmetricFieldGradient
    : public OpCalculateTensor2SymmetricFieldGradient_General<
          Tensor_Dim0, Tensor_Dim1, double, ublas::row_major, DoubleAllocator> {

  OpCalculateTensor2SymmetricFieldGradient(
      const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
      const EntityType zero_type = MBVERTEX)
      : OpCalculateTensor2SymmetricFieldGradient_General<
            Tensor_Dim0, Tensor_Dim1, double, ublas::row_major,
            DoubleAllocator>(field_name, data_ptr, zero_type) {}
};

/**@}*/

/** \name Transform tensors and vectors */

/**@{*/

/**
 * @brief Calculate \f$ \pmb\sigma_{ij} = \mathbf{D}_{ijkl} \pmb\varepsilon_{kl}
 * \f$
 *
 * @tparam DIM
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int DIM_01, int DIM_23, int S = 0>
struct OpTensorTimesSymmetricTensor
    : public ForcesAndSourcesCore::UserDataOperator {

  using EntData = DataForcesAndSourcesCore::EntData;
  using UserOp = ForcesAndSourcesCore::UserDataOperator;

  OpTensorTimesSymmetricTensor(const std::string field_name,
                               boost::shared_ptr<MatrixDouble> in_mat,
                               boost::shared_ptr<MatrixDouble> out_mat,
                               boost::shared_ptr<MatrixDouble> d_mat)
      : UserOp(field_name, OPROW), inMat(in_mat), outMat(out_mat), dMat(d_mat) {
    // Only is run for vertices
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
    if (!inMat)
      THROW_MESSAGE("Pointer for in mat is null");
    if (!outMat)
      THROW_MESSAGE("Pointer for out mat is null");
    if (!dMat)
      THROW_MESSAGE("Pointer for tensor mat is null");
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = getGaussPts().size2();
    auto t_D = getFTensor4DdgFromMat<DIM_01, DIM_23, S>(*(dMat));
    auto t_in = getFTensor2SymmetricFromMat<DIM_01>(*(inMat));
    outMat->resize((DIM_23 * (DIM_23 + 1)) / 2, nb_gauss_pts, false);
    auto t_out = getFTensor2SymmetricFromMat<DIM_23>(*(outMat));
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      t_out(i, j) = t_D(i, j, k, l) * t_in(k, l);
      ++t_in;
      ++t_out;
    }
    MoFEMFunctionReturn(0);
  }

private:
  FTensor::Index<'i', DIM_01> i;
  FTensor::Index<'j', DIM_01> j;
  FTensor::Index<'k', DIM_23> k;
  FTensor::Index<'l', DIM_23> l;

  boost::shared_ptr<MatrixDouble> inMat;
  boost::shared_ptr<MatrixDouble> outMat;
  boost::shared_ptr<MatrixDouble> dMat;
};

template <int DIM>
struct OpSymmetrizeTensor : public ForcesAndSourcesCore::UserDataOperator {

  using EntData = DataForcesAndSourcesCore::EntData;
  using UserOp = ForcesAndSourcesCore::UserDataOperator;

  OpSymmetrizeTensor(const std::string field_name,
                     boost::shared_ptr<MatrixDouble> in_mat,
                     boost::shared_ptr<MatrixDouble> out_mat)
      : UserOp(field_name, OPROW), inMat(in_mat), outMat(out_mat) {
    // Only is run for vertices
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
    if (!inMat)
      THROW_MESSAGE("Pointer not set for in matrix");
    if (!outMat)
      THROW_MESSAGE("Pointer not set for in matrix");
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
    const size_t nb_gauss_pts = getGaussPts().size2();
    auto t_in = getFTensor2FromMat<DIM, DIM>(*(inMat));
    outMat->resize((DIM * (DIM + 1)) / 2, nb_gauss_pts, false);
    auto t_out = getFTensor2SymmetricFromMat<DIM>(*(outMat));
    for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
      t_out(i, j) = (t_in(i, j) || t_in(j, i)) / 2;
      ++t_in;
      ++t_out;
    }
    MoFEMFunctionReturn(0);
  }

private:
  FTensor::Index<'i', DIM> i;
  FTensor::Index<'j', DIM> j;
  boost::shared_ptr<MatrixDouble> inMat;
  boost::shared_ptr<MatrixDouble> outMat;
};

struct OpScaleMatrix : public ForcesAndSourcesCore::UserDataOperator {

  using EntData = DataForcesAndSourcesCore::EntData;
  using UserOp = ForcesAndSourcesCore::UserDataOperator;

  OpScaleMatrix(const std::string field_name, const double scale,
                boost::shared_ptr<MatrixDouble> in_mat,
                boost::shared_ptr<MatrixDouble> out_mat)
      : UserOp(field_name, OPROW), scale(scale), inMat(in_mat),
        outMat(out_mat) {
    // Only is run for vertices
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
    if (!inMat)
      THROW_MESSAGE("Pointer not set for in matrix");
    if (!outMat)
      THROW_MESSAGE("Pointer not set for in matrix");
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
    outMat->resize(inMat->size1(), inMat->size2(), false);
    noalias(*outMat) = scale * (*inMat);
    MoFEMFunctionReturn(0);
  }

private:
  const double scale;
  boost::shared_ptr<MatrixDouble> inMat;
  boost::shared_ptr<MatrixDouble> outMat;
};

/**@}*/

/** \name H-div/H-curls (Vectorial bases) values at integration points */

/**@{*/

/** \brief Get vector field for H-div approximation
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim0, class T, class L, class A>
struct OpCalculateHVecVectorField_General
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHVecVectorField_General(
      const std::string field_name,
      boost::shared_ptr<ublas::matrix<T, L, A>> data_ptr,
      const EntityType zero_type = MBEDGE, const int zero_side = 0)
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

private:
  boost::shared_ptr<ublas::matrix<T, L, A>> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

template <int Tensor_Dim, class T, class L, class A>
MoFEMErrorCode OpCalculateHVecVectorField_General<Tensor_Dim, T, L, A>::doWork(
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
struct OpCalculateHVecVectorField_General<Tensor_Dim, double, ublas::row_major,
                                          DoubleAllocator>
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHVecVectorField_General(const std::string field_name,
                                     boost::shared_ptr<MatrixDouble> data_ptr,
                                     const EntityType zero_type = MBEDGE,
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

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

template <int Tensor_Dim>
MoFEMErrorCode OpCalculateHVecVectorField_General<
    Tensor_Dim, double, ublas::row_major,
    DoubleAllocator>::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_integration_points = this->getGaussPts().size2();
  if (type == zeroType && side == zeroSide) {
    dataPtr->resize(Tensor_Dim, nb_integration_points, false);
    dataPtr->clear();
  }
  const size_t nb_dofs = data.getFieldData().size();
  if (!nb_dofs)
    MoFEMFunctionReturnHot(0);
  const size_t nb_base_functions = data.getN().size2() / Tensor_Dim;
  FTensor::Index<'i', Tensor_Dim> i;
  auto t_n_hdiv = data.getFTensor1N<Tensor_Dim>();
  auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
  for (size_t gg = 0; gg != nb_integration_points; ++gg) {
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
 *
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateHVecVectorField
    : public OpCalculateHVecVectorField_General<
          Tensor_Dim, double, ublas::row_major, DoubleAllocator> {
  using OpCalculateHVecVectorField_General<
      Tensor_Dim, double, ublas::row_major,
      DoubleAllocator>::OpCalculateHVecVectorField_General;
};

/** \brief Get vector field for H-div approximation
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template <int Tensor_Dim>
struct OpCalculateHVecVectorFieldDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHVecVectorFieldDot(const std::string field_name,
                                boost::shared_ptr<MatrixDouble> data_ptr,
                                const EntityType zero_type = MBEDGE,
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

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

template <int Tensor_Dim>
MoFEMErrorCode OpCalculateHVecVectorFieldDot<Tensor_Dim>::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_integration_points = this->getGaussPts().size2();
  if (type == zeroType && side == zeroSide) {
    dataPtr->resize(Tensor_Dim, nb_integration_points, false);
    dataPtr->clear();
  }

  auto &local_indices = data.getIndices();
  const size_t nb_dofs = local_indices.size();
  if (nb_dofs) {

    std::array<double, MAX_DOFS_ON_ENTITY> dot_dofs_vector;
    const double *array;
    CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
    for (size_t i = 0; i != nb_dofs; ++i)
      if (local_indices[i] != -1)
        dot_dofs_vector[i] = array[local_indices[i]];
      else
        dot_dofs_vector[i] = 0;
    CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

    const size_t nb_base_functions = data.getN().size2() / 3;
    FTensor::Index<'i', Tensor_Dim> i;
    auto t_n_hdiv = data.getFTensor1N<Tensor_Dim>();
    auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
    for (size_t gg = 0; gg != nb_integration_points; ++gg) {
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_data(i) += t_n_hdiv(i) * dot_dofs_vector[bb];
        ++t_n_hdiv;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_hdiv;
      ++t_data;
    }
  }
  MoFEMFunctionReturn(0);
}

/**
 * @brief Calculate divergence of vector field
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam BASE_DIM
 * @tparam SPACE_DIM
 */
template <int BASE_DIM, int SPACE_DIM>
struct OpCalculateHdivVectorDivergence
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHdivVectorDivergence(const std::string field_name,
                                  boost::shared_ptr<VectorDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const size_t nb_base_functions = data.getN().size2() / BASE_DIM;
    FTensor::Index<'i', BASE_DIM> i;
    FTensor::Index<'j', SPACE_DIM> j;
    auto t_n_diff_hdiv = data.getFTensor2DiffN<BASE_DIM, SPACE_DIM>();
    auto t_data = getFTensor0FromVec(*dataPtr);
    for (size_t gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_data += t_dof * t_n_diff_hdiv(j, j);
        ++t_n_diff_hdiv;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_n_diff_hdiv;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

/**
 * @brief Calculate gradient of vector field
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam BASE_DIM
 * @tparam SPACE_DIM
 */
template <int BASE_DIM, int SPACE_DIM>
struct OpCalculateHVecVectorGradient
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHVecVectorGradient(const std::string field_name,
                                boost::shared_ptr<MatrixDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(BASE_DIM * SPACE_DIM, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const size_t nb_base_functions = data.getN().size2() / BASE_DIM;
    FTensor::Index<'i', BASE_DIM> i;
    FTensor::Index<'j', SPACE_DIM> j;
    auto t_base_diff = data.getFTensor2DiffN<BASE_DIM, SPACE_DIM>();
    auto t_data = getFTensor2FromMat<BASE_DIM, SPACE_DIM>(*dataPtr);
    for (size_t gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_data(i, j) += t_dof * t_base_diff(i, j);
        ++t_base_diff;
        ++t_dof;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_base_diff;
      ++t_data;
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

/**
 * @brief Calculate divergence of vector field dot
 * @ingroup mofem_forces_and_sources_user_data_operators
 *
 * @tparam Tensor_Dim dimension of space
 */
template <int Tensor_Dim1, int Tensor_Dim2>
struct OpCalculateHdivVectorDivergenceDot
    : public ForcesAndSourcesCore::UserDataOperator {

  OpCalculateHdivVectorDivergenceDot(const std::string field_name,
                                     boost::shared_ptr<VectorDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(nb_integration_points, false);
      dataPtr->clear();
    }

    const auto &local_indices = data.getLocalIndices();
    const int nb_dofs = local_indices.size();
    if (nb_dofs) {

      std::array<double, MAX_DOFS_ON_ENTITY> dot_dofs_vector;
      const double *array;
      CHKERR VecGetArrayRead(getFEMethod()->ts_u_t, &array);
      for (size_t i = 0; i != local_indices.size(); ++i)
        if (local_indices[i] != -1)
          dot_dofs_vector[i] = array[local_indices[i]];
        else
          dot_dofs_vector[i] = 0;
      CHKERR VecRestoreArrayRead(getFEMethod()->ts_u_t, &array);

      const size_t nb_base_functions = data.getN().size2() / Tensor_Dim1;
      FTensor::Index<'i', Tensor_Dim1> i;
      auto t_n_diff_hdiv = data.getFTensor2DiffN<Tensor_Dim1, Tensor_Dim2>();
      auto t_data = getFTensor0FromVec(*dataPtr);
      for (size_t gg = 0; gg != nb_integration_points; ++gg) {
        int bb = 0;
        for (; bb != nb_dofs; ++bb) {
          double div = 0;
          for (auto ii = 0; ii != Tensor_Dim2; ++ii)
            div += t_n_diff_hdiv(ii, ii);
          t_data += dot_dofs_vector[bb] * div;
          ++t_n_diff_hdiv;
        }
        for (; bb != nb_base_functions; ++bb)
          ++t_n_diff_hdiv;
        ++t_data;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
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

  OpCalculateHcurlVectorCurl(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);

    MatrixDouble curl_mat;

    const size_t nb_base_functions = data.getN().size2() / Tensor_Dim;
    FTensor::Index<'i', Tensor_Dim> i;
    FTensor::Index<'j', Tensor_Dim> j;
    FTensor::Index<'k', Tensor_Dim> k;
    auto t_n_diff_hcurl = data.getFTensor2DiffN<Tensor_Dim, Tensor_Dim>();
    auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
    for (int gg = 0; gg != nb_integration_points; ++gg) {

      auto t_dof = data.getFTensor0FieldData();
      int bb = 0;
      for (; bb != nb_dofs; ++bb) {

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

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
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

  OpCalculateHVecTensorField(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (nb_dofs) {
      const size_t nb_base_functions = data.getN().size2() / 3;
      FTensor::Index<'i', Tensor_Dim0> i;
      FTensor::Index<'j', Tensor_Dim1> j;
      auto t_n_hvec = data.getFTensor1N<3>();
      auto t_data = getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(*dataPtr);
      for (size_t gg = 0; gg != nb_integration_points; ++gg) {
        auto t_dof = data.getFTensor1FieldData<Tensor_Dim0>();
        size_t bb = 0;
        for (; bb != nb_dofs / Tensor_Dim0; ++bb) {
          t_data(i, j) += t_dof(i) * t_n_hvec(j);
          ++t_n_hvec;
          ++t_dof;
        }
        for (; bb != nb_base_functions; ++bb)
          ++t_n_hvec;
        ++t_data;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
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

  OpCalculateHTensorTensorField(const std::string field_name,
                                boost::shared_ptr<MatrixDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == zeroSide) {
      dataPtr->resize(Tensor_Dim0 * Tensor_Dim1, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (!nb_dofs)
      MoFEMFunctionReturnHot(0);
    const size_t nb_base_functions =
        data.getN().size2() / (Tensor_Dim0 * Tensor_Dim1);
    FTensor::Index<'i', Tensor_Dim0> i;
    FTensor::Index<'j', Tensor_Dim1> j;
    auto t_n_hten = data.getFTensor2N<Tensor_Dim0, Tensor_Dim1>();
    auto t_data = getFTensor2FromMat<Tensor_Dim0, Tensor_Dim1>(*dataPtr);
    for (size_t gg = 0; gg != nb_integration_points; ++gg) {
      auto t_dof = data.getFTensor0FieldData();
      size_t bb = 0;
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

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
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

  OpCalculateHVecTensorDivergence(const std::string field_name,
                                  boost::shared_ptr<MatrixDouble> data_ptr,
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
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == 0) {
      dataPtr->resize(Tensor_Dim0, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (nb_dofs) {
      const size_t nb_base_functions = data.getN().size2() / 3;
      FTensor::Index<'i', Tensor_Dim0> i;
      FTensor::Index<'j', Tensor_Dim1> j;
      auto t_n_diff_hvec = data.getFTensor2DiffN<3, Tensor_Dim1>();
      auto t_data = getFTensor1FromMat<Tensor_Dim0>(*dataPtr);
      for (size_t gg = 0; gg != nb_integration_points; ++gg) {
        auto t_dof = data.getFTensor1FieldData<Tensor_Dim0>();
        size_t bb = 0;
        for (; bb != nb_dofs / Tensor_Dim0; ++bb) {
          double div = t_n_diff_hvec(j, j);
          t_data(i) += t_dof(i) * div;
          ++t_n_diff_hvec;
          ++t_dof;
        }
        for (; bb < nb_base_functions; ++bb)
          ++t_n_diff_hvec;
        ++t_data;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
};

/**
 * @brief Calculate trace of vector (Hdiv/Hcurl) space
 *
 * @tparam Tensor_Dim
 * @tparam OpBase
 */
template <int Tensor_Dim, typename OpBase>
struct OpCalculateHVecTensorTrace : public OpBase {

  OpCalculateHVecTensorTrace(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> data_ptr,
                             const EntityType zero_type = MBEDGE,
                             const int zero_side = 0)
      : OpBase(field_name, OpBase::OPROW), dataPtr(data_ptr),
        zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const size_t nb_integration_points = OpBase::getGaussPts().size2();
    if (type == zeroType && side == 0) {
      dataPtr->resize(Tensor_Dim, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (nb_dofs) {
      auto t_normal = OpBase::getFTensor1Normal();
      t_normal(i) /= sqrt(t_normal(j) * t_normal(j));
      const size_t nb_base_functions = data.getN().size2() / 3;
      auto t_base = data.getFTensor1N<3>();
      auto t_data = getFTensor1FromMat<Tensor_Dim>(*dataPtr);
      for (size_t gg = 0; gg != nb_integration_points; ++gg) {
        auto t_dof = data.getFTensor1FieldData<Tensor_Dim>();
        size_t bb = 0;
        for (; bb != nb_dofs / Tensor_Dim; ++bb) {
          t_data(i) += t_dof(i) * (t_base(j) * t_normal(j));
          ++t_base;
          ++t_dof;
        }
        for (; bb < nb_base_functions; ++bb)
          ++t_base;
        ++t_data;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
  FTensor::Index<'i', Tensor_Dim> i;
  FTensor::Index<'j', Tensor_Dim> j;
};

template <>
struct OpCalculateHVecTensorTrace<
    3, FaceElementForcesAndSourcesCore::UserDataOperator>
    : public FaceElementForcesAndSourcesCore::UserDataOperator {

  OpCalculateHVecTensorTrace(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> data_ptr,
                             const EntityType zero_type = MBEDGE,
                             const int zero_side = 0)
      : UserDataOperator(field_name, OPROW), dataPtr(data_ptr),
        zeroType(zero_type), zeroSide(zero_side) {
    if (!dataPtr)
      THROW_MESSAGE("Pointer is not set");
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const size_t nb_integration_points = getGaussPts().size2();
    if (type == zeroType && side == 0) {
      dataPtr->resize(3, nb_integration_points, false);
      dataPtr->clear();
    }
    const size_t nb_dofs = data.getFieldData().size();
    if (nb_dofs) {
      auto t_normal = getFTensor1Normal();
      t_normal(i) /= sqrt(t_normal(j) * t_normal(j));
      const size_t nb_base_functions = data.getN().size2() / 3;
      auto t_base = data.getFTensor1N<3>();
      auto t_data = getFTensor1FromMat<3>(*dataPtr);
      for (size_t gg = 0; gg != nb_integration_points; ++gg) {
        auto t_dof = data.getFTensor1FieldData<3>();
        if (getNormalsAtGaussPts().size1() == nb_integration_points) {
          VectorDouble n = getNormalsAtGaussPts(gg);
          auto t_n = getFTensor1FromPtr<3>(&*n.data().begin());
          t_normal(i) = t_n(i) / sqrt(t_n(j) * t_n(j));
        }
        size_t bb = 0;
        for (; bb != nb_dofs / 3; ++bb) {
          t_data(i) += t_dof(i) * (t_base(j) * t_normal(j));
          ++t_base;
          ++t_dof;
        }
        for (; bb < nb_base_functions; ++bb)
          ++t_base;
        ++t_data;
      }
    }
    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> dataPtr;
  const EntityHandle zeroType;
  const int zeroSide;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
};

/**@}*/

/** \name Other operators */

/**@{*/

/**@}*/

inline auto getFaceJac(MatrixDouble &jac, const FTensor::Number<2> &) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 2, 2>{
      &jac(0, 0), &jac(1, 0), &jac(2, 0), &jac(3, 0)};
}

inline auto getFaceJac(MatrixDouble &jac, const FTensor::Number<3> &) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>{
      &jac(0, 0), &jac(1, 0), &jac(4, 0), &jac(2, 0), &jac(3, 0),
      &jac(5, 0), &jac(6, 0), &jac(7, 0), &jac(8, 0)};
}

/** \name Operators for faces */

/**@{*/

/** \brief Calculate jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

  \ingroup mofem_forces_and_sources_tri_element

*/
template <int DIM> struct OpCalculateJacForFaceImpl;

template <>
struct OpCalculateJacForFaceImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpCalculateJacForFaceImpl(MatrixDouble &jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE),
        jac(jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  MatrixDouble &jac;
};

template <>
struct OpCalculateJacForFaceImpl<3> : public OpCalculateJacForFaceImpl<2> {

  using OpCalculateJacForFaceImpl<2>::OpCalculateJacForFaceImpl;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

using OpCalculateJacForFace = OpCalculateJacForFaceImpl<2>;

using OpCalculateJacForFaceEmbeddedIn3DSpace = OpCalculateJacForFaceImpl<3>;

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

  \ingroup mofem_forces_and_sources_tri_element

*/
template <int DIM> struct OpCalculateInvJacForFaceImpl;

template <>
struct OpCalculateInvJacForFaceImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpCalculateInvJacForFaceImpl(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(NOSPACE),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  MatrixDouble &invJac;
};

template <>
struct OpCalculateInvJacForFaceImpl<3>
    : public OpCalculateInvJacForFaceImpl<2> {

  using OpCalculateInvJacForFaceImpl<2>::OpCalculateInvJacForFaceImpl;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

using OpCalculateInvJacForFace = OpCalculateInvJacForFaceImpl<2>;

using OpCalculateInvJacForFaceEmbeddedIn3DSpace =
    OpCalculateInvJacForFaceImpl<3>;

/** \brief Transform local reference derivatives of shape functions to global
derivatives

\ingroup mofem_forces_and_sources_tri_element

*/
template <int DIM> struct OpSetInvJacSpaceForFaceImpl;

template <>
struct OpSetInvJacSpaceForFaceImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetInvJacSpaceForFaceImpl(MatrixDouble &inv_jac, FieldSpace space)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(space),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  MatrixDouble &invJac;
  MatrixDouble diffNinvJac;
};

template <>
struct OpSetInvJacSpaceForFaceImpl<3> : public OpSetInvJacSpaceForFaceImpl<2> {

  using OpSetInvJacSpaceForFaceImpl<2>::OpSetInvJacSpaceForFaceImpl;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

struct OpSetInvJacH1ForFace : public OpSetInvJacSpaceForFaceImpl<2> {
  OpSetInvJacH1ForFace(MatrixDouble &inv_jac)
      : OpSetInvJacSpaceForFaceImpl(inv_jac, H1) {}
};

struct OpSetInvJacL2ForFace : public OpSetInvJacSpaceForFaceImpl<2> {
  OpSetInvJacL2ForFace(MatrixDouble &inv_jac)
      : OpSetInvJacSpaceForFaceImpl(inv_jac, L2) {}
};

struct OpSetInvJacH1ForFaceEmbeddedIn3DSpace
    : public OpSetInvJacSpaceForFaceImpl<3> {
  OpSetInvJacH1ForFaceEmbeddedIn3DSpace(MatrixDouble &inv_jac)
      : OpSetInvJacSpaceForFaceImpl(inv_jac, H1) {}
};

struct OpSetInvJacL2ForFaceEmbeddedIn3DSpace
    : public OpSetInvJacSpaceForFaceImpl<3> {
  OpSetInvJacL2ForFaceEmbeddedIn3DSpace(MatrixDouble &inv_jac)
      : OpSetInvJacSpaceForFaceImpl(inv_jac, L2) {}
};

/**
 * \brief Transform local reference derivatives of shape function to
 global derivatives for face

 * \ingroup mofem_forces_and_sources_tri_element
 */
template <int DIM> struct OpSetInvJacHcurlFaceImpl;

template <>
struct OpSetInvJacHcurlFaceImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetInvJacHcurlFaceImpl(MatrixDouble &inv_jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(HCURL),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  MatrixDouble &invJac;
  MatrixDouble diffHcurlInvJac;
};

template <>
struct OpSetInvJacHcurlFaceImpl<3> : public OpSetInvJacHcurlFaceImpl<2> {
  using OpSetInvJacHcurlFaceImpl<2>::OpSetInvJacHcurlFaceImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

using OpSetInvJacHcurlFace = OpSetInvJacHcurlFaceImpl<2>;
using OpSetInvJacHcurlFaceEmbeddedIn3DSpace = OpSetInvJacHcurlFaceImpl<3>;

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
template <int DIM> struct OpSetContravariantPiolaTransformOnFace2DImpl;

template <>
struct OpSetContravariantPiolaTransformOnFace2DImpl<2>
    : public FaceElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetContravariantPiolaTransformOnFace2DImpl(MatrixDouble &jac)
      : FaceElementForcesAndSourcesCoreBase::UserDataOperator(HCURL), jAc(jac) {
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

protected:
  MatrixDouble &jAc;
  MatrixDouble piolaN;
  MatrixDouble piolaDiffN;
};

template <>
struct OpSetContravariantPiolaTransformOnFace2DImpl<3>
    : public OpSetContravariantPiolaTransformOnFace2DImpl<2> {
  using OpSetContravariantPiolaTransformOnFace2DImpl<
      2>::OpSetContravariantPiolaTransformOnFace2DImpl;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

using OpSetContravariantPiolaTransformOnFace2D =
    OpSetContravariantPiolaTransformOnFace2DImpl<2>;
using OpSetContravariantPiolaTransformOnFace2DEmbeddedIn3DSpace =
    OpSetContravariantPiolaTransformOnFace2DImpl<3>;

/**@}*/

/** \name Operators for edges */

/**@{*/

struct OpSetContravariantPiolaTransformOnEdge2D
    : public EdgeElementForcesAndSourcesCoreBase::UserDataOperator {

  OpSetContravariantPiolaTransformOnEdge2D(const FieldSpace space = HCURL)
      : EdgeElementForcesAndSourcesCoreBase::UserDataOperator(space) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  std::vector<double> l1;
};

/**
 * @deprecated Name is deprecated and this is added for back compatibility
 */
using OpSetContrariantPiolaTransformOnEdge =
    OpSetContravariantPiolaTransformOnEdge2D;

/**@}*/

/** \name Operator for fat prisms */

/**@{*/

/**
 * @brief Operator for fat prism element updating integration weights in the
 * volume.
 *
 * Jacobian on the distorted element is nonconstant. This operator updates
 * integration weight on prism to take into account nonconstat jacobian.
 *
 * \f[
 * W_i = w_i \left( \frac{1}{2V} \left\| \frac{\partial \mathbf{x}}{\partial
 * \pmb\xi} \right\| \right)
 * \f]
 * where \f$w_i\f$ is integration weight at integration point \f$i\f$,
 * \f$\mathbf{x}\f$ is physical coordinate, and \f$\pmb\xi\f$ is reference
 * element coordinate.
 *
 */
struct OpMultiplyDeterminantOfJacobianAndWeightsForFatPrisms
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  OpMultiplyDeterminantOfJacobianAndWeightsForFatPrisms()
      : FatPrismElementForcesAndSourcesCore::UserDataOperator(H1) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);
};

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  FIXME Generalize function for arbitrary face orientation in 3d space
  FIXME Calculate to Jacobins for two faces

  \ingroup mofem_forces_and_sources_prism_element

*/
struct OpCalculateInvJacForFatPrism
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  OpCalculateInvJacForFatPrism(boost::shared_ptr<MatrixDouble> inv_jac_ptr)
      : FatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJacPtr(inv_jac_ptr), invJac(*invJacPtr) {}

  OpCalculateInvJacForFatPrism(MatrixDouble &inv_jac)
      : FatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJac(inv_jac) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  const boost::shared_ptr<MatrixDouble> invJacPtr;
  MatrixDouble &invJac;
};

/** \brief Transform local reference derivatives of shape functions to global
derivatives

FIXME Generalize to curved shapes
FIXME Generalize to case that top and bottom face has different shape

\ingroup mofem_forces_and_sources_prism_element

*/
struct OpSetInvJacH1ForFatPrism
    : public FatPrismElementForcesAndSourcesCore::UserDataOperator {

  OpSetInvJacH1ForFatPrism(boost::shared_ptr<MatrixDouble> inv_jac_ptr)
      : FatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJacPtr(inv_jac_ptr), invJac(*invJacPtr) {}

  OpSetInvJacH1ForFatPrism(MatrixDouble &inv_jac)
      : FatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJac(inv_jac) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  const boost::shared_ptr<MatrixDouble> invJacPtr;
  MatrixDouble &invJac;
  MatrixDouble diffNinvJac;
};

// Flat prism

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  FIXME Generalize function for arbitrary face orientation in 3d space
  FIXME Calculate to Jacobins for two faces

  \ingroup mofem_forces_and_sources_prism_element

*/
struct OpCalculateInvJacForFlatPrism
    : public FlatPrismElementForcesAndSourcesCore::UserDataOperator {

  OpCalculateInvJacForFlatPrism(MatrixDouble &inv_jac_f3)
      : FlatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJacF3(inv_jac_f3) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  MatrixDouble &invJacF3;
};

/** \brief Transform local reference derivatives of shape functions to global
derivatives

FIXME Generalize to curved shapes
FIXME Generalize to case that top and bottom face has different shape

\ingroup mofem_forces_and_sources_prism_element

*/
struct OpSetInvJacH1ForFlatPrism
    : public FlatPrismElementForcesAndSourcesCore::UserDataOperator {

  OpSetInvJacH1ForFlatPrism(MatrixDouble &inv_jac_f3)
      : FlatPrismElementForcesAndSourcesCore::UserDataOperator(H1),
        invJacF3(inv_jac_f3) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  MatrixDouble &invJacF3;
  MatrixDouble diffNinvJac;
};

/**@}*/

/** \name Operation on matrices at integration points */

/**@{*/

template <int DIM>
struct OpInvertMatrix : public ForcesAndSourcesCore::UserDataOperator {

  OpInvertMatrix(boost::shared_ptr<MatrixDouble> in_ptr,
                 boost::shared_ptr<VectorDouble> det_ptr,
                 boost::shared_ptr<MatrixDouble> out_ptr)
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST), inPtr(in_ptr),
        outPtr(out_ptr), detPtr(det_ptr) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    return doWorkImpl(side, type, data, FTensor::Number<DIM>());
  }

private:
  boost::shared_ptr<MatrixDouble> inPtr;
  boost::shared_ptr<MatrixDouble> outPtr;
  boost::shared_ptr<VectorDouble> detPtr;

  MoFEMErrorCode doWorkImpl(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data,
                            const FTensor::Number<3> &);
};

template <int DIM>
MoFEMErrorCode
OpInvertMatrix<DIM>::doWorkImpl(int side, EntityType type,
                                DataForcesAndSourcesCore::EntData &data,
                                const FTensor::Number<3> &) {
  MoFEMFunctionBegin;

  if (!inPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for inPtr matrix not allocated");
  if (!detPtr)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Pointer for detPtr matrix not allocated");

  const auto nb_rows = inPtr->size1();
  const auto nb_integration_pts = inPtr->size2();

  // Calculate determinant
  {
    detPtr->resize(nb_integration_pts, false);
    auto t_in = getFTensor2FromMat<3, 3>(*inPtr);
    auto t_det = getFTensor0FromVec(*detPtr);
    for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
      determinantTensor3by3(t_in, t_det);
      ++t_in;
      ++t_det;
    }
  }

  // Invert jacobian
  if (outPtr) {
    outPtr->resize(nb_rows, nb_integration_pts, false);
    auto t_in = getFTensor2FromMat<3, 3>(*inPtr);
    auto t_out = getFTensor2FromMat<3, 3>(*outPtr);
    auto t_det = getFTensor0FromVec(*detPtr);
    for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
      invertTensor3by3(t_in, t_det, t_out);
      ++t_in;
      ++t_out;
      ++t_det;
    }
  }

  MoFEMFunctionReturn(0);
}

/**@}*/

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
