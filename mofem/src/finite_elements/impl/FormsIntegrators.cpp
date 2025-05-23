/** \file FormsIntegrators.cpp

\brief Implementation of Elements on Entities for Forces and Sources

*/



namespace MoFEM {

//! [Storage and set boundary conditions]


EssentialBcStorage::HashVectorStorage EssentialBcStorage::feStorage;

OpSetBc::OpSetBc(std::string field_name, bool yes_set,
                 boost::shared_ptr<std::vector<unsigned char>> boundary_marker)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      yesSet(yes_set), boundaryMarker(boundary_marker) {}

MoFEMErrorCode OpSetBc::doWork(int side, EntityType type,
                               EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  if (boundaryMarker) {
    if (!data.getIndices().empty())
      if (!data.getFieldEntities().empty()) {
        if (auto e_ptr = data.getFieldEntities()[0]) {
          auto indices = data.getIndices();
          for (int r = 0; r != data.getIndices().size(); ++r) {
            const auto loc_index = data.getLocalIndices()[r];
            if (loc_index >= 0) {
              if (yesSet) {
                if ((*boundaryMarker)[loc_index]) {
                  indices[r] = -1;
                }
              } else {
                if (!(*boundaryMarker)[loc_index]) {
                  indices[r] = -1;
                }
              }
            }
          }
          EssentialBcStorage::feStorage[e_ptr->getName()].push_back(
              boost::make_shared<EssentialBcStorage>(indices));
          e_ptr->getWeakStoragePtr() =
              EssentialBcStorage::feStorage[e_ptr->getName()].back();
        }
      }
  }
  MoFEMFunctionReturn(0);
}

OpUnSetBc::OpUnSetBc(std::string field_name)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW) {}

MoFEMErrorCode OpUnSetBc::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  EssentialBcStorage::feStorage[rowFieldName].clear();
  MoFEMFunctionReturn(0);
}

/**
 * @brief Set values to vector in operator
 *
 * @param V
 * @param data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <>
MoFEMErrorCode
VecSetValues<EssentialBcStorage>(Vec V,
                                 const EntitiesFieldData::EntData &data,
                                 const double *ptr, InsertMode iora) {
  MoFEMFunctionBegin;

  if(!V)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "Pointer to PETSc vector is null");

  CHKERR VecSetOption(V, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
  if (!data.getFieldEntities().empty()) {
    if (auto e_ptr = data.getFieldEntities()[0]) {
      if (auto stored_data_ptr =
              e_ptr->getSharedStoragePtr<EssentialBcStorage>()) {
        return ::VecSetValues(V, stored_data_ptr->entityIndices.size(),
                              &*stored_data_ptr->entityIndices.begin(), ptr,
                              iora);
      }
    }
  }
  return ::VecSetValues(V, data.getIndices().size(),
                        &*data.getIndices().begin(), ptr, iora);
  MoFEMFunctionReturn(0);
}

/**
 * @brief Set values to matrix in operator
 *
 * @param M
 * @param row_data
 * @param col_data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <>
MoFEMErrorCode MatSetValues<EssentialBcStorage>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const double *ptr,
    InsertMode iora) {

  if(!M)
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "Pointer to PETSc matrix is null");

  if (!row_data.getFieldEntities().empty()) {
    if (auto e_ptr = row_data.getFieldEntities()[0]) {
      if (auto stored_data_ptr =
              e_ptr->getSharedStoragePtr<EssentialBcStorage>()) {
        return ::MatSetValues(M, stored_data_ptr->entityIndices.size(),
                              &*stored_data_ptr->entityIndices.begin(),
                              col_data.getIndices().size(),
                              &*col_data.getIndices().begin(), ptr, iora);
      }
    }
  }
  return ::MatSetValues(
      M, row_data.getIndices().size(), &*row_data.getIndices().begin(),
      col_data.getIndices().size(), &*col_data.getIndices().begin(), ptr, iora);
}

//! [Storage and set boundary conditions


}