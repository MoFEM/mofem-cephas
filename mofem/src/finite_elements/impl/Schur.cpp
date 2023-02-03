/**
 * @file Schur.cpp
 * @brief Implementation of Schur Complement
 * @date 2023-02-03
 * 
 * @copyright Copyright (c) 2023
 * 
 */

namespace MoFEM {

template <>
MoFEMErrorCode
MatSetValues<SchurComplement<0>>(Mat M,
                                   const EntitiesFieldData::EntData &row_data,
                                   const EntitiesFieldData::EntData &col_data,
                                   const double *ptr, InsertMode iora) {
  MoFEMFunctionBegin;

  if (!row_data.getFieldEntities().empty()) {
    if (auto e_ptr = row_data.getFieldEntities()[0]) {
      e_ptr->getWeakStoragePtr() =
    }
      //   if (auto e_ptr = data.getFieldEntities()[0]) {
      //     if (auto stored_data_ptr = e_ptr->getWeakStoragePtr().lock()) {
      //       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
      //               "Expected no storage on entity. Slot is needed to store
      //               schur " "data");
      //     }

      //   }
    }

  MoFEMFunctionReturn(0);
}
}