/**
 * @file Schur.hpp
 * @brief Assemble Schur complement
 * @date 2023-02-02
 *
 * @copyright Copyright (c) 2023
 *
 * To create nested system of Schur complements, you push sequence of operator,
 * to set up data on entities, and then assemble complements.
 *
 */

#ifndef __SCHUR_HPP__
#define __SCHUR_HPP__

namespace MoFEM {

struct SchurComplement {
  boost::shared_ptr<std::vector<MatrixDouble>> colFieldMat;
};

template <>
MoFEMErrorCode
MatSetValues<SchurComplement<0>>(Mat M,
                                 const EntitiesFieldData::EntData &row_data,
                                 const EntitiesFieldData::EntData &col_data,
                                 const double *ptr, InsertMode iora);

// /**
//  * @brief Assemble cand clear
//  *
//  */
// struct OpSchurAssemble : public ForcesAndSourcesCore::UserDataOperator {

//   OpSchurAssemble()
//       : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE) {}

// private:
// };


// template <>
// MoFEMErrorCode
// MatSetValues<OpSchurComplement<N>>(Mat M,
//                                    const EntitiesFieldData::EntData &row_data,
//                                    const EntitiesFieldData::EntData &col_data,
//                                    const double *ptr, InsertMode iora);

// template <>
// MoFEMErrorCode
// MatSetValues<OpSchurAssemble>(Mat M, const EntitiesFieldData::EntData &row_data,
//                               const EntitiesFieldData::EntData &col_data,
//                               const double *ptr, InsertMode iora);

} // namespace MoFEM

#endif //__SCHUR_HPP__