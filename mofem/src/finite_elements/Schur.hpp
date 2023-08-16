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

/**
 * @brief Clear Schur complement internal data
 * 
 */
struct OpSchurAssembleBegin : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleBegin()
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Assemble Schur complement (Implementation)
 *
 */
struct OpSchurAssembleEndImpl : public ForcesAndSourcesCore::UserDataOperator {

  using MatSetValuesRaw = boost::function<MoFEMErrorCode(
      Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n,
      const PetscInt idxn[], const PetscScalar v[], InsertMode addv)>;

  static MatSetValuesRaw matSetValuesSchurRaw;

  /**
   * @brief Construct a new Op Schur Assemble End object
   *
   * @param fields_name list of fields
   * @param field_ents list of entities on which schur complement is applied (can be empty)
   * @param sequence_of_aos list of maps from base problem to Schur complement matrix
   * @param sequence_of_mats list of Schur complement matrices
   * @param sym_schur true if Schur complement is symmetric
   * @param symm_op true if block diagonal is symmetric
   */
  OpSchurAssembleEndImpl(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur, bool symm_op = true);

  /**
   * @brief Construct a new Op Schur Assemble End object
   *
   * @param fields_name list of fields
   * @param field_ents list of entities on which schur complement is applied (can be empty)
   * @param sequence_of_aos list of maps from base problem to Schur complement matrix
   * @param sequence_of_mats list of Schur complement matrices
   * @param sym_schur true if Schur complement is symmetric
   * @param diag_eps add epsilon on diagonal of inverted matrix 
   * @param symm_op true if block diagonal is symmetric
   */
  OpSchurAssembleEndImpl(std::vector<std::string> fields_name,
                         std::vector<boost::shared_ptr<Range>> field_ents,
                         std::vector<SmartPetscObj<AO>> sequence_of_aos,
                         std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                         std::vector<bool> sym_schur,
                         std::vector<double> diag_eps, bool symm_op = true);

protected:

  template <typename I>
  MoFEMErrorCode doWorkImpl(int side, EntityType type,
                            EntitiesFieldData::EntData &data);

  std::vector<std::string> fieldsName;
  std::vector<boost::shared_ptr<Range>> fieldEnts;
  std::vector<SmartPetscObj<AO>> sequenceOfAOs;
  std::vector<SmartPetscObj<Mat>> sequenceOfMats;
  std::vector<bool> symSchur;
  std::vector<double> diagEps;

  MatrixDouble invMat;
  MatrixDouble invDiagOffMat;
  MatrixDouble offMatInvDiagOffMat;
  MatrixDouble transOffMatInvDiagOffMat;
};

struct SCHUR_DSYSV; ///< SY	symmetric
struct SCHUR_DGESV; ///< GE	general (i.e., nonsymmetric, in some cases
                    ///< rectangular)

/**
 * @brief Assemble Schur complement
 *
 */
template <typename I> struct OpSchurAssembleEnd;

template <>
struct OpSchurAssembleEnd<SCHUR_DSYSV> : public OpSchurAssembleEndImpl {
  using OpSchurAssembleEndImpl::OpSchurAssembleEndImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

template <>
struct OpSchurAssembleEnd<SCHUR_DGESV> : public OpSchurAssembleEndImpl {
  using OpSchurAssembleEndImpl::OpSchurAssembleEndImpl;
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Schur complement data storage
 * 
 */
struct SchurL2Mats : public boost::enable_shared_from_this<SchurL2Mats> {

  SchurL2Mats(const size_t idx, const UId uid_row, const UId uid_col);
  virtual ~SchurL2Mats() = default;

  const UId uidRow;
  const UId uidCol;

  inline auto &getMat() const { return locMats[iDX]; }
  inline auto &getRowInd() const { return rowIndices[iDX]; }
  inline auto &getColInd() const { return colIndices[iDX]; }

  using MatSetValuesPtr = boost::function<MoFEMErrorCode(
      Mat M, const EntitiesFieldData::EntData &row_data,
      const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
      InsertMode iora)>;

  static MatSetValuesPtr matSetValuesPtr; ///< backend assembly function

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);

private:
  const size_t iDX;

  friend OpSchurAssembleBegin;
  friend OpSchurAssembleEndImpl;

  struct idx_mi_tag {};
  struct uid_mi_tag {};
  struct row_mi_tag {};
  struct col_mi_tag {};

  using SchurL2Storage = multi_index_container<
      SchurL2Mats,
      indexed_by<

          ordered_unique<
              tag<uid_mi_tag>,
              composite_key<
                  SchurL2Mats,

                  member<SchurL2Mats, const UId, &SchurL2Mats::uidRow>,
                  member<SchurL2Mats, const UId, &SchurL2Mats::uidCol>

                  >>,

          ordered_non_unique<tag<row_mi_tag>, member<SchurL2Mats, const UId,
                                                     &SchurL2Mats::uidRow>>,

          ordered_non_unique<tag<col_mi_tag>, member<SchurL2Mats, const UId,
                                                     &SchurL2Mats::uidCol>>

          >>;

  static boost::ptr_vector<MatrixDouble> locMats;
  static boost::ptr_vector<VectorInt> rowIndices;
  static boost::ptr_vector<VectorInt> colIndices;
  static SchurL2Storage schurL2Storage;
};

template <>
inline MoFEMErrorCode
MatSetValues<SchurL2Mats>(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora) {
  return SchurL2Mats::MatSetValues(M, row_data, col_data, mat, iora);
}

template <>
inline MoFEMErrorCode
VecSetValues<SchurL2Mats>(Vec V, const EntitiesFieldData::EntData &data,
                          const VectorDouble &nf, InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

template <>
inline MoFEMErrorCode MatSetValues<AssemblyTypeSelector<SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<SchurL2Mats>(M, row_data, col_data, mat, iora);
}

template <>
inline MoFEMErrorCode VecSetValues<AssemblyTypeSelector<SCHUR>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return VecSetValues<SchurL2Mats>(V, data, nf, iora);
}

} // namespace MoFEM

#endif //__SCHUR_HPP__