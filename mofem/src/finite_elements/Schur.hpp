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
struct OpSchurAssembleBegin : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleBegin()
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};

/**
 * @brief Assemble cand clear
 *
 */
struct OpSchurAssembleEnd : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleEnd(std::vector<std::string> fields_name,
                     std::vector<boost::shared_ptr<Range>> field_ents,
                     std::vector<SmartPetscObj<AO>> sequence_of_aos,
                     std::vector<SmartPetscObj<Mat>> sequence_of_mats,
                     std::vector<bool> sym_schur)
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPSPACE),
        fieldsName(fields_name), fieldEnts(field_ents),
        sequenceOfAOs(sequence_of_aos), sequenceOfMats(sequence_of_mats),
        symSchur(sym_schur) {}

protected:
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

  std::vector<std::string> fieldsName;
  std::vector<boost::shared_ptr<Range>> fieldEnts;
  std::vector<SmartPetscObj<AO>> sequenceOfAOs;
  std::vector<SmartPetscObj<Mat>> sequenceOfMats;
  std::vector<bool> symSchur;

private:

  MatrixDouble invMat;
  MatrixDouble invDiagOffMat;
  MatrixDouble offMatInvDiagOffMat;
  MatrixDouble transOffMatInvDiagOffMat;


};

struct SchurL2Mats : public boost::enable_shared_from_this<SchurL2Mats> {

  SchurL2Mats(const size_t idx, const UId uid_row, const UId uid_col);
  virtual ~SchurL2Mats() = default;

  const size_t iDX;
  const UId uidRow;
  const UId uidCol;

  inline auto &getMat() const { return locMats[iDX]; }
  inline auto &getRowInd() const { return rowIndices[iDX]; }
  inline auto &getColInd() const { return colIndices[iDX]; }

  static MoFEMErrorCode MatSetValues(Mat M,
                                     const EntitiesFieldData::EntData &row_data,
                                     const EntitiesFieldData::EntData &col_data,
                                     const MatrixDouble &mat, InsertMode iora);

private:
  friend OpSchurAssembleBegin;
  friend OpSchurAssembleEnd;

  struct idx_mi_tag {};
  struct uid_mi_tag {};
  struct row_mi_tag {};
  struct col_mi_tag {};

  using SchurL2Storage = multi_index_container<
      SchurL2Mats,
      indexed_by<

          ordered_unique<tag<idx_mi_tag>,
                         member<SchurL2Mats, const size_t, &SchurL2Mats::iDX>>,

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