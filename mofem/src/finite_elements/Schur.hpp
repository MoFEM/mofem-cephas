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

struct SchurEvents {
  static PetscLogEvent MOFEM_EVENT_schurL2MatsMatSetValues;
  static PetscLogEvent MOFEM_EVENT_opSchurAssembleEnd;
  static PetscLogEvent MOFEM_EVENT_diagBlockStrutureSetValues;
  static PetscLogEvent MOFEM_EVENT_diagBlockStrutureMult;
  SchurEvents();
};

/**
 * @brief Clear Schur complement internal data
 * 
 */
struct OpSchurAssembleBegin : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleBegin();

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

struct SchurL2Mats;

template <>
MoFEMErrorCode
MatSetValues<SchurL2Mats>(Mat M, const EntitiesFieldData::EntData &row_data,
                          const EntitiesFieldData::EntData &col_data,
                          const MatrixDouble &mat, InsertMode iora);

template <>
inline MoFEMErrorCode
VecSetValues<SchurL2Mats>(Vec V, const EntitiesFieldData::EntData &data,
                          const VectorDouble &nf, InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

template <>
MoFEMErrorCode MatSetValues<AssemblyTypeSelector<SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora);

template <>
inline MoFEMErrorCode VecSetValues<AssemblyTypeSelector<SCHUR>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return VecSetValues<SchurL2Mats>(V, data, nf, iora);
}

/**
 * @brief Create Schur complement
 * 
 * @param dm 
 * @param schur_is 
 * @return std::tuple<SmartPetscObj<IS>, SmartPetscObj<Mat>> 
 */
std::tuple<SmartPetscObj<IS>, SmartPetscObj<IS>> createSchurISDiff(DM dm,
                                                                   IS schur_is);
struct DiagBlockStruture;

/**
 * @brief Create a Mat Diag Blocks object
 * 
 * @return Mat 
 */
boost::shared_ptr<DiagBlockStruture> createSchurBlockMatStructure(

    DM dm,                                 //< dm
    std::vector<std::string> fields_names, //< block field
    std::vector<std::string> fe_names,     //< block fes
    std::vector<boost::shared_ptr<ForcesAndSourcesCore>>
        fe_ptrs //< block elements

);

using SchurShellMatData =
    std::pair<SmartPetscObj<Mat>, boost::shared_ptr<DiagBlockStruture>>;

/**
 * @brief Create a Schur Mat object
 * 
 * @param dm 
 * @param data 
 * @return std::pair<SmartPetscObj<Mat>, boost::shared_ptr<DiagBlockStruture>>
 */
SchurShellMatData
createSchurBlockMat(DM dm, boost::shared_ptr<DiagBlockStruture> data);

template <>
MoFEMErrorCode
MatSetValues<DiagBlockStruture>(Mat M,
                                const EntitiesFieldData::EntData &row_data,
                                const EntitiesFieldData::EntData &col_data,
                                const MatrixDouble &mat, InsertMode iora);

template <>
inline MoFEMErrorCode MatSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<DiagBlockStruture>(M, row_data, col_data, mat, iora);
}

template <>
inline MoFEMErrorCode VecSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

using SchurNestMatrixData =
    std::pair<std::array<SmartPetscObj<Mat>, 4>,
              std::array<boost::shared_ptr<DiagBlockStruture>, 4>>;

SchurNestMatrixData
getSchurNestMatArray(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                     SchurShellMatData A);

/**
 * @brief Create a Mat Diag Blocks object
 *
 * @return Mat
 */
std::pair<SmartPetscObj<Mat>, SchurNestMatrixData>
createSchurNestedMatrix(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                        SchurNestMatrixData schur_net_dat);

} // namespace MoFEM

#endif //__SCHUR_HPP__