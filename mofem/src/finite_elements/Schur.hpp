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
 * \note Try septate floating points operations from book keeping. Also, align
 * memory that blocks follow floating point operations.
 *
 */

#ifndef __SCHUR_HPP__
#define __SCHUR_HPP__

namespace MoFEM {

/**
 * @brief Structure to register events for Schur block assembly and solver
 */
struct SchurEvents {
  static PetscLogEvent MOFEM_EVENT_schurL2MatsMatSetValues;
  static PetscLogEvent MOFEM_EVENT_opSchurAssembleEnd;
  static PetscLogEvent MOFEM_EVENT_diagBlockStrutureSetValues;
  static PetscLogEvent MOFEM_EVENT_diagBlockStrutureMult;
  static PetscLogEvent MOFEM_EVENT_diagBlockStrutureSolve;
  static PetscLogEvent MOFEM_EVENT_zeroRowsAndCols;
  SchurEvents();
};

struct SchurL2Mats;
struct DiagBlockStruture;

struct OpSchurAssembleBase
    : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleBase() = delete;

  using MatSetValuesRaw = boost::function<MoFEMErrorCode(
      Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n,
      const PetscInt idxn[], const PetscScalar v[], InsertMode addv)>;
  static MatSetValuesRaw matSetValuesSchurRaw;

private:
  using UserDataOperator::UserDataOperator;
};

OpSchurAssembleBase *createOpSchurAssembleBegin();

/**
 * @brief Construct a new Op Schur Assemble End object
 *
 * @param fields_name list of fields
 * @param field_ents list of entities on which schur complement is applied
 * (can be empty)
 * @param sequence_of_aos list of maps from base problem to Schur complement
 * matrix
 * @param sequence_of_mats list of Schur complement matrices
 * @param sym_schur true if Schur complement is symmetric
 * @param symm_op true if block diagonal is symmetric
 */
OpSchurAssembleBase *createOpSchurAssembleEnd(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, bool symm_op = false,
    boost::shared_ptr<DiagBlockStruture> diag_blocks = nullptr);

/**
 * @brief Construct a new Op Schur Assemble End object
 *
 * @param fields_name list of fields
 * @param field_ents list of entities on which schur complement is applied
 * (can be empty)
 * @param sequence_of_aos list of maps from base problem to Schur complement
 * matrix
 * @param sequence_of_mats list of Schur complement matrices
 * @param sym_schur true if Schur complement is symmetric
 * @param diag_eps add epsilon on diagonal of inverted matrix
 * @param symm_op true if block diagonal is symmetric
 */
OpSchurAssembleBase *createOpSchurAssembleEnd(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, std::vector<double> diag_eps,
    bool symm_op = false,
    boost::shared_ptr<DiagBlockStruture> diag_blocks = nullptr);

struct SchurBackendMatSetValuesPtr {
  SchurBackendMatSetValuesPtr() = delete;
  using MatSetValuesPtr = boost::function<MoFEMErrorCode(
      Mat M, const EntitiesFieldData::EntData &row_data,
      const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
      InsertMode iora)>;
  static MatSetValuesPtr matSetValuesPtr;      ///< backend assembly function
  static MatSetValuesPtr matSetValuesBlockPtr; ///< backend assembly function
};

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

using SchurFieldPair = std::pair<std::string, std::string>;

using SchurFEOpsFEandFields = std::vector<

    std::pair<std::string, std::vector<SchurFieldPair>>

    >;

/**
 * @brief Create a Mat Diag Blocks object
 *
 * @return Mat
 */
boost::shared_ptr<DiagBlockStruture> createSchurBlockMatStructure(

    DM dm,                                 //< dm
    SchurFEOpsFEandFields schur_fe_op_vec //< block elements


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

/***
 * @brief Specialization of MatSetValues for DiagBlockStruture
 */
template <>
MoFEMErrorCode
MatSetValues<DiagBlockStruture>(Mat M,
                                const EntitiesFieldData::EntData &row_data,
                                const EntitiesFieldData::EntData &col_data,
                                const MatrixDouble &mat, InsertMode iora);

/***
 * @brief Specialisation of MatSetValues for AssemblyTypeSelector<BLOCK_MAT>
 */
template <>
inline MoFEMErrorCode MatSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<DiagBlockStruture>(M, row_data, col_data, mat, iora);
}

/***
 * @brief Specialisation of VecSetValues for AssemblyTypeSelector<BLOCK_MAT>
 */
template <>
inline MoFEMErrorCode VecSetValues<AssemblyTypeSelector<BLOCK_MAT>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

using SchurNestMatrixData =
    std::pair<std::array<SmartPetscObj<Mat>, 4>,
              std::array<boost::shared_ptr<DiagBlockStruture>, 4>>;

/**
 * @brief Get the Schur Nest Mat Array object
 *
 * \code {.cpp}
 *
 * auto nested_data = getSchurNestMatArray(
 *
 *       {schur_dm, block_dm}, shell_data,
 *
 *       {"TENSOR"}, {nullptr}
 *
 *   );
 *
 * auto [nested_mat, nested_data_] =
 *  createSchurNestedMatrix({schur_dm, block_dm}, nested_data);
 * \endcode
 * 
 * @param dms schur dm, and block dm
 * @param block
 * @param fields_name
 * @param field_ents
 * @return SchurNestMatrixData
 */
SchurNestMatrixData
getSchurNestMatArray(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                     SchurShellMatData A,

                     std::vector<std::string> fields_name, //< a00 fields
                     std::vector<boost::shared_ptr<Range>>
                         field_ents //< a00 ranges (can be null)

);

/**
 * @brief Create a Mat Diag Blocks object
 *
 * \code {.cpp}
 *
 * auto [nested_mat, nested_data] = createSchurNestedMatrix(
 *
 *       {schur_dm, block_dm},
 *
 *       getSchurNestMatArray(
 *
 *           {schur_dm, block_dm}, shell_data,
 *
 *           {"TENSOR"}, {nullptr}
 *
 *       )
 *
 *  );
 *
 * \endcode
 *
 * @return Mat
 */
std::pair<SmartPetscObj<Mat>, SchurNestMatrixData>
createSchurNestedMatrix(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                        SchurNestMatrixData schur_net_dat);


/**
 * @brief Set PC for Schur block
 * 
 * @param pc 
 * @return MoFEMErrorCode 
 */
MoFEMErrorCode setSchurMatSolvePC(SmartPetscObj<PC> pc);

struct SchurL2MatsBlock;

/***
 * @brief Specialization of MatSetValues for SchurL2MatsBlock
*/
template <>
MoFEMErrorCode
MatSetValues<SchurL2MatsBlock>(Mat M,
                               const EntitiesFieldData::EntData &row_data,
                               const EntitiesFieldData::EntData &col_data,
                               const MatrixDouble &mat, InsertMode iora);

/***
 * @brief Specialization of VecSetValues for SchurL2MatsBlock
*/
template <>
inline MoFEMErrorCode
VecSetValues<SchurL2MatsBlock>(Vec V, const EntitiesFieldData::EntData &data,
                               const VectorDouble &nf, InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

/***
 * @brief Specialisation of MatSetValues for AssemblyTypeSelector<BLOCK_SCHUR>
*/
template <>
inline MoFEMErrorCode MatSetValues<AssemblyTypeSelector<BLOCK_SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<SchurL2MatsBlock>(M, row_data, col_data, mat, iora);
}

/***
 * @brief Specialisation of VecSetValues for AssemblyTypeSelector<BLOCK_SCHUR>
*/
template <>
inline MoFEMErrorCode VecSetValues<AssemblyTypeSelector<BLOCK_SCHUR>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return VecSetValues<EssentialBcStorage>(V, data, nf, iora);
}

} // namespace MoFEM

#endif //__SCHUR_HPP__