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
 * \note Try separate floating points operations from book keeping. Also, align
 * memory that blocks follow floating point operations.
 *
 */

#ifndef __SCHUR_HPP__
#define __SCHUR_HPP__

namespace MoFEM {

constexpr const char MoFEM_BLOCK_MAT[] = "mofem_block_mat";

/**
 * @brief Structure to register events for Schur block assembly and solver
 */
struct SchurEvents {
  static PetscLogEvent MOFEM_EVENT_schurMatSetValues;
  static PetscLogEvent MOFEM_EVENT_opSchurAssembleEnd;
  static PetscLogEvent MOFEM_EVENT_BlockStructureSetValues;
  static PetscLogEvent MOFEM_EVENT_BlockStructureMult;
  static PetscLogEvent MOFEM_EVENT_BlockStructureSolve;
  static PetscLogEvent MOFEM_EVENT_zeroRowsAndCols;
  SchurEvents();
};

struct SchurElemMats;
struct SchurElemMatsBlock;
struct SchurElemMatsPreconditionedBlock;

struct OpSchurAssembleBase : public ForcesAndSourcesCore::UserDataOperator {

  OpSchurAssembleBase() = delete;

private:
  using UserDataOperator::UserDataOperator;
};

OpSchurAssembleBase *createOpSchurAssembleBegin();

/**
 * @brief Construct a new Op Schur Assemble End object
 *
 * @param fields_name list of fields (can be empty)
 * @param field_ents list of entities on which schur complement is applied
 * @param schur_aos maps dofs indices from main problem to schur complement
 * @param schur_mat schur matrix
 * @param sym_schur true if schur (matrix) complement is symmetric
 * @param symm_op true if block diagonal is symmetric
 */
OpSchurAssembleBase *createOpSchurAssembleEnd(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    SmartPetscObj<AO> ao = SmartPetscObj<AO>(),
    SmartPetscObj<Mat> schur = SmartPetscObj<Mat>(), bool sym_schur = false,
    bool symm_op = false,
    boost::shared_ptr<BlockStructure> diag_blocks = nullptr);

using SchurFieldPair = std::pair<std::string, std::string>;

using SchurFEOpsFEandFields = std::vector<

    std::pair<std::string, std::vector<SchurFieldPair>>

    >;

/**
 * @brief Create a Mat Diag Blocks object
 *
 * @return Mat
 */
boost::shared_ptr<BlockStructure> createBlockMatStructure(

    DM dm,                                //< dm
    SchurFEOpsFEandFields schur_fe_op_vec //< block elements

);

using SchurShellMatData =
    std::pair<SmartPetscObj<Mat>, boost::shared_ptr<BlockStructure>>;

/**
 * @brief Create a Schur Mat object
 *
 * @param dm
 * @param data
 * @return std::pair<SmartPetscObj<Mat>, boost::shared_ptr<BlockStructure>>
 */
SchurShellMatData createBlockMat(DM dm, boost::shared_ptr<BlockStructure> data);

using NestSchurData = std::tuple<

    std::array<SmartPetscObj<Mat>, 4>,
    std::array<boost::shared_ptr<BlockStructure>, 4>,
    boost::shared_ptr<BlockStructure>,
    std::pair<SmartPetscObj<IS>, SmartPetscObj<IS>>

    >;

/**
 * @brief Get the Schur Nest Mat Array object
 *
 * @param dms schur dm, and block dm
 * @param block mat A data
 * @param fields_name list of fields
 * @param field_ents list of entities on which schur complement is applied
 * @param add_preconditioner_block add block for preconditioner
 * @return boost::shared_ptr<NestSchurData>
 */
boost::shared_ptr<NestSchurData>
getNestSchurData(std::pair<SmartPetscObj<DM>, SmartPetscObj<DM>> dms,
                 boost::shared_ptr<BlockStructure> block_mat_data,

                 std::vector<std::string> fields_name, //< a00 fields
                 std::vector<boost::shared_ptr<Range>> field_ents, //< a00 ents
                 bool add_preconditioner_block = false

);

/**
 * @brief Switch preconditioner
 * 
 * @param block_mat_data 
 * @return MoFEMErrorCode 
 */
MoFEMErrorCode
schurSwitchPreconditioner(boost::shared_ptr<BlockStructure> block_mat_data);

/**
 * @brief Save block matrix as a mesh
 * 
 * @param block_mat_data 
 * @param filename 
 * @return MoFEMErrorCode 
 */
MoFEMErrorCode
schurSaveBlockMesh(boost::shared_ptr<BlockStructure> block_mat_data,
                   std::string filename);

/**
 * @brief Create a Mat Diag Blocks object
 *
 * \code {.cpp}
 *
 * auto [nested_mat, nested_data_ptr] = createSchurNestedMatrix(
 *
 *       getNestSchurData(
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
std::pair<SmartPetscObj<Mat>, boost::shared_ptr<NestSchurData>>
createSchurNestedMatrix(boost::shared_ptr<NestSchurData> schur_net_data_ptr);

template <>
MoFEMErrorCode DMMoFEMSetNestSchurData(DM dm, boost::shared_ptr<NestSchurData>);

/**
 * @brief Set PC for A00 block
 *
 * @param pc
 * @return MoFEMErrorCode
 */
MoFEMErrorCode setSchurA00MatSolvePC(SmartPetscObj<PC> pc);

DEPRECATED MoFEMErrorCode setSchurMatSolvePC(SmartPetscObj<PC> pc);

struct SchurBackendMatSetValuesPtr {
  SchurBackendMatSetValuesPtr() = delete;
  using MatSetValuesPtr = boost::function<MoFEMErrorCode(
      Mat M, const EntitiesFieldData::EntData &row_data,
      const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
      InsertMode iora)>;
  static MatSetValuesPtr matSetValuesPtr; ///< backend assembly function
  static MatSetValuesPtr
      matSetValuesBlockPtr; ///< backend assembly block mat function
  static MatSetValuesPtr
      matSetValuesPreconditionedBlockPtr; ///< backend assembly block
                                          ///< preconditioner mat function
};

template <>
MoFEMErrorCode
MatSetValues<SchurElemMats>(Mat M, const EntitiesFieldData::EntData &row_data,
                            const EntitiesFieldData::EntData &col_data,
                            const MatrixDouble &mat, InsertMode iora);

template <>
inline MoFEMErrorCode
VecSetValues<SchurElemMats>(Vec V, const EntitiesFieldData::EntData &data,
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
  return VecSetValues<SchurElemMats>(V, data, nf, iora);
}

/***
 * @brief Specialization of MatSetValues for BlockStructure
 */
template <>
MoFEMErrorCode
MatSetValues<BlockStructure>(Mat M, const EntitiesFieldData::EntData &row_data,
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
  return MatSetValues<BlockStructure>(M, row_data, col_data, mat, iora);
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

/***
 * @brief Specialization of MatSetValues for SchurElemMatsBlock
 */
template <>
MoFEMErrorCode
MatSetValues<SchurElemMatsBlock>(Mat M,
                                 const EntitiesFieldData::EntData &row_data,
                                 const EntitiesFieldData::EntData &col_data,
                                 const MatrixDouble &mat, InsertMode iora);

/***
 * @brief Specialization of VecSetValues for SchurElemMatsBlock
 */
template <>
inline MoFEMErrorCode
VecSetValues<SchurElemMatsBlock>(Vec V, const EntitiesFieldData::EntData &data,
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
  return MatSetValues<SchurElemMatsBlock>(M, row_data, col_data, mat, iora);
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

/***
 * @brief Specialization of MatSetValues for SchurElemMatsPreconditionedBlock
 */
template <>
MoFEMErrorCode MatSetValues<SchurElemMatsPreconditionedBlock>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora);

/***
 * @brief Specialisation of MatSetValues for
 * AssemblyTypeSelector<BLOCK_PRECONDITIONER_SCHUR>
 */
template <>
inline MoFEMErrorCode
MatSetValues<AssemblyTypeSelector<BLOCK_PRECONDITIONER_SCHUR>>(
    Mat M, const EntitiesFieldData::EntData &row_data,
    const EntitiesFieldData::EntData &col_data, const MatrixDouble &mat,
    InsertMode iora) {
  return MatSetValues<SchurElemMatsPreconditionedBlock>(M, row_data, col_data,
                                                        mat, iora);
}

/***
 * @brief Specialisation of VecSetValues for AssemblyTypeSelector<BLOCK_SCHUR>
 */
template <>
inline MoFEMErrorCode
VecSetValues<AssemblyTypeSelector<BLOCK_PRECONDITIONER_SCHUR>>(
    Vec V, const EntitiesFieldData::EntData &data, const VectorDouble &nf,
    InsertMode iora) {
  return ::VecSetValues(V, data.getIndices().size(),
                        &*data.getIndices().begin(), &*nf.begin(), iora);
}

/**
 * @deprecated do not use
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
DEPRECATED OpSchurAssembleBase *createOpSchurAssembleEnd(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, bool symm_op,
    boost::shared_ptr<BlockStructure> diag_blocks = nullptr);

/**
 * @deprecated do not use
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
DEPRECATED OpSchurAssembleBase *createOpSchurAssembleEnd(
    std::vector<std::string> fields_name,
    std::vector<boost::shared_ptr<Range>> field_ents,
    std::vector<SmartPetscObj<AO>> sequence_of_aos,
    std::vector<SmartPetscObj<Mat>> sequence_of_mats,
    std::vector<bool> sym_schur, std::vector<double> diag_eps, bool symm_op,
    boost::shared_ptr<BlockStructure> diag_blocks = nullptr);

} // namespace MoFEM

#endif //__SCHUR_HPP__