/** \file definitions.h
 * \brief useful compiler directives and definitions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __DEFINITONS_H__
#define __DEFINITONS_H__

// taken from http://stackoverflow.com/questions/295120/c-mark-as-deprecated
#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

/** \brief Interfaces IDs
 *
 * To manage different complexities related to field, finite elements mesh
 * refinements, etc. a appropriate interfaces related to each complexities are
 * created. Interfaces by itself could vary by functionality or the same
 * function can me managed with two interfaces with waring level of abstraction.
 *
 */
enum Interfaces {
  UNKNOWNINTERFACE = 1 << 0,

  CORE_INTERFACE = 1 << 0 | 1 << 1,
  DEPRECATED_CORE_INTERFACE = 1 << 0 | 1 << 2,
  PROBLEMSMANAGER_INTERFACE = 1 << 0 | 1 << 3,
  MATRIX_MANAGER_INTERFACE = 1 << 0 | 1 << 4,
  SIMPLE_INTERFACE = 1 << 0 | 1 << 5,
  BASIC_INTERFACE = 1 << 0 | 1 << 6,
  COMM_INTERFACE = 1 << 0 | 1 << 7,
  LOGMANAGER_INTERFACE = 1 << 0 | 1 << 8,

  MESH_REFINE = 1 << 1 | 1 << 2,
  PRISM_INTEFACE = 1 << 1 | 1 << 3,
  SERIES_RECORDER = 1 << 1 | 1 << 4,
  ISMANAGER_INTERFACE = 1 << 1 | 1 << 5,
  VECMANAGER_INTERFACE = 1 << 1 | 1 << 6,
  FIELDBLAS_INTERFACE = 1 << 1 | 1 << 7,
  BITREFMANAGER_INTERFACE = 1 << 1 | 1 << 8,

  TOOLS = 1 << 1 | 1 << 10,

  // Independent Interfaces
  TETGEN_INTERFACE = 1 << 2 | 1 << 3,
  MED_INTERFACE = 1 << 2 | 1 << 4,
  NODEMERGER_INTERFACE = 1 << 2 | 1 << 5,
  BITLEVELCOUPLER_INTERFACE = 1 << 2 | 1 << 6,
  PRISMSFROMSURFACE_INTERFACE = 1 << 2 | 1 << 7,
  MESHSETSMANAGER_INTERFACE = 1 << 2 | 1 << 8,
  COORDSSYSTEMMANAGER_INTERFACE = 1 << 2 | 1 << 9,
  CUTMESH_INTERFACE = 1 << 2 | 1 << 10,
  FIELDEVALUATOR_INTERFACE = 1 << 2 | 1 << 11
};

enum LoopInterfaces {
  // Loop Methods
  PETSC_DATA_METHOD = 1 << 3,
  KSP_METHOD = 1 << 3 | 1 << 4,
  SNES_METHOD = 1 << 3 | 1 << 5,
  TS_METHOD = 1 << 3 | 1 << 6,
  BASIC_METHOD = 1 << 3 | 1 << 4 | 1 << 5 | 1 << 6,
  FE_METHOD = 1 << 3 | 1 << 4 | 1 << 5 | 1 << 6 | 1 << 7,
  ENTITY_METHOD = 1 << 3 | 1 << 4 | 1 << 5 | 1 << 6 | 1 << 9,
  DOF_METHOD = 1 << 3 | 1 << 4 | 1 << 5 | 1 << 6 | 1 << 9,
};

/**
 * \brief interfaces for PETSc DM interfaces
 */
enum DMInterfaces {
  UNKNOWN_DM_INTERFACE = 1 << 4 | 1 << 5,
  DMCTX_INTERFACE = 1 << 4 | 1 << 6
};

/**
 * \brief Interfaces uses to manage base functions
 */
enum BaseIntefaces {
  UNKNOWN_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6,
  LEGENDRE_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 7,
  LOBATTO_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 8,
  KERNEL_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 9,
  JACOBI_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 10,
  INTEGRATED_JACOBI_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 11,
  ENT_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 7,
  TET_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 8,
  TRI_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 9,
  EDGE_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 10,
  FATPRISM_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 11,
  FLATPRISM_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 12,
  QUAD_BASE_FUNCTION_INTERFACE = 1 << 5 | 1 << 6 | 1 << 13
};

/** \brief Error handling
 *
 * This is complementary to PETSC error codes. The numerical values for
 * these are defined in include/petscerror.h. The names are defined in err.c
 *
 * MoAB error messages are defined in moab/Types.hpp
 *
 */
enum MoFEMErrorCodes {
  MOFEM_SUCCESS = 0,
  MOFEM_DATA_INCONSISTENCY = 100,
  MOFEM_NOT_IMPLEMENTED = 101,
  MOFEM_NOT_FOUND = 102,
  MOFEM_OPERATION_UNSUCCESSFUL = 103,
  MOFEM_IMPOSIBLE_CASE = 104,
  MOFEM_INVALID_DATA = 105,
  MOFEM_NOT_INSTALLED = 106,
  MOFEM_MOFEMEXCEPTION_THROW = 107,
  MOFEM_STD_EXCEPTION_THROW = 108,
  MOFEM_ATOM_TEST_INVALID = 109,
  MOFEM_MOAB_ERROR = 110
};

const static char *const MoFEMErrorCodesNames[] = {
    "MOFEM_SUCCESS",
    "MOFEM_DATA_INCONSISTENCY",
    "MOFEM_NOT_IMPLEMENTED",
    "MOFEM_NOT_FOUND",
    "MOFEM_OPERATION_UNSUCCESSFUL",
    "MOFEM_IMPOSIBLE_CASE",
    "MOFEM_INVALID_DATA",
    "MOFEM_MOFEMEXCEPTION_THROW",
    "MOFEM_STD_EXCEPTION_THROW",
    "MOFEM_ATOM_TEST_INVALID",
    "MOFEM_MOAB_ERROR"};

/// \brief approximation base
enum FieldApproximationBase {
  NOBASE = 0,
  AINSWORTH_LEGENDRE_BASE =
      1, ///< Ainsworth Cole (Legendre) approx. base \cite NME:NME847
  AINSWORTH_LOBATTO_BASE, ///< Like AINSWORTH_LEGENDRE_BASE but with Lobatto
                          ///< base instead Legendre \cite beriot2015efficient
  AINSWORTH_BERNSTEIN_BEZIER_BASE, ///< See \cite ainsworth2011bernstein and
                                   ///< \cite ainsworth2018bernstein
  DEMKOWICZ_JACOBI_BASE, ///< Construction of base is by Demkowicz \cite
                         ///< fuentes2015orientation
  USER_BASE,             ///< user implemented approximation base
  LASTBASE
};

const static char *const ApproximationBaseNames[] = {
    "NOBASE",
    "AINSWORTH_LEGENDRE_BASE",
    "AINSWORTH_LOBATTO_BASE",
    "AINSWORTH_BERNSTEIN_BEZIER_BASE",
    "DEMKOWICZ_JACOBI_BASE",
    "USER_BASE",
    "LASTBASE"};

/// \brief approximation spaces
enum FieldSpace {
  NOSPACE = 0,
  NOFIELD = 1, ///< scalar or vector of scalars describe (no true field)
  H1,          ///< continuous field
  HCURL,       ///< field with continuous tangents
  HDIV,        ///< field with continuous normal traction
  L2,          ///< field with C-1 continuity
  LASTSPACE    ///< FieldSpace in [ 0, LASTSPACE )
};

const static char *const FieldSpaceNames[] = {
    "NOSPACE", "NOFIELD", "H1", "HCURL", "HDIV", "L2", "LASTSPACE"};

/// \brief Those types control how functions respond on arguments, f.e. error
/// handling
enum MoFEMTypes { MF_ZERO = 0, MF_EXCL = 1 << 0, MF_EXIST = 1 << 1 };

/// \brief RowColData
enum RowColData { ROW = 0, COL, DATA, LASTROWCOLDATA };

/**
 *  Controls adjency multi_index container (e.g. BYROW is adjacenciecy by field
 * on on rows), see \ref MoFEM::FieldEntityEntFiniteElementAdjacencyMap
 *
 */
enum ByWhat {
  BYROW = 1 << 0,
  BYCOL = 1 << 1,
  BYDATA = 1 << 2,
  BYROWDATA = 1 << 0 | 1 << 2,
  BYCOLDATA = 1 << 1 | 1 << 2,
  BYROWCOL = 1 << 0 | 1 << 1,
  BYALL = 1 << 0 | 1 << 1 | 1 << 2
};

/**
 * \brief Types of sets and boundary conditions
 *
 */
enum CubitBC {
  UNKNOWNSET = 0,
  NODESET = 1 << 0,
  SIDESET = 1 << 1,
  BLOCKSET = 1 << 2,
  MATERIALSET = 1 << 3,
  DISPLACEMENTSET = 1 << 4,
  FORCESET = 1 << 5,
  PRESSURESET = 1 << 6,
  VELOCITYSET = 1 << 7,
  ACCELERATIONSET = 1 << 8,
  TEMPERATURESET = 1 << 9,
  HEATFLUXSET = 1 << 10,
  INTERFACESET = 1 << 11,
  UNKNOWNNAME = 1 << 12,
  MAT_ELASTICSET = 1 << 13, ///< block name is "MAT_ELASTIC"
  MAT_INTERFSET = 1 << 14,
  MAT_THERMALSET = 1 << 15,  ///< block name is "MAT_THERMAL"
  BODYFORCESSET = 1 << 16,   ///< block name is "BODY_FORCES"
  MAT_MOISTURESET = 1 << 17, ///< block name is "MAT_MOISTURE"
  DIRICHLET_BC = 1 << 18,
  NEUMANN_BC = 1 << 19,
  LASTSET_BC = 1 << 20
};

// OFF_DEPRECATED static const unsigned int UNKNOWNCUBITNAME = UNKNOWNNAME;
// OFF_DEPRECATED static const unsigned int LASTCUBITSET = LASTSET_BC;

/**
 * \brief Names of types of sets and boundary conditions
 */
const static char *const CubitBCNames[] = {
    "UNKNOWNSET",     "NODESET",         "SIDESET",         "BLOCKSET",
    "MATERIALSET",    "DISPLACEMENTSET", "FORCESET",        "PRESSURESET",
    "VELOCITYSET",    "ACCELERATIONSET", "TEMPERATURESET",  "HEATFLUXSET",
    "INTERFACESET",   "UNKNOWNNAME",     "MAT_ELASTICSET",  "MAT_INTERFSET",
    "MAT_THERMALSET", "BODYFORCESSET",   "MAT_MOISTURESET", "DIRICHLET_BC",
    "NEUMANN_BC",     "LASTSET_BC"};

/**
 * \brief Format in rows of vectorial base functions
 */
enum HVecFormatting { HVEC0 = 0, HVEC1, HVEC2 };

/**
 * \brief Format in rows of vectorial base gradients of base functions
 */
enum HVecDiffFormatting {
  HVEC0_0 = 0,
  HVEC1_0,
  HVEC2_0,
  HVEC0_1,
  HVEC1_1,
  HVEC2_1,
  HVEC0_2,
  HVEC1_2,
  HVEC2_2
};

/**
 * \brief Verbosity levels
 */
enum VERBOSITY_LEVELS {
  DEFAULT_VERBOSITY = -1,
  QUIET = 0,
  VERBOSE,
  VERY_VERBOSE,
  NOISY,
  VERY_NOISY
};

#define MAX_CORE_TMP 1       ///< maximal number of cores
#define BITREFEDGES_SIZE 6   ///< number of edges on tetrahedral
#define BITREFLEVEL_SIZE 128 ///< max number of refinements
#define BITFIELDID_SIZE 32   ///< max number of fields
#define BITFEID_SIZE 32      ///< max number of finite elements
#define BITPROBLEMID_SIZE 32 ///< max number of problems
#define BITINTERFACEUID_SIZE 32

#define MYPCOMM_INDEX 0 ///< default communicator number PCOMM

// This Is form MOAB
#define MB_TYPE_WIDTH 4
#define MB_ID_WIDTH (8 * sizeof(EntityHandle) - MB_TYPE_WIDTH)
#define MB_TYPE_MASK ((EntityHandle)0xF << MB_ID_WIDTH)
//             2^MB_TYPE_WIDTH-1 ------^

#define MB_START_ID ((EntityID)1) ///< All entity id's currently start at 1
#define MB_END_ID                                                              \
  ((EntityID)MB_ID_MASK) ///< Last id is the complement of the MASK
#define MB_ID_MASK (~MB_TYPE_MASK)

#define MAX_DOFS_ON_ENTITY 512 ///< Maximal number of DOFs on entity
#define MAX_PROCESSORS_NUMBER 1024 ///< Maximal number of processors
#define DOF_UID_MASK                                                           \
  (MAX_DOFS_ON_ENTITY - 1) ///< Mask for DOF number on entity form UId
#define ENTITY_UID_MASK (~DOF_UID_MASK)

#define NOT_USED(x) ((void)(x))

/** \brief set barrier start
 * Run code in sequence, starting from process 0, and ends on last process.
 *
 * It can be only used for testing. Do not use that function as a part of these
 * code.
 *
 */
#define BARRIER_PCOMM_RANK_START(PCMB)                                         \
  {                                                                            \
    for (unsigned int i = 0; i < PCMB->proc_config().proc_rank(); i++)         \
      MPI_Barrier(PCMB->proc_config().proc_comm());                            \
  };

/** \deprecated Do use this macro, instead use BARRIER_PCOMM_RANK_START
 */
#define BARRIER_RANK_START(PCMB)                                               \
  {                                                                            \
    macro_is_deprecated_using_deprecated_function();                           \
    for (unsigned int i = 0; i < PCMB->proc_config().proc_rank(); i++)         \
      MPI_Barrier(PCMB->proc_config().proc_comm());                            \
  };

/** \brief set barrier start
 * Run code in sequence, starting from process 0, and ends on last process.
 *
 * It can be only used for testing. Do not use that function as a part of these
 * code.
 *
 */
#define BARRIER_PCOMM_RANK_END(PCMB)                                           \
  {                                                                            \
    for (unsigned int i = PCMB->proc_config().proc_rank();                     \
         i < PCMB->proc_config().proc_size(); i++)                             \
      MPI_Barrier(PCMB->proc_config().proc_comm());                            \
  };

/** \deprecated Do use this macro, instead use BARRIER_PCOMM_RANK_START
 */
#define BARRIER_RANK_END(PCMB)                                                 \
  {                                                                            \
    macro_is_deprecated_using_deprecated_function();                           \
    for (unsigned int i = PCMB->proc_config().proc_rank();                     \
         i < PCMB->proc_config().proc_size(); i++)                             \
      MPI_Barrier(PCMB->proc_config().proc_comm());                            \
  };

/** \brief set barrier start
 * Run code in sequence, starting from process 0, and ends on last process.
 *
 * It can be only used for testing. Do not use that function as a part of these
 * code.
 *
 */
#define BARRIER_MOFEM_RANK_START(MOFEM)                                        \
  {                                                                            \
    for (int i = 0; i < (MOFEM)->get_comm_rank(); i++)                         \
      MPI_Barrier((MOFEM)->get_comm());                                        \
  };

/** \brief set barrier start
 * Run code in sequence, starting from process 0, and ends on last process.
 *
 * It can be only used for testing. Do not use that function as a part of these
 * code.
 *
 */
#define BARRIER_MOFEM_RANK_END(MOFEM)                                          \
  {                                                                            \
    for (int i = (MOFEM)->get_comm_rank(); i < (MOFEM)->get_comm_size(); i++)  \
      MPI_Barrier((MOFEM)->get_comm());                                        \
  };

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Is used to indicate that macro is deprecated
 * Do nothing just triggers error at the compilation
 */
DEPRECATED void macro_is_deprecated_using_deprecated_function();

#ifdef __cplusplus
}
#endif

/**
 * \brief First executable line of each MoFEM function, used for error handling.
 Final line of MoFEM functions should be MoFEMFunctionReturn(0);

   \node Not collective

   Example
   \code
   PetscErrorCode fun()  {
    int something;
    MoFEMFunctionBegin;
    MoFEMFunctionReturn(0);
   }
   \endcode

 */
#define MoFEMFunctionBegin                                                     \
  PetscFunctionBegin;                                                          \
  try {

/**
 * @brief \brief Catch errors
 *
 * Usage in main functions
 * \code
 * int main(int argc, char *argv[]) {
 *
 * MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
 *
 * try {
 *
 * // More code here
 *
 * }
 * CATCH_ERRORS;
 *
 * return MoFEM::Core::Finalize();
 *
 * }
 * \endcode
 *
 */
#define CATCH_ERRORS                                                           \
  catch (MoFEMExceptionInitial const &ex) {                                    \
    return PetscError(PETSC_COMM_SELF, ex.lINE, PETSC_FUNCTION_NAME, __FILE__, \
                      ex.errorCode, PETSC_ERROR_INITIAL, ex.what());           \
  }                                                                            \
  catch (MoFEMExceptionRepeat const &ex) {                                     \
    return PetscError(PETSC_COMM_SELF, ex.lINE, PETSC_FUNCTION_NAME, __FILE__, \
                      ex.errorCode, PETSC_ERROR_REPEAT, " ");                  \
  }                                                                            \
  catch (MoFEMException const &ex) {                                           \
    SETERRQ(PETSC_COMM_SELF, ex.errorCode, ex.errorMessage);                   \
  }                                                                            \
  catch (boost::bad_weak_ptr & ex) {                                           \
    std::string message("Boost bad weak ptr: " + std::string(ex.what()) +      \
                        " at " + boost::lexical_cast<std::string>(__LINE__) +  \
                        " : " + std::string(__FILE__) + " in " +               \
                        std::string(PETSC_FUNCTION_NAME));                     \
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, message.c_str());      \
  }                                                                            \
  catch (std::out_of_range & ex) {                                             \
    std::string message("Std out of range error: " + std::string(ex.what()) +  \
                        " at " + boost::lexical_cast<std::string>(__LINE__) +  \
                        " : " + std::string(__FILE__) + " in " +               \
                        std::string(PETSC_FUNCTION_NAME));                     \
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, message.c_str());      \
  }                                                                            \
  catch (std::exception const &ex) {                                           \
    std::string message("Std error: " + std::string(ex.what()) + " at " +      \
                        boost::lexical_cast<std::string>(__LINE__) + " : " +   \
                        std::string(__FILE__) + " in " +                       \
                        std::string(PETSC_FUNCTION_NAME));                     \
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, message.c_str());      \
  }

/**
 * \brief Last executable line of each PETSc function used for error handling.
 * Replaces return()
 * @param  a error code
 *
 * \note MoFEMFunctionReturn has to be used with MoFEMFunctionBegin and can be
 * used only at the end of the function. If is need to return function in
 * earlier use MoFEMFunctionReturnHot
 *
 */
#define MoFEMFunctionReturn(a)                                                 \
  }                                                                            \
  CATCH_ERRORS                                                                 \
  PetscFunctionReturn(a)

/**
  * \brief First executable line of each MoFEM function, used for error
  handling. Final line of MoFEM functions should be MoFEMFunctionReturn(0); Use
  of this function allows for lighter profiling by default.

  \node Not collective

  Example:
  \code
  PetscErrorCode fun()  {
  int something;
  MoFEMFunctionBeginHot;

  // some work here

  MoFEMFunctionReturnHot(0);
  }
  \endcode
*/
#define MoFEMFunctionBeginHot PetscFunctionBeginHot

/**
 * \brief Last executable line of each PETSc function used for error handling.
 * Replaces return()
 * @param  a error code
 */
#define MoFEMFunctionReturnHot(a) PetscFunctionReturn(a)

#define CHKERRQ_PETSC(n) CHKERRQ(n)
/**
 * \brief check error code of MoAB function
 * @param  a MoABErrorCode
 */
#define CHKERRQ_MOAB(a)                                                        \
  if (PetscUnlikely(MB_SUCCESS != (a))) {                                      \
    std::string error_str = (unsigned)(a) <= (unsigned)MB_FAILURE              \
                                ? moab::ErrorCodeStr[a]                        \
                                : "INVALID ERROR CODE";                        \
    std::string str("MOAB error (" + boost::lexical_cast<std::string>((a)) +   \
                    ") " + error_str + " at line " +                           \
                    boost::lexical_cast<std::string>(__LINE__) + " : " +       \
                    std::string(__FILE__));                                    \
    SETERRQ(PETSC_COMM_SELF, MOFEM_MOAB_ERROR, str.c_str());                   \
  }

/**
 * \brief Check error code of MoFEM/MOAB/PETSc function
 * @param  a MoFEMErrorCode
 *
 * \code
 * MoFEMErrorCode fun() {
 * MoFEMFunctionBeginHot;
 * rval = fun_moab(); CHKERRG(rval);
 * ierr = fun_petsc(); CHKERRG(ierr);
 * merr = fun_mofem(); CHKERRG(merr);
 * MoFEMFunctionReturnHot(0);
 * \endcode
 *
 * \note Function detect type of errocode using specialized template function
 * getErrorType, i.e. condition is evaluated at compilation time.
 *
 */
#define CHKERRG(n)                                                             \
  if ((boost::is_same<BOOST_TYPEOF((n)),                                       \
                      MoFEMErrorCodeGeneric<PetscErrorCode>>::value)) {        \
    CHKERRQ_PETSC((n));                                                        \
  } else if (boost::is_same<BOOST_TYPEOF((n)),                                 \
                            MoFEMErrorCodeGeneric<moab::ErrorCode>>::value) {  \
    CHKERRQ_MOAB((n));                                                         \
  }

/**
 * @brief Inline error check
 *
 * \code
 *
 * MoFEMErrorCode foo() {
 *   MoFEMFunctionBegin;
 *
 *   // Call other functions
 *   CHKERR fun_moab();
 *   CHKERR fun_petsc();
 *   CHKERR fun_mofem();
 *
 *   // Throw error
 *   SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "Some error message");
 *
 *   MoFEMFunctionReturn(0);
 * }
 *
 * int main(int argc, char *argv[]) {
 *
 * // Initailise MoFEM and Petsc
 * MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
 *
 * try {
 *
 *   moab::Core mb_instance; // MoAB database
 *   moab::Interface &moab = mb_instance;
 *   MoFEM::Core core(moab); // MOFEM database
 *   MoFEM::CoreInterface &m_field = core;
 *
 *   CHKERR foo(); // Call function
 *
 * }
 * CATCH_ERRORS;
 *
 * return MoFEM::Core::Finalize();
 *
 * }
 *
 * \endcode
 *
 */
#define CHKERR MoFEM::ErrorChecker<__LINE__>() <<

/**
 * \brief Check error code of MoAB function and throw MoFEM exception
 * @param  a MoABErrorCode
 */
#define MOAB_THROW(a)                                                          \
  {                                                                            \
    if (PetscUnlikely(MB_SUCCESS != (a))) {                                    \
      std::string error_str = (unsigned)(a) <= (unsigned)MB_FAILURE            \
                                  ? moab::ErrorCodeStr[a]                      \
                                  : "INVALID ERROR CODE";                      \
      std::string str("MOAB error (" + boost::lexical_cast<std::string>((a)) + \
                      ") " + error_str + " at line " +                         \
                      boost::lexical_cast<std::string>(__LINE__) + " : " +     \
                      std::string(__FILE__));                                  \
      throw MoFEMException(MOFEM_MOAB_ERROR, str.c_str());                     \
    }                                                                          \
  }

/**
 * \brief Throw MoFEM exception
 * @param  a message
 */
#define THROW_MESSAGE(a)                                                       \
  {                                                                            \
    throw MoFEMException(                                                      \
        MOFEM_MOFEMEXCEPTION_THROW,                                            \
        ("MoFEM error " + boost::lexical_cast<std::string>((a)) +              \
         " at line " + boost::lexical_cast<std::string>(__LINE__) + " : " +    \
         std::string(__FILE__))                                                \
            .c_str());                                                         \
  }

/**
 * \brief Convert number to string
 * @param  x number
 */
#define SSTR(x) toString(x)

#define TENSOR1_VEC_PTR(VEC) &VEC[0], &VEC[1], &VEC[2]

#define SYMMETRIC_TENSOR4_MAT_PTR(MAT)                                         \
  &MAT(0, 0), &MAT(0, 1), &MAT(0, 2), &MAT(0, 3), &MAT(0, 4), &MAT(0, 5),      \
      &MAT(1, 0), &MAT(1, 1), &MAT(1, 2), &MAT(1, 3), &MAT(1, 4), &MAT(1, 5),  \
      &MAT(2, 0), &MAT(2, 1), &MAT(2, 2), &MAT(2, 3), &MAT(2, 4), &MAT(2, 5),  \
      &MAT(3, 0), &MAT(3, 1), &MAT(3, 2), &MAT(3, 3), &MAT(3, 4), &MAT(3, 5),  \
      &MAT(4, 0), &MAT(4, 1), &MAT(4, 2), &MAT(4, 3), &MAT(4, 4), &MAT(4, 5),  \
      &MAT(5, 0), &MAT(5, 1), &MAT(5, 2), &MAT(5, 3), &MAT(5, 4), &MAT(5, 5)

#define TENSOR4_MAT_PTR(MAT) &MAT(0, 0), MAT.size2()

#define TENSOR2_MAT_PTR(MAT)                                                   \
  &MAT(0, 0), &MAT(1, 0), &MAT(2, 0), &MAT(3, 0), &MAT(4, 0), &MAT(5, 0),      \
      &MAT(6, 0), &MAT(7, 0), &MAT(8, 0)

#define SYMMETRIC_TENSOR2_MAT_PTR(MAT)                                         \
  &MAT(0, 0), &MAT(0, 1), &MAT(0, 2), &MAT(0, 3), &MAT(0, 4), &MAT(0, 5)

#define SYMMETRIC_TENSOR2_VEC_PTR(VEC)                                         \
  &VEC[0], &VEC[1], &VEC[2], &VEC[3], &VEC[4], &VEC[5]

#endif //__DEFINITONS_H__
