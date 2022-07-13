/** \file Exceptions.hpp
 * \brief Exceptions and error handlers
 */



#ifndef __EXCEPTIONS_HPP__
#define __EXCEPTIONS_HPP__

namespace MoFEM {

/**
 * @brief Exceptions and handling errors data structures
 *
 */
namespace Exceptions {
/**
 * \brief Exception to catch
 */
struct MoFEMException : public std::exception {
  const int errorCode;
  char errorMessage[1024];
  MoFEMException(const MoFEMErrorCodes error_code)
      : MoFEMException(static_cast<int>(error_code)) {}
  MoFEMException(const MoFEMErrorCodes error_code, const char error_message[])
      : errorCode(error_code) {
    strncpy(errorMessage, error_message, sizeof(errorMessage));
    errorMessage[sizeof(errorMessage) - 1] = '\0';
  }
  const char *what() const throw() { return errorMessage; }

protected:
  MoFEMException(const int error_code) : errorCode(error_code) {
    strcpy(errorMessage, "Houston we have a problem, something is wrong");
  }
};

struct MoFEMExceptionRepeat : public MoFEMException {
  const int lINE;
  MoFEMExceptionRepeat(const int error_code, const int line)
      : MoFEMException(error_code), lINE(line) {
    strcpy(errorMessage, " ");
  }
};

struct MoFEMExceptionInitial : public MoFEMExceptionRepeat {
  MoFEMExceptionInitial(const int error_code, const char error_message[],
                        const int line)
      : MoFEMExceptionRepeat(error_code, line) {
    strncpy(errorMessage, error_message, sizeof(errorMessage));
    errorMessage[sizeof(errorMessage) - 1] = '\0';
  }
};

typedef moab::ErrorCode MoABErrorCode; ///< MoAB error code
typedef PetscErrorCode MoFEMErrorCode; ///< MoFEM/PETSc error code

template <typename TYPE> struct MoFEMErrorCodeGeneric {
  MoFEMErrorCodeGeneric(const TYPE) {}
};

template <> struct MoFEMErrorCodeGeneric<PetscErrorCode> {
  PetscErrorCode iERR;
  MoFEMErrorCodeGeneric(const PetscErrorCode ierr) : iERR(ierr) {}
  inline operator PetscErrorCode() const { return iERR; }
};

template <> struct MoFEMErrorCodeGeneric<moab::ErrorCode> {
  moab::ErrorCode rVAL;
  MoFEMErrorCodeGeneric(const moab::ErrorCode rval) : rVAL(rval) {}
  inline operator moab::ErrorCode() const { return rVAL; }
};

static MoFEMErrorCodeGeneric<moab::ErrorCode> rval =
    MoFEMErrorCodeGeneric<moab::ErrorCode>(MB_SUCCESS);
static MoFEMErrorCodeGeneric<PetscErrorCode> ierr =
    MoFEMErrorCodeGeneric<PetscErrorCode>(0);

/**
 * \brief Error check for inline function check.
 *
 * This class is not used directly, it is called in CHKERR. In case of the error
 * pass line number and that is catch at the end of the function. Information is
 * enriched by function name and file name. Then error is pushed to PETSc error
 * stack.
 *
 * \note This class has no variables and line number is set at compilation.
 * Adding variables to this function will reduce efficiency of the code. Do
 * not do that.
 *
 */
template <int LINE> struct ErrorChecker {

  /**
   * @brief Operator for handling PetscErrorCode and MoFEMErrorCode
   *
   */
  inline void operator<<(const MoFEMErrorCode err) {
    if (PetscUnlikely(err)) {
      throw MoFEMExceptionRepeat(err, LINE);
    }
    return;
  }

  /**
   * @brief Operator for handling moab::ErrorCode
   *
   */
  inline void operator<<(const moab::ErrorCode err) {
    if (PetscLikely(MB_SUCCESS != err)) {
      std::string error_str = (unsigned)err <= (unsigned)MB_FAILURE
                                  ? moab::ErrorCodeStr[err]
                                  : "INVALID ERROR CODE";
      std::string str("MOAB error (" + boost::lexical_cast<std::string>(err) +
                      ") " + error_str);
      throw MoFEMExceptionInitial(MOFEM_MOAB_ERROR, str.c_str(), LINE);
    }
    return;
  }
};

} // namespace Exceptions

using namespace Exceptions;

} // namespace MoFEM

#endif // __EXCEPTIONS_HPP__