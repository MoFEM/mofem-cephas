/** \file DMCtxImpl.hpp
  \brief Implementation of DM context. You should not use it directly
  */

#ifndef __DMCTX_IMPL_H
#define __DMCTX_IMPL_H

namespace MoFEM {

struct DMCtxImpl : public DMCtx {

  DMCtxImpl();

  int useCount() { return referenceNumber; }
  int incrementReference() { return ++referenceNumber; }

  int rAnk = -1; //< processor rank
  int sIze = -1; //< communication size

  int verbosity = VERBOSE; ///< verbosity
  int referenceNumber = 0; //< reference number

  // sub problem
  PetscBool isSubDM = PETSC_FALSE;
  std::vector<std::string> rowSubFields;
  std::vector<std::string> colSubFields;
  const Problem *problemMainOfSubPtr; ///< pointer to main problem to sub-problem

  PetscBool isCompDM = PETSC_FALSE;
  std::vector<std::string> rowCompPrb;
  std::vector<std::string> colCompPrb;
  std::map<std::string, boost::shared_ptr<Range>> mapTypeRow;
  std::map<std::string, boost::shared_ptr<Range>> mapTypeCol;

  // Options
  PetscBool isPartitioned = PETSC_FALSE;  ///< true if read mesh is on parts
  PetscBool isSquareMatrix = PETSC_TRUE;  ///< true if rows equals to cols
  PetscBool destroyProblem = PETSC_FALSE; ///< If true destroy problem with DM
  PetscBool isProblemBuild = PETSC_FALSE; ///< True if problem is build

  Interface *mField_ptr = nullptr; ///< MoFEM interface

  // pointer to data structures
  const Problem *problemPtr = nullptr; ///< pointer to problem data structure
  std::string problemName;             ///< Problem name

  // schur block matrix
  boost::shared_ptr<BlockStruture> blocMatDataPtr;
  boost::shared_ptr<NestSchurData> nestedSchurDataPtr;
};

} // namespace MoFEM

#endif //__DMCTX_IMPL_H