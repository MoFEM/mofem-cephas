/** \file VertexElementForcesAndSourcesCore.cpp

\brief Implementation of vertex element

*/



#ifdef __cplusplus
extern "C" {
#endif
#include <cblas.h>
#include <lapack_wrap.h>
#include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

VertexElementForcesAndSourcesCore::VertexElementForcesAndSourcesCore(
    Interface &m_field)
    : ForcesAndSourcesCore(m_field) {
  CHK_THROW_MESSAGE(createDataOnElement(MBVERTEX),
                 "Problem with creation data on element");
};

MoFEMErrorCode VertexElementForcesAndSourcesCore::operator()() {
  MoFEMFunctionBegin;

  if (numeredEntFiniteElementPtr->getEntType() != MBVERTEX)
    MoFEMFunctionReturnHot(0);


  EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
  coords.resize(3, false);
  CHKERR mField.get_moab().get_coords(&ent, 1, &*coords.data().begin());

  // Iterate over operators
  CHKERR loopOverOperators();

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode VertexElementForcesAndSourcesCore::UserDataOperator::setPtrFE(
    ForcesAndSourcesCore *ptr) {
  MoFEMFunctionBeginHot;
  if (!(ptrFE = dynamic_cast<VertexElementForcesAndSourcesCore *>(ptr)))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "User operator and finite element do not work together");
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
