/** \file SaveVertexDofOnTag.hpp

Save field DOFS on vertices/tags
This is another example how to use MoFEM::DofMethod when some operator for each
node need to be applied.

*/



#ifndef __SAVEVERTEXDOFONTAG_HPP__
#define __SAVEVERTEXDOFONTAG_HPP__

#include <MoFEM.hpp>
using namespace MoFEM;

namespace BasicFiniteElements {

/** \brief Save field DOFS on vertices/tags
 */
struct SaveVertexDofOnTag : public MoFEM::DofMethod {

  MoFEM::Interface &mField;
  std::string tagName;
  SaveVertexDofOnTag(MoFEM::Interface &m_field, std::string tag_name)
      : mField(m_field), tagName(tag_name) {}

  Tag tH;

  MoFEMErrorCode preProcess() {
    MoFEMFunctionBegin;
    if (!fieldPtr) {
      SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Null pointer, probably field not found");
    }
    if (fieldPtr->getSpace() != H1) {
      SETERRQ(mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
              "Field must be in H1 space");
    }
    std::vector<double> def_vals(fieldPtr->getNbOfCoeffs(), 0);
    rval = mField.get_moab().tag_get_handle(tagName.c_str(), tH);
    if (rval != MB_SUCCESS) {
      CHKERR mField.get_moab().tag_get_handle(
          tagName.c_str(), fieldPtr->getNbOfCoeffs(), MB_TYPE_DOUBLE, tH,
          MB_TAG_CREAT | MB_TAG_SPARSE, &def_vals[0]);
    }

    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode postProcess() {
    MoFEMFunctionBegin;
    MoFEMFunctionReturn(0);
  }

  MoFEMErrorCode operator()() {
    MoFEMFunctionBegin;
    if (dofPtr->getEntType() != MBVERTEX)
      MoFEMFunctionReturnHot(0);
    EntityHandle ent = dofPtr->getEnt();
    int rank = dofPtr->getNbOfCoeffs();
    double tag_val[rank];

    CHKERR mField.get_moab().tag_get_data(tH, &ent, 1, tag_val);
    tag_val[dofPtr->getDofCoeffIdx()] = dofPtr->getFieldData();
    CHKERR mField.get_moab().tag_set_data(tH, &ent, 1, &tag_val);
    MoFEMFunctionReturn(0);
  }
};

} // BasicFiniteElements

#endif // __SAVEVERTEXDOFONTAG_HPP__
