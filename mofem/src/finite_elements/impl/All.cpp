#include <MoFEM.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  // #include <gm_rule.h>
  #include <quad.h>
#ifdef __cplusplus
}
#endif

#include "impl/EntitiesFieldData.cpp"
#include "impl/DataOperators.cpp"
#include "impl/ForcesAndSourcesCore.cpp"
#include "impl/UserDataOperators.cpp"
#include "impl/HODataOperators.cpp"
#include "impl/MeshProjectionDataOperators.cpp"
#include "impl/BaseDerivativesDataOperators.cpp"
#include "impl/FormsIntegrators.cpp"
#include "impl/VolumeElementForcesAndSourcesCore.cpp"
#include "impl/FaceElementForcesAndSourcesCore.cpp"
#include "impl/EdgeElementForcesAndSourcesCore.cpp"
#include "impl/EdgeElementForcesAndSourcesCoreOnParent.cpp"
#include "impl/VertexElementForcesAndSourcesCore.cpp"
#include "impl/FlatPrismElementForcesAndSourcesCore.cpp"
#include "impl/ContactPrismElementForcesAndSourcesCore.cpp"
#include "impl/FatPrismElementForcesAndSourcesCore.cpp"
#include "impl/VolumeElementForcesAndSourcesCoreOnSide.cpp"
#include "impl/FaceElementForcesAndSourcesCoreOnSide.cpp"
#include "impl/FaceElementForcesAndSourcesCoreOnParent.cpp"
#include "impl/VolumeElementForcesAndSourcesCoreOnContactPrismSide.cpp"
#include "impl/Schur.cpp"