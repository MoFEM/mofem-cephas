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

#include "impl/DataStructures.cpp"
#include "impl/DataOperators.cpp"
#include "impl/ForcesAndSourcesCore.cpp"
#include "impl/UserDataOperators.cpp"
#include "impl/VolumeElementForcesAndSourcesCore.cpp"
#include "impl/FaceElementForcesAndSourcesCore.cpp"
#include "impl/EdgeElementForcesAndSourcesCore.cpp"
#include "impl/VertexElementForcesAndSourcesCore.cpp"
#include "impl/FlatPrismElementForcesAndSourcesCore.cpp"
#include "impl/ContactPrismElementForcesAndSourcesCore.cpp"
#include "impl/FatPrismElementForcesAndSourcesCore.cpp"
#include "impl/VolumeElementForcesAndSourcesCoreOnSide.cpp"
#include "impl/FaceElementForcesAndSourcesCoreOnSide.cpp"
#include "impl/VolumeElementForcesAndSourcesCoreOnVolumeSide.cpp"