/** \file FaceElementForcesAndSourcesCoreOnParent.hpp
  \brief Implementation of face element integrating parent

*/



#ifndef __FACEELEMENTFORCESANDSOURCESCORE_ONPARENT__HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_ONPARENT__HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base face element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct FaceElementForcesAndSourcesCoreOnChildParent
    : public FaceElementForcesAndSourcesCore {

  using FaceElementForcesAndSourcesCore::FaceElementForcesAndSourcesCore;

protected:
  int getRule(int order);
  MoFEMErrorCode setGaussPts(int order);

private:
};

/**
 * @deprecated do not use needed for back compatibility
 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreOnChildParentSwitch
    : public FaceElementForcesAndSourcesCoreOnChildParent {
  using FaceElementForcesAndSourcesCoreOnChildParent::
      FaceElementForcesAndSourcesCoreOnChildParent;
  using UserDataOperator =
      FaceElementForcesAndSourcesCoreOnChildParent::UserDataOperator;
};

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__
