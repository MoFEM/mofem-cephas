/** \file EdgeElementForcesAndSourcesCoreOnParent.hpp
  \brief Implementation of edge element integrating parent

*/



#ifndef __EDGEELEMENTFORCESANDSOURCESCORE_ONPARENT__HPP__
#define __EDGEELEMENTFORCESANDSOURCESCORE_ONPARENT__HPP__

using namespace boost::numeric;

namespace MoFEM {

/**
 * \brief Base face element used to integrate on skeleton
 * \ingroup mofem_forces_and_sources_volume_element
 */
struct EdgeElementForcesAndSourcesCoreOnChildParent
    : public EdgeElementForcesAndSourcesCore {

  using EdgeElementForcesAndSourcesCore::EdgeElementForcesAndSourcesCore;

protected:
  int getRule(int order);
  MoFEMErrorCode setGaussPts(int order);

private:
};

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__