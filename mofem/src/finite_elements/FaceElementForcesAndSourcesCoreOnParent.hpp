/** \file FaceElementForcesAndSourcesCoreOnParent.hpp
  \brief Implementation of face element integrating parent

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

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

  int getRule(int order);

protected:
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
};

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__
