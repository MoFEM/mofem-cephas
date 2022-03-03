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
struct FaceElementForcesAndSourcesCoreOnChildParentBase
    : public FaceElementForcesAndSourcesCoreBase {

  using FaceElementForcesAndSourcesCoreBase::
      FaceElementForcesAndSourcesCoreBase;

  int getRule(int order);

protected:
  MoFEMErrorCode setGaussPts(int order);

private:
};

/**
 * @brief Face side finite element with switches
 *
 * Using SWITCH to off functions
 *
 * @tparam SWITCH
 */
template <int SWITCH>
struct FaceElementForcesAndSourcesCoreOnChildParentSwitch
    : public FaceElementForcesAndSourcesCoreOnChildParentBase {

  // using FaceElementForcesAndSourcesCoreOnChildParentBase::
  //     FaceElementForcesAndSourcesCoreOnChildParentBase;

  FaceElementForcesAndSourcesCoreOnChildParentSwitch(MoFEM::Interface &m_field)
      : FaceElementForcesAndSourcesCoreOnChildParentBase(m_field) {}


  MoFEMErrorCode operator()();
};

template <int SWITCH>
MoFEMErrorCode
FaceElementForcesAndSourcesCoreOnChildParentSwitch<SWITCH>::operator()() {
  return opSwitch<SWITCH>();
}

} // namespace MoFEM

#endif //__FACEELEMENTFORCESANDSOURCESCORE_ONSIDE___HPP__
