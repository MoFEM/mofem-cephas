/** \file FaceElementForcesAndSourcesCoreOnParent.hpp
  \brief Implementation of face element integrating parent

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
