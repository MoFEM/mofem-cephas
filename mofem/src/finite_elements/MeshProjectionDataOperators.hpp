/** \file MeshProjectionDataOperators.hpp
  * \brief Mesh projection operators

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

#ifndef __MESH_PROJECTION_DATA_OPERATORS_HPP__
#define __MESH_PROJECTION_DATA_OPERATORS_HPP__

namespace MoFEM {

struct OpRunParent : public ForcesAndSurcesCore::UserDataOperator {
  OpRunParent(boost::shared_ptr<ForcesAndSourcesCore> parent_ele_ptr,
              BitRefLevel bit_parent, BitRefLevel bit_parent_mask,
              boost::shared_ptr<ForcesAndSourcesCore> this_ele_ptr,
              BitRefLevel bit_this, BitRefLevel bit_this_mask, int verb = QUIET,
              Sev sev = Sev::verbose);

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data);

private:
  boost::shared_ptr<ForcesAndSourcesCore> parentElePtr;
  boost::shared_ptr<ForcesAndSourcesCore> thisElePtr;
  const BitRefLevel bitParent;
  const BitRefLevel bitParentMask;
  const BitRefLevel bitThis;
  const BitRefLevel bitThisMask;
  int verbosity;
  Sev severityLevel;
};
}

#endif //__MESH_PROJECTION_DATA_OPERATORS_HPP__