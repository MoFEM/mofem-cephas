/** \file DeprecatedCoreInterface.cpp
 * \brief implementation of deprecated functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#include <MoFEM.hpp>

namespace MoFEM {

MoFEMErrorCode
DeprecatedCoreInterface::seed_ref_level_2D(const EntityHandle meshset,
                                           const BitRefLevel &bit, int verb) {
  return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset, 2, bit,
                                                            verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::seed_ref_level_3D(const EntityHandle meshset,
                                           const BitRefLevel &bit, int verb) {
  return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset, 3, bit,
                                                            verb);
}

MoFEMErrorCode DeprecatedCoreInterface::seed_ref_level(const Range &ents,
                                                       const BitRefLevel &bit,
                                                       const bool only_tets,
                                                       int verb) {
  return getInterface<BitRefManager>()->setBitRefLevel(ents, bit, only_tets,
                                                       verb);
}

MoFEMErrorCode DeprecatedCoreInterface::partition_mesh(const Range &ents,
                                                       const int dim,
                                                       const int adj_dim,
                                                       const int n_parts,
                                                       int verb) {
  return getInterface<ProblemsManager>()->partitionMesh(
      ents, dim, adj_dim, n_parts, NULL, NULL, NULL, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::synchronise_entities(Range &ent,
                                                             int verb) {
  return getInterface<CommInterface>()->synchroniseEntities(ent, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::synchronise_field_entities(const std::string &name,
                                                    int verb) {
  return getInterface<CommInterface>()->synchroniseFieldEntities(name, verb);
}

} // namespace MoFEM
