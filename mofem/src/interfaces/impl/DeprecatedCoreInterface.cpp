/** \file DeprecatedCoreInterface.cpp
 * \brief implementation of deprecated functions
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
