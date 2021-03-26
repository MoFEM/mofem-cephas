/** \file CoordSysMultiIndices.cpp
 * \brief Tensor coordinate system
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

namespace MoFEM {

CoordSys::CoordSys(const moab::Interface &moab, const EntityHandle meshset)
    : meshSet(meshset), tagCoordSysName(NULL) {
  // Change those tags only by modifiers
  // dim
  Tag th_coord_sys_dim;
  rval = moab.tag_get_handle("_CoordSysDim", th_coord_sys_dim);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_coord_sys_dim, &meshSet, 1,
                             (const void **)&tagCoordSysDim);
  MOAB_THROW(rval);
  // Coord Sys Name
  Tag th_coord_sys_name;
  rval = moab.tag_get_handle("_CoordSysName", th_coord_sys_name);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_coord_sys_name, &meshset, 1,
                             (const void **)&tagCoordSysName,
                             &tagCoordSysNameSize);
  MOAB_THROW(rval);
}

} // namespace MoFEM
