/** \file CoordSysMultiIndices.cpp
 * \brief MoFEM Coordinate system
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>

namespace MoFEM {

  CoordSys::CoordSys(Interface &moab,const EntityHandle meshset):
  meshSet(meshset),
  tagIdData(NULL),
  tagCoordSysName(NULL) {
    // Change those tags only by modifiers
    ErrorCode rval;
    // Id
    Tag th_coord_sys_id;
    rval = moab.tag_get_handle("_CoordSysId",th_coord_sys_id); CHKERR_THROW(rval);
    rval = moab.tag_get_by_ptr(th_coord_sys_id,&meshSet,1,(const void **)&tagIdData); CHKERR_THROW(rval);
    // dim
    Tag th_coord_sys_dim;
    rval = moab.tag_get_handle("_CoordSysDim",th_coord_sys_dim); CHKERR_THROW(rval);
    rval = moab.tag_get_by_ptr(th_coord_sys_dim,&meshSet,1,(const void **)&tagCoordSysDim); CHKERR_THROW(rval);
    // Coord Sys Name
    Tag th_coord_sys_name;
    rval = moab.tag_get_handle("_CoordSysName",th_coord_sys_name); CHKERR(rval);
    rval = moab.tag_get_by_ptr(th_coord_sys_name,&meshset,1,(const void **)&tagCoordSysName,&tagCoordSysNameSize); CHKERR_THROW(rval);

  }

}
