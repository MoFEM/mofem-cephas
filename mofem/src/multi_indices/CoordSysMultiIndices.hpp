/** \file CordSysMultiIndices.hpp
 * \ingroup coordsys_multi_indices
 * \brief Coordinate systems attached to DOFs
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __COORDSYSMULTIINDICES_HPP__
#define __COORDSYSMULTIINDICES_HPP__

namespace MoFEM {

  struct CoordSys {
    EntityHandle meshSet; 		      ///< keeps entities for this meshset
    const int* tagIdData;
    const char* tagCoordSysName; 		///< tag keeps name of the field
    int tagCoordSysNameSize;
    CoordSys(Interface &moab,const EntityHandle meshset);
    inline int getId() const { return *tagIdData; };
    inline EntityHandle getMeshSet() const { return meshSet; };
    inline boost::string_ref getNameRef() const { return boost::string_ref((char *)tagCoordSysName,tagCoordSysNameSize); };
    inline string getName() const { return string((char *)tagCoordSysName,tagCoordSysNameSize); };
  };

  typedef multi_index_container<
    CoordSys,
    indexed_by<
      ordered_non_unique<
        tag<CoordSysID_mi_tag>, const_mem_fun<CoordSys,int,&CoordSys::getId> >,
      ordered_unique<
        tag<Meshset_mi_tag>, member<CoordSys,EntityHandle,&CoordSys::meshSet>
      >,
      ordered_unique<
        tag<CoordSysName_mi_tag >, const_mem_fun<CoordSys,boost::string_ref,&CoordSys::getNameRef>
      >
  > > CoordSys_multiIndex;

}

#endif //__COORDSYSMULTIINDICES_HPP__

/***************************************************************************//**
 * \defgroup coordsys_multi_indices MoFEM coordinate system for multi indices
 * \ingroup mofem
 ******************************************************************************/
