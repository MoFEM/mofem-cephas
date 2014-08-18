/** \file TetGebInterface.hpp
 * \brief TetGen interface 
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __TETGENINTERFACE_HPP__
#define __TETGENINTERFACE_HPP__

#include "FieldUnknownInterface.hpp"

struct tetgenio;

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMTetGegInterface = MOFEMuuid( BitIntefaceId(TETGEN_INTERFACE) );

struct TetGenInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  TetGenInterface(MoFEM::Core& core): cOre(core) {};

  PetscErrorCode inData(
    Range& ents,tetgenio& in,
    map<EntityHandle,unsigned long>& moab_tetgen_map,
    map<unsigned long,EntityHandle>& tetgen_moab_map);

  PetscErrorCode outData(
    Range& ents,tetgenio& in,tetgenio& out,
    map<EntityHandle,unsigned long>& moab_tetgen_map,
    map<unsigned long,EntityHandle>& tetgen_moab_map);

  PetscErrorCode setFacetMarkers(
    Raneg ents,int marker,tetgenio& in);

};

}

#endif //__TETGENINTERFACE_HPP__
