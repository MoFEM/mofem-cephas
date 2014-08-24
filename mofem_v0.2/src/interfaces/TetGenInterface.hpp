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

  typedef map<EntityHandle,unsigned long> moabTetGen_Map;
  typedef map<unsigned long,EntityHandle> tetGenMoab_Map;
  typedef map<int,Range> idxRange_Map;

  PetscErrorCode inData(
    Range& ents,tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  PetscErrorCode outData(
    tetgenio& in,tetgenio& out,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    Range *ents = NULL);

  PetscErrorCode outData(
    tetgenio& in,tetgenio& out,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    BitRefLevel bit);

  PetscErrorCode setFaceData(
    vector<pair<Range,int> >& markers,
    tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  PetscErrorCode getTiangleAttributes(
    tetGenMoab_Map& tetgen_moab_map,tetgenio& out,
    Range *ents = NULL,idxRange_Map *ents_map = NULL);

  PetscErrorCode setReginData(vector<pair<EntityHandle,int> >& regions,tetgenio& in);
  PetscErrorCode getReginData(
    tetGenMoab_Map& tetgen_moab_map,tetgenio& out,
    Range *ents = NULL,idxRange_Map *ents_map = NULL);

  PetscErrorCode tetRahedralize(char switches[],tetgenio& in,tetgenio& out);
  PetscErrorCode loadPoly(char file_name[],tetgenio& in);


  //Tools for TetGen

  PetscErrorCode checkPlanar_Trinagle(double coords[],bool *result,const double eps = 1e-9);
  PetscErrorCode groupPlanar_Triangle(Range &tris,vector<Range> &sorted,const double eps = 1e-9);
  PetscErrorCode groupRegion_Triangle(Range &tris,vector<vector<Range> > &sorted,const double eps = 1e-9);
  PetscErrorCode makePolygonFacet(Range &ents,Range &polygons);

};

}

#endif //__TETGENINTERFACE_HPP__
