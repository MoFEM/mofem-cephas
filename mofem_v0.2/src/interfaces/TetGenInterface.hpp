/** \file TetGenInterface.hpp
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

/** \brief use TetGen to generate mesh
  * \ingroup mofem
  */
struct TetGenInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  TetGenInterface(MoFEM::Core& core): cOre(core) {};

  typedef map<EntityHandle,unsigned long> moabTetGen_Map;
  typedef map<unsigned long,EntityHandle> tetGenMoab_Map;
  typedef map<int,Range> idxRange_Map;

  
  /** \brief create tetgen data structure form range of moab entities
    
    \param ents range of entitities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map maping moab to tegen enties
    \param tetgen_moab_map maping tegen to moab entities

    */
  PetscErrorCode inData(
    Range& ents,tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  enum tetGenNodesTypes { RIDGEVERTEX = 0, FREESEGVERTEX = 1, FREEFACETVERTEX = 2, FREEVOLVERTEX = 3 };

  
  /** \brief set point tags and type
    */
  PetscErrorCode setGeomData(
    tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    map<int,Range> &type_ents);

  /** \brief get entities for TetGen data structure

    \param ents range of entitities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map maping moab to tegen enties
    \param tetgen_moab_map maping tegen to moab entities
    \param ents rerun entities which are in TetGen dara strucure
    \param id_in_tags use tags as entity handles, if that is a case use tag to find moab vertex id
    \param error_if_created throw error if node need to be created

    */
  PetscErrorCode outData(
    tetgenio& in,tetgenio& out,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    Range *ents = NULL,
    bool id_in_tags = false,
    bool error_if_created = false);

  /** \brief get entities for TetGen data structure

    \param ents range of entitities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map maping moab to tegen enties
    \param tetgen_moab_map maping tegen to moab entities
    \param ents rerun entities which are in TetGen dara strucure
    \param bit set level to created entities 
    \param error_if_created throw error if node need to be created

    */
  PetscErrorCode outData(
    tetgenio& in,tetgenio& out,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    BitRefLevel bit,
    bool id_in_tags = false,
    bool error_if_created = false);
  
  /** \brief set markers to faces

    \param markers data structure with markers
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map maping moab to tegen enties
    \param tetgen_moab_map maping tegen to moab entities

    */
  PetscErrorCode setFaceData(
    vector<pair<Range,int> >& markers,
    tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  /** \brief get markers to faces

    \param markers data structure with markers
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map maping moab to tegen enties
    \param tetgen_moab_map maping tegen to moab entities

    */
  PetscErrorCode getTriangleMarkers(
    tetGenMoab_Map& tetgen_moab_map,tetgenio& out,
    Range *ents = NULL,idxRange_Map *ents_map = NULL,bool only_non_zero = true);


  /** \brief set region data to tetrahedral
    */
  PetscErrorCode setReginData(vector<pair<EntityHandle,int> >& regions,tetgenio& in);


  /** \brief get region data to tetrahedral
    */
  PetscErrorCode getReginData(
    tetGenMoab_Map& tetgen_moab_map,tetgenio& out,
    Range *ents = NULL,idxRange_Map *ents_map = NULL);

  /** \brief run tetgen
    */
  PetscErrorCode tetRahedralize(char switches[],tetgenio& in,tetgenio& out);

  /** \brief load poly file
    */
  PetscErrorCode loadPoly(char file_name[],tetgenio& in);

  //Tools for TetGen, i.e. geometry reconstruction from mesh

  PetscErrorCode checkPlanar_Trinagle(double coords[],bool *result,const double eps = 1e-9);
  PetscErrorCode groupPlanar_Triangle(Range &tris,vector<Range> &sorted,const double eps = 1e-9);
  PetscErrorCode groupRegion_Triangle(Range &tris,vector<vector<Range> > &sorted,const double eps = 1e-9);

  //FIXME: assumes that are no holes
  PetscErrorCode makePolygonFacet(Range &ents,Range &polygons,
    bool reduce_edges = false,Range *not_reducable_nodes = NULL,const double eps = 1e-9); 

};

}

#endif //__TETGENINTERFACE_HPP__
