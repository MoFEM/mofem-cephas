/** \file TetGenInterface.hpp
 * \brief TetGen interface
 *
 * MoFEM TetGen interface
 *
 */

/*
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

class tetgenio;

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


  /** \brief create TetGen data structure form range of moab entities

    \param ents range of entities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping moab to TetGen entities
    \param tetgen_moab_map mapping tegen to moab entities

    */
  PetscErrorCode inData(
    Range& ents,tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  enum tetGenNodesTypes {
    RIDGEVERTEX = 0,
    FREESEGVERTEX = 1,
    FREEFACETVERTEX = 2,
    FREEVOLVERTEX = 3 };

  /** \brief set point tags and type

  Set type of entity, look in TetGen manual for details

  \code
  map<int,Range> types_ents;
  //RIDGEVERTEX
  types_ents[TetGenInterface::RIDGEVERTEX].merge(region_tets_skin_without_boundary_nodes);
  //FREESEGVERTEX
  types_ents[TetGenInterface::FREESEGVERTEX].merge(crack_surface_tris_skin_nodes);
  //FREEFACETVERTEX
  types_ents[TetGenInterface::FREEFACETVERTEX].merge(region_tets_skin_nodes);
  types_ents[TetGenInterface::FREEFACETVERTEX] =  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEFACETVERTEX] =  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  //FREEVOLVERTEX
  types_ents[TetGenInterface::FREEVOLVERTEX].merge(region_nodes);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREEFACETVERTEX]);
  \endcode

  */
  PetscErrorCode setGeomData(
    tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map,
    map<int,Range> &type_ents);

  /** \brief get entities for TetGen data structure

    \param ents range of entities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to moab entities
    \param ents rerun entities which are in TetGen Dara structure
    \param id_in_tags use tags as entity handles, if that is a case use tag to find MoAB vertex id
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

    \param ents range of entities (tetrahedrons or nodes)
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities
    \param ents rerun entities which are in TetGen data structure
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
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities

    */
  PetscErrorCode setFaceData(
    vector<pair<Range,int> >& markers,
    tetgenio& in,
    moabTetGen_Map& moab_tetgen_map,
    tetGenMoab_Map& tetgen_moab_map);

  /** \brief get markers to faces

    \param markers data structure with markers
    \param in tegen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities

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

  /** \brief Group surface strangles in planar regions

    \param tris input triangles
    \param sorted output sorted planar faces
    \param eps tolerance

  */
  PetscErrorCode groupRegion_Triangle(Range &tris,vector<vector<Range> > &sorted,const double eps = 1e-9);

  /** make planar polygon facet

    \param ents surface triangles
    \param plygons output list of polygons
    \param reduce_edges reduce edges if on the line
    \param not_reducable_nodes do not reduce node on edge if in this range
    \param eps tolerance

    \bug assumes that are no holes

    */
  PetscErrorCode makePolygonFacet(Range &ents,Range &polygons,
    bool reduce_edges = false,Range *not_reducable_nodes = NULL,const double eps = 1e-9);

};

}

#endif //__TETGENINTERFACE_HPP__
