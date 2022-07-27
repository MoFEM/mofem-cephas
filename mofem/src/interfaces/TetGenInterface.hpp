/** \file TetGenInterface.hpp
 * \brief TetGen interface
 *
 * MoFEM TetGen interface
 *
 *
 * \ingroup mesh_tetgen
 */

#ifndef __TETGENINTERFACE_HPP__
#define __TETGENINTERFACE_HPP__

#include "UnknownInterface.hpp"

class tetgenio;

namespace MoFEM {

/** \brief TetGen interface

  * \ingroup mesh_tetgen
  */
struct TetGenInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  TetGenInterface(const MoFEM::Core &core)
      : cOre(const_cast<MoFEM::Core &>(core)) {}

  typedef std::map<EntityHandle, unsigned long> moabTetGen_Map;
  typedef std::map<unsigned long, EntityHandle> tetGenMoab_Map;
  typedef std::map<int, Range> idxRange_Map;

  /** \brief create TetGen data structure form range of moab entities
    *
    * Move mesh to TetGen data structures
    *
    * \param ents range of entities (tetrahedrons or nodes)
    * \param in tetgen data structure (look to TetGen user manual)
    * \param moab_tetgen_map mapping moab to TetGen entities
    * \param tetgen_moab_map mapping tetgen to moab entities
    */
  MoFEMErrorCode inData(Range &ents, tetgenio &in,
                        moabTetGen_Map &moab_tetgen_map,
                        tetGenMoab_Map &tetgen_moab_map, Tag th = NULL);

  enum tetGenNodesTypes {
    RIDGEVERTEX = 0,
    FREESEGVERTEX = 1,
    FREEFACETVERTEX = 2,
    FREEVOLVERTEX = 3
  };

  /** \brief set point tags and type

  Set type of entity, look in TetGen manual for details

  \code
  std::map<int,Range> types_ents;
  //RIDGEVERTEX
  types_ents[TetGenInterface::RIDGEVERTEX].merge(region_tets_skin_without_boundary_nodes);
  //FREESEGVERTEX
  types_ents[TetGenInterface::FREESEGVERTEX].merge(crack_surface_tris_skin_nodes);
  //FREEFACETVERTEX
  types_ents[TetGenInterface::FREEFACETVERTEX].merge(region_tets_skin_nodes);
  types_ents[TetGenInterface::FREEFACETVERTEX] =
  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEFACETVERTEX] =
  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  //FREEVOLVERTEX
  types_ents[TetGenInterface::FREEVOLVERTEX].merge(region_nodes);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
  subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
  subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] =
  subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREEFACETVERTEX]);
  \endcode

  */
  MoFEMErrorCode setGeomData(tetgenio &in, moabTetGen_Map &moab_tetgen_map,
                             tetGenMoab_Map &tetgen_moab_map,
                             std::map<int, Range> &type_ents);

  /** \brief get entities for TetGen data structure

    \param ents range of entities (tetrahedrons or nodes)
    \param in tetgen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to moab entities
    \param ents rerun entities which are in TetGen Data structure
    \param id_in_tags use tags as entity handles, if that is a case use tag to
    find MoAB vertex id
    \param error_if_created throw error if node need to be created

    */
  MoFEMErrorCode outData(tetgenio &in, tetgenio &out,
                         moabTetGen_Map &moab_tetgen_map,
                         tetGenMoab_Map &tetgen_moab_map, Range *ents = NULL,
                         bool id_in_tags = false, bool error_if_created = false,
                         bool assume_first_nodes_the_same = false,
                         Tag th = nullptr);

  /** \brief get entities for TetGen data structure

    \param ents range of entities (tetrahedrons or nodes)
    \param in tetgen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities
    \param ents rerun entities which are in TetGen data structure
    \param bit set level to created entities
    \param error_if_created throw error if node need to be created

    */
  MoFEMErrorCode outData(tetgenio &in, tetgenio &out,
                         moabTetGen_Map &moab_tetgen_map,
                         tetGenMoab_Map &tetgen_moab_map, BitRefLevel bit,
                         bool id_in_tags = false, bool error_if_created = false,
                         bool assume_first_nodes_the_same = false,
                         Tag th = nullptr);

  /** \brief set markers to faces

    \param markers data structure with markers
    \param in tetgen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities

    */
  MoFEMErrorCode setFaceData(std::vector<std::pair<Range, int> > &markers,
                             tetgenio &in, moabTetGen_Map &moab_tetgen_map,
                             tetGenMoab_Map &tetgen_moab_map);

  /** \brief get markers to faces

    \param markers data structure with markers
    \param in tetgen data structure (look to TetGen user manual)
    \param moab_tetgen_map mapping MoAB to TetGen entities
    \param tetgen_moab_map mapping TetGen to MoAB entities

    */
  MoFEMErrorCode getTriangleMarkers(tetGenMoab_Map &tetgen_moab_map,
                                    tetgenio &out, Range *ents = NULL,
                                    idxRange_Map *ents_map = NULL,
                                    bool only_non_zero = true);

  /** \brief set region data to tetrahedral
    */
  MoFEMErrorCode
  setRegionData(std::vector<std::pair<EntityHandle, int> > &regions, tetgenio &in,
               Tag th = NULL);

  /** \brief get region data to tetrahedral
    */
  MoFEMErrorCode getRegionData(tetGenMoab_Map &tetgen_moab_map, tetgenio &out,
                              Range *ents = NULL,
                              idxRange_Map *ents_map = NULL);

  /** \brief run tetgen
    */
  MoFEMErrorCode tetRahedralize(char switches[], tetgenio &in, tetgenio &out);

  /** \brief load poly file
    */
  MoFEMErrorCode loadPoly(char file_name[], tetgenio &in);

  // Tools for TetGen, i.e. geometry reconstruction from mesh

  MoFEMErrorCode checkPlanar_Trinagle(double coords[], bool *result,
                                      const double eps = 1e-9);
  MoFEMErrorCode groupPlanar_Triangle(Range &tris, std::vector<Range> &sorted,
                                      const double eps = 1e-9, Tag th = NULL);

  /** \brief Group surface triangles in planar regions

    \param tris input triangles
    \param sorted output sorted planar faces
    \param eps tolerance

  */
  MoFEMErrorCode groupRegion_Triangle(Range &tris,
                                      std::vector<std::vector<Range> > &sorted,
                                      const double eps = 1e-9);

  /** make planar polygon facet

    \param ents surface triangles
    \param polygons output list of polygons
    \param reduce_edges reduce edges if on the line
    \param not_reducable_nodes do not reduce node on edge if in this range
    \param eps tolerance

    \bug assumes that are no holes

    */
  MoFEMErrorCode makePolygonFacet(Range &ents, Range &polygons,
                                  bool reduce_edges = false,
                                  Range *not_reducable_nodes = NULL,
                                  const double eps = 1e-9,Tag th = NULL);
};
}

#endif //__TETGENINTERFACE_HPP__

/**
 * \defgroup mesh_tetgen TetGen interface
 * \brief Interface to run TetGen
 *
 * \ingroup mofem
 */
