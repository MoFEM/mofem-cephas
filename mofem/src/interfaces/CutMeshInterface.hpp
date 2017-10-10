/** \file CutMeshInterface.hpp
 * \brief Cut mesh interface
 *
 * \ingroup mesh_cut
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

#ifndef __CUTMESHINTERFACE_HPP__
#define __CUTMESHINTERFACE_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

  static const MOFEMuuid IDD_MOFEMCutMesh = MOFEMuuid( BitIntefaceId(CUTMESH_INTERFACE) );

  /**
   *  \brief Interface to cut meshes
   *
   * \ingroup mesh_cut
   */
  struct CutMeshInterface: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MoFEM::Core& cOre;
    CutMeshInterface(const MoFEM::Core &core);
    ~CutMeshInterface() {}

    /**
     * \brief set surface entities
     * @param  surface entities which going to be added
     * @return         error code
     */
    PetscErrorCode setSurface(const Range &surface);

    /**
     * \brief copy surface entities
     * @param  surface entities which going to be added
     * @return         error code
     */
    PetscErrorCode copySurface(
      const Range &surface,
      Tag th = NULL,
      double *shift = NULL,
      double *origin = NULL,
      double *transform = NULL
    );

    /**
     * \brief set volume entities
     * @param  volume entities which going to be added
     * @return         error code
     */
    PetscErrorCode setVolume(const Range &volume);

    /**
     * \brief merge surface entities
     * @param  surface entities which going to be added
     * @return         error code
     */
    PetscErrorCode mergeSurface(const Range &surface);

    /**
     * \brief merge volume entities
     * @param  volume entities which going to be added
     * @return         error code
     */
    PetscErrorCode mergeVolume(const Range &volume);

    /**
     * \brief build tree
     * @return error code
     */
    PetscErrorCode buildTree();

    /**
     * \brief find edges to cut
     * @param  verb verbosity level
     * @return      error code
     */
    PetscErrorCode findEdgesToCut(const double low_tol = 0,int verb = 0);

    PetscErrorCode getEntsOnCutSurface(
      const double low_tol = 0,int verb  = 0
    );

    /**
     * \brief cut edges
     *
     * For edges to cut (calculated by findEdgesToCut), edges are split in the
     * middle and then using MoFEM::MeshRefinement interface, tetrahedra mesh
     * are cut.
     *
     * @param  bit BitRefLevel of new mesh created by cutting edges
     * @return     error code
     */
    PetscErrorCode cutEdgesInMiddle(const BitRefLevel bit);

    /**
     * \brief projecting of mid edge nodes on new mesh on surface
     * @return error code
     */
    PetscErrorCode moveMidNodesOnCutEdges(Tag th = NULL);

    /**
     * \brief Find edges to trimEdges

     * To make this work, you need to find edges to cut (findEdgesToCut), then
     * cut edges in the middle (cutEdgesInMiddle) and finally project edges on
     * the surface (moveMidNodesOnCutEdges)

     * @param  verb verbosity level
     * @return      error code
     */
    PetscErrorCode findEdgesToTrim(Tag th = NULL,const double tol = 1e-4,int verb = 0);

    /**
     * \brief trim edges
     * @param  bit bit level of the trimed mesh
     * @return     error code
     */
    PetscErrorCode trimEdgesInTheMiddle(const BitRefLevel bit,Tag th = NULL,const double tol = 1e-4);

    /**
     * \brief move trimed edges mid nodes
     * @return error code
     */
    PetscErrorCode moveMidNodesOnTrimedEdges(Tag th = NULL);

    /**
     * \brief split sides
     * @param  split_bit split bit level
     * @param  bit       bit level of split mesh
     * @param  ents      ents on the surface which is going to be split
     * @return           error code
     */
    PetscErrorCode splitSides(
      const BitRefLevel split_bit,
      const BitRefLevel bit,
      const Range &ents,
      Tag th = NULL
    );

    /**
     * \brief split sides of trimmed surface
     * @param  split_bit split bit level
     * @param  bit       bit level of split mesh
     * @return           error code
     */
    PetscErrorCode splitTrimSides(
      const BitRefLevel split_bit,
      const BitRefLevel bit,
      Tag th = NULL
    );

    PetscErrorCode mergeBadEdgesOnSurface(
      const int fraction_level,
      const Range& tets,const Range& surface,const Range& fixed_verts,
      Tag th_quality,Tag th_position,
      Range& out_tets,Range& new_surf
    );

    inline const Range& getCutEdges() const { return cutEdges; }
    inline const Range& getCutVolumes() const { return cutVolumes; }
    inline const Range& getNewCutVolumes() const { return cutNewVolumes; }
    inline const Range& getNewCutSurfaces() const { return cutNewSurfaces; }
    inline const Range& getNewCutVertices() const { return cutNewVertices; }
    inline const Range& getZeroDistanceEnts() const { return zeroDistanseEnts; }


    inline const Range& getTrimEdges() const { return trimEdges; }

    inline const Range& getNewTrimVolumes() const { return trimNewVolumes; }
    inline const Range& getNewTrimSurfaces() const { return trimNewSurfaces; }
    inline const Range& getNewTrimVertices() const { return trimNewVertices; }

  private:

    Range sUrface;
    Range vOlume;

    boost::shared_ptr<OrientedBoxTreeTool> treeSurfPtr;
    EntityHandle rootSetSurf;

    Range cutEdges;
    Range cutVolumes;
    Range cutNewVolumes;
    Range cutNewSurfaces;
    Range zeroDistanseEnts;
    Range zeroDistanseVerts;
    Range cutNewVertices;

    Range trimNewVolumes;
    Range trimNewVertices;
    Range trimNewSurfaces;

    Range trimEdges;

    struct TreeData {
      double dIst;
      double lEngth;
      VectorDouble3 unitRayDir;
      VectorDouble3 rayPoint;
    };

    map<EntityHandle,TreeData> edgesToCut;
    map<EntityHandle,TreeData> verticecOnCutEdges;
    map<EntityHandle,TreeData> edgesToTrim;
    map<EntityHandle,TreeData> verticecOnTrimEdges;

    // #ifdef WITH_TETGEN
    //
    // map<EntityHandle,unsigned long> moabTetGenMap;
    // map<unsigned long,EntityHandle> tetGenMoabMap;
    // boost::ptr_vector<tetgenio> tetGenData;
    //
    // #endif

    PetscErrorCode getRayForEdge(
      const EntityHandle ent,
      VectorAdaptor ray_point,
      VectorAdaptor unit_ray_dir,
      double &ray_length
    ) const;

    // /**
    //  * Find if segment in on the plain
    //  * @param  s0 segemnt fisrt point
    //  * @param  s1 segment second point
    //  * @param  x0 point on the plain
    //  * @param  n  normal on the plain
    //  * @param  s  intersect point
    //  * @return    1 - intersect, 2 - segment on the plain, 0 - no intersect
    //  */
    // int segmentPlane(
    //   VectorAdaptor s0,
    //   VectorAdaptor s1,
    //   VectorAdaptor x0,
    //   VectorAdaptor n,
    //   double &s
    // ) const;

    PetscErrorCode mergeNodes(
      EntityHandle father,
      EntityHandle mother,
      Range& proc_tets,
      Range &new_surf,
      const bool only_if_improve_quality = false,
      const double move = 0,
      const int line_search = 0,
      Tag th = NULL,
      const int verb = 0
    );

    double aveLength;

  };

}

#endif //__CUTMESHINTERFACE_HPP__

/***************************************************************************//**
 * \defgroup mesh_cut Mech cutter
 * \brief Interface to mesh cut mesh
 *
 * \ingroup mofem
 ******************************************************************************/
