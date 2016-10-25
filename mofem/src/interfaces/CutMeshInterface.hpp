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
    PetscErrorCode findEdgesToCut(int verb = 0);

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
    PetscErrorCode moveMidNodesOnCutEdges();

    /**
     * \brief Find edges to trimEdges

     * To make this work, you need to find edges to cut (findEdgesToCut), then
     * cut edges in the middle (cutEdgesInMiddle) and finally project edges on
     * the surface (moveMidNodesOnCutEdges)

     * @param  verb verbosity level
     * @return      error code
     */
    PetscErrorCode findEdgesToTrim(int verb = 0);

    PetscErrorCode trimEdgesInTheMiddle(const BitRefLevel bit);



    // PetscErrorCode findTetOnTheFront(int verb = 0);
    // #ifdef WITH_TETGEN
    // PetscErrorCode imprintFront(
    //   const BitRefLevel bit,
    //   const double tetgen_face_angle,
    //   int verb = 0
    // );
    // #endif //WITH_TETGEN


    inline const Range& getCutEdges() const { return cutEdges; }
    inline const Range& getCutEdgesOutside() const { return cutEdgesOutside; }
    inline const Range& getCutVolumes() const { return cutVolumes; }
    inline const Range& getNewCutVolumes() const { return cutNewVolumes; }
    inline const Range& getNewCutSurfaces() const { return cutNewSurfaces; }
    inline const Range& getNewCutVertices() const { return cutNewVertices; }
    // inline const Range& getFrontTets() const { return frontTets; }

    inline const Range& getTrimEdges() const { return trimEdges; }
    inline const Range& getOutsideEdges() const { return outsideEdges; }

    inline const Range& getNewTrimVolumes() const { return trimNewVolumes; }
    inline const Range& getNewTrimSurfaces() const { return trimNewSurfaces; }
    inline const Range& getNewTrimVertices() const { return trimNewVertices; }

  private:

    Range sUrface;
    Range vOlume;

    boost::shared_ptr<OrientedBoxTreeTool> treeSurfPtr;
    boost::shared_ptr<OrientedBoxTreeTool> treeVolPtr;
    EntityHandle rootSetSurf;
    EntityHandle rootSetVol;

    Range verticesOnSurface;

    Range cutEdges;
    Range cutEdgesOutside;
    Range cutVolumes;
    Range cutNewVolumes;
    Range cutNewSurfaces;
    Range cutNewVertices;

    Range trimNewVolumes;
    Range trimNewVertices;
    Range trimNewSurfaces;


    Range trimEdges;
    Range outsideEdges;

    // Range frontTets;

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

    int segmentPlane(
      VectorAdaptor s0,
      VectorAdaptor s1,
      VectorAdaptor x0,
      VectorAdaptor n,
      double &s
    ) const;


  };

}

#endif //__CUTMESHINTERFACE_HPP__

/***************************************************************************//**
 * \defgroup mesh_cut Mech cutter
 * \brief Interface to mesh cut mesh
 *
 * \ingroup mofem
 ******************************************************************************/
