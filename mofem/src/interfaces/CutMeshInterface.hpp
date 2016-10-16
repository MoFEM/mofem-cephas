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

    PetscErrorCode setSurface(const Range &surface);
    PetscErrorCode setVolume(const Range &volume);
    PetscErrorCode mergeSurface(const Range &surface);
    PetscErrorCode mergeVolume(const Range &volume);

    PetscErrorCode buildTree();
    PetscErrorCode findToCut(int verb = 0);
    PetscErrorCode findTetOnTheFront(int verb = 0);

    inline const Range& getCutEdges() const { return cutEdges; }
    inline const Range& getCutVolumes() const { return cutVolumes; }

    #ifdef WITH_TETGEN
    PetscErrorCode imprintFront(
      const BitRefLevel bit,
      const double tetgen_face_angle,
      int verb = 0
    );
    #endif //WITH_TETGEN

    PetscErrorCode cutTets(const BitRefLevel bit);

    PetscErrorCode moveNodes();

    inline const Range& getNewCutVolumes() const { return cutNewVolumes; }
    inline const Range& getNewCutSurfaces() const { return cutNewSurfaces; }
    inline const Range& getNewCutVertices() const { return cutNewVertices; }
    inline const Range& getFrontTets() const { return frontTets; }

  private:

    Range sUrface;
    Range vOlume;

    boost::shared_ptr<OrientedBoxTreeTool> obTreeSurfPtr;
    EntityHandle rootSetSurf;
    boost::shared_ptr<OrientedBoxTreeTool> obTreeVolPtr;
    EntityHandle rootSetVol;

    Range cutEdges;
    Range cutVolumes;
    Range cutNewVolumes;
    Range cutNewSurfaces;
    Range cutNewVertices;

    Range frontTets;

    struct EdgeData {
      double dIst;
      double lEngth;
      VectorDouble3 unitRayDir;
      VectorDouble3 rayPoint;
    };

    map<EntityHandle,EdgeData> edgesToCut;
    map<EntityHandle,EdgeData> verticecOnCutEdges;
    map<EntityHandle,EdgeData> facesFrontToRemesh;

    PetscErrorCode getRayForEdge(
      const EntityHandle ent,
      double *ray_point,
      double *unit_ray_dir,
      double *ray_length
    ) const;

    #ifdef WITH_TETGEN

    map<EntityHandle,unsigned long> moabTetGenMap;
    map<unsigned long,EntityHandle> tetGenMoabMap;
    boost::ptr_vector<tetgenio> tetGenData;

    #endif


  };

}

#endif //__CUTMESHINTERFACE_HPP__

/***************************************************************************//**
 * \defgroup mesh_cut Mech cutter
 * \brief Interface to mesh cut mesh
 *
 * \ingroup mofem
 ******************************************************************************/
