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

static const MOFEMuuid IDD_MOFEMCutMesh =
    MOFEMuuid(BitIntefaceId(CUTMESH_INTERFACE));

/**
 *  \brief Interface to cut meshes
 *
 * \ingroup mesh_cut
 */
struct CutMeshInterface : public UnknownInterface {

  PetscErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  CutMeshInterface(const MoFEM::Core &core);
  ~CutMeshInterface() {}

  int lineSearchSteps;
  int nbMaxMergingCycles;
  int nbMaxTrimSearchIterations;

  /**
   * \brief Get options from command line
   * @return error cdoe
   */
  PetscErrorCode getOptions() {
    MoFEMFunctionBeginHot;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "MOFEM Cut mesh options",
                             "none");
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-cut_lineserach_steps",
                           "number of bisection steps wich line search do to "
                           "find optimal merged nodes position",
                           "", lineSearchSteps, &lineSearchSteps, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt("-cut_max_merging_cycles",
                           "number of maximal merging cycles", "",
                           nbMaxMergingCycles, &nbMaxMergingCycles, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsInt(
        "-cut_max_trim_iterations", "number of maximal merging cycles", "",
        nbMaxTrimSearchIterations, &nbMaxTrimSearchIterations, PETSC_NULL);
    CHKERRQ(ierr);
    ierr = PetscOptionsEnd();
    CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }

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
  PetscErrorCode copySurface(const Range &surface, Tag th = NULL,
                             double *shift = NULL, double *origin = NULL,
                             double *transform = NULL);

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
  PetscErrorCode mergeVolumes(const Range &volume);

  /**
   * \brief build tree
   * @return error code
   */
  PetscErrorCode buildTree();

  PetscErrorCode cutAndTrim(const BitRefLevel &bit_level1,
                            const BitRefLevel &bit_level2, Tag th,
                            const double tol_cut, const double tol_cut_close,
                            const double tol_trim, const double tol_trim_close,
                            Range *fixed_edges = NULL,
                            Range *corner_nodes = NULL,
                            const bool update_meshsets = false);

  PetscErrorCode
  cutTrimAndMerge(const int fraction_level, const BitRefLevel &bit_level1,
                  const BitRefLevel &bit_level2, const BitRefLevel &bit_level3,
                  Tag th, const double tol_cut, const double tol_cut_close,
                  const double tol_trim, const double tol_trim_close,
                  Range &fixed_edges, Range &corner_nodes,
                  const bool update_meshsets = false);

  /**
   * \brief find edges to cut
   * @param  verb verbosity level
   * @return      error code
   */
  PetscErrorCode findEdgesToCut(const double low_tol = 0, int verb = 0);

  PetscErrorCode getEntsOnCutSurface(const double low_tol = 0, int verb = 0);

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
  PetscErrorCode findEdgesToTrim(Tag th = NULL, const double tol = 1e-4,
                                 int verb = 0);

  /**
   * \brief trim edges
   * @param  bit bit level of the trimed mesh
   * @return     error code
   */
  PetscErrorCode trimEdgesInTheMiddle(const BitRefLevel bit, Tag th = NULL,
                                      const double tol = 1e-4);

  /**
   * \brief move trimed edges mid nodes
   * @return error code
   */
  PetscErrorCode moveMidNodesOnTrimedEdges(Tag th = NULL);

  /**
   * \brief Remove patalogical elements on surface internal front
   *
   * Internal surface skin is a set of edges in interia of the body on boundary
   * of surface. This set of edges is called surface front. If surface face has
   * three nodes on surface front, non of the face nodes is split and should be
   * removed from surface if it is going to be split.
   *
   * @param  split_bit split bit level
   * @param  bit       bit level of split mesh
   * @param  ents      ents on the surface which is going to be split
   * @return           error code
   */
  PetscErrorCode removePathologicalFrontTris(const BitRefLevel split_bit,
                                             Range &ents);

  /**
   * \brief split sides
   * @param  split_bit split bit level
   * @param  bit       bit level of split mesh
   * @param  ents      ents on the surface which is going to be split
   * @return           error code
   */
  PetscErrorCode splitSides(const BitRefLevel split_bit, const BitRefLevel bit,
                            const Range &ents, Tag th = NULL);

  PetscErrorCode mergeBadEdges(const int fraction_level, const Range &tets,
                               const Range &surface, const Range &fixed_edges,
                               const Range &corner_nodes, Range &merged_nodes,
                               Range &out_tets, Range &new_surf, Tag th,
                               const bool update_meshsets = false,
                               const BitRefLevel *bit_ptr = NULL);

  PetscErrorCode mergeBadEdges(const int fraction_level,
                               const BitRefLevel cut_bit,
                               const BitRefLevel trim_bit,
                               const BitRefLevel bit, const Range &surface,
                               const Range &fixed_edges,
                               const Range &corner_nodes, Tag th);

#ifdef WITH_TETGEN

  PetscErrorCode
  rebuildMeshWithTetGen(vector<string> &switches, const BitRefLevel &mesh_bit,
                        const BitRefLevel &bit, const Range &surface,
                        const Range &fixed_edges, const Range &corner_nodes,
                        Tag th = NULL, const bool debug = false);

#endif

  /**
   * \brief set coords to tag
   * @param  th tag handle
   * @return    error code
   */
  PetscErrorCode setTagData(Tag th);

  /**
   * \brief set coords from tag
   * @param  th tag handle
   * @return    error code
   */
  PetscErrorCode setCoords(Tag th);

  inline const Range &getCutEdges() const { return cutEdges; }
  inline const Range &getCutVolumes() const { return cutVolumes; }
  inline const Range &getNewCutVolumes() const { return cutNewVolumes; }
  inline const Range &getNewCutSurfaces() const { return cutNewSurfaces; }
  inline const Range &getNewCutVertices() const { return cutNewVertices; }
  inline const Range &getZeroDistanceEnts() const { return zeroDistanseEnts; }

  inline const Range &getTrimEdges() const { return trimEdges; }
  inline const Range &getNewTrimVolumes() const { return trimNewVolumes; }
  inline const Range &getNewTrimSurfaces() const { return trimNewSurfaces; }
  inline const Range &getNewTrimVertices() const { return trimNewVertices; }

  inline const Range &getMergedVolumes() const { return mergedVolumes; }
  inline const Range &getMergedSurfaces() const { return mergedSurfaces; }

  inline const Range &getTetgenSurfaces() const { return tetgenSurfaces; }

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

  Range mergedVolumes;
  Range mergedSurfaces;

  Range tetgenSurfaces;

  struct TreeData {
    double dIst;
    double lEngth;
    VectorDouble3 unitRayDir;
    VectorDouble3 rayPoint;
  };

  map<EntityHandle, TreeData> edgesToCut;
  map<EntityHandle, TreeData> verticecOnCutEdges;
  map<EntityHandle, TreeData> edgesToTrim;
  map<EntityHandle, TreeData> verticecOnTrimEdges;

#ifdef WITH_TETGEN

  map<EntityHandle, unsigned long> moabTetGenMap;
  map<unsigned long, EntityHandle> tetGenMoabMap;
  boost::ptr_vector<tetgenio> tetGenData;

#endif

  PetscErrorCode getRayForEdge(const EntityHandle ent, VectorAdaptor ray_point,
                               VectorAdaptor unit_ray_dir,
                               double &ray_length) const;

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

  double aveLength;
};
} // namespace MoFEM

#endif //__CUTMESHINTERFACE_HPP__

/**
 * \defgroup mesh_cut Mech cutter
 * \brief Interface to mesh cut mesh
 *
 * \ingroup mofem
 */
