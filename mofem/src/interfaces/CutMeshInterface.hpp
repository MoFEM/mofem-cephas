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

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  CutMeshInterface(const MoFEM::Core &core);
  ~CutMeshInterface() {}

  int lineSearchSteps;
  int nbMaxMergingCycles;
  int nbMaxTrimSearchIterations;

  /**
   * \brief Get options from command line
   * @return error code
   */
  MoFEMErrorCode getOptions() {
    MoFEMFunctionBegin;
    CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "MOFEM Cut mesh options",
                             "none");
    
    CHKERR PetscOptionsInt("-cut_linesearch_steps",
                           "number of bisection steps which line search do to "
                           "find optimal merged nodes position",
                           "", lineSearchSteps, &lineSearchSteps, PETSC_NULL);
    
    CHKERR PetscOptionsInt("-cut_max_merging_cycles",
                           "number of maximal merging cycles", "",
                           nbMaxMergingCycles, &nbMaxMergingCycles, PETSC_NULL);
    
    CHKERR PetscOptionsInt(
        "-cut_max_trim_iterations", "number of maximal merging cycles", "",
        nbMaxTrimSearchIterations, &nbMaxTrimSearchIterations, PETSC_NULL);

    ierr = PetscOptionsEnd();
    CHKERRG(ierr);
    MoFEMFunctionReturn(0);
  }

  /**
   * \brief set surface entities
   * @param  surface entities which going to be added
   * @return         error code
   */
  MoFEMErrorCode setSurface(const Range &surface);

  /**
   * \brief copy surface entities
   * @param  surface entities which going to be added
   * @return         error code
   */
  MoFEMErrorCode copySurface(const Range &surface, Tag th = NULL,
                             double *shift = NULL, double *origin = NULL,
                             double *transform           = NULL,
                             const std::string save_mesh = "");

  /**
   * \brief set volume entities
   * @param  volume entities which going to be added
   * @return         error code
   */
  MoFEMErrorCode setVolume(const Range &volume);

  /**
   * \brief merge surface entities
   * @param  surface entities which going to be added
   * @return         error code
   */
  MoFEMErrorCode mergeSurface(const Range &surface);

  /**
   * \brief merge volume entities
   * @param  volume entities which going to be added
   * @return         error code
   */
  MoFEMErrorCode mergeVolumes(const Range &volume);

  /**
   * \brief build tree
   * @return error code
   */
  MoFEMErrorCode buildTree();

  MoFEMErrorCode
  cutAndTrim(const BitRefLevel &bit_level1, const BitRefLevel &bit_level2,
             Tag th, const double tol_cut, const double tol_cut_close,
             const double tol_trim, const double tol_trim_close,
             Range *fixed_edges = NULL, Range *corner_nodes = NULL,
             const bool update_meshsets = false, const bool debug = true);

  MoFEMErrorCode
  cutTrimAndMerge(const int fraction_level, const BitRefLevel &bit_level1,
                  const BitRefLevel &bit_level2, const BitRefLevel &bit_level3,
                  Tag th, const double tol_cut, const double tol_cut_close,
                  const double tol_trim, const double tol_trim_close,
                  Range &fixed_edges, Range &corner_nodes,
                  const bool update_meshsets = false, const bool debug = false);

  /**
   * \brief find edges to cut
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode findEdgesToCut(Range *fixed_edges, Range *corner_nodes,
                                const double low_tol = 0, int verb = 0,
                                const bool debug = false);

  MoFEMErrorCode projectZeroDistanceEnts(Range *fixed_edges,
                                         Range *corner_nodes,
                                         const double low_tol = 0,
                                         const int verb = QUIET,
                                         const bool debug = false);

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
  MoFEMErrorCode cutEdgesInMiddle(const BitRefLevel bit,const bool debug = false);

  /**
   * \brief projecting of mid edge nodes on new mesh on surface
   * @return error code
   */
  MoFEMErrorCode moveMidNodesOnCutEdges(Tag th = NULL);

  /**
   * \brief Find edges to trimEdges

   * To make this work, you need to find edges to cut (findEdgesToCut), then
   * cut edges in the middle (cutEdgesInMiddle) and finally project edges on
   * the surface (moveMidNodesOnCutEdges)

   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode findEdgesToTrim(Range *fixed_edges, Range *corner_nodes,
                                 Tag th = NULL, const double tol = 1e-4,
                                 int verb = 0);

  /**
   * \brief trim edges
   * @param  bit bit level of the trimmed mesh
   * @return     error code
   */
  MoFEMErrorCode trimEdgesInTheMiddle(const BitRefLevel bit, Tag th = NULL,
                                      const double tol = 1e-4,
                                      const bool debug = false);

  /**
   * \brief move trimmed edges mid nodes
   * @return error code
   */
  MoFEMErrorCode moveMidNodesOnTrimmedEdges(Tag th = NULL);

  /**
   * \brief Remove pathological elements on surface internal front
   *
   * Internal surface skin is a set of edges in iterior of the body on boundary
   * of surface. This set of edges is called surface front. If surface face has
   * three nodes on surface front, non of the face nodes is split and should be
   * removed from surface if it is going to be split.
   *
   * @param  split_bit split bit level
   * @param  bit       bit level of split mesh
   * @param  ents      ents on the surface which is going to be split
   * @return           error code
   */
  MoFEMErrorCode removePathologicalFrontTris(const BitRefLevel split_bit,
                                             Range &ents);

  /**
   * \brief split sides
   * @param  split_bit split bit level
   * @param  bit       bit level of split mesh
   * @param  ents      ents on the surface which is going to be split
   * @return           error code
   */
  MoFEMErrorCode splitSides(const BitRefLevel split_bit, const BitRefLevel bit,
                            const Range &ents, Tag th = NULL);

  /**
   * @brief Merge edges
   *
   * Sort all edges, where sorting is by quality calculated as edge length times
   * quality of tets adjacent to the edge. Edge is merged if quality if the mesh
   * is improved.
   *
   * @param fraction_level Fraction of edges attemt to be merged at iteration
   * @param tets Tets of the mesh which edges are merged
   * @param surface Surface created by edge spliting
   * @param fixed_edges edges which are geometrical corners of the body
   * @param corner_nodes vertices on the corners
   * @param merged_nodes  merged nodes
   * @param out_tets  returned test after merge
   * @param new_surf  new surface without merged edges
   * @param th  tag with nodal positons
   * @param bit_ptr set bit ref level to mesh without merged edges
   * @param debug
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode mergeBadEdges(const int fraction_level, const Range &tets,
                               const Range &surface, const Range &fixed_edges,
                               const Range &corner_nodes, Range &merged_nodes,
                               Range &out_tets, Range &new_surf, Tag th,
                               const bool update_meshsets = false,
                               const BitRefLevel *bit_ptr = NULL,
                               const bool debug = false);

  /**
   * @brief Merge edges
   *
   * Sort all edges, where sorting is by quality calculated as edge length times
   * quality of tets adjacent to the edge. Edge is merged if quality if the mesh
   * is improved.
   */
  MoFEMErrorCode mergeBadEdges(const int fraction_level,
                               const BitRefLevel cut_bit,
                               const BitRefLevel trim_bit,
                               const BitRefLevel bit, const Range &surface,
                               const Range &fixed_edges,
                               const Range &corner_nodes, Tag th = NULL,
                                const bool update_meshsets = false,
                               const bool debug = false);

#ifdef WITH_TETGEN

  MoFEMErrorCode
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
  MoFEMErrorCode setTagData(Tag th,const BitRefLevel bit = BitRefLevel());

  /**
   * \brief set coords from tag
   * @param  th tag handle
   * @return    error code
   */
  MoFEMErrorCode setCoords(Tag th, const BitRefLevel bit = BitRefLevel(),
                           const BitRefLevel mask = BitRefLevel().set());

  inline const Range &getVolume() const { return vOlume; }
  inline const Range &getSurface() const { return sUrface; }

  inline const Range &getCutEdges() const { return cutEdges; }
  inline const Range &getCutVolumes() const { return cutVolumes; }
  inline const Range &getNewCutVolumes() const { return cutNewVolumes; }
  inline const Range &getNewCutSurfaces() const { return cutNewSurfaces; }
  inline const Range &getNewCutVertices() const { return cutNewVertices; }
  inline const Range &projectZeroDistanceEnts() const { return zeroDistanceEnts; }

  inline const Range &getTrimEdges() const { return trimEdges; }
  inline const Range &getNewTrimVolumes() const { return trimNewVolumes; }
  inline const Range &getNewTrimSurfaces() const { return trimNewSurfaces; }
  inline const Range &getNewTrimVertices() const { return trimNewVertices; }

  inline const Range &getMergedVolumes() const { return mergedVolumes; }
  inline const Range &getMergedSurfaces() const { return mergedSurfaces; }

  inline const Range &getTetgenSurfaces() const { return tetgenSurfaces; }

  MoFEMErrorCode saveCutEdges();

  MoFEMErrorCode saveTrimEdges();

  inline boost::shared_ptr<OrientedBoxTreeTool>& getTreeSurfPtr() {
    return treeSurfPtr;
  }

private:
  Range sUrface;
  Range vOlume;

  boost::shared_ptr<OrientedBoxTreeTool> treeSurfPtr;
  EntityHandle rootSetSurf;

  Range cutEdges;
  Range cutVolumes;
  Range cutNewVolumes;
  Range cutNewSurfaces;
  Range zeroDistanceEnts;
  Range zeroDistanceVerts;
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
  map<EntityHandle, TreeData> verticesOnCutEdges;
  map<EntityHandle, TreeData> edgesToTrim;
  map<EntityHandle, TreeData> verticesOnTrimEdges;

#ifdef WITH_TETGEN

  map<EntityHandle, unsigned long> moabTetGenMap;
  map<unsigned long, EntityHandle> tetGenMoabMap;
  boost::ptr_vector<tetgenio> tetGenData;

#endif

  MoFEMErrorCode getRayForEdge(const EntityHandle ent, VectorAdaptor &ray_point,
                               VectorAdaptor &unit_ray_dir,
                               double &ray_length) const;

  // /**
  //  * Find if segment in on the plain
  //  * @param  s0 segment first point
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

  double aveLength; ///< Average edge length
  double maxLength; ///< Maximal edge length
  
};
} // namespace MoFEM

#endif //__CUTMESHINTERFACE_HPP__

/**
 * \defgroup mesh_cut CutMeshInterface
 * \brief Interface to mesh cut mesh
 *
 * \ingroup mofem
 */
