/** \file PrismsFromSurface.hpp
 * \brief PrismsFromSurface interface
 *
 * Create prisms from surface triangle elements
 */

/* MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#ifndef __PRISMS_FORM_SURFACE_HPP__
#define __PRISMS_FORM_SURFACE_HPP__

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMPrismsFromSurface =
    MOFEMuuid(BitIntefaceId(PRISMSFROMSURFACE_INTERFACE));

/** \brief merge node from two bit levels
 * \ingroup mofem
 */
struct PrismsFromSurfaceInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  PrismsFromSurfaceInterface(const MoFEM::Core &core)
      : cOre(const_cast<MoFEM::Core &>(core)) {}

  std::map<EntityHandle, EntityHandle> createdVertices;

  /**
   * \brief Make prisms from triangles
   * @param  ents       Range of triangles
   * @param  swap_nodes If set true prism's nodes are swapped to satisfy the
   * canonical ordering (required if surface normal is pointing inwards)
   * @param  prisms     Returned range of prisms
   * @param  verb       Verbosity level
   * @return            Error code
   */
  MoFEMErrorCode createPrisms(const Range &ents, const bool swap_nodes,
                              Range &prisms, int verb = -1);

  /// \deprecated Use the function with the same name and a bool parameter
  /// *swap_nodes*, if set true prism's nodes are swapped to satisfy the
  /// canonical ordering (required if surface normal is pointing inwards)
  DEPRECATED MoFEMErrorCode createPrisms(const Range &ents, Range &prisms,
                                         int verb = -1);
  /**
   * \brief Seed prism entities by bit level
   * @param  prisms Range of entities
   * @param  bit    BitRefLevel
   * @param  verb   Verbosity level
   * @return        Error code
   */
  MoFEMErrorCode seedPrismsEntities(Range &prisms, const BitRefLevel &bit,
                                    int verb = -1);

  /**
   * \brief Make prisms by extruding top or bottom prisms
   * @param  prisms      Input prisms
   * @param  from_down  Use top or down face, if true from f3
   * @param  out_prisms  Returned prisms entities
   * @param  verb        Verbosity level
   * @return             Error code
   */
  MoFEMErrorCode createPrismsFromPrisms(const Range &prisms, bool from_down,
                                        Range &out_prisms, int verb = -1);

  /**
   * Set uniform thickness
   * @param  prisms   Range of prisms
   * @param  director3 Displacement of face 3
   * @param  director4 Displacement of face 4
   * @return
   */
  MoFEMErrorCode setThickness(const Range &prisms, const double director3[],
                              const double director4[]);

  /**
   * Set normal thickness
   * @param prisms   Range of prisms
   * @param thickness normal thickness
   * @return
   */
  MoFEMErrorCode setNormalThickness(const Range &prisms, double thickness3,
                                    double thickness4);

  /**
   * @brief Add quads to bockset
   *
   * If quad is adjacent to extruded edge, is added to given blockset
   *
   * @param prisms
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode updateMeshestByEdgeBlock(const Range &prisms);

  /**
   * @brief Add prism to bockset
   *
   * If quad is adjacent to extruded triangle, is added to given blockset
   *
   * @param prisms
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode updateMeshestByTriBlock(const Range &prisms);
};

} // namespace MoFEM

#endif //__PRISMS_FORM_SURFACE_HPP__
