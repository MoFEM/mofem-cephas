/** \file HangingNodes.hpp
 * \brief Mesh refinemnt with hanging nodes
 */

/*
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

#ifndef __HANGINGNODES_HPP__
#define __HANGINGNODES_HPP__

namespace MoFEM {

struct HangingNodes : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  HangingNodes(const MoFEM::Core &core);

  virtual ~HangingNodes() = default;

  /**
   * \brief Make vertices in the middle of edges in meshset and add them to
   * Refinement levels defined by bit
   *
   * \param Range consisting edges for refine
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried
   * recursively.  Returns the contents of meshsets, but not the meshsets
   * themselves if true.
   */
  MoFEMErrorCode addVerticesInTheMiddleOfEdges(const Range &edges,
                                               const BitRefLevel &bit,
                                               EntityHandle start_v = 0);

	/**
	 * @brief Refine elements
	 * 
	 * @param ents 
	 * @param bit 
	 * @param verb 
	 * @param start_v 
	 * @return MoFEMErrorCode 
	 */
  MoFEMErrorCode refineMesh(const Range &ents, const BitRefLevel &bit,
                            EntityHandle start_v = 0);

private:
  MoFEM::Core &cOre;
}

} // namespace MoFEM

#endif // __HANGINGNODES_HPP__

/**
 * \defgroup mofem_refiner MeshRefinement
 * \brief Refine mesh by splitting edges
 *
 */