/** \file GetAdjacancies.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
 *
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

#include <version.h>
#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

  PetscErrorCode Core::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return Tools(*this).getAdjacenciesEquality(from_entiti,to_dimension,adj_entities);
  }
  PetscErrorCode Core::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return Tools(*this).getAdjacenciesAny(from_entiti,to_dimension,adj_entities);
  }
  PetscErrorCode Core::get_adjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return Tools(*this).getAdjacencies(problem_ptr,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }
  PetscErrorCode Core::get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return Tools(*this).getAdjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

}
