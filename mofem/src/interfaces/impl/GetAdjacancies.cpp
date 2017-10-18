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

namespace MoFEM {

  PetscErrorCode Core::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return BitRefManager(*this).getAdjacenciesEquality(from_entiti,to_dimension,adj_entities);
  }
  PetscErrorCode Core::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return BitRefManager(*this).getAdjacenciesAny(from_entiti,to_dimension,adj_entities);
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
    return BitRefManager(*this).getAdjacencies(problem_ptr,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
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
    return BitRefManager(*this).getAdjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

}
