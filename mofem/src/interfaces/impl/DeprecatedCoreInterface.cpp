/** \file DeprecatedCoreInterface.cpp
 * \brief implementation of deprecated functions
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

MoFEMErrorCode
DeprecatedCoreInterface::seed_ref_level_2D(const EntityHandle meshset,
                                           const BitRefLevel &bit, int verb) {
  return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset, 2, bit,
                                                            verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::seed_ref_level_3D(const EntityHandle meshset,
                                           const BitRefLevel &bit, int verb) {
  return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset, 3, bit,
                                                            verb);
}

MoFEMErrorCode DeprecatedCoreInterface::seed_ref_level(const Range &ents,
                                                       const BitRefLevel &bit,
                                                       const bool only_tets,
                                                       int verb) {
  return getInterface<BitRefManager>()->setBitRefLevel(ents, bit, only_tets,
                                                       verb);
}

MoFEMErrorCode DeprecatedCoreInterface::seed_ref_level_MESHSET(
    const EntityHandle meshset, const BitRefLevel &bit, int verb) {
  return getInterface<BitRefManager>()->setBitLevelToMeshset(meshset, bit,
                                                             verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::seed_finite_elements(const EntityHandle meshset,
                                              int verb) {

  MoFEMFunctionBegin;
  Range entities;
  CHKERR getInterface<CoreInterface>()->get_moab().get_entities_by_handle(
      meshset, entities, true);
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(
      entities, BitRefLevel(), verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
DeprecatedCoreInterface::seed_finite_elements(const Range &entities, int verb) {
  return getInterface<BitRefManager>()->setElementsBitRefLevel(
      entities, BitRefLevel(), verb);
}

MoFEMErrorCode DeprecatedCoreInterface::field_axpy(
    const double alpha, const std::string &field_name_x,
    const std::string &field_name_y, bool error_if_missing,
    bool creat_if_missing) {
  return getInterface<FieldBlas>()->fieldAxpy(
      alpha, field_name_x, field_name_y, error_if_missing, creat_if_missing);
}
MoFEMErrorCode
DeprecatedCoreInterface::set_field(const double val, const EntityType type,
                                   const std::string &field_name) {
  return getInterface<FieldBlas>()->setField(val, type, field_name);
}
MoFEMErrorCode
DeprecatedCoreInterface::set_field(const double val, const EntityType type,
                                   const Range &ents,
                                   const std::string &field_name) {
  return getInterface<FieldBlas>()->setField(val, type, ents, field_name);
}
MoFEMErrorCode
DeprecatedCoreInterface::field_scale(const double alpha,
                                     const std::string &field_name) {
  return getInterface<FieldBlas>()->fieldScale(alpha, field_name);
}

MoFEMErrorCode DeprecatedCoreInterface::get_adjacencies_equality(
    const EntityHandle from_entiti, const int to_dimension,
    Range &adj_entities) const {
  return getInterface<BitRefManager>()->getAdjacenciesEquality(
      from_entiti, to_dimension, adj_entities);
}

MoFEMErrorCode
DeprecatedCoreInterface::get_adjacencies_any(const EntityHandle from_entiti,
                                             const int to_dimension,
                                             Range &adj_entities) const {
  return getInterface<BitRefManager>()->getAdjacenciesAny(
      from_entiti, to_dimension, adj_entities);
}

MoFEMErrorCode DeprecatedCoreInterface::get_adjacencies(
    const Problem *problem_ptr, const EntityHandle *from_entities,
    const int num_netities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  return getInterface<BitRefManager>()->getAdjacencies(
      problem_ptr, from_entities, num_netities, to_dimension, adj_entities,
      operation_type, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::get_adjacencies(
    const BitRefLevel &bit, const EntityHandle *from_entities,
    const int num_netities, const int to_dimension, Range &adj_entities,
    const int operation_type, const int verb) const {
  return getInterface<BitRefManager>()->getAdjacencies(
      bit, from_entities, num_netities, to_dimension, adj_entities,
      operation_type, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::VecCreateSeq(const std::string &name,
                                                     RowColData rc,
                                                     Vec *V) const {
  return getInterface<VecManager>()->vecCreateSeq(name, rc, V);
}
MoFEMErrorCode DeprecatedCoreInterface::VecCreateGhost(const std::string &name,
                                                       RowColData rc,
                                                       Vec *V) const {
  return getInterface<VecManager>()->vecCreateGhost(name, rc, V);
}
MoFEMErrorCode DeprecatedCoreInterface::ISCreateProblemOrder(
    const std::string &problem, RowColData rc, int min_order, int max_order,
    IS *is, int verb) const {
  return getInterface<ISManager>()->isCreateProblemOrder(problem, rc, min_order,
                                                         max_order, is);
}
MoFEMErrorCode DeprecatedCoreInterface::ISCreateProblemFieldAndRank(
    const std::string &problem, RowColData rc, const std::string &field,
    int min_coeff_idx, int max_coeff_idx, IS *is, int verb) const {
  return getInterface<ISManager>()->isCreateProblemFieldAndRank(
      problem, rc, field, min_coeff_idx, max_coeff_idx, is);
}
MoFEMErrorCode
DeprecatedCoreInterface::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem, const std::string &x_field_name,
    RowColData x_rc, const std::string &y_problem,
    const std::string &y_field_name, RowColData y_rc, std::vector<int> &idx,
    std::vector<int> &idy, int verb) const {
  return getInterface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem, x_field_name, x_rc, y_problem, y_field_name, y_rc, idx, idy);
}

MoFEMErrorCode
DeprecatedCoreInterface::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem, const std::string &x_field_name,
    RowColData x_rc, const std::string &y_problem,
    const std::string &y_field_name, RowColData y_rc, IS *ix, IS *iy,
    int verb) const {
  return getInterface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem, x_field_name, x_rc, y_problem, y_field_name, y_rc, ix, iy);
}
MoFEMErrorCode DeprecatedCoreInterface::VecScatterCreate(
    Vec xin, const std::string &x_problem, const std::string &x_field_name,
    RowColData x_rc, Vec yin, const std::string &y_problem,
    const std::string &y_field_name, RowColData y_rc, VecScatter *newctx,
    int verb) const {
  return getInterface<VecManager>()->vecScatterCreate(
      xin, x_problem, x_field_name, x_rc, yin, y_problem, y_field_name, y_rc,
      newctx);
}
MoFEMErrorCode DeprecatedCoreInterface::VecScatterCreate(
    Vec xin, const std::string &x_problem, RowColData x_rc, Vec yin,
    const std::string &y_problem, RowColData y_rc, VecScatter *newctx,
    int verb) const {
  return getInterface<VecManager>()->vecScatterCreate(xin, x_problem, x_rc, yin,
                                                      y_problem, y_rc, newctx);
}
MoFEMErrorCode DeprecatedCoreInterface::set_local_ghost_vector(
    const Problem *problem_ptr, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  return getInterface<VecManager>()->setLocalGhostVector(problem_ptr, rc, V,
                                                         mode, scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_local_ghost_vector(
    const std::string &name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  return getInterface<VecManager>()->setLocalGhostVector(name, rc, V, mode,
                                                         scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_global_ghost_vector(
    const Problem *problem_ptr, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  return getInterface<VecManager>()->setGlobalGhostVector(problem_ptr, rc, V,
                                                          mode, scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_global_ghost_vector(
    const std::string &name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode) const {
  return getInterface<VecManager>()->setGlobalGhostVector(name, rc, V, mode,
                                                          scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_other_local_ghost_vector(
    const Problem *problem_ptr, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode, int verb) {
  return getInterface<VecManager>()->setOtherLocalGhostVector(
      problem_ptr, field_name, cpy_field_name, rc, V, mode, scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_other_local_ghost_vector(
    const std::string &name, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode, int verb) {
  return getInterface<VecManager>()->setOtherLocalGhostVector(
      name, field_name, cpy_field_name, rc, V, mode, scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_other_global_ghost_vector(
    const Problem *problem_ptr, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode, int verb) {
  return getInterface<VecManager>()->setOtherGlobalGhostVector(
      problem_ptr, field_name, cpy_field_name, rc, V, mode, scatter_mode);
}
MoFEMErrorCode DeprecatedCoreInterface::set_other_global_ghost_vector(
    const std::string &name, const std::string &field_name,
    const std::string &cpy_field_name, RowColData rc, Vec V, InsertMode mode,
    ScatterMode scatter_mode, int verb) {
  return getInterface<VecManager>()->setOtherGlobalGhostVector(
      name, field_name, cpy_field_name, rc, V, mode, scatter_mode);
}

MoFEMErrorCode DeprecatedCoreInterface::shift_right_bit_ref(const int shift,
                                                            int verb) {
  return getInterface<BitRefManager>()->shiftRightBitRef(
      shift, BitRefLevel().set(), verb);
}

MoFEMErrorCode DeprecatedCoreInterface::get_entities_by_type_and_ref_level(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    const EntityHandle meshset, int verb) {
  return getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, mask, type, meshset, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::get_entities_by_type_and_ref_level(
    const BitRefLevel &bit, const BitRefLevel &mask, const EntityType type,
    Range &ents, int verb) {
  return getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      bit, mask, type, ents, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::get_entities_by_ref_level(const BitRefLevel &bit,
                                                   const BitRefLevel &mask,
                                                   const EntityHandle meshset) {
  return getInterface<BitRefManager>()->getEntitiesByRefLevel(bit, mask,
                                                              meshset);
}

MoFEMErrorCode DeprecatedCoreInterface::get_entities_by_ref_level(
    const BitRefLevel &bit, const BitRefLevel &mask, Range &ents) {
  return getInterface<BitRefManager>()->getEntitiesByRefLevel(bit, mask, ents);
}

MoFEMErrorCode DeprecatedCoreInterface::update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,
    const EntityHandle child, EntityType child_type, const bool recursive,
    int verb) {
  return getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(
      parent, child_bit, child, child_type, recursive, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::update_field_meshset_by_entities_children(
    const BitRefLevel &child_bit, int verb) {
  return getInterface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(
      child_bit, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::update_field_meshset_by_entities_children(
    const std::string name, const BitRefLevel &child_bit, int verb) {
  return getInterface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(
      name, child_bit, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::update_finite_element_meshset_by_entities_children(
    const std::string name, const BitRefLevel &child_bit,
    const EntityType fe_ent_type, int verb) {
  return getInterface<BitRefManager>()
      ->updateFiniteElementMeshsetByEntitiesChildren(name, child_bit,
                                                     fe_ent_type, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::print_cubit_displacement_set() const {
  return getInterface<MeshsetsManager>()->printDisplacementSet();
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode DeprecatedCoreInterface::print_cubit_pressure_set() const {
  return getInterface<MeshsetsManager>()->printPressureSet();
}

MoFEMErrorCode DeprecatedCoreInterface::print_cubit_force_set() const {
  return getInterface<MeshsetsManager>()->printForceSet();
}

MoFEMErrorCode DeprecatedCoreInterface::print_cubit_materials_set() const {
  return getInterface<MeshsetsManager>()->printMaterialsSet();
}

bool DeprecatedCoreInterface::check_msId_meshset(
    const int ms_id, const CubitBCType cubit_bc_type) {
  return getInterface<MeshsetsManager>()->checkMeshset(ms_id, cubit_bc_type);
}

MoFEMErrorCode DeprecatedCoreInterface::add_cubit_msId(
    const CubitBCType cubit_bc_type, const int ms_id, const std::string name) {
  return getInterface<MeshsetsManager>()->addMeshset(cubit_bc_type, ms_id,
                                                     name);
}

MoFEMErrorCode DeprecatedCoreInterface::set_cubit_msId_attribites(
    const CubitBCType cubit_bc_type, const int ms_id,
    const std::vector<double> &attributes, const std::string name) {
  return getInterface<MeshsetsManager>()->setAtributes(cubit_bc_type, ms_id,
                                                       attributes, name);
}
MoFEMErrorCode
DeprecatedCoreInterface::set_cubit_msId_attribites_data_structure(
    const CubitBCType cubit_bc_type, const int ms_id,
    const GenericAttributeData &data, const std::string name) {
  return getInterface<MeshsetsManager>()->setAtributesByDataStructure(
      cubit_bc_type, ms_id, data, name);
}
MoFEMErrorCode DeprecatedCoreInterface::set_cubit_msId_bc_data_structure(
    const CubitBCType cubit_bc_type, const int ms_id,
    const GenericCubitBcData &data) {
  return getInterface<MeshsetsManager>()->setBcData(cubit_bc_type, ms_id, data);
}
MoFEMErrorCode
DeprecatedCoreInterface::delete_cubit_msId(const CubitBCType cubit_bc_type,
                                           const int ms_id) {
  return getInterface<MeshsetsManager>()->deleteMeshset(cubit_bc_type, ms_id);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId(
    const int ms_id, const CubitBCType cubit_bc_type,
    const CubitMeshSets **cubit_meshset_ptr) {
  return getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
      ms_id, cubit_bc_type, cubit_meshset_ptr);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int ms_id, const CubitBCType cubit_bc_type, const int dimension,
    Range &entities, const bool recursive) {
  return getInterface<MeshsetsManager>()->getEntitiesByDimension(
      ms_id, cubit_bc_type.to_ulong(), dimension, entities, recursive);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int ms_id, const CubitBCType cubit_bc_type, Range &entities,
    const bool recursive) {
  return getInterface<MeshsetsManager>()->getEntitiesByDimension(
      ms_id, cubit_bc_type.to_ulong(), entities, recursive);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int ms_id, const unsigned int cubit_bc_type, const int dimension,
    Range &entities, const bool recursive) {
  return getInterface<MeshsetsManager>()->getEntitiesByDimension(
      ms_id, cubit_bc_type, dimension, entities, recursive);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int ms_id, const unsigned int cubit_bc_type, Range &entities,
    const bool recursive) {
  return getInterface<MeshsetsManager>()->getEntitiesByDimension(
      ms_id, cubit_bc_type, entities, recursive);
}
MoFEMErrorCode DeprecatedCoreInterface::get_cubit_msId_meshset(
    const int ms_id, const unsigned int cubit_bc_type, EntityHandle &meshset) {
  return getInterface<MeshsetsManager>()->getMeshset(ms_id, cubit_bc_type,
                                                     meshset);
}

MoFEMErrorCode
DeprecatedCoreInterface::get_cubit_meshsets(const unsigned int cubit_bc_type,
                                            Range &meshsets) {
  return getInterface<MeshsetsManager>()->getMeshsetsByType(cubit_bc_type,
                                                            meshsets);
}

MoFEMErrorCode DeprecatedCoreInterface::build_problem_on_distributed_mesh(
    const std::string &name, const bool square_matrix, int verb) {
  return getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(
      name, square_matrix, verb);
}
MoFEMErrorCode DeprecatedCoreInterface::build_problem_on_distributed_mesh(
    Problem *problem_ptr, const bool square_matrix, int verb) {
  return getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(
      problem_ptr, square_matrix, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::build_problem(Problem *problem_ptr,
                                                      const bool square_matrix,
                                                      int verb) {
  return getInterface<ProblemsManager>()->buildProblem(problem_ptr,
                                                       square_matrix, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::build_problem(const std::string &problem_name,
                                       const bool square_matrix, int verb) {
  return getInterface<ProblemsManager>()->buildProblem(problem_name,
                                                       square_matrix, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::partition_simple_problem(const std::string &name,
                                                  int verb) {
  return getInterface<ProblemsManager>()->partitionSimpleProblem(name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::partition_compose_problem(
    const std::string &name, const std::string &problem_for_rows,
    bool copy_rows, const std::string &problem_for_cols, bool copy_cols,
    int verb) {
  return getInterface<ProblemsManager>()->inheritPartition(
      name, problem_for_rows, copy_rows, problem_for_cols, copy_cols, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::build_sub_problem(
    const std::string &out_name, const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col, const std::string &main_problem,
    const bool square_matrix, int verb) {
  return getInterface<ProblemsManager>()->buildSubProblem(
      out_name, fields_row, fields_col, main_problem, square_matrix, nullptr,
      nullptr, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::partition_problem(const std::string &name, int verb) {
  return getInterface<ProblemsManager>()->partitionProblem(name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::partition_finite_elements(
    const std::string &name, bool part_from_moab, int low_proc, int hi_proc,
    int verb) {
  return getInterface<ProblemsManager>()->partitionFiniteElements(
      name, part_from_moab, low_proc, hi_proc, verb);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
DeprecatedCoreInterface::partition_ghost_dofs(const std::string &name,
                                              int verb) {
  return getInterface<ProblemsManager>()->partitionGhostDofs(name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::get_problem_elements_layout(
    const std::string &name, const std::string &fe_name, PetscLayout *layout,
    int verb) {
  return getInterface<ProblemsManager>()->getProblemElementsLayout(
      name, fe_name, layout);
}

MoFEMErrorCode DeprecatedCoreInterface::partition_mesh(const Range &ents,
                                                       const int dim,
                                                       const int adj_dim,
                                                       const int n_parts,
                                                       int verb) {
  return getInterface<ProblemsManager>()->partitionMesh(
      ents, dim, adj_dim, n_parts, NULL, NULL, NULL, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::clear_dofs_fields(const BitRefLevel &bit,
                                           const BitRefLevel &mask, int verb) {
  return clear_dofs_fields_by_bit_ref(bit, mask, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::clear_ents_fields(const BitRefLevel &bit,
                                           const BitRefLevel &mask, int verb) {
  return clear_ents_fields_by_bit_ref(bit, mask, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::clear_finite_elements(
    const BitRefLevel &bit, const BitRefLevel &mask, int verb) {
  return clear_finite_elements_by_bit_ref(bit, mask, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_EDGEs(
    const Range &edges, const std::string &name, int verb) {
  return add_ents_to_field_by_type(edges, MBEDGE, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_EDGEs(
    const EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBEDGE, name, true, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_TRIs(
    const EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBTRI, name, true, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_TRIs(
    const Range &tris, const std::string &name, int verb) {
  return add_ents_to_field_by_type(tris, MBTRI, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_VERTICEs(
    const Range &nodes, const std::string &name, int verb) {
  return add_ents_to_field_by_type(nodes, MBVERTEX, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_VERTICEs(
    const EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBVERTEX, name, true, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_TETs(
    const Range &tets, const std::string &name, int verb) {
  return add_ents_to_field_by_type(tets, MBTET, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_TETs(
    const EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBTET, name, true, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_QUADs(
    const Range &quads, const std::string &name, int verb) {
  return add_ents_to_field_by_type(quads, MBQUAD, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_QUADs(
    EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBQUAD, name, true, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_PRISMs(
    const Range &prisms, const std::string &name, int verb) {
  return add_ents_to_field_by_type(prisms, MBPRISM, name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::add_ents_to_field_by_PRISMs(
    EntityHandle meshset, const std::string &name, int verb) {
  return add_ents_to_field_by_type(meshset, MBPRISM, name, true, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::MatCreateMPIAIJWithArrays(const std::string &name,
                                                   Mat *Aij, int verb) {
  return getInterface<MatrixManager>()
      ->createMPIAIJWithArrays<PetscGlobalIdx_mi_tag>(name, Aij, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::MatCreateMPIAdj_with_Idx_mi_tag(
    const std::string &name, Mat *Adj, int verb) {
  return getInterface<MatrixManager>()->createMPIAdjWithArrays<Idx_mi_tag>(
      name, Adj, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::MatCreateSeqAIJWithArrays(
    const std::string &name, Mat *Aij, PetscInt **i, PetscInt **j,
    PetscScalar **v, int verb) {
  MoFEMFunctionBegin;
  if (i || j || v)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "Not implemented");
  CHKERR getInterface<MatrixManager>()
      ->createSeqAIJWithArrays<PetscLocalIdx_mi_tag>(name, Aij, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DeprecatedCoreInterface::loop_finite_elements(
    const Problem *problem_ptr, const std::string &fe_name, FEMethod &method,
    int lower_rank, int upper_rank, MoFEMTypes bh, int verb) {
  return getInterface<CoreInterface>()->loop_finite_elements(
      problem_ptr, fe_name, method, lower_rank, upper_rank, nullptr, bh,
      CacheTupleSharedPtr(), verb);
}

MoFEMErrorCode DeprecatedCoreInterface::loop_finite_elements(
    const std::string &problem_name, const std::string &fe_name,
    FEMethod &method, int lower_rank, int upper_rank, MoFEMTypes bh, int verb) {
  return getInterface<CoreInterface>()->loop_finite_elements(
      problem_name, fe_name, method, lower_rank, upper_rank, nullptr, bh,
      CacheTupleSharedPtr(), verb);
}

MoFEMErrorCode DeprecatedCoreInterface::loop_finite_elements(
    const std::string &problem_name, const std::string &fe_name,
    FEMethod &method, MoFEMTypes bh, int verb) {
  return getInterface<CoreInterface>()->loop_finite_elements(
      problem_name, fe_name, method, get_comm_rank(), get_comm_rank(), nullptr,
      bh, CacheTupleSharedPtr(), verb);
}

MoFEMErrorCode DeprecatedCoreInterface::make_entities_multishared(
    const EntityHandle *entities, const int num_entities, const int my_proc,
    int verb) {
  return getInterface<CommInterface>()->makeEntitiesMultishared(
      entities, num_entities, my_proc, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::make_entities_multishared(
    Range &entities, const int my_proc, int verb) {
  return getInterface<CommInterface>()->makeEntitiesMultishared(entities,
                                                                my_proc, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::make_field_entities_multishared(
    const std::string field_name, const int owner_proc, int verb) {
  return getInterface<CommInterface>()->makeFieldEntitiesMultishared(
      field_name, owner_proc, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::exchange_field_data(const std::string field_name,
                                             int verb) {
  return getInterface<CommInterface>()->exchangeFieldData(field_name, verb);
}

MoFEMErrorCode DeprecatedCoreInterface::synchronise_entities(Range &ent,
                                                             int verb) {
  return getInterface<CommInterface>()->synchroniseEntities(ent, verb);
}

MoFEMErrorCode
DeprecatedCoreInterface::synchronise_field_entities(const std::string &name,
                                                    int verb) {
  return getInterface<CommInterface>()->synchroniseFieldEntities(name, verb);
}

DEPRECATED MoFEMErrorCode DeprecatedCoreInterface::remove_parents_by_by_bit_ref(
    const BitRefLevel bit, const BitRefLevel mask, int verb) {
  return remove_parents_by_bit_ref(bit, mask, verb);
}

} // namespace MoFEM
