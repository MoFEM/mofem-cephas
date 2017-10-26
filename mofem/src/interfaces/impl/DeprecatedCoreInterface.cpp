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

  PetscErrorCode DeprecatedCoreInterface::seed_ref_level_2D(
    const EntityHandle meshset,const BitRefLevel &bit,int verb
  ) {
    return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset,2,bit,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
    return getInterface<BitRefManager>()->setBitRefLevelByDim(meshset,3,bit,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::seed_ref_level(const Range &ents,const BitRefLevel &bit,const bool only_tets,int verb) {
    return getInterface<BitRefManager>()->setBitRefLevel(ents,bit,only_tets,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
    return getInterface<BitRefManager>()->setBitLevelToMeshset(meshset,bit,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::field_axpy(
    const double alpha,const std::string& field_name_x,const std::string& field_name_y,
    bool error_if_missing,bool creat_if_missing) {
    return getInterface<FieldBlas>()->fieldAxpy(alpha,field_name_x,field_name_y,error_if_missing,creat_if_missing);
  }
  PetscErrorCode DeprecatedCoreInterface::set_field(const double val,const EntityType type,const std::string& field_name) {
    return getInterface<FieldBlas>()->setField(val,type,field_name);
  }
  PetscErrorCode DeprecatedCoreInterface::set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name) {
    return getInterface<FieldBlas>()->setField(val,type,ents,field_name);
  }
  PetscErrorCode DeprecatedCoreInterface::field_scale(const double alpha,const std::string& field_name) {
    return getInterface<FieldBlas>()->fieldScale(alpha,field_name);
  }

  PetscErrorCode DeprecatedCoreInterface::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return getInterface<BitRefManager>()->getAdjacenciesEquality(from_entiti,to_dimension,adj_entities);
  }

  PetscErrorCode DeprecatedCoreInterface::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return getInterface<BitRefManager>()->getAdjacenciesAny(from_entiti,to_dimension,adj_entities);
  }

  PetscErrorCode DeprecatedCoreInterface::get_adjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return getInterface<BitRefManager>()->getAdjacencies(problem_ptr,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return getInterface<BitRefManager>()->getAdjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const {
    return getInterface<VecManager>()->vecCreateSeq(name,rc,V);
  }
  PetscErrorCode DeprecatedCoreInterface::VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const {
    return getInterface<VecManager>()->vecCreateGhost(name,rc,V);
  }
  PetscErrorCode DeprecatedCoreInterface::ISCreateProblemOrder(
    const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb
  ) const {
    return getInterface<ISManager>()->isCreateProblemOrder(problem,rc,min_order,max_order,is);
  }
  PetscErrorCode DeprecatedCoreInterface::ISCreateProblemFieldAndRank(
    const std::string &problem,
    RowColData rc,
    const std::string &field,
    int min_coeff_idx,
    int max_coeff_idx,
    IS *is,
    int verb
  ) const {
    return getInterface<ISManager>()->isCreateProblemFieldAndRank(problem,rc,field,min_coeff_idx,max_coeff_idx,is);
  }
  PetscErrorCode DeprecatedCoreInterface::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy,int verb
  ) const {
    return getInterface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,idx,idy
    );
  }

  PetscErrorCode DeprecatedCoreInterface::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy,int verb
  ) const {
    return getInterface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,ix,iy
    );
  }
  PetscErrorCode DeprecatedCoreInterface::VecScatterCreate(
    Vec xin,const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    Vec yin,const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    VecScatter *newctx,int verb
  ) const {
    return getInterface<VecManager>()->vecScatterCreate(
      xin,x_problem,x_field_name,x_rc,
      yin,y_problem,y_field_name,y_rc,
      newctx
    );
  }
  PetscErrorCode DeprecatedCoreInterface::VecScatterCreate(
    Vec xin,
    const std::string &x_problem,
    RowColData x_rc,
    Vec yin,
    const std::string &y_problem,
    RowColData y_rc,
    VecScatter *newctx,
    int verb
  ) const {
    return getInterface<VecManager>()->vecScatterCreate(
      xin,x_problem,x_rc,yin,y_problem,y_rc,newctx
    );
  }
  PetscErrorCode DeprecatedCoreInterface::set_local_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return getInterface<VecManager>()->setLocalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
  }
  PetscErrorCode DeprecatedCoreInterface::set_local_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return getInterface<VecManager>()->setLocalGhostVector(name,rc,V,mode,scatter_mode);
  }
  PetscErrorCode DeprecatedCoreInterface::set_global_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return getInterface<VecManager>()->setGlobalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
  }
  PetscErrorCode DeprecatedCoreInterface::set_global_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return getInterface<VecManager>()->setGlobalGhostVector(name,rc,V,mode,scatter_mode);
  }
  PetscErrorCode DeprecatedCoreInterface::set_other_local_ghost_vector(
    const Problem *problem_ptr,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
    Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
  ) {
    return getInterface<VecManager>()->setOtherLocalGhostVector(
      problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode DeprecatedCoreInterface::set_other_local_ghost_vector(
    const std::string &name,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
    Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
  ) {
    return getInterface<VecManager>()->setOtherLocalGhostVector(
      name,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode DeprecatedCoreInterface::set_other_global_ghost_vector(
    const Problem *problem_ptr,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb
  ) {
    return getInterface<VecManager>()->setOtherGlobalGhostVector(
      problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode DeprecatedCoreInterface::set_other_global_ghost_vector(
    const std::string &name,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb
  ) {
    return getInterface<VecManager>()->setOtherGlobalGhostVector(
      name,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }

  PetscErrorCode DeprecatedCoreInterface::shift_right_bit_ref(const int shift,int verb) {
    return getInterface<BitRefManager>()->shiftRightBitRef(shift,BitRefLevel().set(),verb);
  }

  PetscErrorCode DeprecatedCoreInterface::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
    return getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit,mask,type,meshset,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
    return getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit,mask,type,ents,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
    return getInterface<BitRefManager>()->getEntitiesByRefLevel(bit,mask,meshset);
  }

  PetscErrorCode DeprecatedCoreInterface::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) {
    return getInterface<BitRefManager>()->getEntitiesByRefLevel(bit,mask,ents);
  }

  PetscErrorCode DeprecatedCoreInterface::update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child,
    EntityType child_type,const bool recursive,int verb
  ) {
    return getInterface<BitRefManager>()->updateMeshsetByEntitiesChildren(parent,child_bit,child,child_type,recursive,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb) {
    return getInterface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(child_bit,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::update_field_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,int verb) {
    return getInterface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(name,child_bit,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::update_finite_element_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb) {
    return getInterface<BitRefManager>()->updateFiniteElementMeshsetByEntitiesChildren(name,child_bit,fe_ent_type,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::print_cubit_displacement_set() const {
    return getInterface<MeshsetsManager>()->printDisplacementSet();
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode DeprecatedCoreInterface::print_cubit_pressure_set() const {
    return getInterface<MeshsetsManager>()->printPressureSet();
  }

  PetscErrorCode DeprecatedCoreInterface::print_cubit_force_set() const {
    return getInterface<MeshsetsManager>()->printForceSet();
  }

  PetscErrorCode DeprecatedCoreInterface::print_cubit_materials_set() const {
    return getInterface<MeshsetsManager>()->printMaterialsSet();
  }

  bool DeprecatedCoreInterface::check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type) {
    return getInterface<MeshsetsManager>()->checkMeshset(ms_id,cubit_bc_type);
  }

  PetscErrorCode DeprecatedCoreInterface::add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
    return getInterface<MeshsetsManager>()->addMeshset(cubit_bc_type,ms_id,name);
  }

  PetscErrorCode DeprecatedCoreInterface::set_cubit_msId_attribites(
    const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name
  ) {
    return getInterface<MeshsetsManager>()->setAtributes(cubit_bc_type,ms_id,attributes,name);
  }
  PetscErrorCode DeprecatedCoreInterface::set_cubit_msId_attribites_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name
  ) {
    return getInterface<MeshsetsManager>()->setAtributesByDataStructure(cubit_bc_type,ms_id,data,name);
  }
  PetscErrorCode DeprecatedCoreInterface::set_cubit_msId_bc_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
  ) {
    return getInterface<MeshsetsManager>()->setBcData(cubit_bc_type,ms_id,data);
  }
  PetscErrorCode DeprecatedCoreInterface::delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id) {
    return getInterface<MeshsetsManager>()->deleteMeshset(cubit_bc_type,ms_id);
  }
  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
    return getInterface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,cubit_bc_type,cubit_meshset_ptr);
  }
  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int msId,const CubitBCType cubit_bc_type,const int dimension,Range &entities,const bool recursive
  ) {
    return getInterface<MeshsetsManager>()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),dimension,entities,recursive);
  }
  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,Range &entities,const bool recursive) {
    return getInterface<MeshsetsManager>()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),entities,recursive);
  }
  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(
    const int ms_id,const unsigned int cubit_bc_type,const int dimension,Range &entities,const bool recursive
  ) {
    MoFEMFunctionBeginHot;
    ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),dimension,entities,recursive); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }
  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive) {
    MoFEMFunctionBeginHot;
    ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),entities,recursive); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode DeprecatedCoreInterface::get_cubit_msId_meshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
    return getInterface<MeshsetsManager>()->getMeshset(ms_id,cubit_bc_type,meshset);
  }

  PetscErrorCode DeprecatedCoreInterface::get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) {
    return getInterface<MeshsetsManager>()->getMeshsetsByType(cubit_bc_type,meshsets);
  }

  PetscErrorCode DeprecatedCoreInterface::build_problem_on_distributed_mesh(
    const std::string &name,const bool square_matrix,int verb
  ) {
    return getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(name,square_matrix,verb);
  }
  PetscErrorCode DeprecatedCoreInterface::build_problem_on_distributed_mesh(
    Problem *problem_ptr,const bool square_matrix,int verb
  ) {
    return getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(problem_ptr,square_matrix,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::build_problem(Problem *problem_ptr,const bool square_matrix,int verb) {
    return getInterface<ProblemsManager>()->buildProblem(problem_ptr,square_matrix,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::build_problem(const std::string &problem_name,const bool square_matrix,int verb) {
    return getInterface<ProblemsManager>()->buildProblem(problem_name,square_matrix,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::partition_simple_problem(const std::string &name,int verb) {
    return getInterface<ProblemsManager>()->partitionSimpleProblem(name,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::partition_compose_problem(
    const std::string &name,
    const std::string &problem_for_rows,
    bool copy_rows,
    const std::string &problem_for_cols,
    bool copy_cols,
    int verb
  ) {
    return getInterface<ProblemsManager>()->inheretPartition(
      name,problem_for_rows,copy_rows,problem_for_cols,copy_cols,verb
    );
  }

  PetscErrorCode DeprecatedCoreInterface::build_sub_problem(
    const std::string &out_name,
    const std::vector<std::string> &fields_row,
    const std::vector<std::string> &fields_col,
    const std::string &main_problem,
    const bool square_matrix,
    int verb
  ) {
    return getInterface<ProblemsManager>()->buildSubProblem(
      out_name,fields_row,fields_col,main_problem,square_matrix,verb
    );
  }

  PetscErrorCode DeprecatedCoreInterface::partition_problem(const std::string &name,int verb) {
    return getInterface<ProblemsManager>()->partitionProblem(name,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::partition_finite_elements(
    const std::string &name,
    bool part_from_moab,
    int low_proc,
    int hi_proc,
    int verb
  ) {
    return getInterface<ProblemsManager>()->partitionFiniteElements(
      name,part_from_moab,low_proc,hi_proc,verb
    );
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode DeprecatedCoreInterface::partition_ghost_dofs(const std::string &name,int verb) {
    return getInterface<ProblemsManager>()->partitionGhostDofs(name,verb);
  }

  PetscErrorCode DeprecatedCoreInterface::get_problem_elements_layout(const std::string &name,const std::string &fe_name,PetscLayout *layout,int verb) {
    return getInterface<ProblemsManager>()->getProblemElementsLayout(name,fe_name,layout);
  }

  PetscErrorCode DeprecatedCoreInterface::partition_mesh(
    const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb
  ) {
    return getInterface<ProblemsManager>()->partitionMesh(ents,dim,adj_dim,n_parts,NULL,NULL,NULL,verb);
  }



}
