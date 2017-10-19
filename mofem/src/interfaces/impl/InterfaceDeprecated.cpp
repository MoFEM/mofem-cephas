/** \file InterfaceDeprecated.cpp
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

  PetscErrorCode Interface_DEPRECATED_VERSION::seed_ref_level_2D(
    const EntityHandle meshset,const BitRefLevel &bit,int verb
  ) {
    return query_interface<BitRefManager>()->setBitRefLevelByDim(meshset,2,bit,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
    return query_interface<BitRefManager>()->setBitRefLevelByDim(meshset,3,bit,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::seed_ref_level(const Range &ents,const BitRefLevel &bit,const bool only_tets,int verb) {
    return query_interface<BitRefManager>()->setBitRefLevel(ents,bit,only_tets,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
    return query_interface<BitRefManager>()->setBitLevelToMeshset(meshset,bit,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::field_axpy(
    const double alpha,const std::string& field_name_x,const std::string& field_name_y,
    bool error_if_missing,bool creat_if_missing) {
    return query_interface<FieldBlas>()->fieldAxpy(alpha,field_name_x,field_name_y,error_if_missing,creat_if_missing);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_field(const double val,const EntityType type,const std::string& field_name) {
    return query_interface<FieldBlas>()->setField(val,type,field_name);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_field(const double val,const EntityType type,const Range &ents,const std::string& field_name) {
    return query_interface<FieldBlas>()->setField(val,type,ents,field_name);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::field_scale(const double alpha,const std::string& field_name) {
    return query_interface<FieldBlas>()->fieldScale(alpha,field_name);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return query_interface<BitRefManager>()->getAdjacenciesEquality(from_entiti,to_dimension,adj_entities);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) const {
    return query_interface<BitRefManager>()->getAdjacenciesAny(from_entiti,to_dimension,adj_entities);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_adjacencies(
    const Problem *problem_ptr,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return query_interface<BitRefManager>()->getAdjacencies(problem_ptr,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,
    const int num_netities,
    const int to_dimension,
    Range &adj_entities,
    const int operation_type,
    const int verb
  ) const {
    return query_interface<BitRefManager>()->getAdjacencies(bit,from_entities,num_netities,to_dimension,adj_entities,operation_type,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::VecCreateSeq(const std::string &name,RowColData rc,Vec *V) const {
    return query_interface<VecManager>()->vecCreateSeq(name,rc,V);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::VecCreateGhost(const std::string &name,RowColData rc,Vec *V) const {
    return query_interface<VecManager>()->vecCreateGhost(name,rc,V);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::ISCreateProblemOrder(
    const std::string &problem,RowColData rc,int min_order,int max_order,IS *is,int verb
  ) const {
    return query_interface<ISManager>()->isCreateProblemOrder(problem,rc,min_order,max_order,is);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::ISCreateProblemFieldAndRank(
    const std::string &problem,
    RowColData rc,
    const std::string &field,
    int min_coeff_idx,
    int max_coeff_idx,
    IS *is,
    int verb
  ) const {
    return query_interface<ISManager>()->isCreateProblemFieldAndRank(problem,rc,field,min_coeff_idx,max_coeff_idx,is);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    std::vector<int> &idx,std::vector<int> &idy,int verb
  ) const {
    return query_interface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,idx,idy
    );
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::ISCreateFromProblemFieldToOtherProblemField(
    const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    IS *ix,IS *iy,int verb
  ) const {
    return query_interface<ISManager>()->isCreateFromProblemFieldToOtherProblemField(
      x_problem,x_field_name,x_rc,y_problem,y_field_name,y_rc,ix,iy
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::VecScatterCreate(
    Vec xin,const std::string &x_problem,const std::string &x_field_name,RowColData x_rc,
    Vec yin,const std::string &y_problem,const std::string &y_field_name,RowColData y_rc,
    VecScatter *newctx,int verb
  ) const {
    return query_interface<VecManager>()->vecScatterCreate(
      xin,x_problem,x_field_name,x_rc,
      yin,y_problem,y_field_name,y_rc,
      newctx
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::VecScatterCreate(
    Vec xin,
    const std::string &x_problem,
    RowColData x_rc,
    Vec yin,
    const std::string &y_problem,
    RowColData y_rc,
    VecScatter *newctx,
    int verb
  ) const {
    return query_interface<VecManager>()->vecScatterCreate(
      xin,x_problem,x_rc,yin,y_problem,y_rc,newctx
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_local_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return query_interface<VecManager>()->setLocalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_local_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return query_interface<VecManager>()->setLocalGhostVector(name,rc,V,mode,scatter_mode);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_global_ghost_vector(
    const Problem *problem_ptr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return query_interface<VecManager>()->setGlobalGhostVector(problem_ptr,rc,V,mode,scatter_mode);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_global_ghost_vector(
    const std::string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode
  ) const {
    return query_interface<VecManager>()->setGlobalGhostVector(name,rc,V,mode,scatter_mode);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_other_local_ghost_vector(
    const Problem *problem_ptr,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
    Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
  ) {
    return query_interface<VecManager>()->setOtherLocalGhostVector(
      problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_other_local_ghost_vector(
    const std::string &name,const std::string& field_name,const std::string& cpy_field_name,RowColData rc,
    Vec V,InsertMode mode,ScatterMode scatter_mode,int verb
  ) {
    return query_interface<VecManager>()->setOtherLocalGhostVector(
      name,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_other_global_ghost_vector(
    const Problem *problem_ptr,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb
  ) {
    return query_interface<VecManager>()->setOtherGlobalGhostVector(
      problem_ptr,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_other_global_ghost_vector(
    const std::string &name,
    const std::string& field_name,
    const std::string& cpy_field_name,
    RowColData rc,
    Vec V,
    InsertMode mode,
    ScatterMode scatter_mode,
    int verb
  ) {
    return query_interface<VecManager>()->setOtherGlobalGhostVector(
      name,field_name,cpy_field_name,rc,V,mode,scatter_mode
    );
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::shift_right_bit_ref(const int shift,int verb) {
    return query_interface<BitRefManager>()->shiftRightBitRef(shift,BitRefLevel().set(),verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb) {
    return query_interface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit,mask,type,meshset,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb) {
    return query_interface<BitRefManager>()->getEntitiesByTypeAndRefLevel(bit,mask,type,ents,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) {
    return query_interface<BitRefManager>()->getEntitiesByRefLevel(bit,mask,meshset);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) {
    return query_interface<BitRefManager>()->getEntitiesByRefLevel(bit,mask,ents);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child,
    EntityType child_type,const bool recursive,int verb
  ) {
    return query_interface<BitRefManager>()->updateMeshsetByEntitiesChildren(parent,child_bit,child,child_type,recursive,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb) {
    return query_interface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(child_bit,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::update_field_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,int verb) {
    return query_interface<BitRefManager>()->updateFieldMeshsetByEntitiesChildren(name,child_bit,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::update_finite_element_meshset_by_entities_children(const std::string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb) {
    return query_interface<BitRefManager>()->updateFiniteElementMeshsetByEntitiesChildren(name,child_bit,fe_ent_type,verb);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::print_cubit_displacement_set() const {
    return query_interface_type<MeshsetsManager>()->printDisplacementSet(); 
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::print_cubit_pressure_set() const {
    return query_interface_type<MeshsetsManager>()->printPressureSet();
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::print_cubit_force_set() const {
    return query_interface_type<MeshsetsManager>()->printForceSet();
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::print_cubit_materials_set() const {
    return query_interface_type<MeshsetsManager>()->printMaterialsSet();
  }

  bool Interface_DEPRECATED_VERSION::check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type) {
    return query_interface<MeshsetsManager>()->checkMeshset(ms_id,cubit_bc_type);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
    return query_interface<MeshsetsManager>()->addMeshset(cubit_bc_type,ms_id,name);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::set_cubit_msId_attribites(
    const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name
  ) {
    return query_interface<MeshsetsManager>()->setAttribites(cubit_bc_type,ms_id,attributes,name);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_cubit_msId_attribites_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name
  ) {
    return query_interface<MeshsetsManager>()->setAttribitesByDataStructure(cubit_bc_type,ms_id,data,name);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::set_cubit_msId_bc_data_structure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
  ) {
    return query_interface<MeshsetsManager>()->setBcData(cubit_bc_type,ms_id,data);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id) {
    return query_interface<MeshsetsManager>()->deleteMeshset(cubit_bc_type,ms_id);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
    return query_interface<MeshsetsManager>()->getCubitMeshsetPtr(ms_id,cubit_bc_type,cubit_meshset_ptr);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId_entities_by_dimension(
    const int msId,const CubitBCType cubit_bc_type,const int dimension,Range &entities,const bool recursive
  ) {
    return query_interface<MeshsetsManager>()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),dimension,entities,recursive);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,Range &entities,const bool recursive) {
    return query_interface<MeshsetsManager>()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),entities,recursive);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId_entities_by_dimension(
    const int ms_id,const unsigned int cubit_bc_type,const int dimension,Range &entities,const bool recursive
  ) {
    MoFEMFunctionBeginHot;
    ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),dimension,entities,recursive); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }
  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive) {
    MoFEMFunctionBeginHot;
    ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),entities,recursive); CHKERRQ(ierr);
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_msId_meshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
    return query_interface<MeshsetsManager>()->getMeshset(ms_id,cubit_bc_type,meshset);
  }

  PetscErrorCode Interface_DEPRECATED_VERSION::get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) {
    return query_interface<MeshsetsManager>()->getMeshsetsByType(cubit_bc_type,meshsets);
  }


}
