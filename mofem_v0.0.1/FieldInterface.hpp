/** \file FieldInterface.hpp
 * \brief Myltindex containes, data structures and other low-level functions 
*/

/*! \mainpage Index Page
 * 
 * \section Introduction
 *
 * This is database for Finite Element System build on the MOAB
 * <https://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>.  It is tailored to keep
 * data for multiphisc problem with arbitray level of approximation and
 * diffrent levels of mesh refinments.  
 *
 * Finite Element Package for hp-adaptivity of MultiPhysics Problems for
 * H1,Hdiv,Hcurl and L2 approximation spaces.
 *
 * \section License
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted solely by Lukasz Kaczmarczyk.  It can be
 * freely used for educational and research purposes by other institutions. If
 * you use this software pleas cite my work. 
 *
 * MoFEM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 *
 * This is the introduction.
 *
 */

#ifndef __MOABFIELD_HPP__
#define __MOABFIELD_HPP__

#include "common.hpp"

namespace MoFEM {

/**
 * \brief The user interface
 * 
 * This class is used by user to: <br>
 *  (*) creat approximation fields,  <br>
 *  (*) define elements, <br>
 *  (*) define problems, <br>
 *  (*) manage refined meshses
 */
struct FieldInterface {
  /// get moab interface
  virtual Interface& get_moab() = 0; 

  /** 
    * \brief check data consistency in ents_moabfield
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_field(const string& name) = 0;

  /** 
    * \brief get entities form CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: from CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 - Volume(tetrahedral))
    * \param Range containing the retreived entities
    * \param recursive If true, meshsets containing meshsets are queried recursively. Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get all entities types from CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: from CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param entities form meshset
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get meshset from CUBIT Id and CUBIT type
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: from CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param meshset where to store the retreived entities
    */
  virtual PetscErrorCode get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) = 0;

  /** 
    * \brief get all CUBIT meshsets by CUBIT type
    *
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more). 
    * \param meshsets is range of meshsets
    */
  virtual PetscErrorCode get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets) = 0;

   /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
    * 	...
    * }
    *
    */
  virtual moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_begin() = 0;

   /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
    * 	...
    * }
    *
    */
  virtual moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_end() = 0;
    
    /**
     * \brief Iterator that loops over all the Cubit MeshSets in a moFEM field
     *
     * \param mField moFEM Field
     * \param iterator 
     */

  #define _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,IT) \
    moabCubitMeshSet_multiIndex::iterator IT = MFIELD.get_CubitMeshSets_begin(); IT!=MFIELD.get_CubitMeshSets_end(); IT++

  /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_end(const unsigned int CubitBCType) = 0;
    
    /**
    * \brief Iterator that loops over a specific Cubit MeshSet in a moFEM field
    *
    * \param mField moFEM Field
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param iterator 
    */

  #define _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_end(CUBITBCTYPE); IT++

  /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_end(const unsigned int CubitBCType) = 0;
    
    /**
     * \brief Iterator that loops over a specific Cubit MeshSet having a particular BC meshset in a moFEM field
     *
     * \param mField moFEM Field
     * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
     * \param iterator 
     *
     * Example: \code
       for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NodeSet|DisplacementSet,it) {
      	...
     * } \endcode
     */

  #define _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_bySetType_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_bySetType_end(CUBITBCTYPE); IT++


  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_begin(const string& name) = 0;
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_end(const string& name) = 0;

  #define _IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    moabCubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator IT = MFIELD.get_CubitMeshSets_byName_begin(NAME); \
    IT!=MFIELD.get_CubitMeshSets_byName_end(NAME); IT++

  virtual PetscErrorCode printCubitDisplacementSet() = 0;
  virtual PetscErrorCode printCubitPressureSet() = 0;
  virtual PetscErrorCode printCubitForceSet() = 0;
  virtual PetscErrorCode printCubitMaterials() = 0;

  /**
  * \brief seed 3D entities (Volume entities only) in the meshset and their adjacencies (only TETs adjencies) in a particular BitRefLevel
  * 
  * \param EntityHandle MeshSet
  * \param BitRefLevel bitLevel
  *
  * \brief Example:\code
  EntityHandle meshset1; //contains ent1,ent2,ent3
  BitRefLevel myLevel0;
  myLevel0.set(0);
  seed_ref_level_3D(meshset1,myLevel0); 
  //refine meshset1 into meshset2 and get new ents which are ent4, ent5
  EntityHandle meshset2; //contains ent1,ent2,ent3,ent4,ent5
  BitRefLevel myLevel1;
  myLevel1.set(1);
  seed_ref_level_3D(meshset2,myLevel1);  \endcode
   
  * So entities 1,2,3 would be assigned to bit level 0 and 1 <br>
  * ent1[1,1,0,0,0,0,0], ent2[1,1,0,0,0,0,0], ent3[1,1,0,0,0,0,0], <br> 
  * and entities 4 and 5 are assigned to bit level 1 only <br>
  * ent4[0,1,0,0,0,0,0], ent5[0,1,0,0,0,0,0] <br>
  * 
  */
  virtual PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /** brief seed ref level by MESHSET that contains entities other than volumes
   * 
   * \param EntityHandle MeshSet
   * \param BitRefLevel bitLevel
   */
  virtual PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit) = 0;

  /**
   * \brief make vertices in the middle of edges in meshset and add them to refinment levels defined by bit
   *
   * Takes entities fromm meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies.
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief make vertices in the middle of edges in meshset and add them to refinment levels defined by bit
   *
   * Takes entities fromm meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies.
   *
   * \param Range consisting edges for refine 
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(const Range &edges,const BitRefLevel &bit,int verb = -1) = 0;

  /**\brief refine TET in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = true) = 0;

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const Range &tets,const BitRefLevel &bit,const bool respect_interface = true) = 0;


  /**\brief refine PRISM in the meshset
   *
   * \param EntityHandle meshset
   * \param BitRefLevel bitLevel
   */

  virtual PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /**\brief refinem meshset, i.e. add child of refined entities to meshset
   *
   * \param EntityHandle meshset where to save the child refined entities
   * \param BitRefLevel bitLevel
   * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
   */
    
  virtual PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,
    const bool recursive = false,int verb = -1) = 0;

  /**\brief add FEs entities (Tets/Prisms) from ref level given by bit to meshset
   *
   * \param BitRefLevel bitLevel
   * \param EntityHandle meshset
   */
  virtual PetscErrorCode refine_get_finite_elements(const BitRefLevel &bit,const EntityHandle meshset) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityHandle meshset   
   *
   */
  virtual PetscErrorCode refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param Range   
   *
   *
   */
  virtual PetscErrorCode refine_get_ents(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) = 0;


  /** \brief Get childed entities form meshset containing parent entities 
    * 
    * Search for refined entities of given type whose parent are entities in the
    * parent meshset. It can be used for example to transfer information about
    * boundary conditions to refined mesh or splited mesh by interface
    * elements. It is used by function refine_MESHSET, to update MESHSET finite elements.
    * 
    * \param parent meshset
    * \param child_bit refinment level
    * \param type of refined entity
    * \param child_type meshset where child entities are stored (if the child meshset is set to be the parent meshset, the parent would be updated with the refined entities)
    * \param recursive if true parent meshset is searched recurively
    *
   **/
  virtual PetscErrorCode refine_get_childern(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1) = 0;

  /** 
    * \brief add approximation field
    *
    * \param name of the field
    * \param space approximation space (H1, Hdiv, Hcurl, L2 and NoField (dofs adjacent to meshset) 
    * \prama rank of the field, f.e. temperature has rank 1, displacement in 3d has rank 3
    */
  virtual PetscErrorCode add_field(const string& name,const FieldSpace space,const ApproximationRank rank,int verb = -1) = 0;

  /** 
    * \brief set field entities on vertices
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set vertices
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities form adjacencies of edges
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities form adjacencies of triangles
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities from adjacencies of tetrahedrals
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    *
    * \param meshset containing set of the entities (use 0 for all the entities in the meshset)
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order 
    */
  virtual PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order,int verb = -1) = 0;

  /// \brief list entities in the field
  virtual PetscErrorCode list_field() const = 0;

  /** \brief create field meshset
   *
   * \param name of Field
   * Example:\code
   EntityHandle disp_files_meshset = mField.get_field_meshset("DISPLACEMENT");
   * \endcode
   */
    
  virtual EntityHandle get_field_meshset(const string& name) const = 0;

  /**
    * \brief add finite element
    * \param name finite element name
    *
    * Example \code
      ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
      ierr = mField.add_finite_element("PLASTIC"); CHKERRQ(ierr);
   \endcode
    */
  virtual PetscErrorCode add_finite_element(const string &MoFEMFiniteElement_name) = 0;

  /** \brief set field data which finite element use
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual PetscErrorCode modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_filed) = 0;

    /** \brief set field row which finite element use
     *
     * \param name finite element name
     * \param name field name
     */
    
  virtual PetscErrorCode modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) = 0;

    /** \brief set field col which finite element use
     *
     * \param name finite element name
     * \param name field name
     */
    
    virtual PetscErrorCode modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_row) = 0;

  /** \brief add TET entities fromm meshset to finite element database given by name
   *
   * \param meshset contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive = false) = 0;

  /** \brief add TET elements from given refinment level to finite element database given by name 
   *
   * \param BitRefLevel BitLevel
   * \param name Finite Element name
   */

  virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type,int verb = -1) = 0;

  /** \brief add MESHSET element to finite element database given by name 
   *
   * \param meshset contains all entities that could be used for finite element
   * \param name Finite Element name
   */
    
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name) = 0;
    
    /** \brief add MESHSETs contained in meshset to finite element database given by name 
     *
     * \param meshset contains all meshsets with entities that could be used for finite element
     * \param name Finite Element name
     */
    
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSETs(const EntityHandle meshset,const string& name) = 0;

  /// list finite elements in database
  virtual PetscErrorCode list_finite_elements() const = 0;

  /// list adjacencies
  virtual PetscErrorCode list_adjacencies() const = 0;

  /// add Finite Element Problem
  virtual PetscErrorCode add_problem(const string& name) = 0;

  /* \brief add finite element to problem, this add entities assigned to finite element to a particular problem
   *
   * \param name Problem name
   * \param name Finite Element name
   */
    
  virtual PetscErrorCode modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name) = 0;

  /** \brief add ref level to problem
    *
   * if same finite element is solved using different level of refinements, than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC"); CHKERRQ(ierr);
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC"); CHKERRQ(ierr);
   
   ierr = mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1); CHKERRQ(ierr);
   ierr = mField.modify_problem_ref_level_add_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2); CHKERRQ(ierr);
   *\endcode
   * Two Problems exist and solved independently, both are elastic, but solved using different mesh refinement <br>
    */
    
  virtual PetscErrorCode modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) = 0;

  /// list problems
  virtual PetscErrorCode list_problem() const = 0;

  /** build fields
   */
  virtual PetscErrorCode build_fields(int verb = -1) = 0;

  /** build finite elements
    */
  virtual PetscErrorCode build_finite_elements(int verb = -1) = 0;

  /**
    * \brief build ajacencies 
    * \param bit adjacencies for refine level
    * This function will get information of adjacent finite elements and fields of all entities
    * If this is not perfomed, partitioning the problem is not possible
    */
  virtual PetscErrorCode build_adjacencies(const BitRefLevel bit) = 0;

  /** \brief build problem data structures
   */
    
  virtual PetscErrorCode build_problems(int verb = -1) = 0;

  /** \brief partition problem dofs
   *
   * \param name problem name
   */
  virtual PetscErrorCode simple_partition_problem(const string &name,int verb = -1) = 0;


  /** \brief partition problem dofs
   *
   * \param name problem name
   */
  virtual PetscErrorCode partition_problem(const string &name,int verb = -1) = 0;

  /**
    * \brief build indexing and partition problem inhereting indexing and partitioning from other problem
    *
    * \param name problem name
    * \param problem_for_rows problem used to index rows
    * \param problem_for_cols problem used to index cols
    */
  virtual PetscErrorCode compose_problem(const string &name,const string &problem_for_rows,const string &problem_for_cols,int var = -1) = 0;

  /**
    * \brief build indexing and partition problem inhereting indexing and partitioning from other problem
    *
    * \param name problem name
    * \param problem_for_rows problem used to index rows
    * \param copy_rows just copy rows dofs
    * \param problem_for_cols problem used to index cols
    * \param copy_cols just copy cols dofs
    *
    */
  virtual PetscErrorCode compose_problem(const string &name,const string &problem_for_rows,bool copy_rows,const string &problem_for_cols,bool copy_cols,int verb = -1) = 0;

  /** \brief determine ghost nodes
   *
   * \param name problem name 
   */
   
  virtual PetscErrorCode partition_ghost_dofs(const string &name,int verb = -1) = 0;

  /** \brief partition finite elements
   *
   * Function which partition finite elements based on dofs partitioning.<br>
   * In addition it sets information about local row and cols dofs at given element on partition. 
   *
   * \param name problem name 
   * \param do_skip if true, MultiIndices are set for finite element only on entities own by given partition
   */
  virtual PetscErrorCode partition_finite_elements(const string &name,bool do_skip = true,int verb = -1) = 0;

  /**
    * \brief add finite elements to the meshset
    *
    * Add finite elements to do meshset. 
    * \param name is problem name
    * \param fe_name
    * \param meshset
    */
  virtual PetscErrorCode problem_get_FE(const string &name,const string &fe_name,const EntityHandle meshset) = 0;

  /** \brief create ghost vector for problem
   *
   * \param name problem name
   * \param RowColData specify what data is taken from Row, Col or Data
   * \param Vec the vector where data is stored
   */
  virtual PetscErrorCode VecCreateGhost(const string &name,RowColData rc,Vec *V) = 0;

  /**
    * \brief create Mat (MPIAIJ) for problem
    *
    * \param name of the problem
    */
  virtual PetscErrorCode MatCreateMPIAIJWithArrays(const string &name,Mat *Aij,int verb = -1) = 0;

  /**
    * \brief create Mat (AIJ) for problem
    *
    * \param name of the problem
    */
  virtual PetscErrorCode MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,int verb = -1) = 0;

  /**
    * \brief create scatter for vectors form one to another problem
    *
    * \param xin vector
    * \param x_proble problem name
    * \param yin vector
    * \param y_problem problem name
    * \param newctx scatter
    */
  virtual PetscErrorCode VecScatterCreate(Vec xin,string &x_problem,RowColData x_rc,Vec yin,string &y_problem,RowColData y_rc,VecScatter *newctx,int verb = -1) = 0;

  /** 
    * \brief set values of vector from/to meshdatabase
    *
    * \param name of the problem
    * \param RowColData for row or column:e (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    * 
    * SCATTER_REVERSE set data to field entities from V vector.
    *
    * SCATTER_FORWARD set vector V from data field entities
    *
    */
  virtual PetscErrorCode set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;

  /** 
    * \brief set values of vector from/to meshdatabase
    *
    * \param name of the problem
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  virtual PetscErrorCode set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;

  /** \brief Copy vector to field which is not part of the problem
    *
    * \param name problem name
    * \param field_name field name used for indexing petsc vectors used in the problem
    * \param cpy_field field name where data from vector are stored
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  virtual PetscErrorCode set_other_global_VecCreateGhost(
    const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1) = 0;

  //topologity
  /** \brief create two children meshsets in the meshset containing terahedrals on two sides of faces
   *
   * \param msId Id of meshset 
   * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,
    const bool recursive = false,int verb = -1) = 0;

  /** \brief create two children meshsets in the meshset conta9ning terahedrals on two sides of faces
   *
   * Get tets adj to faces. Take skin form tets and get edges from that skin. Take skin form triangles (the face).
   * Subtrac skin faces edges form skin edges in order to get eges on the boundary of the face which is in the
   * voulume of the body, but is not on the boundary.
   * Each child set has a child contaning nodes which can be split and skin edges.
   * After that simply iterate under all tets on one side which are adjacent to the face are found.
   * Side tets are stored in to children meshsets of the SideSet meshset.
   */
  virtual PetscErrorCode get_msId_3dENTS_sides(const EntityHandle SideSet,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   * \param meshset meshset to get entities from
   * \param BitRefLevel new level where refinement would be stored
   * \param msId meshset ID imported from cubit 
   * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
   * \param add_intefece_entities meshset which contain the interface
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const Cubit_BC_bitset CubitBCType,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SideSet,const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;

  struct SnesMethod {
    enum snes_context { ctx_SNESSetFunction, ctx_SNESSetJacobian, ctx_SNESNone };
    //
    snes_context snes_ctx;
    SnesMethod(): snes_ctx(ctx_SNESNone) {};
    //
    PetscErrorCode set_snes_ctx(const snes_context ctx_);
    //
    SNES snes;
    PetscErrorCode set_snes(SNES _snes);
    Vec snes_x,snes_f;
    Mat *snes_A,*snes_B;
    MatStructure *snes_flag;
  };
  struct TSMethod {
    enum ts_context { ctx_TSSetRHSFunction, ctx_TSSetRHSJacobian, ctx_TSSetIFunction, ctx_TSSetIJacobian, ctx_TSTSMonitorSet, ctx_TSNone };
    //
    ts_context ts_ctx;
    TSMethod(): ts_ctx(ctx_TSNone) {};
    //
    PetscErrorCode set_ts_ctx(const ts_context ctx_);
    //
    TS ts;
    PetscErrorCode set_ts(TS _ts);
    Vec ts_u,ts_u_t,ts_F;
    Mat *ts_A,*ts_B;
    MatStructure *ts_flag;
    //
    PetscInt ts_step;
    PetscReal ts_a,ts_t;
  };

  struct BasicMethod: public SnesMethod,TSMethod {
    BasicMethod();    
    //
    PetscErrorCode set_problem(const MoFEMProblem *_problem_ptr);
    PetscErrorCode set_fields(const MoFEMField_multiIndex *_moabfields);
    PetscErrorCode set_ents_multiIndex(const MoFEMEntity_multiIndex *_ents_moabfield);
    PetscErrorCode set_dofs_multiIndex(const DofMoFEMEntity_multiIndex *_dofs_moabfield);
    PetscErrorCode set_fes_multiIndex(const MoFEMFiniteElement_multiIndex *_finite_elements);
    PetscErrorCode set_fes_data_multiIndex(const EntMoFEMFiniteElement_multiIndex *_finite_elements_moabents);
    PetscErrorCode set_adjacencies(const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *_fem_adjacencies);
    //
    virtual PetscErrorCode preProcess() = 0;
    virtual PetscErrorCode operator()() = 0;
    virtual PetscErrorCode postProcess() = 0;
    //
    const MoFEMProblem *problem_ptr;
    const MoFEMField_multiIndex *moabfields;
    const MoFEMEntity_multiIndex *ents_moabfield;
    const DofMoFEMEntity_multiIndex *dofs_moabfield;
    const MoFEMFiniteElement_multiIndex *finite_elements;
    const EntMoFEMFiniteElement_multiIndex *finite_elements_moabents;
    const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *fem_adjacencies;
  };

  /**
    * \brief structure for finite element method
    *
    * It can be used to calulate stiffnes matrices, residuals, load vectors etc.
    */  
  struct FEMethod: public BasicMethod {
    FEMethod();

    /** \brief function is run at the beginig of looop
     *
     * It is used to zeroing matrices and vectors, calculation of shape
     * functions on reference element, preporocessing boundary conditions, etc.
     */
    PetscErrorCode preProcess();

    /** \brief function is run for every finite element 
     *
     * It is used to calulate element local matrices and assembly. It can be
     * used for post-processing.
     */
    PetscErrorCode operator()();

    /** \brief function is run at the end of looop
     *
     * It is used to assembly matrices and vectors, calulating global variables,
     * f.e. total internal energy, ect.
     * 
     * Iterating over dofs:
     * Example1 iterating over dofs in row by name of the field
     * for(_IT_GET_FEROW_DOFS_FOR_LOOP_("DISPLACEMENT")) { ... } 
     * 
     * Example2 iterating over dofs in row by name of the field and type of the entity
     * for(_IT_GET_FEROW_DOFS_FOR_LOOP_("DISPLACEMENT",MBVERTEX)) { ... } 
     * 
     * Example2 iterating over dofs in row by name of the field, type of the entity and side number
     * for(_IT_GET_FEROW_DOFS_FOR_LOOP_("DISPLACEMENT",MBEDGE,0)) { ... }
     * 
     */
    PetscErrorCode postProcess();
    //
    PetscErrorCode set_fe(const NumeredMoFEMFiniteElement *_fe_ptr); 
    PetscErrorCode set_data_multIndex(const FEDofMoFEMEntity_multiIndex *_data_multiIndex);
    PetscErrorCode set_row_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_row_multiIndex);
    PetscErrorCode set_col_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_col_multiIndex);
    string fe_name;
    const NumeredMoFEMFiniteElement *fe_ptr;
    const FEDofMoFEMEntity_multiIndex *data_multiIndex;
    const FENumeredDofMoFEMEntity_multiIndex *row_multiIndex;
    const FENumeredDofMoFEMEntity_multiIndex *col_multiIndex;

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,  
      const string &field_name,const EntityType type,const int side_number) const {
      return index.lower_bound(boost::make_tuple(field_name,type,side_number));
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,  
      const string &field_name,const EntityType type,const int side_number) const {
      return index.upper_bound(boost::make_tuple(field_name,type,side_number));
    } 

      /** \brief loop over all dofs which are on a particular FE row, field, entity type and canonical side number
       *
       * \param FE finite elements
       * \param Name field name
       * \param Type moab entity type (MBVERTEX, MBEDGE etc)
       * \param Side side canonical number
       * \param IT the interator in use
       */
    #define _IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->row_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->row_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    ///loop over all dofs which are on a particular FE column, field, entity type and canonical side number

    #define _IT_GET_FECOL_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->col_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->col_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

      ///loop over all dofs which are on a particular FE data, field, entity type and canonical side number

    #define _IT_GET_FEDATA_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->data_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->data_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
      return index.lower_bound(boost::make_tuple(field_name,type));
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
      return index.upper_bound(boost::make_tuple(field_name,type));
    } 

      ///loop over all dofs which are on a particular FE row, field and entity type
    #define _IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->row_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->row_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); IT++
      ///loop over all dofs which are on a particular FE column, field and entity type
    #define _IT_GET_FECOL_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->col_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->col_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); IT++
      ///loop over all dofs which are on a particular FE data, field and entity type
    #define _IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->data_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->data_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name) const {
      return index.lower_bound(field_name);
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name) const {
      return index.upper_bound(field_name);
    } 
      ///loop over all dofs which are on a particular FE row and field
    #define _IT_GET_FEROW_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->row_multiIndex->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->row_multiIndex->get<FieldName_mi_tag>(),NAME); IT++
      ///loop over all dofs which are on a particular FE column and field
    #define _IT_GET_FECOL_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->col_multiIndex->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->col_multiIndex->get<FieldName_mi_tag>(),NAME); IT++
      ///loop over all dofs which are on a particular FE data and field
    #define _IT_GET_FEDATA_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->data_multiIndex->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->data_multiIndex->get<FieldName_mi_tag>(),NAME); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const EntityHandle ent) const {
      return index.lower_bound(ent);
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const EntityHandle ent) const {
      return index.upper_bound(ent);
    } 
      ///loop over all dofs which are on a particular FE row and given element entity (handle from moab)
                #define _IT_GET_FEROW_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->row_multiIndex->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->row_multiIndex->get<MoABEnt_mi_tag>(),ENT); IT++
      ///loop over all dofs which are on a particular FE column and given element entity (handle from moab)
    #define _IT_GET_FECOL_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->col_multiIndex->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->col_multiIndex->get<MoABEnt_mi_tag>(),ENT); IT++
      ///loop over all dofs which are on a particular FE data and given element entity (handle from moab)
    #define _IT_GET_FEDATA_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->data_multiIndex->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->data_multiIndex->get<MoABEnt_mi_tag>(),ENT); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
      return index.lower_bound(boost::make_tuple(field_name,ent));
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
      return index.upper_bound(boost::make_tuple(field_name,ent));
    } 
      ///loop over all dofs which are on a particular FE row, field and given element entity (handle from moab)
    #define _IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->row_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->row_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); IT++
      ///loop over all dofs which are on a particular FE column, field and given element entity (handle from moab)
    #define _IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->col_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->col_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); IT++
      ///loop over all dofs which are on a particular FE data, field and given element entity (handle from moab)
    #define _IT_GET_FEDATA_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->data_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->data_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); IT++

  };

  struct EntMethod: public BasicMethod {
    EntMethod();
    
    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();
    
    PetscErrorCode set_dof(const NumeredDofMoFEMEntity *_dof_ptr);
    const NumeredDofMoFEMEntity *dof_ptr;
  };

  /** \brief Set data for BasicMethod 
    *
    * This function set data about problem, adjacencies and other MultIindices
    * in database. This function can be used a special case when user need to
    * do some pre- and post-processing before matrix or vector is initiated, or
    * to assemble matrix for group of FEMethods. Is used by calsses classes
    * SnesCtx and TsCtx. Look for more details there.
    *
    * FIXME: Here we need example
    *
    * \param problem_name name of the problem
    * \param method user method derived from BasicMethod
    *
  **/
  virtual PetscErrorCode problem_basic_method_preProcess(const string &problem_name,BasicMethod &method,int verb = -1) = 0;

  /** \brief Set data for BasicMethod 
    *
    * This function set data about problem, adjacencies and other MultIindices
    * in database. This function can be used a special case when user need to
    * do some pre- and post-processing before matrix or vector is initiated, or
    * to assemble matrix for group of FEMethods. Is used by calsses classes
    * SnesCtx and TsCtx. Look for more details there.
    *
    * FIXME: Here we need example
    *
    * \param problem_name name of the problem
    * \param method user method derived from BasicMethod
    *
  **/
  virtual PetscErrorCode problem_basic_method_postProcess(const string &problem_name,BasicMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over finite elements. 
   *
   * This function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods exectuted on each element are given by
   * class derived form FieldInterface::FEMethod. At beginig of each loop user definded
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * Methods are executed only for local elements at given processor.
   *
   * For more details pleas look to examples.
   *
   * \param problem_name fe_name \param method is class derived form
   * FieldInterface::FEMethod
  **/ 
  virtual PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to lower rank. 
   *
   * This function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods exectuted on each element are given by
   * class derived form FieldInterface::FEMethod. At beginig of each loop user definded
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * For more details please look to examples.
   *
   * \param problem_name fe_name \param method is class derived form
   * FieldInterface::FEMethod
  **/ 
  virtual PetscErrorCode loop_finite_elements(
    const string &problem_name,const string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,int verb = -1) = 0;

  /** \brief Make a loop over entities
    *
    */
  virtual PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1) = 0;

  /** \brief Get problem database (datastructure) 
    *
    */
  virtual PetscErrorCode get_problem(const string &problem_name,const MoFEMProblem **problem_ptr) = 0;

  /** \brief Get dofs multi index
    *
    */
  virtual PetscErrorCode get_dofs(const DofMoFEMEntity_multiIndex **dofs_moabfield_ptr) = 0;

  /** 
    * \brief get begin iterator of filed ents of given name (instead you can use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_begin(const string &field_name) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_end(const string &field_name) = 0;

  ///loop over all dofs from a moFEM field and particular field
  #define _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MFIELD.get_ent_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_ent_moabfield_by_name_end(NAME); IT++

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_begin(const string &field_name) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_end(const string &field_name) = 0;

    ///loop over all dofs from a moFEM field and particular field
  #define _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MFIELD.get_dofs_by_name_begin(NAME); \
      IT != MFIELD.get_dofs_by_name_end(NAME); IT++

  /** \brief Get finite elements multi index
    *
    */
  virtual PetscErrorCode get_finite_elements(const MoFEMFiniteElement_multiIndex **finite_elements_ptr) = 0;

  /** 
    * \brief get begin iterator of finite elements of given name (instead you can use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_begin(const string &fe_name) = 0;

  /** 
    * \brief get end iterator of finite elements of given name (instead you can use _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_end(const string &fe_name) = 0;

    ///loop over all finite elements from a moFEM field and FE
  #define _IT_GET_FES_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator IT = MFIELD.get_fes_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_fes_moabfield_by_name_end(NAME); IT++


};

}

#endif // __MOABFIELD_HPP__
