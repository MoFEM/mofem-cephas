/** \file moabField.hpp
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
 *  (*) define elenets, <br>
 *  (*) define problems, <br>
 *  (*) manage refined meshses
 */
struct moabField {
  /// get moab interface
  virtual Interface& get_moab() = 0; 

  /** 
    * \brief check data consistency in ents_moabfield
    *
    */
  virtual PetscErrorCode check_NumbetOfEnts_in_ents_moabfield(const string& name) = 0;

  /** 
    * \brief get entities form CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param entities form meshset
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get entities form CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param entities form meshset
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get meshset form CUBIT Id and CUBIT type
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more) 
    * \param meshset 
    */
  virtual PetscErrorCode get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) = 0;

  /** 
    * \brief get all CUBIT meshset by CUBIT type
    *
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more). 
    * \param meshsets is range of meshsets
    */
  virtual PetscErrorCode get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets) = 0;

   /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BYFOR_LOOP_(mFiled,it) {
    * 	...
    * }
    *
    */
  virtual moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_begin() = 0;

   /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BYFOR_LOOP_(mFiled,it) {
    * 	...
    * }
    *
    */
  virtual moabCubitMeshSet_multiIndex::iterator get_CubitMeshSets_end() = 0;

  #define _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,IT) \
    moabCubitMeshSet_multiIndex::iterator IT = MFIELD.get_CubitMeshSets_begin(); IT!=MFIELD.get_CubitMeshSets_end(); IT++

  /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mFiled,NodeSet|DisplacementSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mFiled,NodeSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_end(const unsigned int CubitBCType) = 0;

  #define _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    moabCubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_end(CUBITBCTYPE); IT++

  /** 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mFiled,NodeSet|DisplacementSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mFiled,NodeSet,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NodeSet, SideSet or BlockSet and more)
    */
  virtual moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_end(const unsigned int CubitBCType) = 0;

  #define _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    moabCubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_bySetType_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_bySetType_end(CUBITBCTYPE); IT++


  /// seed ref level by 3D ents in the meshset and its adjacencies (only TETs adjencies)
  virtual PetscErrorCode seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit) = 0;

  /// seed ref level by MESHSET
  virtual PetscErrorCode seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit) = 0;

  /**
   * \brief make vetices in the middle of edges in meshset and add them to refinment levels defined by bit
   *
   * Takes entities form meshsets and queried recursively (get entities from meshsets in meshsets, usually have to be used for CUBIT meshset).
   * If meshset does not contain any edges, get entities in dimension 3 and get edge adjacencies. 
   */
  virtual PetscErrorCode add_verices_in_the_middel_of_edges(
    const EntityHandle meshset,const BitRefLevel &bit,const bool recursive = false,int verb = -1) = 0;

  /// refine TET in the meshset
  virtual PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = true) = 0;

  /// refine PRISM in the meshset
  virtual PetscErrorCode refine_PRISM(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /// refinem meshset, i.e. add child of refined entities to meshset
  virtual PetscErrorCode refine_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,
    const bool recursive = false,int verb = -1) = 0;

  /// add FEs form ref level given by bit to meshset
  virtual PetscErrorCode refine_get_finite_elements(const BitRefLevel &bit,const EntityHandle meshset) = 0;

  /// add ents form ref level given by bit to meshset
  virtual PetscErrorCode refine_get_ents(const BitRefLevel &bit,const EntityHandle meshset) = 0;

  /** \brief Get childed entities form meshste containing parent entities 
    * 
    * Serch for refained entities of given type whose paren are entities in the
    * parent meshset. It can be used for example to transfer information about
    * boundary conditions to refined mesh or splited mesh by interface
    * elements. It is used by function refine_MESHSET, to update MESHSET finite elements.
    * 
    * \param paren meshset
    * \param child_bit refinment level
    * \param type of refined entity
    * \param child_type meshset where child entities are stored
    * \param recursive if true parent meshset is searched recurively
    *
   **/
  virtual PetscErrorCode refine_get_childern(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1) = 0;

  /** 
    * \brief add appeoximation field
    *
    * \param name of the field
    * \param space approximation space (H1, Hdiv, Hcurl, L2 and NoField (dofs adjacent to meshset) 
    * \prama rank of the field, f.e. temeraure has rank 1, displacement in 3d has rank 3
    */
  virtual PetscErrorCode add_field(const string& name,const FieldSpace space,const ApproximationRank rank,int verb = -1) = 0;

  /** 
    * \brief set field entities on veryices
    *
    * The lower dimension entities are added depending on the space type
    * \param meshet contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities form adjacencies of triangles
    *
    * The lower dimension entities are added depending on the space type
    * \param meshet contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities from adjacencies of tetrahedrals
    *
    * The lower dimension entities are added depending on the space type
    * \param meshet contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /**
    * \brief Set order of the entities in the field
    *
    * \param meshset containg set of the entities
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order 
    */
  virtual PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order) = 0;

  /// \brief list entities in the field
  virtual PetscErrorCode list_field() const = 0;

  /// \brief get field meshset
  virtual EntityHandle get_field_meshset(const string& name) const = 0;

  /**
    * \brief add finite element
    * \param name finite elenent name
    */
  virtual PetscErrorCode add_finite_element(const string &MoFEMFiniteElement_name) = 0;

  /// \brief set field data which element use
  virtual PetscErrorCode modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_filed) = 0;

  /// \brief set field row
  virtual PetscErrorCode modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) = 0;

  /// \brief set field col
  virtual PetscErrorCode modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_row) = 0;

  /// add TET elements form meshset to finite element database given by name 
  virtual PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive = false) = 0;

  /// add TET elements from given refinment level to finite element database given by name 
  virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type) = 0;

  /// add MESHSET element to finite element database given by name 
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name) = 0;

  /// add MESHSETs contained in meshset to finite element database given by name 
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSETs(const EntityHandle meshset,const string& name) = 0;

  /// list finite elements in database
  virtual PetscErrorCode list_finite_elements() const = 0;

  /// list adjacencies
  virtual PetscErrorCode list_adjacencies() const = 0;

  /// add problem
  virtual PetscErrorCode add_problem(const string& name) = 0;

  /// \brief add finite element to problem
  virtual PetscErrorCode modify_problem_add_finite_element(const string &name_problem,const string &MoFEMFiniteElement_name) = 0;

  /// \brief add ref level to problem
  virtual PetscErrorCode modify_problem_ref_level_add_bit(const string &name_problem,const BitRefLevel &bit) = 0;

  /// list problems
  virtual PetscErrorCode list_problem() const = 0;

  /// build fields
  virtual PetscErrorCode build_fields(int verb = -1) = 0;

  /// build finite elements
  virtual PetscErrorCode build_finite_elements(int verb = -1) = 0;

  /**
    * \brief build ajacencies 
    * \param bit adjacencies for refine level
    */
  virtual PetscErrorCode build_adjacencies(const BitRefLevel bit) = 0;

  /// \brief build proble data structures
  virtual PetscErrorCode build_problems(int verb = -1) = 0;

  /// partition problem dofs
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

  /// determine ghost nodes
  virtual PetscErrorCode partition_ghost_dofs(const string &name,int verb = -1) = 0;

  /// partition finite elements
  virtual PetscErrorCode partition_finite_elements(const string &name,bool do_skip = true,int verb = -1) = 0;

  /// erase inactive dofs form field
  virtual PetscErrorCode erase_inactive_dofs_moabfield() = 0;

  /**
    * \brief add finite elements to the meshset
    *
    * Add finite elements to de meshset. 
    * \param name is problem name
    * \param fe_name
    * \param meshset
    */
  virtual PetscErrorCode problem_get_FE(const string &name,const string &fe_name,const EntityHandle meshset) = 0;

  /// create ghost vector for problem
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
    * \brief set values of vector form/to meshdatabase
    *
    * \param name of the problem
    * \param RowColData for row or column:e (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    * 
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    * SCATTER_FORWARD set vector V from data field entities
    *
    */
  virtual PetscErrorCode set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;

  /** 
    * \brief set values of vector form/to meshdatabase
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



  //topoltgy
  /// create two cheldren meshsets in the meshset contaning terahedrals on two sides of faces
  virtual PetscErrorCode get_msId_3dENTS_sides(const int msId,const Cubit_BC_bitset CubitBCType,
    const bool recursive = false,int verb = -1) = 0;

  /** 
   * \brief create two cheldren meshsets in the meshset contaning terahedrals on two sides of faces
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
   * \brif split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const Cubit_BC_bitset CubitBCType,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brif split nodes and other entities of tetrahedrals in children sets and add prism elements
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
    PetscErrorCode set_moabfields(const MoFEMField_multiIndex *_moabfields);
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
     * It is used to zeroing matrices and vectors, clalualtion of shape
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
     * It is used to assembly matrices and vectors, calulating global varibles,
     * f.e. total internal energy, ect.
     * 
     * Iterating over dofs:
     * Example1 iteratong over dofs in row by name of the field
     * for(_IT_GET_FEROW_DOFS_FOR_LOOP_("DISPLACEMENT")) { ... } 
     * 
     * Example2 iteratong over dofs in row by name of the field and type of the entity
     * for(_IT_GET_FEROW_DOFS_FOR_LOOP_("DISPLACEMENT",MBVERTEX)) { ... } 
     * 
     * Example2 iteratong over dofs in row by name of the field, type of the entity and side number
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

    #define _IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->row_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->row_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    #define _IT_GET_FECOL_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->col_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->col_multiIndex->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

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

    #define _IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->row_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->row_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); IT++
    #define _IT_GET_FECOL_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->col_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag2>::type>(FE->col_multiIndex->get<Composite_mi_tag2>(),NAME,TYPE); IT++
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

    #define _IT_GET_FEROW_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->row_multiIndex->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->row_multiIndex->get<FieldName_mi_tag>(),NAME); IT++
    #define _IT_GET_FECOL_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->col_multiIndex->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->col_multiIndex->get<FieldName_mi_tag>(),NAME); IT++
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

    #define _IT_GET_FEROW_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->row_multiIndex->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->row_multiIndex->get<MoABEnt_mi_tag>(),ENT); IT++
    #define _IT_GET_FECOL_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->col_multiIndex->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->col_multiIndex->get<MoABEnt_mi_tag>(),ENT); IT++
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

    #define _IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->row_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->row_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); IT++
    #define _IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->col_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag3>::type>(FE->col_multiIndex->get<Composite_mi_tag3>(),NAME,ENT); IT++
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
    * moabSnes and moabTs. Look for more details there.
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
    * moabSnes and moabTs. Look for more details there.
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
   * Thsis function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods exectuted on each element are given by
   * class derived form moabField::FEMethod. At beginig of each loop user definded
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * Methods are executed only for local elements at given processor.
   *
   * For more details pleas look to examples.
   *
   * \param problem_name fe_name \param method is class derived form
   * moabField::FEMethod
  **/ 
  virtual PetscErrorCode loop_finite_elements(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over finite elements on partitions from upper to lower rank. 
   *
   * Thsis function is like swiss knife, is can be used to post-processing or matrix
   * and vectors assembly. It makes loop over given finite element for given
   * problem. The particular methods exectuted on each element are given by
   * class derived form moabField::FEMethod. At beginig of each loop user definded
   * function (method)  preProcess() is called, for each element operator() is
   * executed, at the end loop finalizes with user defined function (method)
   * postProcess().
   *
   * For more details pleas look to examples.
   *
   * \param problem_name fe_name \param method is class derived form
   * moabField::FEMethod
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
  virtual PetscErrorCode get_problems_database(const string &problem_name,const MoFEMProblem **problem_ptr) = 0;

  /** \brief Get dofs multi index
    *
    */
  virtual PetscErrorCode get_dofs_moabfield(const DofMoFEMEntity_multiIndex **dofs_moabfield_ptr) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_moabfield_by_name_begin(const string &field_name) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_moabfield_by_name_end(const string &field_name) = 0;

  #define _IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MFIELD.get_dofs_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_dofs_moabfield_by_name_end(NAME); IT++

  /** \brief Get finite elements multi index
    *
    */
  virtual PetscErrorCode get_finite_elements(const MoFEMFiniteElement_multiIndex **finite_elements_ptr) = 0;

  /** 
    * \brief get begin iterator of finite elements of given name (instead you can use _IT_GET_FES_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_FES_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_begin(const string &fe_name) = 0;

  /** 
    * \brief get end iterator of finite elements of given name (instead you can use _IT_GET_FES_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    *
    * for(_IT_GET_FES_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param fe_name
    */
  virtual EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator get_fes_moabfield_by_name_end(const string &fe_name) = 0;

  #define _IT_GET_FES_MOABFIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    EntMoFEMFiniteElement_multiIndex::index<MoFEMFiniteElement_name_mi_tag>::type::iterator IT = MFIELD.get_fes_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_fes_moabfield_by_name_end(NAME); IT++


};

}

#endif // __MOABFIELD_HPP__
