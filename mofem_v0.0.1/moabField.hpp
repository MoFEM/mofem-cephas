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
    * \brief get entities form CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more. 
    * \param entities form meshset
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get entities form CUBIT/meshset 
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more. 
    * \param entities form meshset
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false) = 0;

  /** 
    * \brief get meshset form CUBIT Id and CUBIT type
    *
    * \param msId id of the BlockSet/SideSet/BlockSet: form CUBIT
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more. 
    * \param meshsset 
    */
  virtual PetscErrorCode get_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) = 0;

  /** 
    * \brief get all CUBIT meshset by CUBIT type
    *
    * \param CubitBCType see Cubit_BC (NodeSet, SideSet or BlockSet and more. 
    * \param meshsets is range of meshsets
    */
  virtual PetscErrorCode get_CubitBCType_meshsets(const unsigned int CubitBCType,Range &meshsets) = 0;

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
  virtual PetscErrorCode refine_get_FE(const BitRefLevel &bit,const EntityHandle meshset) = 0;

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
    * \param space approximation space
    * \prama rank of the field, f.e. temeraure has rank 1, displacement in 3d has rank 3
    */
  virtual PetscErrorCode add_BitFieldId(const string& name,const FieldSpace space,const ApproximationRank rank,int verb = -1) = 0;

  /** 
    * \brief set tetrahedrals part of the given field 
    *
    * The lower dimension entities are added depending on the space type
    * \param meshet contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name) = 0;

  /** 
    * \brief set prisms part of the given field (works only for L2 space)
    *
    * The lower dimension entities are added depending on the space type (works only for L2 space)
    * \param meshet contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_PRISMs(const EntityHandle meshset,const string& name) = 0;

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
  virtual EntityHandle get_meshset_by_BitFieldId(const string& name) const = 0;

  /**
    * \brief add finite element
    * \param name finite elenent name
    */
  virtual PetscErrorCode add_MoFEMFE(const string &MoFEMFE_name) = 0;

  /// \brief set field data which element use
  virtual PetscErrorCode modify_MoFEMFE_data_add_bit(const string &MoFEMFE_name,const string &name_filed) = 0;

  /// \brief set field row
  virtual PetscErrorCode modify_MoFEMFE_row_add_bit(const string &MoFEMFE_name,const string &name_row) = 0;

  /// \brief set field col
  virtual PetscErrorCode modify_MoFEMFE_col_add_bit(const string &MoFEMFE_name,const string &name_row) = 0;

  /// add TET elements form meshset to finite element database given by name 
  virtual PetscErrorCode add_ents_to_MoFEMFE_by_TETs(const EntityHandle meshset,const string &name) = 0;

  /// add TET elements to the refinment level to finite element database given by name 
  virtual PetscErrorCode add_ents_to_MoFEMFE_EntType_by_bit_ref(const BitRefLevel &bit_ref,const string &name,EntityType type) = 0;

  /// add MESHSET element to finite element database given by name 
  virtual PetscErrorCode add_ents_to_MoFEMFE_by_MESHSET(const EntityHandle meshset,const string& name) = 0;

  /// add MESHSETs contained in meshset to finite element database given by name 
  virtual PetscErrorCode add_ents_to_MoFEMFE_by_MESHSETs(const EntityHandle meshset,const string& name) = 0;

  /// list finite elements in database
  virtual PetscErrorCode list_MoFEMFE() const = 0;

  /// list adjacencies
  virtual PetscErrorCode list_adjacencies() const = 0;

  /// add problem
  virtual PetscErrorCode add_BitProblemId(const string& name) = 0;

  /// \brief add finite element to problem
  virtual PetscErrorCode modify_problem_MoFEMFE_add_bit(const string &name_problem,const string &MoFEMFE_name) = 0;

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
  virtual PetscErrorCode partition_problems(const string &name,int verb = -1) = 0;

  /// determine ghost nodes
  virtual PetscErrorCode partition_ghost_dofs(const string &name) = 0;

  /// partition finite elements
  virtual PetscErrorCode partition_finite_elements(const string &name,int verb = -1) = 0;

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
    * \brief set values of vector form/to meshdatabase
    *
    * \param name of the problem
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    */
  virtual PetscErrorCode set_local_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;

  /** 
    * \brief set values of vector form/to meshdatabase
    *
    * \param name of the problem
    * \param RowColData for row or column
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    */
  virtual PetscErrorCode set_global_VecCreateGhost(const string &name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;


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
    SNES snes;
    PetscErrorCode set_snes(SNES _snes) { 
      PetscFunctionBegin;
      snes = _snes;
      PetscFunctionReturn(0);
    }
    Vec snes_x,snes_f;
    PetscErrorCode set_x(Vec _x) {
      PetscFunctionBegin;
      snes_x = _x;
      PetscFunctionReturn(0);
    }
    PetscErrorCode set_f(Vec _f) {
      PetscFunctionBegin;
      snes_f = _f;
      PetscFunctionReturn(0);
    }
    Mat *snes_A,*snes_B;
    MatStructure *snes_flag;
    PetscErrorCode set_A(Mat *_A) {
      PetscFunctionBegin;
      snes_A = _A;
      PetscFunctionReturn(0);
    }
    PetscErrorCode set_B(Mat *_B) {
      PetscFunctionBegin;
      snes_A = _B;
      PetscFunctionReturn(0);
    }
    PetscErrorCode set_flag(MatStructure *_flag) {
      PetscFunctionBegin;
      snes_flag = _flag;
      PetscFunctionReturn(0);
    }
  };

  struct BasicMethod: public SnesMethod {
    BasicMethod();    

    PetscErrorCode set_moabfields(const MoFEMField_multiIndex *_moabfields);
    PetscErrorCode set_dofs_multiIndex(const DofMoFEMEntity_multiIndex *_dofs_moabfield);
    PetscErrorCode set_fes_multiIndex(const MoFEMFE_multiIndex *_finite_elements);
    PetscErrorCode set_fes_data_multiIndex(const EntMoFEMFE_multiIndex *_finite_elements_data);
    PetscErrorCode set_adjacencies(const MoFEMAdjacencies_multiIndex *_fem_adjacencies);
    //
    const MoFEMField_multiIndex *moabfields;
    const DofMoFEMEntity_multiIndex *dofs_moabfield;
    const MoFEMFE_multiIndex *finite_elements;
    const EntMoFEMFE_multiIndex *finite_elements_data;
    const MoFEMAdjacencies_multiIndex *fem_adjacencies;
  };

  /**
    * \brief structure for finite element method
    *
    * It can be used to calulate stiffnes matrices, residuals, load vectors etc.
    */  
  struct FEMethod: public BasicMethod {
    Interface& moab;
    FEMethod(Interface& _moab);

    /** \brief function is run at the beginig of looop
     *
     * It is used to zeroing matrices and vectors, clalualtion of shape
     * functions on reference element, preporocessing boundary conditions, etc.
     */
    virtual PetscErrorCode preProcess();

    /** \brief function is run for every finite element 
     *
     * It is used to calulate element local matrices and assembly. It can be
     * used for post-processing.
     */
    virtual PetscErrorCode operator()();

    /** \brief function is run at the end of looop
     *
     * It is used to assembly matrices and vectors, calulating global varibles,
     * f.e. total internal energy, ect.
     */
    virtual PetscErrorCode postProcess();

    PetscErrorCode set_problem(const MoFEMProblem *_problem_ptr);
    PetscErrorCode set_fe(const NumeredMoFEMFE *_fe_ptr); 
    PetscErrorCode set_data_multIndex(const FEDofMoFEMEntity_multiIndex *_data_multiIndex);
    PetscErrorCode set_row_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_row_multiIndex);
    PetscErrorCode set_col_multIndex(const FENumeredDofMoFEMEntity_multiIndex *_col_multiIndex);
    const MoFEMProblem *problem_ptr;
    const NumeredMoFEMFE *fe_ptr;
    const FEDofMoFEMEntity_multiIndex *data_multiIndex;
    const FENumeredDofMoFEMEntity_multiIndex *row_multiIndex;
    const FENumeredDofMoFEMEntity_multiIndex *col_multiIndex;
  };

  struct EntMethod: public BasicMethod {
    Interface& moab;
    EntMethod(Interface& _moab);
    
    virtual PetscErrorCode preProcess();
    virtual PetscErrorCode operator()();
    virtual PetscErrorCode postProcess();
    
    PetscErrorCode set_problem(const MoFEMProblem *_problem_ptr);
    PetscErrorCode set_dof(const NumeredDofMoFEMEntity *_dof_ptr);
    const MoFEMProblem *problem_ptr;
    const NumeredDofMoFEMEntity *dof_ptr;
  };

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
   * \param problem_name \param fe_name \param method is class derived form
   * moabField::FEMethod
  **/ virtual PetscErrorCode loop_finite_elements(const string
&problem_name,const string &fe_name,FEMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over entities
    *
    */
  virtual PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1) = 0;

};

}

#endif // __MOABFIELD_HPP__
