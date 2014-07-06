/** \file FieldInterface.hpp
 * \brief MoFEM interface 
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
    * \brief check data consistency in entitiesPtr
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_field(const string& name) = 0;

  /** 
    * \brief check data consistency in entitiesPtr
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_field() = 0;

  /** 
    * \brief check data consistency in finiteElementsMoFEMEnts
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_finite_element(const string& name) = 0;

  /** 
    * \brief check data consistency in finiteElementsMoFEMEnts
    *
    */
  virtual PetscErrorCode check_number_of_ents_in_ents_finite_element() = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief check for CUBIT Id and CUBIT type
    *
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    */
  virtual bool check_msId_meshset(const int msId,const CubitBC_BitSet CubitBCType) = 0;

  /**
    * \ingroup mofem_bc 
    * \brief add cubit meshset
    *
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    *
    */
  virtual PetscErrorCode add_Cubit_msId(const CubitBC_BitSet CubitBCType,const int msId) = 0;

  /**
    * \ingroup mopfem_bc
    * \brief delete cubit meshset
    *
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    *
    */
  virtual PetscErrorCode delete_Cubit_msId(const CubitBC_BitSet CubitBCType,const int msId) = 0;

  /** 
    * \ingroup mofem_bc
    * \brief get cubit meshset
    */
  virtual PetscErrorCode get_Cubit_msId(const int msId,const CubitBC_BitSet CubitBCType,const CubitMeshSets **cubit_meshset_ptr) = 0;

  /** 
    * \ingroup mofem_bc
    * \brief get entities from CUBIT/meshset of a particular entity dimension \n
	  * Nodeset can contain nodes, edges, triangles and tets. This applies to other CubitBCType meshsets too. \n
	  * The nodeset's meshset contain the nodes in the MIDDLE of the surface or volume which is done by default in Cubit,\n
		* Hence if all nodes on a particular nodeset are required,\n
	  * one should get all triangles or tetrahedrals for which the nodeset was create in Cubit,\n
	  * and get all the connectivities of tris/tets.
		*
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    * \param dimensions (0 - Nodes, 1 - Edges, 2 - Faces, 3 - Volume(tetrahedral))
    * \param Range containing the retreived entities
    * \param recursive If true, meshsets containing meshsets are queried recursively. Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, const int dimension,Range &entities,const bool recursive = false) = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief get entities related to CUBIT/meshset, \n
	  * NODESET will get Vertices only, even if the NODESET contains egdes, tris and tets\n
    * SIDESET will get Tris, BLOCKSET will get Tets, DISPLACEMENTSET and FORCESET are stored in NODESET, PRESSURESET is stored in Sideset.
	  *
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    * \param Range containing the retreived entities related to the CubitBCType
    * \param recursive If true, meshsets containing meshsets are queried recursively.  Returns the contents of meshsets, but not the meshsets themselves if true.
    */
  virtual PetscErrorCode get_Cubit_msId_entities_by_dimension(const int msId,const unsigned int CubitBCType, Range &entities,const bool recursive = false) = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief get meshset from CUBIT Id and CUBIT type
    *
    * \param msId id of the BLOCKSET/SIDESET/BLOCKSET: from CUBIT
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
    * \param meshset where to store the retreived entities
    */
  virtual PetscErrorCode get_Cubit_msId_meshset(const int msId,const unsigned int CubitBCType,EntityHandle &meshset) = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief get all CUBIT meshsets by CUBIT type
    *
    * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more). 
    * \param meshsets is range of meshsets
    */
  virtual PetscErrorCode get_Cubit_meshsets(const unsigned int CubitBCType,Range &meshsets) = 0;

   /** 
    * \ingroup mofem_bc 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
    * 	...
    * }
    *
    */
  virtual CubitMeshSet_multiIndex::iterator get_CubitMeshSets_begin() = 0;

   /** 
    * \ingroup mofem_bc 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,it) {
    * 	...
    * }
    *
    */
  virtual CubitMeshSet_multiIndex::iterator get_CubitMeshSets_end() = 0;
    
    /**
     * \ingroup mofem_bc 
     * \brief Iterator that loops over all the Cubit MeshSets in a moFEM field
     *
     * \param mField moFEM Field
     * \param iterator 
     */

  #define _IT_CUBITMESHSETS_FOR_LOOP_(MFIELD,IT) \
    CubitMeshSet_multiIndex::iterator IT = MFIELD.get_CubitMeshSets_begin(); IT!=MFIELD.get_CubitMeshSets_end(); IT++

  /**
    * \ingroup mofem_bc 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator get_CubitMeshSets_end(const unsigned int CubitBCType) = 0;
    
    /**
      * \ingroup mofem_bc 
      * \brief Iterator that loops over a specific Cubit MeshSet in a moFEM field
      *
      * \param mField moFEM Field
      * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
      * \param iterator 
      */

  #define _IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_end(CUBITBCTYPE); IT++

  /** 
    * \ingroup mofem_bc 
    * \brief get begin iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_begin(const unsigned int CubitBCType) = 0;

  /** 
    * \ingroup mofem_bc 
    * \brief get end iterator of cubit mehset of given type (instead you can use _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT)
    *
    * for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
    * 	...
    * }
    *
    * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
    */
  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator get_CubitMeshSets_bySetType_end(const unsigned int CubitBCType) = 0;
    
  /**
   * \ingroup mofem_bc 
   * \brief Iterator that loops over a specific Cubit MeshSet having a particular BC meshset in a moFEM field
   *
   * \param mField moFEM Field
   * \param CubitBCType see CubitBC (NODESET, SIDESET or BLOCKSET and more) 
   * \param iterator 
   *
   * Example: \code
     for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,NODESET|DISPLACEMENTSET,it) {
    	...
   * } \endcode
   */
  #define _IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(MFIELD,CUBITBCTYPE,IT) \
    CubitMeshSet_multiIndex::index<CubitMeshSets_mask_meshset_mi_tag>::type::iterator IT = MFIELD.get_CubitMeshSets_bySetType_begin(CUBITBCTYPE); \
    IT!=MFIELD.get_CubitMeshSets_bySetType_end(CUBITBCTYPE); IT++


  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_begin(const string& name) = 0;
  virtual CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator get_CubitMeshSets_byName_end(const string& name) = 0;

  /**
    * \ingroup mofem_bc 
    * \brief Iterator that loops over Cubit BlockSet having a particular name
    *
    * \param MFIELD mField
    * \param NAME name
    * \param IT iterator 
    *
    * Example: \code
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,"SOME_BLOCK_NANE",it) {
    	...
    * } \endcode
    */
  #define _IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    CubitMeshSet_multiIndex::index<CubitMeshSets_name>::type::iterator IT = MFIELD.get_CubitMeshSets_byName_begin(NAME); \
    IT!=MFIELD.get_CubitMeshSets_byName_end(NAME); IT++

  virtual PetscErrorCode print_cubit_displacement_set() = 0;
  virtual PetscErrorCode print_cubit_pressure_set() = 0;
  virtual PetscErrorCode print_cubit_force_set() = 0;
  virtual PetscErrorCode printCubitTEMPERATURESET() = 0;
  virtual PetscErrorCode printCubitHeatFluxSet() = 0;
  virtual PetscErrorCode print_cubit_materials_set() = 0;

  virtual PetscErrorCode rebuild_database(int verb = -1) = 0;

  /**
    * \ingroup mofem_series
    * Add series recorder 
    *
    * \param name of series
    */
  virtual PetscErrorCode add_series_recorder(const string& serie_name) = 0;

  /**
    * \ingroup mofem_series
    * initialize sries for recording
    *
    * \param series name
    */
  virtual PetscErrorCode initialize_series_recorder(const string& serie_name) = 0;

  /**
    * \ingroup mofem_series
    * finalize series for recording, recorded data are not accessible until finializd
    *
    * \param series name
    */
  virtual PetscErrorCode finalize_series_recorder(const string& serie_name) = 0;

  /**
    * \ingroup mofem_series
    * begin series recording
    * 
    * \param series name
    */
  virtual PetscErrorCode record_begin(const string& serie_name) = 0;

  /**
    * \ingroup mofem_series
    * record problem 
    *
    * \param series name
    * \param problem pointer
    * \param rc could be Row or Col
    */
  virtual PetscErrorCode record_problem(const string& serie_name,const MoFEMProblem *problemPtr,RowColData rc) = 0;

  /**
    * \ingroup mofem_series
    * record problem
    *
    * \param series name
    * \param problem name
    * \param rc could be Row or Col
    */
  virtual PetscErrorCode record_problem(const string& serie_name,const string& problem_name,RowColData rc) = 0;

  /**
    * \ingroup mofem_series
    * record field
    * 
    * \param field name
    * \param bit ref level
    * \param mask for bir ref level
    */
  virtual PetscErrorCode record_field(const string& serie_name,const string& field_name,const BitRefLevel &bit,const BitRefLevel &mask) = 0;

  /**
    * \ingroup mofem_series
    * end series recording 
    *
    * \param series name
    */
  virtual PetscErrorCode record_end(const string& serie_name) = 0;
  
  /**
    * \ingroup mofem_series
    * load data from series into dofs database
    * 
    * \param series name 
    * \param step number 
    */
  virtual PetscErrorCode load_series_data(const string& serie_name,const int step_number) = 0;

  /**
    * \ingroup mofem_series
    * print series
    */
  virtual PetscErrorCode print_series_steps() = 0;

  /** \brief check if series is in database
   * \ingroup mofem_series
   *
   * \param name field name
   * \return true if field exist
   *
   */
  virtual bool check_series(const string& name) const = 0;

  virtual SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_begin(const string& name) = 0;
  virtual SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_end(const string& name) = 0;

  /** \brief loop over recorded series step
    * \ingroup mofem_series
    *
    * \param NAME series name
    * \parma IT iterator variable
    *
    * Example: \code
      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(mField,"TEST_SERIES1",sit)) {

	ierr = mField.load_series_data("TEST_SERIES1",sit->get_step_number()); CHKERRQ(ierr);
    * } \endcode
    *
    */
  #define _IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator IT = MFIELD.get_series_steps_byName_begin(NAME); \
    IT!=MFIELD.get_series_steps_byName_end(NAME); IT++

  /**
  * Create finite elements based from eneties in meshses. Throw error if entity is not in database
  * 
  * \param EntityHandle meshset
  *
  */
  virtual PetscErrorCode seed_finiteElementsPtr(const EntityHandle meshset,int verb = -1) = 0;

  /**
  * Create finite elements based from eneties in meshses. Throw error if entity is not in database
  * 
  * \param Range entities
  *
  */
  virtual PetscErrorCode seed_finiteElementsPtr(const Range &entities,int verb = -1) = 0;

  /**
  * \brief seed 2D entities (Triangles entities only) in the meshset and their adjacencies (only TRIs adjencies) in a particular BitRefLevel
  * 
  * \param EntityHandle MeshSet
  * \param BitRefLevel bitLevel
  * 
  */
  virtual PetscErrorCode seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb = -1) = 0;

  /**
  * \brief seed 2D entities (Triangles entities only) in the meshset and their adjacencies (only TRIs adjencies) in a particular BitRefLevel
  * 
  * \param Range of tris
  * \param BitRefLevel bitLevel
  * 
  */
  virtual PetscErrorCode seed_ref_level_2D(const Range &ents2d,const BitRefLevel &bit,int verb = -1) = 0;

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

  /**
   * \brief seed 3D entities (Volume entities only) in the meshset and their adjacencies (only TETs adjencies) in a particular BitRefLevel
   */ 
  virtual PetscErrorCode seed_ref_level_3D(const Range &ents3d,const BitRefLevel &bit,int verb = -1) = 0;

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
  virtual PetscErrorCode refine_TET(const EntityHandle meshset,const BitRefLevel &bit,const bool respect_interface = false) = 0;

  /**\brief refine TET in the meshset
   *
   * \param Range of tets to refine
   * \param BitRefLevel bitLevel
   * \param If TRUE, interface elements would be refined too
   */
  virtual PetscErrorCode refine_TET(const Range &tets,const BitRefLevel &bit,const bool respect_interface = false) = 0;

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

  /**\brief add all ents from ref level given by bit to meshset
    * \ingroup mofem_ref_ents
    *
    * \param BitRefLevel bitLevel
    * \param BitRefLevel mask
    * \param EntityType type of entities
    * \param EntityHandle meshset   
    *
    */
  virtual PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const EntityHandle meshset,int verb = -1) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityType type of entities
   * \param Range ents   
   *
   */
  virtual PetscErrorCode get_entities_by_type_and_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,Range &ents,int verb = -1) = 0;


  /**\brief add all ents from ref level given by bit to meshset
   * \ingroup mofem_ref_ents
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param EntityHandle meshset   
   *
   */
  virtual PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,const EntityHandle meshset) = 0;

  /**\brief add all ents from ref level given by bit to meshset
   *
   * \param BitRefLevel bitLevel
   * \param BitRefLevel mask
   * \param Range   
   *
   *
   */
  virtual PetscErrorCode get_entities_by_ref_level(const BitRefLevel &bit,const BitRefLevel &mask,Range &ents) = 0;

  /**\brief add ref level to entities
   *
   */
  virtual PetscErrorCode add_ref_level_to_entities(const BitRefLevel &bit,Range &ents) = 0;

  /**\brief add ref level to entities
   *
   */
  virtual PetscErrorCode set_ref_level_to_entities(const BitRefLevel &bit,Range &ents) = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specfied dimension.
    * \ingroup mofem_ref_ents
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies_equality(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specfied dimension.
    * \ingroup mofem_ref_ents
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies_any(const EntityHandle from_entiti,const int to_dimension,Range &adj_entities) = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specfied dimension.
    * \ingroup mofem_ref_ents
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies(
    const MoFEMProblem *problemPtr,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0) = 0;

  /** \brief Get the adjacencies associated with a entity to entities of a specfied dimension.
    * \ingroup mofem_ref_ents
    *
    * bit ref level of adjacent entities is equal to bit ref level of adjacent entities
    */
  virtual PetscErrorCode get_adjacencies(
    const BitRefLevel &bit,
    const EntityHandle *from_entities,const int num_netities,const int to_dimension,Range &adj_entities,const int operation_type = Interface::INTERSECT,const int verb = 0) = 0;


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
  virtual PetscErrorCode update_meshset_by_entities_children(
    const EntityHandle parent, const BitRefLevel &child_bit,const EntityHandle child, EntityType child_type,
    const bool recursive = false, int verb = -1) = 0;

  /** \brief update fields meshesets by child entities
    * \ingroup mofem_field
    */
  virtual PetscErrorCode update_field_meshset_by_entities_children(const BitRefLevel &child_bit,int verb = -1) = 0;

  /** \brief update field mesheset by child entities
    * \ingroup mofem_field
    */
  virtual PetscErrorCode update_field_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,int verb = -1) = 0;

  /** \brief update finite element mesheset by child entities
    */
  virtual PetscErrorCode update_finite_element_meshset_by_entities_children(const string name,const BitRefLevel &child_bit,const EntityType fe_ent_type,int verb = -1) = 0;

  /** \brief delete enttities form mofem and moab database 
    */
  virtual PetscErrorCode delete_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const bool remove_parent = false,int verb = -1) = 0;

  /** \brief remove entities form mofem database
    */
  virtual PetscErrorCode remove_ents_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief remove finite element from mofem database
    */
  virtual PetscErrorCode delete_finiteElementsPtr_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief left shift bit ref level
    * this results of deletion of enetities on far left side
    */
  virtual PetscErrorCode shift_left_bit_ref(const int shif,int verb = -1) = 0;

  /** \brief right shift bit ref level
    *
    */
  virtual PetscErrorCode shift_right_bit_ref(const int shift,int verb = -1) = 0;

  /** 
    * \brief add approximation field
    * \ingroup mofem_field
    *
    * \param name of the field
    * \param space approximation space (H1, Hdiv, Hcurl, L2 and NoField (dofs adjacent to meshset) 
    * \prama rank of the field, f.e. temperature has rank 1, displacement in 3d has rank 3
    */
  virtual PetscErrorCode add_field(const string& name,const FieldSpace space,const ApproximationRank rank,enum MoFEMTypes bh = MF_EXCL,int verb = -1) = 0;

  /** 
    * \brief set field entities on vertices
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param nodes contains set vertices
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const Range &nodes,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities on vertices
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set vertices
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities form adjacencies of edges
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_EDGEs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities form adjacencies of triangles
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set triangles
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TRIs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities from adjacencies of tetrahedrals
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param meshset contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const EntityHandle meshset,const string& name,int verb = -1) = 0;

  /** 
    * \brief set field entities from adjacencies of tetrahedrals
    * \ingroup mofem_field
    *
    * The lower dimension entities are added depending on the space type
    * \param range contains set tetrahedrals
    * \param name of the field
    */
  virtual PetscErrorCode add_ents_to_field_by_TETs(const Range &tets,const string& name,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field(const string& name,const EntityHandle meshset,const EntityType type,int verb = -1) = 0;

  /**
    * \brief remove entities from field
    * \ingroup mofem_field
    *
    */
  virtual PetscErrorCode remove_ents_from_field(const string& name,const Range &ents,int verb = -1) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param meshset containing set of the entities (use 0 for all the entities in the meshset)
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order 
    */
  virtual PetscErrorCode set_field_order(const EntityHandle meshset,const EntityType type,const string& name,const ApproximationOrder order,int verb = -1) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param entities 
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order 
    */
  virtual PetscErrorCode set_field_order(const Range &ents,const string& name,const ApproximationOrder order,int verb = -1) = 0;

  /**
    * \brief Set order approximation of the entities in the field
    * \ingroup mofem_field
    *
    * \param bit refinement level
    * \param mask bit mask
    * \param type selected type of the entities f.e. MBTET, MBTRI, MBEDGE, MBVERTEX, see moab documentation
    * \param order approximation order 
    */
  virtual PetscErrorCode set_field_order_by_entity_type_and_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const EntityType type,const string& name,const ApproximationOrder order,int verb = -1) = 0;


  /** \brief list entities in the field
    * \ingroup mofem_field 
    */
  virtual PetscErrorCode list_fields() const = 0;

  /** \brief get field meshsett
   * \ingroup mofem_field
   *
   * \param name of Field
   * Example:\code
   EntityHandle disp_files_meshset = mField.get_field_meshset("DISPLACEMENT");
   * \endcode
   */    
  virtual EntityHandle get_field_meshset(const string& name) const = 0;

  /** \brief check if field is in database
   *
   * \param name field name
   * \return true if field exist
   *
   */
  virtual bool check_field(const string& name) const = 0;

  /** \brief get field structure
   * \ingroup mofem_field
   *
   * \param name field name
   * \return const MoFEMField*
   *
   */
  virtual const MoFEMField* get_field_structure(const string& name) = 0;

  /**
    * \brief add finite element
    * \param name finite element name
    *
    * Example \code
      ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
      ierr = mField.add_finite_element("PLASTIC"); CHKERRQ(ierr);
   \endcode
    */
  virtual PetscErrorCode add_finite_element(const string &MoFEMFiniteElement_name,enum MoFEMTypes bh = MF_EXCL) = 0;

  /** \brief set field data which finite element use
   *
   * \param name finite element name
   * \param name field name
   *
   * This function will set memory in the form of a vector
   */
  virtual PetscErrorCode modify_finite_element_add_field_data(const string &MoFEMFiniteElement_name,const string &name_filed) = 0;
  virtual PetscErrorCode modify_finite_element_off_field_data(const string &MoFEMFiniteElement_name,const string &name_filed) = 0;

    /** \brief set field row which finite element use
     *
     * \param name finite element name
     * \param name field name
     */
  virtual PetscErrorCode modify_finite_element_add_field_row(const string &MoFEMFiniteElement_name,const string &name_row) = 0;
  virtual PetscErrorCode modify_finite_element_off_field_row(const string &MoFEMFiniteElement_name,const string &name_row) = 0;


    /** \brief set field col which finite element use
     *
     * \param name finite element name
     * \param name field name
     */  
  virtual PetscErrorCode modify_finite_element_add_field_col(const string &MoFEMFiniteElement_name,const string &name_row) = 0;
  virtual PetscErrorCode modify_finite_element_off_field_col(const string &MoFEMFiniteElement_name,const string &name_row) = 0;

  /** \brief add EDGES entities fromm meshset to finite element database given by name
   *
   * \param range contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_EDGEs(const Range& edge,const string &name) = 0;

  /** \brief add VERTICES entities fromm meshset to finite element database given by name
   *
   * \param range contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_VERTICEs(const Range& vert,const string &name) = 0;

  /** \brief add TRI entities fromm meshset to finite element database given by name
   *
   * \param range contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_TRIs(const Range& tris,const string &name) = 0;


  /** \brief add TET entities fromm meshset to finite element database given by name
   *
   * \param range contains tetrahedrals
   * \param name Finite Element name
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_TETs(const Range& tets,const string &name) = 0;

  /** \brief add TET entities fromm meshset to finite element database given by name
   *
   * \param meshset contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_TETs(const EntityHandle meshset,const string &name,const bool recursive = false) = 0;

  /** \brief add PRISM entities fromm meshset to finite element database given by name
   *
   * \param range contains tetrahedrals
   * \param name Finite Element name
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_PRISMs(const Range& prims,const string &name) = 0;

  /** \brief add TET entities fromm meshset to finite element database given by name
   *
   * \param meshset contains tetrahedrals
   * \param name Finite Element name
   * \param recursive if true parent meshset is searched recursively
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_PRISMs(const EntityHandle meshset,const string &name,const bool recursive = false) = 0;

  /** \brief add TET elements from given refinment level to finite element database given by name 
   *
   * \param BitRefLevel bit
   * \param Finite Element name
   * \param Finite Elenent type
   * \param verrbose level
   */
  virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const string &name,EntityType type,int verb = -1) = 0;

  
  /** get finite element meshset
   *
   */
  virtual EntityHandle get_finite_element_meshset(const string& name) const = 0;


  /** \brief remove elements from given refinment level to finite element database
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param verrbose level
   */
  virtual PetscErrorCode remove_ents_from_finite_element_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief remove elements from given refinment level to finite element database
   *
   */
  virtual PetscErrorCode remove_ents_from_finite_element(const string &name,const EntityHandle meshset,const EntityType type,int verb = -1) = 0;

  /** \brief remove elements from given refinment level to finite element database
   *
   */
  virtual PetscErrorCode remove_ents_from_finite_element(const string &name,const Range &ents,int verb = -1) = 0;

  /** \brief add TET elements from given refinment level to finite element database given by name 
   *
   * \param BitRefLevel bit
   * \param BitRefLevel mask
   * \param Finite Element name
   * \param Finite Elenent type
   * \param verrbose level
   */
  virtual PetscErrorCode add_ents_to_finite_element_EntType_by_bit_ref(const BitRefLevel &bit,const BitRefLevel &mask,const string &name,EntityType type,int verb = -1) = 0;

  /** \brief add MESHSET element to finite element database given by name 
   *
   * \param meshset contains all entities that could be used for finite element
   * \param name Finite Element name
   */
  virtual PetscErrorCode add_ents_to_finite_element_by_MESHSET(const EntityHandle meshset,const string& name,const bool recursive = false) = 0;

  /// list finite elements in database
  virtual PetscErrorCode list_finiteElementsPtr() const = 0;

  /// list adjacencies
  virtual PetscErrorCode list_adjacencies() const = 0;

  /// add Finite Element Problem
  virtual PetscErrorCode add_problem(const string& name,enum MoFEMTypes bh = MF_EXCL,int verb = -1) = 0;
  
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

  /** \brief set ref level for problem
   *
   * if same finite element is solved using different level of refinements, than the level of refinement has to be specificied to problem in query
   *
   * \param name Problem name
   * \param BitRefLevel bitLevel
   * Example: \code
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF1","ELASTIC"); CHKERRQ(ierr);
   ierr = mField.modify_problem_add_finite_element("BEAM_BENDING_ON_MESH_REF2","ELASTIC"); CHKERRQ(ierr);
   
   ierr = mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF1",bit_level1); CHKERRQ(ierr);
   ierr = mField.modify_problem_ref_level_set_bit("BEAM_BENDING_ON_MESH_REF2",bit_level2); CHKERRQ(ierr);
   *\endcode
   * Two Problems exist and solved independently, both are elastic, but solved using different mesh refinement <br>
   */
  virtual PetscErrorCode modify_problem_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) = 0;

  /** \brief set dof mask ref level for problem
   *
   */
  virtual PetscErrorCode modify_problem_dof_mask_ref_level_set_bit(const string &name_problem,const BitRefLevel &bit) = 0;

  /// list problems
  virtual PetscErrorCode list_problem() const = 0;

  /** build fields
    * \ingroup mofem_field
   */
  virtual PetscErrorCode build_fields(int verb = -1) = 0;

  /** list dofs
    * \ingroup mofem_dofs
   */
  virtual PetscErrorCode list_dofs_by_field_name(const string &name,bool synchronised = false) const = 0;

  /** clear fields
    * \ingroup mofem_dofs
   */
  virtual PetscErrorCode clear_dofs_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** clear fields
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_ents_fields(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** clear fields
    * \ingroup mofem_dofs
   */
  virtual PetscErrorCode clear_dofs_fields(const string &name,const Range ents,int verb = -1) = 0;

  /** clear fields
    * \ingroup mofem_field
   */
  virtual PetscErrorCode clear_ents_fields(const string &name,const Range enst,int verb = -1) = 0;

  /** build finite elements
    */
  virtual PetscErrorCode build_finiteElementsPtr(int verb = -1) = 0;

  /** clear finite elements
    */
  virtual PetscErrorCode clear_finiteElementsPtr(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** clear finite elements
    */
  virtual PetscErrorCode clear_finiteElementsPtr(const string &name,const Range &ents,int verb = -1) = 0;

  /** \brief build ajacencies 
    *
    * \param bit adjacencies for refine level
    *
    * This function will get information of adjacent finite elements and fields
    * of all entities If this is not perfomed, partitioning the problem is not
    * possible. Adjacency map is based on degrees of freedom adjacent to
    * elements. This linked to gemetric element connectivity. 
    *
    * If new degrees of freedom or new finite elements are added to the
    * database, adjacency map has to be rebuild.
    *
    */
  virtual PetscErrorCode build_adjacencies(const BitRefLevel &bit,int verb = -1) = 0;

  /** \brief clear adjacency map for finite elements on given bit level
    *
    * \param bit
    * \param mask
    */
  virtual PetscErrorCode clear_adjacencies_finiteElementsPtr(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief clear adjacency map for entities on given bit level
    *
    * \param bit
    * \param mask
    */
  virtual PetscErrorCode clear_adjacencies_entities(const BitRefLevel &bit,const BitRefLevel &mask,int verb = -1) = 0;

  /** \brief build problem data structures
   */
  virtual PetscErrorCode build_problems(int verb = -1) = 0;

  /** \brief clear problems
   */
  virtual PetscErrorCode clear_problems(int verb = -1) = 0;


  /** \brief partition problem dofs
   *
   * \param name problem name
   */
  virtual PetscErrorCode simple_partition_problem(const string &name,const int all_on_part = -1,int verb = -1) = 0;


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
  virtual PetscErrorCode partition_finiteElementsPtr(const string &name,bool do_skip = true,int verb = -1) = 0;

  /** \brief check if matrix fill in correspond to finite element indices
    *
    */
  virtual PetscErrorCode partition_check_matrix_fill_in(const string &problem_name,int verb) = 0;

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
  virtual PetscErrorCode MatCreateSeqAIJWithArrays(const string &name,Mat *Aij,PetscInt **i,PetscInt **j,PetscScalar **v,int verb = -1) = 0;

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
    * \param pointer to porblem struture
    * \param RowColData for row or column (i.e. Row,Col)
    * \param V vector
    * \param mode see petsc manual for VecSetValue (ADD_VALUES or INSERT_VALUES)
    * \param scatter_mode see petsc manual for ScatterMode (The available modes are: SCATTER_FORWARD or SCATTER_REVERSE)
    *
    * SCATTER_REVERSE set data to field entities form V vector.
    *
    */
  virtual PetscErrorCode set_global_VecCreateGhost(const MoFEMProblem *problemPtr,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode) = 0;

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
    * \param pointer to poroblem multi_index 
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
  virtual PetscErrorCode set_other_local_VecCreateGhost(
    const MoFEMProblem *problemPtr,const string& fiel_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1) = 0;

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
  virtual PetscErrorCode set_other_local_VecCreateGhost(
    const string &name,const string& field_name,const string& cpy_field_name,RowColData rc,Vec V,InsertMode mode,ScatterMode scatter_mode,int verb = -1) = 0;

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

  /** \brief axpy fields 
    * \ingroup mofem_field_operators
    *
    * field_y = field_y + alpha*field_x
    *
    * \param alpha
    * \param field_name_x name of field_x
    * \param field_name_y name of field_y
    * \param error_if_missing throw error if entity/dof exist in field_x but not on field_y
    * \param create_if_missing creat dof in field_y from fiedl_x if it is not database
    *
    */
  virtual PetscErrorCode field_axpy(const double alpha,const string& fiel_name_x,const string& field_name_y,bool error_if_missing = false,bool creat_if_missing = false) = 0;

  /** \brief scale field
    * \ingroup mofem_field_operators
    * 
    * \param alpha is a scaling factor
    * \field_name  is a field name
    *
    */
  virtual PetscErrorCode field_scale(const double alpha,const string& field_name) = 0;

  /** \brief set field 
    * \ingroup mofem_field_operators
    *
    * field_y = val
    *
    * \param val
    * \param entity type
    * \param field_name
    *
    */
  virtual PetscErrorCode set_field(const double val,const EntityType type,const string& field_name) = 0;

  //topologity
  /** \brief create two children meshsets in the meshset containing terahedrals on two sides of faces
    *
    * \param msId Id of meshset 
    * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
    * \param mesh_bit_level add interface on bit level is bit_level = BitRefLevel.set() then add interfece on all bit levels
    * \param recursive if true parent meshset is searched recursively
    */
  virtual PetscErrorCode get_msId_3dENTS_sides(
    const int msId,
    const CubitBC_BitSet CubitBCType,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1) = 0;

  /** \brief create two children meshsets in the meshset contaning terahedrals on two sides of faces
   *
   * Get tets adj to faces. Take skin form tets and get edges from that skin. Take skin form triangles (the face).
   * Subtrac skin faces edges form skin edges in order to get eges on the boundary of the face which is in the
   * voulume of the body, but is not on the boundary.
   * Each child set has a child contaning nodes which can be split and skin edges.
   * After that simply iterate under all tets on one side which are adjacent to the face are found.
   * Side tets are stored in to children meshsets of the SIDESET meshset.
   */
  virtual PetscErrorCode get_msId_3dENTS_sides(
    const EntityHandle SIDESET,
    const BitRefLevel mesh_bit_level,
    const bool recursive,int verb = -1) = 0;

  /**
   * \brief split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   * \param meshset meshset to get entities from
   * \param BitRefLevel new level where refinement would be stored
   * \param msId meshset ID imported from cubit 
   * \param CubitBCType type of meshset (NODESET, SIDESET or BLOCKSET and more)
   * \param add_intefece_entities meshset which contain the interface
   * \param recursive if true parent meshset is searched recursively
   *
   * Each inteface face has two tages, 
   *    const int def_side[] = {0};
   *	rval = moab.tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
   *	  th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERR_PETSC(rval);
   * 
   *  	const EntityHandle def_node[] = {0};
   *	rval = moab.tag_get_handle("SIDE_INTFACE_ELEMENT",1,MB_TYPE_HANDLE,
   *	  th_side_elem,MB_TAG_CREAT|MB_TAG_SPARSE,def_node); CHKERR_PETSC(rval);
   * 
   * First tag inform abot inteface side, second tag inform about side adjacent 
   * inteface element.
   *
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const int msId,const CubitBC_BitSet CubitBCType,
    const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;

  /**
   * \brief split nodes and other entities of tetrahedrals in children sets and add prism elements
   * 
   * The all new entities (prisms, tets) are added to refinment level given by bit
   *
   * \param meshset 
   * \param refinment bit level of new mesh
   * \param inheret_from_bit_level inheret nodes and other entities form this bit level. 
   * \param add_iterfece_entities add prism elements at interface
   * \param recuslsive do meshesets in the meshset
   * 
   * note inheret_from_bit_level is need to be specidied to some meshset
   * with interfaces. Some nodes on some refinment levels dividing edges but
   * not splitting faces. Inhereteing thise nodes will not split faces.
   *
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;


  /**
   * \brief data structure for snes (nonlinear solver) context
   * \ingroup mofem_loop_methods
   *
   * Struture stores context data which are set in finctions run by PETSc SNES functions.
   *
   */
  struct SnesMethod {
    enum SNESContext { CTX_SNESSETFUNCTION, CTX_SNESSETJACOBIAN, CTX_SNESNONE };
    //
    SNESContext snes_ctx;
    SnesMethod(): snes_ctx(CTX_SNESNONE) {};
    //
    PetscErrorCode set_snes_ctx(const SNESContext ctx_);
    //
    SNES snes;
    PetscErrorCode set_snes(SNES _snes);
    Vec snes_x,snes_f;
    Mat *snes_A,*snes_B;
    MatStructure *snes_flag;
    virtual ~SnesMethod() {};
  };

  /**
   * \brief data structure for ts (time stepping) context
   * \ingroup mofem_loop_methods
   *
   * Struture stores context data which are set in finctions run by PETSc Time Stepping functions.
   */
  struct TSMethod {
    enum TSContext { CTX_TSSETRHSFUNCTION, CTX_TSSETRHSJACOBIAN, CTX_TSSETIFUNCTION, CTX_TSSETIJACOBIAN, CTX_TSTSMONITORSET, CTX_TSNONE };
    //
    TSContext ts_ctx;
    TSMethod(): ts_ctx(CTX_TSNONE),ts_a(0),ts_t(0) {};
    //
    PetscErrorCode set_ts_ctx(const TSContext ctx_);
    //
    TS ts;
    PetscErrorCode set_ts(TS _ts);
    Vec ts_u,ts_u_t,ts_F;
    Mat *ts_A,*ts_B;
    MatStructure *ts_flag;
    //
    PetscInt ts_step;
    PetscReal ts_a,ts_t;
    virtual ~TSMethod() {};
  };

  /**
   * \brief Data strutucture to exchange data between mofem and User Loop Methods.
   * \ingroup mofem_loop_methods
   *
   * It allows to exchange data between MoFEM and user functoions. It stores informaton about multi-indices.
   */
  struct BasicMethod: public SnesMethod,TSMethod {
    BasicMethod();    
    //
    PetscErrorCode setProblem(const MoFEMProblem *_problemPtr);
    PetscErrorCode setFields(const MoFEMField_multiIndex *_fieldsPtr);
    PetscErrorCode setEnts(
      const  RefMoFEMEntity_multiIndex *_refinedEntitiesPtr,
      const MoFEMEntity_multiIndex *_entitiesPtr);
    PetscErrorCode setDofs(const DofMoFEMEntity_multiIndex *_dofsPtr);
    PetscErrorCode setFiniteElements(
      const RefMoFEMElement_multiIndex *_refinedFiniteElementsPtr,
      const MoFEMFiniteElement_multiIndex *_finiteElementsPtr);
    PetscErrorCode setFiniteElementsEntities(const EntMoFEMFiniteElement_multiIndex *_finiteElementsEntitiesPtr);
    PetscErrorCode setAdjacencies(const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *_adjacenciesPtr);
    //
    virtual PetscErrorCode preProcess() = 0;
    virtual PetscErrorCode operator()() = 0;
    virtual PetscErrorCode postProcess() = 0;
    //
    const RefMoFEMEntity_multiIndex *refinedEntitiesPtr;
    const RefMoFEMElement_multiIndex *refinedFiniteElementsPtr;
    const MoFEMProblem *problemPtr;
    const MoFEMField_multiIndex *fieldsPtr;
    const MoFEMEntity_multiIndex *entitiesPtr;
    const DofMoFEMEntity_multiIndex *dofsPtr;
    const MoFEMFiniteElement_multiIndex *finiteElementsPtr;
    const EntMoFEMFiniteElement_multiIndex *finiteElementsEntitiesPtr;
    const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap_multiIndex *adjacenciesPtr;
    virtual ~BasicMethod() {};
  };

  /**
    * \brief structure for User Loop Methods on finite elements
    * \ingroup mofem_loop_methods
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
    PetscErrorCode setFE(const NumeredMoFEMFiniteElement *_fe_ptr); 
    PetscErrorCode setData(const FEDofMoFEMEntity_multiIndex *_data_ptr);
    PetscErrorCode setRowData(const FENumeredDofMoFEMEntity_multiIndex *_row_ptr);
    PetscErrorCode setColData(const FENumeredDofMoFEMEntity_multiIndex *_col_ptr);
    string feName;
    const NumeredMoFEMFiniteElement *fePtr;
    const FEDofMoFEMEntity_multiIndex *dataPtr;
    const FENumeredDofMoFEMEntity_multiIndex *rowPtr;
    const FENumeredDofMoFEMEntity_multiIndex *colPtr;

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
       * \ingroup mofem_loop_methods
       *
       * \param FE finite elements
       * \param Name field name
       * \param Type moab entity type (MBVERTEX, MBEDGE etc)
       * \param Side side canonical number
       * \param IT the interator in use
       */
    #define _IT_GET_FEROW_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->rowPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->rowPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    /** \brief loop over all dofs which are on a particular FE column, field, entity type and canonical side number
      * \ingroup mofem_loop_methods
      */ 
    #define _IT_GET_FECOL_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->colPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->colPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    /** \brief loop over all dofs which are on a particular FE data, field, entity type and canonical side number
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEDATA_BY_SIDE_DOFS_FOR_LOOP_(FE,NAME,TYPE,SIDE,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->dataPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type>(FE->dataPtr->get<Composite_mi_tag>(),NAME,TYPE,SIDE); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
      return index.lower_bound(boost::make_tuple(field_name,type));
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityType type) const {
      return index.upper_bound(boost::make_tuple(field_name,type));
    } 

    /** \brief loop over all dofs which are on a particular FE row, field and entity type
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

    /** \brief loop over all dofs which are on a particular FE column, field and entity type
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FECOL_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

    /** \brief loop over all dofs which are on a particular FE data, field and entity type
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(FE,NAME,TYPE,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Type_mi_tag>(),NAME,TYPE); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name) const {
      return index.lower_bound(field_name);
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name) const {
      return index.upper_bound(field_name);
    } 

    /** \brief loop over all dofs which are on a particular FE row and field
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEROW_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->rowPtr->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->rowPtr->get<FieldName_mi_tag>(),NAME); IT++

    /** \brief loop over all dofs which are on a particular FE column and field
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FECOL_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->colPtr->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->colPtr->get<FieldName_mi_tag>(),NAME); IT++

    /** \brief loop over all dofs which are on a particular FE data and field
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEDATA_DOFS_FOR_LOOP_(FE,NAME,IT) \
    FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->dataPtr->get<FieldName_mi_tag>(),NAME); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type>(FE->dataPtr->get<FieldName_mi_tag>(),NAME); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const EntityHandle ent) const {
      return index.lower_bound(ent);
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const EntityHandle ent) const {
      return index.upper_bound(ent);
    } 

    /** \brief loop over all dofs which are on a particular FE row and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEROW_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
      FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->rowPtr->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->rowPtr->get<MoABEnt_mi_tag>(),ENT); IT++

    /** \brief loop over all dofs which are on a particular FE column and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FECOL_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->colPtr->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->colPtr->get<MoABEnt_mi_tag>(),ENT); IT++

    /** \brief loop over all dofs which are on a particular FE data and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEDATA_DOFS_BY_ENT_FOR_LOOP_(FE,ENT,IT) \
    FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->dataPtr->get<MoABEnt_mi_tag>(),ENT); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type>(FE->dataPtr->get<MoABEnt_mi_tag>(),ENT); IT++

    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_begin(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
      return index.lower_bound(boost::make_tuple(field_name,ent));
    } 
    template<class MULTIINDEX>
    typename MULTIINDEX::iterator get_end(const MULTIINDEX &index,const string &field_name,const EntityHandle ent) const {
      return index.upper_bound(boost::make_tuple(field_name,ent));
    } 

    /** \brief loop over all dofs which are on a particular FE row, field and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEROW_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->rowPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

    /** \brief loop over all dofs which are on a particular FE column, field and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FECOL_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
      IT = FE->get_begin<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
      IT != FE->get_end<FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->colPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

    /** \brief loop over all dofs which are on a particular FE data, field and given element entity (handle from moab)
      * \ingroup mofem_loop_methods
      */
    #define _IT_GET_FEDATA_DOFS_BY_NAME_AND_ENT_FOR_LOOP_(FE,NAME,ENT,IT) \
    FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator \
      IT = FE->get_begin<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); \
      IT != FE->get_end<FEDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type>(FE->dataPtr->get<Composite_Name_And_Ent_mi_tag>(),NAME,ENT); IT++

  };

  /**
   * \brief Data strutucture to exchange data between mofem and User Loop Methods on Entirties.
   * \ingroup mofem_loop_methods
   *
   * It allows to exchange data between MoFEM and user functoions. It stores informaton about multi-indices.
   */
  struct EntMethod: public BasicMethod {
    EntMethod();
    
    PetscErrorCode preProcess();
    PetscErrorCode operator()();
    PetscErrorCode postProcess();
 
    PetscErrorCode setDof(const DofMoFEMEntity *_dof_ptr);
    const DofMoFEMEntity *dof_ptr;
    PetscErrorCode set_numered_dof(const NumeredDofMoFEMEntity *_dof_ptr);
    const NumeredDofMoFEMEntity *dofPtr;
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
  virtual PetscErrorCode loop_finiteElementsPtr(const string &problem_name,const string &fe_name,FEMethod &method,int verb = -1) = 0;

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
  virtual PetscErrorCode loop_finiteElementsPtr(
    const string &problem_name,const string &fe_name,FEMethod &method,
    int lower_rank,int upper_rank,int verb = -1) = 0;

  /** \brief Make a loop over entities
    *
    */
  virtual PetscErrorCode loop_dofs(const string &problem_name,const string &field_name,RowColData rc,EntMethod &method,int verb = -1) = 0;

  /** \brief Make a loop over entities
    * \ingroup mofem_dofs
    *
    */
  virtual PetscErrorCode loop_dofs(const string &field_name,EntMethod &method,int verb = -1) = 0;


  /** \brief Get ref entities from database (datastructure) 
    *
    */
  virtual PetscErrorCode get_ref_ents(const RefMoFEMEntity_multiIndex **refinedEntitiesPtr_ptr) = 0;

  /** \brief Get problem database (datastructure) 
    *
    */
  virtual PetscErrorCode get_problem(const string &problem_name,const MoFEMProblem **problemPtr) = 0;

  /** \brief Get dofs multi index
    * \ingroup mofem_dofs
    *
    */
  virtual PetscErrorCode get_dofs(const DofMoFEMEntity_multiIndex **dofsPtr_ptr) = 0;

  /** 
    * \brief get begin iterator of filed ents of given name (instead you can use _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_dofs
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
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_ent_moabfield_by_name_end(const string &field_name) = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_dofs
    */
  #define _IT_GET_ENT_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MFIELD.get_ent_moabfield_by_name_begin(NAME); \
      IT != MFIELD.get_ent_moabfield_by_name_end(NAME); IT++

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_begin(const string &field_name) const = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator get_dofs_by_name_end(const string &field_name) const = 0;

  /** loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_dofs
    */
  #define _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,IT) \
    DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator IT = MFIELD.get_dofs_by_name_begin(NAME); \
      IT != MFIELD.get_dofs_by_name_end(NAME); IT++

  /** 
    * \brief get begin iterator of filed dofs of given name and ent(instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_begin(const string &field_name,const EntityHandle ent) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name and ent (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,ENT,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator get_dofs_by_name_and_ent_end(const string &field_name,const EntityHandle ent) = 0;

  ///loop over all dofs from a moFEM field and particular field
  #define _IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(MFIELD,NAME,ENT,IT) \
    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator IT = MFIELD.get_dofs_by_name_and_ent_begin(NAME,ENT); \
      IT != MFIELD.get_dofs_by_name_and_ent_end(NAME,ENT); IT++

  /** 
    * \brief get field data from entity and field
    * \ingroup mofem_dofs
    * 
    * this funciont is not recommended to be used in finite elemeny implementation
    *
    */
  template <typename DIT>
  PetscErrorCode get_FielData(const string& name,const EntityHandle *ent,const int num_ents,DIT dit,int *count = NULL) {
    PetscFunctionBegin;
    if(count!=NULL) *count = 0;
    for(int nn = 0;nn<num_ents;nn++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_((*this),name,ent[nn],it)) {
	*(dit++) = it->get_FieldData();
	if(count!=NULL) (*count)++;
      }
    }
    PetscFunctionReturn(0);
  }

  /** 
    * \brief get field data from entity and field
    * \ingroup mofem_dofs
    * 
    * this funciont is not recommended to be used in finite elemeny implementation
    *
    */
  template <typename DIT>
  PetscErrorCode get_FielData(const string& name,const Range &ents,DIT dit,int *count = NULL) {
    PetscFunctionBegin;
    if(count!=NULL) *count = 0;
    for(Range::const_iterator eit = ents.begin();eit!=ents.end();eit++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_((*this),name,*eit,it)) {
	*(dit++) = it->get_FieldData();
	if(count!=NULL) (*count)++;
      }
    }
    PetscFunctionReturn(0);
  }

  /** 
    * \brief get begin iterator of filed dofs of given name and ent type (instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_begin(const string &field_name,const EntityType type) = 0;

  /** 
    * \brief get begin iterator of filed dofs of given name end ent type(instead you can use _IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(MFIELD,NAME,TYPE,IT)
    * \ingroup mofem_dofs
    *
    * for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT)) {
    * 	...
    * }
    *
    * \param field_name  
    */
  virtual DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator get_dofs_by_name_and_type_end(const string &field_name,const EntityType type) = 0;

  /** \brief loop over all dofs from a moFEM field and particular field
    * \ingroup mofem_dofs
    */
  #define _IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(MFIELD,NAME,TYPE,IT) \
    DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag>::type::iterator IT = MFIELD.get_dofs_by_name_and_type_begin(NAME,TYPE); \
      IT != MFIELD.get_dofs_by_name_and_type_end(NAME,TYPE); IT++

  /** \brief Get finite elements multi index
    *
    */
  virtual PetscErrorCode get_finiteElementsPtr(const MoFEMFiniteElement_multiIndex **finiteElementsPtr_ptr) = 0;

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

/***************************************************************************//**
 * \defgroup mofem_bc Handling boundary conditions
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_series Recording and reading series
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_field Fields
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_dofs Dofs
 * \ingroup mofem_field
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_field_operators Oprators
 * \ingroup mofem_field
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_ref_ents Getting entities and adjacencies
 * \ingroup mofem
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup mofem_loop_methods Methods for Loops
 * \ingroup mofem
 ******************************************************************************/








