/** \file PrismInterface.hpp
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

#ifndef __PRISMINTERFACE_HPP__
#define __PRISMINTERFACE_HPP__

#include "FieldUnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMPrismInterface = MOFEMuuid( BitIntefaceId(PRISM_INTEFACE) );

struct PrismInterface: public FieldUnknownInterface {

  ///destructor
  virtual ~PrismInterface() {}

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
   * not splitting faces. Inhereteing those nodes will not split faces.
   *
   */
  virtual PetscErrorCode get_msId_3dENTS_split_sides(
    const EntityHandle meshset,const BitRefLevel &bit,
    const BitRefLevel &inheret_from_bit_level,const BitRefLevel &inheret_from_bit_level_mask,
    const EntityHandle SIDESET,const bool add_iterfece_entities,const bool recursive = false,int verb = -1) = 0;



};

}

#endif // __PRISMINTERFACE_HPP__


