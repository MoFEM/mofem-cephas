/** \file NodeMerger.hpp
 * \brief NodeMerger interface 
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

#ifndef __NODE_MERGER_HPP__
#define __NODE_MERGER_HPP__

namespace MoFEM {


static const MOFEMuuid IDD_MOFENNodeMerger = MOFEMuuid( BitIntefaceId(NODEMERGER_INTERFACE) );

struct NodeMergerInterface: public FieldUnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, FieldUnknownInterface** iface);

  MoFEM::Core& cOre;
  NodeMergerInterface(MoFEM::Core& core): cOre(core) {};

  PetscErrorCode mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,Range *tets = NULL);
  PetscErrorCode mergeNodes(EntityHandle father,EntityHandle mother,BitRefLevel bit,BitRefLevel tets_from_bit_ref_level);

};

}

#endif //__NODE_MERGER_HPP__
