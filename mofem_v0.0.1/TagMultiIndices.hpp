/** \file common.hpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk. 
 * It can be freely used for educational and research purposes 
 * by other institutions. If you use this softwre pleas cite my work. 
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

#ifndef __TAGMULTIINDICES_HPP__
#define __TAGMULTIINDICES_HPP__

namespace MoFEM {

  /// MultiIndex Tag for field id 
  struct CubitMeshSets_mi_tag {};
  struct CubitMeshSets_mask_meshset_mi_tag {};
  struct CubitMeshSets_bc_data_mi_tag {};
  struct CubitMeshSets_name {};
  struct BitFieldId_mi_tag {};
  struct Unique_mi_tag {};
  struct MoABEnt_mi_tag {};
  struct EntType_mi_tag {};
  struct Composite_unique_mi_tag {};
  struct Composite_mi_tag {};
  struct Composite_mi_tag2 {};
  struct Composite_mi_tag3 {};
  struct MoFEMFiniteElement_Meshset_mi_tag {};
  struct BitFEId_mi_tag {};
  struct MoFEMFiniteElement_name_mi_tag {};
  struct SideNumber_mi_tag{};
  struct MoABEnt_MoABEnt_mi_tag {};
  struct MoABEnt_mi_tag2 {};
  struct Idx_mi_tag { 
    static const bool IamNotPartitioned;
    /// extract dof index from iterator 
    template<class IT>
    static DofIdx get_index(const IT &it) { return it->dof_idx; }
  };
  struct PetscGlobalIdx_mi_tag {};
  struct PetscLocalIdx_mi_tag {};
  struct Part_mi_tag {
    static const bool IamNotPartitioned;
    /// extract global dof index from iterator 
    template<class IT>
    static DofIdx get_index(const IT &it) { return it->petsc_gloabl_dof_idx; }
  };
  struct Unique_MoABEnt_mi_tag {};
  struct Unique_MoFEMFiniteElement_mi_tag {};
  struct MoABEnt_MoFEMFiniteElement_mi_tag {};
  struct Meshset_mi_tag {};
  /// MultiIndex Tag for field name
  struct FieldName_mi_tag {};
  struct BitFieldId_space_mi_tag {};
  struct MoFEMFiniteElement_Part_mi_tag {};
  struct BitProblemId_mi_tag {};
  struct MoFEMProblem_mi_tag {};

  struct Composite_EntityType_And_ParentEntityType_mi_tag {};
  struct Composite_EntityHandle_And_ParentEntityType_mi_tag {};
  struct Composite_of_ParentEnt_And_BitsOfRefinedEdges_mi_tag {};
  struct Composite_Name_And_Ent_And_EndDofIdx {};
  struct Composite_Name_And_Ent {};
  struct Composite_Name_And_Type {};
  struct Composite_Name_Type_And_Side_Number {};
  struct Composite_Name_Ent_And_Part {};

}

#endif // __TAGMULTIINDICES_HPP__
