/** \file TagMultiIndices.hpp
 * \brief Tags for Multi-index containers
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef __TAGMULTIINDICES_HPP__
#define __TAGMULTIINDICES_HPP__

namespace MoFEM {

/// MultiIndex Tag for field id
struct CubitMeshsetType_mi_tag {};
struct CubitMeshsetMaskedType_mi_tag {};
struct CubitMeshsets_name {};
struct Composite_Cubit_msId_And_MeshsetType_mi_tag {};

struct BitFieldId_mi_tag {};
struct Unique_mi_tag {};
struct DOF_Unique_mi_tag {};
struct FE_Unique_mi_tag {};
struct Ent_mi_tag {};
struct FEEnt_mi_tag {};
struct EntType_mi_tag {};
struct FiniteElement_Meshset_mi_tag {};
struct BitFEId_mi_tag {};
struct FiniteElement_name_mi_tag {};
struct SideNumber_mi_tag {};
struct EntDofIdx_mi_tag {};
struct Space_mi_tag {};

struct Idx_mi_tag {
  static const bool IamNotPartitioned;
  /// extract dof index from iterator
  template <class IT> static inline DofIdx get_index(const IT &it) {
    return (*it)->getDofIdx();
  }
};
struct PetscGlobalIdx_mi_tag {
  static const bool IamNotPartitioned;
  /// extract global dof index from iterator
  template <class IT> static inline DofIdx get_index(const IT &it) {
    return (*it)->getPetscGlobalDofIdx();
  }
};
struct PetscLocalIdx_mi_tag {
  static const bool IamNotPartitioned;
  /// extract global dof index from iterator
  template <class IT> static inline DofIdx get_index(const IT &it) {
    return (*it)->getPetscLocalDofIdx();
  }
};

struct Part_mi_tag {};

struct Ent_Ent_mi_tag {};
struct Ent_Owner_mi_tag {};

struct Unique_Ent_mi_tag {};
struct Unique_FiniteElement_mi_tag {};
struct Ent_FiniteElement_mi_tag {};
struct Meshset_mi_tag {};

/// MultiIndex Tag for field order
struct Order_mi_tag {};

/// MultiIndex Tag for field name
struct FieldName_mi_tag {};
struct BitFieldId_space_mi_tag {};
struct BitProblemId_mi_tag {};
struct Problem_mi_tag {};

struct Ent_ParallelStatus {};
struct Proc_mi_tag {};

struct Composite_mi_tag {};
struct Composite_Unique_mi_tag {};
struct Composite_EntType_and_ParentEntType_mi_tag {};
struct Composite_ParentEnt_And_EntType_mi_tag {};
struct Composite_ParentEnt_And_BitsOfRefinedEdges_mi_tag {};
struct Composite_Name_And_Ent_And_EntDofIdx_mi_tag {};
struct Composite_Ent_And_EntDofIdx_mi_tag {};
struct Composite_Name_And_Ent_mi_tag {};
struct Composite_Part_And_Order_mi_tag {};
struct Composite_Name_Ent_Order_And_CoeffIdx_mi_tag {};
struct Composite_Ent_Order_And_CoeffIdx_mi_tag {};
struct Composite_Name_Ent_And_Part_mi_tag {};
struct Composite_Name_And_Part_mi_tag {};
struct Composite_Ent_and_ShortId_mi_tag {};
struct Composite_EntType_and_Space_mi_tag {};

struct SeriesID_mi_tag {};
struct SeriesName_mi_tag {};
struct Composite_SeriesID_And_Step_mi_tag {};
struct Composite_SeriesName_And_Step_mi_tag {};
struct Composite_SeriesName_And_Time_mi_tag {};

} // namespace MoFEM

#endif // __TAGMULTIINDICES_HPP__
