/** \file CoreDataStructures.cpp
 * \brief Myltindex containers, data structures and other low-level functions
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

const bool Idx_mi_tag::IamNotPartitioned = true;
const bool PetscGlobalIdx_mi_tag::IamNotPartitioned = false;
const bool PetscLocalIdx_mi_tag::IamNotPartitioned = false;
const bool Part_mi_tag::IamNotPartitioned = false;

//fields
MoFEMField::MoFEMField(Interface &moab,const EntityHandle meshset,const CoordSys *coord_sys_ptr):
  meshSet(meshset),
  coordSysPtr(coord_sys_ptr),
  tag_id_data(NULL),
  tag_space_data(NULL),
  tag_nb_coeff_data(NULL),
  tag_name_data(NULL),
  tag_name_size(0) {
  //Change those tags only by modifiers
  ErrorCode rval;
  //id
  Tag th_field_id;
  rval = moab.tag_get_handle("_FieldId",th_field_id); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_field_id,&meshSet,1,(const void **)&tag_id_data); CHKERR_THROW(rval);
  //space
  Tag th_field_space;
  rval = moab.tag_get_handle("_FieldSpace",th_field_space); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_field_space,&meshSet,1,(const void **)&tag_space_data); CHKERR_THROW(rval);
  //approx. base
  Tag th_field_base;
  rval = moab.tag_get_handle("_FieldBase",th_field_base); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_field_base,&meshSet,1,(const void **)&tag_base_data); CHKERR_THROW(rval);
  //name
  Tag th_field_name;
  rval = moab.tag_get_handle("_FieldName",th_field_name); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_field_name,&meshSet,1,(const void **)&tag_name_data,&tag_name_size); CHKERR_THROW(rval);
  //name prefix
  Tag th_field_name_data_name_prefix;
  rval = moab.tag_get_handle("_FieldName_DataNamePrefix",th_field_name_data_name_prefix); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_field_name_data_name_prefix,&meshSet,1,(const void **)&tag_name_prefix_data,&tag_name_prefix_size); CHKERR_THROW(rval);
  string name_data_prefix((char *)tag_name_prefix_data,tag_name_prefix_size);
  //data
  string tag_data_name = name_data_prefix+get_name();
  rval = moab.tag_get_handle(tag_data_name.c_str(),th_FieldData); CHKERR_THROW(rval);
  //order
  string tag_approximation_order_name = "_App_Order_"+get_name();
  rval = moab.tag_get_handle(tag_approximation_order_name.c_str(),th_AppOrder); CHKERR_THROW(rval);
  //dof order
  string tag_dof_approximation_order_name = "_App_Dof_Order"+get_name();
  rval = moab.tag_get_handle(tag_dof_approximation_order_name.c_str(),th_AppDofOrder); CHKERR_THROW(rval);
  //rank
  Tag th_rank;
  string Tag_rank_name = "_Field_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_rank_name.c_str(),th_rank); CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_rank,&meshSet,1,(const void **)&tag_nb_coeff_data); CHKERR_THROW(rval);
  //dof rank
  string Tag_dof_rank_name = "_Field_Dof_Rank_"+get_name();
  rval = moab.tag_get_handle(Tag_dof_rank_name.c_str(),th_DofRank); CHKERR_THROW(rval);
  for(int tt = 0;tt<MBMAXTYPE;tt++) {
    forder_table[tt] = NULL;
  }
  switch (*tag_space_data) {
    case H1:
      forder_table[MBVERTEX] = fNBVERTEX_H1_AINSWORTH_COLE;
      forder_table[MBEDGE] = fNBEDGE_H1_AINSWORTH_COLE;
      forder_table[MBTRI] = fNBFACETRI_H1_AINSWORTH_COLE;
      forder_table[MBQUAD] = fNBFACEQUAD_H1_AINSWORTH_COLE;
      forder_table[MBTET] = fNBVOLUMETET_H1_AINSWORTH_COLE;
      forder_table[MBPRISM] = fNBVOLUMEPRISM_H1_AINSWORTH_COLE;
      break;
    case HDIV:
      forder_table[MBVERTEX] = fNBVERTEX_HDIV_AINSWORTH_COLE;
      forder_table[MBEDGE] = fNBEDGE_HDIV_AINSWORTH_COLE;
      forder_table[MBTRI] = fNBFACETRI_HDIV_AINSWORTH_COLE;
      forder_table[MBTET] = fNBVOLUMETET_HDIV_AINSWORTH_COLE;
      break;
    case HCURL:
      forder_table[MBVERTEX] = fNBVERTEX_HCURL_AINSWORTH_COLE;
      forder_table[MBEDGE] = fNBEDGE_HCURL_AINSWORTH_COLE;
      forder_table[MBTRI] = fNBFACETRI_HCURL_AINSWORTH_COLE;
      forder_table[MBTET] = fNBVOLUMETET_HCURL_AINSWORTH_COLE;
      break;
    case L2:
      forder_table[MBVERTEX] = fNBVERTEX_L2;
      forder_table[MBEDGE] = fNBEDGE_L2_AINSWORTH_COLE;
      forder_table[MBTRI] = fNBFACETRI_L2_AINSWORTH_COLE;
      forder_table[MBTET] = fNBVOLUMETET_L2_AINSWORTH_COLE;
      break;
    case NOFIELD:
    for(EntityType t = MBVERTEX;t<MBMAXTYPE;t++) {
      // Concept of approximation order make no sense is there is no field
      forder_table[t] = fNBENTITYSET_NOFIELD;
    }
      break;
    default:
      THROW_AT_LINE("not implemented");
  }
}

ostream& operator<<(ostream& os,const MoFEMField& e) {
  os << "name "<<e.get_name_ref()<<" BitFieldId "<< e.get_id().to_ulong() << " bit number " << e.get_bit_number()
    << " space " << FieldSpaceNames[e.get_space()] << " rank " << e.get_nb_of_coeffs() << " meshset " << e.meshSet;
  return os;
}

//MoFEMEntityEntMoFEMFiniteElementAdjacencyMap
MoFEMEntityEntMoFEMFiniteElementAdjacencyMap::MoFEMEntityEntMoFEMFiniteElementAdjacencyMap(
  const MoFEMEntity *_MoFEMEntity_ptr,const EntMoFEMFiniteElement *_EntMoFEMFiniteElement_ptr
):
by_other(0),
MoFEMEntity_ptr(_MoFEMEntity_ptr),
EntMoFEMFiniteElement_ptr(_EntMoFEMFiniteElement_ptr) {}

ostream& operator<<(ostream& os,const MoFEMEntityEntMoFEMFiniteElementAdjacencyMap& e) {
  os << "by_other " << bitset<3>(e.by_other) << " "
    << *e.MoFEMEntity_ptr << endl << *e.EntMoFEMFiniteElement_ptr->fe_ptr;
  return os;
}

//....
PetscErrorCode test_moab(Interface &moab,const EntityHandle ent) {
  PetscFunctionBegin;
  //tets type
  EntityType type = (EntityType)((ent&MB_TYPE_MASK)>>MB_ID_WIDTH);
  if(type != moab.type_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistencies with type_from_handle");
  //tets id
  EntityID id = (EntityType)(ent&MB_ID_MASK);
  if(id != moab.id_from_handle(ent)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"inconsistencies with id_from_handle");
  PetscFunctionReturn(0);
}

}
