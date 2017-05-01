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
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

// Not partitioned
const bool Idx_mi_tag::IamNotPartitioned = true;

// This tag is used for partitioned problems
const bool PetscGlobalIdx_mi_tag::IamNotPartitioned = false;
const bool PetscLocalIdx_mi_tag::IamNotPartitioned = false;

//fields
Field::Field(
  const Interface &moab,
  const EntityHandle meshset,
  const boost::shared_ptr<CoordSys> coord_sys_ptr
):
moab(const_cast<Interface&>(moab)),
meshSet(meshset),
coordSysPtr(coord_sys_ptr),
tag_id_data(NULL),
tag_space_data(NULL),
tag_nb_coeff_data(NULL),
tag_name_data(NULL),
tag_name_size(0),
sequenceEntContainer(
  boost::make_shared<SequenceEntContainer>()
),
sequenceDofContainer(
  boost::make_shared<SequenceDofContainer>()
) {
  //Change those tags only by modifiers
  ErrorCode rval;
  //id
  Tag th_field_id;
  rval = moab.tag_get_handle("_FieldId",th_field_id); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_field_id,&meshSet,1,(const void **)&tag_id_data); MOAB_THROW(rval);
  //space
  Tag th_field_space;
  rval = moab.tag_get_handle("_FieldSpace",th_field_space); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_field_space,&meshSet,1,(const void **)&tag_space_data); MOAB_THROW(rval);
  //approx. base
  Tag th_field_base;
  rval = moab.tag_get_handle("_FieldBase",th_field_base); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_field_base,&meshSet,1,(const void **)&tag_base_data); MOAB_THROW(rval);
  //name
  Tag th_field_name;
  rval = moab.tag_get_handle("_FieldName",th_field_name); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_field_name,&meshSet,1,(const void **)&tag_name_data,&tag_name_size); MOAB_THROW(rval);
  //name prefix
  Tag th_field_name_data_name_prefix;
  rval = moab.tag_get_handle("_FieldName_DataNamePrefix",th_field_name_data_name_prefix); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_field_name_data_name_prefix,&meshSet,1,(const void **)&tag_name_prefix_data,&tag_name_prefix_size); MOAB_THROW(rval);
  std::string name_data_prefix((char *)tag_name_prefix_data,tag_name_prefix_size);
  //data
  std::string tag_data_name = name_data_prefix+getName();
  rval = moab.tag_get_handle(tag_data_name.c_str(),th_FieldData); MOAB_THROW(rval);
  //order
  std::string tag_approximation_order_name = "_App_Order_"+getName();
  rval = moab.tag_get_handle(tag_approximation_order_name.c_str(),th_AppOrder); MOAB_THROW(rval);
  //rank
  Tag th_rank;
  std::string Tag_rank_name = "_Field_Rank_"+getName();
  rval = moab.tag_get_handle(Tag_rank_name.c_str(),th_rank); MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_rank,&meshSet,1,(const void **)&tag_nb_coeff_data); MOAB_THROW(rval);
  bit_number = getBitNumberCalculate();
  for(int tt = 0;tt<MBMAXTYPE;tt++) {
    forder_table[tt] = NULL;
  }
  switch(*tag_base_data) {
    case AINSWORTH_LEGENDRE_BASE:
    case AINSWORTH_LOBBATO_BASE:
    switch(*tag_space_data) {
      case H1:
      forder_table[MBVERTEX] = fNBVERTEX_H1;
      forder_table[MBEDGE] = fNBEDGE_H1;
      forder_table[MBTRI] = fNBFACETRI_H1;
      forder_table[MBQUAD] = fNBFACEQUAD_H1;
      forder_table[MBTET] = fNBVOLUMETET_H1;
      forder_table[MBPRISM] = fNBVOLUMEPRISM_H1;
      break;
      case HDIV:
      forder_table[MBVERTEX] = fNBVERTEX_HDIV;
      forder_table[MBEDGE] = fNBEDGE_HDIV;
      forder_table[MBTRI] = fNBFACETRI_AINSWORTH_HDIV;
      forder_table[MBTET] = fNBVOLUMETET_AINSWORTH_HDIV;
      break;
      case HCURL:
      forder_table[MBVERTEX] = fNBVERTEX_HCURL;
      forder_table[MBEDGE] = fNBEDGE_HCURL;
      forder_table[MBTRI] = fNBFACETRI_HCURL;
      forder_table[MBTET] = fNBVOLUMETET_HCURL;
      break;
      case L2:
      forder_table[MBVERTEX] = fNBVERTEX_L2;
      forder_table[MBEDGE] = fNBEDGE_L2;
      forder_table[MBTRI] = fNBFACETRI_L2;
      forder_table[MBTET] = fNBVOLUMETET_L2;
      break;
      case NOFIELD:
      for(EntityType t = MBVERTEX;t<MBMAXTYPE;t++) {
        // Concept of approximation order make no sense is there is no field
        forder_table[t] = fNBENTITYSET_NOFIELD;
      }
      break;
      default:
      THROW_MESSAGE("unknown approximation space");
    }
    break;
    case DEMKOWICZ_JACOBI_BASE:
    switch(*tag_space_data) {
      case HDIV:
      forder_table[MBVERTEX] = fNBVERTEX_HDIV;
      forder_table[MBEDGE] = fNBEDGE_HDIV;
      forder_table[MBTRI] = fNBFACETRI_DEMKOWICZ_HDIV;
      forder_table[MBTET] = fNBVOLUMETET_DEMKOWICZ_HDIV;
      break;
      default:
      THROW_MESSAGE("unknown approximation space or not yet implemented");
    }
    break;
    case AINSOWRTH_BERNSTEIN_BEZIER_BASE:
      THROW_MESSAGE("AINSOWRTH_BERNSTEIN_BEZIER_BASE not implemented yer")
    break;
    case USER_BASE:
    for(int ee = 0;ee<MBMAXTYPE;ee++) {
      forder_table[ee] = fNBENTITY_GENERIC;
    }
    break;
    default:
    if(*tag_space_data!=NOFIELD) {
      THROW_MESSAGE("unknown approximation base");
    } else {
      for(EntityType t = MBVERTEX;t<MBMAXTYPE;t++) {
        forder_table[t] = fNBENTITYSET_NOFIELD;
      }
    }
  }
  // // Set DOFs orders on entities
  // for(EntityType ee = MBVERTEX;ee!=MBMAXTYPE;ee++) {
  //   getDofOrderMap(ee).resize(MAX_DOFS_ON_ENTITY,-1);
  // }
}

std::ostream& operator<<(std::ostream& os,const Field& e) {
  os
  << "name " <<e.getNameRef()
  << " BitFieldId "<< e.getId().to_ulong()
  << " bit number " << e.getBitNumber()
  << " space " << FieldSpaceNames[e.getSpace()]
  << " approximation base " << ApproximationBaseNames[e.getApproxBase()]
  << " rank " << e.getNbOfCoeffs()
  << " meshset " << e.meshSet;
  return os;
}

//FieldEntityEntFiniteElementAdjacencyMap
FieldEntityEntFiniteElementAdjacencyMap::FieldEntityEntFiniteElementAdjacencyMap(
  const boost::shared_ptr<FieldEntity> mofem_ent_ptr,
  const boost::shared_ptr<EntFiniteElement> ent_fe_ptr
):
by_other(0),
mofemEntPtr(mofem_ent_ptr),
entFePtr(ent_fe_ptr) {}

std::ostream& operator<<(std::ostream& os,const FieldEntityEntFiniteElementAdjacencyMap& e) {
  os << "by_other " << std::bitset<3>(e.by_other) << " "
    << *e.mofemEntPtr << std::endl << *e.entFePtr->sFePtr;
  return os;
}

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
