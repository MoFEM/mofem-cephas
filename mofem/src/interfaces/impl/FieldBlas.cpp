/** \file FieldBlas.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
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

#include <version.h>
#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

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
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::field_axpy(const double alpha,const string& field_name_x,const string& field_name_y,
  bool error_if_missing,bool creat_if_missing) {
  PetscFunctionBegin;
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator x_fit = fIelds.get<FieldName_mi_tag>().find(field_name_x);
  if(x_fit==fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"x field < %s > not found, (top tip: check spelling)",field_name_x.c_str());
  }
  MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator y_fit = fIelds.get<FieldName_mi_tag>().find(field_name_y);
  if(y_fit==fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"y field < %s > not found, (top tip: check spelling)",field_name_y.c_str());
  }
  if(x_fit->get_space() != y_fit->get_space()) {
    SETERRQ2(PETSC_COMM_SELF,1,"space for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
  }
  if(x_fit->get_nb_of_coeffs() != y_fit->get_nb_of_coeffs()) {
    SETERRQ2(PETSC_COMM_SELF,1,"rank for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
  }
  MoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator x_eit;
  x_eit = entsFields.get<FieldName_mi_tag>().lower_bound(field_name_x.c_str());
  for(;x_eit!=entsFields.get<FieldName_mi_tag>().upper_bound(field_name_x.c_str());x_eit++) {
    int nb_dofs_on_x_entity = x_eit->tag_FieldData_size/sizeof(FieldData);
    for(int dd = 0;dd<nb_dofs_on_x_entity;dd++) {
      ApproximationOrder dof_order = x_eit->tag_dof_order_data[dd];
      FieldCoefficientsNumber dof_rank = x_eit->tag_dof_rank_data[dd];
      FieldData data = x_eit->tag_FieldData[dd];
      DofMoFEMEntity_multiIndex::index<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>::type::iterator dit;
      dit = dofsField.get<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>().find(
        boost::make_tuple(field_name_y.c_str(),x_eit->get_ent(),dof_order,dof_rank)
      );
      if(dit == dofsField.get<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>().end()) {
        if(creat_if_missing) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not yet implemented");
        } else {
          if(error_if_missing) {
            ostringstream ss;
            ss << "dof on ent " << x_eit->get_ent() << " order " << dof_order << " rank " << dof_rank << " does not exist";
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
          } else {
            continue;
          }
        }
      }
      (*dit)->get_FieldData() += alpha*data;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field(const double val,const EntityType type,const string& field_name) {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag >::type::iterator dit,hi_dit;
  dit = dofsField.get<Composite_Name_And_Type_mi_tag >().lower_bound(boost::make_tuple(field_name,type));
  hi_dit = dofsField.get<Composite_Name_And_Type_mi_tag >().upper_bound(boost::make_tuple(field_name,type));
  for(;dit!=hi_dit;dit++) {
    (*dit)->get_FieldData() = val;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::set_field(const double val,const EntityType type,const Range &ents,const string& field_name) {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<Composite_Name_And_Type_mi_tag >::type::iterator dit,hi_dit;
  dit = dofsField.get<Composite_Name_And_Type_mi_tag >().lower_bound(boost::make_tuple(field_name,type));
  hi_dit = dofsField.get<Composite_Name_And_Type_mi_tag >().upper_bound(boost::make_tuple(field_name,type));
  EntityHandle ent,last = 0;
  bool cont;
  for(;dit!=hi_dit;dit++) {
    ent = (*dit)->get_ent();
    if(ent != last) {
      if(ents.find(ent)==ents.end()) {
        cont = true;
      } else {
        cont = false;
      }
      last = ent;
    }
    if(cont) continue;
    (*dit)->get_FieldData() = val;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::field_scale(const double alpha,const string& field_name) {
  PetscFunctionBegin;
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit,hi_dit;
  dit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    (*dit)->get_FieldData() *= alpha;
  }
  PetscFunctionReturn(0);
}

}
