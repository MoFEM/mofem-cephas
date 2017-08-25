/** \file FieldBlas.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_section_manager
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
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

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
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

#include <FieldBlas.hpp>

namespace MoFEM {

  PetscErrorCode FieldBlas::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMFieldBlas) {
      *iface = dynamic_cast<FieldBlas*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  FieldBlas::FieldBlas(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)),
  dEbug(false) {
  }
  FieldBlas::~FieldBlas() {
  }

  PetscErrorCode FieldBlas::fieldAxpy(
    const double alpha,const std::string& field_name_x,const std::string& field_name_y,
    bool error_if_missing,bool creat_if_missing
  ) {
    const MoFEM::Interface &m_field = cOre;
    const Field_multiIndex *fields_ptr;
    const FieldEntity_multiIndex *field_ents;
    const DofEntity_multiIndex *dofs_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_fields(&fields_ptr); CHKERRQ(ierr);
    ierr = m_field.get_field_ents(&field_ents); CHKERRQ(ierr);
    ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
    Field_multiIndex::index<FieldName_mi_tag>::type::iterator x_fit = fields_ptr->get<FieldName_mi_tag>().find(field_name_x);
    if(x_fit==fields_ptr->get<FieldName_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"x field < %s > not found, (top tip: check spelling)",field_name_x.c_str());
    }
    Field_multiIndex::index<FieldName_mi_tag>::type::iterator y_fit = fields_ptr->get<FieldName_mi_tag>().find(field_name_y);
    if(y_fit==fields_ptr->get<FieldName_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"y field < %s > not found, (top tip: check spelling)",field_name_y.c_str());
    }
    if((*x_fit)->getSpace() != (*y_fit)->getSpace()) {
      SETERRQ2(PETSC_COMM_SELF,1,"space for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
    }
    if((*x_fit)->getNbOfCoeffs() != (*y_fit)->getNbOfCoeffs()) {
      SETERRQ2(PETSC_COMM_SELF,1,"rank for field < %s > and field <%s> are not compatible",field_name_x.c_str(),field_name_y.c_str());
    }
    FieldEntityByFieldName::iterator x_eit;
    x_eit = field_ents->get<FieldName_mi_tag>().lower_bound(field_name_x.c_str());
    for(;x_eit!=field_ents->get<FieldName_mi_tag>().upper_bound(field_name_x.c_str());x_eit++) {
      int nb_dofs_on_x_entity = (*x_eit)->tag_FieldData_size/sizeof(FieldData);
      for(int dd = 0;dd<nb_dofs_on_x_entity;dd++) {
        ApproximationOrder dof_order = (*x_eit)->getDofOrderMap()[dd];
        FieldCoefficientsNumber dof_rank = dd%(*x_eit)->getNbOfCoeffs();
        FieldData data = (*x_eit)->tag_FieldData[dd];
        DofEntity_multiIndex::index<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>::type::iterator dit;
        dit = dofs_ptr->get<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>().find(
          boost::make_tuple(field_name_y.c_str(),(*x_eit)->getEnt(),dof_order,dof_rank)
        );
        if(dit == dofs_ptr->get<Composite_Name_Ent_Order_And_CoeffIdx_mi_tag>().end()) {
          if(creat_if_missing) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"not yet implemented");
          } else {
            if(error_if_missing) {
              std::ostringstream ss;
              ss << "dof on ent " << (*x_eit)->getEnt() << " order " << dof_order << " rank " << dof_rank << " does not exist";
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
            } else {
              continue;
            }
          }
        }
        (*dit)->getFieldData() += alpha*data;
      }
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode FieldBlas::setField(const double val,const EntityType type,const std::string& field_name) {
    const MoFEM::Interface &m_field = cOre;
    const DofEntity_multiIndex *dofs_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
    DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag >::type::iterator dit,hi_dit;
    dit = dofs_ptr->get<Composite_Name_And_Type_mi_tag >().lower_bound(boost::make_tuple(field_name,type));
    hi_dit = dofs_ptr->get<Composite_Name_And_Type_mi_tag >().upper_bound(boost::make_tuple(field_name,type));
    for(;dit!=hi_dit;dit++) {
      (*dit)->getFieldData() = val;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode FieldBlas::setField(const double val,const EntityType type,const Range &ents,const std::string& field_name) {
    const MoFEM::Interface &m_field = cOre;
    const DofEntity_multiIndex *dofs_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
    DofEntity_multiIndex::index<Composite_Name_And_Type_mi_tag >::type::iterator dit,hi_dit;
    dit = dofs_ptr->get<Composite_Name_And_Type_mi_tag >().lower_bound(boost::make_tuple(field_name,type));
    hi_dit = dofs_ptr->get<Composite_Name_And_Type_mi_tag >().upper_bound(boost::make_tuple(field_name,type));
    EntityHandle ent,last = 0;
    bool cont = true;
    for(;dit!=hi_dit;dit++) {
      ent = (*dit)->getEnt();
      if(ent != last) {
        if(ents.find(ent)==ents.end()) {
          cont = true;
        } else {
          cont = false;
        }
        last = ent;
      }
      if(cont) continue;
      (*dit)->getFieldData() = val;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode FieldBlas::fieldScale(const double alpha,const std::string& field_name) {
    const MoFEM::Interface &m_field = cOre;
    const DofEntity_multiIndex *dofs_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
    DofEntityByFieldName::iterator dit,hi_dit;
    dit = dofs_ptr->get<FieldName_mi_tag>().lower_bound(field_name);
    hi_dit = dofs_ptr->get<FieldName_mi_tag>().upper_bound(field_name);
    for(;dit!=hi_dit;dit++) {
      (*dit)->getFieldData() *= alpha;
    }
    PetscFunctionReturn(0);
  }

}
