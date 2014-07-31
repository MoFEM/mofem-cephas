/** \file CoreDataStructures.cpp
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


#include <petscsys.h>
#include <cblas.h>

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>
#include <CoreDataStructures.hpp>

namespace MoFEM {

//moab dof
DofMoFEMEntity::DofMoFEMEntity(const MoFEMEntity *_MoFEMEntity_ptr,const ApproximationOrder _dof_order,const ApproximationRank _dof_rank,const DofIdx _dof): 
    interface_MoFEMEntity<MoFEMEntity>(_MoFEMEntity_ptr), dof(_dof),active(false) {
  if(field_ptr->tag_dof_order_data==NULL) {
    ostringstream ss;
    ss << "at " << __LINE__ << " in " << __FILE__;
    ss << " field_ptr->tag_dof_order_data==NULL";
    ss << " (top tip: check if order set to vertices is 1)";
    throw(ss.str().c_str());
  }
  assert(field_ptr->tag_dof_rank_data!=NULL);
  ((ApproximationOrder*)field_ptr->tag_dof_order_data)[dof] = _dof_order;
  ((ApproximationRank*)field_ptr->tag_dof_rank_data)[dof] = _dof_rank;
  local_uid = get_local_unique_id_calculate();
  global_uid = get_global_unique_id_calculate();
  short_uid = get_non_nonunique_short_id_calculate();
}
ostream& operator<<(ostream& os,const DofMoFEMEntity& e) {
  os << "dof_uid " << e.get_global_unique_id()
    << " dof_order " << e.get_dof_order()
    << " dof_rank " << e.get_dof_rank()
    << " dof " << e.dof
    << " active " << e.active
    << " " << *(e.field_ptr); /*
    << " Data " << e.get_FieldData();*/
  return os;
}
DofMoFEMEntity_active_change::DofMoFEMEntity_active_change(bool _active): active(_active) {}
void DofMoFEMEntity_active_change::operator()(DofMoFEMEntity &_dof_) {
  _dof_.active = active;
  assert((_dof_.get_dof_order()<=_dof_.get_max_order()));
}

//numered dof
NumeredDofMoFEMEntity::NumeredDofMoFEMEntity(const DofMoFEMEntity* _DofMoFEMEntity_ptr): 
    interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr),
    dof_idx(-1),petsc_gloabl_dof_idx(-1),petsc_local_dof_idx(-1),part(-1) {}
ostream& operator<<(ostream& os,const NumeredDofMoFEMEntity& e) {
  os << "idx " << e.dof_idx << " part " << e.part 
    << " petsc idx " << e.petsc_gloabl_dof_idx 
    << " ( " << e.petsc_local_dof_idx <<  " ) "
    << *e.field_ptr;
  return os;
}

FEDofMoFEMEntity::FEDofMoFEMEntity(boost::tuple<SideNumber *,const DofMoFEMEntity *> t):
  BaseFEDofMoFEMEntity(t.get<0>()), interface_DofMoFEMEntity<DofMoFEMEntity>(t.get<1>()) {}


FEDofMoFEMEntity::FEDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const DofMoFEMEntity *_DofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_DofMoFEMEntity<DofMoFEMEntity>(_DofMoFEMEntity_ptr) {}
ostream& operator<<(ostream& os,const FEDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

FENumeredDofMoFEMEntity::FENumeredDofMoFEMEntity(
  SideNumber *_side_number_ptr,
  const NumeredDofMoFEMEntity *_NumeredDofMoFEMEntity_ptr): 
  BaseFEDofMoFEMEntity(_side_number_ptr), interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity>(_NumeredDofMoFEMEntity_ptr) {}
FENumeredDofMoFEMEntity::FENumeredDofMoFEMEntity(
  boost::tuple<SideNumber *,const NumeredDofMoFEMEntity *> t): 
  BaseFEDofMoFEMEntity(t.get<0>()), interface_NumeredDofMoFEMEntity<NumeredDofMoFEMEntity>(t.get<1>()) {}

ostream& operator<<(ostream& os,const FENumeredDofMoFEMEntity& e) {
  os << "local dof MoFEMFiniteElement idx " 
    << "side_number " << e.side_number_ptr->side_number << " "
    << "sense " << e.side_number_ptr->sense << " "
    << *e.field_ptr;
  return os;
}

}
