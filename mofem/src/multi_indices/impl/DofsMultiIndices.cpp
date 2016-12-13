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

namespace MoFEM {

//moab dof
DofEntity::DofEntity(
  const boost::shared_ptr<MoFEMEntity> entity_ptr,
  const ApproximationOrder dof_order,
  const FieldCoefficientsNumber dof_rank,
  const DofIdx dof
):
interface_MoFEMEntity<MoFEMEntity>(entity_ptr),
active(false),
dof(dof) {

  if(!entity_ptr) {
    THROW_MESSAGE("MoFEMEntity pinter not initialized");
  }
  if(!sPtr) {
    THROW_MESSAGE("MoFEMEntity pinter not initialized");
  }
  if(!getMoFEMEntityPtr()) {
    THROW_MESSAGE("MoFEMEntity pinter not initialized");
  }

  if(sFieldPtr->tag_dof_order_data==NULL) {
    std::ostringstream ss;
    ss << "at " << __LINE__ << " in " << __FILE__;
    ss << " sFieldPtr->tag_dof_order_data==NULL";
    ss << " (top tip: check if order set to vertices is 1)";
    //throw(ss.str().c_str());
    PetscTraceBackErrorHandler(
      PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,ss.str().c_str(),PETSC_NULL
    );
    PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,ss.str().c_str(),PETSC_NULL
    );
  }
  assert(sFieldPtr->tag_dof_rank_data!=NULL);
  ((ApproximationOrder*)sFieldPtr->tag_dof_order_data)[dof] = dof_order;
  ((FieldCoefficientsNumber*)sFieldPtr->tag_dof_rank_data)[dof] = dof_rank;
  // short_uid = get_non_nonunique_short_id_calculate(dof);
}

std::ostream& operator<<(std::ostream& os,const DofEntity& e) {
  os << "dof_uid " << e.getGlobalUniqueId()
  << " dof_order " << e.getDofOrder()
  << " dof_rank " << e.getDofCoeffIdx()
  << " dof " << e.getEntDofIdx()
  << " active " << e.active
  << " " << *(e.sFieldPtr);
  return os;
}

DofEntity_active_change::DofEntity_active_change(bool _active): active(_active) {}
void DofEntity_active_change::operator()(boost::shared_ptr<DofEntity> &_dof_) {
  _dof_->active = active;
  if(active && _dof_->getDofOrder()>_dof_->getMaxOrder()) {
    cerr << *_dof_ << endl;
    THROW_MESSAGE("Set DoF active which has order larger than maximal order set to entity");
  }
}

//numered dof
NumeredDofEntity::NumeredDofEntity(
  const boost::shared_ptr<DofEntity> dof_entity_ptr,
  const int dof_idx,
  const int petsc_gloabl_dof_idx,
  const int petsc_local_dof_idx,
  const int part

):
interface_DofEntity<DofEntity>(dof_entity_ptr),
dofIdx(dof_idx),
petscGloablDofIdx(petsc_gloabl_dof_idx),
petscLocalDofIdx(petsc_local_dof_idx),
pArt(part) {
}

std::ostream& operator<<(std::ostream& os,const NumeredDofEntity& e) {
  os << "idx " << e.dofIdx << " part " << e.pArt
  << " petsc idx " << e.petscGloablDofIdx
  << " ( " << e.petscLocalDofIdx <<  " ) "
  << *e.sFieldPtr;
  return os;
}

FEDofEntity::FEDofEntity(
  boost::tuple<boost::shared_ptr<SideNumber>,const boost::shared_ptr<DofEntity> > t
):
BaseFEDofEntity(t.get<0>()),
interface_DofEntity<DofEntity>(t.get<1>()) {
}


FEDofEntity::FEDofEntity(
  boost::shared_ptr<SideNumber> side_number_ptr,
  const boost::shared_ptr<DofEntity> dof_ptr
):
BaseFEDofEntity(side_number_ptr),
interface_DofEntity<DofEntity>(dof_ptr) {
}

std::ostream& operator<<(std::ostream& os,const FEDofEntity& e) {
  os << "local dof FiniteElement idx "
    << "side_number " << (int)e.sideNumberPtr->side_number << " "
    << "sense " << (int)e.sideNumberPtr->sense << " "
    << *e.sFieldPtr;
  return os;
}

FENumeredDofEntity::FENumeredDofEntity(
  boost::shared_ptr<SideNumber> side_number_ptr,
  const boost::shared_ptr<NumeredDofEntity> dof_ptr
):
BaseFEDofEntity(side_number_ptr),
interface_NumeredDofEntity<NumeredDofEntity>(dof_ptr) {
}

FENumeredDofEntity::FENumeredDofEntity(
  boost::tuple<boost::shared_ptr<SideNumber>,const boost::shared_ptr<NumeredDofEntity> > t
):
BaseFEDofEntity(t.get<0>()),
interface_NumeredDofEntity<NumeredDofEntity>(t.get<1>()) {
}

std::ostream& operator<<(std::ostream& os,const FENumeredDofEntity& e) {
  os << "local dof FiniteElement idx "
    << "side_number " << (int)e.sideNumberPtr->side_number << " "
    << "sense " << (int)e.sideNumberPtr->sense << " "
    << *e.sFieldPtr;
  return os;
}



}
