/** \file SeriesMultiIndices.cpp
 * \brief Data strutures for time steries and load series storage,
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

#include <CoreDataStructures.hpp>

namespace MoFEM {

MoFEMSeries::MoFEMSeries(Interface &moab,const EntityHandle _meshset): 
  meshset(_meshset),tag_name_data(NULL),tag_name_size(0),
  record_begin(false),record_end(false) {
  ErrorCode rval;

  Tag th_SeriesName;
  rval = moab.tag_get_handle("_SeriesName",th_SeriesName); CHKERR(rval);
  rval = moab.tag_get_by_ptr(th_SeriesName,&meshset,1,(const void **)&tag_name_data,&tag_name_size); CHKERR_THROW(rval);
 
  const int def_val_len = 0;

  //uids
  string Tag_DataUIDs_SeriesName = "_SeriesDataUIDs_"+get_name();
  rval = moab.tag_get_handle(Tag_DataUIDs_SeriesName.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_SeriesDataUIDs,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);

  //data
  string Tag_Data_SeriesName = "_SeriesData_"+get_name();
  rval = moab.tag_get_handle(Tag_Data_SeriesName.c_str(),def_val_len,MB_TYPE_OPAQUE,
    th_SeriesData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);

}

PetscErrorCode MoFEMSeries::get_nb_steps(Interface &moab,int &nb_steps) const {
  PetscFunctionBegin;
  ErrorCode rval;
  rval = moab.num_contained_meshsets(meshset,&nb_steps); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEMSeries::push_dofs(UId uid,FieldData val) {
  PetscFunctionBegin;
  if(!record_begin) {
    SETERRQ(PETSC_COMM_SELF,1,"you neet to set recording");
  }
  uids.push_back(uid);
  data.push_back(val);
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEMSeries::begin() { 
  PetscFunctionBegin;
  if(record_begin) {
    SETERRQ(PETSC_COMM_SELF,1,"recording already  begin");
  }
  record_begin = true;
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEMSeries::end() { 
  PetscFunctionBegin;
  if(!record_begin) {
    SETERRQ(PETSC_COMM_SELF,1,"recording not begin it can not be ended");
  }
  record_begin = false;
  record_end = true;
  ia.push_back(uids.size());
  PetscFunctionReturn(0);
}

PetscErrorCode MoFEMSeries::read(Interface &moab) {
  PetscFunctionBegin;
  ErrorCode rval;

  if(record_begin) {
    SETERRQ(PETSC_COMM_SELF,1,"all series data will be lost");
  }
  if(record_end) {
    SETERRQ(PETSC_COMM_SELF,1,"all series data will be lost");
  }

  vector<EntityHandle> contained;
  rval = moab.get_contained_meshsets(meshset,contained); CHKERR_PETSC(rval);
  ia.resize(0);
  uids.resize(0);
  data.resize(0);
  ia.push_back(0);
  for(int mm = 0;mm<contained.size();mm++) {
    //uids
    {
      const UId* tag_data;
      int tag_size;
      rval = moab.tag_get_by_ptr(th_SeriesDataUIDs,&meshset,1,(const void **)&tag_data,&tag_size); CHKERR_PETSC(rval);
      int nb = tag_size/sizeof(UId);
      uids.insert(uids.end(),tag_data,&tag_data[nb]);
    }
    {
      const FieldData* tag_data;
      int tag_size;
      rval = moab.tag_get_by_ptr(th_SeriesData,&meshset,1,(const void **)&tag_data,&tag_size); CHKERR_PETSC(rval);
      int nb = tag_size/sizeof(FieldData);
      data.insert(data.end(),tag_data,&tag_data[nb]);
    }
    if(data.size() != uids.size()) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    ia.push_back(data.size());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode MoFEMSeries::save(Interface &moab) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if(record_begin) {
    SETERRQ(PETSC_COMM_SELF,1,"switch off recording");
  }
  if(!record_end) {
    SETERRQ(PETSC_COMM_SELF,1,"finish recording");
  }

  ErrorCode rval;
  vector<EntityHandle> contained;
  rval = moab.get_contained_meshsets(meshset,contained); CHKERR_PETSC(rval);
  int nb_contained = contained.size();
  if(nb_contained < ia.size()-1) {
    contained.resize(ia.size());
  }
  for(int mm = ia.size()-1;mm<nb_contained;mm++) {
    rval = moab.remove_entities(meshset,&contained[mm],1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&contained[mm],1); CHKERR_PETSC(rval);
  }
  for(int mm = nb_contained;mm<ia.size()-1;mm++) {
    EntityHandle new_meshset;
    rval = moab.create_meshset(MESHSET_SET,new_meshset); CHKERR_PETSC(rval);
    rval = moab.add_entities(meshset,&new_meshset,1); CHKERR_PETSC(rval);
  }
  contained.resize(0);
  rval = moab.get_contained_meshsets(meshset,contained); CHKERR_PETSC(rval);
  if(contained.size() != ia.size()-1) {
    SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency nb_contained != ia.size()-1 %d!=%d",contained.size(),ia.size()-1);
  }

  for(int ii = 1;ii<ia.size();ii++) {
    void const* tag_data[] = { &uids[ia[ii-1]] };
    int tag_sizes[] = { (ia[ii]-ia[ii-1])*sizeof(UId) };
    rval = moab.tag_set_by_ptr(th_SeriesDataUIDs,&contained[ii-1],1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  }
  
  for(int ii = 1;ii<ia.size();ii++) {
    void const* tag_data[] = { &data[ia[ii-1]] };
    int tag_sizes[] = { (ia[ii]-ia[ii-1])*sizeof(FieldData) };
    rval = moab.tag_set_by_ptr(th_SeriesData,&contained[ii-1],1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  }

  PetscFunctionReturn(0);
}

MoFEMSeriesStep::MoFEMSeriesStep(const MoFEMSeries *_MoFEMSeries_ptr,const int _step_number): 
  interface_MoFEMSeries<MoFEMSeries>(_MoFEMSeries_ptr),step_number(_step_number) {}

PetscErrorCode MoFEMSeriesStep::get(Interface &moab,DofMoFEMEntity_multiIndex &dofsMoabField) const {
  PetscFunctionBegin;
  ErrorCode rval;

  vector<EntityHandle> contained;
  rval = moab.get_contained_meshsets(ptr->meshset,contained); CHKERR_PETSC(rval);
  if(contained.size()<=step_number) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }

  UId *uids_ptr;
  int uids_size;
  rval = moab.tag_get_by_ptr(ptr->th_SeriesDataUIDs,&contained[step_number],1,(const void **)&uids_ptr,&uids_size); CHKERR_PETSC(rval);
  uids_size /= sizeof(UId);

  FieldData *data_ptr; 
  int data_size;
  rval = moab.tag_get_by_ptr(ptr->th_SeriesData,&contained[step_number],1,(const void **)&data_ptr,&data_size); CHKERR_PETSC(rval);
  data_size /= sizeof(FieldData);

  if(data_size != uids_size) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }

  for(int ii = 0;ii<uids_size;ii++) {
    if(ii>=uids_size) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    DofMoFEMEntity_multiIndex::index<Unique_mi_tag>::type::iterator dit = dofsMoabField.get<Unique_mi_tag>().find(uids_ptr[ii]);
    if(dit!=dofsMoabField.get<Unique_mi_tag>().end()) {
      dit->get_FieldData() = data_ptr[ii];
    } else {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
  }

  PetscFunctionReturn(0);
}

ostream& operator<<(ostream& os,const MoFEMSeries& e) {
  os << "name " << e.get_name() << " meshset " << e.get_meshset();
  return os;
}

ostream& operator<<(ostream& os,const MoFEMSeriesStep& e) {
  os << *(e.get_MoFEMSeries_ptr()) << " step number " << e.step_number ;
  return os;
}


}
