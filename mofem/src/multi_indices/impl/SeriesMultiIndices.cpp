/** \file SeriesMultiIndices.cpp
 * \brief Data structures for time series and load series storage,
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

namespace MoFEM {

FieldSeries::FieldSeries(moab::Interface &moab, const EntityHandle _meshset)
    : meshset(_meshset), tagName(NULL), tagNameSize(0), record_begin(false),
      record_end(false) {

  Tag th_SeriesName;
  rval = moab.tag_get_handle("_SeriesName", th_SeriesName);
  MOAB_THROW(rval);
  rval = moab.tag_get_by_ptr(th_SeriesName, &meshset, 1,
                             (const void **)&tagName, &tagNameSize);
  MOAB_THROW(rval);

  const int def_val_len = 0;

  // time
  std::string Tag_SeriesTime = "_SeriesTime_" + getName();
  double def_time = 0;
  rval = moab.tag_get_handle(Tag_SeriesTime.c_str(), 1, MB_TYPE_DOUBLE,
                             th_SeriesTime, MB_TAG_CREAT | MB_TAG_SPARSE,
                             &def_time);
  MOAB_THROW(rval);

  // handles
  std::string Tag_DataHandles_SeriesName = "_SeriesDataHandles_" + getName();
  rval = moab.tag_get_handle(
      Tag_DataHandles_SeriesName.c_str(), def_val_len, MB_TYPE_HANDLE,
      th_SeriesDataHandles, MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_VARLEN, NULL);
  MOAB_THROW(rval);

  // uids
  std::string Tag_DataUIDs_SeriesName = "_SeriesDataUIDs_" + getName();
  rval = moab.tag_get_handle(
      Tag_DataUIDs_SeriesName.c_str(), def_val_len, MB_TYPE_OPAQUE,
      th_SeriesDataUIDs,
      MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES | MB_TAG_VARLEN, NULL);
  MOAB_THROW(rval);

  // data
  std::string Tag_Data_SeriesName = "_SeriesData_" + getName();
  rval = moab.tag_get_handle(
      Tag_Data_SeriesName.c_str(), def_val_len, MB_TYPE_OPAQUE, th_SeriesData,
      MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES | MB_TAG_VARLEN, NULL);
  MOAB_THROW(rval);
}

MoFEMErrorCode FieldSeries::get_nb_steps(moab::Interface &moab,
                                         int &nb_steps) const {
  MoFEMFunctionBegin;
  CHKERR moab.num_contained_meshsets(meshset, &nb_steps);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeries::push_dofs(const EntityHandle ent, const ShortId uid,
                                      const FieldData val) {
  MoFEMFunctionBegin;
  if (!record_begin) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "you neet to set recording");
  }
  handles.push_back(ent);
  uids.push_back(uid);
  data.push_back(val);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeries::begin() {
  MoFEMFunctionBegin;
  if (record_begin) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "recording already  begin");
  }
  record_begin = true;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeries::end(double t) {
  MoFEMFunctionBegin;
  if (!record_begin) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "recording not begin it can not be ended");
  }
  record_begin = false;
  record_end = true;
  ia.push_back(uids.size());
  time.push_back(t);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeries::read(moab::Interface &moab) {
  MoFEMFunctionBegin;

  if (record_begin) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "all series data will be lost");
  }
  if (record_end) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "all series data will be lost");
  }

  std::vector<EntityHandle> contained;
  CHKERR moab.get_contained_meshsets(meshset, contained);
  ia.resize(0);
  time.resize(0);
  handles.resize(0);
  uids.resize(0);
  data.resize(0);
  ia.push_back(0);

  for (unsigned int mm = 0; mm < contained.size(); mm++) {
    // time
    {
      double t;
      CHKERR moab.tag_set_data(th_SeriesTime, &meshset, 1, &t);
      time.push_back(t);
    }
    // handles
    {
      const EntityHandle *tag_data;
      int tag_size;
      CHKERR moab.tag_get_by_ptr(th_SeriesDataHandles, &meshset, 1,
                                 (const void **)&tag_data, &tag_size);
      handles.insert(handles.end(), tag_data, &tag_data[tag_size]);
    }
    // uids
    {
      const ShortId *tag_data;
      int tag_size;
      CHKERR moab.tag_get_by_ptr(th_SeriesDataUIDs, &meshset, 1,
                                 (const void **)&tag_data, &tag_size);
      int nb = tag_size / sizeof(ShortId);
      uids.insert(uids.end(), tag_data, &tag_data[nb]);
    }
    if (handles.size() != uids.size()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    // data
    {
      const FieldData *tag_data;
      int tag_size;
      CHKERR moab.tag_get_by_ptr(th_SeriesData, &meshset, 1,
                                 (const void **)&tag_data, &tag_size);
      int nb = tag_size / sizeof(FieldData);
      data.insert(data.end(), tag_data, &tag_data[nb]);
    }
    if (data.size() != uids.size()) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
    }
    ia.push_back(data.size());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeries::save(moab::Interface &moab) const {
  MoFEMFunctionBegin;

  if (record_begin) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "switch off recording");
  }
  if (!record_end) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "finish recording");
  }

  std::vector<EntityHandle> contained;
  CHKERR moab.get_contained_meshsets(meshset, contained);
  unsigned int nb_contained = contained.size();
  if (nb_contained < ia.size() - 1) {
    contained.resize(ia.size());
  }
  for (unsigned int mm = ia.size() - 1; mm < nb_contained; mm++) {
    CHKERR moab.remove_entities(meshset, &contained[mm], 1);
    CHKERR moab.delete_entities(&contained[mm], 1);
  }
  for (unsigned int mm = nb_contained; mm < ia.size() - 1; mm++) {
    EntityHandle new_meshset;
    CHKERR moab.create_meshset(MESHSET_SET, new_meshset);
    CHKERR moab.add_entities(meshset, &new_meshset, 1);
  }
  contained.resize(0);
  CHKERR moab.get_contained_meshsets(meshset, contained);
  if (contained.size() != ia.size() - 1) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency nb_contained != ia.size()-1 %d!=%d",
             contained.size(), ia.size() - 1);
  }

  // time
  for (unsigned int ii = 1; ii < ia.size(); ii++) {
    CHKERR moab.tag_set_data(th_SeriesTime, &contained[ii - 1], 1,
                             &time[ii - 1]);
  }

  // handles
  for (unsigned int ii = 1; ii < ia.size(); ii++) {
    void const *tag_data[] = {&handles[ia[ii - 1]]};
    int tag_sizes[] = {(ia[ii] - ia[ii - 1])};
    CHKERR moab.tag_set_by_ptr(th_SeriesDataHandles, &contained[ii - 1], 1,
                               tag_data, tag_sizes);
  }
  // uids
  for (unsigned int ii = 1; ii < ia.size(); ii++) {
    void const *tag_data[] = {&uids[ia[ii - 1]]};
    int tag_sizes[] = {(ia[ii] - ia[ii - 1]) * (int)sizeof(ShortId)};
    CHKERR moab.tag_set_by_ptr(th_SeriesDataUIDs, &contained[ii - 1], 1,
                               tag_data, tag_sizes);
  }

  // data
  for (unsigned int ii = 1; ii < ia.size(); ii++) {
    void const *tag_data[] = {&data[ia[ii - 1]]};
    int tag_sizes[] = {(ia[ii] - ia[ii - 1]) * (int)sizeof(FieldData)};
    CHKERR moab.tag_set_by_ptr(th_SeriesData, &contained[ii - 1], 1, tag_data,
                               tag_sizes);
  }

  MoFEMFunctionReturn(0);
}

FieldSeriesStep::FieldSeriesStep(moab::Interface &moab,
                                 const FieldSeries *_FieldSeries_ptr,
                                 const int _step_number)
    : interface_FieldSeries<FieldSeries>(_FieldSeries_ptr),
      step_number(_step_number) {

  ierr = get_time_init(moab);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

MoFEMErrorCode FieldSeriesStep::get(moab::Interface &moab,
                                    DofEntity_multiIndex &dofsField) const {
  MoFEMFunctionBegin;

  std::vector<EntityHandle> contained;
  CHKERR moab.get_contained_meshsets(ptr->meshset, contained);
  if (contained.size() <= (unsigned int)step_number) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }

  EntityHandle *handles_ptr;
  int handles_size;
  CHKERR moab.tag_get_by_ptr(ptr->th_SeriesDataHandles, &contained[step_number],
                             1, (const void **)&handles_ptr, &handles_size);

  ShortId *uids_ptr;
  int uids_size;
  CHKERR moab.tag_get_by_ptr(ptr->th_SeriesDataUIDs, &contained[step_number], 1,
                             (const void **)&uids_ptr, &uids_size);
  uids_size /= sizeof(ShortId);

  if (handles_size != uids_size) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }

  FieldData *data_ptr;
  int data_size;
  CHKERR moab.tag_get_by_ptr(ptr->th_SeriesData, &contained[step_number], 1,
                             (const void **)&data_ptr, &data_size);
  data_size /= sizeof(FieldData);

  if (data_size != uids_size) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }

  typedef multi_index_container<
      boost::shared_ptr<DofEntity>,
      indexed_by<

          hashed_non_unique<
              tag<Composite_Ent_and_ShortId_mi_tag>,
              composite_key<
                  DofEntity,
                  const_mem_fun<DofEntity, EntityHandle, &DofEntity::getEnt>,
                  const_mem_fun<DofEntity, ShortId,
                                &DofEntity::getNonNonuniqueShortId>>>

          >>
      DofEntity_multiIndex_short_uid_view;

  DofEntity_multiIndex_short_uid_view short_uid_view;
  short_uid_view.insert(dofsField.begin(), dofsField.end());

  for (int ii = 0; ii < uids_size; ii++) {
    EntityHandle ent = handles_ptr[ii];
    ShortId uid = uids_ptr[ii];
    FieldData val = data_ptr[ii];

    auto dit = short_uid_view.find(boost::make_tuple(ent, uid));
    if (dit != short_uid_view.end()) {
      (*dit)->getFieldData() = val;
    } else {
      SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
              "data inconsistency, getting data series, dof on ENTITY and "
              "ShortId can't be found");
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldSeriesStep::get_time_init(moab::Interface &moab) {
  MoFEMFunctionBegin;

  std::vector<EntityHandle> contained;
  CHKERR moab.get_contained_meshsets(ptr->meshset, contained);
  if (contained.size() <= (unsigned int)step_number) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  double *time_ptr;
  int size;
  CHKERR moab.tag_get_by_ptr(ptr->th_SeriesTime, &contained[step_number], 1,
                             (const void **)&time_ptr, &size);
  time = *time_ptr;
  MoFEMFunctionReturn(0);
}

std::ostream &operator<<(std::ostream &os, const FieldSeries &e) {
  os << "name " << e.getName() << " meshset " << e.getMeshset();
  return os;
}

std::ostream &operator<<(std::ostream &os, const FieldSeriesStep &e) {
  os << *(e.get_FieldSeries_ptr()) << " step number " << e.step_number
     << " time " << e.get_time();
  return os;
}

} // namespace MoFEM
