/** \file SeriesRecorder.cpp
 * \brief Record and save data series in time or load steps
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

#define SeriesRecorderFunctionBegin                                            \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("SYNC", "SeriesRecorder");                                     \
  MOFEM_LOG_TAG("WORLD", "SeriesRecorder")

namespace MoFEM {

// const static int debug = 1;

MoFEMErrorCode SeriesRecorder::query_interface(const MOFEMuuid &uuid,
                                               UnknownInterface **iface) const {
  MoFEMFunctionBegin;
  *iface = NULL;
  if (uuid == IDD_MOFEMSeriesRecorder) {
    *iface = const_cast<SeriesRecorder *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturn(0);
}

SeriesRecorder::SeriesRecorder(const Core &core)
    : cOre(const_cast<Core &>(core)) {}

MoFEMErrorCode SeriesRecorder::getTags(int verb) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  const int def_val_len = 0;
  CHKERR moab.tag_get_handle(
      "_SeriesName", def_val_len, MB_TYPE_OPAQUE, th_SeriesName,
      MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::clearMap() {
  MoFEMFunctionBegin;
  sEries.clear();
  seriesSteps.clear();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::initialiseDatabaseFromMesh(int verb) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  SeriesRecorderFunctionBegin;
  Range meshsets;
  const int def_val_len = 0;
  CHKERR moab.tag_get_handle(
      "_SeriesName", def_val_len, MB_TYPE_OPAQUE, th_SeriesName,
      MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);

  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, false);
  for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
    const void *tagName;
    int tagNameSize;
    rval = moab.tag_get_by_ptr(th_SeriesName, &*mit, 1, (const void **)&tagName,
                               &tagNameSize);
    if (rval == MB_SUCCESS) {
      std::pair<Series_multiIndex::iterator, bool> p =
          sEries.insert(FieldSeries(moab, *mit));
      if (verb > QUIET)
        MOFEM_LOG("SYNC", Sev::inform) << "read series " << *p.first;
    }
  }
  // //build series steps
  for (Series_multiIndex::iterator sit = sEries.begin(); sit != sEries.end();
       sit++) {
    int nb_steps;
    CHKERR sit->get_nb_steps(moab, nb_steps);
    int ss = 0;
    for (; ss < nb_steps; ss++) {
      std::pair<SeriesStep_multiIndex::iterator, bool> p =
          seriesSteps.insert(FieldSeriesStep(moab, &*sit, ss));
      if (verb > QUIET)
        MOFEM_LOG("SYNC", Sev::inform) << "add series step " << *p.first;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
SeriesRecorder::add_series_recorder(const std::string &series_name) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  EntityHandle meshset;
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  void const *tag_data[] = {series_name.c_str()};
  int tag_sizes[1];
  tag_sizes[0] = series_name.size();
  CHKERR moab.tag_set_by_ptr(th_SeriesName, &meshset, 1, tag_data, tag_sizes);
  std::pair<Series_multiIndex::iterator, bool> p =
      sEries.insert(FieldSeries(moab, meshset));
  if (!p.second) {
    SETERRQ1(PETSC_COMM_SELF, 1, "series recorder <%s> is already there",
             series_name.c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
SeriesRecorder::delete_recorder_series(const std::string &series_name) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit, hi_ssit;
  ssit = seriesSteps.get<SeriesName_mi_tag>().lower_bound(series_name);
  hi_ssit = seriesSteps.get<SeriesName_mi_tag>().upper_bound(series_name);
  seriesSteps.get<SeriesName_mi_tag>().erase(ssit, hi_ssit);
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit;
  sit = sEries.get<SeriesName_mi_tag>().find(series_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "series recorder <%s> not exist and can be deleted",
             series_name.c_str());
  }
  EntityHandle series_meshset = sit->getMeshset();
  CHKERR moab.tag_delete(sit->th_SeriesTime);
  CHKERR moab.tag_delete(sit->th_SeriesDataHandles);
  CHKERR moab.tag_delete(sit->th_SeriesDataUIDs);
  CHKERR moab.tag_delete(sit->th_SeriesData);
  sEries.get<SeriesName_mi_tag>().erase(sit);
  std::vector<EntityHandle> contained;
  CHKERR moab.get_contained_meshsets(series_meshset, contained);
  CHKERR moab.remove_entities(series_meshset, &contained[0], contained.size());
  CHKERR moab.delete_entities(&contained[0], contained.size());
  CHKERR moab.delete_entities(&series_meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::record_problem(const std::string &serie_name,
                                              const Problem *problemPtr,
                                              RowColData rc) {
  MoFEMFunctionBegin;

  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "series recorder <%s> not exist", serie_name.c_str());
  }
  switch (rc) {
  case ROW:
    CHKERR const_cast<FieldSeries *>(&*sit)->push_dofs(
        problemPtr->numeredRowDofsPtr->begin(),
        problemPtr->numeredRowDofsPtr->end());
    break;
  case COL:
    CHKERR const_cast<FieldSeries *>(&*sit)->push_dofs(
        problemPtr->numeredColDofsPtr->begin(),
        problemPtr->numeredColDofsPtr->end());
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::record_problem(const std::string &serie_name,
                                              const std::string &problem_name,
                                              RowColData rc) {

  Interface &m_field = cOre;
  MoFEMFunctionBegin;
  const Problem *problem_ptr;
  CHKERR m_field.get_problem(problem_name, &problem_ptr);
  CHKERR record_problem(serie_name, problem_ptr, rc);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::record_field(const std::string &serie_name,
                                            const std::string &field_name,
                                            const BitRefLevel &bit,
                                            const BitRefLevel &mask) {
  Interface &m_field = cOre;
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "serie recorder <%s> not exist",
             serie_name.c_str());
  }

  const auto bit_number = m_field.get_field_bit_number(field_name);

  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId(bit_number));
  if (dit == dofs_ptr->get<Unique_mi_tag>().end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "field <%s> not exist", field_name.c_str());
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(bit_number));

  for (; dit != hi_dit; dit++) {
    const BitRefLevel &dof_bit = (*dit)->getBitRefLevel();
    if ((dof_bit & mask) != dof_bit)
      continue;
    if ((dof_bit & bit).any()) {
      EntityHandle ent = (*dit)->getEnt();
      ShortId uid = (*dit)->getNonNonuniqueShortId();
      FieldData val = (*dit)->getFieldData();
      CHKERR const_cast<FieldSeries *>(&*sit)->push_dofs(ent, uid, val);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::record_begin(const std::string &serie_name) {

  MoFEMFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "serie recorder <%s> not exist",
             serie_name.c_str());
  }
  CHKERR const_cast<FieldSeries *>(&*sit)->begin();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::record_end(const std::string &serie_name,
                                          double time) {

  MoFEMFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "serie recorder <%s> not exist",
             serie_name.c_str());
  }
  CHKERR const_cast<FieldSeries *>(&*sit)->end(time);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
SeriesRecorder::initialize_series_recorder(const std::string &serie_name) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "serie recorder <%s> not exist",
             serie_name.c_str());
  }
  CHKERR const_cast<FieldSeries *>(&*sit)->read(moab);
  const_cast<FieldSeries *>(&*sit)->record_begin = false;
  const_cast<FieldSeries *>(&*sit)->record_end = false;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
SeriesRecorder::finalize_series_recorder(const std::string &serie_name) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(serie_name);
  if (sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, 1, "serie recorder <%s> not exist",
             serie_name.c_str());
  }
  CHKERR sit->save(moab);
  int nb_steps;
  CHKERR sit->get_nb_steps(moab, nb_steps);
  int ss = 0;
  for (; ss < nb_steps; ss++) {
    /*std::pair<SeriesStep_multiIndex::iterator,bool> p =*/seriesSteps.insert(
        FieldSeriesStep(moab, &*sit, ss));
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode SeriesRecorder::print_series_steps() {
  SeriesRecorderFunctionBegin;

  for (auto &sit : sEries.get<SeriesName_mi_tag>())
    MOFEM_LOG("SYNC", Sev::inform) << "series " << sit;

  for (auto &ssit : seriesSteps.get<SeriesName_mi_tag>())
    MOFEM_LOG("SYNC", Sev::inform) << "series steps " << ssit;

  MoFEMFunctionReturn(0);
}

bool SeriesRecorder::check_series(const std::string &name) const {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit =
      sEries.get<SeriesName_mi_tag>().find(name);
  if (sit != sEries.get<SeriesName_mi_tag>().end())
    return true;
  return false;
}

MoFEMErrorCode SeriesRecorder::load_series_data(const std::string &serie_name,
                                                const int step_number) {

  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;
  SeriesStep_multiIndex::index<
      Composite_SeriesName_And_Step_mi_tag>::type::iterator sit;
  sit = seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().find(
      boost::make_tuple(serie_name, step_number));
  if (sit == seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().end()) {
    SETERRQ2(PETSC_COMM_SELF, 1, "series <%s> and step %d not found",
             serie_name.c_str(), step_number);
  }
  CHKERR sit->get(moab, *(const_cast<DofEntity_multiIndex *>(dofs_ptr)));
  MoFEMFunctionReturn(0);
}

SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator
SeriesRecorder::get_series_steps_byName_begin(const std::string &name) {
  return seriesSteps.get<SeriesName_mi_tag>().lower_bound(name);
}

SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator
SeriesRecorder::get_series_steps_byName_end(const std::string &name) {
  return seriesSteps.get<SeriesName_mi_tag>().upper_bound(name);
}

} // namespace MoFEM
