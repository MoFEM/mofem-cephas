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

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode SeriesRecorder::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_MOFEMSeriesRecorder) {
    *iface = dynamic_cast<SeriesRecorder*>(this);
    MoFEMFunctionReturnHot(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<UnknownInterface*>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  MoFEMFunctionReturnHot(0);
}

SeriesRecorder::SeriesRecorder(const MoFEM::Core &core):
cOre(const_cast<MoFEM::Core&>(core)) {
}

PetscErrorCode SeriesRecorder::getTags(int verb) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  const int def_val_len = 0;
  rval = moab.tag_get_handle("_SeriesName",def_val_len,MB_TYPE_OPAQUE,
    th_SeriesName,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,NULL
  ); CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::clearMap() {
  MoFEMFunctionBeginHot;
  sEries.clear();
  seriesSteps.clear();
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::initialiseDatabseInformationFromMesh(int verb) {


  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Range meshsets;
  const int def_val_len = 0;
  rval = moab.tag_get_handle("_SeriesName",def_val_len,MB_TYPE_OPAQUE,
    th_SeriesName,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,NULL
  ); CHKERRQ_MOAB(rval);

  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERRQ_MOAB(rval);
  for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
    const void* tagName;
    int tagNameSize;
    rval = moab.tag_get_by_ptr(
      th_SeriesName,&*mit,1,(const void **)&tagName,&tagNameSize
    );
    if(rval == MB_SUCCESS) {
      std::pair<Series_multiIndex::iterator,bool> p = sEries.insert(FieldSeries(moab,*mit));
      if(verb > 0) {
        std::ostringstream ss;
        ss << "read series " << *p.first << std::endl;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
    }
  }
  // //build series steps
  for(Series_multiIndex::iterator sit = sEries.begin();sit!=sEries.end();sit++) {
    int nb_steps;
    ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
    int ss = 0;
    for(;ss<nb_steps;ss++) {
      std::pair<SeriesStep_multiIndex::iterator,bool> p = seriesSteps.insert(FieldSeriesStep(moab,&*sit,ss));
      if(verb > 0) {
        std::ostringstream ss;
        ss << "add series step " << *p.first << std::endl;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::add_series_recorder(const std::string &series_name) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  void const* tag_data[] = { series_name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = series_name.size();
  rval = moab.tag_set_by_ptr(th_SeriesName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
  std::pair<Series_multiIndex::iterator,bool> p = sEries.insert(FieldSeries(moab,meshset));
  if(!p.second) {
    SETERRQ1(PETSC_COMM_SELF,1,"series recorder <%s> is already there",series_name.c_str());
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::delete_recorder_series(const std::string& series_name) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit,hi_ssit;
  ssit = seriesSteps.get<SeriesName_mi_tag>().lower_bound(series_name);
  hi_ssit = seriesSteps.get<SeriesName_mi_tag>().upper_bound(series_name);
  seriesSteps.get<SeriesName_mi_tag>().erase(ssit,hi_ssit);
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit;
  sit = sEries.get<SeriesName_mi_tag>().find(series_name);
  if(sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(
      m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,
      "series recorder <%s> not exist and can be deleted",series_name.c_str()
    );
  }
  EntityHandle series_meshset = sit->getMeshset();
  rval = moab.tag_delete(sit->th_SeriesTime); CHKERRQ_MOAB(rval);
  rval = moab.tag_delete(sit->th_SeriesDataHandles); CHKERRQ_MOAB(rval);
  rval = moab.tag_delete(sit->th_SeriesDataUIDs); CHKERRQ_MOAB(rval);
  rval = moab.tag_delete(sit->th_SeriesData); CHKERRQ_MOAB(rval);
  sEries.get<SeriesName_mi_tag>().erase(sit);
  std::vector<EntityHandle> contained;
  rval = moab.get_contained_meshsets(series_meshset,contained); CHKERRQ_MOAB(rval);
  rval = moab.remove_entities(series_meshset,&contained[0],contained.size()); CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&contained[0],contained.size()); CHKERRQ_MOAB(rval);
  rval = moab.delete_entities(&series_meshset,1); CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::record_problem(const std::string& serie_name,const Problem *problemPtr,RowColData rc) {
  MoFEMFunctionBeginHot;

  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"series recorder <%s> not exist",serie_name.c_str());
  }
  switch (rc) {
    case ROW:
      ierr = const_cast<FieldSeries*>(&*sit)->push_dofs(
        problemPtr->numeredDofsRows->begin(),problemPtr->numeredDofsRows->end()
      ); CHKERRQ(ierr);
      break;
    case COL:
      ierr = const_cast<FieldSeries*>(&*sit)->push_dofs(
        problemPtr->numeredDofsCols->begin(),problemPtr->numeredDofsCols->end()
      ); CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::record_problem(const std::string& serie_name,const std::string& problem_name,RowColData rc) {

  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  const Problem *problem_ptr;
  ierr = m_field.get_problem(problem_name,&problem_ptr); CHKERRQ(ierr);
  ierr = record_problem(serie_name,problem_ptr,rc); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::record_field(
  const std::string& serie_name,const std::string& field_name,
  const BitRefLevel &bit,const BitRefLevel &mask
) {

  MoFEM::Interface &m_field = cOre;
  const DofEntity_multiIndex *dofs_ptr;
  MoFEMFunctionBeginHot;
  ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  DofEntityByFieldName::iterator dit = dofs_ptr->get<FieldName_mi_tag>().lower_bound(field_name);
  if(dit == dofs_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field <%s> not exist",field_name.c_str());
  }
  DofEntityByFieldName::iterator hi_dit = dofs_ptr->get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    const BitRefLevel &dof_bit = (*dit)->getBitRefLevel();
    if((dof_bit&mask) != dof_bit) continue;
    if((dof_bit&bit).any()) {
      EntityHandle ent = (*dit)->getEnt();
      ShortId uid = (*dit)->getNonNonuniqueShortId();
      FieldData val = (*dit)->getFieldData();
      ierr = const_cast<FieldSeries*>(&*sit)->push_dofs(ent,uid,val); CHKERRQ(ierr);
    }
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::record_begin(const std::string& serie_name) {

  MoFEMFunctionBeginHot;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<FieldSeries*>(&*sit)->begin(); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::record_end(const std::string& serie_name,double time) {

  MoFEMFunctionBeginHot;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<FieldSeries*>(&*sit)->end(time); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::initialize_series_recorder(const std::string& serie_name) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<FieldSeries*>(&*sit)->read(moab); CHKERRQ(ierr);
  const_cast<FieldSeries*>(&*sit)->record_begin = false;
  const_cast<FieldSeries*>(&*sit)->record_end = false;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::finalize_series_recorder(const std::string& serie_name) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = sit->save(moab); CHKERRQ(ierr);
  int nb_steps;
  ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
  int ss = 0;
  for(;ss<nb_steps;ss++) {
    /*std::pair<SeriesStep_multiIndex::iterator,bool> p =*/ seriesSteps.insert(FieldSeriesStep(moab,&*sit,ss));
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode SeriesRecorder::print_series_steps() {
  //
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBeginHot;
  std::ostringstream ss;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().begin();
  for(;sit!=sEries.get<SeriesName_mi_tag>().end();sit++) {
    ss << "series " << *sit << std::endl;
  }
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit = seriesSteps.get<SeriesName_mi_tag>().begin();
  for(;ssit!=seriesSteps.get<SeriesName_mi_tag>().end();ssit++) {
    ss << "serises steps " << *ssit << std::endl;
  }
  PetscPrintf(m_field.get_comm(),ss.str().c_str());
  MoFEMFunctionReturnHot(0);
}

bool SeriesRecorder::check_series(const std::string& name) const {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(name);
  if(sit!=sEries.get<SeriesName_mi_tag>().end()) return true;
  return false;
}

PetscErrorCode SeriesRecorder::load_series_data(const std::string& serie_name,const int step_number) {

  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  const DofEntity_multiIndex *dofs_ptr;
  MoFEMFunctionBeginHot;
  SeriesStep_multiIndex::index<Composite_SeriesName_And_Step_mi_tag>::type::iterator sit;
  sit = seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().find(boost::make_tuple(serie_name,step_number));
  if(sit == seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().end()) {
    SETERRQ2(PETSC_COMM_SELF,1,"series <%s> and step %d not found",serie_name.c_str(),step_number);
  }
  ierr = m_field.get_dofs(&dofs_ptr); CHKERRQ(ierr);
  ierr = sit->get(moab,*(const_cast<DofEntity_multiIndex*>(dofs_ptr))); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator SeriesRecorder::get_series_steps_byName_begin(const std::string& name) {
  return seriesSteps.get<SeriesName_mi_tag>().lower_bound(name);
}

SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator SeriesRecorder::get_series_steps_byName_end(const std::string& name) {
  return seriesSteps.get<SeriesName_mi_tag>().upper_bound(name);
}

}
