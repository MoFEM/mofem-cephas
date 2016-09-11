/** \file SeriesRecorderCore.cpp
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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
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
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

namespace MoFEM {

// const static int debug = 1;

PetscErrorCode Core::add_series_recorder(const std::string &series_name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  void const* tag_data[] = { series_name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = series_name.size();
  rval = moab.tag_set_by_ptr(th_SeriesName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
  std::pair<Series_multiIndex::iterator,bool> p = sEries.insert(MoFEMSeries(moab,meshset));
  if(!p.second) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> is already there",series_name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_recorder_series(const std::string& series_name) {
  PetscFunctionBegin;
  //PetscErrorCode ierr;
  ErrorCode rval;
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit,hi_ssit;
  ssit = seriesSteps.get<SeriesName_mi_tag>().lower_bound(series_name);
  hi_ssit = seriesSteps.get<SeriesName_mi_tag>().upper_bound(series_name);
  seriesSteps.get<SeriesName_mi_tag>().erase(ssit,hi_ssit);
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit;
  sit = sEries.get<SeriesName_mi_tag>().find(series_name);
  if(sit == sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist and can be deleted",series_name.c_str());
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
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_problem(const std::string& serie_name,const MoFEMProblem *problemPtr,RowColData rc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"series recorder <%s> not exist",serie_name.c_str());
  }
  switch (rc) {
    case ROW:
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(
        problemPtr->numered_dofs_rows->begin(),problemPtr->numered_dofs_rows->end()
      ); CHKERRQ(ierr);
      break;
    case COL:
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(
        problemPtr->numered_dofs_cols->begin(),problemPtr->numered_dofs_cols->end()
      ); CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_problem(const std::string& serie_name,const std::string& problem_name,RowColData rc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const MoFEMProblem *problemPtr;
  ierr = get_problem(problem_name,&problemPtr); CHKERRQ(ierr);
  ierr = record_problem(serie_name,problemPtr,rc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_field(const std::string& serie_name,const std::string& field_name,
  const BitRefLevel &bit,const BitRefLevel &mask) {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  PetscFunctionBegin;
  PetscErrorCode ierr;
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  if(dit == dofsField.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field <%s> not exist",field_name.c_str());
  }
  DofEntity_multiIndex::index<FieldName_mi_tag>::type::iterator hi_dit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    const BitRefLevel &dof_bit = (*dit)->getBitRefLevel();
    if((dof_bit&mask) != dof_bit) continue;
    if((dof_bit&bit).any()) {
      EntityHandle ent = (*dit)->getEnt();
      ShortId uid = (*dit)->getNonNonuniqueShortId();
      FieldData val = (*dit)->getFieldData();
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(ent,uid,val); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_begin(const std::string& serie_name) {
  PetscFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->begin(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_end(const std::string& serie_name,double time) {
  PetscFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->end(time); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::initialize_series_recorder(const std::string& serie_name) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->read(moab); CHKERRQ(ierr);
  const_cast<MoFEMSeries*>(&*sit)->record_begin = false;
  const_cast<MoFEMSeries*>(&*sit)->record_end = false;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::finalize_series_recorder(const std::string& serie_name) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==sEries.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = sit->save(moab); CHKERRQ(ierr);
  int nb_steps;
  ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
  int ss = 0;
  for(;ss<nb_steps;ss++) {
    /*std::pair<SeriesStep_multiIndex::iterator,bool> p =*/ seriesSteps.insert(MoFEMSeriesStep(moab,&*sit,ss));
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::print_series_steps() {
  PetscFunctionBegin;
  std::ostringstream ss;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().begin();
  for(;sit!=sEries.get<SeriesName_mi_tag>().end();sit++) {
    ss << "series " << *sit << std::endl;
  }
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit = seriesSteps.get<SeriesName_mi_tag>().begin();
  for(;ssit!=seriesSteps.get<SeriesName_mi_tag>().end();ssit++) {
    ss << "serises steps " << *ssit << std::endl;
  }
  PetscPrintf(comm,ss.str().c_str());
  PetscFunctionReturn(0);
}
bool Core::check_series(const std::string& name) const {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = sEries.get<SeriesName_mi_tag>().find(name);
  if(sit!=sEries.get<SeriesName_mi_tag>().end()) return true;
  return false;
}
PetscErrorCode Core::load_series_data(const std::string& serie_name,const int step_number) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  SeriesStep_multiIndex::index<Composite_SeriesName_And_Step_mi_tag>::type::iterator sit;
  sit = seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().find(boost::make_tuple(serie_name,step_number));
  if(sit == seriesSteps.get<Composite_SeriesName_And_Step_mi_tag>().end()) {
    SETERRQ2(PETSC_COMM_SELF,1,"series <%s> and step %d not found",serie_name.c_str(),step_number);
  }
  ierr = sit->get(moab,dofsField); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator Core::get_series_steps_byName_begin(const std::string& name) {
  return seriesSteps.get<SeriesName_mi_tag>().lower_bound(name);
}
SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator Core::get_series_steps_byName_end(const std::string& name) {
  return seriesSteps.get<SeriesName_mi_tag>().upper_bound(name);
}

}
