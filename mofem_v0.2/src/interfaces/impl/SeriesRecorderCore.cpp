/** \file SeriesRecorderCore.cpp
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

#include <moab/ParallelComm.hpp>

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h> 

#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>

#include <LoopMethods.hpp>

#include <boost/ptr_container/ptr_map.hpp>
#include <Core.hpp>

#include <CoreDataStructures.hpp>

namespace MoFEM {

const static int debug = 1;

PetscErrorCode Core::add_series_recorder(const string &series_name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
  void const* tag_data[] = { series_name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = series_name.size();
  rval = moab.tag_set_by_ptr(th_SeriesName,&meshset,1,tag_data,tag_sizes); CHKERR_PETSC(rval);
  pair<Series_multiIndex::iterator,bool> p = series.insert(MoFEMSeries(moab,meshset));
  if(!p.second) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> is already there",series_name.c_str());
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::delete_recorder_series(const string& series_name) {
  PetscFunctionBegin;
  //PetscErrorCode ierr;
  ErrorCode rval;
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit,hi_ssit;
  ssit = series_steps.get<SeriesName_mi_tag>().lower_bound(series_name);
  hi_ssit = series_steps.get<SeriesName_mi_tag>().upper_bound(series_name);
  series_steps.get<SeriesName_mi_tag>().erase(ssit,hi_ssit);
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit;
  sit = series.get<SeriesName_mi_tag>().find(series_name);
  if(sit == series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist and can be deleted",series_name.c_str());
  }
  EntityHandle series_meshset = sit->get_meshset();
  rval = moab.tag_delete(sit->th_SeriesTime); CHKERR_PETSC(rval);
  rval = moab.tag_delete(sit->th_SeriesDataHandles); CHKERR_PETSC(rval);
  rval = moab.tag_delete(sit->th_SeriesDataUIDs); CHKERR_PETSC(rval);
  rval = moab.tag_delete(sit->th_SeriesData); CHKERR_PETSC(rval);
  series.get<SeriesName_mi_tag>().erase(sit);
  vector<EntityHandle> contained;
  rval = moab.get_contained_meshsets(series_meshset,contained); CHKERR_PETSC(rval);
  rval = moab.remove_entities(series_meshset,&contained[0],contained.size()); CHKERR_PETSC(rval);
  rval = moab.delete_entities(&contained[0],contained.size()); CHKERR_PETSC(rval);
  rval = moab.delete_entities(&series_meshset,1); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_problem(const string& serie_name,const MoFEMProblem *problemPtr,RowColData rc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  switch (rc) {
    case ROW:
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(problemPtr->numered_dofs_rows.begin(),problemPtr->numered_dofs_rows.end()); CHKERRQ(ierr);
      break;
    case COL:
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(problemPtr->numered_dofs_cols.begin(),problemPtr->numered_dofs_cols.end()); CHKERRQ(ierr);
      break;
    default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INSONSISTENCY,"data inconsistency");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_problem(const string& serie_name,const string& problem_name,RowColData rc) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const MoFEMProblem *problemPtr;
  ierr = get_problem(problem_name,&problemPtr); CHKERRQ(ierr);
  ierr = record_problem(serie_name,problemPtr,rc); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_field(const string& serie_name,const string& field_name,
  const BitRefLevel &bit,const BitRefLevel &mask) {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  PetscFunctionBegin;
  PetscErrorCode ierr;
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator dit = dofsMoabField.get<FieldName_mi_tag>().lower_bound(field_name);
  if(dit == dofsMoabField.get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"field <%s> not exist",field_name.c_str());
  }
  DofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator hi_dit = dofsMoabField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    const BitRefLevel &dof_bit = dit->get_BitRefLevel();
    if((dof_bit&mask) != dof_bit) continue;
    if((dof_bit&bit).any()) {
      EntityHandle ent = dit->get_ent();
      ShortId uid = dit->get_non_nonunique_short_id();
      FieldData val = dit->get_FieldData();
      ierr = const_cast<MoFEMSeries*>(&*sit)->push_dofs(ent,uid,val); CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_begin(const string& serie_name) {
  PetscFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->begin(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::record_end(const string& serie_name,double time) {
  PetscFunctionBegin;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->end(time); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::initialize_series_recorder(const string& serie_name) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = const_cast<MoFEMSeries*>(&*sit)->read(moab); CHKERRQ(ierr);
  const_cast<MoFEMSeries*>(&*sit)->record_begin = false;
  const_cast<MoFEMSeries*>(&*sit)->record_end = false;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::finalize_series_recorder(const string& serie_name) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(serie_name);
  if(sit==series.get<SeriesName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"serie recorder <%s> not exist",serie_name.c_str());
  }
  ierr = sit->save(moab); CHKERRQ(ierr);
  int nb_steps;
  ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
  int ss = 0;
  for(;ss<nb_steps;ss++) {
    /*pair<SeriesStep_multiIndex::iterator,bool> p =*/ series_steps.insert(MoFEMSeriesStep(moab,&*sit,ss));
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::print_series_steps() {
  PetscFunctionBegin;
  ostringstream ss;
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().begin();
  for(;sit!=series.get<SeriesName_mi_tag>().end();sit++) {
    ss << "series " << *sit << endl;
  }
  SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator ssit = series_steps.get<SeriesName_mi_tag>().begin();
  for(;ssit!=series_steps.get<SeriesName_mi_tag>().end();ssit++) {
    ss << "serises steps " << *ssit << endl;
  }
  PetscPrintf(comm,ss.str().c_str());
  PetscFunctionReturn(0);
}
bool Core::check_series(const string& name) const {
  Series_multiIndex::index<SeriesName_mi_tag>::type::iterator sit = series.get<SeriesName_mi_tag>().find(name);
  if(sit!=series.get<SeriesName_mi_tag>().end()) return true;
  return false;
}
PetscErrorCode Core::load_series_data(const string& serie_name,const int step_number) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  SeriesStep_multiIndex::index<Composite_SeriesName_And_Step_mi_tag>::type::iterator sit;
  sit = series_steps.get<Composite_SeriesName_And_Step_mi_tag>().find(boost::make_tuple(serie_name,step_number));
  if(sit == series_steps.get<Composite_SeriesName_And_Step_mi_tag>().end()) {
    SETERRQ2(PETSC_COMM_SELF,1,"series <%s> and step %d not found",serie_name.c_str(),step_number);
  }
  ierr = sit->get(moab,dofsMoabField); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator Core::get_series_steps_byName_begin(const string& name) {
  return series_steps.get<SeriesName_mi_tag>().lower_bound(name);
}
SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator Core::get_series_steps_byName_end(const string& name) {
  return series_steps.get<SeriesName_mi_tag>().upper_bound(name);
}

}
