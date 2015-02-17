/** \file SeriesMultiIndices.hpp 
 * \brief Myltindex containes, for mofem fields data structures and other low-level functions 
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

#ifndef __SERIESMULTIINDICES_HPP__
#define __SERIESMULTIINDICES_HPP__

namespace MoFEM {

struct MoFEMSeriesStep;

struct MoFEMSeries {

  EntityHandle meshset;
  const void* tag_name_data;		///< tag keeps name of the series
  int tag_name_size; 			///< number of bits necessery to keep field series

  bool record_begin;
  bool record_end;

  MoFEMSeries(Interface &moab,const EntityHandle _meshset);

  /// get meshset
  inline EntityHandle get_meshset() const { return meshset; }
  inline EntityID get_meshset_id() const { return (EntityID)(meshset&MB_ID_MASK); }
  /// get string_ref of series
  inline boost::string_ref get_name_ref() const { return boost::string_ref((char *)tag_name_data,tag_name_size); }
  /// get series name
  inline string get_name() const { return string((char *)tag_name_data,tag_name_size); }
  
  Tag th_SeriesData;
  Tag th_SeriesDataUIDs;
  Tag th_SeriesDataHandles;
  Tag th_SeriesTime;

  PetscErrorCode get_nb_steps(Interface &moab,int &nb_setps) const;

  vector<int> ia;
  vector<double> time;
  vector<EntityHandle> handles;
  vector<ShortId> uids;
  vector<FieldData> data;

  PetscErrorCode set_time(double time);
  PetscErrorCode push_dofs(const EntityHandle ent,const ShortId uid,const FieldData val);
 
  template<typename IT>
  PetscErrorCode push_dofs(IT it,IT hi_it) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for(;it!=hi_it;it++) {
      ierr = push_dofs(it->get_ent(),it->get_non_nonunique_short_id(),it->get_FieldData()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode begin();
  PetscErrorCode end(double time = 0);
  PetscErrorCode read(Interface &moab);
  PetscErrorCode save(Interface &moab) const;

  inline const MoFEMSeries* get_MoFEMSeries_ptr() const { return const_cast<MoFEMSeries*>(this); };

  friend ostream& operator<<(ostream& os,const MoFEMSeries& e);


};

template<typename T>
struct interface_MoFEMSeries {
  const T *ptr;
  interface_MoFEMSeries(const T *_ptr): ptr(_ptr) {}

  /// get meshset
  inline EntityHandle get_meshset() const { return ptr->get_meshset(); }
  inline EntityID get_meshset_id() const { return ptr->get_meshset_id(); }
  /// get string_ref of series
  inline boost::string_ref get_name_ref() const { return ptr->get_name_ref(); }
  /// get series name
  inline string get_name() const { return ptr->get_name(); }

  inline const MoFEMSeries* get_MoFEMSeries_ptr() const { return ptr->get_MoFEMSeries_ptr(); };

};

struct MoFEMSeriesStep: public interface_MoFEMSeries<MoFEMSeries> {

  typedef interface_MoFEMSeries<MoFEMSeries> interface_type_MoFEMSeries;

  int step_number;
  MoFEMSeriesStep(Interface &moab,const MoFEMSeries *_MoFEMSeries_ptr,const int _step_number);

  inline int get_step_number() const { return step_number; };
  PetscErrorCode get(Interface &moab,DofMoFEMEntity_multiIndex &dofsMoabField) const;
 
  double time;
  PetscErrorCode get_time_init(Interface &moab);
  inline double get_time() const { return time; }

  friend ostream& operator<<(ostream& os,const MoFEMSeriesStep& e);

};

typedef multi_index_container<
  MoFEMSeries,
  indexed_by<
  ordered_unique<
    tag<SeriesID_mi_tag>, const_mem_fun<MoFEMSeries,EntityID,&MoFEMSeries::get_meshset_id> >,
  ordered_unique<
    tag<SeriesName_mi_tag>, const_mem_fun<MoFEMSeries,boost::string_ref,&MoFEMSeries::get_name_ref> >
  > > Series_multiIndex;

typedef multi_index_container<
  MoFEMSeriesStep,
  indexed_by<
   ordered_unique<
      tag<Composite_SeriesID_And_Step_mi_tag>, 
      composite_key<
	MoFEMSeriesStep,
	const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,EntityID,&MoFEMSeriesStep::get_meshset_id>,
	member<MoFEMSeriesStep,int,&MoFEMSeriesStep::step_number>
      > >,
   ordered_unique<
      tag<Composite_SeriesName_And_Step_mi_tag>, 
      composite_key<
	MoFEMSeriesStep,
	const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::get_name_ref>,
	member<MoFEMSeriesStep,int,&MoFEMSeriesStep::step_number>
      > >,
    ordered_non_unique<
      tag<SeriesName_mi_tag>, const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::get_name_ref> >,
   ordered_non_unique<
      tag<Composite_SeriesName_And_Time_mi_tag>, 
      composite_key<
	MoFEMSeriesStep,
	const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::get_name_ref>,
	const_mem_fun<MoFEMSeriesStep,double,&MoFEMSeriesStep::get_time>
      > >
  > > SeriesStep_multiIndex;

}

/***************************************************************************//**
 * \defgroup series_multi_indices Series structures and multi-indices
 * \ingroup mofem
 ******************************************************************************/


#endif // __SERIESMULTIINDICES_HPP__
