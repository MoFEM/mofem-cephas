/** \file SeriesMultiIndices.hpp
 * \brief Myltindex containers, for mofem fields data structures and other low-level functions
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

/**
 * \brief Structure for recording (time) series
 *
 * \ingroup series_multi_indices
 */
struct MoFEMSeries {

  EntityHandle meshset;
  const void* tag_name_data;		///< tag keeps name of the series
  int tag_name_size; 			      ///< number of bits necessary to keep field series

  bool record_begin;
  bool record_end;

  MoFEMSeries(Interface &moab,const EntityHandle _meshset);

  /// get meshset
  inline EntityHandle getMeshset() const { return meshset; }
  inline EntityID get_meshset_id() const { return (EntityID)(meshset&MB_ID_MASK); }
  /// get string_ref of series
  inline boost::string_ref getNameRef() const { return boost::string_ref((char *)tag_name_data,tag_name_size); }
  /// get series name
  inline std::string getName() const { return std::string((char *)tag_name_data,tag_name_size); }

  Tag th_SeriesData;
  Tag th_SeriesDataUIDs;
  Tag th_SeriesDataHandles;
  Tag th_SeriesTime;

  PetscErrorCode get_nb_steps(Interface &moab,int &nb_setps) const;

  std::vector<int> ia;
  std::vector<double> time;
  std::vector<EntityHandle> handles;
  std::vector<ShortId> uids;
  std::vector<FieldData> data;

  PetscErrorCode set_time(double time);
  PetscErrorCode push_dofs(const EntityHandle ent,const ShortId uid,const FieldData val);

  template<typename IT>
  PetscErrorCode push_dofs(IT it,IT hi_it) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for(;it!=hi_it;it++) {
      ierr = push_dofs((*it)->getEnt(),(*it)->getNonNonuniqueShortId(),(*it)->getFieldData()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode begin();
  PetscErrorCode end(double time = 0);
  PetscErrorCode read(Interface &moab);
  PetscErrorCode save(Interface &moab) const;

  inline const MoFEMSeries* get_MoFEMSeries_ptr() const { return const_cast<MoFEMSeries*>(this); };

  friend std::ostream& operator<<(std::ostream& os,const MoFEMSeries& e);


};

template<typename T>
struct interface_MoFEMSeries {
  const T *ptr;
  interface_MoFEMSeries(const T *_ptr): ptr(_ptr) {}

  /// get meshset
  inline EntityHandle getMeshset() const { return ptr->getMeshset(); }
  inline EntityID get_meshset_id() const { return ptr->get_meshset_id(); }
  /// get string_ref of series
  inline boost::string_ref getNameRef() const { return ptr->getNameRef(); }
  /// get series name
  inline std::string getName() const { return ptr->getName(); }

  inline const MoFEMSeries* get_MoFEMSeries_ptr() const { return ptr->get_MoFEMSeries_ptr(); };

};

/**
 * \brief Structure for keeping time and step
 *
 * \ingroup series_multi_indices
 */
struct MoFEMSeriesStep: public interface_MoFEMSeries<MoFEMSeries> {

  typedef interface_MoFEMSeries<MoFEMSeries> interface_type_MoFEMSeries;

  int step_number;
  MoFEMSeriesStep(Interface &moab,const MoFEMSeries *_MoFEMSeries_ptr,const int _step_number);

  inline int get_step_number() const { return step_number; };
  PetscErrorCode get(Interface &moab,DofEntity_multiIndex &dofsField) const;

  double time;
  PetscErrorCode get_time_init(Interface &moab);
  inline double get_time() const { return time; }

  friend std::ostream& operator<<(std::ostream& os,const MoFEMSeriesStep& e);

};

/**
 * \brief Series multi index
 *
 * \ingroup series_multi_indices
 */
typedef multi_index_container<
  MoFEMSeries,
  indexed_by<
  ordered_unique<
    tag<SeriesID_mi_tag>, const_mem_fun<MoFEMSeries,EntityID,&MoFEMSeries::get_meshset_id> >,
  ordered_unique<
    tag<SeriesName_mi_tag>, const_mem_fun<MoFEMSeries,boost::string_ref,&MoFEMSeries::getNameRef> >
  > > Series_multiIndex;

/**
 * \brief Step multi index
 *
 * \ingroup series_multi_indices
 */
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
	     const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::getNameRef>,
	     member<MoFEMSeriesStep,int,&MoFEMSeriesStep::step_number>
      > >,
    ordered_non_unique<
      tag<SeriesName_mi_tag>, const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::getNameRef> >,
        ordered_non_unique<
        tag<Composite_SeriesName_And_Time_mi_tag>,
      composite_key<
	     MoFEMSeriesStep,
	     const_mem_fun<MoFEMSeriesStep::interface_type_MoFEMSeries,boost::string_ref,&MoFEMSeriesStep::getNameRef>,
	     const_mem_fun<MoFEMSeriesStep,double,&MoFEMSeriesStep::get_time>
      > >
  > > SeriesStep_multiIndex;

}

/***************************************************************************//**
 * \defgroup series_multi_indices Series structures and multi-indices
 * \brief Interface used to record fields and play them back (for example response of structure subjected to earthquake)
 *
 * The idea of this interface is taken from Opensees
 * <http://opensees.berkeley.edu>, tailored and generalised for MoFEM code. In
 * principle, one can record tape in for a sequence of points (for example steps
 * in time) and play it back. One can create several tapes and record all fields
 * or some of them.
 *
 * \ingroup mofem
 ******************************************************************************/


#endif // __SERIESMULTIINDICES_HPP__
