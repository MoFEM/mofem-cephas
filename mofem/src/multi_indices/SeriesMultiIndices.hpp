/** \file SeriesMultiIndices.hpp
 * \brief Multindex containers, for mofem fields data structures and other
 low-level functions
 *
 */

#ifndef __SERIESMULTIINDICES_HPP__
#define __SERIESMULTIINDICES_HPP__

namespace MoFEM {

struct FieldSeriesStep;

/**
 * \brief Structure for recording (time) series
 *
 * \ingroup series_multi_indices
 *
 * \bug Fix member functions names. Not follow naming convention.
 * \bug Member functions not documented
 *
 */
struct FieldSeries {

  EntityHandle meshset;
  const void *tagName; ///< tag keeps name of the series
  int tagNameSize;     ///< number of bits necessary to keep field series

  bool record_begin;
  bool record_end;

  FieldSeries(Interface &moab, const EntityHandle _meshset);

  /// get meshset
  inline EntityHandle getMeshset() const { return meshset; }
  inline EntityID get_meshset_id() const {
    return (EntityID)(meshset & MB_ID_MASK);
  }
  /// get string_ref of series
  inline boost::string_ref getNameRef() const {
    return boost::string_ref((char *)tagName, tagNameSize);
  }
  /// get series name
  inline std::string getName() const {
    return std::string((char *)tagName, tagNameSize);
  }

  Tag th_SeriesData;
  Tag th_SeriesDataUIDs;
  Tag th_SeriesDataHandles;
  Tag th_SeriesTime;

  MoFEMErrorCode get_nb_steps(Interface &moab, int &nb_setps) const;

  std::vector<int> ia;
  std::vector<double> time;
  std::vector<EntityHandle> handles;
  std::vector<ShortId> uids;
  std::vector<FieldData> data;

  MoFEMErrorCode set_time(double time);
  MoFEMErrorCode push_dofs(const EntityHandle ent, const ShortId uid,
                           const FieldData val);

  template <typename IT> MoFEMErrorCode push_dofs(IT it, IT hi_it) {
    MoFEMFunctionBeginHot;

    for (; it != hi_it; it++) {
      ierr = push_dofs((*it)->getEnt(), (*it)->getNonNonuniqueShortId(),
                       (*it)->getFieldData());
      CHKERRG(ierr);
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode begin();
  MoFEMErrorCode end(double time = 0);
  MoFEMErrorCode read(Interface &moab);
  MoFEMErrorCode save(Interface &moab) const;

  inline const FieldSeries *get_FieldSeries_ptr() const {
    return const_cast<FieldSeries *>(this);
  };

  friend std::ostream &operator<<(std::ostream &os, const FieldSeries &e);
};

template <typename T> struct interface_FieldSeries {
  const T *ptr;
  interface_FieldSeries(const T *_ptr) : ptr(_ptr) {}

  /// get meshset
  inline EntityHandle getMeshset() const { return ptr->getMeshset(); }
  inline EntityID get_meshset_id() const { return ptr->get_meshset_id(); }
  /// get string_ref of series
  inline boost::string_ref getNameRef() const { return ptr->getNameRef(); }
  /// get series name
  inline std::string getName() const { return ptr->getName(); }

  inline const FieldSeries *get_FieldSeries_ptr() const {
    return ptr->get_FieldSeries_ptr();
  };
};

/**
 * \brief Structure for keeping time and step
 *
 * \ingroup series_multi_indices
 *
 * \bug Fix member functions names. Not follow namining convention.
 */
struct FieldSeriesStep : public interface_FieldSeries<FieldSeries> {

  typedef interface_FieldSeries<FieldSeries> interface_type_FieldSeries;

  int step_number;
  FieldSeriesStep(Interface &moab, const FieldSeries *_FieldSeries_ptr,
                  const int _step_number);

  inline int get_step_number() const { return step_number; };
  MoFEMErrorCode get(Interface &moab, DofEntity_multiIndex &dofsField) const;

  double time;
  MoFEMErrorCode get_time_init(Interface &moab);
  inline double get_time() const { return time; }

  friend std::ostream &operator<<(std::ostream &os, const FieldSeriesStep &e);
};

/**
 * \brief Series multi index
 *
 * \ingroup series_multi_indices
 */
typedef multi_index_container<
    FieldSeries,
    indexed_by<ordered_unique<tag<SeriesID_mi_tag>,
                              const_mem_fun<FieldSeries, EntityID,
                                            &FieldSeries::get_meshset_id>>,
               ordered_unique<tag<SeriesName_mi_tag>,
                              const_mem_fun<FieldSeries, boost::string_ref,
                                            &FieldSeries::getNameRef>>>>
    Series_multiIndex;

/**
 * \brief Step multi index
 *
 * \ingroup series_multi_indices
 */
typedef multi_index_container<
    FieldSeriesStep,
    indexed_by<
        ordered_unique<
            tag<Composite_SeriesID_And_Step_mi_tag>,
            composite_key<
                FieldSeriesStep,
                const_mem_fun<FieldSeriesStep::interface_type_FieldSeries,
                              EntityID, &FieldSeriesStep::get_meshset_id>,
                member<FieldSeriesStep, int, &FieldSeriesStep::step_number>>>,
        ordered_unique<
            tag<Composite_SeriesName_And_Step_mi_tag>,
            composite_key<
                FieldSeriesStep,
                const_mem_fun<FieldSeriesStep::interface_type_FieldSeries,
                              boost::string_ref, &FieldSeriesStep::getNameRef>,
                member<FieldSeriesStep, int, &FieldSeriesStep::step_number>>>,
        ordered_non_unique<
            tag<SeriesName_mi_tag>,
            const_mem_fun<FieldSeriesStep::interface_type_FieldSeries,
                          boost::string_ref, &FieldSeriesStep::getNameRef>>,
        ordered_non_unique<
            tag<Composite_SeriesName_And_Time_mi_tag>,
            composite_key<
                FieldSeriesStep,
                const_mem_fun<FieldSeriesStep::interface_type_FieldSeries,
                              boost::string_ref, &FieldSeriesStep::getNameRef>,
                const_mem_fun<FieldSeriesStep, double,
                              &FieldSeriesStep::get_time>>>>>
    SeriesStep_multiIndex;

} // namespace MoFEM

/***************************************************************************/ /**
* \defgroup series_multi_indices Series structures and multi-indices
* \brief Interface used to record fields and play them back (for example response of structure subjected to earthquake)
*
* The idea of this interface is taken from Opensees
* <http://opensees.berkeley.edu>, tailored and generalised for MoFEM code. In
* principle, one can record tape in for a sequence of points (for example steps
* in time) and play it back. One can create several tapes and record all fields
* or some of them.
*
* User do not use multi-indices data structures directly, recording and playing
* data is by MoFEM::SeriesRecorder interface, more in module \ref mofem_series
*
* \ingroup mofem
******************************************************************************/

#endif // __SERIESMULTIINDICES_HPP__
