/** \file SeriesRecorder.hpp
 * \brief MoFEM interface
 *
 * Interface for recording and saving data series, for example in time or load stepping problems.
 *
 * \ingroup mofem_series
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __SERIESRECORDER_HPP__
#define __SERIESRECORDER_HPP__

#include "UnknownInterface.hpp"

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMSeriesRecorder = MOFEMuuid( BitIntefaceId(SERIES_RECORDER) );

/** Record (time) data series
 * \ingroup mofem_series

  Is abstraction of Core interface.

  \bug fix names of this interface to follow naming convention

 */
struct SeriesRecorder: public UnknownInterface {

  PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

  MoFEM::Core& cOre;
  SeriesRecorder(const MoFEM::Core &core);

  ///destructor
  ~SeriesRecorder() {}

  /**
   * \brief get tags handlers used on meshsets

   */
  PetscErrorCode getTags(int verb = -1);

  inline Tag get_th_SeriesName() { return th_SeriesName; }

  /**
   * \brief clear multi-index container
   * @return error code
   */
  PetscErrorCode clearMap();

  /**
   * \brier initialize container form data on mesh
   * @return error code
   */
  PetscErrorCode initialiseDatabseFromMesh(int verb = 0);

  /**
   * \brief return pointer to meshset manager
   */
  SeriesRecorder* get_series_recorder_ptr() { return this; }

  /**
   * \brief return pointer to meshset manager
   */
  const SeriesRecorder* get_series_recorder_ptr() const { return this; }

 /**
    * \ingroup mofem_series
    * Add series recorder
    *
    * \param name of series
    */
  virtual PetscErrorCode add_series_recorder(const std::string& series_name);

 /**
    * \ingroup mofem_series
    * delete recorded series
    *
    * \param name of series
    */
  virtual PetscErrorCode delete_recorder_series(const std::string& series_name);

  /**
    * \ingroup mofem_series
    * initialize sries for recording
    *
    * \param series name
    */
  virtual PetscErrorCode initialize_series_recorder(const std::string& serie_name);

  /**
    * \ingroup mofem_series
    * finalize series for recording, recorded data are not accessible until finalize
    *
    * \param series name
    */
  virtual PetscErrorCode finalize_series_recorder(const std::string& serie_name);

  /**
    * \ingroup mofem_series
    * begin series recording
    *
    * \param series name
    */
  virtual PetscErrorCode record_begin(const std::string& serie_name);

  /**
    * \ingroup mofem_series
    * record problem
    *
    * \param series name
    * \param problem pointer
    * \param rc could be Row or Col
    */
  virtual PetscErrorCode record_problem(const std::string& serie_name,const Problem *problemPtr,RowColData rc);

  /**
    * \ingroup mofem_series
    * record problem
    *
    * \param series name
    * \param problem name
    * \param rc could be Row or Col
    */
  virtual PetscErrorCode record_problem(const std::string& serie_name,const std::string& problem_name,RowColData rc);

  /**
    * \ingroup mofem_series
    * record field
    *
    * \param field name
    * \param bit ref level
    * \param mask for bit ref level
    */
  virtual PetscErrorCode record_field(const std::string& serie_name,const std::string& field_name,const BitRefLevel &bit,const BitRefLevel &mask);

  /**
    * \ingroup mofem_series
    * end series recording
    *
    * \param series name
    */
  virtual PetscErrorCode record_end(const std::string& serie_name,double time = 0);

  /**
    * \ingroup mofem_series
    * load data from series into dofs database
    *
    * \param series name
    * \param step number
    */
  virtual PetscErrorCode load_series_data(const std::string& serie_name,const int step_number);

  /**
    * \ingroup mofem_series
    * print series
    */
  virtual PetscErrorCode print_series_steps();

  /** \brief check if series is in database
   * \ingroup mofem_series
   *
   * \param name field name
   * \return true if field exist
   *
   */
  virtual bool check_series(const std::string& name) const;

  virtual SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_begin(const std::string& name);
  virtual SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator get_series_steps_byName_end(const std::string& name);

protected:

  Tag th_SeriesName;                    ///< Recorded series name

  Series_multiIndex sEries;							///< recorded series
  SeriesStep_multiIndex seriesSteps;		///< recorded series steps


  /** \brief loop over recorded series step
    * \ingroup mofem_series
    *
    * \param NAME series name
    * \parma IT iterator variable
    *
    * Example: \code
      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"TEST_SERIES1",sit)) {

	ierr = mField.load_series_data("TEST_SERIES1",sit->get_step_number()); CHKERRQ(ierr);
    * } \endcode
    *
    */
  #define _IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(RECORDER,NAME,IT) \
    SeriesStep_multiIndex::index<SeriesName_mi_tag>::type::iterator IT = (RECORDER)->get_series_steps_byName_begin(NAME); \
    IT!=(RECORDER)->get_series_steps_byName_end(NAME); IT++




};

}

#endif // __SERIESRECORDER_HPP__

/***************************************************************************//**
 * \defgroup mofem_series Recording and reading series
 * Recorder for time steps and solution sequences
 *
 * \ingroup mofem
 ******************************************************************************/
