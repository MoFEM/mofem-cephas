/** \file MedInterface.hpp
 * \brief Med file interface interface
 *
 * Interface loading mesh and data on mesh directly to mofem & moab
 *
 * \todo Reading higher order entities
 * \todo Reading fields data
 *
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

#ifdef WITH_MED

#ifndef __MED_INTERFACE_HPP__
#define __MED_INTERFACE_HPP__

namespace MoFEM {

/** \brief Interface for load MED files

  * \ingroup mofem_med_files

*/
struct MedInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MedInterface(const MoFEM::Core &core);

  /**
   * \brief Get MED file name from command line
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode getFileNameFromCommandLine(int verb = 1);

  /** \brief Check if file name is given in command line
   */
  inline PetscBool getFlgFile() const { return flgFile; }

  /**
   * \brief Get field names in MED file
   * @param  file file name
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode medGetFieldNames(const string &file, int verb = 1);

  /**
   * \brief get field names in MED file
   *
   * File name is get form command line
   *
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode medGetFieldNames(int verb = 1);

  /**
   * \brief read MED file
   * @param  file filed name
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode readMed(const string &file, int verb = 1);

  /**
   * \brief read MED file
   *
   * File name is form command line
   *
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode readMed(int verb = 1);

  /**
   * \brief write MED file
   * @param  file file name
   * @param  verb verbosity level
   * @return      error code
   */
  MoFEMErrorCode writeMed(const string &file, int verb = 1);

  /**
   * Read fields
   * @param  file_name   file name
   * @param  file_index  mesh index
   * @param  loadSeries  load field into series
   * @param  onlyStep    read only one step
   * @return             error code
   */
  MoFEMErrorCode readFields(const std::string &file_name,
                            const std::string &field_name,
                            const bool load_series = false,
                            const int only_step = -1, int verb = 1);

  struct FieldData {
    std::string fieldName;
    std::string meshName;
    bool localMesh;
    // Med_field_type fieldType;
    std::vector<std::string> componentNames;
    std::vector<std::string> componentUnits;
    std::string dtUnit;
    int ncSteps;
  };

  std::vector<std::string> meshNames;     ///< list of meshes in MED file
  std::vector<EntityHandle> meshMeshsets; ///< meshset for each mesh
  std::map<std::string, FieldData> fieldNames;
  std::string medFileName; ///< MED file name

private:
  MoFEM::Core &cOre; ///< core database
  // std::vector<std::string> fieldNames;  ///< list of fields in MED file
  PetscBool flgFile; ///< true if file name given in command line

  /**
   * \brief read mesh from MED file
   * @param  file            file name
   * @param  index           index of mesh
   * @param  family_elem_map map of families and elements
   * @param  verb            verbosity level
   * @return                 error code
   */
  MoFEMErrorCode readMesh(const string &file, const int index,
                          std::map<int, Range> &family_elem_map, int verb = 1);

  /**
   * \brief read family and groups
   * @param  file               file name
   * @param  index              mesh index
   * @param  family_elem_mapint map of families and elements
   * @param  group_elem_map     map of groups and elements
   * @param  verb               verbosity level
   * @return                    error code
   */
  MoFEMErrorCode readFamily(const string &file, const int index,
                            const std::map<int, Range> &family_elem_mapint,
                            std::map<string, Range> &group_elem_map,
                            int verb = 1);

  /**
   * \brief make from groups meshsets
   * @param  group_elem_map map of groups and elements
   * @param  verb           verbosity level
   * @return                error code
   */
  MoFEMErrorCode makeBlockSets(const std::map<string, Range> &group_elem_map,
                               int verb = 1);
};

std::ostream &operator<<(std::ostream &os,
                         const MedInterface::FieldData &field_data);

} // namespace MoFEM

#endif // __MED_INTERFACE_HPP__
#endif // WITH_MED

/***************************************************************************/ /**
* \defgroup mofem_med_files Reading and writing med files
* \ingroup mofem
*
* \brief Interface for reading and writing med files
******************************************************************************/
