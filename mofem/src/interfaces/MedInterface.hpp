/** \file MedInterface.hpp
 * \brief Med file interface interface
 *
 * Interface loading mesh and data on mesh directly to mofem & moab
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

  static const MOFEMuuid IDD_MOFEMMedInterface = MOFEMuuid( BitIntefaceId(MED_INTERFACE) );

  /** \brief Interface for load MED files
  * \ingroup mofem
  */
  struct MedInterface: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface);

    MedInterface(const MoFEM::Core& core);

    PetscErrorCode getFileNameFromCommandLine(int verb = 1);

    /** \brief Check if file name is given in command line
    */
    inline PetscBool getFlgFile() const { return flgFile; }

    PetscErrorCode medGetFieldNames(const string &file,int verb = 1);

    PetscErrorCode medGetFieldNames(int verb = 1);

    PetscErrorCode readMed(const string &file,int verb = 1);

    PetscErrorCode readMed(int verb = 1);

    PetscErrorCode writeMed(const string &file,int verb = 1);

  private:

    MoFEM::Core& cOre;
    std::vector<std::string> fieldNames;
    std::vector<std::string> meshNames;
    std::string medFileName;
    PetscBool flgFile;

    PetscErrorCode readMesh(
      const string &file,
      const int index,
      std::map<int,Range> &family_elem_map,
      int verb = 1
    );

    PetscErrorCode readFamily(
      const string &file,
      const int index,
      const std::map<int,Range> &family_elem_mapint,
      std::map<string,Range> &group_elem_map,
      int verb = 1
    );

    PetscErrorCode makeBlockSets(
      const std::map<string,Range> &group_elem_map,
      int verb = 1
    );

  };

}

#endif // __MED_INTERFACE_HPP__
#endif // WITH_MED
