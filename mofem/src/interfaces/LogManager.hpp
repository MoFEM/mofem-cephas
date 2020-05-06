/**
 * @file LogManager.hpp
 * @brief Log and register warnings
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

#pragma once

#define __LOGMANAGER_HPP__

namespace MoFEM {

// static const MOFEMuuid IDD_MOFEMLogManager =
//     MOFEMuuid(BitIntefaceId(LogManager_INTERFACE));

/**
 * \brief Problem manager is used to build and partition problems
 * \ingroup mofem_warring_manager
 *
 */
struct LogManager : public UnknownInterface {

  // MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
  //                                UnknownInterface **iface) const;

  LogManager(const MoFEM::Core &core);
  virtual ~LogManager() = default;

  enum WarringType { SELF, WORLD, SYNCHRONIZED };

  MoFEMErrorCode getSubInterfaceOptions();

  /**
   * \brief Get options from command line
   * @return error code
   */
  MoFEMErrorCode getOptions();

  class SelfStreamBuf : public std::stringbuf {
    virtual int sync() {
      if (!this->str().empty()) {
        PetscPrintf(PETSC_COMM_SELF, "%s", this->str().c_str());
        this->str("");
      }
      return 0;
    }
  };

  struct WorldStreamBuf : public std::stringbuf  {
    WorldStreamBuf(MPI_Comm comm) : cOmm(comm) {}
    virtual int sync() {
      if (!this->str().empty()) {
        PetscPrintf(cOmm, "%s", this->str().c_str());
        this->str("");
      }
      return 0;
    }
  private:
    MPI_Comm cOmm;
  };

  struct SynchronizedStreamBuf : public std::stringbuf  {
    SynchronizedStreamBuf(MPI_Comm comm) : cOmm(comm) {}
    virtual int sync() {
      if (!this->str().empty()) {
        PetscSynchronizedPrintf(cOmm, "%s", this->str().c_str());
        PetscSynchronizedFlush(cOmm, PETSC_STDOUT);
        this->str("");
      }
      return 0;
    }

  private:
    MPI_Comm cOmm;
  };

private:
  MoFEM::Core &cOre;
  MPI_Comm cOmm; ///< MoFEM communicator
};
}

/**
 * \defgroup mofem_warring_manager
 * \brief Warning manager
 *
 * \ingroup mofem
 */