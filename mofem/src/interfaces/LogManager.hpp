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

  struct SelfStream {
    SelfStream() = default;

    template<typename T>
    SelfStream &operator<<(const std::ostringstream &os) {
      PetscPrintf(PETSC_COMM_SELF, "%s", os.str().c_str());
      return *this;
    };
  };

  struct WorldStream {
    WorldStream(MPI_Comm comm) : cOmm(comm) {}
    WorldStream &operator<<(const std::ostringstream &os) {
      PetscPrintf(cOmm, "%s", os.str().c_str());
      return *this;
    }

  private:
    MPI_Comm cOmm;
  };

  struct SynchronizedStream {
    SynchronizedStream(MPI_Comm comm) : cOmm(comm) {}
    SynchronizedStream &operator<<(const std::ostringstream &os) {
      PetscSynchronizedPrintf(cOmm, "%s", os.str().c_str());
      return *this;
    }
    void flush() { PetscSynchronizedFlush(cOmm, PETSC_STDOUT); }

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