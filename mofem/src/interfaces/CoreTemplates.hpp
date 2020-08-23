/** \file CoreTemplates.hpp
 * \brief Core interface class for user interface
 *
 * Low level data structures not used directly by user
 *
 * FIXME It is a mess with names of core cpp files need better organization
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

#include <CoordSysMultiIndices.hpp>
#include <CoordSystemsManager.hpp>
#include <LogManager.hpp>
#include <BitRefManager.hpp>
#include <MeshsetsManager.hpp>
#include <SeriesRecorder.hpp>

namespace MoFEM {

using Sev = MoFEM::LogManager::SeverityLevel;

template <int N> MoFEMErrorCode CoreTmp<N>::rebuild_database(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = this->verbose;
  CHKERR this->clearMap();
  CHKERR this->getTags(verb);
  CHKERR this->initialiseDatabaseFromMesh(verb);
  MoFEMFunctionReturn(0);
}

template <int N>
MoFEMErrorCode CoreTmp<N>::set_moab_interface(moab::Interface &new_moab,
                                              int verb) {
  return this->setMoabInterface(new_moab, verb);
};

template <int N>
boost::shared_ptr<RefEntityTmp<0>>
CoreTmp<N>::make_shared_ref_entity(const EntityHandle ent) {
  return Core::makeSharedRefEntity(*this, ent);
}

template <int V>
CoreTmp<0>::CoreTmp(moab::Interface &moab, MPI_Comm comm, const int verbose,
                    CoreValue<V>)
    : moab(moab), cOmm(0), verbose(verbose),
      initaliseAndBuildField(PETSC_FALSE),
      initaliseAndBuildFiniteElements(PETSC_FALSE) {

  ierr = coreGenericConstructor(moab, comm, verbose);
  CHKERRABORT(comm, ierr);

  if (verbose > QUIET) {
    MOFEM_LOG_CHANNEL("WORLD");
    MOFEM_LOG_C("WORLD", Sev::verbose, "Core number < %d >", V);
  }
}

template <int N>
CoreTmp<N>::CoreTmp(moab::Interface &moab, ///< MoAB interface
                    MPI_Comm comm,         ///< MPI communicator
                    const int verbose      ///< Verbosity level
                    )
    : CoreTmp<N - 1>(moab, comm, verbose, CoreValue<N>()) {

  // Register sub-interfaces
  ierr = this->registerSubInterfaces();
  CHKERRABORT(comm, ierr);
  ierr = this->clearMap();
  CHKERRABORT(comm, ierr);
  ierr = this->getTags();
  CHKERRABORT(comm, ierr);
  ierr = this->getOptions(verbose);
  CHKERRABORT(comm, ierr);

  this->basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  this->setRefEntBasicDataPtr(*this, this->basicEntityDataPtr);

  ierr = this->initialiseDatabaseFromMesh(verbose);
  CHKERRABORT(comm, ierr);
}

} // namespace MoFEM
