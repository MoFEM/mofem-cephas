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

#ifndef __LOGMANAGER_HPP__
#define __LOGMANAGER_HPP__




namespace MoFEM {

static const MOFEMuuid IDD_MOFEMLogManager =
    MOFEMuuid(BitIntefaceId(LOGMANAGER_INTERFACE));

/**
 * \brief Problem manager is used to build and partition problems
 * \ingroup mofem_warring_manager
 *
 */
struct LogManager : public UnknownInterface {

  enum SeverityLevel { normal, notification, warning, error, critical };

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  LogManager(const MoFEM::Core &core);
  virtual ~LogManager() = default;

  enum WarringType { SELF, WORLD, SYNCHRONIZED };

  MoFEMErrorCode getSubInterfaceOptions();

  /**
   * \brief Get options from command line
   * @return error code
   */
  MoFEMErrorCode getOptions();

  boost::shared_ptr<std::ostream> getStrmSelf();
  boost::shared_ptr<std::ostream> getStrmWorld();
  boost::shared_ptr<std::ostream> getStrmSync();

private:
  MoFEM::Core &cOre;

  struct InternalData;
  boost::shared_ptr<InternalData> internalDataPtr;

  MoFEMErrorCode setUpLog();
};
}

#endif //__LOGMANAGER_HPP__

/**
 * \defgroup mofem_warring_manager
 * \brief Warning manager
 *
 * \ingroup mofem
 */