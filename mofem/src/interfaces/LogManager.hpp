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

#define BOOST_LOG_DYN_LINK
#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/expressions.hpp>

namespace MoFEM {}

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMLogManager =
    MOFEMuuid(BitIntefaceId(LOGMANAGER_INTERFACE));

/**
 * \brief Problem manager is used to build and partition problems
 * \ingroup mofem_warring_manager
 *
 */
struct LogManager : public UnknownInterface {

  enum SeverityLevel {
    very_noisy,
    noisy,
    very_verbose,
    verbose,
    warning,
    fault,
    critical
  };

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  LogManager(const MoFEM::Core &core);
  virtual ~LogManager() = default;

  enum WarringType { SELF, WORLD, SYNCHRONIZED };

  MoFEMErrorCode getSubInterfaceOptions();

  typedef boost::log::sources::severity_channel_logger<SeverityLevel,
                                                       std::string>
      LoggerType;

  inline LoggerType getLogSelf();
  inline LoggerType getLogWorld();
  inline LoggerType getLogSync();

  /**
   * \brief Get options from command line
   * @return error code
   */
  MoFEMErrorCode getOptions();

private:
  MoFEM::Core &cOre;

  struct InternalData;
  boost::shared_ptr<InternalData> internalDataPtr;

  MoFEMErrorCode setUpLog();
};

// The operator puts a human-friendly representation of the severity level to
// the stream
std::ostream &operator<<(std::ostream &strm,
                         const LogManager::SeverityLevel &level);

namespace LogKeywords {

BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", LogManager::SeverityLevel)
BOOST_LOG_ATTRIBUTE_KEYWORD(channel, "Channel", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope",
                            boost::log::attributes::named_scope::value_type)
BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline",
                            boost::log::attributes::timer::value_type)

} // namespace LogKeywords

LogManager::LoggerType LogManager::getLogSelf() {
  return LoggerType(boost::log::keywords::channel = "SELF");
}
LogManager::LoggerType LogManager::getLogWorld() {
  return LoggerType(boost::log::keywords::channel = "WORLD");
}
LogManager::LoggerType LogManager::getLogSync() {
  return LoggerType(boost::log::keywords::channel = "WORLD");
}

} // namespace MoFEM

#endif //__LOGMANAGER_HPP__

/**
 * \defgroup mofem_warring_manager
 * \brief Warning manager
 *
 * \ingroup mofem
 */