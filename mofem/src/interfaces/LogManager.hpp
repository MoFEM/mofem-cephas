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
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace attrs = boost::log::attributes;
namespace logging = boost::log;

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
    inform,
    warning,
    error,
    critical
  };

  static constexpr std::array<char *const, critical + 1> severityStrings = {

      (char *)"very_noisy", (char *)"noisy",   (char *)"very_verbose",
      (char *)"verbose",    (char *)"inform",  (char *)"warning",
      (char *)"error",      (char *)"critical"

  };

  enum LogAttributesBits {
    BitLineID = 1 << 0,
    BitScope = 1 << 1,
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

  static void addAttributes(LogManager::LoggerType &lg, const int bit = 0);

  static void addAttributes(const std::string channel, const int bit = 0);

  static LoggerType &setLog(const std::string channel, const int bit = 0);

  static LoggerType &getLog(const std::string channel);

  static void addTag(LogManager::LoggerType &lg, const std::string tag);

  static void addTag(const std::string channel, const std::string tag);

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

BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", LogManager::SeverityLevel)
BOOST_LOG_ATTRIBUTE_KEYWORD(channel, "Channel", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope",
                            boost::log::attributes::named_scope::value_type)
BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline",
                            boost::log::attributes::timer::value_type)

} // namespace LogKeywords

} // namespace MoFEM

#define MOFEM_LOG(channel, severity)                                           \
  BOOST_LOG_SEV(MoFEM::LogManager::getLog(channel), severity)

#endif //__LOGMANAGER_HPP__

/**
 * \defgroup mofem_warring_manager
 * \brief Warning manager
 *
 * \ingroup mofem
 */