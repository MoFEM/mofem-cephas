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

namespace attrs = boost::log::attributes;
namespace logging = boost::log;
namespace keywords = boost::log::keywords;
namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace attrs = boost::log::attributes;
namespace expr = boost::log::expressions;

namespace MoFEM {

static const MOFEMuuid IDD_MOFEMLogManager =
    MOFEMuuid(BitIntefaceId(LOGMANAGER_INTERFACE));

/**
 * \brief Log manager is used to build and partition problems
 * \ingroup mofem_log_manager
 *
 */
struct LogManager : public UnknownInterface {

  /**
   * @brief Severity levels
   * \ingroup mofem_log_manager
   * 
   */
  enum SeverityLevel { noisy, verbose, inform, warning, error };

  static constexpr std::array<char *const, error + 1> severityStrings = {

      (char *)"noisy", (char *)"verbose", (char *)"inform", (char *)"warning",
      (char *)"error"

  };

  /**
   * @brief Tag attributes switches
   * \ingroup mofem_log_manager
   * 
   */
  enum LogAttributesBits {
    BitLineID = 1 << 0,
    BitScope = 1 << 1,
  };

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const;

  LogManager(const MoFEM::Core &core);
  virtual ~LogManager() = default;

  /**
   * @brief Definition of the channel logger
   * 
   */
  typedef boost::log::sources::severity_channel_logger<SeverityLevel,
                                                       std::string>
      LoggerType;

  typedef sinks::synchronous_sink<sinks::text_ostream_backend> SinkType;

  /**
   * @brief Add attributes to logger
   * \ingroup mofem_log_manager
   * 
   * @param lg 
   * @param bit 
   */
  static void addAttributes(LogManager::LoggerType &lg, const int bit = 0);

  /**
   * @brief Add attributes to channel
   * \ingroup mofem_log_manager
   * 
   * @param channel 
   * @param bit 
   */
  static void addAttributes(const std::string channel, const int bit = 0);

  /**
   * @brief Set ans resset chanel logger
   * \ingroup mofem_log_manager
   * 
   * @param channel 
   * @return LoggerType& 
   */
  static LoggerType &setLog(const std::string channel);

  /**
   * @brief Get logger by channel
   * \ingroup mofem_log_manager
   * 
   * @param channel 
   * @return LoggerType& 
   */
  static LoggerType &getLog(const std::string channel);

  /**
   * @brief Add tag to logger
   * \ingroup mofem_log_manager
   * 
   * @param lg 
   * @param tag 
   */
  static void addTag(LogManager::LoggerType &lg, const std::string tag);

  /**
   * @brief Add tag to channel
   * \ingroup mofem_log_manager
   * 
   * @param channel 
   * @param tag 
   */
  static void addTag(const std::string channel, const std::string tag);

  /**
   * @brief Get the strm self object
   * 
   * @return boost::shared_ptr<std::ostream> 
   */
  static boost::shared_ptr<std::ostream> getStrmSelf(); 

  /**
   * @brief Get the strm world object
   * 
   * @return boost::shared_ptr<std::ostream> 
   */
  static boost::shared_ptr<std::ostream> getStrmWorld();

  /**
   * @brief Get the strm sync object
   * 
   * @return boost::shared_ptr<std::ostream> 
   */
  static boost::shared_ptr<std::ostream> getStrmSync();

  /**
   * @brief Create a sink object
   *
   * @param stream_ptr
   * @param comm_filter
   * @return boost::shared_ptr<SinkType>
   */
  static boost::shared_ptr<SinkType>
  createSink(boost::shared_ptr<std::ostream> stream_ptr,
             std::string comm_filter);

  /**
   * @brief Create default sinks
   *
   */
  static void createDefaultSinks(MPI_Comm comm);

  /** @brief Get logger option 
   *
   * This function is called by MoFEM core when this interface is registred
   * into database.
   *
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode getOptions();

  /**
   * @brief Dummy file pointer (DO NOT USE)
   * 
   * \note This is for internal use only/
   * 
   */
  static FILE *dummy_mofem_fd;

  /**
   * @brief Use to handle PETSc output
   * 
   * \note This is for internal use only/
   * 
   * @param fd 
   * @param format 
   * @param Argp 
   * @return PetscErrorCode 
   */
  static PetscErrorCode logPetscFPrintf(FILE *fd, const char format[],
                                        va_list Argp);

  /**
   * @brief Converts formatted output to string
   * 
   * @param fmt 
   * @param args 
   * @return std::string 
   */
  static std::string getVLikeFormatedString(const char *fmt, va_list args);

  /**
   * @brief Converts formatted output to string
   * 
   * @param fmt 
   * @param args 
   * @return std::string 
   */
  static std::string getCLikeFormatedString(const char *fmt, ...);

  /**
   * @brief Default record formatter
   * 
   * @param rec 
   * @param strm 
   */
  static void recordFormatterDefault(logging::record_view const &rec,
                                     logging::formatting_ostream &strm);

private:
  MoFEM::Core &cOre;

  struct InternalData;
  static boost::shared_ptr<InternalData> internalDataPtr;

  static std::string petscStringCache;

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
BOOST_LOG_ATTRIBUTE_KEYWORD(proc_attr, "Proc", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope",
                            boost::log::attributes::named_scope::value_type)
BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline",
                            boost::log::attributes::timer::value_type)

} // namespace LogKeywords

} // namespace MoFEM

extern "C" {
PetscErrorCode PetscVFPrintfDefault(FILE *fd, const char *format, va_list Argp);
}

/**
 * @brief Set and reset channel
 * \ingroup mofem_log_manager
 *
 * \code
 * MOFEM_LOG_CHANNEL("WORLD");
 * \endcode
 *
 * Are three default type of channels, SELF, ech processor prints to the
 * standard output, WORLD, only processor one prints, and SYNC all processors
 * prints synchronously.
 *
 *
 */
#define MOFEM_LOG_CHANNEL(channel)                                             \
  { MoFEM::LogManager::setLog(channel); }

/**
 * @brief Add attributes to channel
 * \ingroup mofem_log_manager
 * 
 * \code
 * MOFEM_LOG_ATTRIBUTES("SYNC", LogManager::BitLineID | LogManager::BitScope);
 * \endcode
 * 
 */
#define MOFEM_LOG_ATTRIBUTES(channel, bit)                                     \
  { MoFEM::LogManager::addAttributes(channel, bit); }

/**
 * @brief Log
 * \ingroup mofem_log_manager
 *
 * \code
 * MOFEM_LOG("WORLD", LogManager::SeverityLevel::inform) << "Hello world";
 * \endcode
 *
 */
#define MOFEM_LOG(channel, severity)                                           \
  BOOST_LOG_SEV(MoFEM::LogManager::getLog(channel), severity)

#define MOFEM_LOG_C(channel, severity, format, ...)                            \
  MOFEM_LOG(channel, severity)                                                 \
      << MoFEM::LogManager::getCLikeFormatedString(format, __VA_ARGS__)

/** \brief Set scope
 * \ingroup mofem_log_manager
 * 
 * Macro for function scope markup. The scope name is constructed with help of
 * compiler and contains the current function signature. The scope name is
 * pushed to the end of the current thread scope list.
 *
 * Not all compilers have support for this macro. The exact form of the scope
 * name may vary from one compiler to another.
 */
#define MOFEM_LOG_FUNCTION()                                                   \
  BOOST_LOG_NAMED_SCOPE_INTERNAL(                                              \
      BOOST_LOG_UNIQUE_IDENTIFIER_NAME(_boost_log_named_scope_sentry_),        \
      PETSC_FUNCTION_NAME, __FILE__, __LINE__,                                 \
      ::boost::log::attributes::named_scope_entry::function)

/**
 * @brief Tag channel
 * \ingroup mofem_log_manager
 *
 * Tag channel tag is set until MOFEM_LOG_CHANNEL is called, then new tag can be
 * set.
 *
 */
#define MOFEM_LOG_TAG(channel, tag) MoFEM::LogManager::addTag(channel, tag);

/**
 * @brief Synchronise "SYNC" channel
 * 
 */
#define MOFEM_LOG_SYNCHORMISE(comm)                                            \
  PetscSynchronizedFlush(comm, MoFEM::LogManager::dummy_mofem_fd);

#endif //__LOGMANAGER_HPP__

/**
 * \defgroup mofem_log_manager
 * \brief Log manager
 *
 * Logging manager based on Boost.Log
 * (<a href="https://www.boost.org/doc/libs/1_63_0/libs/log/doc/html/index.html">Boost.Log v2</a>)
 *
 * \ingroup mofem
 */