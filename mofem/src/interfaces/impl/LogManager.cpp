/**
 * @file LogManager.cpp
 * @brief Log and register warnings
 * 
 */

#include <MoFEM.hpp>

#undef likely

#define BOOST_LOG_DYN_LINK 

#include <cstddef>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/core/null_deleter.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;


namespace boost {
namespace log {
namespace expressions {
BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity",
                            MoFEM::LogManager::SeverityLevel)
BOOST_LOG_ATTRIBUTE_KEYWORD(tag_attr, "Tag", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(scope, "Scope", attrs::named_scope::value_type)
BOOST_LOG_ATTRIBUTE_KEYWORD(timeline, "Timeline", attrs::timer::value_type)
} // namespace expressions
} // namespace log
} // namespace boost


// The operator puts a human-friendly representation of the severity level to
// the stream
std::ostream &operator<<(std::ostream &strm, MoFEM::LogManager::SeverityLevel level) {
  static const char *strings[] = {"normal", "notification", "warning", "error",
                                  "critical"};

  if (static_cast<std::size_t>(level) < sizeof(strings) / sizeof(*strings))
    strm << strings[level];
  else
    strm << static_cast<int>(level);

  return strm;
}

namespace MoFEM {

class SelfStreamBuf : public std::stringbuf {
  virtual int sync() {
    if (!this->str().empty()) {
      PetscPrintf(PETSC_COMM_SELF, "%s", this->str().c_str());
      this->str("");
    }
    return 0;
  }
};

struct WorldStreamBuf : public std::stringbuf {
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

struct SynchronizedStreamBuf : public std::stringbuf {
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

struct LogManager::InternalData
    : public boost::enable_shared_from_this<LogManager::InternalData> {

  SelfStreamBuf selfBuf;
  WorldStreamBuf worldBuf;
  SynchronizedStreamBuf syncBuf;

  std::ostream strmSelf;
  std::ostream strmWorld;
  std::ostream strmSync;

  InternalData(MPI_Comm comm)
      : worldBuf(comm), syncBuf(comm), strmSelf(&selfBuf), strmWorld(&worldBuf),
        strmSync(&syncBuf) {}
    };

    LogManager::LogManager(const MoFEM::Core &core)
        : cOre(const_cast<MoFEM::Core &>(core)),
          internalDataPtr(new InternalData(
              static_cast<MoFEM::Interface &>(cOre).get_comm())) {}

    MoFEMErrorCode LogManager::query_interface(const MOFEMuuid &uuid,
                                               UnknownInterface **iface) const {
      MoFEMFunctionBeginHot;
      *iface = NULL;
      if (uuid == IDD_MOFEMLogManager) {
        *iface = const_cast<LogManager *>(this);
        MoFEMFunctionReturnHot(0);
      }
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
      MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode LogManager::getSubInterfaceOptions() { return getOptions(); }

MoFEMErrorCode LogManager::getOptions() {
  MoFEMFunctionBegin;
  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "log_",
                           "Warning interface options", "none");
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);

  CHKERR setUpLog();

  MoFEMFunctionReturn(0);
}

boost::shared_ptr<std::ostream> LogManager::getStrmSelf() {
  return boost::shared_ptr<std::ostream>(internalDataPtr,
                                         &internalDataPtr->strmSelf);
}
boost::shared_ptr<std::ostream> LogManager::getStrmWorld() {
  return boost::shared_ptr<std::ostream>(internalDataPtr,
                                         &internalDataPtr->strmWorld);
}
boost::shared_ptr<std::ostream> LogManager::getStrmSync() {
  return boost::shared_ptr<std::ostream>(internalDataPtr,
                                         &internalDataPtr->strmSync);
}

MoFEMErrorCode LogManager::setUpLog() {
  MoFEMFunctionBegin;

  auto stream_ptr = getStrmWorld();

  auto core_log = logging::core::get();
  auto backend = boost::make_shared<sinks::text_ostream_backend>();
  backend->add_stream(stream_ptr);

  typedef sinks::synchronous_sink<sinks::text_ostream_backend> sink_t;
  auto sink = boost::make_shared<sink_t>(backend);

  sink->set_formatter(

      expr::stream
      << std::hex << std::setw(8) << std::setfill('0')
      << boost::log::expressions::line_id << std::dec << std::setfill(' ')
      << ": <" << boost::log::expressions::severity << ">\t"
      << boost::log::expressions::format_named_scope("Scope", keywords::format =
                                                                  "[%f:%l]")
      << "(" << boost::log::expressions::scope << ") "
      << expr::if_(expr::has_attr(boost::log::expressions::tag_attr))
             [expr::stream << "[" << boost::log::expressions::tag_attr << "] "]
      << expr::if_(expr::has_attr(boost::log::expressions::timeline))
             [expr::stream << "[" << boost::log::expressions::timeline << "] "]
      << expr::smessage

  );

  core_log->add_sink(sink);

  logging::add_common_attributes();
  core_log->add_global_attribute("LineID", attrs::counter<unsigned int>(1));
  core_log->add_global_attribute("TimeStamp", attrs::local_clock());
  core_log->add_global_attribute("Scope", attrs::named_scope());

  backend->auto_flush(true);

  MoFEMFunctionReturn(0);
}

} // MOFEM namespace
