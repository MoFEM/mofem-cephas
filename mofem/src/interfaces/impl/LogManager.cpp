/**
 * @file LogManager.cpp
 * @brief Log and register warnings
 *
 */

#include <MoFEM.hpp>

#undef likely

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
namespace sinks = boost::log::sinks;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;

namespace MoFEM {

using namespace MoFEM::LogKeywords;

std::ostream &operator<<(std::ostream &strm,
                         const LogManager::SeverityLevel &level) {
  static const char *strings[] = {

      "very_noisy", "noisy", "very_verbose", "verbose",
      "warning",    "error", "critical"

  };

  // if (static_cast<std::size_t>(level) < sizeof(strings) / sizeof(*strings))
  strm << strings[level];
  // else
  // strm << static_cast<int>(level);

  return strm;
}

struct LogManager::InternalData
    : public boost::enable_shared_from_this<LogManager::InternalData> {

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

  SelfStreamBuf selfBuf;
  WorldStreamBuf worldBuf;
  SynchronizedStreamBuf syncBuf;

  std::ostream strmSelf;
  std::ostream strmWorld;
  std::ostream strmSync;

  boost::shared_ptr<std::ostream> getStrmSelf() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmSelf);
  }
  boost::shared_ptr<std::ostream> getStrmWorld() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmWorld);
  }
  boost::shared_ptr<std::ostream> getStrmSync() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmSync);
  }

  LoggerType lgSelf;
  LoggerType lgWorld;
  LoggerType lgSync;

  InternalData(MPI_Comm comm)
      : worldBuf(comm), syncBuf(comm), strmSelf(&selfBuf), strmWorld(&worldBuf),
        strmSync(&syncBuf), lgSelf(keywords::channel = "SELF"),
        lgWorld(keywords::channel = "WORLD"),
        lgSync(keywords::channel = "SYNC") {}

  virtual ~InternalData() = default;
};

LogManager::LogManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)),
      internalDataPtr(
          new InternalData(static_cast<MoFEM::Interface &>(cOre).get_comm())) {}

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

MoFEMErrorCode LogManager::setUpLog() {
  MoFEMFunctionBegin;

  auto create_sink = [&](auto stream_ptr, auto comm_filter) {
    auto backend = boost::make_shared<sinks::text_ostream_backend>();
    backend->add_stream(stream_ptr);
    backend->auto_flush(true);

    typedef sinks::synchronous_sink<sinks::text_ostream_backend> sink_t;
    auto sink = boost::make_shared<sink_t>(backend);
    sink->set_filter((expr::has_attr(channel) && channel == comm_filter));

    sink->set_formatter(

        expr::stream

        << std::hex << std::setw(8) << std::setfill('0') << line_id

        << std::dec << std::setfill(' ') << ": <" << severity << ">\t"

        << boost::log::expressions::format_named_scope(
               "Scope", keywords::format = "[%f:%l]")
        << "(" << scope << ") "

        << expr::if_(expr::has_attr(
               tag_attr))[expr::stream << "[" << tag_attr << "] "]

        << expr::if_(expr::has_attr(
               timeline))[expr::stream << "[" << timeline << "] "]
        << expr::smessage

    );

    return sink;
  };

  auto core_log = logging::core::get();
  core_log->add_sink(create_sink(internalDataPtr->getStrmSelf(), "SELF"));
  core_log->add_sink(create_sink(internalDataPtr->getStrmWorld(), "WORLD"));
  core_log->add_sink(create_sink(internalDataPtr->getStrmSync(), "SYNC"));

  logging::add_common_attributes();
  core_log->add_global_attribute("LineID", attrs::counter<unsigned int>(1));
  core_log->add_global_attribute("TimeStamp", attrs::local_clock());
  core_log->add_global_attribute("Scope", attrs::named_scope());

  auto &strm = *(internalDataPtr->getStrmSelf());
  strm << SeverityLevel::very_noisy << endl; 
  strm.flush();
  // cout << strm;

  MoFEMFunctionReturn(0);
}

LogManager::LoggerType &LogManager::getLogSelf() {
  return internalDataPtr->lgSelf;
}
LogManager::LoggerType &LogManager::getLogWorld() {
  return internalDataPtr->lgWorld;
}
LogManager::LoggerType &LogManager::getLogSync() {
  return internalDataPtr->lgSync;
}

} // namespace MoFEM
