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

static int dummy_file_ptr = 1;

namespace MoFEM {

using namespace MoFEM::LogKeywords;

constexpr std::array<char *const, LogManager::SeverityLevel::error + 1>
    LogManager::severityStrings;

std::ostream &operator<<(std::ostream &strm,
                         const LogManager::SeverityLevel &level) {

  strm << "<" << LogManager::severityStrings[level] << "> ";

  return strm;
}

struct LogManager::InternalData
    : public boost::enable_shared_from_this<LogManager::InternalData> {

  class SelfStreamBuf : public std::stringbuf {
    virtual int sync() {
      if (!this->str().empty()) {
        PetscFPrintf(PETSC_COMM_SELF, LogManager::dummy_mofem_fd, "%s",
                     this->str().c_str());
        this->str("");
      }
      return 0;
    }
  };

  struct WorldStreamBuf : public std::stringbuf {
    WorldStreamBuf(MPI_Comm comm) : cOmm(comm) {}
    virtual int sync() {
      if (!this->str().empty()) {
        PetscFPrintf(cOmm, LogManager::dummy_mofem_fd, "%s",
                     this->str().c_str());
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
        PetscSynchronizedFPrintf(cOmm, LogManager::dummy_mofem_fd, "%s",
                                 this->str().c_str());
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

  static bool logQuiet;
  static std::map<std::string, LoggerType> logChannels;

  boost::shared_ptr<std::ostream> getStrmSelf() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmSelf);
  }
  boost::shared_ptr<std::ostream> getStrmWorld() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmWorld);
  }
  boost::shared_ptr<std::ostream> getStrmSync() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &strmSync);
  }

  InternalData(MPI_Comm comm)
      : worldBuf(comm), syncBuf(comm), strmSelf(&selfBuf), strmWorld(&worldBuf),
        strmSync(&syncBuf) {}

  virtual ~InternalData() = default;
};

bool LogManager::InternalData::logQuiet = false;
std::map<std::string, LogManager::LoggerType>
    LogManager::InternalData::logChannels;

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
  PetscInt sev_level = SeverityLevel::inform;
  PetscBool log_scope = PETSC_FALSE;
  PetscBool log_quiet = PETSC_FALSE;

  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "log_",
                           "Logging interface options", "none");

  CHKERR PetscOptionsEList("-severity_level", "Scope level", "",
                           severityStrings.data(), SeverityLevel::error + 1,
                           severityStrings[sev_level], &sev_level, PETSC_NULL);

  CHKERR PetscOptionsBool("-scope", "Log scope", "", log_scope, &log_scope,
                          NULL);

  CHKERR PetscOptionsBool("-quiet", "Quiet log attributes", "", log_quiet,
                          &log_quiet, NULL);

  ierr = PetscOptionsEnd();
  CHKERRG(ierr);

  CHKERR setUpLog();

  logging::core::get()->set_filter(MoFEM::LogKeywords::severity >= sev_level);

  if (log_scope)
    logging::core::get()->add_global_attribute("Scope", attrs::named_scope());

  if(log_quiet)
    LogManager::InternalData::logQuiet = true;

  MoFEMFunctionReturn(0);
}

void LogManager::recordFormatterDefault(logging::record_view const &rec,
                                        logging::formatting_ostream &strm) {


  if (!LogManager::InternalData::logQuiet) {

    auto sev = rec[severity];
    auto p = rec[proc_attr];
    auto l = rec[line_id];
    auto s = rec[scope];
    auto tg = rec[tag_attr];
    auto tl = rec[timeline];

    auto set_color = [&](const auto str) {
#if defined(PETSC_HAVE_UNISTD_H) && defined(PETSC_USE_ISATTY)
      if (isatty(fileno(stdout)))
        strm << str;
#endif
    };

    if (!p.empty()) {
      strm << "[";
      set_color("\033[32m");
      strm << p;
      set_color("\033[0m");
      strm << "] ";
    }

    if (sev > SeverityLevel::inform) {
      set_color("\033[31m");
      if (sev > SeverityLevel::warning)
        set_color("\033[1m");
    } else
      set_color("\033[34m");

    strm << sev;

    set_color("\033[0m");

    if (!l.empty())
      strm << std::hex << std::setw(8) << std::setfill('0') << l.get()
           << std::dec << std::setfill(' ') << ": ";

    if (!s.empty()) {
      for (::boost::log::attributes::named_scope_list::const_iterator iter =
               s->begin();
           iter != s->end(); ++iter) {
        const auto path = std::string(iter->file_name.data());
        const auto file = path.substr(path.find_last_of("/\\") + 1);
        strm << "(" << file << ":" << iter->line << ">" << iter->scope_name
             << ")";
      }
      strm << " ";
    }

    if (!tg.empty()) {

      set_color("\033[1m");
      strm << "[" << tg.get() << "] ";
      set_color("\033[0m");
    }

    if (!tl.empty())
      strm << "[" << tl.get() << "] ";
  }

  auto msg = rec[logging::expressions::smessage];

  strm << msg;
}

boost::shared_ptr<LogManager::SinkType>
LogManager::createSink(boost::shared_ptr<std::ostream> stream_ptr,
                       std::string comm_filter) {

  auto backend = boost::make_shared<sinks::text_ostream_backend>();
  if (stream_ptr)
    backend->add_stream(stream_ptr);
  backend->auto_flush(true);

  auto sink = boost::make_shared<SinkType>(backend);
  sink->set_filter((expr::has_attr(channel) && channel == comm_filter));
  sink->set_formatter(&recordFormatterDefault);

  return sink;
}

static char dummy_file;
FILE *LogManager::dummy_mofem_fd = (FILE *)&dummy_file;

MoFEMErrorCode LogManager::setUpLog() {
  MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  auto core_log = logging::core::get();
  core_log->add_sink(createSink(internalDataPtr->getStrmSelf(), "SELF"));
  core_log->add_sink(createSink(internalDataPtr->getStrmWorld(), "WORLD"));
  core_log->add_sink(createSink(internalDataPtr->getStrmSync(), "SYNC"));

  MoFEMFunctionReturn(0);
}

void LogManager::addAttributes(LogManager::LoggerType &lg, const int bit) {

  if (bit == 0)
    return;

  if (bit & (BitLineID | BitScope)) {

    if (bit & BitLineID)
      lg.add_attribute("LineID", attrs::counter<unsigned int>(1));

    if (bit & BitScope)
      lg.add_attribute("Scope", attrs::named_scope());

  } else {
    THROW_MESSAGE("Wrong cast");
  }
}

void LogManager::addAttributes(const std::string channel, const int bit) {
  addAttributes(getLog(channel), bit);
}

void LogManager::addTag(LogManager::LoggerType &lg, const std::string tag) {
  lg.add_attribute("Tag", attrs::constant<std::string>(tag));
}

void LogManager::addTag(const std::string channel, const std::string tag) {
  getLog(channel).add_attribute("Tag", attrs::constant<std::string>(tag));
}

LogManager::LoggerType &LogManager::setLog(const std::string channel) {
  InternalData::logChannels[channel] =
      LoggerType(boost::log::keywords::channel = channel);
  return InternalData::logChannels[channel];
}

LogManager::LoggerType &LogManager::getLog(const std::string channel) {
  return InternalData::logChannels.at(channel);
}

PetscErrorCode LogManager::logPetscFPrintf(FILE *fd, const char format[],
                                           va_list Argp) {
  PetscErrorCode ierr;
  PetscFunctionBegin;
  if (fd != stdout && fd != stderr && fd != dummy_mofem_fd) {
    ierr = PetscVFPrintfDefault(fd, format, Argp);
    CHKERR(ierr);

  } else {
    std::array<char, 1024> buff;
    size_t length;
    ierr = PetscVSNPrintf(buff.data(), 1024, format, &length, Argp);
    CHKERRQ(ierr);

    auto remove_line_break = [](auto &&msg) {
      if (!msg.empty() && msg.back() == '\n')
        msg = std::string_view(msg.data(), msg.size() - 1);
      return msg;
    };

    std::ostringstream ss;

    const std::string str(buff.data());
    if (!str.empty()) {
      if (fd != dummy_mofem_fd) {
        MOFEM_LOG("PETSC", MoFEM::LogManager::SeverityLevel::inform)
            << ss.str() << remove_line_break(std::string(buff.data()));
      } else {
        std::clog << ss.str() << std::string(buff.data());
      }
    }
  }
  PetscFunctionReturn(0);
}

std::string LogManager::getVLikeFormatedString(const char *fmt, va_list args) {
  std::array<char, 1024> buf;
  vsprintf(buf.data(), fmt, args);
  return std::string(buf.data());
}

std::string LogManager::getCLikeFormatedString(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  auto str = getVLikeFormatedString(fmt, args);
  va_end(args);
  return str;
}

} // namespace MoFEM
