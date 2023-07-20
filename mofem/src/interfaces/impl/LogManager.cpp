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

  class NulStreambuf : public std::streambuf {
    char dummyBuffer[64];

  protected:
    virtual int overflow(int c) {
      setp(dummyBuffer, dummyBuffer + sizeof(dummyBuffer));
      return (c == traits_type::eof()) ? '\0' : c;
    }
  };

  SelfStreamBuf selfBuf;
  WorldStreamBuf worldBuf;
  SynchronizedStreamBuf syncBuf;

  std::ostream strmSelf;
  std::ostream strmWorld;
  std::ostream strmSync;

  class NulOStream : private NulStreambuf, public std::ostream {
  public:
    NulOStream() : std::ostream(this) {}
    NulStreambuf *rdbuf() { return this; }
  };

  NulOStream nullSelf;

  static bool logQuiet;
  static bool noColors;
  static bool sinksAdd;
  static bool logTime;

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
  boost::shared_ptr<std::ostream> getStrmNull() {
    return boost::shared_ptr<std::ostream>(shared_from_this(), &nullSelf);
  }

  InternalData(MPI_Comm comm)
      : worldBuf(comm), syncBuf(comm), strmSelf(&selfBuf), strmWorld(&worldBuf),
        strmSync(&syncBuf) {}

  virtual ~InternalData() = default;
};

bool LogManager::InternalData::logQuiet = false;
bool LogManager::InternalData::noColors = false;
bool LogManager::InternalData::sinksAdd = true;
bool LogManager::InternalData::logTime = false;

std::map<std::string, LogManager::LoggerType>
    LogManager::InternalData::logChannels;

boost::shared_ptr<LogManager::InternalData> LogManager::internalDataPtr;

LogManager::LogManager(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)) {}

MoFEMErrorCode
LogManager::query_interface(boost::typeindex::type_index type_index,
                            UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<LogManager *>(this);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode LogManager::getOptions() {
  MoFEMFunctionBegin;
  PetscInt sev_level = SeverityLevel::inform;
  PetscBool log_scope = PETSC_FALSE;
  PetscBool log_quiet = PETSC_FALSE;
  PetscBool log_no_colors = PETSC_FALSE;
  PetscBool log_time = PETSC_FALSE;

  CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "log_",
                           "Logging interface options", "none");

  CHKERR PetscOptionsEList("-severity_level", "Severity level", "",
                           severityStrings.data(), SeverityLevel::error + 1,
                           severityStrings[sev_level], &sev_level, PETSC_NULL);
  CHKERR PetscOptionsEList("-sl", "Severity level", "",
                           severityStrings.data(), SeverityLevel::error + 1,
                           severityStrings[sev_level], &sev_level, PETSC_NULL);

  CHKERR PetscOptionsBool("-scope", "Log scope", "", log_scope, &log_scope,
                          NULL);

  CHKERR PetscOptionsBool("-quiet", "Quiet log attributes", "", log_quiet,
                          &log_quiet, NULL);

  CHKERR PetscOptionsBool("-no_color", "Terminal with no colors", "",
                          log_no_colors, &log_no_colors, NULL);
  CHKERR PetscOptionsBool("-nc", "Terminal with no colors", "", log_no_colors,
                          &log_no_colors, NULL);

  CHKERR PetscOptionsBool("-time", "Log time", "",
                          log_time, &log_time, NULL);

  ierr = PetscOptionsEnd();
  CHKERRG(ierr);

  logging::core::get()->set_filter(MoFEM::LogKeywords::severity >= sev_level);

  if (log_scope)
    logging::core::get()->add_global_attribute("Scope", attrs::named_scope());

  if (log_quiet)
    LogManager::InternalData::logQuiet = true;

  if (log_no_colors)
    LogManager::InternalData::noColors = true;

  if (log_time)
    LogManager::InternalData::logTime = true;

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
      if (isatty(fileno(stdout)) && !LogManager::InternalData::noColors)
        strm << str;
#endif
    };

    if(LogManager::InternalData::logTime) {

      auto local_time = boost::posix_time::second_clock::local_time();
      strm << "(Local time ";
      strm << local_time.date().year() << "-" << local_time.date().month()
           << "-" << local_time.date().day() << " "
           << local_time.time_of_day().hours() << ":"
           << local_time.time_of_day().minutes() << ":"
           << local_time.time_of_day().seconds();
      strm << ") ";

    }

    if (!p.empty()) {
      strm << "[";
      set_color("\033[32m");
      strm << p;
      set_color("\033[0m");
      strm << "] ";
    }

    switch (sev.get()) {
    case SeverityLevel::error:
      set_color("\033[1m");
    case SeverityLevel::warning:
      set_color("\033[31m");
      break;
    case SeverityLevel::inform:
      set_color("\033[34m");
      break;
    case SeverityLevel::verbose:
      set_color("\033[35m");
      break;
    case SeverityLevel::noisy:
      set_color("\033[36m");
      break;
    }

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

void LogManager::createDefaultSinks(MPI_Comm comm) {
  
  internalDataPtr = boost::make_shared<InternalData>(comm);

  auto core_log = logging::core::get();
  core_log->remove_all_sinks();
  core_log->add_sink(LogManager::createSink(
      boost::shared_ptr<std::ostream>(&std::clog, boost::null_deleter()),
      "PETSC"));
  core_log->add_sink(createSink(getStrmSelf(), "SELF"));
  core_log->add_sink(createSink(getStrmWorld(), "WORLD"));
  core_log->add_sink(createSink(getStrmSync(), "SYNC"));
  core_log->add_sink(createSink(getStrmNull(), "NULL"));

  LogManager::setLog("PETSC");
  LogManager::setLog("SELF");
  LogManager::setLog("WORLD");
  LogManager::setLog("SYNC");
  LogManager::setLog("NULL");

  MOFEM_LOG_TAG("PETSC", "petsc");

  int rank;
  MPI_Comm_rank(comm, &rank);
  core_log->add_global_attribute("Proc", attrs::constant<unsigned int>(rank));
}

boost::shared_ptr<std::ostream> LogManager::getStrmSelf() {
  return internalDataPtr->getStrmSelf();
}

boost::shared_ptr<std::ostream> LogManager::getStrmWorld() {
  return internalDataPtr->getStrmWorld();
}

boost::shared_ptr<std::ostream> LogManager::getStrmSync() {
  return internalDataPtr->getStrmSync();
}

boost::shared_ptr<std::ostream> LogManager::getStrmNull() {
  return internalDataPtr->getStrmNull();
}

static char dummy_file;
FILE *LogManager::dummy_mofem_fd = (FILE *)&dummy_file;

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

bool LogManager::checkIfChannelExist(const std::string channel) {
  return InternalData::logChannels.find(channel) !=
         InternalData::logChannels.end();
}

std::string LogManager::petscStringCache = std::string();

PetscErrorCode LogManager::logPetscFPrintf(FILE *fd, const char format[],
                                           va_list Argp) {
  MoFEMFunctionBegin;
  if (fd != stdout && fd != stderr && fd != dummy_mofem_fd) {
    CHKERR PetscVFPrintfDefault(fd, format, Argp);

  } else {

    std::array<char, 1024> buff;
    size_t length;
    CHKERR PetscVSNPrintf(buff.data(), 1024, format, &length, Argp);

    auto get_str = [&buff]() {
      std::string str;
      if (!petscStringCache.empty())
        str = petscStringCache + std::string(buff.data());
      else
        str = std::string(buff.data());
      return str;
    };
    const auto str = get_str();

    if (!str.empty()) {
      if (fd != dummy_mofem_fd) {
        
        MoFEM::LogManager::SeverityLevel sev =
            MoFEM::LogManager::SeverityLevel::inform;
        if (str.find("WARNING") != std::string::npos)
          sev = MoFEM::LogManager::SeverityLevel::warning;

        std::istringstream is(str);
        std::string line;
        std::vector<std::string> log_list;

        while (getline(is, line, '\n')) 
          log_list.push_back(line);

        if (str.back() != '\n') {
          petscStringCache = log_list.back();
          log_list.pop_back();
        } else
          petscStringCache.clear();

        for(auto &line : log_list) 
          MOFEM_LOG("PETSC", sev) << line;

      } else {
        std::clog << str;
      }
    }
  }
  MoFEMFunctionReturn(0);
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
