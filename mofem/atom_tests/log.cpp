

#include <MoFEM.hpp>
#include <LogManager.hpp>

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
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace attrs = boost::log::attributes;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;

using namespace MoFEM;

void log_fun1(MoFEM::Interface &m_field) {
  BOOST_LOG_FUNCTION();
  auto world_log = m_field.getInterface<LogManager>()->getLogWorld();
  BOOST_LOG_SEV(world_log, LogManager::SeverityLevel::verbose)
      << "Hello, world!";
}

void log_fun2(MoFEM::Interface &m_field) {
  BOOST_LOG_FUNCTION();
  auto sync_log = m_field.getInterface<LogManager>()->getLogSync();
  BOOST_LOG_SEV(sync_log, LogManager::SeverityLevel::warning) << "Hello, sync!";
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    logging::core::get()->set_filter(MoFEM::LogKeywords::severity >=
                                     LogManager::SeverityLevel::noisy);

    auto self_log = m_field.getInterface<LogManager>()->getLogSelf();

    {
      BOOST_LOG_NAMED_SCOPE("log test");
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::noisy)
          << "Hello, self!";
    }

    {
      BOOST_LOG_NAMED_SCOPE("more log test")
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::verbose)
          << "Hello, self!";
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::very_noisy)
          << "Hello, self!";
    }

    {
      BOOST_LOG_NAMED_SCOPE("test functions")
      log_fun1(m_field);
      log_fun2(m_field);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}