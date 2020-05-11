/**
 * @file log.cpp
 * @brief Example and test how to log 
 * 
 * This is an example of how to use the logger.
 * 
 */

#include <MoFEM.hpp>

#include <thread>
#include <chrono> 

using namespace MoFEM;

void log_fun1(MoFEM::Interface &m_field) {

  auto world_log = LogManager::setLogWorld();
  LogManager::addTag(world_log, "My tag");

  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  BOOST_LOG_SEV(world_log, LogManager::SeverityLevel::verbose)
      << "Hello, world!";

  // sleep for half a second
  std::this_thread::sleep_for(std::chrono::milliseconds(300));

  BOOST_LOG_SEV(world_log, LogManager::SeverityLevel::verbose)
      << "Hello, second time world!";
}

void log_fun2(MoFEM::Interface &m_field) {

  auto sync_log =
      LogManager::setLogSync(LogManager::BitLineID | LogManager::BitScope);

  BOOST_LOG_FUNCTION();
  BOOST_LOG_SEV(sync_log, LogManager::SeverityLevel::warning) << "Hello, sync!";
  BOOST_LOG_SEV(sync_log, LogManager::SeverityLevel::warning)
      << "Hello again, sync!";
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);


  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    // logging::core::get()->set_filter(MoFEM::LogKeywords::severity >=
    //                                  LogManager::SeverityLevel::noisy);


    {
      auto self_log = LogManager::setLogSelf(LogManager::BitScope);
      BOOST_LOG_NAMED_SCOPE("log test");
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::noisy)
          << "Hello, self!";
    }

    {
      auto self_log = LogManager::setLogSelf();
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::verbose)
          << "Hello, self!";
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::noisy)
          << "Hello, self!";
      BOOST_LOG_SEV(self_log, LogManager::SeverityLevel::very_noisy)
          << "Hello, self!";
    }

    {
      auto self_inform = LogManager::setLogSelf();
      BOOST_LOG_SEV(self_inform, LogManager::SeverityLevel::inform)
          << "Hello, work inform";
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