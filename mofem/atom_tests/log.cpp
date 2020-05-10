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

  auto world_log = m_field.getInterface<LogManager>()->getLogWorld();
  // world_log.add_attribute("LineID", attrs::counter<unsigned int>(1));
  world_log.add_attribute("TimeStamp", attrs::local_clock());

  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  BOOST_LOG_FUNCTION();
  BOOST_LOG_SEV(world_log, LogManager::SeverityLevel::verbose)
      << "Hello, world!";
  std::this_thread::sleep_for(std::chrono::seconds(1));
  BOOST_LOG_SEV(world_log, LogManager::SeverityLevel::verbose)
      << "Hello, world!";
}

void log_fun2(MoFEM::Interface &m_field) {

  auto sync_log = m_field.getInterface<LogManager>()->getLogSync();
  sync_log.add_attribute("LineID", attrs::counter<unsigned int>(1));
  sync_log.add_attribute("Scope", attrs::named_scope());
  sync_log.add_attribute("TimeStamp", attrs::local_clock());

  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
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

    logging::core::get()->set_filter(MoFEM::LogKeywords::severity >=
                                     LogManager::SeverityLevel::very_noisy);

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