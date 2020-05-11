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

MoFEMErrorCode log_fun1(MoFEM::Interface &m_field) {
  MoFEMFunctionBegin;

  MOFEM_LOG_CHANNEL("WORLD");
  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  MOFEM_LOG_TAG("WORLD", "Tag this output");

  MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose) << "Hello, world!";

  // sleep for half a second
  std::this_thread::sleep_for(std::chrono::milliseconds(300));

  MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose)
      << "Hello, second time world!";

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode log_fun2(MoFEM::Interface &m_field) {
  MoFEMFunctionBegin;

  MOFEM_LOG_FUNCTION();
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_ATTRIBUTES("SYNC", LogManager::BitLineID | LogManager::BitScope);
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::warning) << "Hello, sync!";
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::warning) << "Hello again, sync!";
  
  MoFEMFunctionReturn(0);
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

    MOFEM_LOG_CHANNEL("SELF");
    {
      MOFEM_LOG("SELF", LogManager::SeverityLevel::critical)
          << "Hello, self critical!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::error)
          << "Hello, self error!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::warning)
          << "Hello, self warning!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::inform)
          << "Hello, self inform!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::verbose)
          << "Hello, self verbose!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::noisy)
          << "Hello, self noisy!";
      MOFEM_LOG("SELF", LogManager::SeverityLevel::very_noisy)
          << "Hello, self very noisy!";
    }

    MOFEM_LOG_ATTRIBUTES("SELF", LogManager::BitScope);
    {
      MOFEM_LOG("SELF", LogManager::SeverityLevel::noisy)
          << "Hello, self with scope!";
    }

    {
      CHKERR log_fun1(m_field);
      CHKERR log_fun2(m_field);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}