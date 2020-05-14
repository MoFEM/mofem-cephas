/**
 * @file log.cpp
 * @example log.cpp
 * @brief Example and test how to log
 *
 * This is an example of how to use the logger.
 *
 */

#include <MoFEM.hpp>

#include <thread>
#include <chrono>

using namespace MoFEM;

MoFEMErrorCode log_fun1() {
  MoFEMFunctionBegin;

  MOFEM_LOG_CHANNEL("WORLD");
  // Log time
  BOOST_LOG_SCOPED_THREAD_ATTR("Timeline", attrs::timer());
  // Log tag
  MOFEM_LOG_TAG("WORLD", "Tag this output");

  MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose) << "Hello, world!";

  // sleep for half a second
  std::this_thread::sleep_for(std::chrono::milliseconds(300));

  MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose)
      << "Hello, second time world!";

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode log_fun3() {
  MoFEMFunctionBegin;
  MOFEM_LOG_FUNCTION();

  // Log scope
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose) << "Hello, sync!";
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose) << "Hello again, sync!";

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode log_fun2() {
  MoFEMFunctionBegin;

  // Log scope
  MOFEM_LOG_FUNCTION();
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_ATTRIBUTES("SYNC", LogManager::BitLineID | LogManager::BitScope);
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose) << "Hello, sync!";
  MOFEM_LOG("SYNC", LogManager::SeverityLevel::verbose) << "Hello again, sync!";

  CHKERR log_fun3();

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

    CHKERR PetscPrintf(PETSC_COMM_WORLD,
                       "Testing logging for obsolete way of printing messages");

    // Set "WORLD channel" 
    MOFEM_LOG_CHANNEL("WORLD");
    {
      MOFEM_LOG("WORLD", LogManager::SeverityLevel::error)
          << "Hello, self error!";
      MOFEM_LOG("WORLD", LogManager::SeverityLevel::warning)
          << "Hello, self warning!";
      MOFEM_LOG("WORLD", LogManager::SeverityLevel::inform)
          << "Hello, self inform!";
      MOFEM_LOG("WORLD", LogManager::SeverityLevel::verbose)
          << "Hello, self verbose!";
      MOFEM_LOG("WORLD", LogManager::SeverityLevel::noisy)
          << "Hello, self noisy!";
    }

    {
      MOFEM_C_LOG("WORLD", LogManager::SeverityLevel::inform, "%s %d %d %d",
                  "Hello C, self error!", 1, 2, 3);
    }

    {
      CHKERR log_fun1();
      CHKERR log_fun2();
      MOFEM_LOG_SYNCHORMISE(m_field.get_comm());
    }


    // SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "Trigger error");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}