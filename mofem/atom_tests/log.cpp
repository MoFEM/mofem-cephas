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

  MOFEM_LOG("WORLD", Sev::verbose) << "Hello, world!";

  // sleep for half a second
  std::this_thread::sleep_for(std::chrono::milliseconds(300));

  MOFEM_LOG("WORLD", Sev::verbose) << "Hello, second time world!";

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode log_fun3() {
  MoFEMFunctionBegin;
  MOFEM_LOG_FUNCTION();

  // Log scope
  MOFEM_LOG("SYNC", Sev::verbose) << "Hello, sync!";
  MOFEM_LOG("SYNC", Sev::verbose) << "Hello again, sync!";

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode log_fun2() {
  MoFEMFunctionBegin;

  // Log scope
  MOFEM_LOG_FUNCTION();
  MOFEM_LOG_CHANNEL("SYNC");
  MOFEM_LOG_ATTRIBUTES("SYNC", LogManager::BitLineID | LogManager::BitScope);
  MOFEM_LOG("SYNC", Sev::verbose) << "Hello, sync!";
  MOFEM_LOG("SYNC", Sev::verbose) << "Hello again, sync!";

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
                       "\nTesting logging for obsolete way of printing "
                       "messages\nnext line\nnext line\n\n");

    CHKERR PetscPrintf(PETSC_COMM_WORLD,
                       "Ala have ");
    CHKERR PetscPrintf(PETSC_COMM_WORLD, "a ");
    CHKERR PetscPrintf(PETSC_COMM_WORLD,
                       "cat\n");

    CHKERR PetscPrintf(PETSC_COMM_WORLD, "WARNING\n");

    // Set "WORLD channel"
    MOFEM_LOG_CHANNEL("WORLD");
    {
      MOFEM_LOG("WORLD", Sev::error) << "Hello, self error!";
      MOFEM_LOG("WORLD", Sev::warning) << "Hello, self warning!";
      MOFEM_LOG("WORLD", Sev::inform) << "Hello, self inform!";
      MOFEM_LOG("WORLD", Sev::verbose) << "Hello, self verbose!";
      MOFEM_LOG("WORLD", Sev::noisy) << "Hello, self noisy!";
    }

    {
      MOFEM_LOG_C("WORLD", Sev::inform, "%s %d %d %d", "Hello C, self error!",
                  1, 2, 3);
    }

    {
      CHKERR log_fun1();
      CHKERR log_fun2();
      MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
    }

    // Create channel
    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
    LogManager::setLog("ATOM_TEST");

    // add sink to channel
    core_log->add_sink(LogManager::createSink(
        boost::make_shared<std::ofstream>("log0.log"), "ATOM_TEST"));

    // add sink to channel other way
    logging::add_file_log(keywords::file_name = "log1.log",
                          keywords::filter =
                              MoFEM::LogKeywords::channel == "ATOM_TEST");

    // add skink to channel third way
    auto backend = boost::make_shared<sinks::text_ostream_backend>();
    backend->add_stream(boost::make_shared<std::ofstream>("log2.log"));
    auto sink = boost::make_shared<LogManager::SinkType>(backend);
    sink->set_filter((expr::has_attr(MoFEM::LogKeywords::channel) &&
                      MoFEM::LogKeywords::channel == "ATOM_TEST"));
    core_log->add_sink(sink);


    MOFEM_LOG_TAG("ATOM_TEST", "atom test");
    MOFEM_LOG("ATOM_TEST", Sev::inform) << "Test atom test channel";


    // SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY, "Trigger error");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
}