/**
 * @file matrix_function.cpp
 * @example matrix_function.cpp
 * @brief Test and example for matrix function
 *
 */

#define FTENSOR_DEBUG
#include <FTensor.hpp>
#include <MatrixFunction.hpp>

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
  LogManager::setLog("ATOM_TEST");

  try {
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
