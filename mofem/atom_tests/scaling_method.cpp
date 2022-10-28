/** \file scaling_method.cpp

  \brief Testing interface for reading and writing CSV files containing time series data.

*/

#include <MoFEM.hpp>
#include <iostream>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {
    auto scaling_method = boost::make_shared<TimeScale>("scalar_data.csv", ',');
    if(scaling_method->getScale(2.0) != 2.4)
    {
      std::cout << scaling_method->getScale(2.0) << std::endl;
      // SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
      //   "*** Invalid scalar value for time < %d >", 2.0);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
