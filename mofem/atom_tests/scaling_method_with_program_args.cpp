/** \file scaling_method.cpp

  \brief Testing interface for reading and writing CSV files containing time
  series data.

*/

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  const char delimiter = ',';
  std::vector<double> scalar_values = {1.1, 2.4, 3.6,  4.1,  3.1,
                                       5.1, 9.1, 10.5, 11.2, 15.3};
  try {
    auto time_scale = std::make_shared<TimeScale>();
    for (int i = 1; i <= scalar_values.size(); i++) {
      if (time_scale->getScale(double(i)) != scalar_values[i - 1]) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Validation for data scaling from csv "
                 "failed for time: %d value: %d",
                 double(i), time_scale->getScale(i));
      }
    }
    double time1 = 3.0;
    double time0 = 2.0;
    double scale1 = scalar_values[2];
    double scale0 = scalar_values[1];
    double input_time = 2.5;
    double interp_t = (input_time - time0) / (time1 - time0);
    double expected_scale = scale0 + (scale1 - scale0) * interp_t;
    double actual_scale = time_scale->getScale(2.5);
    if (expected_scale != actual_scale) {
      SETERRQ2(
          PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
          "Validation for data scaling from csv failed for time: %f value: %f",
          2.5, time_scale->getScale(2.5));
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
