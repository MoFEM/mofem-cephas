/** \file scaling_method.cpp

  \brief Testing interface for reading and writing CSV files containing time
  series data.

*/

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";
double lerp(double a, double b, double t) { return a + t * (b - a); }
int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  const char delimiter = ',';
  std::string fileName = "scalar_data.csv";
  std::vector<double> scalarValues = {1.1, 2.4, 3.6,  4.1,  3.1,
                                      5.1, 9.1, 10.5, 11.2, 15.3};
  try {
    auto timeScale = std::make_shared<TimeScale>();
    for (int i = 1; i <= scalarValues.size(); i++) {
      if (timeScale->getScale(double(i)) != scalarValues[i - 1]) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                 "Validation for data scaling from csv "
                 "failed for time: %d value: %d",
                 double(i), timeScale->getScale(i));
      }
    }
    double time1 = 3.0;
    double time0 = 2.0;
    double scale1 = scalarValues[2];
    double scale0 = scalarValues[1];
    double input_time = 2.5;
    double interp_t = (input_time - time0) / (time1 - time0);
    double expected_scale = scale0 + (scale1 - scale0) * interp_t;
    double actual_scale = timeScale->getScale(2.5);
    if (expected_scale != actual_scale) {
      SETERRQ2(
          PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
          "Validation for data scaling from csv failed for time: %d value: %d",
          2.5, timeScale->getScale(2.5));
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
