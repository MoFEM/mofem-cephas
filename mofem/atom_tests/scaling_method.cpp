/** \file scaling_method.cpp

  \brief Testing interface for reading and writing CSV files containing time series data.

*/

#include <MoFEM.hpp>
#include <iostream>
using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  const char delimiter = ',';
  std::string fileName = "scalar_data.csv";
  std::vector<double> scalarValues = {1.1, 2.4, 3.6, 4.1, 3.1, 5.1, 9.1, 10.5, 11.2, 15.3};
  try {
    auto timeScale = std::make_shared<TimeScale>(fileName, delimiter);
    auto timeScaleLinearScaling = std::make_shared<TimeScale>();
    for(int i = 1; i <= scalarValues.size(); i++) {
      if(timeScale->getScale(double(i)) != scalarValues[i])
      {
        std::cout << timeScale->getScale(2.0) << std::endl;
        std::cout << "Error: Time: " <<  double(i) << " Value: " << timeScale->getScale(double(i)) << "\n";
      }
    }
    for(int i = 1; i <= scalarValues.size(); i++) {
      if(timeScaleLinearScaling->getScale(double(i)) != double(i))
      {
        std::cout << timeScaleLinearScaling->getScale(2.0) << std::endl;
        std::cout << "Error: Time: " <<  double(i) << " Value: " << timeScaleLinearScaling->getScale(double(i)) << "\n";
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
