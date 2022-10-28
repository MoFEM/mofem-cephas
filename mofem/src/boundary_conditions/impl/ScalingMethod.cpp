/**
 * @file TimeScaling.cpp
 * @brief
 * @version 0.1
 * @date 2022-08-12
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>
#include <fstream>
#include <regex>
namespace MoFEM {

double ScalingMethod::getScale(const double time) {
  THROW_MESSAGE("getScale not implemented");
}

TimeScale::TimeScale(string file_name, bool error_if_file_not_given) 
  :
  TimeScale(file_name, ',', error_if_file_not_given) {
  
}


TimeScale::TimeScale(string file_name, char delimiter, bool error_if_file_not_given)
    : readFile(0), debug(0), fileName(file_name), delimiter(delimiter),
      errorIfFileNotGiven(error_if_file_not_given),fLg(PETSC_TRUE) {
  CHK_THROW_MESSAGE(timeData(), "Error in reading time data");
}

MoFEMErrorCode TimeScale::timeData() {
  MoFEMFunctionBegin;
  // char time_file_name[255];
  // CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, fileName.c_str(),
  //                              time_file_name, 255, &fLg);
  // if (!fLg && errorIfFileNotGiven) {
  //   SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
  //            "*** ERROR %s (time_data FILE NEEDED)", fileName.c_str());
  // }
  // if (!fLg) {
  //   MOFEM_LOG_C("WORLD", Sev::warning,
  //               "The %s file not provided. Loading scaled with time.",
  //               nAme.c_str());
  //   MoFEMFunctionReturnHot(0);
  // }
  std::ifstream in_file_stream(fileName);
  
  if(!in_file_stream.is_open()) {
        SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", fileName);
  }
  in_file_stream.seekg(0);
  std::string line;
  double time = 0.0, value = 0.0;
  tSeries[time] = value;
  std::regex rgx(",");
  std::sregex_token_iterator end;
  while (std::getline(in_file_stream, line)) {
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    auto time_str = iter;
    auto value_str = ++iter;
    if(time_str == end || value_str == end) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "*** ERROR read data file error (check input time data file) ");
    }
    time = std::stod(time_str->str());
    value = std::stod(value_str->str());
    tSeries[time] = value;
  }
  in_file_stream.close();
  if(in_file_stream.is_open()){
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "*** ERROR file close unsuccessful");
  }
  if (debug) {
    std::map<double, double>::iterator tit = tSeries.begin();
    for (; tit != tSeries.end(); tit++) {
      PetscPrintf(PETSC_COMM_WORLD, "** read time series %3.2e time %3.2e\n",
                  tit->first, tit->second);
    }
  }
  readFile = 1;
  MoFEMFunctionReturn(0);
}

double TimeScale::getScale(const double time) {
  if (readFile == 0) {
    CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Data file not read");
  }
  double scale = 0;
  auto it = tSeries.find(time);
  if(it != tSeries.end()){
    scale = it->second;
  } else {
    auto upper = tSeries.upper_bound(time);
    auto lower = tSeries.lower_bound(time);
    double time1 = upper->first;
    double scale1 = upper->second;
    double time0 = lower->first;
    double scale0 = lower->second;
    double dt = time - time0;
    scale = scale0 + ((scale1 - scale0)/(time1 - time0)) * dt;
  }
  return scale;
}
} // namespace MoFEM
