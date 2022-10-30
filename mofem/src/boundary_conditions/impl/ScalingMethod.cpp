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
  TimeScale(file_name, defaultDelimiter, error_if_file_not_given) {
  
}


TimeScale::TimeScale(string file_name, char delimiter, bool error_if_file_not_given)
    : debug(0), fileName(file_name), delimiter(delimiter),
      errorIfFileNotGiven(error_if_file_not_given) {
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
  if(!in_file_stream.is_open() && errorIfFileNotGiven) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          "*** ERROR data file < %s > open unsuccessful", fileName);
  }
  else if(!in_file_stream.is_open() && !errorIfFileNotGiven) {
      MOFEM_LOG("WORLD", Sev::warning) << "*** Warning dadta file "<< fileName << " open unsuccessful. Using linear time scaling." ;
        scalingMethod = [this](double time){ return this->getLinearScale(time);};
  }
  in_file_stream.seekg(0);
  std::string line;
  double time = 0.0, value = 0.0;
  tSeries[time] = value;
  char * p_delimiter = &delimiter;
  std::regex rgx(static_cast<const char*>(p_delimiter));
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
  scalingMethod = [this](double time){ return this->getScaleFromData(time);};
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
  MoFEMFunctionReturn(0);
}

double TimeScale::getLinearScale(const double time) {
  return time;
}

double TimeScale::getScaleFromData(const double time) {
  double scale = 0;
  auto it = tSeries.find(time);
  if(it != tSeries.end()){
    scale = it->second;
  } else {
    auto upper = tSeries.upper_bound(time);
    auto lower = tSeries.lower_bound(time);
    double t = (time - lower->first)/(upper->first - lower->first);
    double scale1 = upper->second;
    double scale0 = lower->second;
    scale = lerp(scale0, scale1, t);
  }
  return scale;
}

double TimeScale::getScale(const double time) {
  return scalingMethod(time);
}

double TimeScale::lerp(const double a, const double b, const double t) {
  return a + t*(b-a);
}
} // namespace MoFEM
