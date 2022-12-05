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

namespace MoFEM {

double ScalingMethod::getScale(const double time) {
  THROW_MESSAGE("getScale not implemented");
}

const std::string TimeScale::defaultDelimiter =
    "(\\s*,\\s*|\\s+)"; ///< comma or space

TimeScale::TimeScale(std::string file_name, bool error_if_file_not_given)
    : TimeScale(file_name, defaultDelimiter, error_if_file_not_given) {}

TimeScale::TimeScale(std::string file_name, std::string delimiter,
                     bool error_if_file_not_given)
    : fileName(file_name), errorIfFileNotGiven(error_if_file_not_given) {
  CHK_THROW_MESSAGE(timeData(delimiter), "Error in reading time data");
}

MoFEMErrorCode TimeScale::timeData(std::string delimiter) {
  MoFEMFunctionBegin;
  PetscBool arg_found = PETSC_FALSE;
  char time_file_name[255] = {'\0'};
  CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, fileNameFlag.c_str(),
                               time_file_name, 255, &arg_found);
  if (arg_found) {
    fileName = std::string(time_file_name);
  }
  if (!arg_found && fileName.empty() && errorIfFileNotGiven) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR %s (time_data FILE NEEDED)", fileName.c_str());
  } else if (!arg_found && fileName.empty() && !errorIfFileNotGiven) {
    MOFEM_LOG_C("WORLD", Sev::warning,
                "The %s file not provided. Loading scaled with time.",
                fileName.c_str());
    scalingMethod = [this](double time) { return this->getLinearScale(time); };
    MoFEMFunctionReturnHot(0);
  }

  std::ifstream in_file_stream(fileName);
  if (!in_file_stream.is_open() && errorIfFileNotGiven) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", fileName.c_str());
  } else if (!in_file_stream.is_open() && !errorIfFileNotGiven) {
    MOFEM_LOG("WORLD", Sev::warning)
        << "*** Warning data file " << fileName
        << " open unsuccessful. Using linear time scaling.";
    scalingMethod = [this](double time) { return this->getLinearScale(time); };
  }
  in_file_stream.seekg(0);
  std::string line;
  double time = 0.0, value = 0.0;
  tSeries[time] = value;

  std::regex rgx(delimiter.c_str());
  std::sregex_token_iterator end;
  while (std::getline(in_file_stream, line)) {
    std::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
    auto time_str = iter;
    auto value_str = ++iter;
    if (time_str == end || value_str == end) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "*** ERROR read data file error (check input time data file) ");
    }
    time = std::stod(time_str->str());
    value = std::stod(value_str->str());
    tSeries[time] = value;
  }
  in_file_stream.close();
  scalingMethod = [this](double time) { return this->getScaleFromData(time); };
  if (in_file_stream.is_open()) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "*** ERROR file close unsuccessful");
  }
  MoFEMFunctionReturn(0);
}

double TimeScale::getLinearScale(const double time) { return time; }

double TimeScale::getScaleFromData(const double time) {
  if (tSeries.empty())
    return 0.;

  auto it = tSeries.find(time);
  if (it != tSeries.end()) {
    return it->second;
  } else {
    auto it = tSeries.lower_bound(time);
    if (it == tSeries.end()) {
      return (--it)->second;
    }
    if (it == tSeries.begin()) {
      return it->second;
    }
    auto upper = *(it);
    it--;
    if (it == tSeries.end()) {
      return upper.second;
    }
    auto lower = *it;
    double t = (time - lower.first) / (upper.first - lower.first);
    double scale1 = upper.second;
    double scale0 = lower.second;
    return scale0 + t * (scale1 - scale0);
  }
}

TimeScaleVector::TimeScaleVector(string name, bool error_if_file_not_given)
    : readFile(0), debug(0), nAme(name),
      errorIfFileNotGiven(error_if_file_not_given) {
  CHK_THROW_MESSAGE(timeData(), "Error in reading time data");
}

MoFEMErrorCode TimeScaleVector::timeData() {

  MoFEMFunctionBeginHot;

  char time_file_name[255];
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, nAme.c_str(),
                               time_file_name, 255, &flg);
  CHKERRG(ierr);
  if (flg != PETSC_TRUE) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR %s (time_data FILE NEEDED)", nAme.c_str());
  }
  FILE *time_data = fopen(time_file_name, "r");
  if (time_data == NULL) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", time_file_name);
  }
  double no1 = 0.0;
  std::array<double, 3> no2{0, 0, 0};
  tSeries[no1] = no2;
  while (!feof(time_data)) {
    int n =
        fscanf(time_data, "%lf %lf %lf %lf", &no1, &no2[0], &no2[1], &no2[2]);
    if (n < 0) {
      fgetc(time_data);
      continue;
    }
    if (n != 4) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "*** ERROR read data file error (check input time data file) "
               "{ n = %d }",
               n);
    }
    tSeries[no1] = no2;
  }
  int r = fclose(time_data);
  if (debug) {

    for (auto &[ts, vec] : tSeries) {
      PetscPrintf(PETSC_COMM_WORLD,
                  "** read accelerogram %3.2e time %3.2e %3.2e %3.2e\n",
                  ts, vec[0], vec[1], vec[2]);
    }
  }
  if (r != 0) {
    SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR file close unsuccessful");
  }
  readFile = 1;

  MoFEMFunctionReturnHot(0);
}

std::array<double, 3> TimeScaleVector::getVector(const double time) {

  if (readFile == 0) {
    CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Data file not read");
  }

  std::array<double, 3> acc{0, 0, 0};
  std::array<double, 3> acc0 = tSeries.begin()->second;
  std::array<double, 3> Nf{0, 0, 0};

  double t0 = 0, t1, dt;
  for (auto &[ts, vec] : tSeries) {
    if (ts > time) {
      t1 = ts;
      dt = time - t0;
      for (int i = 0; i != 3; ++i)
        acc[i] = acc0[i] + ((vec[i] - acc0[i]) / (t1 - t0)) * dt;
      break;
    }
    t0 = ts;
    acc0 = vec;
    acc = acc0;
  }
  for (int i = 0; i != 3; ++i)
    Nf[i] += acc[i];
  return Nf;
}

} // namespace MoFEM
