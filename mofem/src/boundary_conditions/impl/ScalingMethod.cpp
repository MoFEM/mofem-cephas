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

TimeScale::TimeScale(string name, bool error_if_file_not_given)
    : readFile(0), debug(0), nAme(name),
      errorIfFileNotGiven(error_if_file_not_given) {
  CHK_THROW_MESSAGE(timeData(), "Error in reading time data");
}

MoFEMErrorCode TimeScale::timeData() {
  MoFEMFunctionBegin;
  char time_file_name[255];
  CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, nAme.c_str(),
                               time_file_name, 255, &fLg);
  if (!fLg && errorIfFileNotGiven) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR %s (time_data FILE NEEDED)", nAme.c_str());
  }
  if (!fLg) {
    MOFEM_LOG_C("WORLD", Sev::warning,
                "The %s file not provided. Loading scaled with time.",
                nAme.c_str());
    MoFEMFunctionReturnHot(0);
  }
  FILE *time_data = fopen(time_file_name, "r");
  if (time_data == NULL) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", time_file_name);
  }
  double no1 = 0.0, no2 = 0.0;
  tSeries[no1] = no2;
  while (!feof(time_data)) {
    int n = fscanf(time_data, "%lf %lf", &no1, &no2);
    if ((n <= 0) || ((no1 == 0) && (no2 == 0))) {
      fgetc(time_data);
      continue;
    }
    if (n != 2) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "*** ERROR read data file error (check input time data file) "
               "{ n = %d }",
               n);
    }
    tSeries[no1] = no2;
  }
  int r = fclose(time_data);
  if (debug) {
    std::map<double, double>::iterator tit = tSeries.begin();
    for (; tit != tSeries.end(); tit++) {
      PetscPrintf(PETSC_COMM_WORLD, "** read time series %3.2e time %3.2e\n",
                  tit->first, tit->second);
    }
  }
  if (r != 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "*** ERROR file close unsuccessful");
  }
  readFile = 1;
  MoFEMFunctionReturn(0);
}

double TimeScale::getScale(const double time) {
  if (!fLg) {
    return time; // scale with time, by default
  }
  if (readFile == 0) {
    CHK_THROW_MESSAGE(MOFEM_OPERATION_UNSUCCESSFUL, "Data file not read");
  }

  double scale = 0;
  double t0 = 0, t1, s0 = tSeries[0], s1, dt;
  for (auto &t : tSeries) {
    if (t.first > time) {
      t1 = t.first;
      s1 = t.second;
      dt = time - t0;
      scale = s0 + ((s1 - s0) / (t1 - t0)) * dt;
      break;
    }
    t0 = t.first;
    s0 = t.second;
    scale = s0;
  }
  return scale;
}




} // namespace MoFEM
