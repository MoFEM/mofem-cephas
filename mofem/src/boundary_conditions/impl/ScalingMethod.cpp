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

MoFEMErrorCode
ScalingMethod::applyScalingSeries(const FEMethod *fe,
                                  boost::ptr_vector<ScalingMethod> &methods_op,
                                  VectorDouble &nf) {
  MoFEMFunctionBegin;
  for (auto &v : methods_op) {
    CHKERR v.scaleNf(fe, nf);
  }
  MoFEMFunctionReturn(0);
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

/**
 * @brief Scale force the right hand vector
 *
 * @param fe
 * @param Nf
 * @return MoFEMErrorCode
 */
MoFEMErrorCode TimeScale::scaleNf(const FEMethod *fe, VectorDouble &Nf) {
  MoFEMFunctionBegin;
  Nf *= getScale(fe->ts_t);
  MoFEMFunctionReturn(0);
}

TimeVector::TimeVector(std::string name)
    : readFile(0), debug(0), nAme(name) {
  CHK_THROW_MESSAGE(timeData(), "Error reading time data");
}

MoFEMErrorCode TimeVector::timeData() {
  MoFEMFunctionBegin;
  char time_file_name[255];
  PetscBool flg = PETSC_TRUE;
  CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, nAme.c_str(),
                               time_file_name, 255, &flg);
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
  VectorDouble no2(3);
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
    std::map<double, VectorDouble>::iterator tit = tSeries.begin();
    for (; tit != tSeries.end(); tit++) {
      PetscPrintf(PETSC_COMM_WORLD,
                  "** read accelerogram %3.2e time %3.2e %3.2e %3.2e\n",
                  tit->first, tit->second[0], tit->second[1], tit->second[2]);
    }
  }
  if (r != 0) {
    SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR file close unsuccessful");
  }
  readFile = 1;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode TimeVector::scaleNf(const FEMethod *fe, VectorDouble &Nf) {
  MoFEMFunctionBegin;
  if (readFile == 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "Acceleration data file not read");
  }

  double ts_t = fe->ts_t;
  VectorDouble acc(3);
  VectorDouble acc0 = tSeries[0], acc1(3);
  double t0 = 0, t1, dt;
  for (auto &t : tSeries) {
    if (t.first > ts_t) {
      t1 = t.first;
      acc1 = t.second;
      dt = ts_t - t0;
      acc = acc0 + ((acc1 - acc0) / (t1 - t0)) * dt;
      break;
    }
    t0 = t.first;
    acc0 = t.second;
    acc = acc0;
  }
  Nf += acc;

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
