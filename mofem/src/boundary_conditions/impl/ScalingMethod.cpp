/**
 * @file ScalingMethod.cpp
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

TimeScale::TimeScale(std::string file_name, bool error_if_file_not_given,
                     ScalingFun def_scaling_fun)
    : TimeScale(file_name, defaultDelimiter, error_if_file_not_given,
                def_scaling_fun) {}

TimeScale::TimeScale(std::string file_name, std::string delimiter,
                     bool error_if_file_not_given, ScalingFun def_scaling_fun)
    : fileName(file_name), errorIfFileNotGiven(error_if_file_not_given),
      defScalingMethod(def_scaling_fun) {
  CHK_THROW_MESSAGE(timeData(file_name, delimiter),
                    "Error in reading time data");
}

MoFEMErrorCode TimeScale::timeData(std::string fileName,
                                   std::string delimiter) {
  MoFEMFunctionBegin;
  MOFEM_LOG_CHANNEL("WORLD");
  // Set the argument found flag to false as default
  char time_file_name[255] = {'\0'};
  // Check to see if a filename has been provided
  if (fileName.empty()) {
    // If no filename, look for command line argument
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, fileNameFlag.c_str(),
                                 time_file_name, 255, &argFileScale);
    if (argFileScale) {
      fileName = std::string(time_file_name);
    }
  } else {
    // Set the command line flag to true for correct flow control using provided
    // filename
    argFileScale = PETSC_TRUE;
  }
  if (!argFileScale && fileName.empty() && errorIfFileNotGiven) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR %s (time_data FILE NEEDED)", fileName.c_str());
  } else if (!argFileScale && fileName.empty() && !errorIfFileNotGiven) {
    MOFEM_LOG_C("WORLD", Sev::warning,
                "The %s file not provided. Loading scaled with time.",
                fileName.c_str());
    scalingMethod = [this](double time) {
      return this->defScalingMethod(time);
    };
    MoFEMFunctionReturnHot(0);
  }

  std::ifstream in_file_stream(fileName);
  if (!in_file_stream.is_open() && errorIfFileNotGiven) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", fileName.c_str());
  } else if (!in_file_stream.is_open() && !errorIfFileNotGiven) {
    MOFEM_LOG("WORLD", Sev::warning)
        << "*** Warning data file " << fileName
        << " open unsuccessful. Using default time scaling.";
    scalingMethod = [this](double time) {
      return this->defScalingMethod(time);
    };
    MoFEMFunctionReturnHot(0);
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

double TimeScale::getScale(const double time) { return scalingMethod(time); }

// Vector scale

// template <int SPACE_DIM> TimeScaleVector<SPACE_DIM>::TimeScaleVector() {}
template <int SPACE_DIM>
TimeScaleVector<SPACE_DIM>::TimeScaleVector(string name,
                                            bool error_if_file_not_given)
    : readFile(0), debug(0), nAme(name),
      errorIfFileNotGiven(error_if_file_not_given) {
  CHK_THROW_MESSAGE(timeData(), "Error in reading time data");
}

template <int SPACE_DIM>
TimeScaleVector<SPACE_DIM>::TimeScaleVector(std::string name, int ms_id,
                                            bool error_if_file_not_given)
    : readFile(0), debug(0), errorIfFileNotGiven(error_if_file_not_given) {
  nAme = name + std::to_string(ms_id);
  CHK_THROW_MESSAGE(timeData(), "Error in reading time data");
}

template <int SPACE_DIM> MoFEMErrorCode TimeScaleVector<SPACE_DIM>::timeData() {
  MoFEMFunctionBegin;

  char time_file_name[255];
  PetscBool flg = PETSC_FALSE;
  if (!nAme.empty())
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, nAme.c_str(),
                                 time_file_name, 255, &flg);

  if (!flg) {
    if (errorIfFileNotGiven)
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "*** ERROR %s (time_data FILE NEEDED)", nAme.c_str());
    MoFEMFunctionReturnHot(0);
  }

  FILE *time_data = fopen(time_file_name, "r");
  if (time_data == NULL) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "*** ERROR data file < %s > open unsuccessful", time_file_name);
  }
  double no1 = 0.0;
  FTensor::Index<'i', 3> i;
  (tSeries[no1])(i) = 0.;
  while (!feof(time_data)) {
    FTensor::Tensor1<double, 3> no2{0., 0., 0.};
    int n =
        fscanf(time_data, "%lf %lf %lf %lf", &no1, &no2(0), &no2(1), &no2(2));
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
    (tSeries[no1])(i) = no2(i);
  }
  int r = fclose(time_data);

  MOFEM_LOG_CHANNEL("WORLD");
  for (auto &[ts, vec] : tSeries) {
    MOFEM_TAG_AND_LOG_C("WORLD", Sev::verbose, "TimeScaleVector",
                        "** read vector %3.2e time %3.2e %3.2e %3.2e", ts,
                        vec(0), vec(1), vec(2));
  }
  MOFEM_LOG_CHANNEL("WORLD");

  if (r != 0) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
            "*** ERROR file close unsuccessful");
  }
  readFile = 1;

  if (readFile == 1)
    scalingMethod = [this](double time) {
      return this->getVectorFromData(time);
    };

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM>
FTensor::Tensor1<double, SPACE_DIM>
TimeScaleVector<SPACE_DIM>::getVector(const double time) {
  return scalingMethod(time);
}

template <int SPACE_DIM>
FTensor::Tensor1<double, SPACE_DIM>
TimeScaleVector<SPACE_DIM>::getVectorFromData(const double time) {

  FTensor::Tensor1<double, SPACE_DIM> Nf;
  FTensor::Tensor1<double, 3> acc;
  FTensor::Tensor1<double, 3> acc0 = tSeries.begin()->second;

  FTensor::Index<'I', SPACE_DIM> I;
  FTensor::Index<'i', SPACE_DIM> i;

  double t0 = 0, t1, dt;
  for (auto &[ts, vec] : tSeries) {
    if (ts > time) {
      t1 = ts;
      dt = time - t0;
      acc(I) = acc0(I) + ((vec(I) - acc0(I)) / (t1 - t0)) * dt;
      break;
    }
    t0 = ts;
    acc0 = vec;
    acc = acc0;
  }
  Nf(i) = acc(i);
  return Nf;
}

template class TimeScaleVector<3>;
template class TimeScaleVector<2>;

} // namespace MoFEM
