// /** \file NormsOperators.cpp

// \brief Generic user data operators for evaluate fields, and other common
// purposes.

// */

template struct MoFEM::OpCalcNormL2Tensor1<2>;
template struct MoFEM::OpCalcNormL2Tensor1<3>;
template struct MoFEM::OpCalcNormL2Tensor2<2, 2>;
template struct MoFEM::OpCalcNormL2Tensor2<3, 3>;
template struct MoFEM::OpCalcDifferenceTensor1<2>;
template struct MoFEM::OpCalcDifferenceTensor1<3>;
template struct MoFEM::OpCalcDifferenceTensor2<2, 2>;
template struct MoFEM::OpCalcDifferenceTensor2<3, 3>;

namespace MoFEM {

OpCalcNormL2Tensor0::OpCalcNormL2Tensor0(
    const std::string field_name, boost::shared_ptr<VectorDouble> data_ptr,
    SmartPetscObj<Vec> data_vec, const int index)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
}

MoFEMErrorCode OpCalcNormL2Tensor0::doWork(int side, EntityType type,
                                           EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get values
  auto t_data = getFTensor0FromVec(*dataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol;
    // add to element norm
    norm_on_element += alpha * t_data * t_data;
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }

  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

template <int DIM>
OpCalcNormL2Tensor1<DIM>::OpCalcNormL2Tensor1(
    const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
    SmartPetscObj<Vec> data_vec, const int index)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
}

template <int DIM>
MoFEMErrorCode
OpCalcNormL2Tensor1<DIM>::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index
  FTensor::Index<'i', DIM> i;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get vector values
  auto t_data = getFTensor1FromMat<DIM>(*dataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol;
    // add to element norm
    norm_on_element += alpha * t_data(i) * t_data(i);
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }

  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

template <int DIM_1, int DIM_2>
OpCalcNormL2Tensor2<DIM_1, DIM_2>::OpCalcNormL2Tensor2(
    const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
    SmartPetscObj<Vec> data_vec, const int index)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
}

template <int DIM_1, int DIM_2>
MoFEMErrorCode
OpCalcNormL2Tensor2<DIM_1, DIM_2>::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index
  FTensor::Index<'i', DIM_1> i;
  FTensor::Index<'j', DIM_2> j;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get vector values
  auto t_data = getFTensor2FromMat<DIM_1, DIM_2>(*dataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // take into account Jacobian
    const double alpha = t_w * vol;
    // add to element norm
    norm_on_element += alpha * t_data(i, j) * t_data(i, j);
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }

  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

OpCalcDifferenceTensor0::OpCalcDifferenceTensor0(
    const std::string field_name, boost::shared_ptr<VectorDouble> data_1_ptr,
    boost::shared_ptr<VectorDouble> data_2_ptr,
    boost::shared_ptr<VectorDouble> result_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      data1Ptr(data_1_ptr), data2Ptr(data_2_ptr), resultPtr(result_ptr) {
  if (!data1Ptr || !data2Ptr)
    THROW_MESSAGE("Pointer is not set");
}

MoFEMErrorCode
OpCalcDifferenceTensor0::doWork(int side, EntityType type,
                                EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // resize resulting pointer
  resultPtr->resize(nb_integration_points);
  // get values
  auto t_data_1 = getFTensor0FromVec(*data1Ptr);
  auto t_data_2 = getFTensor0FromVec(*data2Ptr);
  auto t_result = getFTensor0FromVec(*resultPtr);

  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // get difference
    t_result = t_data_1 - t_data_2;
    // move to another data values
    ++t_data_1;
    ++t_data_2;
    ++t_result;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM>
OpCalcDifferenceTensor1<DIM>::OpCalcDifferenceTensor1(
    const std::string field_name, boost::shared_ptr<MatrixDouble> data_1_ptr,
    boost::shared_ptr<MatrixDouble> data_2_ptr,
    boost::shared_ptr<MatrixDouble> result_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      data1Ptr(data_1_ptr), data2Ptr(data_2_ptr), resultPtr(result_ptr) {
  if (!data1Ptr || !data2Ptr)
    THROW_MESSAGE("Pointer is not set");
}

template <int DIM>
MoFEMErrorCode
OpCalcDifferenceTensor1<DIM>::doWork(int side, EntityType type,
                                     EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index
  FTensor::Index<'i', DIM> i;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // resize resulting pointer
  resultPtr->resize(DIM, nb_integration_points);
  // get values
  auto t_data_1 = getFTensor1FromMat<DIM>(*data1Ptr);
  auto t_data_2 = getFTensor1FromMat<DIM>(*data2Ptr);
  auto t_result = getFTensor1FromMat<DIM>(*resultPtr);

  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // get difference
    t_result(i) = t_data_1(i) - t_data_2(i);
    // move to another data values
    ++t_data_1;
    ++t_data_2;
    ++t_result;
  }

  MoFEMFunctionReturn(0);
}

template <int DIM_1, int DIM_2>
OpCalcDifferenceTensor2<DIM_1, DIM_2>::OpCalcDifferenceTensor2(
    const std::string field_name, boost::shared_ptr<MatrixDouble> data_1_ptr,
    boost::shared_ptr<MatrixDouble> data_2_ptr,
    boost::shared_ptr<MatrixDouble> result_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      data1Ptr(data_1_ptr), data2Ptr(data_2_ptr), resultPtr(result_ptr) {
  if (!data1Ptr || !data2Ptr)
    THROW_MESSAGE("Pointer is not set");
}

template <int DIM_1, int DIM_2>
MoFEMErrorCode OpCalcDifferenceTensor2<DIM_1, DIM_2>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index
  FTensor::Index<'i', DIM_1> i;
  FTensor::Index<'j', DIM_2> j;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // resize resulting pointer
  resultPtr->resize(DIM_1 * DIM_2, nb_integration_points);
  // get values
  auto t_data_1 = getFTensor2FromMat<DIM_1, DIM_2>(*data1Ptr);
  auto t_data_2 = getFTensor2FromMat<DIM_1, DIM_2>(*data2Ptr);
  auto t_result = getFTensor2FromMat<DIM_1, DIM_2>(*resultPtr);

  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // get difference
    t_result(i, j) = t_data_1(i, j) - t_data_2(i, j);
    // move to another data values
    ++t_data_1;
    ++t_data_2;
    ++t_result;
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
