// /** \file NormsOperators.cpp

// \brief Generic user data operators for evaluate fields, and other common
// purposes.

// */

template struct MoFEM::OpCalcNormL2Tesnosr1<2>;
template struct MoFEM::OpCalcNormL2Tesnosr1<3>;
template struct MoFEM::OpCalcNormL2Tesnosr2<2, 2>;
template struct MoFEM::OpCalcNormL2Tesnosr2<3, 3>;

namespace MoFEM {

OpCalcNormL2Tesnosr0::OpCalcNormL2Tesnosr0(
    const std::string field_name, boost::shared_ptr<VectorDouble> data_ptr,
    SmartPetscObj<Vec> data_vec, const int index)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
}

MoFEMErrorCode OpCalcNormL2Tesnosr0::doWork(int side, EntityType type,
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
  double norm_on_element;
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
OpCalcNormL2Tesnosr1<DIM>::OpCalcNormL2Tesnosr1(
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
OpCalcNormL2Tesnosr1<DIM>::doWork(int side, EntityType type,
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
  double norm_on_element;
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
OpCalcNormL2Tesnosr2<DIM_1, DIM_2>::OpCalcNormL2Tesnosr2(
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
OpCalcNormL2Tesnosr2<DIM_1, DIM_2>::doWork(int side, EntityType type,
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
  double norm_on_element;
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

} // namespace MoFEM
