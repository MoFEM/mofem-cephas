/** \file NormsOperators.cpp

\brief User data operators for calculating norms and differences between
fields

*/

namespace MoFEM {

OpCalcNormL2Tensor0::OpCalcNormL2Tensor0(
    boost::shared_ptr<VectorDouble> data_ptr, SmartPetscObj<Vec> data_vec,
    const int index, boost::shared_ptr<VectorDouble> diff_data_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index),
      diffDataPtr(diff_data_ptr) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
  if (!diffDataPtr)
    diffDataPtr = dataPtr;
}

MoFEMErrorCode OpCalcNormL2Tensor0::doWork(int side, EntityType type,
                                           EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // calculate the difference between data pointers and save them in diffDataPtr
  if (dataPtr != diffDataPtr)
    *diffDataPtr -= *dataPtr;

  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get values
  auto t_data = getFTensor0FromVec(*diffDataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // add to element norm
    norm_on_element += t_w * t_data * t_data;
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }
  // scale with volume of the element
  norm_on_element *= vol;
  // add to dataVec at iNdex position
  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

template <int DIM>
OpCalcNormL2Tensor1<DIM>::OpCalcNormL2Tensor1(
    boost::shared_ptr<MatrixDouble> data_ptr, SmartPetscObj<Vec> data_vec,
    const int index, boost::shared_ptr<MatrixDouble> diff_data_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index),
      diffDataPtr(diff_data_ptr) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
  if (!diffDataPtr)
    diffDataPtr = dataPtr;
}

template <int DIM>
MoFEMErrorCode
OpCalcNormL2Tensor1<DIM>::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // calculate the difference between data pointers and save them in diffDataPtr
  if (dataPtr != diffDataPtr)
    *diffDataPtr -= *dataPtr;

  // Declare FTensor index
  FTensor::Index<'i', DIM> i;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get vector values
  auto t_data = getFTensor1FromMat<DIM>(*diffDataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; ++gg) {
    // add to element norm
    norm_on_element += t_w * (t_data(i) * t_data(i));
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }
  // scale with volume of the element
  norm_on_element *= vol;
  // add to dataVec at iNdex position
  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

template <int DIM_1, int DIM_2>
OpCalcNormL2Tensor2<DIM_1, DIM_2>::OpCalcNormL2Tensor2(
    boost::shared_ptr<MatrixDouble> data_ptr, SmartPetscObj<Vec> data_vec,
    const int index, boost::shared_ptr<MatrixDouble> diff_data_ptr)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index),
      diffDataPtr(diff_data_ptr) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
  if (!diffDataPtr)
    diffDataPtr = dataPtr;
}

template <int DIM_1, int DIM_2>
MoFEMErrorCode
OpCalcNormL2Tensor2<DIM_1, DIM_2>::doWork(int side, EntityType type,
                                          EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // calculate the difference between data pointers and save them in diffDataPtr
  if (dataPtr != diffDataPtr)
    *diffDataPtr -= *dataPtr;

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
  auto t_data = getFTensor2FromMat<DIM_1, DIM_2>(*diffDataPtr);
  // initialise double to store norm values
  double norm_on_element = 0.;
  // loop over integration points
  for (int gg = 0; gg != nb_integration_points; gg++) {
    // add to element norm
    norm_on_element += t_w * (t_data(i, j) * t_data(i, j));
    // move to another integration weight
    ++t_w;
    // move to another data values
    ++t_data;
  }
  // scale with volume of the element
  norm_on_element *= vol;
  // add to dataVec at iNdex position
  CHKERR VecSetValue(dataVec, iNdex, norm_on_element, ADD_VALUES);

  MoFEMFunctionReturn(0);
}

OpGetTensor0fromFunc::OpGetTensor0fromFunc(
    boost::shared_ptr<VectorDouble> data_ptr, ScalarFunc scalar_function)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      dataPtr(data_ptr), sFunc(scalar_function) {}

MoFEMErrorCode OpGetTensor0fromFunc::doWork(int side, EntityType type,
                                            EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // get number of integration points
  const int nb_integration_pts = getGaussPts().size2();
  // resize dataPtr to store values at integration points
  dataPtr->resize(nb_integration_pts);
  // get values
  auto t_val = getFTensor0FromVec(*dataPtr);

  // get coordinates at integration point
  auto t_coords = getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != nb_integration_pts; ++gg) {
    // set value at integration point to value of scalar function at the
    // coordinates
    t_val = sFunc(t_coords(0), t_coords(1), t_coords(2));
    // move to another value
    ++t_val;
    // move to another coordinates
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

template <int SPACE_DIM, int BASE_DIM>
OpGetTensor1fromFunc<SPACE_DIM, BASE_DIM>::OpGetTensor1fromFunc(
    boost::shared_ptr<MatrixDouble> data_ptr, VectorFunc vector_function)
    : ForcesAndSourcesCore::UserDataOperator(
          NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
      dataPtr(data_ptr), vFunc(vector_function) {}

template <int SPACE_DIM, int BASE_DIM>
MoFEMErrorCode OpGetTensor1fromFunc<SPACE_DIM, BASE_DIM>::doWork(
    int side, EntityType type, EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index set to space dimension
  FTensor::Index<'i', SPACE_DIM> i;
  // get number of integration points
  const int nb_integration_pts = getGaussPts().size2();
  // resize dataPtr to store values at integration points
  dataPtr->resize(BASE_DIM, nb_integration_pts);
  // get values
  auto t_val = getFTensor1FromMat<BASE_DIM>(*(dataPtr));

  // get coordinates at integration point
  auto t_coords = getFTensor1CoordsAtGaussPts();
  // loop over integration points
  for (int gg = 0; gg != nb_integration_pts; ++gg) {
    // get function values of vector function at the coordinates
    auto func_val = vFunc(t_coords(0), t_coords(1), t_coords(2));
    // translate function values to tensor
    auto t_func_val = getFTensor1FromArray<SPACE_DIM, SPACE_DIM>(func_val);
    // set values at integration point to function values at the coordinates
    t_val(i) = t_func_val(i);

    // move to another value
    ++t_val;
    // move to another coordinates
    ++t_coords;
  }

  MoFEMFunctionReturn(0);
}

template struct OpCalcNormL2Tensor1<2>;
template struct OpCalcNormL2Tensor1<3>;
template struct OpCalcNormL2Tensor2<2, 2>;
template struct OpCalcNormL2Tensor2<3, 3>;
template struct OpGetTensor1fromFunc<2, 2>;
template struct OpGetTensor1fromFunc<2, 3>;
template struct OpGetTensor1fromFunc<3, 3>;

} // namespace MoFEM
