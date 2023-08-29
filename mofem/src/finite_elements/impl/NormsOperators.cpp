// /** \file NormsOperators.cpp

// \brief Generic user data operators for evaluate fields, and other common
// purposes.

// */

template struct MoFEM::OpCalcNormL2Tesnosr1<2>;
template struct MoFEM::OpCalcNormL2Tesnosr1<3>;

namespace MoFEM {

template <int FIELD_DIM>
OpCalcNormL2Tesnosr1<FIELD_DIM>::OpCalcNormL2Tesnosr1(
    const std::string field_name, boost::shared_ptr<MatrixDouble> data_ptr,
    SmartPetscObj<Vec> data_vec, const int index)
    : ForcesAndSourcesCore::UserDataOperator(
          field_name, ForcesAndSourcesCore::UserDataOperator::OPROW),
      dataPtr(data_ptr), dataVec(data_vec), iNdex(index) {
  if (!dataPtr)
    THROW_MESSAGE("Pointer is not set");
}

template <int FIELD_DIM>
MoFEMErrorCode
OpCalcNormL2Tesnosr1<FIELD_DIM>::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  // Declare FTensor index
  FTensor::Index<'i', FIELD_DIM> i;
  // get number of integration points
  const auto nb_integration_points = getGaussPts().size2();
  // get element volume
  const double vol = getMeasure();
  // get integration weights
  auto t_w = getFTensor0IntegrationWeight();
  // get vector values
  auto t_data = getFTensor1FromMat<FIELD_DIM>(*dataPtr);
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

} // namespace MoFEM
