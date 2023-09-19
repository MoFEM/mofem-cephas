/**
 * @file DGProjection.cpp
 */

namespace MoFEM {

OpDGProjectionMassMatrix::OpDGProjectionMassMatrix(
    int order, boost::shared_ptr<MatrixDouble> mass_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s, int verb, Sev sev)
    : OpBaseDerivativesBase(mass_ptr, data_l2, b, s, verb, sev),
      baseOrder(order) {}

MoFEMErrorCode
OpDGProjectionMassMatrix::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (sPace != L2) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Space should be set to L2");
  }

  CHKERR calculateBase(

      [this]() { return this->baseOrder; }

  );
  CHKERR calculateMass();

  MoFEMFunctionReturn(0);
}

OpDGProjectionCoefficients::OpDGProjectionCoefficients(
    boost::shared_ptr<MatrixDouble> data_ptr,
    boost::shared_ptr<MatrixDouble> coeffs_ptr,
    boost::shared_ptr<MatrixDouble> mass_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s,
    const LogManager::SeverityLevel sev)
    : OpBaseDerivativesBase(mass_ptr, data_l2, b, s, VERBOSE, sev),
      dataPtr(data_ptr), coeffsPtr(coeffs_ptr) {}

MoFEMErrorCode
OpDGProjectionCoefficients::doWork(int side, EntityType type,
                                   EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  auto nb_integration_pts = getGaussPts().size2();

  const auto fe_type = getFEType();
  constexpr auto side_number = 0;
  auto &ent_data = dataL2->dataOnEntities[fe_type][side_number];
  auto nb_base_functions = ent_data.getN(base).size2();

  auto nb_data_coeffs = dataPtr->size1();
  auto &rhs = *coeffsPtr;
  auto &data = *dataPtr;

  rhs.resize(nb_base_functions, nb_data_coeffs, false);
  rhs.clear();

  using Tensor0 = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>;
  MatrixDouble trans_data = trans(data);

  // assemble rhs
  auto t_base = ent_data.getFTensor0N(base, 0, 0);
  auto t_w = getFTensor0IntegrationWeight();
  for (auto gg = 0; gg != nb_integration_pts; ++gg) {
    Tensor0 t_rhs(&*rhs.data().begin());
    for (auto bb = 0; bb != nb_base_functions; ++bb) {
      double alpha = t_w * t_base;
      Tensor0 t_data(&trans_data(gg, 0));
      for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
        t_rhs += alpha * t_data;
        ++t_data;
        ++t_rhs;
      }
      ++t_base;
    }
    ++t_w;
  }

  // solve for coefficients
  for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
    ublas::matrix_column<MatrixDouble> mc(rhs, cc);
    cholesky_solve(*baseMassPtr, mc, ublas::lower());
  }

  MoFEMFunctionReturn(0);
}

OpDGProjectionEvaluation::OpDGProjectionEvaluation(
    boost::shared_ptr<MatrixDouble> data_ptr,
    boost::shared_ptr<MatrixDouble> coeffs_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s,
    const LogManager::SeverityLevel sev)
    : OpDGProjectionCoefficients(data_ptr, coeffs_ptr,
                                 boost::make_shared<MatrixDouble>(), data_l2, b,
                                 s, sev) {}

MoFEMErrorCode
OpDGProjectionEvaluation::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

  if (sPace != L2) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Space should be set to L2");
  }

  auto fe_type = getFEType();
  constexpr auto side_number = 0;
  auto order = dataL2->dataOnEntities[fe_type][side_number].getOrder();
  CHKERR calculateBase(

      [order]() { return order; }

  );

  auto &ent_data = dataL2->dataOnEntities[fe_type][side_number];
  const auto nb_base_functions = ent_data.getN(base).size2();
  auto nb_data_coeffs = dataPtr->size1();
  auto &data = *dataPtr;
  auto &coeffs = *coeffsPtr;

  auto nb_integration_pts = getGaussPts().size2();
  data.resize(nb_integration_pts, nb_data_coeffs);
  data.clear();

  using Tensor0 = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>;

  auto t_base = ent_data.getFTensor0N(base, 0, 0);
  for (auto gg = 0; gg != nb_integration_pts; ++gg) {
    Tensor0 t_coeffs(&*coeffs.data().begin());
    for (auto bb = 0; bb != nb_base_functions; ++bb) {
      double alpha = t_base;
      Tensor0 t_data(&data(gg, 0));
      for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
        t_data += alpha * t_coeffs;
        ++t_data;
        ++t_coeffs;
      }
      ++t_base;
    }
  }

  *dataPtr = trans(*dataPtr);  

  MoFEMFunctionReturn(0);
}
} // namespace MoFEM