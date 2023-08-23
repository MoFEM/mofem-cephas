/**
 * @file DGProjection.cpp
 */

namespace MoFEM {

OpAliceProjector::OpAliceProjector(
    int order, boost::shared_ptr<MatrixDouble> mass_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s, int verb, Sev sev)
    : OpBaseDerivativesMass<1>(mass_ptr, data_l2, b, s, verb, sev),
      baseOrder(order) {
  getOrder = [this](ForcesAndSourcesCore *fe_ptr) { return baseOrder; };
}

OpAliceMapping::OpAliceMapping(boost::shared_ptr<MatrixDouble> data_ptr,
                               boost::shared_ptr<MatrixDouble> coeffs_ptr,
                               boost::shared_ptr<MatrixDouble> mass_ptr,
                               boost::shared_ptr<EntitiesFieldData> data_l2,
                               const FieldApproximationBase b,
                               const FieldSpace s,
                               const LogManager::SeverityLevel sev)
    : OpAliceProjector(0, mass_ptr, data_l2, b, s, VERBOSE, sev),
      dataPtr(data_ptr), coeffsPtr(coeffs_ptr) {}

MoFEMErrorCode OpAliceMapping::doWork(int side, EntityType type,
                                      EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  auto nb_integration_pts = getGaussPts().size2();
  auto t_w = getFTensor0IntegrationWeight();

  const auto fe_type = getFEType();
  auto &ent_data = dataL2->dataOnEntities[fe_type][0];
  auto nb_base_functions = ent_data.getN(base).size2();
  auto t_base = ent_data.getFTensor0N(base, 0, 0);
  auto nb_data_coeffs = dataPtr->size2();

  coeffsPtr->resize(nb_data_coeffs, nb_base_functions, false);

  for (auto gg = 0; gg != nb_integration_pts; ++gg) {
    for (auto bb = 0; bb != nb_base_functions; ++bb) {
      double alpha = t_w * t_base;
      auto t_f =
          FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(&(*coeffsPtr)(bb, 0));
      FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_data(
          &(*dataPtr)(gg, 0));
      for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
        t_f += alpha * t_data;
        ++t_f;
        ++t_data;
      }
      ++t_base;
    }
    ++t_w;
  }

  for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
    ublas::matrix_row<MatrixDouble> mc(*coeffsPtr, cc);
    cholesky_solve(*baseMassPtr, mc, ublas::lower());
  }

  MoFEMFunctionReturn(0);
}

OpBobMapping::OpBobMapping(boost::shared_ptr<MatrixDouble> data_ptr,
                           boost::shared_ptr<MatrixDouble> coeffs_ptr,
                           boost::shared_ptr<MatrixDouble> mass_ptr,
                           boost::shared_ptr<EntitiesFieldData> data_l2,
                           const FieldApproximationBase b, const FieldSpace s,
                           const LogManager::SeverityLevel sev)
    : OpAliceMapping(data_ptr, coeffs_ptr, mass_ptr, data_l2, b, s, sev) {}

MoFEMErrorCode OpBobMapping::doWork(int side, EntityType type,
                                    EntitiesFieldData::EntData &data) {

  MoFEMFunctionBegin;

  if (sPace != L2) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Space should be set to L2");
  }

  auto fe_type = getFEType();
  auto &ent_data = dataL2->dataOnEntities[fe_type][0];
  // baseOrder = ent_data.getOrder(base);
  CHKERR calculateBase(fe_type);

  auto &base_functions = ent_data.getN(base);
  const auto nb_base_functions = base_functions.size2();
  auto nb_data_coeffs = dataPtr->size2();

  if (nb_base_functions) {
    auto nb_integration_pts = getGaussPts().size2();
    dataPtr->resize(nb_integration_pts, nb_data_coeffs);
    auto t_base = ent_data.getFTensor0N(base);
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      for (int bb = 0; bb != nb_base_functions; ++bb) {
        auto t_coeffs = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
            &(*coeffsPtr)(bb, 0));
        FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_data(
            &(*dataPtr)(gg, 0));
        for (auto cc = 0; cc != nb_data_coeffs; ++cc) {
          t_data += t_base * t_coeffs;
          ++t_coeffs;
          ++t_data;
        }
        ++t_base;
      }
    }
  }
  MoFEMFunctionReturn(0);
}
} // namespace MoFEM