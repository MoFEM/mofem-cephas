/**
 * @file FormsBrokenSpaceConstraintImpl.cpp
 * @date 2024-07-04
 * 
 * @copyright Copyright (c) 2024
 * 
 */

namespace MoFEM {

OpGetBrokenBaseSideData::OpGetBrokenBaseSideData(
    const std::string field_name,
    boost::shared_ptr<BrokenBaseSideData> broken_base_side_data)
    : OP(field_name, field_name, OP::OPROW),
      brokenBaseSideData(broken_base_side_data) {}

MoFEMErrorCode
OpGetBrokenBaseSideData::doWork(int row_side, EntityType row_type,
                                EntitiesFieldData::EntData &row_data) {
  MoFEMFunctionBegin;

  brokenBaseSideData->getSide() = row_side;
  brokenBaseSideData->getType() = row_type;

  brokenBaseSideData->getData().sEnse = row_data.sEnse;
  brokenBaseSideData->getData().sPace = row_data.sPace;
  brokenBaseSideData->getData().bAse = row_data.bAse;
  brokenBaseSideData->getData().iNdices = row_data.iNdices;
  brokenBaseSideData->getData().localIndices = row_data.localIndices;
  brokenBaseSideData->getData().dOfs = row_data.dOfs;
  brokenBaseSideData->getData().fieldEntities = row_data.fieldEntities;
  brokenBaseSideData->getData().fieldData = row_data.fieldData;

  auto base = brokenBaseSideData->getData().bAse;
  for (auto dd = 0; dd != BaseDerivatives::LastDerivative; ++dd) {
    brokenBaseSideData->getData().baseFunctionsAndBaseDerivatives[dd][base] =
        row_data.baseFunctionsAndBaseDerivatives[dd][base];
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM