/**
 * @file EssentialTemperatureCubitBcData.cpp
 * @brief Essential boundary conditions
 * @version 13.1
 * @date 2022-09-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

EssentialPreProc<TemperatureCubitBcData>::EssentialPreProc(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
    std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : mField(m_field), fePtr(fe_ptr), vecOfTimeScalingMethods(smv) {}

MoFEMErrorCode EssentialPreProc<TemperatureCubitBcData>::operator()() {
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto fb = mField.getInterface<FieldBlas>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto temp_bc = bc.second->tempBcPtr) {

        auto &bc_id = bc.first;
        std::regex field_rgx("^(.*)_(.*)_(.*)$");
        std::smatch match_field_name;
        std::string field_name;
        std::string block_name;

        if (std::regex_search(bc_id, match_field_name, field_rgx)) {
          field_name = match_field_name[2];
          block_name = match_field_name[3];
        } else {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Field name and block name can not be resolved");
        }

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply EssentialPreProc<TemperatureCubitBcData>: "
              << problem_name << "_" << field_name << "_" << block_name;

          double v;

          auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
            MoFEMFunctionBegin;
            std::fill(field_entity_ptr->getEntFieldData().begin(),
                      field_entity_ptr->getEntFieldData().end(), v);
            MoFEMFunctionReturn(0);
          };

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
          v = temp_bc->data.value1;
          for (auto s : vecOfTimeScalingMethods) {
            v *= s->getScale(fe_method_ptr->ts_t);
          }
          CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
        }
      }
    }

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM