/**
 * @file EssentialDisplacementCubitBcData.cpp
 * @brief Essential boundary conditions
 * @version 13.1
 * @date 2022-09-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

EssentialPreProc<DisplacementCubitBcData>::EssentialPreProc(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
    std::vector<boost::shared_ptr<ScalingMethod>> smv, bool get_coords)
    : mField(m_field), fePtr(fe_ptr), vecOfTimeScalingMethods(smv),
      getCoords(get_coords) {}

MoFEMErrorCode EssentialPreProc<DisplacementCubitBcData>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto fb = mField.getInterface<FieldBlas>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto disp_bc = bc.second->dispBcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          auto get_field_coeffs = [&](auto field_name) {
            auto field_ptr = mField.get_field_structure(field_name);
            return field_ptr->getNbOfCoeffs();
          };
          const auto nb_field_coeffs = get_field_coeffs(field_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply EssentialPreProc<DisplacementCubitBcData>: "
              << problem_name << "_" << field_name << "_" << block_name;

          double v;
          int coeff;
          std::array<std::vector<double>, 3> coords;
          int idx;

          auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
            MoFEMFunctionBegin;
            if (getCoords) {
              field_entity_ptr->getEntFieldData()[coeff] =
                  coords[coeff][idx] + v;
              ++idx;
            } else {
              field_entity_ptr->getEntFieldData()[coeff] = v;
            }
            MoFEMFunctionReturn(0);
          };

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
          if (getCoords) {
            for (auto d : {0, 1, 2})
              coords[d].resize(verts.size());
            CHKERR mField.get_moab().get_coords(verts, &*coords[0].begin(),
                                                &*coords[1].begin(),
                                                &*coords[2].begin());
          }

          if (disp_bc->data.flag1) {
            v = disp_bc->data.value1;
            for (auto s : vecOfTimeScalingMethods) {
              v *= s->getScale(fe_method_ptr->ts_t);
            }
            idx = 0;
            coeff = 0;
            CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
          } else if (disp_bc->data.flag2 && nb_field_coeffs > 1) {
            v = disp_bc->data.value2;
            for (auto s : vecOfTimeScalingMethods) {
              v *= s->getScale(fe_method_ptr->ts_t);
            }
            idx = 0;
            coeff = 1;
            CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
          } else if (disp_bc->data.flag3 && nb_field_coeffs > 2) {
            v = disp_bc->data.value3;
            for (auto s : vecOfTimeScalingMethods) {
              v *= s->getScale(fe_method_ptr->ts_t);
            }
            idx = 0;
            coeff = 2;
            CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
          }
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