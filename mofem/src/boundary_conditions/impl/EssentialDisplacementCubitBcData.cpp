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

          FTensor::Tensor1<double, 3> t_angles(0., 0., 0.);
          FTensor::Tensor1<double, 3> t_vals(0., 0., 0.);

          auto scale_value = [&](const double &c) {
            double val = c;
            for (auto s : vecOfTimeScalingMethods) {
              val *= s->getScale(fe_method_ptr->ts_t);
            }
            return val;
          };

          if (disp_bc->data.flag1 == 1)
            t_vals(0) = scale_value(-disp_bc->data.value1);
          if (disp_bc->data.flag2 == 1)
            t_vals(1) = scale_value(-disp_bc->data.value2);
          if (disp_bc->data.flag3 == 1)
            t_vals(2) = scale_value(-disp_bc->data.value3);
          if (disp_bc->data.flag4 == 1)
            t_angles(0) = scale_value(-disp_bc->data.value4);
          if (disp_bc->data.flag5 == 1)
            t_angles(1) = scale_value(-disp_bc->data.value5);
          if (disp_bc->data.flag6 == 1)
            t_angles(2) = scale_value(-disp_bc->data.value6);

          int coeff;
          std::array<std::vector<double>, 3> coords;
          int idx;

          const bool is_rotation = disp_bc->data.flag4 || disp_bc->data.flag5 || disp_bc->data.flag6;

          auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
            MoFEMFunctionBegin;
            auto v = t_vals(coeff);
            if (is_rotation) {
              FTensor::Tensor1<double, 3> t_coords(
                  coords[0][idx], coords[1][idx], coords[2][idx]);
              v += _getRotDisp(t_angles, t_coords)(coeff);
            }
            if (getCoords)
              v += coords[coeff][idx];

            field_entity_ptr->getEntFieldData()[coeff] = v;
            ++idx;

            MoFEMFunctionReturn(0);
          };

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
          if (getCoords || is_rotation) {
            for (auto d : {0, 1, 2})
              coords[d].resize(verts.size());
            CHKERR mField.get_moab().get_coords(verts, &*coords[0].begin(),
                                                &*coords[1].begin(),
                                                &*coords[2].begin());
          }

          if (disp_bc->data.flag1 || disp_bc->data.flag5 ||
              disp_bc->data.flag6) {
            idx = 0;
            coeff = 0;
            CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
          }
          if (disp_bc->data.flag2 || disp_bc->data.flag4 ||
              disp_bc->data.flag6 && nb_field_coeffs > 1) {
            idx = 0;
            coeff = 1;
            CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
          }
          if ((disp_bc->data.flag3 || disp_bc->data.flag4 ||
               disp_bc->data.flag5 && nb_field_coeffs > 2) ||
              is_rotation && nb_field_coeffs > 1) {
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