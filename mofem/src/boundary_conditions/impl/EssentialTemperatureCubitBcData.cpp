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

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply EssentialPreProc<TemperatureCubitBcData>: "
              << problem_name << "_" << field_name << "_" << block_name;

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
          auto v = temp_bc->data.value1;
          for (auto s : vecOfTimeScalingMethods) {
            v *= s->getScale(fe_method_ptr->ts_t);
          }

          //   std::array<std::vector<double>, 3> coords;
          //   auto verts_check = bc.second->bcEnts.subset_by_type(MBVERTEX);
            
          //  for (auto d : {0, 1, 2})
          //     coords[d].resize(verts_check.size());

          //   CHKERR mField.get_moab().get_coords(verts_check, &*coords[0].begin(),
          //                                       &*coords[1].begin(),
          //                                       &*coords[2].begin());
          // CHKERR PetscPrintf(PETSC_COMM_WORLD, "coords z: ");                                                                                                
          //   for (int ii = 0; ii != verts_check.size();  ++ii)
          //     CHKERR PetscPrintf(PETSC_COMM_WORLD, " %e ", coords[2][ii]);
              
          //   CHKERR PetscPrintf(PETSC_COMM_WORLD, "\n");   

          auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
            MoFEMFunctionBegin;
            for (auto &d : field_entity_ptr->getEntFieldData())
              d = v;
            MoFEMFunctionReturn(0);
          };
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