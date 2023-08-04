/**
 * @file EssentialMPCsData.cpp
 * @brief MPC boundary conditions
 * @version 13.1
 * @date 2022-09-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

MoFEMErrorCode setMPCParentAdjacency() {
  MoFEMFunctionBegin;

  MoFEMFunctionReturn(0);
}


EssentialPreProc<MPCsType>::EssentialPreProc(MoFEM::Interface &m_field,
                                   boost::shared_ptr<FEMethod> fe_ptr,
                                   bool is_spatial_positions)
    : m_field(m_field), fePtr(fe_ptr), isSpatialPositions(is_spatial_positions) {
}

MoFEMErrorCode EssentialPreProc<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = m_field.getInterface<BcManager>();
    // auto fb = m_field.getInterface<FieldBlas>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto disp_bc = bc.second->dispBcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          auto get_field_coeffs = [&](auto field_name) {
            auto field_ptr = m_field.get_field_structure(field_name);
            return field_ptr->getNbOfCoeffs();
          };
          // const auto nb_field_coeffs = get_field_coeffs(field_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply MultiPointConstraints PreProc<MPCsType>: "
              << problem_name << "_" << field_name << "_" << block_name;
          //TODO: here will we insert BcManager functionality
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