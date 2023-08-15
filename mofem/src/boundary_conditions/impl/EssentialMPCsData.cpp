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

// for coupling
// remove slave nodes from problem, set lambda funcs
// nothing for master nodes
// for rigid body
// set MASTER to no field, remove from problem, assign BCs (later)
// for tie
// basically the same as for rigid body, but dont calculate the distance

EssentialPreProc<MPCsType>::EssentialPreProc(MoFEM::Interface &m_field,
                                             boost::shared_ptr<FEMethod> fe_ptr,
                                             bool is_spatial_positions)
    : mField(m_field), fePtr(fe_ptr), isSpatialPositions(is_spatial_positions) {
}

MoFEMErrorCode EssentialPreProc<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (isSpatialPositions) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  MOFEM_LOG("BcMngWorld", Sev::warning)
      << "EssentialPreProc<MPCsType> has no effect here. ";

  // auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
  //   MoFEMFunctionBegin;
  //   auto nb_coeffs = field_entity_ptr->getNbDofsOnEnt();
  //   for (int i = 0; i != nb_coeffs; i++)
  //     field_entity_ptr->getEntFieldData()[i] = t_vals(i);
  //   MoFEMFunctionReturn(0);
  // };

  // if (!mpc_bc->isMaster) {
  //   auto master_ptr = mpc_bc->mpcMasterPtr;
  //   Range master_ents = *(mpc_bc->mpcMasterEntsPtr);
  //   Range slave_ents = bc.second->bcEnts;

  //   CHKERR fb->fieldLambdaOnEntities(lambda_master_val, field_name,
  //                                    &master_ents);
  //   t_vals(I) /= master_ents.size();
  //   CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &slave_ents);
  // }
  // FIXME: handle the cases for different directions
  // mpc_bc->data.flag1;

  MoFEMFunctionReturn(0);
}

EssentialPreProcLhs<MPCsType>::EssentialPreProcLhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Mat> lhs, SmartPetscObj<AO> ao)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vLhs(lhs), vAO(ao) {}

MoFEMErrorCode EssentialPreProcLhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    // auto fb = mField.getInterface<FieldBlas>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          // auto get_field_coeffs = [&](auto field_name) {
          //   auto field_ptr = mField.get_field_structure(field_name);
          //   return field_ptr->getNbOfCoeffs();
          // };
          // const auto nb_field_coeffs = get_field_coeffs(field_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply MultiPointConstraints PreProc<MPCsType>: "
              << problem_name << "_" << field_name << "_" << block_name;

          auto mpc_type = mpc_bc->mpcType;
          switch (mpc_type) {
          case MPC::COUPLING:
          case MPC::TIE:
          case MPC::RIGID_BODY:
            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "MPC type not implemented");
          }

          // FIXME: there handle the vertices differently (block level)
          auto master_verts = mpc_bc->mpcMasterEnts.subset_by_type(MBVERTEX);
          auto slave_verts = mpc_bc->mpcSlaveEnts.subset_by_type(MBVERTEX);

          // FIXME: for now only ALL dofs
          auto prb_name = fe_method_ptr->problemPtr->getName();

          SmartPetscObj<IS> is_xyz_m[2];
          SmartPetscObj<IS> is_xyz_s[2];
          RowColData rc[] = {ROW, COL};
          for (int i = 0; i != 2; i++) {
            CHKERR is_mng->isCreateProblemFieldAndRank(
                prb_name, rc[i], field_name, 0, MAX_DOFS_ON_ENTITY, is_xyz_m[i],
                &master_verts);
            CHKERR is_mng->isCreateProblemFieldAndRank(
                prb_name, rc[i], field_name, 0, MAX_DOFS_ON_ENTITY, is_xyz_s[i],
                &slave_verts);
          }
          // PetscInt size_check;
          // ISGetSize(is_xyz_m[0], P&size_check);

          // MOFEM_LOG("WORLD", Sev::noisy)
          //     << "is_xyz_m[0] size: " << size_check;

          if (is_xyz_m[0] && is_xyz_s[0]) {

            if (auto fe_ptr = fePtr.lock()) {
              SmartPetscObj<Mat> B =
                  vLhs ? vLhs : SmartPetscObj<Mat>(fe_ptr->B, true);
              // The user is responsible for assembly if vLhs is provided
              if (fe_ptr->matAssembleSwitch && !vLhs) {
                if (*fe_ptr->matAssembleSwitch) {
                  CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                  CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
                  *fe_ptr->matAssembleSwitch = false;
                }
              }

              // if (vAO) {
              //   MOFEM_LOG("WORLD", Sev::noisy) << "Apply AO to IS";
              //   CHKERR AOApplicationToPetscIS(vAO, is_xyz_m[0]);
              //   CHKERR AOApplicationToPetscIS(vAO, is_xyz_s[0]);
              // }

              // CHKERR MatSetOption(B, MAT_USE_HASH_TABLE, PETSC_TRUE);
              CHKERR MatSetOption(B, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);

              // CHKERR MatZeroRowsIS(B, is_xyz_m[0], 0, PETSC_NULL, PETSC_NULL);
              CHKERR MatZeroRowsIS(B, is_xyz_s[0], 0, PETSC_NULL, PETSC_NULL);

              // CHKERR MatSetValuesIS(B, is_xyz_m[0], is_xyz_m[1],
              // &master_val,
              //                       INSERT_VALUES);
              // CHKERR MatSetValuesIS(B, is_xyz_s[0], is_xyz_s[1],
              // &slave_val,
              //                       INSERT_VALUES);

              auto set_mat_values = [&](auto row_is, auto col_is, double val) {
                // we need to use MatSetValues here
                MoFEMFunctionBeginHot;
                const int *row_index_ptr;
                CHKERR ISGetIndices(row_is, &row_index_ptr);
                const int *col_index_ptr;
                CHKERR ISGetIndices(col_is, &col_index_ptr);
                int size;
                CHKERR ISGetLocalSize(row_is, &size);
                MatrixDouble locMat(size, size);
                fill(locMat.data().begin(), locMat.data().end(), 0.0);
                for (int i = 0; i != size; i++)
                  locMat(i, i) = val;

                CHKERR ::MatSetValues(B, size, row_index_ptr, size,
                                      col_index_ptr, &*locMat.data().begin(),
                                      INSERT_VALUES);

                CHKERR ISRestoreIndices(row_is, &row_index_ptr);
                CHKERR ISRestoreIndices(col_is, &col_index_ptr);
                // MatView(Aij, PETSC_VIEWER_STDOUT_SELF);

                MoFEMFunctionReturnHot(0);
              };
              // CHKERR set_mat_values(is_xyz_m[0], is_xyz_m[1], 1);
              CHKERR set_mat_values(is_xyz_s[0], is_xyz_s[1], vDiag);
              CHKERR set_mat_values(is_xyz_s[0], is_xyz_m[1], -vDiag);
            };

            // User is responsible for assembly if vLhs is provided

            // ISView(is_sum, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERR MatZeroRowsColumnsIS(B, is_sum, vDiag, PETSC_NULL,
            //                             PETSC_NULL);
          } else {
            MOFEM_LOG("WORLD", Sev::error)
                << "Cannot create ISs for MPCs for field: " << field_name;
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

EssentialPreProcRhs<MPCsType>::EssentialPreProcRhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Vec> rhs)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vRhs(rhs) {}

MoFEMErrorCode EssentialPreProcRhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    auto vec_mng = mField.getInterface<VecManager>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          // auto get_field_coeffs = [&](auto field_name) {
            // auto field_ptr = mField.get_field_structure(field_name);
            // return field_ptr->getNbOfCoeffs();
          // };

          // const auto nb_field_coeffs = get_field_coeffs(field_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply MultiPointConstraints PreProc<MPCsType>: "
              << problem_name << "_" << field_name << "_" << block_name;

          auto mpc_type = mpc_bc->mpcType;
          switch (mpc_type) {
          case MPC::COUPLING:
          case MPC::TIE:
          case MPC::RIGID_BODY:
            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "MPC type not implemented");
          }

          // FIXME: there handle the vertices differently (block level)
          auto master_verts = mpc_bc->mpcMasterEnts.subset_by_type(MBVERTEX);
          auto slave_verts = mpc_bc->mpcSlaveEnts.subset_by_type(MBVERTEX);

          // FIXME: for now only ALL dofs
          auto prb_name = fe_method_ptr->problemPtr->getName();

          SmartPetscObj<IS> is_xyz_m;
          SmartPetscObj<IS> is_xyz_s;
          CHKERR is_mng->isCreateProblemFieldAndRank(prb_name, ROW, field_name,
                                                     0, MAX_DOFS_ON_ENTITY,
                                                     is_xyz_m, &master_verts);
          CHKERR is_mng->isCreateProblemFieldAndRank(prb_name, ROW, field_name,
                                                     0, MAX_DOFS_ON_ENTITY,
                                                     is_xyz_s, &slave_verts);

          // PetscInt size_check;
          // ISGetSize(is_xyz_m[0], P&size_check);

          // MOFEM_LOG("WORLD", Sev::noisy)
          //     << "is_xyz_m[0] size: " << size_check;

          // ISView(is_xyz_m[0], PETSC_VIEWER_STDOUT_WORLD);
          // ISView(is_xyz_s[0], PETSC_VIEWER_STDOUT_WORLD);
          // ISView(is_xyz_m[1], PETSC_VIEWER_STDOUT_WORLD);
          // ISView(is_xyz_s[1], PETSC_VIEWER_STDOUT_WORLD);

          if (is_xyz_m && is_xyz_s) {

            if (auto fe_ptr = fePtr.lock()) {
              auto snes_ctx = fe_ptr->snes_ctx;
              auto ts_ctx = fe_ptr->ts_ctx;
              const bool is_nonlinear = snes_ctx != FEMethod::CTX_SNESNONE ||
                                         ts_ctx != FEMethod::CTX_TSNONE;
              SmartPetscObj<Vec> F =
                  vRhs ? vRhs : SmartPetscObj<Vec>(fe_ptr->f, true);

              if (fe_ptr->vecAssembleSwitch && !vRhs) {
                CHKERR VecGhostUpdateBegin(F, ADD_VALUES, SCATTER_REVERSE);
                CHKERR VecGhostUpdateEnd(F, ADD_VALUES, SCATTER_REVERSE);
                CHKERR VecAssemblyBegin(F);
                CHKERR VecAssemblyEnd(F);
                *fe_ptr->vecAssembleSwitch = false;
              }

              auto set_vec_values = [&]() {

                MoFEMFunctionBeginHot;
                const int *m_index_ptr;
                CHKERR ISGetIndices(is_xyz_m, &m_index_ptr);
                const int *s_index_ptr;
                CHKERR ISGetIndices(is_xyz_s, &s_index_ptr);
                int size;
                CHKERR ISGetLocalSize(is_xyz_m, &size);
                // for testing
                int size2;
                CHKERR ISGetLocalSize(is_xyz_s, &size2);
                double *f;
                CHKERR VecGetArray(F, &f);

                if (is_nonlinear) {
                  auto x = fe_ptr->x;
                  auto tmp_x = vectorDuplicate(F);
                  
                  CHKERR vec_mng->setLocalGhostVector(
                      problem_name, ROW, tmp_x, INSERT_VALUES, SCATTER_FORWARD);
                  const double *v;
                  const double *u;
                  CHKERR VecGetArrayRead(tmp_x, &u);
                  CHKERR VecGetArrayRead(x, &v);

                  for (auto i = 0; i != size; ++i) {
                    f[s_index_ptr[i]] =
                        vDiag * (v[s_index_ptr[i]] - v[m_index_ptr[i]]);
                  }

                  CHKERR VecRestoreArrayRead(x, &v);
                  CHKERR VecRestoreArrayRead(tmp_x, &u);
                } else {
                  for (auto i = 0; i != size; ++i) {
                    f[s_index_ptr[i]] = 0;
                  }
                }

                CHKERR VecRestoreArray(F, &f);
                CHKERR ISRestoreIndices(is_xyz_m, &m_index_ptr);
                CHKERR ISRestoreIndices(is_xyz_s, &s_index_ptr);

                MoFEMFunctionReturnHot(0);
              };

              CHKERR set_vec_values();
          
            };

            // User is responsible for assembly if vLhs is provided

            // ISView(is_sum, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERR MatZeroRowsColumnsIS(B, is_sum, vDiag, PETSC_NULL,
            //                             PETSC_NULL);
          } else {
            MOFEM_LOG("WORLD", Sev::error)
                << "Cannot create ISs for MPCs for field: " << field_name;
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