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

EssentialPostProcLhs<MPCsType>::EssentialPostProcLhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Mat> lhs, SmartPetscObj<AO> ao)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vLhs(lhs), vAO(ao) {}

MoFEMErrorCode EssentialPostProcLhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    const auto problem_name = fe_ptr->problemPtr->getName();

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

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
          constexpr auto max_nb_dofs_per_node = 6;

          if (nb_field_coeffs > max_nb_dofs_per_node)
            MOFEM_LOG("WORLD", Sev::error)
                << "MultiPointConstraints PreProcLhs<MPCsType>: support only "
                   "up to 6 dofs per node.";
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

          auto prb_name = fe_ptr->problemPtr->getName();
          auto get_flag = [&](int idx) { return (&mpc_bc->data.flag1)[idx]; };

          if (mpc_bc->mpcMasterEnts.size() != mpc_bc->mpcSlaveEnts.size()) {
            MOFEM_LOG("WORLD", Sev::error)
                << "Size of master nodes != Size of slave nodes : "
                << mpc_bc->mpcMasterEnts.size() << " != " << mpc_bc->mpcSlaveEnts.size();
            // SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
            //         "data inconsistency");
          }

          auto add_is = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          SmartPetscObj<IS> is_xyz_row_sum;
          SmartPetscObj<IS> is_m_row[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_m_col[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s_row[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s_col[max_nb_dofs_per_node];

          for (int dd = 0; dd != nb_field_coeffs; dd++) {
            if (get_flag(dd)) {
              // SmartPetscObj<IS> is_xyz_m;
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, ROW, field_name, dd, dd, is_m_row[dd],
                  &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, COL, field_name, dd, dd, is_m_col[dd],
                  &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, COL, field_name, dd, dd, is_s_col[dd],
                  &slave_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, ROW, field_name, dd, dd, is_s_row[dd],
                  &slave_verts);
              // ISView(is_s_row[dd], PETSC_VIEWER_STDOUT_WORLD);
              // ISView(is_s_col[dd], PETSC_VIEWER_STDOUT_WORLD);

              if (!mpc_bc->isReprocitical) {
                if (!is_xyz_row_sum) {
                  is_xyz_row_sum = is_s_row[dd];
                } else {
                  is_xyz_row_sum = add_is(is_xyz_row_sum, is_s_row[dd]);
                }
              }
            }
          }

          // if (is_xyz_row_sum) {
            SmartPetscObj<Mat> B =
                vLhs ? vLhs : SmartPetscObj<Mat>(fe_ptr->B, true);
            // The user is responsible for assembly if vLhs is provided
            MatType typem;
            CHKERR MatGetType(B, &typem);
            CHKERR MatSetOption(B, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
            if (typem == MATMPIBAIJ)
              CHKERR MatSetOption(B, MAT_USE_HASH_TABLE, PETSC_TRUE);

            if ((*fe_ptr->matAssembleSwitch) && !vLhs) {
              if (*fe_ptr->matAssembleSwitch) {
                CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
                *fe_ptr->matAssembleSwitch = false;
              }
            }

            if (vAO) {
              MOFEM_LOG("WORLD", Sev::error) << "No support for AO yet";
              //   MOFEM_LOG("WORLD", Sev::noisy) << "Apply AO to IS";
              //   CHKERR AOApplicationToPetscIS(vAO, is_xyz_m);
              //   CHKERR AOApplicationToPetscIS(vAO, is_xyz_s_c);
              //   CHKERR AOApplicationToPetscIS(vAO, is_xyz_s);
            }

            CHKERR MatZeroRowsIS(B, is_xyz_row_sum, 0, PETSC_NULL,
                                 PETSC_NULL);
          // }
            auto set_mat_values = [&](auto row_is, auto col_is, double val, double perturb = 0) {
              MoFEMFunctionBeginHot;
              const int *row_index_ptr;
              CHKERR ISGetIndices(row_is, &row_index_ptr);
              const int *col_index_ptr;
              CHKERR ISGetIndices(col_is, &col_index_ptr);
              int size, size2;
              CHKERR ISGetLocalSize(row_is, &size);
              CHKERR ISGetLocalSize(col_is, &size2);

              MatrixDouble locMat(size, size2);
              fill(locMat.data().begin(), locMat.data().end(), 0.0);
              for (int i = 0; i != size; i++)
                for (int j = 0; j != size2; j++)
                  if (i == j)
                    locMat(i, j) = val;

              if (mpc_bc->isReprocitical) {
                locMat *= perturb;
                CHKERR ::MatSetValues(B, size, row_index_ptr, size2,
                                      col_index_ptr, &*locMat.data().begin(),
                                      ADD_VALUES);
              } else {
                CHKERR ::MatSetValues(B, size, row_index_ptr, size2,
                                      col_index_ptr, &*locMat.data().begin(),
                                      INSERT_VALUES);
              }

              CHKERR ISRestoreIndices(row_is, &row_index_ptr);
              CHKERR ISRestoreIndices(col_is, &col_index_ptr);

              MoFEMFunctionReturnHot(0);
            };

            for (int dd = 0; dd != nb_field_coeffs; dd++) {
              if (get_flag(dd)) {
                if (!mpc_bc->isReprocitical) {
                  CHKERR set_mat_values(is_s_row[dd], is_s_col[dd], vDiag);
                  CHKERR set_mat_values(is_s_row[dd], is_m_col[dd], -vDiag);
                } 
                else {
                  auto &pn = mpc_bc->mPenalty;
                  CHKERR set_mat_values(is_s_row[dd], is_s_col[dd], vDiag, pn);
                  CHKERR set_mat_values(is_s_row[dd], is_m_col[dd], -vDiag, pn);
                  CHKERR set_mat_values(is_m_row[dd], is_m_col[dd], vDiag, pn);
                  CHKERR set_mat_values(is_m_row[dd], is_s_col[dd], -vDiag, pn);
                }
              }
            }

        } // if regex
      } // mpc loop
    } // bc loop
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  } // if fe_ptr
  MoFEMFunctionReturn(0);
}

EssentialPostProcRhs<MPCsType>::EssentialPostProcRhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Vec> rhs)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vRhs(rhs) {}

MoFEMErrorCode EssentialPostProcRhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    auto vec_mng = mField.getInterface<VecManager>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();
    // get connectivity from edge and assign master, slave ents
    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

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
          constexpr auto max_nb_dofs_per_node = 6;

          if (nb_field_coeffs > max_nb_dofs_per_node)
            MOFEM_LOG("WORLD", Sev::error)
                << "MultiPointConstraints PreProcLhs<MPCsType>: support only "
                   "up to 6 dofs per node for now.";

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

          auto master_verts = mpc_bc->mpcMasterEnts.subset_by_type(MBVERTEX);
          auto slave_verts = mpc_bc->mpcSlaveEnts.subset_by_type(MBVERTEX);

          auto prb_name = fe_method_ptr->problemPtr->getName();

          auto get_flag = [&](int idx) { return (&mpc_bc->data.flag1)[idx]; };

          auto add_is = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          SmartPetscObj<IS> is_m[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_sum;

          for (int dd = 0; dd != nb_field_coeffs; dd++) {
            if (get_flag(dd)) {
        
              CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                  prb_name, ROW, field_name, dd, dd, is_m[dd], &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                  prb_name, ROW, field_name, dd, dd, is_s[dd], &slave_verts);
              // slave
              if (!is_sum) {
                is_sum = is_s[dd];
              } else
                is_sum = add_is(is_sum, is_s[dd]);
  
            }
          }

          if (is_sum) {

            if (auto fe_ptr = fePtr.lock()) {
              auto snes_ctx = fe_ptr->snes_ctx;
              auto ts_ctx = fe_ptr->ts_ctx;
              bool is_nonlinear = snes_ctx != FEMethod::CTX_SNESNONE ||
                                        ts_ctx != FEMethod::CTX_TSNONE;
              is_nonlinear = is_nonlinear || mpc_bc->isReprocitical;
              SmartPetscObj<Vec> F =
                  vRhs ? vRhs : SmartPetscObj<Vec>(fe_ptr->f, true);

              if (fe_ptr->vecAssembleSwitch) {
                if ((*fe_ptr->vecAssembleSwitch) && !vRhs) {
                  CHKERR VecGhostUpdateBegin(F, ADD_VALUES, SCATTER_REVERSE);
                  CHKERR VecGhostUpdateEnd(F, ADD_VALUES, SCATTER_REVERSE);
                  CHKERR VecAssemblyBegin(F);
                  CHKERR VecAssemblyEnd(F);
                  *fe_ptr->vecAssembleSwitch = false;
                }
              }

              auto set_vec_values = [&](auto is_xyz_m, auto is_xyz_s, double val) {
                MoFEMFunctionBeginHot;
                const int *m_index_ptr;
                CHKERR ISGetIndices(is_xyz_m, &m_index_ptr);
                const int *s_index_ptr;
                CHKERR ISGetIndices(is_xyz_s, &s_index_ptr);
                int size_m;
                CHKERR ISGetLocalSize(is_xyz_m, &size_m);
                int size_s;
                CHKERR ISGetLocalSize(is_xyz_s, &size_s);

                // ISView(is_xyz_m, PETSC_VIEWER_STDOUT_WORLD);
                double *f;
                CHKERR VecGetArray(F, &f);
                if (size_m != size_s)
                  MOFEM_LOG("WORLD", Sev::error)
                      << "Size of master IS != Size of slave IS : " << size_m
                      << " != " << size_s;
                
                  if (is_nonlinear) {
                    auto x = fe_ptr->x;
                    auto tmp_x = vectorDuplicate(F);

                    CHKERR vec_mng->setLocalGhostVector(problem_name, ROW,
                                                        tmp_x, INSERT_VALUES,
                                                        SCATTER_FORWARD);
                    const double *v;
                    const double *u;

                    CHKERR VecGetArrayRead(tmp_x, &u);
                    CHKERR VecGetArrayRead(x, &v);
                    const int contrb = mpc_bc->isReprocitical ? 1 : 0;
                    if (size_m && size_s) {
                      for (auto i = 0; i != size_s; ++i) {
                        f[s_index_ptr[i]] =
                            val * (v[s_index_ptr[i]] - v[m_index_ptr[i]]) +
                            f[s_index_ptr[i]] * contrb;
                      }
                    }
                    CHKERR VecRestoreArrayRead(x, &v);
                    CHKERR VecRestoreArrayRead(tmp_x, &u);
                  } else {
                    if (size_m && size_s) 
                    for (auto i = 0; i != size_s; ++i) {
                      f[s_index_ptr[i]] = 0;
                    }
                  }
                
                CHKERR VecRestoreArray(F, &f);
                CHKERR ISRestoreIndices(is_xyz_m, &m_index_ptr);
                CHKERR ISRestoreIndices(is_xyz_s, &s_index_ptr);

                MoFEMFunctionReturnHot(0);
              };

              for (int dd = 0; dd != nb_field_coeffs; dd++)
                if (get_flag(dd)) {

                  if (!mpc_bc->isReprocitical) {
                    CHKERR set_vec_values(is_m[dd], is_s[dd], vDiag);
                  } else {
                    auto &pn = mpc_bc->mPenalty;
                    CHKERR set_vec_values(is_m[dd], is_s[dd], pn);
                    CHKERR set_vec_values(is_s[dd], is_m[dd], pn);
                  }
                }
            };

            // User is responsible for assembly if vLhs is provided

            // ISView(is_xyz_m, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERR MatZeroRowsColumnsIS(B, is_xyz_m, vDiag, PETSC_NULL,
            //                             PETSC_NULL);
          } else {
            MOFEM_LOG("WORLD", Sev::error)
                << "Cannot create ISs for MPCs for field (Rhs): " << field_name;
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