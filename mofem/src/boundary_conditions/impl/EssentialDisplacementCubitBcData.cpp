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

          FTensor::Tensor1<double, 3> t_angles{0., 0., 0.};
          FTensor::Tensor1<double, 3> t_vals{0., 0., 0.};
          FTensor::Tensor1<double, 3> t_off{0., 0., 0.};

          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation const *>(
                      disp_bc.get())) {
            for (int a = 0; a != 3; ++a)
              t_off(a) = ext_disp_bc->rotOffset[a];
          }

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

          const bool is_rotation =
              disp_bc->data.flag4 || disp_bc->data.flag5 || disp_bc->data.flag6;

          auto lambda = [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
            MoFEMFunctionBegin;

            auto v = t_vals(coeff);
            if (is_rotation) {
              FTensor::Tensor1<double, 3> t_coords(
                  coords[0][idx], coords[1][idx], coords[2][idx]);
              v += DisplacementCubitBcDataWithRotation::GetRotDisp(
                  t_angles, t_coords, t_off)(coeff);
            }
            if (getCoords) {
              v += coords[coeff][idx];
            }

            field_entity_ptr->getEntFieldData()[coeff] = v;
            ++idx;

            MoFEMFunctionReturn(0);
          };

          auto zero_lambda =
              [&](boost::shared_ptr<FieldEntity> field_entity_ptr) {
                MoFEMFunctionBegin;
                auto size = field_entity_ptr->getEntFieldData().size();
                for (int i = coeff; i < size; i += nb_field_coeffs)
                  field_entity_ptr->getEntFieldData()[i] = 0;
                MoFEMFunctionReturn(0);
              };

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);
          auto not_verts = subtract(bc.second->bcEnts, verts);

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
            CHKERR fb->fieldLambdaOnEntities(zero_lambda, field_name,
                                             &not_verts);
          }
          if (disp_bc->data.flag2 || disp_bc->data.flag4 ||
              disp_bc->data.flag6) {
            idx = 0;
            coeff = 1;
            if (nb_field_coeffs > 1) {
              CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
              CHKERR fb->fieldLambdaOnEntities(zero_lambda, field_name,
                                               &not_verts);
            }
          }
          if (disp_bc->data.flag3 || disp_bc->data.flag4 ||
              disp_bc->data.flag5 || is_rotation) {
            idx = 0;
            coeff = 2;
            if (nb_field_coeffs > 2) {
              CHKERR fb->fieldLambdaOnEntities(lambda, field_name, &verts);
              CHKERR fb->fieldLambdaOnEntities(zero_lambda, field_name,
                                               &not_verts);
            }
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

EssentialPostProcRhs<DisplacementCubitBcData>::EssentialPostProcRhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Vec> rhs)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vRhs(rhs) {}

MoFEMErrorCode EssentialPostProcRhs<DisplacementCubitBcData>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto vec_mng = mField.getInterface<VecManager>();
    auto is_mng = mField.getInterface<ISManager>();

    const auto problem_name = fe_method_ptr->problemPtr->getName();

    SmartPetscObj<IS> is_sum;

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto disp_bc = bc.second->dispBcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply EssentialPreProc<DisplacementCubitBcData>: "
              << problem_name << "_" << field_name << "_" << block_name;

          const bool is_rotation =
              disp_bc->data.flag4 || disp_bc->data.flag5 || disp_bc->data.flag6;

          auto ents = bc.second->bcEnts;

          std::array<SmartPetscObj<IS>, 3> is_xyz;
          auto prb_name = fe_method_ptr->problemPtr->getName();

          if (disp_bc->data.flag1 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                prb_name, ROW, field_name, 0, 0, is_xyz[0], &ents);
          }
          if (disp_bc->data.flag2 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                prb_name, ROW, field_name, 1, 1, is_xyz[1], &ents);
          }
          if (disp_bc->data.flag3 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                prb_name, ROW, field_name, 2, 2, is_xyz[2], &ents);
          }

          auto get_is_sum = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          for (auto &is : is_xyz) {
            if (is) {
              if (!is_sum) {
                is_sum = is;
              } else {
                is_sum = get_is_sum(is_sum, is);
              }
            }
          }
        }
      }
    }

    if (is_sum) {
      if (auto fe_ptr = fePtr.lock()) {

        auto snes_ctx = fe_ptr->snes_ctx;
        auto ts_ctx = fe_ptr->ts_ctx;

        SmartPetscObj<Vec> f =
            vRhs ? vRhs : SmartPetscObj<Vec>(fe_ptr->f, true);

        if (fe_ptr->vecAssembleSwitch) {
          if ((*fe_ptr->vecAssembleSwitch) && !vRhs) {
            CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
            CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
            CHKERR VecAssemblyBegin(f);
            CHKERR VecAssemblyEnd(f);
            *fe_ptr->vecAssembleSwitch = false;
          }
        }

        const int *index_ptr;
        CHKERR ISGetIndices(is_sum, &index_ptr);
        int size;
        CHKERR ISGetLocalSize(is_sum, &size);
        double *a;
        CHKERR VecGetArray(f, &a);

        auto tmp_x = vectorDuplicate(f);
        CHKERR vec_mng->setLocalGhostVector(problem_name, ROW, tmp_x,
                                            INSERT_VALUES, SCATTER_FORWARD);
        const double *u;
        CHKERR VecGetArrayRead(tmp_x, &u);

        if (snes_ctx != FEMethod::CTX_SNESNONE ||
            ts_ctx != FEMethod::CTX_TSNONE) {

          auto x = fe_ptr->x;
          CHKERR VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
          CHKERR VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

          const double *v;
          CHKERR VecGetArrayRead(x, &v);

          for (auto i = 0; i != size; ++i) {
            a[index_ptr[i]] = vDiag * (v[index_ptr[i]] - u[index_ptr[i]]);
          }

          CHKERR VecRestoreArrayRead(x, &v);

        } else {
          for (auto i = 0; i != size; ++i) {
            a[index_ptr[i]] = vDiag * u[index_ptr[i]];
          }
        }

        CHKERR VecRestoreArrayRead(tmp_x, &u);
        CHKERR VecRestoreArray(f, &a);
        CHKERR ISRestoreIndices(is_sum, &index_ptr);
      }
    }

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}

EssentialPostProcLhs<DisplacementCubitBcData>::EssentialPostProcLhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Mat> lhs, SmartPetscObj<AO> ao)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vLhs(lhs), vAO(ao) {}

MoFEMErrorCode EssentialPostProcLhs<DisplacementCubitBcData>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();

    const auto problem_name = fe_method_ptr->problemPtr->getName();

    SmartPetscObj<IS> is_sum;

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto disp_bc = bc.second->dispBcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply EssentialPreProc<DisplacementCubitBcData>: "
              << problem_name << "_" << field_name << "_" << block_name;

          const bool is_rotation =
              disp_bc->data.flag4 || disp_bc->data.flag5 || disp_bc->data.flag6;

          auto ents = bc.second->bcEnts;

          std::array<SmartPetscObj<IS>, 3> is_xyz;
          auto prb_name = fe_method_ptr->problemPtr->getName();

          if (disp_bc->data.flag1 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRank(
                prb_name, ROW, field_name, 0, 0, is_xyz[0], &ents);
          }
          if (disp_bc->data.flag2 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRank(
                prb_name, ROW, field_name, 1, 1, is_xyz[1], &ents);
          }
          if (disp_bc->data.flag3 || is_rotation) {
            CHKERR is_mng->isCreateProblemFieldAndRank(
                prb_name, ROW, field_name, 2, 2, is_xyz[2], &ents);
          }

          auto get_is_sum = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          for (auto &is : is_xyz) {
            if (is) {
              if (!is_sum) {
                is_sum = is;
              } else {
                is_sum = get_is_sum(is_sum, is);
              }
            }
          }
        }
      }
    }

    if (is_sum) {
      if (auto fe_ptr = fePtr.lock()) {
        SmartPetscObj<Mat> B =
            vLhs ? vLhs : SmartPetscObj<Mat>(fe_ptr->B, true);
        // User is responsible for assembly if vLhs is provided
        if ((*fe_ptr->matAssembleSwitch) && !vLhs) {
          if (*fe_ptr->matAssembleSwitch) {
            CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
            CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
            *fe_ptr->matAssembleSwitch = false;
          }
        }
        if (vAO) {
          MOFEM_LOG("WORLD", Sev::noisy) << "Apply AO to IS";
          CHKERR AOApplicationToPetscIS(vAO, is_sum);
        }
        // ISView(is_sum, PETSC_VIEWER_STDOUT_WORLD);
        CHKERR MatZeroRowsColumnsIS(B, is_sum, vDiag, PETSC_NULL, PETSC_NULL);
      }
    }

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}

EssentialPreProcReaction<DisplacementCubitBcData>::EssentialPreProcReaction(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
    SmartPetscObj<Vec> rhs, LogManager::SeverityLevel sev)
    : mField(m_field), fePtr(fe_ptr), vRhs(rhs), sevLevel(sev),
      printBlockName(PETSC_FALSE) {

  CHK_THROW_MESSAGE(PetscOptionsGetBool(PETSC_NULL, PETSC_NULL,
                                        "-reaction_print_block_name",
                                        &printBlockName, PETSC_NULL),
                    "can not get option");
}

MoFEMErrorCode EssentialPreProcReaction<DisplacementCubitBcData>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  enum { X = 0, Y, Z, MX, MY, MZ, LAST };

  if (auto fe_ptr = fePtr.lock()) {

    SmartPetscObj<Vec> f = vRhs ? vRhs : SmartPetscObj<Vec>(fe_ptr->f, true);

    if (fe_ptr->vecAssembleSwitch) {
      if ((*fe_ptr->vecAssembleSwitch) && !vRhs) {
        CHKERR VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE);
        CHKERR VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE);
        CHKERR VecAssemblyBegin(f);
        CHKERR VecAssemblyEnd(f);
        *fe_ptr->vecAssembleSwitch = false;
      }
    }

    auto get_low_hi_uid_by_entities = [](auto bit_number, auto f, auto s) {
      return std::make_pair(DofEntity::getLoFieldEntityUId(bit_number, f),
                            DofEntity::getHiFieldEntityUId(bit_number, s));
    };

    auto get_low_hi = [fe_ptr](auto lo_uid, auto hi_uid) {
      auto it = fe_ptr->problemPtr->numeredRowDofsPtr->get<Unique_mi_tag>()
                    .lower_bound(lo_uid);
      auto hi_it = fe_ptr->problemPtr->numeredRowDofsPtr->get<Unique_mi_tag>()
                       .upper_bound(hi_uid);
      return std::make_pair(it, hi_it);
    };

    auto mpi_array_reduce = [this](auto &array) {
      std::array<double, LAST> array_sum{0, 0, 0, 0, 0, 0};
      MPI_Allreduce(&array[0], &array_sum[0], LAST, MPI_DOUBLE, MPI_SUM,
                    mField.get_comm());
      return array_sum;
    };

    const double *a;
    CHKERR VecGetArrayRead(f, &a);

    auto bc_mng = mField.getInterface<BcManager>();
    const auto problem_name = fe_ptr->problemPtr->getName();
    const auto nb_local_dofs = fe_ptr->problemPtr->nbLocDofsRow;

    std::array<double, LAST> total_reactions{0, 0, 0, 0, 0, 0};

    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto disp_bc = bc.second->dispBcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          MOFEM_TAG_AND_LOG("WORLD", sevLevel, "Essential")
              << "EssentialPreProc<DisplacementCubitBcData>: " << problem_name
              << "_" << field_name << "_" << block_name;
          auto bit_number = mField.get_field_bit_number(field_name);

          FTensor::Tensor1<double, 3> t_off{0., 0., 0.};
          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation const *>(
                      disp_bc.get())) {
            for (int a = 0; a != 3; ++a)
              t_off(a) = ext_disp_bc->rotOffset[a];
          }

          auto verts = bc.second->bcEnts.subset_by_type(MBVERTEX);

          FTensor::Index<'i', 3> i;
          FTensor::Index<'j', 3> j;
          FTensor::Index<'k', 3> k;

          auto get_coords_vec = [&]() {
            VectorDouble coords_vec(3 * verts.size());
            if (verts.size()) {
              CHKERR mField.get_moab().get_coords(verts, &*coords_vec.begin());
            }
            return coords_vec;
          };

          auto coords_vec = get_coords_vec();
          std::array<double, LAST> reactions{0, 0, 0, 0, 0, 0};

          for (auto pit = verts.const_pair_begin();
               pit != verts.const_pair_end(); ++pit) {
            auto [lo_uid, hi_uid] =
                get_low_hi_uid_by_entities(bit_number, pit->first, pit->second);
            auto [lo, hi] = get_low_hi(lo_uid, hi_uid);

            for (; lo != hi; ++lo) {
              const auto loc_dof = (*lo)->getPetscLocalDofIdx();
              if (loc_dof < nb_local_dofs) {

                const auto coeff = (*lo)->getDofCoeffIdx();

                if (

                    ((disp_bc->data.flag1 || disp_bc->data.flag4) &&
                     coeff == 0) ||
                    ((disp_bc->data.flag2 || disp_bc->data.flag5) &&
                     coeff == 1) ||
                    ((disp_bc->data.flag3 || disp_bc->data.flag6) &&
                     coeff == 2)) {

                  const auto ent = (*lo)->getEnt();
                  reactions[coeff] += a[loc_dof];

                  auto force = [&]() {
                    FTensor::Tensor1<double, 3> t_force{0., 0., 0.};
                    t_force(coeff) = a[loc_dof];
                    return t_force;
                  };

                  auto coord = [&]() {
                    const auto idx = verts.index(ent);
                    FTensor::Tensor1<double, 3> t_coords{
                        coords_vec[3 * idx + X], coords_vec[3 * idx + Y],
                        coords_vec[3 * idx + Z]};
                    t_coords(i) -= t_off(i);
                    return t_coords;
                  };

                  auto moment = [&](auto &&t_force, auto &&t_coords) {
                    FTensor::Tensor1<double, 3> t_moment;
                    t_moment(i) =
                        (FTensor::levi_civita<double>(i, j, k) * t_coords(k)) *
                        t_force(j);
                    return t_moment;
                  };

                  auto t_moment = moment(force(), coord());
                  reactions[MX] += t_moment(X);
                  reactions[MY] += t_moment(Y);
                  reactions[MZ] += t_moment(Z);
                }
              }
            }
          }

          FTensor::Tensor1<double, 3> t_force{reactions[X], reactions[Y],
                                              reactions[Z]};
          FTensor::Tensor1<double, 3> t_moment{reactions[MX], reactions[MY],
                                               reactions[MZ]};
          FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_total_force{
              &total_reactions[X], &total_reactions[Y], &total_reactions[Z]};
          FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_total_moment{
              &total_reactions[MX], &total_reactions[MY], &total_reactions[MZ]};
          t_total_force(i) += t_force(i);
          t_total_moment(i) +=
              t_moment(i) +
              (FTensor::levi_civita<double>(i, j, k) * t_off(k)) * t_force(j);

          auto mpi_reactions = mpi_array_reduce(reactions);
          if (printBlockName) {
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Block %s Offset: %4.3e %4.3e %4.3e",
                                block_name.c_str(), t_off(X), t_off(Y),
                                t_off(Z));
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Block %s Force: %4.3e %4.3e %4.3e",
                                block_name.c_str(), mpi_reactions[X],
                                mpi_reactions[Y], mpi_reactions[Z]);
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Block %s Moment: %4.3e %4.3e %4.3e",
                                block_name.c_str(), mpi_reactions[MX],
                                mpi_reactions[MY], mpi_reactions[MZ]);
          } else {
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Offset: %4.3e %4.3e %4.3e", t_off(X), t_off(Y),
                                t_off(Z));
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Force: %4.3e %4.3e %4.3e", mpi_reactions[X],
                                mpi_reactions[Y], mpi_reactions[Z]);
            MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                                "Moment: %4.3e %4.3e %4.3e", mpi_reactions[MX],
                                mpi_reactions[MY], mpi_reactions[MZ]);
          }
        }
      }
    }

    CHKERR VecRestoreArrayRead(f, &a);

    auto mpi_total_reactions = mpi_array_reduce(total_reactions);
    MOFEM_TAG_AND_LOG_C(
        "WORLD", sevLevel, "Essential", "Total force: %4.3e %4.3e %4.3e",
        mpi_total_reactions[X], mpi_total_reactions[Y], mpi_total_reactions[Z]);
    MOFEM_TAG_AND_LOG_C("WORLD", sevLevel, "Essential",
                        "Total moment: %4.3e %4.3e %4.3e",
                        mpi_total_reactions[MX], mpi_total_reactions[MY],
                        mpi_total_reactions[MZ]);

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM