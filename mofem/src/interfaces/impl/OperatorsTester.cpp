/**
 * @file OperatorsTester.cpp
 * @brief
 * @version 0.13.2
 * @date 2022-12-04
 *
 * @copyright Copyright (c) 2022
 *
 */

namespace MoFEM {

OperatorsTester::OperatorsTester(const Core &core)
    : cOre(const_cast<Core &>(core)) {}

MoFEMErrorCode
OperatorsTester::query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const {
  *iface = const_cast<OperatorsTester *>(this);
  return 0;
}

SmartPetscObj<Vec>
OperatorsTester::setRandomFields(SmartPetscObj<DM> dm,
                                 std::vector<RandomFieldData> random_fields,
                                 boost::shared_ptr<Range> ents_ptr) {

  MoFEM::Interface &m_field = cOre;

  auto prb_ptr = getProblemPtr(dm);

  auto get_random_number = [](auto &range) {
    const auto r = static_cast<double>(std::rand()) / RAND_MAX;
    return range[0] + r * (range[1] - range[0]);
  };

  auto v = smartCreateDMVector(dm);
  double *array;
  CHKERR VecGetArray(v, &array);

  MOFEM_LOG_CHANNEL("WORLD");

  for (auto &rf : random_fields) {
    const auto field_id = m_field.get_field_bit_number(rf.first);
    MOFEM_TAG_AND_LOG("WORLD", Sev::noisy, "OperatorsTester")
        << "Set random field " << rf.first << " " << static_cast<int>(field_id);

    if (!ents_ptr) {
      auto lo = prb_ptr->getNumeredRowDofsPtr()->lower_bound(
          FieldEntity::getLoBitNumberUId(field_id));
      const auto hi = prb_ptr->getNumeredRowDofsPtr()->upper_bound(
          FieldEntity::getHiBitNumberUId(field_id));
      for (; lo != hi; ++lo) {
        const auto idx = (*lo)->getPetscLocalDofIdx();
        array[idx] = get_random_number(rf.second);
      }
    } else {
      for (auto pit = ents_ptr->const_pair_begin();
           pit != ents_ptr->const_pair_end(); ++pit) {
        auto lo =
            prb_ptr->getNumeredRowDofsPtr()->get<Unique_mi_tag>().lower_bound(
                DofEntity::getLoFieldEntityUId(field_id, pit->first));
        auto hi =
            prb_ptr->getNumeredRowDofsPtr()->get<Unique_mi_tag>().upper_bound(
                DofEntity::getHiFieldEntityUId(field_id, pit->second));
        for (; lo != hi; ++lo) {
          const auto idx = (*lo)->getPetscLocalDofIdx();
          array[idx] = get_random_number(rf.second);
        }
      }
    }
  }

  CHKERR VecRestoreArray(v, &array);

  CHKERR VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
  CHKERR VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);

  return v;
};

SmartPetscObj<Vec>
OperatorsTester::assembleVec(SmartPetscObj<DM> dm, std::string fe_name,
                             boost::shared_ptr<FEMethod> pipeline,
                             SmartPetscObj<Vec> x, SmartPetscObj<Vec> delta_x,
                             SmartPetscObj<Vec> delta2_x, double time,
                             double delta_t, CacheTupleWeakPtr cache_ptr) {

  pipeline->data_ctx = PetscData::CtxSetF | PetscData::CtxSetTime;
  auto f = smartCreateDMVector(dm);
  pipeline->f = f;
  pipeline->ts_t = time;
  pipeline->ts_dt = delta_t;

  std::ignore = setPipelineX(pipeline, x, delta_x, delta2_x, delta_t);

  CHK_THROW_MESSAGE(
      DMoFEMMeshToLocalVector(dm, x, INSERT_VALUES, SCATTER_REVERSE),
      "loc to vec");
  CHK_THROW_MESSAGE(DMoFEMLoopFiniteElements(dm, fe_name, pipeline, cache_ptr),
                    "Run assemble pipeline");
  CHK_THROW_MESSAGE(VecAssemblyBegin(f), "Asmb vector");
  CHK_THROW_MESSAGE(VecAssemblyEnd(f), "Asmb vector");

  CHK_THROW_MESSAGE(VecGhostUpdateBegin(f, ADD_VALUES, SCATTER_REVERSE),
                    "Update vec");
  CHK_THROW_MESSAGE(VecGhostUpdateEnd(f, ADD_VALUES, SCATTER_REVERSE),
                    "Update vec");
  CHK_THROW_MESSAGE(VecGhostUpdateBegin(f, INSERT_VALUES, SCATTER_FORWARD),
                    "Update vec");
  CHK_THROW_MESSAGE(VecGhostUpdateEnd(f, INSERT_VALUES, SCATTER_FORWARD),
                    "Update vec");

  return f;
}

SmartPetscObj<Mat>
OperatorsTester::assembleMat(SmartPetscObj<DM> dm, std::string fe_name,
                             boost::shared_ptr<FEMethod> pipeline,
                             SmartPetscObj<Vec> x, SmartPetscObj<Vec> delta_x,
                             SmartPetscObj<Vec> delta2_x, double time,
                             double delta_t, CacheTupleWeakPtr cache_ptr) {
  auto m = smartCreateDMMatrix(dm);

  pipeline->data_ctx =
      PetscData::CtxSetA | PetscData::CtxSetB | PetscData::CtxSetTime;
  pipeline->A = m;
  pipeline->B = m;

  pipeline->ts_t = time;

  pipeline->ts_dt = delta_t;
  pipeline->ts_a = 1. / delta_t;
  pipeline->ts_aa = 1. / pow(delta_t, 2);

  std::ignore = setPipelineX(pipeline, x, delta_x, delta2_x, delta_t);

  CHK_THROW_MESSAGE(
      DMoFEMMeshToLocalVector(dm, x, INSERT_VALUES, SCATTER_REVERSE),
      "loc to vec");
  CHK_THROW_MESSAGE(DMoFEMLoopFiniteElements(dm, fe_name, pipeline, cache_ptr),
                    "Run assemble pipeline");
  CHK_THROW_MESSAGE(MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY), "Asmb mat");
  CHK_THROW_MESSAGE(MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY), "Asmb mat");

  return m;
}

SmartPetscObj<Vec> OperatorsTester::directionalCentralFiniteDifference(
    SmartPetscObj<DM> dm, std::string fe_name,
    boost::shared_ptr<FEMethod> pipeline, SmartPetscObj<Vec> x,
    SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
    SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
    CacheTupleWeakPtr cache_ptr) {

  auto axpy = [&](auto v0, auto diff_v, auto eps) {
    SmartPetscObj<Vec> v;
    if (v0.use_count()) {
      v = smartVectorDuplicate(v0);
      CHK_THROW_MESSAGE(VecCopy(v0, v), "Cpy");
      CHK_THROW_MESSAGE(VecAXPY(v, eps, diff_v), "Add");
      CHK_THROW_MESSAGE(VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD),
                        "ghost");
      CHK_THROW_MESSAGE(VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD),
                        "ghost");
    }
   return v;
  };

  auto fm1 = assembleVec(dm, fe_name, pipeline,
                         //
                         axpy(x, diff_x, -eps),
                         //
                         axpy(delta_x, diff_x, -eps),
                         //
                         axpy(delta2_x, diff_x, -eps),
                         //
                         time, delta_t, cache_ptr);

  auto fp1 = assembleVec(dm, fe_name, pipeline,
                         //
                         axpy(x, diff_x, eps),
                         //
                         axpy(delta_x, diff_x, eps),
                         //
                         axpy(delta2_x, diff_x, eps),
                         //
                         time, delta_t, cache_ptr);

  auto calculate_derivative = [eps](auto fm1, auto fp1) {
    int size;
    VecGetLocalSize(fp1, &size);
    VecAXPY(fp1, -1, fm1);
    double *array;
    VecGetArray(fp1, &array);
    for (int j = 0; j != size; ++j) {
      array[j] /= 2 * eps;
    }
    VecRestoreArray(fp1, &array);
    return fp1;
  };

  return calculate_derivative(fm1, fp1);
}

SmartPetscObj<Vec> OperatorsTester::checkCentralFiniteDifference(
    SmartPetscObj<DM> dm, std::string fe_name,
    boost::shared_ptr<FEMethod> pipeline_rhs,
    boost::shared_ptr<FEMethod> pipeline_lhs, SmartPetscObj<Vec> x,
    SmartPetscObj<Vec> delta_x, SmartPetscObj<Vec> delta2_x,
    SmartPetscObj<Vec> diff_x, double time, double delta_t, double eps,
    CacheTupleWeakPtr cache_ptr) {

  auto fd_diff = directionalCentralFiniteDifference(
      dm, fe_name, pipeline_rhs, x, delta_x, delta2_x, diff_x, time, delta_t,
      eps, cache_ptr);

  auto m = assembleMat(dm, fe_name, pipeline_lhs, x, delta_x, delta2_x, time,
                       delta_t, cache_ptr);

  auto fm = smartVectorDuplicate(fd_diff);
  CHK_THROW_MESSAGE(MatMult(m, diff_x, fm), "Mat mult");
  CHK_THROW_MESSAGE(VecAXPY(fd_diff, -1., fm), "Add");

  return fd_diff;
}

std::pair<SmartPetscObj<Vec>, SmartPetscObj<Vec>>
OperatorsTester::setPipelineX(boost::shared_ptr<FEMethod> pipeline,
                              SmartPetscObj<Vec> x, SmartPetscObj<Vec> delta_x,
                              SmartPetscObj<Vec> delta2_x, double delta_t) {

  // Set dofs vector to finite element instance

  pipeline->x = x;
  pipeline->data_ctx |= PetscData::CTX_SET_X;

  SmartPetscObj<Vec> x_t, x_tt;

  // Set velocity dofs vector to finite element instance

  if (delta_x.use_count()) {
    x_t = smartVectorDuplicate(x);
    VecCopy(delta_x, x_t);
    VecScale(x_t, 1. / delta_t);
    CHK_THROW_MESSAGE(VecGhostUpdateBegin(x_t, INSERT_VALUES, SCATTER_FORWARD),
                      "ghost");
    CHK_THROW_MESSAGE(VecGhostUpdateEnd(x_t, INSERT_VALUES, SCATTER_FORWARD),
                      "ghost");
    pipeline->x_t = x_t;
    pipeline->data_ctx |= PetscData::CTX_SET_X_T;
  } else {
    pipeline->x_t = PETSC_NULL;
  }

  // Set acceleration dofs vector to finite element instance

  if (delta2_x.use_count()) {
    x_tt = smartVectorDuplicate(x);
    VecCopy(delta2_x, x_tt);
    VecScale(x_tt, 1. / pow(delta_t, 2));
    CHK_THROW_MESSAGE(VecGhostUpdateBegin(x_tt, INSERT_VALUES, SCATTER_FORWARD),
                      "ghost");
    CHK_THROW_MESSAGE(VecGhostUpdateEnd(x_tt, INSERT_VALUES, SCATTER_FORWARD),
                      "ghost");
    pipeline->x_tt = x_tt;
    pipeline->data_ctx |= PetscData::CTX_SET_X_TT;
  } else {
    pipeline->x_tt = PETSC_NULL;
  }

  return std::make_pair(x_t, x_tt);
}
} // namespace MoFEM