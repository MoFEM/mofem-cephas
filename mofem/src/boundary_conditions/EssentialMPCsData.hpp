/**
 * @file EssentialMPCsData.hpp
 * @brief Specialisation for MPCsTypes
 * @version 0.13.2
 * @date 2022-09-18
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef _ESSENTIAL_MPCS_DATA_HPP_
#define _ESSENTIAL_MPCS_DATA_HPP_

namespace MoFEM {

enum class MPC { TIE, RIGID_BODY, COUPLING, EMBEDDED_REGION, EQUATION, LAST };

// Base struct for all multipoint constraints data.
struct MPCsType : public DisplacementCubitBcData {

  
  Range mpcMasterEnts;
  Range mpcSlaveEnts;
  std::vector<unsigned char> bcMasterMarkers;
  std::vector<unsigned char> bcSlaveMarkers;

  MPC mpcType;
  MPCsType() : mpcType(MPC::TIE) {}
  // virtual ~MPCsType() {}
  using InitConditionsMap = std::map<EntityHandle, VectorDouble>;
  InitConditionsMap conditionsMap; 
  
  enum class CouplingPairs {
    POINT_TO_POINT,
    POINT_TO_SURFACE,
    SURFACE_TO_SURFACE
  };
};

// struct Kinematic : public MPCsType {};
// struct Distributing : public MPCsType {};

/**
 * @brief Type generating multipoint constraints.
 *
 */
template <> struct EssentialPreProc<MPCsType> {
  EssentialPreProc(MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr,
              bool is_spatial_positions = false);

  MoFEMErrorCode operator()();

  MoFEMErrorCode setMPCParentAdjacency();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  bool isSpatialPositions;

  boost::shared_ptr<ParentFiniteElementAdjacencyFunction<3>>
      parentAdjFunctionDim3;
  boost::shared_ptr<ParentFiniteElementAdjacencyFunction<2>>
      parentAdjFunctionDim2;
  boost::shared_ptr<ParentFiniteElementAdjacencyFunction<1>>
      parentAdjFunctionDim1;
};

/**
 * @brief Specialization for DisplacementCubitBcData
 *
 * @tparam
 */
template <> struct EssentialPostProcRhs<MPCsType> {
  EssentialPostProcRhs(MoFEM::Interface &m_field,
                      boost::shared_ptr<FEMethod> fe_ptr,
                      double diag = 1, SmartPetscObj<Vec> rhs = nullptr);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  double vDiag;
  SmartPetscObj<Vec> vRhs;
};

/**
 * @brief Specialization for MPCsType
 *
 * @tparam
 */
template <> struct EssentialPostProcLhs<MPCsType> {
  EssentialPostProcLhs(MoFEM::Interface &m_field,
                      boost::shared_ptr<FEMethod> fe_ptr,
                      double diag = 1,
                      SmartPetscObj<Mat> lhs = nullptr,
                      SmartPetscObj<AO> ao = nullptr);

  MoFEMErrorCode operator()();

protected:
  MoFEM::Interface &mField;
  boost::weak_ptr<FEMethod> fePtr;
  double vDiag;
  SmartPetscObj<Mat> vLhs;
  SmartPetscObj<AO> vAO;
};


template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
struct OpEssentialRhsImpl<MPCsType, 1, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template LinearForm<
          I>::template OpSource<1, FIELD_DIM> {

  using OpSource = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template LinearForm<I>::template OpSource<1, FIELD_DIM>;

  OpEssentialRhsImpl(const std::string field_name,
                     boost::shared_ptr<MPCsType> bc_data, // FIXME:
                     boost::shared_ptr<Range> ents_ptr,
                     std::vector<boost::shared_ptr<ScalingMethod>> smv);

  VectorFun<FIELD_DIM>
  initConstraintFunction(boost::shared_ptr<MPCsType> bc_data) {
    // FIXME: here implement functions depending on bc_data (MPC data)
    return [this](double x, double y, double z) { return tVal; };
  }

private:
  FTensor::Tensor1<double, FIELD_DIM> tVal;
  VecOfTimeScalingMethods vecOfTimeScalingMethods;
};

template <int FIELD_DIM, AssemblyType A, IntegrationType I, typename OpBase>
OpEssentialRhsImpl<MPCsType, 1, FIELD_DIM, A, I, OpBase>::OpEssentialRhsImpl(
    const std::string field_name,
    boost::shared_ptr<MPCsType> bc_data, // FIXME:
    boost::shared_ptr<Range> ents_ptr,
    std::vector<boost::shared_ptr<ScalingMethod>> smv)
    : OpSource(field_name, initConstraintFunction(bc_data), ents_ptr),
      vecOfTimeScalingMethods(smv) {
  // static_assert(FIELD_DIM > 1, "Is not implemented for scalar field");
}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct OpEssentialLhsImpl<MPCsType, BASE_DIM, FIELD_DIM, A, I, OpBase>
    : FormsIntegrators<OpBase>::template Assembly<A>::template BiLinearForm<
          I>::template OpMass<BASE_DIM, FIELD_DIM> {

  using OpMass = typename FormsIntegrators<OpBase>::template Assembly<
      A>::template BiLinearForm<I>::template OpMass<BASE_DIM, FIELD_DIM>;

  OpEssentialLhsImpl(const std::string field_name,
                     boost::shared_ptr<Range> ents_ptr);
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
OpEssentialLhsImpl<MPCsType, BASE_DIM, FIELD_DIM, A, I,
                   OpBase>::OpEssentialLhsImpl(const std::string field_name,
                                               boost::shared_ptr<Range>
                                                   ents_ptr)
    : OpMass(
          field_name, field_name,
          [](double, double, double) constexpr { return 1; }, ents_ptr) {}

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToRhsPipelineImpl<

    OpEssentialRhsImpl<MPCsType, BASE_DIM, FIELD_DIM, A, I, OpBase>, A, I,
    OpBase

    > {

  AddEssentialToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name,
      boost::shared_ptr<MatrixDouble> field_mat_ptr,
      std::vector<boost::shared_ptr<ScalingMethod>> smv

  ) {
    MoFEMFunctionBegin;

    using OP = typename EssentialBC<OpBase>::template Assembly<A>::
        template LinearForm<I>::template OpEssentialRhs<
            MPCsType, BASE_DIM,
            FIELD_DIM>; // TODO: we might need to implement OpEssentialRhsMaster
                        // and OpEssentialRhsSlave
    using OpInternal = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpBaseTimesVector<BASE_DIM,
                                                               FIELD_DIM, 1>;

    auto add_op = [&](auto &bcs) {
      MoFEMFunctionBeginHot;
      for (auto &m : bcs) {
        if (auto bc = m.second->mpcPtr) {
          auto &bc_id = m.first;
          auto regex_str =
              (boost::format("%s_%s_(.*)") % problem_name %
              field_name).str();
          if (std::regex_match(bc_id, std::regex(regex_str))) {
            MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpMPCsRhs") << *bc;
              pipeline.push_back(
                  new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
              pipeline.push_back(
                  new OP(field_name, bc, m.second->getBcEntsPtr(), smv));
              pipeline.push_back(new OpInternal(
                  field_name, field_mat_ptr,
                  [](double, double, double) constexpr { return 1.; },
                  m.second->getBcEntsPtr()));
              pipeline.push_back(new OpUnSetBc(field_name));

            // else 
            // {
            //   // specific for master
            // }
            // switch (bc->mpcType) {
            // case MPCsType::MPC::TIE:
            // case MPCsType::MPC::RIGID_BODY:
            // case MPCsType::MPC::COUPLING:
            //   // in case for tie, rigid body and coupling for MASTER, we dont
            //   // have to do anything
            //   break;
            // }
            //
          }
        }
      }
      MOFEM_LOG_CHANNEL("SELF");
      MoFEMFunctionReturnHot(0);
    };

    // CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Not implemented for MPCs");
    MoFEMFunctionReturn(0);
  }
};

template <int BASE_DIM, int FIELD_DIM, AssemblyType A, IntegrationType I,
          typename OpBase>
struct AddEssentialToLhsPipelineImpl<

    OpEssentialLhsImpl<MPCsType, BASE_DIM, FIELD_DIM, A, I, OpBase>, A, I,
    OpBase

    > {

  AddEssentialToLhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      MoFEM::Interface &m_field,
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      const std::string problem_name, std::string field_name

  ) {
    MoFEMFunctionBegin;

    // using OP = typename PCsBC<OpBase>::template Assembly<A>::
    //     template BiLinearForm<I>::template OpPCsLhs<
    //         MPCsType, BASE_DIM, FIELD_DIM>;

    // auto add_op = [&](auto &bcs) {
    //   MoFEMFunctionBeginHot;
    //   for (auto &m : bcs) {
    //     if (auto bc = m.second->dispBcPtr) {
    //       auto &bc_id = m.first;
    //       auto regex_str =
    //           (boost::format("%s_%s_(.*)") % problem_name %
    //           field_name).str();
    //       if (std::regex_match(bc_id, std::regex(regex_str))) {
    //         MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "OpPCsLhs") << *bc;
    //         pipeline.push_back(
    //             new OpSetBc(field_name, false, m.second->getBcMarkersPtr()));
    //         pipeline.push_back(new OP(field_name, m.second->getBcEntsPtr()));
    //         pipeline.push_back(new OpUnSetBc(field_name));
    //       }
    //     }
    //   }
    //   MOFEM_LOG_CHANNEL("SELF");
    // MoFEMFunctionReturnHot(0);
    // };

    // CHKERR add_op(m_field.getInterface<BcManager>()->getBcMapByBlockName());

    MoFEMFunctionReturn(0);
  }
};

} // namespace MoFEM

#endif //_ESSENTIAL_MPCS_DATA_HPP_