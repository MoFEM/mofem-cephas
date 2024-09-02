

/** \file ThermoElasticOps.hpp
 * \example ThermoElasticOps.hpp
 */

namespace ThermoElasticOps {

struct BlockedParameters
    : public boost::enable_shared_from_this<BlockedParameters> {
  double coeffExpansion;
  double refTemp;
  double heatConductivity;
  double heatCapacity;

  inline auto getCoeffExpansionPtr() {
    return boost::shared_ptr<double>(shared_from_this(), &coeffExpansion);
  }

  inline auto getRefTempPtr() {
    return boost::shared_ptr<double>(shared_from_this(), &refTemp);
  }

  inline auto getHeatConductivityPtr() {
    return boost::shared_ptr<double>(shared_from_this(), &heatConductivity);
  }

  inline auto getHeatCapacityPtr() {
    return boost::shared_ptr<double>(shared_from_this(), &heatCapacity);
  }
};

MoFEMErrorCode addMatThermalBlockOps(
    MoFEM::Interface &m_field,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
    std::string block_name,
    boost::shared_ptr<BlockedParameters> blockedParamsPtr, Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatThermalBlocks : public DomainEleOp {
    OpMatThermalBlocks(boost::shared_ptr<double> expansion_ptr,
                       boost::shared_ptr<double> ref_temp_ptr,
                       boost::shared_ptr<double> conductivity_ptr,
                       boost::shared_ptr<double> capacity_ptr,
                       MoFEM::Interface &m_field, Sev sev,
                       std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE),
          expansionPtr(expansion_ptr), refTempPtr(ref_temp_ptr),
          conductivityPtr(conductivity_ptr), capacityPtr(capacity_ptr) {
      CHK_THROW_MESSAGE(extractThermalBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      for (auto &b : blockData) {

        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          *expansionPtr = b.expansion;
          *refTempPtr = b.ref_temp;
          *conductivityPtr = b.conductivity;
          *capacityPtr = b.capacity;
          MoFEMFunctionReturnHot(0);
        }
      }

      *expansionPtr = coeff_expansion;
      *refTempPtr = ref_temp;
      *conductivityPtr = heat_conductivity;
      *capacityPtr = heat_capacity;

      MoFEMFunctionReturn(0);
    }

  private:
    struct BlockData {
      double expansion;
      double ref_temp;
      double conductivity;
      double capacity;
      Range blockEnts;
    };

    std::vector<BlockData> blockData;

    MoFEMErrorCode
    extractThermalBlockData(MoFEM::Interface &m_field,
                            std::vector<const CubitMeshSets *> meshset_vec_ptr,
                            Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermal Block") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() < 3) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Expected that block has two attributes");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR
          m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
          return ents;
        };

        MOFEM_TAG_AND_LOG("WORLD", sev, "Mat Thermal Block")
            << m->getName() << ": expansion = " << block_data[0]
            << " ref_temp = " << block_data[1]
            << " conductivity = " << block_data[2] << " capacity "
            << block_data[3];

        blockData.push_back({block_data[0], block_data[1], block_data[2],
                             block_data[3], get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    boost::shared_ptr<double> expansionPtr;
    boost::shared_ptr<double> refTempPtr;
    boost::shared_ptr<double> conductivityPtr;
    boost::shared_ptr<double> capacityPtr;
  };

  pipeline.push_back(new OpMatThermalBlocks(
      blockedParamsPtr->getCoeffExpansionPtr(),
      blockedParamsPtr->getRefTempPtr(),
      blockedParamsPtr->getHeatConductivityPtr(),
      blockedParamsPtr->getHeatCapacityPtr(), m_field, sev,

      // Get blockset using regular expression
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

struct OpKCauchyThermoElasticity : public AssemblyDomainEleOp {
  OpKCauchyThermoElasticity(const std::string row_field_name,
                            const std::string col_field_name,
                            boost::shared_ptr<MatrixDouble> mDptr,
                            boost::shared_ptr<double> coeff_expansion_ptr);

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data);

private:
  boost::shared_ptr<MatrixDouble> mDPtr;
  boost::shared_ptr<double> coeffExpansionPtr;
};

struct OpStressThermal : public DomainEleOp {
  /**
   * @deprecated do not use this constructor
   */
  DEPRECATED OpStressThermal(const std::string field_name,
                             boost::shared_ptr<MatrixDouble> strain_ptr,
                             boost::shared_ptr<VectorDouble> temp_ptr,
                             boost::shared_ptr<MatrixDouble> m_D_ptr,
                             boost::shared_ptr<double> coeff_expansion_ptr,
                             boost::shared_ptr<MatrixDouble> stress_ptr);

  OpStressThermal(boost::shared_ptr<MatrixDouble> strain_ptr,
                  boost::shared_ptr<VectorDouble> temp_ptr,
                  boost::shared_ptr<MatrixDouble> m_D_ptr,
                  boost::shared_ptr<double> coeff_expansion_ptr,
                  boost::shared_ptr<MatrixDouble> stress_ptr);

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

private:
  boost::shared_ptr<MatrixDouble> strainPtr;
  boost::shared_ptr<VectorDouble> tempPtr;
  boost::shared_ptr<MatrixDouble> mDPtr;
  boost::shared_ptr<double> coeffExpansionPtr;
  boost::shared_ptr<MatrixDouble> stressPtr;
};

OpKCauchyThermoElasticity::OpKCauchyThermoElasticity(
    const std::string row_field_name, const std::string col_field_name,
    boost::shared_ptr<MatrixDouble> mDptr,
    boost::shared_ptr<double> coeff_expansion_ptr)
    : AssemblyDomainEleOp(row_field_name, col_field_name,
                          DomainEleOp::OPROWCOL),
      mDPtr(mDptr), coeffExpansionPtr(coeff_expansion_ptr) {
  sYmm = false;
}

MoFEMErrorCode
OpKCauchyThermoElasticity::iNtegrate(EntitiesFieldData::EntData &row_data,
                                     EntitiesFieldData::EntData &col_data) {
  MoFEMFunctionBegin;

  auto &locMat = AssemblyDomainEleOp::locMat;

  const auto nb_integration_pts = row_data.getN().size1();
  const auto nb_row_base_functions = row_data.getN().size2();
  auto t_w = getFTensor0IntegrationWeight();

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
  auto t_row_diff_base = row_data.getFTensor1DiffN<SPACE_DIM>();
  auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM, 0>(*mDPtr);

  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'k', SPACE_DIM> k;
  FTensor::Index<'l', SPACE_DIM> l;
  FTensor::Tensor2_symmetric<double, SPACE_DIM> t_eigen_strain;
  t_eigen_strain(i, j) = (t_D(i, j, k, l) * t_kd(k, l)) * (*coeffExpansionPtr);

  for (auto gg = 0; gg != nb_integration_pts; ++gg) {

    double alpha = getMeasure() * t_w;
    auto rr = 0;
    for (; rr != AssemblyDomainEleOp::nbRows / SPACE_DIM; ++rr) {
      auto t_mat = getFTensor1FromMat<SPACE_DIM, 1>(locMat, rr * SPACE_DIM);
      auto t_col_base = col_data.getFTensor0N(gg, 0);
      for (auto cc = 0; cc != AssemblyDomainEleOp::nbCols; cc++) {

        t_mat(i) -=
            (t_row_diff_base(j) * t_eigen_strain(i, j)) * (t_col_base * alpha);

        ++t_mat;
        ++t_col_base;
      }

      ++t_row_diff_base;
    }
    for (; rr != nb_row_base_functions; ++rr)
      ++t_row_diff_base;

    ++t_w;
  }

  MoFEMFunctionReturn(0);
}

OpStressThermal::OpStressThermal(const std::string field_name,
                                 boost::shared_ptr<MatrixDouble> strain_ptr,
                                 boost::shared_ptr<VectorDouble> temp_ptr,
                                 boost::shared_ptr<MatrixDouble> m_D_ptr,
                                 boost::shared_ptr<double> coeff_expansion_ptr,
                                 boost::shared_ptr<MatrixDouble> stress_ptr)
    : DomainEleOp(field_name, DomainEleOp::OPROW), strainPtr(strain_ptr),
      tempPtr(temp_ptr), mDPtr(m_D_ptr), coeffExpansionPtr(coeff_expansion_ptr),
      stressPtr(stress_ptr) {
  // Operator is only executed for vertices
  std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
}

OpStressThermal::OpStressThermal(boost::shared_ptr<MatrixDouble> strain_ptr,
                                 boost::shared_ptr<VectorDouble> temp_ptr,
                                 boost::shared_ptr<MatrixDouble> m_D_ptr,
                                 boost::shared_ptr<double> coeff_expansion_ptr,
                                 boost::shared_ptr<MatrixDouble> stress_ptr)
    : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE), strainPtr(strain_ptr),
      tempPtr(temp_ptr), mDPtr(m_D_ptr), coeffExpansionPtr(coeff_expansion_ptr),
      stressPtr(stress_ptr) {}

//! [Calculate stress]
MoFEMErrorCode OpStressThermal::doWork(int side, EntityType type,
                                       EntData &data) {
  MoFEMFunctionBegin;
  const auto nb_gauss_pts = getGaussPts().size2();
  stressPtr->resize((SPACE_DIM * (SPACE_DIM + 1)) / 2, nb_gauss_pts);

  constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
  auto t_D = getFTensor4DdgFromMat<SPACE_DIM, SPACE_DIM, 0>(*mDPtr);
  auto t_strain = getFTensor2SymmetricFromMat<SPACE_DIM>(*strainPtr);
  auto t_stress = getFTensor2SymmetricFromMat<SPACE_DIM>(*stressPtr);
  auto t_temp = getFTensor0FromVec(*tempPtr);
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'k', SPACE_DIM> k;
  FTensor::Index<'l', SPACE_DIM> l;
  for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
    t_stress(i, j) =
        t_D(i, j, k, l) * (t_strain(k, l) - t_kd(k, l) * (t_temp - ref_temp) *
                                                (*coeffExpansionPtr));
    ++t_strain;
    ++t_stress;
    ++t_temp;
  }
  MoFEMFunctionReturn(0);
}
//! [Calculate stress]

struct SetTargetTemperature;

} // namespace ThermoElasticOps

//! [Target temperature]

template <AssemblyType A, IntegrationType I, typename OpBase>
struct MoFEM::OpFluxRhsImpl<ThermoElasticOps::SetTargetTemperature, 1, 1, A, I,
                            OpBase>;

template <AssemblyType A, IntegrationType I, typename OpBase>
struct MoFEM::OpFluxLhsImpl<ThermoElasticOps::SetTargetTemperature, 1, 1, A, I,
                            OpBase>;

template <AssemblyType A, IntegrationType I, typename OpBase>
struct MoFEM::AddFluxToRhsPipelineImpl<

    MoFEM::OpFluxRhsImpl<ThermoElasticOps::SetTargetTemperature, 1, 1, A, I,
                         OpBase>,
    A, I, OpBase

    > {

  AddFluxToRhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      boost::shared_ptr<VectorDouble> temp_ptr, std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP_SOURCE = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpSource<1, 1>;
    using OP_TEMP = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template LinearForm<I>::template OpBaseTimesScalar<1>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        std::vector<double> block_data;
        m->getAttributes(block_data);
        if (block_data.size() != 2) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                  "Expected two parameters");
        }
        double target_temperature = block_data[0];
        double beta =
            block_data[1]; // Set temperature parameter [ W/K * (1/m^3)]
        auto ents_ptr = boost::make_shared<Range>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->meshset,
                                                         *(ents_ptr), true);

        MOFEM_TAG_AND_LOG("WORLD", sev, "SetTargetTemperature")
            << "Add " << *m << " target temperature " << target_temperature
            << " penalty " << beta;

        pipeline.push_back(new OP_SOURCE(
            field_name,
            [target_temperature, beta](double, double, double) {
              return target_temperature * beta;
            },
            ents_ptr));
        pipeline.push_back(new OP_TEMP(
            field_name, temp_ptr,
            [beta](double, double, double) { return -beta; }, ents_ptr));
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR add_op(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

            (boost::format("%s(.*)") % block_name).str()

                ))

    );

    MOFEM_LOG_CHANNEL("WORLD");

    MoFEMFunctionReturn(0);
  }
};

template <AssemblyType A, IntegrationType I, typename OpBase>
struct AddFluxToLhsPipelineImpl<

    OpFluxLhsImpl<ThermoElasticOps::SetTargetTemperature, 1, 1, A, I, OpBase>,
    A, I, OpBase

    > {

  AddFluxToLhsPipelineImpl() = delete;

  static MoFEMErrorCode add(

      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pipeline,
      MoFEM::Interface &m_field, const std::string field_name,
      boost::shared_ptr<VectorDouble> temp_ptr, std::string block_name, Sev sev

  ) {
    MoFEMFunctionBegin;

    using OP_MASS = typename FormsIntegrators<OpBase>::template Assembly<
        A>::template BiLinearForm<I>::template OpMass<1, 1>;

    auto add_op = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        std::vector<double> block_data;
        m->getAttributes(block_data);
        if (block_data.size() != 2) {
          SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
                  "Expected two parameters");
        }
        double beta =
            block_data[1]; // Set temperature parameter [ W/K * (1/m^3)]
        auto ents_ptr = boost::make_shared<Range>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->meshset,
                                                         *(ents_ptr), true);

        MOFEM_TAG_AND_LOG("WORLD", sev, "SetTargetTemperature")
            << "Add " << *m << " penalty " << beta;

        pipeline.push_back(new OP_MASS(
            field_name, field_name,
            [beta](double, double, double) { return -beta; }, ents_ptr));
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR add_op(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

            (boost::format("%s(.*)") % block_name).str()

                ))

    );

    MOFEM_LOG_CHANNEL("WORLD");

    MoFEMFunctionReturn(0);
  }
};

/**
 * @brief Rhs for setting initial conditions
 *
 */
template <AssemblyType A, IntegrationType I>
struct OpRhsSetInitT : public AssemblyDomainEleOp {

  OpRhsSetInitT(const std::string field_name,
                boost::shared_ptr<VectorDouble> dot_T_ptr,
                boost::shared_ptr<VectorDouble> T_ptr,
                boost::shared_ptr<MatrixDouble> grad_T_ptr,
                double *initial_T_ptr, double *peak_T_ptr)
      : AssemblyDomainEleOp(field_name, field_name, AssemblyDomainEleOp::OPROW),
        dotTPtr(dot_T_ptr), TPtr(T_ptr), gradTPtr(grad_T_ptr),
        initialTPtr(initial_T_ptr), peakTPtr(peak_T_ptr) {}

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;

    const double vol = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    auto t_base = data.getFTensor0N();
    auto t_diff_base = data.getFTensor1DiffN<SPACE_DIM>();

#ifndef NDEBUG
    if (data.getDiffN().size1() != data.getN().size1())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 1");
    if (data.getDiffN().size2() != data.getN().size2() * SPACE_DIM) {
      MOFEM_LOG("SELF", Sev::error)
          << "Side " << rowSide << " " << CN::EntityTypeName(rowType);
      MOFEM_LOG("SELF", Sev::error) << data.getN();
      MOFEM_LOG("SELF", Sev::error) << data.getDiffN();
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 2");
    }
#endif

    auto t_T = getFTensor0FromVec(*TPtr);
    // auto t_grad_temp = getFTensor1FromMat<SPACE_DIM>(*gradTPtr);

    for (int gg = 0; gg != nbIntegrationPts; ++gg) {

      const double alpha = t_w * vol;

      const double set_T = init_T(*initialTPtr, *peakTPtr, t_coords(0),
                                  t_coords(1), t_coords(2));
      // const double m = get_M(set_T) * alpha;

      int bb = 0;
      for (; bb != nbRows; ++bb) {
        locF[bb] += (t_base * alpha) * (t_T - set_T);
        // locF[bb] += (t_diff_base(i) * m) * t_grad_g(i);
        ++t_base;
        ++t_diff_base;
      }

      for (; bb < nbRowBaseFunctions; ++bb) {
        ++t_base;
        ++t_diff_base;
      }

      ++t_T;
      // ++t_grad_g;

      ++t_coords;
      ++t_w;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> dotTPtr;
  boost::shared_ptr<VectorDouble> TPtr;
  boost::shared_ptr<MatrixDouble> gradTPtr;
  boost::shared_ptr<MatrixDouble> gradQPtr;
  double *initialTPtr;
  double *peakTPtr;
};

/**
 * @brief Lhs for setting initial conditions wrt T
 *
 */
template <AssemblyType A, IntegrationType I>
struct OpLhsSetInitT_dT : public AssemblyDomainEleOp {

  OpLhsSetInitT_dT(const std::string field_name,
                   boost::shared_ptr<VectorDouble> T_ptr)
      : AssemblyDomainEleOp(field_name, field_name,
                            AssemblyDomainEleOp::OPROWCOL),
        TPtr(T_ptr) {
    sYmm = false;
  }

  MoFEMErrorCode iNtegrate(EntitiesFieldData::EntData &row_data,
                           EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;

    const double vol = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    auto t_row_base = row_data.getFTensor0N();
    auto t_row_diff_base = row_data.getFTensor1DiffN<SPACE_DIM>();

#ifndef NDEBUG
    if (row_data.getDiffN().size1() != row_data.getN().size1())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 1");
    if (row_data.getDiffN().size2() != row_data.getN().size2() * SPACE_DIM) {
      MOFEM_LOG("SELF", Sev::error)
          << "Side " << rowSide << " " << CN::EntityTypeName(rowType);
      MOFEM_LOG("SELF", Sev::error) << row_data.getN();
      MOFEM_LOG("SELF", Sev::error) << row_data.getDiffN();
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 2");
    }

    if (col_data.getDiffN().size1() != col_data.getN().size1())
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 1");
    if (col_data.getDiffN().size2() != col_data.getN().size2() * SPACE_DIM) {
      MOFEM_LOG("SELF", Sev::error)
          << "Side " << rowSide << " " << CN::EntityTypeName(rowType);
      MOFEM_LOG("SELF", Sev::error) << col_data.getN();
      MOFEM_LOG("SELF", Sev::error) << col_data.getDiffN();
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong size 2");
    }
#endif

    auto t_T = getFTensor0FromVec(*TPtr);
    // auto t_grad_g = getFTensor1FromMat<SPACE_DIM>(*gradGPtr);

    for (int gg = 0; gg != nbIntegrationPts; gg++) {

      const double alpha = t_w * vol;

      int rr = 0;
      for (; rr != nbRows; ++rr) {

        auto t_col_base = col_data.getFTensor0N(gg, 0);
        auto t_col_diff_base = col_data.getFTensor1DiffN<SPACE_DIM>(gg, 0);

        for (int cc = 0; cc != nbCols; ++cc) {

          locMat(rr, cc) += (t_row_base * t_col_base * alpha);

          ++t_col_base;
          ++t_col_diff_base;
        }

        ++t_row_base;
        ++t_row_diff_base;
      }

      for (; rr < nbRowBaseFunctions; ++rr) {
        ++t_row_base;
        ++t_row_diff_base;
      }

      ++t_T;
      // ++t_grad_g;
      ++t_w;
      ++t_coords;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<VectorDouble> TPtr;
  // boost::shared_ptr<MatrixDouble> gradGPtr;
};

//! [Target temperature]
