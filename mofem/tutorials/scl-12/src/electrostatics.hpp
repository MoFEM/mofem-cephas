
#ifndef __ELECTROSTATICS_HPP__
#define __ELECTROSTATICS_HPP__

#include <MoFEM.hpp>
using namespace MoFEM;
constexpr auto domainField = "POTENTIAL";
constexpr int BASE_DIM = 1;
constexpr int FIELD_DIM = 1;
constexpr int SPACE_DIM = EXECUTABLE_DIMENSION;
const double bodySource = 0.0;

// aliases for elements types from the pipeline
using DomainEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::DomainEle;
using IntEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::BoundaryEle;
using DomainEleOp = DomainEle::UserDataOperator;
using IntEleOp = IntEle::UserDataOperator;
using EntData = EntitiesFieldData::EntData;
using PostProcEle = PostProcBrokenMeshInMoab<DomainEle>;
using SideEle = PipelineManager::ElementsAndOpsByDim<SPACE_DIM>::FaceSideEle;
using SideEleOp = SideEle::UserDataOperator;
using PostProcFaceEle =
    PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>;

template <int SPACE_DIM> struct intPostProc {};

template <> struct intPostProc<2> {
  using intEle = MoFEM::EdgeElementForcesAndSourcesCore;
};

template <> struct intPostProc<3> {
  using intEle = MoFEM::FaceElementForcesAndSourcesCore;
};

using intElementForcesAndSourcesCore = intPostProc<SPACE_DIM>::intEle;

// forms integrators to calculate the LHS and RHS
using OpDomainLhsMatrixK =
    FormsIntegrators<DomainEleOp>::Assembly<PETSC>::BiLinearForm<
        GAUSS>::OpGradGrad<BASE_DIM, FIELD_DIM, SPACE_DIM>; // lhs: K Matrix

using OpInterfaceRhsVectorF = FormsIntegrators<IntEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<BASE_DIM, FIELD_DIM>; // rhs: F Vector
using OpBodySourceVectorb = FormsIntegrators<DomainEleOp>::Assembly<
    PETSC>::LinearForm<GAUSS>::OpSource<BASE_DIM, FIELD_DIM>; // rhs: b Vector

static char help[] = "...\n\n";

// Common data structure for blocks
struct BlockData {
  int iD;
  double chargeDensity;
  double epsPermit;
  Range interfaceEnts;
  Range domainEnts;
  Range electrodeEnts;
  Range internalDomainEnts;
};

struct DataAtIntegrationPts {
  SmartPetscObj<Vec> petscVec;
  double blockPermittivity;
  double blockChrgDens;
  DataAtIntegrationPts(MoFEM::Interface &m_field) {
    blockPermittivity = 0.0;
    blockChrgDens = 0.0;
    PetscInt ghosts[2] = {0, 1};
    if (!m_field.get_comm_rank())
      petscVec = createGhostVector(m_field.get_comm(), 2, 2, 0, ghosts);
    else
      petscVec = createGhostVector(m_field.get_comm(), 0, 2, 2, ghosts);
  }
};

// Operator to get the charge density on the interface
struct OpBlockChargeDensity : public IntEleOp {
  OpBlockChargeDensity(
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr,
      boost::shared_ptr<std::map<int, BlockData>> int_block_sets_ptr,
      const std::string &field_name)
      : IntEleOp(field_name, field_name, OPROWCOL, false),
        commonDataPtr(common_data_ptr), intBlockSetsPtr(int_block_sets_ptr) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type, EntData &row_data,
                        EntData &col_data) {
    MoFEMFunctionBegin;
    for (const auto &m : *intBlockSetsPtr) {
      if (m.second.interfaceEnts.find(getFEEntityHandle()) !=
          m.second.interfaceEnts.end()) {
        commonDataPtr->blockChrgDens = m.second.chargeDensity;
      }
    }
    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
  boost::shared_ptr<std::map<int, BlockData>> intBlockSetsPtr;
};
// Operator to get the permittivity of the domain/block
struct OpBlockPermittivity : public DomainEleOp {

  OpBlockPermittivity(
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr,
      boost::shared_ptr<map<int, BlockData>> perm_block_sets_ptr,
      const std::string &field_name)
      : DomainEleOp(field_name, field_name, OPROWCOL, false),
        commonDataPtr(common_data_ptr), permBlockSetsPtr(perm_block_sets_ptr) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;
    for (auto &m : (*permBlockSetsPtr)) {
      if (m.second.domainEnts.find(getFEEntityHandle()) !=
          m.second.domainEnts.end()) {
        commonDataPtr->blockPermittivity = m.second.epsPermit;
      }
    }
    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
};

// Operator to calculate Total Electrgy density in Post Porocessing
struct OpTotalEnergyDensity : public DomainEleOp {
  OpTotalEnergyDensity(
      const std::string &field_name,
      boost::shared_ptr<MatrixDouble> grad_u_negative,
      boost::shared_ptr<std::map<int, BlockData>> perm_block_sets_ptr,
      boost::shared_ptr<std::map<int, BlockData>> int_domain_block_sets_ptr,
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr,
      SmartPetscObj<Vec> petscVec_Energy)
      : DomainEleOp(field_name, DomainEleOp::OPROW),
        gradUNegative(grad_u_negative), permBlockSetsPtr(perm_block_sets_ptr),
        intDomainBlocksetsPtr(int_domain_block_sets_ptr),
        commonDataPtr(common_data_ptr), PetsctotalEnergy(petscVec_Energy) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    auto t_negative_grad_u = getFTensor1FromMat<SPACE_DIM>(*gradUNegative);
    FTensor::Index<'I', SPACE_DIM> I;
    const double area = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    const auto nb_gauss_pts = getGaussPts().size2();
    double domainArea = 0.0;
    double totalEnergy = 0.0;
    int indexx = 0;
    // Total Energy Density = 0.5 * (E * E) * Perm * Area
    for (const auto &m : *intDomainBlocksetsPtr) {

      double blockPermittivity = 0.0;
      domainArea += area;

      if (m.second.internalDomainEnts.find(getFEEntityHandle()) !=
          m.second.internalDomainEnts.end()) {
        for (const auto &n : *permBlockSetsPtr) {
          if (n.second.domainEnts.find(getFEEntityHandle()) !=
              n.second.domainEnts.end()) {
            blockPermittivity = n.second.epsPermit;
          }
        }

        for (int gg = 0; gg != nb_gauss_pts; gg++) {
          totalEnergy += 0.5 *
                         std::abs(t_negative_grad_u(I) * t_negative_grad_u(I)) *
                         blockPermittivity * t_w * area;
          ++t_negative_grad_u;
          ++t_w;
        }
      }
    }
    CHKERR VecSetValues(PetsctotalEnergy, 1, &indexx, &totalEnergy, ADD_VALUES);

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> gradUNegative;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<std::map<int, BlockData>> intDomainBlocksetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
  SmartPetscObj<Vec> PetsctotalEnergy;
};

struct OpEnergyDensity : public DomainEleOp {

  OpEnergyDensity(
      const std::string &field_name,
      boost::shared_ptr<MatrixDouble> grad_u_negative,
      boost::shared_ptr<MatrixDouble> energy_densiity_ptr,
      boost::shared_ptr<std::map<int, BlockData>> perm_block_sets_ptr,
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr)
      : DomainEleOp(field_name, field_name, OPROWCOL, false),
        energyDensiity(energy_densiity_ptr), gradUNegative(grad_u_negative),
        permBlockSetsPtr(perm_block_sets_ptr), commonDataPtr(common_data_ptr) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;

    const size_t nb_gauss_pts = getGaussPts().size2();
    energyDensiity->resize(SPACE_DIM, nb_gauss_pts, false);
    energyDensiity->clear();
    auto t_energy_densiity = getFTensor1FromMat<SPACE_DIM>(*energyDensiity);
    auto t_negative_grad_u = getFTensor1FromMat<SPACE_DIM>(*gradUNegative);
    FTensor::Index<'I', SPACE_DIM> I;
    double blockPermittivity = 0.0;
    for (const auto &n : *permBlockSetsPtr) {
      if (n.second.domainEnts.find(getFEEntityHandle()) !=
          n.second.domainEnts.end()) {
        blockPermittivity = n.second.epsPermit;
      }
    }
    // E = 0.5 * (E * E) * Permittivity
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      t_energy_densiity(I) =
          0.5 * std::abs(t_negative_grad_u(I) * t_negative_grad_u(I)) *
          blockPermittivity;
      ++t_negative_grad_u;
      ++t_energy_densiity;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> gradUNegative;
  boost::shared_ptr<MatrixDouble> energyDensiity;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
};

struct OpGradTimesperm : public DomainEleOp {

  OpGradTimesperm(
      const std::string &field_name,
      boost::shared_ptr<MatrixDouble> grad_u_negative,
      boost::shared_ptr<MatrixDouble> grad_times_perm,
      boost::shared_ptr<std::map<int, BlockData>> perm_block_sets_ptr,
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr)
      : DomainEleOp(field_name, field_name, OPROWCOL, false),
        gradtimesperm(grad_times_perm), gradUNegative(grad_u_negative),
        permBlockSetsPtr(perm_block_sets_ptr), commonDataPtr(common_data_ptr) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }

  MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                        EntityType col_type,
                        EntitiesFieldData::EntData &row_data,
                        EntitiesFieldData::EntData &col_data) {
    MoFEMFunctionBegin;

    const size_t nb_gauss_pts = getGaussPts().size2();
    gradtimesperm->resize(SPACE_DIM, nb_gauss_pts, false);
    gradtimesperm->clear();
    auto t_grad_times_perm = getFTensor1FromMat<SPACE_DIM>(*gradtimesperm);
    auto t_negative_grad_u = getFTensor1FromMat<SPACE_DIM>(*gradUNegative);
    FTensor::Index<'I', SPACE_DIM> I;
    double blockPermittivity = 0.0;
    for (const auto &n : *permBlockSetsPtr) {
      if (n.second.domainEnts.find(getFEEntityHandle()) !=
          n.second.domainEnts.end()) {
        blockPermittivity = n.second.epsPermit;
      }
    }
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      t_grad_times_perm(I) = t_negative_grad_u(I) * blockPermittivity;

      ++t_negative_grad_u;
      ++t_grad_times_perm;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> gradUNegative;
  boost::shared_ptr<MatrixDouble> gradtimesperm;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
};

// operator to get the electric jump on the all side elecment of the domain: 
// d = E1 * perm1 - E2 * perm1
template <int SPACE_DIM> struct OpELectricJump : public SideEleOp {
  OpELectricJump(
      const std::string field_name, boost::shared_ptr<MatrixDouble> grad_ptr,
      boost::shared_ptr<MatrixDouble> d_jump,
      boost::shared_ptr<DataAtIntegrationPts> common_data_ptr,
      boost::shared_ptr<std::map<int, BlockData>> perm_block_sets_ptr)
      : SideEleOp(field_name, SideEleOp::OPROW, false), gradPtr(grad_ptr),
        djump(d_jump), commonDataPtr(common_data_ptr),
        permBlockSetsPtr(perm_block_sets_ptr)

  {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }

  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;

    FTensor::Index<'i', SPACE_DIM> i;
    auto t_field_grad = getFTensor1FromMat<SPACE_DIM>(*gradPtr);

    double blockPermittivity = 0.0;

    for (const auto &n : *permBlockSetsPtr) {
      if (n.second.domainEnts.find(getFEEntityHandle()) !=
          n.second.domainEnts.end()) {
        blockPermittivity = n.second.epsPermit;
      }
    }
    auto N_InLoop = getNinTheLoop();
    auto sensee = getSkeletonSense();
    auto nb_gauss_pts = getGaussPts().size2();

    if (N_InLoop == 0) {
      djump->resize(SPACE_DIM, nb_gauss_pts, false);
      djump->clear();
    }
    auto t_jump = getFTensor1FromMat<SPACE_DIM>(*djump);

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      t_jump(i) -= t_field_grad(i) * blockPermittivity * sensee;
      ++t_jump;
      ++t_field_grad;
    }

    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<MatrixDouble> gradPtr;
  boost::shared_ptr<MatrixDouble> djump;
  boost::shared_ptr<std::map<int, BlockData>> permBlockSetsPtr;
  boost::shared_ptr<DataAtIntegrationPts> commonDataPtr;
};

// operator to get the alpha value on the inerface: Alpha = d * n

template <int SPACE_DIM> struct OpAlpha : public IntEleOp {
  OpAlpha(const std::string field_name, boost::shared_ptr<MatrixDouble> d_jump,
          SmartPetscObj<Vec> alpha_vec,
          boost::shared_ptr<std::map<int, BlockData>> electrode_block_sets_ptr)
      : IntEleOp(field_name, IntEleOp::OPROW, false), djumpptr(d_jump),
        petscVec(alpha_vec), elecBlockSetsPtr(electrode_block_sets_ptr) {
    std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
    doEntities[MBVERTEX] = true;
  }
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
    MoFEMFunctionBegin;
    FTensor::Index<'i', SPACE_DIM> i;
    int index = 0; // electrode index no.
    const auto fe_ent = getFEEntityHandle();
    auto t_jump = getFTensor1FromMat<SPACE_DIM>(*djumpptr);
    auto t_normal = getFTensor1NormalsAtGaussPts();
    const auto nb_gauss_pts = getGaussPts().size2();
    const double area = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    for (const auto &m : *elecBlockSetsPtr) {
      double alphaPart = 0.0;
      if (m.second.electrodeEnts.find(fe_ent) != m.second.electrodeEnts.end()) {

        for (int gg = 0; gg != nb_gauss_pts; gg++) {
          FTensor::Tensor1<double, SPACE_DIM> t_r;
          t_r(i) = t_normal(i);
          auto t_r_normalized = t_r.normalize();
          alphaPart += (t_jump(i) * t_r(i)) * t_w * area;
          ++t_jump;
          ++t_normal;
          ++t_w;
        }
        CHKERR ::VecSetValues(petscVec, 1, &index, &alphaPart, ADD_VALUES); // add the alpha value to the vector

      }
      ++index;
    }

    MoFEMFunctionReturn(0);
  }

protected:
  boost::shared_ptr<MatrixDouble> djumpptr;
  boost::shared_ptr<std::map<int, BlockData>> elecBlockSetsPtr;
  SmartPetscObj<Vec> petscVec;
};

struct OpElectricField : public ForcesAndSourcesCore::UserDataOperator {
  OpElectricField(boost::shared_ptr<MatrixDouble> grad_u_negative,
                  boost::shared_ptr<MatrixDouble> grad_u)
      : ForcesAndSourcesCore::UserDataOperator(NOSPACE, OPLAST),
        gradUNegative(grad_u_negative), gradU(grad_u) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;

    const size_t nb_gauss_pts = getGaussPts().size2();
    gradUNegative->resize(SPACE_DIM, nb_gauss_pts, false);
    gradUNegative->clear();
    auto t_grad_u = getFTensor1FromMat<SPACE_DIM>(*gradU);
    auto t_negative_grad_u = getFTensor1FromMat<SPACE_DIM>(*gradUNegative);
    FTensor::Index<'I', SPACE_DIM> I;

    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      t_negative_grad_u(I) = -t_grad_u(I);
      ++t_grad_u;
      ++t_negative_grad_u;
    }

    MoFEMFunctionReturn(0);
  }

private:
  boost::shared_ptr<MatrixDouble> gradUNegative;
  boost::shared_ptr<MatrixDouble> gradU;
};

#endif // __ELECTROSTATICS_HPP__