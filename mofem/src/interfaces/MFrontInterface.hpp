/** \file MFrontInterface.hpp
 * \brief MFrontInterface
 *
 * MFrontInterface
 *
 *
 * \ingroup mfront_interface
 */

#ifndef __MFRONT_INTERFACE_HPP__
#define __MFRONT_INTERFACE_HPP__

#include "UnknownInterface.hpp"

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>

namespace MoFEM {

class MFrontInterfaceBase : public UnknownInterface {
public:
  virtual ~MFrontInterfaceBase() = default;
  virtual MoFEMErrorCode getCommandLineParameters() = 0;
  virtual MoFEMErrorCode opFactoryDomainRhs(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name) = 0;
  virtual MoFEMErrorCode opFactoryDomainLhs(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name) = 0;
  virtual MoFEMErrorCode
  setUpdateElementVariablesOperators(ForcesAndSourcesCore::RuleHookFun rule,
                                     std::string field_name) = 0;
  virtual MoFEMErrorCode updateElementVariables(SmartPetscObj<DM> dm,
                                                std::string fe_name) = 0;

  enum DataTags { RHS = 0, LHS };

  struct BlockData {

    BlockData()
        : oRder(-1), isFiniteStrain(false),
          behaviourPath("src/libBehaviour.so"),
          behaviourName("LinearElasticity") {
      dIssipation = 0;
      storedEnergy = 0;
      externalVariable = 0;
    }

    inline MoFEMErrorCode setTag(DataTags tag) {
      MoFEMFunctionBeginHot;
      if (tag == RHS) {
        behDataPtr->K[0] = 0;
      } else {
        behDataPtr->K[0] = 5;
      }
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode setBlockBehaviourData(bool set_params_from_blocks);

    int iD;
    int oRder;

    bool isFiniteStrain;
    string behaviourPath;
    string behaviourName;

    boost::shared_ptr<mgis::behaviour::Behaviour> mGisBehaviour;
    mgis::behaviour::BehaviourDataView bView;
    boost::shared_ptr<mgis::behaviour::BehaviourData> behDataPtr;

    int sizeIntVar;
    int sizeExtVar;
    int sizeGradVar;
    int sizeStressVar;

    vector<double> params;

    double dIssipation;
    double storedEnergy;
    double externalVariable;

    Range eNts;
  };

  struct CommonData {
    CommonData(MoFEM::Interface &m_field) : mField(m_field) {};

    MoFEMErrorCode setBlocks(int dim);
    MoFEMErrorCode getInternalVar(const EntityHandle fe_ent,
                                  const int nb_gauss_pts, const int var_size,
                                  const int grad_size, const int stress_size,
                                  bool is_large_strain = false);

    MoFEMErrorCode setInternalVar(const EntityHandle fe_ent);
    MoFEMErrorCode createTags();
    MoFEMErrorCode clearTags();

    MoFEM::Interface &mField;
    boost::shared_ptr<MatrixDouble> mGradPtr;
    boost::shared_ptr<MatrixDouble> mStressPtr;
    boost::shared_ptr<MatrixDouble> mFullStrainPtr;
    boost::shared_ptr<MatrixDouble> mFullStressPtr;

    boost::shared_ptr<MatrixDouble> mPrevGradPtr;
    boost::shared_ptr<MatrixDouble> mPrevStressPtr;

    boost::shared_ptr<MatrixDouble> mDispPtr;
    boost::shared_ptr<MatrixDouble> materialTangentPtr;
    boost::shared_ptr<MatrixDouble> mFullTangentPtr;
    boost::shared_ptr<MatrixDouble> internalVariablePtr;

    std::map<int, BlockData> setOfBlocksData;
    std::map<EntityHandle, int> blocksIDmap;

    MatrixDouble lsMat;
    MatrixDouble gaussPts;
    MatrixDouble myN;

    Tag internalVariableTag;
    Tag stressTag;
    Tag gradientTag;
  };

  // template <bool IS_LARGE_STRAIN, ModelHypothesis MH>
  // struct OpSaveStress : public MFrontEleType<MH>::DomainEleOp {
  //   static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;

  //   OpSaveStress(const std::string field_name,
  //                boost::shared_ptr<CommonData> common_data_ptr)
  //       : MFrontEleType<MH>::DomainEleOp(field_name,
  //                                       MFrontEleType<MH>::DomainEleOp::OPROW),
  //         commonDataPtr(common_data_ptr) {
  //     std::fill(&MFrontEleType<MH>::DomainEleOp::doEntities[MBEDGE],
  //               &MFrontEleType<MH>::DomainEleOp::doEntities[MBMAXTYPE],
  //               false);
  //   }
  //   MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

  // private:
  //   boost::shared_ptr<CommonData> commonDataPtr;
  // };
};

using EntData = EntitiesFieldData::EntData;

enum ModelHypothesis { TRIDIMENSIONAL, PLANESTRAIN, AXISYMMETRICAL };

template <ModelHypothesis MH> struct MFrontEleType;

template <> struct MFrontEleType<TRIDIMENSIONAL> {

  MFrontEleType() = delete;
  ~MFrontEleType() = delete;

  using DomainEle = MoFEM::VolumeElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
  // using PostProcDomainOnRefinedMesh = PostProcVolumeOnRefinedMesh;

  static constexpr int SPACE_DIM = 3;
};

template <> struct MFrontEleType<PLANESTRAIN> {

  MFrontEleType() = delete;
  ~MFrontEleType() = delete;

  using DomainEle = MoFEM::FaceElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
  // using PostProcDomainOnRefinedMesh = PostProcFaceOnRefinedMesh;

  static constexpr int SPACE_DIM = 2;
};

template <> struct MFrontEleType<AXISYMMETRICAL> {

  MFrontEleType() = delete;
  ~MFrontEleType() = delete;

  using DomainEle = FaceElementForcesAndSourcesCore;
  using DomainEleOp = DomainEle::UserDataOperator;
  // using PostProcDomainOnRefinedMesh = PostProcFaceOnRefinedMesh;

  static constexpr int SPACE_DIM = 2;
};

template <ModelHypothesis MH, AssemblyType AT = AssemblyType::PETSC>
struct MFrontInterface : public MFrontInterfaceBase {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MFrontInterface(const MoFEM::Core &core);

  using DomainEle = typename MFrontEleType<MH>::DomainEle;
  using DomainEleOp = typename MFrontEleType<MH>::DomainEleOp;
  // using PostProcDomainOnRefinedMesh =
  //     typename MFrontEleType<MH>::PostProcDomainOnRefinedMesh;

  static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;

  using OpInternalForce =
      typename FormsIntegrators<DomainEleOp>::template Assembly<AT>::
          template LinearForm<GAUSS>::template OpGradTimesTensor<1, DIM, DIM>;
  using OpAssembleLhsFiniteStrains =
      typename FormsIntegrators<DomainEleOp>::template Assembly<
          AT>::template BiLinearForm<GAUSS>::template OpGradTensorGrad<1, DIM,
                                                                       DIM, 1>;
  using OpAssembleLhsSmallStrains =
      typename FormsIntegrators<DomainEleOp>::template Assembly<AT>::
          template BiLinearForm<GAUSS>::template OpGradSymTensorGrad<1, DIM,
                                                                     DIM, 0>;

  MoFEMErrorCode getCommandLineParameters() override;

  MoFEMErrorCode opFactoryDomainRhs(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name) override;

  MoFEMErrorCode opFactoryDomainLhs(
      boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
      std::string field_name) override;

  // // MoFEMErrorCode testOperators();

  MoFEMErrorCode
  setUpdateElementVariablesOperators(ForcesAndSourcesCore::RuleHookFun rule,
                                     std::string field_name) override;

  MoFEMErrorCode updateElementVariables(SmartPetscObj<DM> dm,
                                        std::string fe_name) override;

  // MoFEMErrorCode postProcessElement(int step, SmartPetscObj<DM> dm,
  //                                   string fe_name);

  MoFEMErrorCode setMonitorPtr(boost::shared_ptr<MoFEM::FEMethod> monitor_ptr) {
    MoFEMFunctionBeginHot;
    monitorPtr = monitor_ptr;
    MoFEMFunctionReturnHot(0);
  }

private:
  MoFEM::Core &cOre;

  string optionsPrefix;

  SmartPetscObj<DM> dM;

  // PetscBool isQuasiStatic;
  PetscBool saveGauss;
  PetscBool saveVolume;

  PetscBool testJacobian;
  PetscReal randomFieldScale;

  // PetscInt oRder;
  // bool isDisplacementField;
  bool isFiniteKinematics;
  // BitRefLevel bIt;

  FieldApproximationBase fieldBase;

  // boost::shared_ptr<PostProcDomainOnRefinedMesh> postProcFe;
  boost::shared_ptr<DomainEle> updateIntVariablesElePtr;

  //   moab::Core mb_postGauss;
  boost::shared_ptr<moab::Interface> moabGaussIntPtr;

  boost::shared_ptr<MoFEM::FEMethod> monitorPtr;

  boost::shared_ptr<CommonData> commonDataPtr;
};

// template <ModelHypothesis MH>
// struct OpSaveGaussPts : public MFrontEleType<MH>::DomainEleOp {
//   static constexpr int DIM = MFrontEleType<MH>::SPACE_DIM;

//   OpSaveGaussPts(const std::string field_name, moab::Interface &moab_mesh,
//                  boost::shared_ptr<CommonData> common_data_ptr)
//       : MFrontEleType<MH>::DomainEleOp(field_name,
//                                       MFrontEleType<MH>::DomainEleOp::OPROW),
//         internalVarMesh(moab_mesh), commonDataPtr(common_data_ptr) {
//     // Operator is only executed for vertices
//     std::fill(&MFrontEleType<MH>::DomainEleOp::doEntities[MBEDGE],
//               &MFrontEleType<MH>::DomainEleOp::doEntities[MBMAXTYPE], false);
//   }
//   MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

// private:
//   boost::shared_ptr<CommonData> commonDataPtr;
//   moab::Interface &internalVarMesh;
// };

} // namespace MoFEM

#endif // __MFRONT_INTERFACE_HPP__
