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

// #include <MoFEM.hpp>
#include "UnknownInterface.hpp"

namespace MoFEM {

// using EntData = MoFEM::EntitiesFieldData::EntData;

enum ModelHypothesis { TRIDIMENSIONAL, PLANESTRAIN, AXISYMMETRICAL };

template <ModelHypothesis H> struct MFrontEleType;

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

  #ifndef SCHUR_ASSEMBLE
    #define SCHUR_ASSEMBLE 0
  #endif

constexpr AssemblyType MFRONT_AT =
    (SCHUR_ASSEMBLE) ? AssemblyType::BLOCK_SCHUR
                     : AssemblyType::PETSC; //< selected assembly type

template <ModelHypothesis H>
struct MFrontInterface : public UnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  MoFEM::Core &cOre;
  MFrontInterface(const MoFEM::Core &core);

  using DomainEle = typename MFrontEleType<H>::DomainEle;
  using DomainEleOp = typename MFrontEleType<H>::DomainEleOp;
  // using PostProcDomainOnRefinedMesh =
  //     typename MFrontEleType<H>::PostProcDomainOnRefinedMesh;

  static constexpr int DIM = MFrontEleType<H>::SPACE_DIM;

  using OpInternalForce =
      typename FormsIntegrators<DomainEleOp>::template Assembly<MFRONT_AT>::
          template LinearForm<GAUSS>::template OpGradTimesTensor<1, DIM, DIM>;
  using OpAssembleLhsFiniteStrains =
      typename FormsIntegrators<DomainEleOp>::template Assembly<MFRONT_AT>::
          template BiLinearForm<GAUSS>::template OpGradTensorGrad<1, DIM, DIM,
                                                                  1>;
  using OpAssembleLhsSmallStrains =
      typename FormsIntegrators<DomainEleOp>::template Assembly<MFRONT_AT>::
          template BiLinearForm<GAUSS>::template OpGradSymTensorGrad<1, DIM,
                                                                     DIM, 0>;

  // MoFEMErrorCode getCommandLineParameters();

  // MoFEMErrorCode setUpdateElementVariablesOperators();

  // MoFEMErrorCode opFactoryDomainRhs(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);
  // MoFEMErrorCode opFactoryDomainLhs(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);

  // MoFEMErrorCode testOperators();

  // MoFEMErrorCode updateElementVariables(SmartPetscObj<DM> dm, string fe_name);

  // MoFEMErrorCode postProcessElement(int step, SmartPetscObj<DM> dm,
  //                                   string fe_name);

  // MoFEMErrorCode setMonitorPtr(boost::shared_ptr<MoFEM::FEMethod> monitor_ptr) {
  //   MoFEMFunctionBeginHot;
  //   monitorPtr = monitor_ptr;
  //   MoFEMFunctionReturnHot(0);
  // }

private:
  // MoFEM::Interface &mField;
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
  // boost::shared_ptr<DomainEle> updateIntVariablesElePtr;

  string positionField;
  string meshNodeField;

  //   moab::Core mb_postGauss;
  // boost::shared_ptr<moab::Interface> moabGaussIntPtr;

  // boost::shared_ptr<MoFEM::FEMethod> monitorPtr;

  // MoFEMErrorCode
  // setIntegrationRule(boost::shared_ptr<ForcesAndSourcesCore> fe_ptr);
  // MoFEMErrorCode setBaseOperators(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);
};

} // namespace MoFEM

#endif // __MFRONT_INTERFACE_HPP__
