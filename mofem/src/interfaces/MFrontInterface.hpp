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

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>
// #include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/LibrariesManager.hxx"

using namespace mgis;
using namespace mgis::behaviour;

namespace MoFEM {

class MFrontInterfaceBase : public UnknownInterface {
public:
  virtual ~MFrontInterfaceBase() = default;
  virtual MoFEMErrorCode getCommandLineParameters() = 0;
};

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
struct MFrontInterface : public MFrontInterfaceBase {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

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

  MoFEMErrorCode getCommandLineParameters() override;

  // MoFEMErrorCode setUpdateElementVariablesOperators();

  // MoFEMErrorCode opFactoryDomainRhs(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);
  // MoFEMErrorCode opFactoryDomainLhs(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);

  // MoFEMErrorCode testOperators();

  // MoFEMErrorCode updateElementVariables(SmartPetscObj<DM> dm, string
  // fe_name);

  // MoFEMErrorCode postProcessElement(int step, SmartPetscObj<DM> dm,
  //                                   string fe_name);

  // MoFEMErrorCode setMonitorPtr(boost::shared_ptr<MoFEM::FEMethod>
  // monitor_ptr) {
  //   MoFEMFunctionBeginHot;
  //   monitorPtr = monitor_ptr;
  //   MoFEMFunctionReturnHot(0);
  // }

enum DataTags { RHS = 0, LHS };

struct BlockData {
  int iD;
  int oRder;

  bool isFiniteStrain;
  string behaviourPath;
  string behaviourName;

  boost::shared_ptr<Behaviour> mGisBehaviour;
  BehaviourDataView bView;
  boost::shared_ptr<BehaviourData> behDataPtr;

  int sizeIntVar;
  int sizeExtVar;
  int sizeGradVar;
  int sizeStressVar;

  vector<double> params;

  double dIssipation;
  double storedEnergy;
  double externalVariable;

  Range eNts;

  BlockData()
      : oRder(-1), isFiniteStrain(false), behaviourPath("src/libBehaviour.so"),
        behaviourName("IsotropicLinearHardeningPlasticity") {
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

  MoFEMErrorCode setBlockBehaviourData(bool set_params_from_blocks) {
    MoFEMFunctionBeginHot;
    if (mGisBehaviour) {

      auto &mgis_bv = *mGisBehaviour;

      sizeIntVar = getArraySize(mgis_bv.isvs, mgis_bv.hypothesis);
      sizeExtVar = getArraySize(mgis_bv.esvs, mgis_bv.hypothesis);
      sizeGradVar = getArraySize(mgis_bv.gradients, mgis_bv.hypothesis);
      sizeStressVar =
          getArraySize(mgis_bv.thermodynamic_forces, mgis_bv.hypothesis);

      behDataPtr = boost::make_shared<BehaviourData>(BehaviourData{mgis_bv});
      bView = make_view(*behDataPtr);
      const int total_number_of_params = mgis_bv.mps.size();
      // const int total_number_of_params = mgis_bv.mps.size() +
      // mgis_bv.params.size() + mgis_bv.iparams.size() +
      // mgis_bv.usparams.size();

      if (set_params_from_blocks) {

        if (params.size() < total_number_of_params)
          SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Not enough parameters supplied for this block. We have %d "
                   "provided where %d are necessary for this block",
                   params.size(), total_number_of_params);

        for (int dd = 0; dd < total_number_of_params; ++dd) {
          setMaterialProperty(behDataPtr->s0, dd, params[dd]);
          setMaterialProperty(behDataPtr->s1, dd, params[dd]);
        }
      }

      if (isFiniteStrain) {
        behDataPtr->K[0] = 0; // no tangent
        behDataPtr->K[1] = 2; // PK1
        behDataPtr->K[2] = 2; // PK1
      } else {
        behDataPtr->K[0] = 0; // no tangent
        behDataPtr->K[1] = 0; // cauchy
      }

      for (auto &mb : {&behDataPtr->s0, &behDataPtr->s1}) {
        mb->dissipated_energy = dIssipation;
        mb->stored_energy = storedEnergy;
        setExternalStateVariable(*mb, 0, externalVariable);
      }
    }

    MoFEMFunctionReturnHot(0);
  }
};

  struct CommonData {

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

    CommonData(MoFEM::Interface &m_field) : mField(m_field) {}

    MoFEMErrorCode setBlocks(int dim) {
      MoFEMFunctionBegin;
      string block_name = "MFRONT";
      for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField, BLOCKSET, it)) {
        if (it->getName().compare(0, block_name.size(), block_name) == 0) {
          std::vector<double> block_data;
          // FIXME: TODO: maybe this should be set only from the command line!!!
          CHKERR it->getAttributes(block_data);
          const int id = it->getMeshsetId();
          EntityHandle meshset = it->getMeshset();
          CHKERR mField.get_moab().get_entities_by_dimension(
              meshset, dim, setOfBlocksData[id].eNts, true);
          for (auto ent : setOfBlocksData[id].eNts)
            blocksIDmap[ent] = id;

          setOfBlocksData[id].iD = id;
          setOfBlocksData[id].params.resize(block_data.size());

          for (int n = 0; n != block_data.size(); n++)
            setOfBlocksData[id].params[n] = block_data[n];
        }
      }

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode getInternalVar(const EntityHandle fe_ent,
                                  const int nb_gauss_pts, const int var_size,
                                  const int grad_size, const int stress_size,
                                  bool is_large_strain = false) {
      MoFEMFunctionBegin;

      auto mget_tag_data = [&](Tag &m_tag,
                               boost::shared_ptr<MatrixDouble> &m_mat,
                               const int &m_size, bool is_def_grad = false) {
        MoFEMFunctionBeginHot;

        double *tag_data;
        int tag_size;
        rval = mField.get_moab().tag_get_by_ptr(
            m_tag, &fe_ent, 1, (const void **)&tag_data, &tag_size);

        if (rval != MB_SUCCESS || tag_size != m_size * nb_gauss_pts) {
          m_mat->resize(nb_gauss_pts, m_size, false);
          m_mat->clear();
          // initialize deformation gradient properly
          if (is_def_grad && is_large_strain)
            for (int gg = 0; gg != nb_gauss_pts; ++gg) {
              (*m_mat)(gg, 0) = 1;
              (*m_mat)(gg, 1) = 1;
              (*m_mat)(gg, 2) = 1;
            }
          void const *tag_data2[] = {&*m_mat->data().begin()};
          const int tag_size2 = m_mat->data().size();
          CHKERR mField.get_moab().tag_set_by_ptr(m_tag, &fe_ent, 1, tag_data2,
                                                  &tag_size2);
        } else {
          MatrixAdaptor tag_vec = MatrixAdaptor(
              nb_gauss_pts, m_size,
              ublas::shallow_array_adaptor<double>(tag_size, tag_data));

          *m_mat = tag_vec;
        }

        MoFEMFunctionReturnHot(0);
      };

      CHKERR mget_tag_data(internalVariableTag, internalVariablePtr, var_size);
      CHKERR mget_tag_data(stressTag, mPrevStressPtr, grad_size);
      CHKERR mget_tag_data(gradientTag, mPrevGradPtr, stress_size, true);

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode setInternalVar(const EntityHandle fe_ent) {
      MoFEMFunctionBegin;

      auto mset_tag_data = [&](Tag &m_tag,
                               boost::shared_ptr<MatrixDouble> &m_mat) {
        MoFEMFunctionBeginHot;
        void const *tag_data[] = {&*m_mat->data().begin()};
        const int tag_size = m_mat->data().size();
        CHKERR mField.get_moab().tag_set_by_ptr(m_tag, &fe_ent, 1, tag_data,
                                                &tag_size);
        MoFEMFunctionReturnHot(0);
      };

      CHKERR mset_tag_data(internalVariableTag, internalVariablePtr);
      CHKERR mset_tag_data(stressTag, mPrevStressPtr);
      CHKERR mset_tag_data(gradientTag, mPrevGradPtr);

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode createTags() {
      MoFEMFunctionBegin;
      double def_val = 0.0;
      const int default_length = 0;
      CHKERR mField.get_moab().tag_get_handle(
          "_INTERNAL_VAR", default_length, MB_TYPE_DOUBLE, internalVariableTag,
          MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);
      CHKERR mField.get_moab().tag_get_handle(
          "_STRESS_TAG", default_length, MB_TYPE_DOUBLE, stressTag,
          MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);
      CHKERR mField.get_moab().tag_get_handle(
          "_GRAD_TAG", default_length, MB_TYPE_DOUBLE, gradientTag,
          MB_TAG_CREAT | MB_TAG_VARLEN | MB_TAG_SPARSE, PETSC_NULL);

      MoFEMFunctionReturn(0);
    }

    MoFEMErrorCode clearTags() {
      MoFEMFunctionBegin;
      double zero = 0;
      for (auto &[id, data] : setOfBlocksData) {
        CHKERR mField.get_moab().tag_clear_data(internalVariableTag, data.eNts,
                                                &zero);
        CHKERR mField.get_moab().tag_clear_data(stressTag, data.eNts, &zero);
        CHKERR mField.get_moab().tag_clear_data(gradientTag, data.eNts, &zero);
      }
      MoFEMFunctionReturn(0);
    }
  };

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
  // boost::shared_ptr<DomainEle> updateIntVariablesElePtr;

  // string positionField;
  // string meshNodeField;

  //   moab::Core mb_postGauss;
  boost::shared_ptr<moab::Interface> moabGaussIntPtr;

  boost::shared_ptr<MoFEM::FEMethod> monitorPtr;

  boost::shared_ptr<CommonData> commonDataPtr;

  // MoFEMErrorCode
  // setIntegrationRule(boost::shared_ptr<ForcesAndSourcesCore> fe_ptr);
  // MoFEMErrorCode setBaseOperators(
  //     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip);
};

} // namespace MoFEM

#endif // __MFRONT_INTERFACE_HPP__
