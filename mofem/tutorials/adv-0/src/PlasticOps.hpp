

/** \file PlasticOps.hpp
 * \example PlasticOps.hpp

\f[
\left\{
\begin{array}{ll}
\frac{\partial \sigma_{ij}}{\partial x_j} - b_i = 0 & \forall x \in \Omega \\
\varepsilon_{ij} = \frac{1}{2}\left( \frac{\partial u_i}{\partial x_j} +
\frac{\partial u_j}{\partial x_i} \right)\\
\sigma_{ij} = D_{ijkl}\left(\varepsilon_{kl}-\varepsilon^p_{kl}\right) \\
\dot{\varepsilon}^p_{kl} - \dot{\tau} \left( \left. \frac{\partial f}{\partial
\sigma_{kl}} \right|_{(\sigma,\tau) } \right) = 0 \\
f(\sigma, \tau) \leq 0,\; \dot{\tau} \geq 0,\;\dot{\tau}f(\sigma, \tau)=0\\
u_i = \overline{u}_i & \forall x \in \partial\Omega_u \\
\sigma_{ij}n_j = \overline{t}_i & \forall x \in \partial\Omega_\sigma \\
\Omega_u \cup \Omega_\sigma = \Omega \\
\Omega_u \cap \Omega_\sigma = \emptyset
\end{array}
\right.
\f]

\f[
\left\{
\begin{array}{ll}
\left(\frac{\partial \delta u_i}{\partial x_j},\sigma_{ij}\right)_\Omega-(\delta
u_i,b_i)_\Omega -(\delta u_i,\overline{t}_i)_{\partial\Omega_\sigma}=0 & \forall
\delta u_i \in H^1(\Omega)\\ \left(\delta\varepsilon^p_{kl} ,D_{ijkl}\left(
\dot{\varepsilon}^p_{kl} - \dot{\tau} A_{kl} \right)\right) = 0
& \forall \delta\varepsilon^p_{ij} \in L^2(\Omega) \cap \mathcal{S} \\
\left(\delta\tau,c_n\dot{\tau} - \frac{1}{2}\left\{c_n \dot{\tau} +
(f(\pmb\sigma,\tau) - \sigma_y) +
\| c_n \dot{\tau} + (f(\pmb\sigma,\tau) - \sigma_y) \|\right\}\right) = 0 &
\forall \delta\tau \in L^2(\Omega) \end{array} \right.
\f]

*/

namespace PlasticOps {

//! [Common data]
struct CommonData : public boost::enable_shared_from_this<CommonData> {

  enum ParamsIndexes {
    YOUNG_MODULUS,
    POISSON_RATIO,
    SIGMA_Y,
    H,
    VIS_H,
    QINF,
    BISO,
    C1_k,
    LAST_PARAM
  };

  using BlockParams = std::array<double, LAST_PARAM>;
  BlockParams blockParams;

  inline auto getParamsPtr() {
    return boost::shared_ptr<BlockParams>(shared_from_this(), &blockParams);
  };

  //! [Common data set externally]
  boost::shared_ptr<MatrixDouble> mDPtr;
  boost::shared_ptr<MatrixDouble> mGradPtr;
  boost::shared_ptr<MatrixDouble> mStrainPtr;
  boost::shared_ptr<MatrixDouble> mStressPtr;
  // Advection
  boost::shared_ptr<MatrixDouble> guidingVelocityPtr;
  boost::shared_ptr<VectorDouble> velocityDotNormalPtr;
  //
  //! [Common data set externally]

  VectorDouble plasticSurface;
  MatrixDouble plasticFlow;
  VectorDouble plasticTau;
  VectorDouble plasticTauDotPartial;
  VectorDouble plasticTauDot;
  MatrixDouble plasticStrain;
  MatrixDouble plasticStrainDot;

  MatrixDouble plasticGradTau;
  MatrixDouble plasticGradStrain;

  VectorDouble resC;
  VectorDouble resCdTau;
  VectorDouble resCdTauAle;
  MatrixDouble resCdStrain;
  MatrixDouble resCdPlasticStrain;
  MatrixDouble resFlow;
  MatrixDouble resFlowDtau;
  MatrixDouble resFlowDtauAle;
  MatrixDouble resFlowDstrain;
  MatrixDouble resFlowDstrainAle;
  MatrixDouble resFlowDstrainDot;

  inline auto getPlasticSurfacePtr() {
    return boost::shared_ptr<VectorDouble>(shared_from_this(), &plasticSurface);
  }
  inline auto getPlasticTauPtr() {
    return boost::shared_ptr<VectorDouble>(shared_from_this(), &plasticTau);
  }
  inline auto getPlasticGradTauPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &plasticGradTau);
  }
  inline auto getPlasticTauDotPtr() {
    return boost::shared_ptr<VectorDouble>(shared_from_this(), &plasticTauDot);
  }
  inline auto getPlasticStrainPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &plasticStrain);
  }
  inline auto getPlasticGradStrainPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &plasticGradStrain);
  }
  inline auto getPlasticStrainDotPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                           &plasticStrainDot);
  }
  inline auto getPlasticFlowPtr() {
    return boost::shared_ptr<MatrixDouble>(shared_from_this(), &plasticFlow);
  }

  static std::array<int, 5> activityData;
};

std::array<int, 5> CommonData::activityData = {0, 0, 0, 0, 0};

//! [Common data]

FTensor::Index<'I', 3> I;
FTensor::Index<'J', 3> J;
FTensor::Index<'M', 3> M;
FTensor::Index<'N', 3> N;

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculatePlasticSurfaceImpl;

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpCalculatePlasticityImpl;

template <int DIM, IntegrationType I, typename DomainEleOp>
struct OpPlasticStressImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticFlowRhsImpl;

template <IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateConstraintsRhsImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticInternalForceLhs_dEPImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticInternalForceLhs_LogStrain_dEPImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticFlowLhs_dUImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticFlowLhs_LogStrain_dUImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticFlowLhs_dEPImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculatePlasticFlowLhs_dTAUImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateConstraintsLhs_dUImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateConstraintsLhs_LogStrain_dUImpl;

template <int DIM, IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateConstraintsLhs_dEPImpl;

template <IntegrationType I, typename AssemblyDomainEleOp>
struct OpCalculateConstraintsLhs_dTAUImpl;

struct SideData {
  // data for skeleton computation
  std::array<EntityHandle, 2> feSideHandle;

  // plastic constraint
  std::array<VectorInt, 2>
      indicesRowSideMapTau; ///< indices on rows for left hand-side
  std::array<VectorInt, 2>
      indicesColSideMapTau; ///< indices on columns for left hand-side
  std::array<MatrixDouble, 2> rowBaseSideMapTau; // base functions on rows
  std::array<MatrixDouble, 2> colBaseSideMapTau; // base function  on columns

  // plastic Strain -- sizes?
  std::array<VectorInt, 2>
      indicesRowSideMapEp; ///< indices on rows for left hand-side
  std::array<VectorInt, 2>
      indicesColSideMapEp; ///< indices on columns for left hand-side
  std::array<MatrixDouble, 2> rowBaseSideMapEp; // base functions on rows
  std::array<MatrixDouble, 2> colBaseSideMapEp; // base function  on columns
  std::array<int, 2> senseMap; // orientation of local element edge/face in
                               // respect to global orientation of edge/face
  std::array<MatrixDouble, 2> epMat;  //< Values of plastic strain
  std::array<VectorDouble, 2> tauVec; //< Values of plastic flow
  // std::array<MatrixDouble, 2> velMat; //< Values of velocity field

  int currentFESide; ///< current side counter
};

boost::shared_ptr<SideEle>
getSideFE(boost::shared_ptr<SideData> side_data_ptr);

struct OpCalculatePlasticConvRotatingFrame : public DomainEleOp {
  OpCalculatePlasticConvRotatingFrame(
      const std::string field_name,
      boost::shared_ptr<CommonData> common_data_ptr);
  MoFEMErrorCode doWork(int side, EntityType type, EntData &data);

protected:
  boost::shared_ptr<CommonData> commonDataPtr;
  boost::shared_ptr<SideData> SideDataPtr;
};
struct OpRhsSkeleton : public BoundaryEleOp {

  OpRhsSkeleton(boost::shared_ptr<SideData> side_data_ptr,
                boost::shared_ptr<SideEle> side_fe_ptr);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<SideData> sideDataPtr;
  boost::shared_ptr<SideEle> sideFEPtr;

  VectorDouble resSkelton;
};

struct OpLhsSkeleton : public BoundaryEleOp {
  OpLhsSkeleton(boost::shared_ptr<SideData> side_data_ptr,
                boost::shared_ptr<SideEle> side_fe_ptr);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<SideData> sideDataPtr;
  boost::shared_ptr<SideEle> sideFEPtr;

  MatrixDouble matSkeleton;
};

struct OpRhsSkeletonEp : public BoundaryEleOp {

  OpRhsSkeletonEp(boost::shared_ptr<SideData> side_data_ptr,
                boost::shared_ptr<SideEle> side_fe_ptr);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<SideData> sideDataPtr;
  boost::shared_ptr<SideEle> sideFEPtr;

  VectorDouble resSkelton;
};

struct OpLhsSkeletonEp : public BoundaryEleOp {
  OpLhsSkeletonEp(boost::shared_ptr<SideData> side_data_ptr,
                boost::shared_ptr<SideEle> side_fe_ptr);
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  boost::shared_ptr<SideData> sideDataPtr;
  boost::shared_ptr<SideEle> sideFEPtr;

  MatrixDouble matSkeleton;
};

template <typename DomainEleOp> struct PlasticityIntegrators {
  template <int DIM, IntegrationType I>
  using OpCalculatePlasticSurface =
      OpCalculatePlasticSurfaceImpl<DIM, I, DomainEleOp>;

  template <int DIM, IntegrationType I>
  using OpCalculatePlasticity = OpCalculatePlasticityImpl<DIM, I, DomainEleOp>;

  template <int DIM, IntegrationType I>
  using OpPlasticStress = OpPlasticStressImpl<DIM, I, DomainEleOp>;

  template <AssemblyType A> struct Assembly {

    using AssemblyDomainEleOp =
        typename FormsIntegrators<DomainEleOp>::template Assembly<A>::OpBase;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticFlowRhs =
        OpCalculatePlasticFlowRhsImpl<DIM, I, AssemblyDomainEleOp>;

    template <IntegrationType I>
    using OpCalculateConstraintsRhs =
        OpCalculateConstraintsRhsImpl<I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticInternalForceLhs_dEP =
        OpCalculatePlasticInternalForceLhs_dEPImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticInternalForceLhs_LogStrain_dEP =
        OpCalculatePlasticInternalForceLhs_LogStrain_dEPImpl<
            DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticFlowLhs_dU =
        OpCalculatePlasticFlowLhs_dUImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticFlowLhs_LogStrain_dU =
        OpCalculatePlasticFlowLhs_LogStrain_dUImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticFlowLhs_dEP =
        OpCalculatePlasticFlowLhs_dEPImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculatePlasticFlowLhs_dTAU =
        OpCalculatePlasticFlowLhs_dTAUImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculateConstraintsLhs_dU =
        OpCalculateConstraintsLhs_dUImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculateConstraintsLhs_LogStrain_dU =
        OpCalculateConstraintsLhs_LogStrain_dUImpl<DIM, I, AssemblyDomainEleOp>;

    template <int DIM, IntegrationType I>
    using OpCalculateConstraintsLhs_dEP =
        OpCalculateConstraintsLhs_dEPImpl<DIM, I, AssemblyDomainEleOp>;

    template <IntegrationType I>
    using OpCalculateConstraintsLhs_dTAU =
        OpCalculateConstraintsLhs_dTAUImpl<I, AssemblyDomainEleOp>;
  };
};

}; // namespace PlasticOps

#include <PlasticOpsGeneric.hpp>
#include <PlasticOpsSmallStrains.hpp>
#include <PlasticOpsLargeStrains.hpp>
#include <PlasticOpsMonitor.hpp>

namespace PlasticOps {

using Pip = boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator>;
using CommonPlasticPtr = boost::shared_ptr<PlasticOps::CommonData>;
using CommonHenckyPtr = boost::shared_ptr<HenckyOps::CommonData>;

template <int DIM>
MoFEMErrorCode
addMatBlockOps(MoFEM::Interface &m_field, std::string block_name, Pip &pip,
               boost::shared_ptr<MatrixDouble> mat_D_Ptr,
               boost::shared_ptr<CommonData::BlockParams> mat_params_ptr,
               double scale_value, Sev sev) {
  MoFEMFunctionBegin;

  struct OpMatBlocks : public DomainEleOp {
    OpMatBlocks(boost::shared_ptr<MatrixDouble> m_D_ptr,
                boost::shared_ptr<CommonData::BlockParams> mat_params_ptr,
                double scale_value, MoFEM::Interface &m_field, Sev sev,
                std::vector<const CubitMeshSets *> meshset_vec_ptr)
        : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE), matDPtr(m_D_ptr),
          matParamsPtr(mat_params_ptr), scaleVal(scale_value) {
      CHK_THROW_MESSAGE(extractBlockData(m_field, meshset_vec_ptr, sev),
                        "Can not get data from block");
    }

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBegin;

      auto getK = [](auto &p) {
        auto young_modulus = p[CommonData::YOUNG_MODULUS];
        auto poisson_ratio = p[CommonData::POISSON_RATIO];
        return young_modulus / (3 * (1 - 2 * poisson_ratio));
      };

      auto getG = [](auto &p) {
        auto young_modulus = p[CommonData::YOUNG_MODULUS];
        auto poisson_ratio = p[CommonData::POISSON_RATIO];
        return young_modulus / (2 * (1 + poisson_ratio));
      };

      auto scale_fun = [this](auto &p) {
        p[CommonData::YOUNG_MODULUS] *= scaleVal;
        p[CommonData::SIGMA_Y] *= scaleVal;
        p[CommonData::H] *= scaleVal;
        p[CommonData::VIS_H] *= scaleVal;
        p[CommonData::QINF] *= scaleVal;
        p[CommonData::C1_k] *= scaleVal;
      };

      for (auto &b : blockData) {
        if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
          *matParamsPtr = b.bParams;
          scale_fun(*matParamsPtr);
          CHKERR getMatDPtr(matDPtr, getK(*matParamsPtr), getG(*matParamsPtr));
          MoFEMFunctionReturnHot(0);
        }
      }

      (*matParamsPtr) = {young_modulus, poisson_ratio, sigmaY, H,
                         visH,          Qinf,          b_iso,  C1_k};
      scale_fun(*matParamsPtr);
      CHKERR getMatDPtr(matDPtr, getK(*matParamsPtr), getG(*matParamsPtr));

      MoFEMFunctionReturn(0);
    }

  private:
    boost::shared_ptr<MatrixDouble> matDPtr;
    boost::shared_ptr<CommonData::BlockParams> matParamsPtr;
    const double scaleVal;

    struct BlockData {
      std::array<double, CommonData::LAST_PARAM> bParams;
      Range blockEnts;
    };
    std::vector<BlockData> blockData;

    /**
     * @brief Extract block data from meshsets
     *
     * @param m_field
     * @param meshset_vec_ptr
     * @param sev
     * @return MoFEMErrorCode
     */
    MoFEMErrorCode
    extractBlockData(MoFEM::Interface &m_field,
                     std::vector<const CubitMeshSets *> meshset_vec_ptr,
                     Sev sev) {
      MoFEMFunctionBegin;

      for (auto m : meshset_vec_ptr) {
        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock") << *m;
        std::vector<double> block_data;
        CHKERR m->getAttributes(block_data);
        if (block_data.size() != CommonData::LAST_PARAM) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Wrong number of block attribute");
        }
        auto get_block_ents = [&]() {
          Range ents;
          CHKERR m_field.get_moab().get_entities_by_handle(m->meshset, ents,
                                                           true);
          return ents;
        };

        CommonData::BlockParams block_params;
        for (auto i = 0; i != CommonData::LAST_PARAM; ++i) {
          block_params[i] = block_data[i];
        }

        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock")
            << "E = " << block_params[CommonData::YOUNG_MODULUS]
            << " nu = " << block_params[CommonData::POISSON_RATIO];
        MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock")
            << std::endl
            << "sigma_y = " << block_params[CommonData::SIGMA_Y] << std::endl
            << "h = " << block_params[CommonData::H] << std::endl
            << "vis_h = " << block_params[CommonData::VIS_H] << std::endl
            << "qinf = " << block_params[CommonData::QINF] << std::endl
            << "biso = " << block_params[CommonData::BISO] << std::endl
            << "C1_k = " << block_params[CommonData::C1_k] << std::endl;

        blockData.push_back({block_params, get_block_ents()});
      }
      MOFEM_LOG_CHANNEL("WORLD");
      MoFEMFunctionReturn(0);
    }

    /**
     * @brief Get elasticity tensor
     *
     * Calculate elasticity tensor for given material parameters
     *
     * @param mat_D_ptr
     * @param bulk_modulus_K
     * @param shear_modulus_G
     * @return MoFEMErrorCode
     *
     */
    MoFEMErrorCode getMatDPtr(boost::shared_ptr<MatrixDouble> mat_D_ptr,
                              double bulk_modulus_K, double shear_modulus_G) {
      MoFEMFunctionBegin;
      //! [Calculate elasticity tensor]
      auto set_material_stiffness = [&]() {
        FTensor::Index<'i', DIM> i;
        FTensor::Index<'j', DIM> j;
        FTensor::Index<'k', DIM> k;
        FTensor::Index<'l', DIM> l;
        constexpr auto t_kd = FTensor::Kronecker_Delta_symmetric<int>();
        double A = (DIM == 2)
                       ? 2 * shear_modulus_G /
                             (bulk_modulus_K + (4. / 3.) * shear_modulus_G)
                       : 1;
        auto t_D = getFTensor4DdgFromMat<DIM, DIM, 0>(*mat_D_ptr);
        t_D(i, j, k, l) =
            2 * shear_modulus_G * ((t_kd(i, k) ^ t_kd(j, l)) / 4.) +
            A * (bulk_modulus_K - (2. / 3.) * shear_modulus_G) * t_kd(i, j) *
                t_kd(k, l);
      };
      //! [Calculate elasticity tensor]
      constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
      mat_D_ptr->resize(size_symm * size_symm, 1);
      set_material_stiffness();
      MoFEMFunctionReturn(0);
    }
  };

  // push operator to calculate material stiffness matrix for each block
  pip.push_back(new OpMatBlocks(
      mat_D_Ptr, mat_params_ptr, scale_value, m_field, sev,

      // Get blockset using regular expression
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

          (boost::format("%s(.*)") % block_name).str()

              ))

          ));

  MoFEMFunctionReturn(0);
}

template <int DIM, IntegrationType I, typename DomainEleOp>
auto createCommonPlasticOps(
    MoFEM::Interface &m_field, std::string block_name,
    boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
    std::string u, std::string ep, std::string tau, double scale, Sev sev) {

  using P = PlasticityIntegrators<DomainEleOp>;

  auto common_plastic_ptr = boost::make_shared<PlasticOps::CommonData>();
  common_plastic_ptr = boost::make_shared<PlasticOps::CommonData>();

  constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
  auto make_d_mat = []() {
    return boost::make_shared<MatrixDouble>(size_symm * size_symm, 1);
  };

  common_plastic_ptr->mDPtr = make_d_mat();
  common_plastic_ptr->mGradPtr = boost::make_shared<MatrixDouble>();
  common_plastic_ptr->mStrainPtr = boost::make_shared<MatrixDouble>();
  common_plastic_ptr->mStressPtr = boost::make_shared<MatrixDouble>();
  common_plastic_ptr->guidingVelocityPtr = boost::make_shared<MatrixDouble>();
  common_plastic_ptr->velocityDotNormalPtr = boost::make_shared<VectorDouble>();

  auto m_D_ptr = common_plastic_ptr->mDPtr;

  CHK_THROW_MESSAGE(addMatBlockOps<DIM>(m_field, block_name, pip, m_D_ptr,
                                        common_plastic_ptr->getParamsPtr(),
                                        scale, sev),
                    "add mat block ops");

  pip.push_back(new OpCalculateScalarFieldValues(
      tau, common_plastic_ptr->getPlasticTauPtr()));
  pip.push_back(new OpCalculateTensor2SymmetricFieldValues<DIM>(
      ep, common_plastic_ptr->getPlasticStrainPtr()));
  pip.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(
      u, common_plastic_ptr->mGradPtr));

  CommonHenckyPtr common_hencky_ptr;

  if (is_large_strains) {
    common_hencky_ptr = boost::make_shared<HenckyOps::CommonData>();
    common_hencky_ptr->matGradPtr = common_plastic_ptr->mGradPtr;
    common_hencky_ptr->matDPtr = common_plastic_ptr->mDPtr;
    common_hencky_ptr->matLogCPlastic =
        common_plastic_ptr->getPlasticStrainPtr();
    common_plastic_ptr->mStrainPtr = common_hencky_ptr->getMatLogC();
    common_plastic_ptr->mStressPtr = common_hencky_ptr->getMatHenckyStress();

    using H = HenckyOps::HenckyIntegrators<DomainEleOp>;

    pip.push_back(new typename H::template OpCalculateEigenVals<DIM, I>(
        u, common_hencky_ptr));
    pip.push_back(
        new typename H::template OpCalculateLogC<DIM, I>(u, common_hencky_ptr));
    pip.push_back(new typename H::template OpCalculateLogC_dC<DIM, I>(
        u, common_hencky_ptr));
    pip.push_back(
        new typename H::template OpCalculateHenckyPlasticStress<DIM, I, 0>(
            u, common_hencky_ptr, m_D_ptr));
    pip.push_back(new typename H::template OpCalculatePiolaStress<DIM, I, 0>(
        u, common_hencky_ptr));
  } else {

    pip.push_back(new OpSymmetrizeTensor<SPACE_DIM>(
        common_plastic_ptr->mGradPtr, common_plastic_ptr->mStrainPtr));
    pip.push_back(new typename P::template OpPlasticStress<DIM, I>(
        common_plastic_ptr, m_D_ptr));
  }

  pip.push_back(new typename P::template OpCalculatePlasticSurface<DIM, I>(
      u, common_plastic_ptr));

  return std::make_tuple(common_plastic_ptr, common_hencky_ptr);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode
opFactoryDomainRhs(MoFEM::Interface &m_field, std::string block_name, Pip &pip,
                   std::string u, std::string ep, std::string tau) {
  MoFEMFunctionBegin;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template LinearForm<I>;
  using OpInternalForceCauchy =
      typename B::template OpGradTimesSymTensor<1, DIM, DIM>;
  using OpInternalForcePiola =
      typename B::template OpGradTimesTensor<1, DIM, DIM>;

  using P = PlasticityIntegrators<DomainEleOp>;

  auto [common_plastic_ptr, common_hencky_ptr] =
      createCommonPlasticOps<DIM, I, DomainEleOp>(m_field, block_name, pip, u,
                                                  ep, tau, scale, Sev::inform);

  auto m_D_ptr = common_plastic_ptr->mDPtr;

  pip.push_back(new OpCalculateTensor2SymmetricFieldValuesDot<DIM>(
      ep, common_plastic_ptr->getPlasticStrainDotPtr()));
  pip.push_back(new OpCalculateScalarFieldValuesDot(
      tau, common_plastic_ptr->getPlasticTauDotPtr()));
  // if rotation not 0 calculate material derivative
  // if (angular_velocity[2] != 0.0) {
  pip.push_back(
      new OpCalculateTensor2SymmetricFieldGradient<SPACE_DIM, SPACE_DIM>(
          ep, common_plastic_ptr->getPlasticGradStrainPtr()));
  pip.push_back(new OpCalculateScalarFieldGradient<SPACE_DIM>(
      tau, common_plastic_ptr->getPlasticGradTauPtr()));
  pip.push_back(
      new OpCalculatePlasticConvRotatingFrame(ep, common_plastic_ptr));
  //}
  pip.push_back(new typename P::template OpCalculatePlasticity<DIM, I>(
      u, common_plastic_ptr, m_D_ptr));

  // Calculate internal forces
  if (common_hencky_ptr) {
    pip.push_back(new OpInternalForcePiola(
        u, common_hencky_ptr->getMatFirstPiolaStress()));
  } else {
    pip.push_back(new OpInternalForceCauchy(u, common_plastic_ptr->mStressPtr));
  }

  pip.push_back(
      new
      typename P::template Assembly<A>::template OpCalculateConstraintsRhs<I>(
          tau, common_plastic_ptr, m_D_ptr));
  pip.push_back(
      new typename P::template Assembly<A>::template OpCalculatePlasticFlowRhs<
          DIM, I>(ep, common_plastic_ptr, m_D_ptr));

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode
opFactoryDomainLhs(MoFEM::Interface &m_field, std::string block_name, Pip &pip,
                   std::string u, std::string ep, std::string tau) {
  MoFEMFunctionBegin;

  using namespace HenckyOps;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template BiLinearForm<I>;
  using OpKPiola = typename B::template OpGradTensorGrad<1, DIM, DIM, 1>;
  using OpKCauchy = typename B::template OpGradSymTensorGrad<1, DIM, DIM, 0>;

  using P = PlasticityIntegrators<DomainEleOp>;

  auto [common_plastic_ptr, common_hencky_ptr] =
      createCommonPlasticOps<DIM, I, DomainEleOp>(m_field, block_name, pip, u,
                                                  ep, tau, scale, Sev::verbose);

  auto m_D_ptr = common_plastic_ptr->mDPtr;

  pip.push_back(new OpCalculateTensor2SymmetricFieldValuesDot<DIM>(
      ep, common_plastic_ptr->getPlasticStrainDotPtr()));
  pip.push_back(new OpCalculateScalarFieldValuesDot(
      tau, common_plastic_ptr->getPlasticTauDotPtr()));
  // if rotation not 0 calculate material derivative
  // if (angular_velocity[2] != 0.0) {
  pip.push_back(
      new OpCalculateTensor2SymmetricFieldGradient<SPACE_DIM, SPACE_DIM>(
          ep, common_plastic_ptr->getPlasticGradStrainPtr()));
  pip.push_back(new OpCalculateScalarFieldGradient<SPACE_DIM>(
      tau, common_plastic_ptr->getPlasticGradTauPtr()));
  pip.push_back(
      new OpCalculatePlasticConvRotatingFrame(ep, common_plastic_ptr));
  //}

  pip.push_back(new typename P::template OpCalculatePlasticity<DIM, I>(
      u, common_plastic_ptr, m_D_ptr));

  if (common_hencky_ptr) {
    using H = HenckyOps::HenckyIntegrators<DomainEleOp>;
    pip.push_back(new typename H::template OpHenckyTangent<DIM, I, 0>(
        u, common_hencky_ptr, m_D_ptr));
    pip.push_back(new OpKPiola(u, u, common_hencky_ptr->getMatTangent()));
    pip.push_back(
        new typename P::template Assembly<A>::
            template OpCalculatePlasticInternalForceLhs_LogStrain_dEP<DIM, I>(
                u, ep, common_plastic_ptr, common_hencky_ptr, m_D_ptr));
  } else {
    pip.push_back(new OpKCauchy(u, u, m_D_ptr));
    pip.push_back(new typename P::template Assembly<A>::
                      template OpCalculatePlasticInternalForceLhs_dEP<DIM, I>(
                          u, ep, common_plastic_ptr, m_D_ptr));
  }

  if (common_hencky_ptr) {
    pip.push_back(
        new typename P::template Assembly<A>::
            template OpCalculateConstraintsLhs_LogStrain_dU<DIM, I>(
                tau, u, common_plastic_ptr, common_hencky_ptr, m_D_ptr));
    pip.push_back(
        new typename P::template Assembly<A>::
            template OpCalculatePlasticFlowLhs_LogStrain_dU<DIM, I>(
                ep, u, common_plastic_ptr, common_hencky_ptr, m_D_ptr));
  } else {
    pip.push_back(
        new
        typename P::template Assembly<A>::template OpCalculateConstraintsLhs_dU<
            DIM, I>(tau, u, common_plastic_ptr, m_D_ptr));
    pip.push_back(
        new
        typename P::template Assembly<A>::template OpCalculatePlasticFlowLhs_dU<
            DIM, I>(ep, u, common_plastic_ptr, m_D_ptr));
  }

  pip.push_back(
      new
      typename P::template Assembly<A>::template OpCalculatePlasticFlowLhs_dEP<
          DIM, I>(ep, ep, common_plastic_ptr, m_D_ptr));
  pip.push_back(
      new
      typename P::template Assembly<A>::template OpCalculatePlasticFlowLhs_dTAU<
          DIM, I>(ep, tau, common_plastic_ptr, m_D_ptr));
  pip.push_back(
      new
      typename P::template Assembly<A>::template OpCalculateConstraintsLhs_dEP<
          DIM, I>(tau, ep, common_plastic_ptr, m_D_ptr));
  pip.push_back(
      new
      typename P::template Assembly<A>::template OpCalculateConstraintsLhs_dTAU<
          I>(tau, tau, common_plastic_ptr));

  MoFEMFunctionReturn(0);
}

template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
MoFEMErrorCode opFactoryDomainReactions(MoFEM::Interface &m_field,
                                        std::string block_name, Pip &pip,
                                        std::string u, std::string ep,
                                        std::string tau) {
  MoFEMFunctionBegin;

  using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
      A>::template LinearForm<I>;
  using OpInternalForceCauchy =
      typename B::template OpGradTimesSymTensor<1, DIM, DIM>;
  using OpInternalForcePiola =
      typename B::template OpGradTimesTensor<1, DIM, DIM>;

  auto [common_plastic_ptr, common_hencky_ptr] =
      createCommonPlasticOps<DIM, I, DomainEleOp>(m_field, block_name, pip, u,
                                                  ep, tau, 1., Sev::inform);

  // Calculate internal forces
  if (common_hencky_ptr) {
    pip.push_back(new OpInternalForcePiola(
        u, common_hencky_ptr->getMatFirstPiolaStress()));
  } else {
    pip.push_back(new OpInternalForceCauchy(u, common_plastic_ptr->mStressPtr));
  }

  MoFEMFunctionReturn(0);
}


OpCalculatePlasticConvRotatingFrame::OpCalculatePlasticConvRotatingFrame(
    const std::string field_name, boost::shared_ptr<CommonData> common_data_ptr)
    : DomainEleOp(field_name, field_name, DomainEleOp::OPROW),
      commonDataPtr(common_data_ptr) {}

MoFEMErrorCode OpCalculatePlasticConvRotatingFrame::doWork(int side,
                                                           EntityType type,
                                                           EntData &data) {
  MoFEMFunctionBegin;

  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'k', SPACE_DIM> k;
  FTensor::Index<'I', 3> I;
  FTensor::Index<'J', 3> J;
  FTensor::Index<'K', 3> K;

  const size_t nb_dofs = data.getIndices().size();
  if (nb_dofs) {
    const size_t nb_integration_pts = data.getN().size1();

    auto t_ep_dot = getFTensor2SymmetricFromMat<SPACE_DIM>(
        *(commonDataPtr->getPlasticStrainDotPtr()));
    auto t_tau_dot =
        getFTensor0FromVec(*(commonDataPtr->getPlasticTauDotPtr()));
    auto t_coords = getFTensor1CoordsAtGaussPts();
    auto t_grad_tau =
        getFTensor1FromMat<SPACE_DIM>(*(commonDataPtr->getPlasticGradTauPtr()));

    // OpCalculateTensor2SymmetricFieldGradient
    auto t_grad_plastic_strain = getFTensor3DgFromMat<SPACE_DIM, SPACE_DIM>(
        *(commonDataPtr->getPlasticGradStrainPtr()));
    commonDataPtr->guidingVelocityPtr->resize(SPACE_DIM, nb_integration_pts,
                                              false);
    commonDataPtr->guidingVelocityPtr->clear();
    auto t_omega =
        getFTensor1FromMat<SPACE_DIM>(*commonDataPtr->guidingVelocityPtr);

    FTensor::Tensor1<double, 3> t1_omega(
        angular_velocity[0], angular_velocity[1], angular_velocity[2]);
    FTensor::Tensor2<double, 3, 3> Omega;
    Omega(I, K) = FTensor::levi_civita<double>(I, J, K) * t1_omega(J);
    auto t_tau = getFTensor0FromVec(commonDataPtr->plasticTau);

    for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
      // FIXME: maybe here we should take cooordinates that are rotated, X
      // but rotation will not change the velocity, right?
      t_omega(i) = Omega(i, j) * t_coords(j);

      t_tau_dot += t_grad_tau(i) * t_omega(i);
      t_ep_dot(i, j) += t_grad_plastic_strain(i, j, k) * t_omega(k);

      ++t_coords;
      ++t_omega;
      ++t_grad_plastic_strain;
      ++t_ep_dot;
      ++t_tau_dot;
      ++t_tau;
      ++t_grad_tau;
    }
  }

  MoFEMFunctionReturn(0);
}

OpRhsSkeleton::OpRhsSkeleton(boost::shared_ptr<SideData> side_data_ptr,
                             boost::shared_ptr<SideEle> side_fe_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      sideDataPtr(side_data_ptr), sideFEPtr(side_fe_ptr) {}

OpLhsSkeleton::OpLhsSkeleton(boost::shared_ptr<SideData> side_data_ptr,
                             boost::shared_ptr<SideEle> side_fe_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      sideDataPtr(side_data_ptr), sideFEPtr(side_fe_ptr) {}

MoFEMErrorCode OpRhsSkeleton::doWork(int side, EntityType type,
                                     EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'I', 3> I;
  FTensor::Index<'J', 3> J;
  FTensor::Index<'K', 3> K;

  const auto in_the_loop =
      sideFEPtr->nInTheLoop; // return number of elements on the side

  auto not_side = [](auto s) {
    return s == LEFT_SIDE ? RIGHT_SIDE : LEFT_SIDE;
  };

  auto get_ntensor = [](auto &base_mat) {
    return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
        &*base_mat.data().begin());
  };

  if (in_the_loop > 0) {

    // get normal of the face or edge
    auto t_normal = getFTensor1Normal();
    const auto nb_gauss_pts = getGaussPts().size2();

    for (auto s0 : {LEFT_SIDE, RIGHT_SIDE}) {

      // gent number of DOFs on the right side.
      const auto nb_rows = sideDataPtr->indicesRowSideMapTau[s0].size();

      if (nb_rows) {

        resSkelton.resize(nb_rows, false);
        resSkelton.clear();

        // get orientation of the local element edge
        const auto opposite_s0 = not_side(s0);
        const auto sense_row = sideDataPtr->senseMap[s0];
#ifndef NDEBUG
        const auto opposite_sense_row = sideDataPtr->senseMap[opposite_s0];
        if (sense_row * opposite_sense_row > 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Should be opposite sign");
#endif

        // iterate the side cols
        const auto nb_row_base_functions =
            sideDataPtr->rowBaseSideMapTau[s0].size2();

        auto t_w = getFTensor0IntegrationWeight();
        auto arr_t_l =
            make_array(getFTensor0FromVec(sideDataPtr->tauVec[LEFT_SIDE]),
                       getFTensor0FromVec(sideDataPtr->tauVec[RIGHT_SIDE]));
        // auto arr_t_vel = make_array(
        //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[LEFT_SIDE]),
        //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[RIGHT_SIDE]));

        auto t_coords = getFTensor1CoordsAtGaussPts();

        FTensor::Tensor1<double, 3> t1_omega(
            angular_velocity[0], angular_velocity[1], angular_velocity[2]);
        FTensor::Tensor2<double, 3, 3> Omega;
        Omega(I, K) = FTensor::levi_civita<double>(I, J, K) * t1_omega(J);

        auto next = [&]() {
          for (auto &t_l : arr_t_l)
            ++t_l;
          // for (auto &t_vel : arr_t_vel)
          //   ++t_vel;
          ++t_coords;
        };

#ifndef NDEBUG
        if (nb_gauss_pts != sideDataPtr->rowBaseSideMapTau[s0].size1())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Inconsistent number of DOFs");
#endif

        auto t_row_base = get_ntensor(sideDataPtr->rowBaseSideMapTau[s0]);
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          FTensor::Tensor1<double, 3> t_vel;
          t_vel(i) = Omega(i, j) * t_coords(j);
          // FTensor::Tensor1<double, SPACE_DIM> t_vel;
          // t_vel(i) = (arr_t_vel[LEFT_SIDE](i) + arr_t_vel[RIGHT_SIDE](i))
          // / 2.;
          const auto dot = sense_row * (t_normal(i) * t_vel(i));
          const auto l_upwind_side = (dot > 0) ? s0 : opposite_s0;
          const auto l_upwind = arr_t_l[l_upwind_side];
          const auto res = t_w * dot * l_upwind;
          next();
          ++t_w;

          auto rr = 0;
          for (; rr != nb_rows; ++rr) {
            resSkelton[rr] += t_row_base * res;
            ++t_row_base;
          }
          for (; rr < nb_row_base_functions; ++rr) {
            ++t_row_base;
          }
        }
        // assemble local operator vector to global vector
        CHKERR ::VecSetValues(getTSf(),
                              sideDataPtr->indicesRowSideMapTau[s0].size(),
                              &*sideDataPtr->indicesRowSideMapTau[s0].begin(),
                              &*resSkelton.begin(), ADD_VALUES);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpLhsSkeleton::doWork(int side, EntityType type,
                                     EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'I', 3> I;
  FTensor::Index<'J', 3> J;
  FTensor::Index<'K', 3> K;

  const auto in_the_loop =
      sideFEPtr->nInTheLoop; // return number of elements on the side

  auto not_side = [](auto s) {
    return s == LEFT_SIDE ? RIGHT_SIDE : LEFT_SIDE;
  };

  auto get_ntensor = [](auto &base_mat) {
    return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
        &*base_mat.data().begin());
  };

  if (in_the_loop > 0) {

    // get normal of the face or edge
    auto t_normal = getFTensor1Normal();
    const auto nb_gauss_pts = getGaussPts().size2();

    for (auto s0 : {LEFT_SIDE, RIGHT_SIDE}) {

      // gent number of DOFs on the right side.
      const auto nb_rows = sideDataPtr->indicesRowSideMapTau[s0].size();

      if (nb_rows) {

        // get orientation of the local element edge
        const auto opposite_s0 = not_side(s0);
        const auto sense_row = sideDataPtr->senseMap[s0];

        // iterate the side cols
        const auto nb_row_base_functions =
            sideDataPtr->rowBaseSideMapTau[s0].size2();

        for (auto s1 : {LEFT_SIDE, RIGHT_SIDE}) {

          // gent number of DOFs on the right side.
          const auto nb_cols = sideDataPtr->indicesColSideMapTau[s1].size();
          const auto sense_col = sideDataPtr->senseMap[s1];

          // resize local element matrix
          matSkeleton.resize(nb_rows, nb_cols, false);
          matSkeleton.clear();

          auto t_w = getFTensor0IntegrationWeight();

          auto t_coords = getFTensor1CoordsAtGaussPts();

          FTensor::Tensor1<double, 3> t1_omega(
              angular_velocity[0], angular_velocity[1], angular_velocity[2]);
          FTensor::Tensor2<double, 3, 3> Omega;
          Omega(I, K) = FTensor::levi_civita<double>(I, J, K) * t1_omega(J);

          // auto arr_t_vel = make_array(
          //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[LEFT_SIDE]),
          //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[RIGHT_SIDE]));

          auto next = [&]() {
            //  for (auto &t_vel : arr_t_vel)
            //    ++t_vel;
            ++t_coords;
          };

          auto t_row_base = get_ntensor(sideDataPtr->rowBaseSideMapTau[s0]);
          for (int gg = 0; gg != nb_gauss_pts; ++gg) {
            FTensor::Tensor1<double, 3> t_vel;
            t_vel(i) = Omega(i, j) * t_coords(j);

            // FTensor::Tensor1<double, SPACE_DIM> t_vel;
            // t_vel(i) =
            //(arr_t_vel[LEFT_SIDE](i) + arr_t_vel[RIGHT_SIDE](i)) / 2.;
            const auto dot = sense_row * (t_normal(i) * t_vel(i));
            const auto l_upwind_side = (dot > 0) ? s0 : opposite_s0;
            const auto sense_upwind = sideDataPtr->senseMap[l_upwind_side];
            auto res = t_w * dot; // * sense_row * sense_upwind;
            next();
            ++t_w;
            auto rr = 0;
            if (s1 == l_upwind_side) {
              for (; rr != nb_rows; ++rr) {
                auto get_ntensor = [](auto &base_mat, auto gg, auto bb) {
                  double *ptr = &base_mat(gg, bb);
                  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(ptr);
                };
                auto t_col_base =
                    get_ntensor(sideDataPtr->colBaseSideMapTau[s1], gg, 0);
                const auto res_row = res * t_row_base;
                ++t_row_base;
                // iterate columns
                for (size_t cc = 0; cc != nb_cols; ++cc) {
                  matSkeleton(rr, cc) += res_row * t_col_base;
                  ++t_col_base;
                }
              }
            }
            for (; rr < nb_row_base_functions; ++rr) {
              ++t_row_base;
            }
          }
          // assemble system
          CHKERR MatSetValues(getTSB(),
                              sideDataPtr->indicesRowSideMapTau[s0].size(),
                              &*sideDataPtr->indicesRowSideMapTau[s0].begin(),
                              sideDataPtr->indicesColSideMapTau[s1].size(),
                              &*sideDataPtr->indicesColSideMapTau[s1].begin(),
                              &*matSkeleton.data().begin(), ADD_VALUES);
        }
      }
    }
  }
  MoFEMFunctionReturn(0);
}

OpRhsSkeletonEp::OpRhsSkeletonEp(boost::shared_ptr<SideData> side_data_ptr,
                                 boost::shared_ptr<SideEle> side_fe_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      sideDataPtr(side_data_ptr), sideFEPtr(side_fe_ptr) {}

OpLhsSkeletonEp::OpLhsSkeletonEp(boost::shared_ptr<SideData> side_data_ptr,
                                 boost::shared_ptr<SideEle> side_fe_ptr)
    : BoundaryEleOp(NOSPACE, BoundaryEleOp::OPSPACE),
      sideDataPtr(side_data_ptr), sideFEPtr(side_fe_ptr) {}

MoFEMErrorCode OpRhsSkeletonEp::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'I', 3> I;
  FTensor::Index<'J', 3> J;
  FTensor::Index<'K', 3> K;
  FTensor::Index<'L', SPACE_DIM> L;

  const auto in_the_loop =
      sideFEPtr->nInTheLoop; // return number of elements on the side

  auto not_side = [](auto s) {
    return s == LEFT_SIDE ? RIGHT_SIDE : LEFT_SIDE;
  };

  auto get_ntensor = [](auto &base_mat) {
    return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
        &*base_mat.data().begin());
  };

  if (in_the_loop > 0) {

    // get normal of the face or edge
    auto t_normal = getFTensor1Normal();
    const auto nb_gauss_pts = getGaussPts().size2();

    for (auto s0 : {LEFT_SIDE, RIGHT_SIDE}) {

      // gent number of DOFs on the right side.
      const auto nb_rows = sideDataPtr->indicesRowSideMapEp[s0].size();

      if (nb_rows) {

        resSkelton.resize(nb_rows, false);
        resSkelton.clear();

        // get orientation of the local element edge
        const auto opposite_s0 = not_side(s0);
        const auto sense_row = sideDataPtr->senseMap[s0];
#ifndef NDEBUG
        const auto opposite_sense_row = sideDataPtr->senseMap[opposite_s0];
        if (sense_row * opposite_sense_row > 0)
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Should be opposite sign");
#endif

        // iterate the side cols
        const auto nb_row_base_functions =
            sideDataPtr->rowBaseSideMapEp[s0].size2();

        auto t_w = getFTensor0IntegrationWeight();
        auto arr_t_l = make_array(getFTensor2SymmetricFromMat<SPACE_DIM>(
                                      sideDataPtr->epMat[LEFT_SIDE]),
                                  getFTensor2SymmetricFromMat<SPACE_DIM>(
                                      sideDataPtr->epMat[RIGHT_SIDE]));
        // auto arr_t_vel = make_array(
        //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[LEFT_SIDE]),
        //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[RIGHT_SIDE]));

        auto t_coords = getFTensor1CoordsAtGaussPts();

        FTensor::Tensor1<double, 3> t1_omega(
            angular_velocity[0], angular_velocity[1], angular_velocity[2]);
        FTensor::Tensor2<double, 3, 3> Omega;
        Omega(I, K) = FTensor::levi_civita<double>(I, J, K) * t1_omega(J);

        auto t_L = symm_L_tensor(FTensor::Number<SPACE_DIM>());

        auto next = [&]() {
          for (auto &t_l : arr_t_l)
            ++t_l;
          // for (auto &t_vel : arr_t_vel)
          //   ++t_vel;
          ++t_coords;
        };

#ifndef NDEBUG
        if (nb_gauss_pts != sideDataPtr->rowBaseSideMapEp[s0].size1())
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Inconsistent number of DOFs");
#endif

        auto t_row_base = get_ntensor(sideDataPtr->rowBaseSideMapEp[s0]);
        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          FTensor::Tensor1<double, 3> t_vel;
          t_vel(i) = Omega(i, j) * t_coords(j);
          // FTensor::Tensor1<double, SPACE_DIM> t_vel;
          // t_vel(i) = (arr_t_vel[LEFT_SIDE](i) + arr_t_vel[RIGHT_SIDE](i))
          // / 2.;
          const auto dot = sense_row * (t_normal(i) * t_vel(i));
          const auto l_upwind_side = (dot > 0) ? s0 : opposite_s0;
          const auto l_upwind = arr_t_l[l_upwind_side];
          FTensor::Tensor2_symmetric<double, SPACE_DIM> t_res;
          t_res(i, j) = t_w * dot * l_upwind(i, j);
          next();
          ++t_w;
          FTensor::Tensor1<double, size_symm> t_rhs;
          t_rhs(L) = t_res(i, j) * t_L(i, j, L);

          auto t_nf = getFTensor1FromArray<size_symm, size_symm>(resSkelton);

          auto rr = 0;
          for (; rr != nb_rows / size_symm; ++rr) {
            t_nf(L) += t_row_base * t_rhs(L);
            ++t_row_base;
            ++t_nf;
          }
          for (; rr < nb_row_base_functions; ++rr) {
            ++t_row_base;
          }
        }
        // assemble local operator vector to global vector
        CHKERR ::VecSetValues(getTSf(),
                              sideDataPtr->indicesRowSideMapEp[s0].size(),
                              &*sideDataPtr->indicesRowSideMapEp[s0].begin(),
                              &*resSkelton.begin(), ADD_VALUES);
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpLhsSkeletonEp::doWork(int side, EntityType type,
                                       EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', SPACE_DIM> i;
  FTensor::Index<'j', SPACE_DIM> j;
  FTensor::Index<'I', 3> I;
  FTensor::Index<'J', 3> J;
  FTensor::Index<'K', 3> K;

  const auto in_the_loop =
      sideFEPtr->nInTheLoop; // return number of elements on the side

  auto not_side = [](auto s) {
    return s == LEFT_SIDE ? RIGHT_SIDE : LEFT_SIDE;
  };

  auto get_ntensor = [](auto &base_mat) {
    return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
        &*base_mat.data().begin());
  };

  if (in_the_loop > 0) {

    // get normal of the face or edge
    auto t_normal = getFTensor1Normal();
    const auto nb_gauss_pts = getGaussPts().size2();

    for (auto s0 : {LEFT_SIDE, RIGHT_SIDE}) {

      // gent number of DOFs on the right side.
      const auto nb_rows = sideDataPtr->indicesRowSideMapEp[s0].size();

      if (nb_rows) {

        // get orientation of the local element edge
        const auto opposite_s0 = not_side(s0);
        const auto sense_row = sideDataPtr->senseMap[s0];

        // iterate the side cols
        const auto nb_row_base_functions =
            sideDataPtr->rowBaseSideMapEp[s0].size2();

        for (auto s1 : {LEFT_SIDE, RIGHT_SIDE}) {

          // gent number of DOFs on the right side.
          const auto nb_cols = sideDataPtr->indicesColSideMapEp[s1].size();
          const auto sense_col = sideDataPtr->senseMap[s1];

          // resize local element matrix
          matSkeleton.resize(nb_rows, nb_cols, false);
          matSkeleton.clear();

          auto t_w = getFTensor0IntegrationWeight();

          auto t_coords = getFTensor1CoordsAtGaussPts();

          FTensor::Tensor1<double, 3> t1_omega(
              angular_velocity[0], angular_velocity[1], angular_velocity[2]);
          FTensor::Tensor2<double, 3, 3> Omega;
          Omega(I, K) = FTensor::levi_civita<double>(I, J, K) * t1_omega(J);

          // auto arr_t_vel = make_array(
          //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[LEFT_SIDE]),
          //     getFTensor1FromMat<SPACE_DIM>(sideDataPtr->velMat[RIGHT_SIDE]));

          auto next = [&]() {
            //  for (auto &t_vel : arr_t_vel)
            //    ++t_vel;
            ++t_coords;
          };

          auto t_row_base = get_ntensor(sideDataPtr->rowBaseSideMapEp[s0]);
          for (int gg = 0; gg != nb_gauss_pts; ++gg) {
            FTensor::Tensor1<double, 3> t_vel;
            t_vel(i) = Omega(i, j) * t_coords(j);

            // FTensor::Tensor1<double, SPACE_DIM> t_vel;
            // t_vel(i) =
            //(arr_t_vel[LEFT_SIDE](i) + arr_t_vel[RIGHT_SIDE](i)) / 2.;
            const auto dot = sense_row * (t_normal(i) * t_vel(i));
            const auto l_upwind_side = (dot > 0) ? s0 : opposite_s0;
            const auto sense_upwind = sideDataPtr->senseMap[l_upwind_side];
            auto res = t_w * dot; // * sense_row * sense_upwind;
            next();
            ++t_w;
            auto rr = 0;
            if (s1 == l_upwind_side) {
              for (; rr != nb_rows / size_symm; ++rr) {
                auto get_ntensor = [](auto &base_mat, auto gg, auto bb) {
                  double *ptr = &base_mat(gg, bb);
                  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(ptr);
                };
                auto t_col_base =
                    get_ntensor(sideDataPtr->colBaseSideMapEp[s1], gg, 0);
                const auto res_row = res * t_row_base;
                ++t_row_base;
                // iterate columns
                for (size_t cc = 0; cc != nb_cols / size_symm; ++cc) {
                  matSkeleton(rr, cc) += res_row * t_col_base;
                  ++t_col_base;
                }
              }
            }
            for (; rr < nb_row_base_functions; ++rr) {
              ++t_row_base;
            }
          }
          // assemble system
          CHKERR ::MatSetValues(getTSB(),
                                sideDataPtr->indicesRowSideMapEp[s0].size(),
                                &*sideDataPtr->indicesRowSideMapEp[s0].begin(),
                                sideDataPtr->indicesColSideMapEp[s1].size(),
                                &*sideDataPtr->indicesColSideMapEp[s1].begin(),
                                &*matSkeleton.data().begin(), ADD_VALUES);
        }
      }
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode getSideFE(MoFEM::Interface &m_field,
                         boost::shared_ptr<SideData> side_data_ptr,
                         boost::shared_ptr<SideEle> side_fe_ptr) {
  MoFEMFunctionBegin;

  auto simple = m_field.getInterface<Simple>();

  auto tau_ptr = boost::make_shared<VectorDouble>();
  auto ep_ptr = boost::make_shared<MatrixDouble>();

  struct OpSideData : public SideEleOp {
    OpSideData(boost::shared_ptr<SideData> side_data_ptr)
        : SideEleOp("TAU", "TAU", SideEleOp::OPROWCOL),
          sideDataPtr(side_data_ptr) {
      std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
      for (auto t = moab::CN::TypeDimensionMap[SPACE_DIM].first;
           t <= moab::CN::TypeDimensionMap[SPACE_DIM].second; ++t)
        doEntities[t] = true;
      sYmm = false;
    }

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type, EntData &row_data,
                          EntData &col_data) {
      MoFEMFunctionBegin;
      if ((CN::Dimension(row_type) == SPACE_DIM) &&
          (CN::Dimension(col_type) == SPACE_DIM)) {

        auto reset = [&](auto nb_in_loop) {
          sideDataPtr->feSideHandle[nb_in_loop] = 0;
          sideDataPtr->indicesRowSideMapTau[nb_in_loop].clear();
          sideDataPtr->indicesColSideMapTau[nb_in_loop].clear();
          sideDataPtr->rowBaseSideMapTau[nb_in_loop].clear();
          sideDataPtr->colBaseSideMapTau[nb_in_loop].clear();
          sideDataPtr->senseMap[nb_in_loop] = 0;
        };

        const auto nb_in_loop = getFEMethod()->nInTheLoop;
        if (nb_in_loop == 0)
          for (auto s : {0, 1})
            reset(s);

        sideDataPtr->currentFESide = nb_in_loop;
        sideDataPtr->senseMap[nb_in_loop] = getSkeletonSense();

      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Should not happen");
      }

      MoFEMFunctionReturn(0);
    };

  private:
    boost::shared_ptr<SideData> sideDataPtr;
  };

  struct OpCalculateVolumeOnSideTau : public SideEleOp {

    OpCalculateVolumeOnSideTau(boost::shared_ptr<SideData> side_data_ptr,
                               boost::shared_ptr<VectorDouble> tau_ptr)
        : SideEleOp("TAU", "TAU", SideEleOp::OPROWCOL),
          sideDataPtr(side_data_ptr), tauPtr(tau_ptr) {
      std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
      for (auto t = moab::CN::TypeDimensionMap[SPACE_DIM].first;
           t <= moab::CN::TypeDimensionMap[SPACE_DIM].second; ++t)
        doEntities[t] = true;
      sYmm = false;
    }

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type, EntData &row_data,
                          EntData &col_data) {
      MoFEMFunctionBegin;

      if ((CN::Dimension(row_type) == SPACE_DIM) &&
          (CN::Dimension(col_type) == SPACE_DIM)) {

        // t_tau = getFTensor0FromVec(*tauPtr);

        const auto nb_in_loop = sideDataPtr->currentFESide;
        sideDataPtr->feSideHandle[nb_in_loop] = getFEEntityHandle();
        sideDataPtr->indicesRowSideMapTau[nb_in_loop] = row_data.getIndices();
        sideDataPtr->indicesColSideMapTau[nb_in_loop] = col_data.getIndices();
        sideDataPtr->rowBaseSideMapTau[nb_in_loop] = row_data.getN();
        sideDataPtr->colBaseSideMapTau[nb_in_loop] = col_data.getN();
        (sideDataPtr->tauVec)[nb_in_loop] = *tauPtr;
        //(sideDataPtr->velMat)[nb_in_loop] = *velPtr;

        // #ifndef NDEBUG
        //         if ((sideDataPtr->tauVec)[nb_in_loop].size() !=
        //             (sideDataPtr->velMat)[nb_in_loop].size2())
        //           SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        //                    "Wrong number of integaration pts %d != %d",
        //                    (sideDataPtr->tauVec)[nb_in_loop].size(),
        //                    (sideDataPtr->velMat)[nb_in_loop].size2());
        //         if ((sideDataPtr->velMat)[nb_in_loop].size1() != SPACE_DIM)
        //           SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        //                    "Wrong size of velocity vector size = %d",
        //                    (sideDataPtr->velMat)[nb_in_loop].size1());
        // #endif

        if (!nb_in_loop) {
          (sideDataPtr->tauVec)[1] = sideDataPtr->tauVec[0];
          //(sideDataPtr->velMat)[1] = (sideDataPtr->velMat)[0];
        } else {
          // #ifndef NDEBUG
          //           if (sideDataPtr->rowBaseSideMapTau[0].size1() !=
          //               sideDataPtr->rowBaseSideMapTau[1].size1()) {
          //             SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          //                      "Wrong number of integration pt %d != %d",
          //                      sideDataPtr->rowBaseSideMapTau[0].size1(),
          //                      sideDataPtr->rowBaseSideMapTau[1].size1());
          //           }
          //           if (sideDataPtr->colBaseSideMapTau[0].size1() !=
          //               sideDataPtr->colBaseSideMapTau[1].size1()) {
          //             SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          //                     "Wrong number of integration pt");
          //           }
          // #endif
        }

      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Should not happen");
      }

      MoFEMFunctionReturn(0);
    };

  private:
    boost::shared_ptr<SideData> sideDataPtr;
    boost::shared_ptr<VectorDouble> tauPtr;
  };

  struct OpCalculateVolumeOnSideEp : public SideEleOp {

    OpCalculateVolumeOnSideEp(boost::shared_ptr<SideData> side_data_ptr,
                              boost::shared_ptr<MatrixDouble> ep_ptr)
        : SideEleOp("EP", "EP", SideEleOp::OPROWCOL),
          sideDataPtr(side_data_ptr), epPtr(ep_ptr) {
      std::fill(&doEntities[MBVERTEX], &doEntities[MBMAXTYPE], false);
      for (auto t = moab::CN::TypeDimensionMap[SPACE_DIM].first;
           t <= moab::CN::TypeDimensionMap[SPACE_DIM].second; ++t)
        doEntities[t] = true;
      sYmm = false;
    }

    MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                          EntityType col_type, EntData &row_data,
                          EntData &col_data) {
      MoFEMFunctionBegin;

      if ((CN::Dimension(row_type) == SPACE_DIM) &&
          (CN::Dimension(col_type) == SPACE_DIM)) {

        // t_tau = getFTensor0FromVec(*tauPtr);

        const auto nb_in_loop = sideDataPtr->currentFESide;
        sideDataPtr->feSideHandle[nb_in_loop] = getFEEntityHandle();
        sideDataPtr->indicesRowSideMapEp[nb_in_loop] = row_data.getIndices();
        sideDataPtr->indicesColSideMapEp[nb_in_loop] = col_data.getIndices();
        sideDataPtr->rowBaseSideMapEp[nb_in_loop] = row_data.getN();
        sideDataPtr->colBaseSideMapEp[nb_in_loop] = col_data.getN();
        (sideDataPtr->epMat)[nb_in_loop] = *epPtr;
        //(sideDataPtr->velMat)[nb_in_loop] = *velPtr;

        // #ifndef NDEBUG
        //         if ((sideDataPtr->tauVec)[nb_in_loop].size() !=
        //             (sideDataPtr->velMat)[nb_in_loop].size2())
        //           SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        //                    "Wrong number of integaration pts %d != %d",
        //                    (sideDataPtr->tauVec)[nb_in_loop].size(),
        //                    (sideDataPtr->velMat)[nb_in_loop].size2());
        //         if ((sideDataPtr->velMat)[nb_in_loop].size1() != SPACE_DIM)
        //           SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        //                    "Wrong size of velocity vector size = %d",
        //                    (sideDataPtr->velMat)[nb_in_loop].size1());
        // #endif

        if (!nb_in_loop) {
          (sideDataPtr->epMat)[1] = sideDataPtr->epMat[0];
          //(sideDataPtr->velMat)[1] = (sideDataPtr->velMat)[0];
        } else {
          // #ifndef NDEBUG
          //           if (sideDataPtr->rowBaseSideMapTau[0].size1() !=
          //               sideDataPtr->rowBaseSideMapTau[1].size1()) {
          //             SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          //                      "Wrong number of integration pt %d != %d",
          //                      sideDataPtr->rowBaseSideMapTau[0].size1(),
          //                      sideDataPtr->rowBaseSideMapTau[1].size1());
          //           }
          //           if (sideDataPtr->colBaseSideMapTau[0].size1() !=
          //               sideDataPtr->colBaseSideMapTau[1].size1()) {
          //             SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
          //                     "Wrong number of integration pt");
          //           }
          // #endif
        }

      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Should not happen");
      }

      MoFEMFunctionReturn(0);
    };

  private:
    boost::shared_ptr<SideData> sideDataPtr;
    boost::shared_ptr<MatrixDouble> epPtr;
  };

  // Create aliased shared pointers, all elements are destroyed if side_fe_ptr
  // is destroyed
  CHKERR AddHOOps<SPACE_DIM, SPACE_DIM, SPACE_DIM>::add(
      side_fe_ptr->getOpPtrVector(), {H1, L2, HDIV}, "GEOMETRY");
  side_fe_ptr->getOpPtrVector().push_back(new OpSideData(side_data_ptr));
  side_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateScalarFieldValues("TAU", tau_ptr));
  side_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateTensor2SymmetricFieldValues<SPACE_DIM>("EP", ep_ptr));

  side_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateVolumeOnSideTau(side_data_ptr, tau_ptr));
  side_fe_ptr->getOpPtrVector().push_back(
      new OpCalculateVolumeOnSideEp(side_data_ptr, ep_ptr));

  MoFEMFunctionReturn(0);
};

} // namespace PlasticOps
