/**
 * @file HookeOps.hpp
 * @author MD Tanzib Ehsan Sanglap (m.sanglap.1@research.gla.ac.uk)
 * @brief Implementation of Hooke's law for linear elasticity
 * @version 0.14.0
 * @date 2025-04-02
 */

#ifndef __HOOKE_OPS_HPP__
#define __HOOKE_OPS_HPP__

namespace HookeOps {


 
 
 
 struct CommonData : public boost::enable_shared_from_this<CommonData> {
   boost::shared_ptr<MatrixDouble> matGradPtr;
   boost::shared_ptr<MatrixDouble> matDPtr;
   boost::shared_ptr<MatrixDouble> matLogCPlastic;
 
   MatrixDouble matCauchyStress;
   MatrixDouble matTangent;

 
   inline auto getMatCauchyStress() {
     return boost::shared_ptr<MatrixDouble>(shared_from_this(),
                                            &matCauchyStress);
   }

 
   inline auto getMatTangent() {
     return boost::shared_ptr<MatrixDouble>(shared_from_this(), &matTangent);
   }
 };


 
 template <int DIM, IntegrationType I, typename DomainEleOp, int S>
 struct OpCalculateCauchyStressImpl;
 

 

 
 template <int DIM, IntegrationType I, typename DomainEleOp, int S>
 struct OpHookeTangentImpl;
 
 
 
 
 
 
 template <int DIM, typename DomainEleOp, int S>
 struct OpCalculateCauchyStressImpl<DIM, GAUSS, DomainEleOp, S>
     : public DomainEleOp {
 
        OpCalculateCauchyStressImpl(const std::string field_name,
                               boost::shared_ptr<CommonData> common_data)
       : DomainEleOp(field_name, DomainEleOp::OPROW),
         commonDataPtr(common_data) {
     std::fill(&DomainEleOp::doEntities[MBEDGE],
               &DomainEleOp::doEntities[MBMAXTYPE], false);
   }
 
   MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
     MoFEMFunctionBegin;
 
     FTensor::Index<'i', DIM> i;
     FTensor::Index<'j', DIM> j;
     FTensor::Index<'k', DIM> k;
     FTensor::Index<'l', DIM> l;
 
     constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
 
     // const size_t nb_gauss_pts = matGradPtr->size2();
     const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
     auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*commonDataPtr->matDPtr);
     auto t_logC = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matLogC);
     constexpr auto size_symm = (DIM * (DIM + 1)) / 2;
     commonDataPtr->matHenckyStress.resize(size_symm, nb_gauss_pts, false);
     auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
 
     for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
       t_T(i, j) = t_D(i, j, k, l) * t_logC(k, l);
       ++t_logC;
       ++t_T;
       ++t_D;
     }
 
     MoFEMFunctionReturn(0);
   }
 
 private:
   boost::shared_ptr<CommonData> commonDataPtr;
 };
 
 
 
 
 
 template <int DIM, typename DomainEleOp, int S>
 struct OpHookeTangentImpl<DIM, GAUSS, DomainEleOp, S> : public DomainEleOp {
    OpHookeTangentImpl(const std::string field_name,
                       boost::shared_ptr<CommonData> common_data,
                       boost::shared_ptr<MatrixDouble> mat_D_ptr = nullptr)
       : DomainEleOp(field_name, DomainEleOp::OPROW),
         commonDataPtr(common_data) {
     std::fill(&DomainEleOp::doEntities[MBEDGE],
               &DomainEleOp::doEntities[MBMAXTYPE], false);
     if (mat_D_ptr)
       matDPtr = mat_D_ptr;
     else
       matDPtr = commonDataPtr->matDPtr;
   }
 
   MoFEMErrorCode doWork(int side, EntityType type, EntData &data) {
     MoFEMFunctionBegin;
 
     FTensor::Index<'i', DIM> i;
     FTensor::Index<'j', DIM> j;
     FTensor::Index<'k', DIM> k;
     FTensor::Index<'l', DIM> l;
     FTensor::Index<'m', DIM> m;
     FTensor::Index<'n', DIM> n;
     FTensor::Index<'o', DIM> o;
     FTensor::Index<'p', DIM> p;
 
     constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
     // const size_t nb_gauss_pts = matGradPtr->size2();
     const size_t nb_gauss_pts = DomainEleOp::getGaussPts().size2();
     commonDataPtr->matTangent.resize(DIM * DIM * DIM * DIM, nb_gauss_pts);
     auto dP_dF =
         getFTensor4FromMat<DIM, DIM, DIM, DIM, 1>(commonDataPtr->matTangent);
 
     auto t_D = getFTensor4DdgFromMat<DIM, DIM, S>(*matDPtr);
     auto t_eig_val = getFTensor1FromMat<DIM>(commonDataPtr->matEigVal);
     auto t_eig_vec = getFTensor2FromMat<DIM, DIM>(commonDataPtr->matEigVec);
     auto t_T = getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matHenckyStress);
     auto t_S =
         getFTensor2SymmetricFromMat<DIM>(commonDataPtr->matSecondPiolaStress);
     auto t_grad = getFTensor2FromMat<DIM, DIM>(*(commonDataPtr->matGradPtr));
     auto t_logC_dC = getFTensor4DdgFromMat<DIM, DIM>(commonDataPtr->matLogCdC);
 
     for (size_t gg = 0; gg != nb_gauss_pts; ++gg) {
 
 #ifdef HENCKY_SMALL_STRAIN
       dP_dF(i, j, k, l) = t_D(i, j, k, l);
 #else
 
       FTensor::Tensor2<double, DIM, DIM> t_F;
       t_F(i, j) = t_grad(i, j) + t_kd(i, j);
 
       FTensor::Tensor1<double, DIM> eig;
       FTensor::Tensor2<double, DIM, DIM> eigen_vec;
       FTensor::Tensor2_symmetric<double, DIM> T;
       eig(i) = t_eig_val(i);
       eigen_vec(i, j) = t_eig_vec(i, j);
       T(i, j) = t_T(i, j);
 
       // rare case when two eigen values are equal
       auto nb_uniq = get_uniq_nb<DIM>(&eig(0));
 
       FTensor::Tensor4<double, DIM, DIM, DIM, DIM> dC_dF;
       dC_dF(i, j, k, l) = (t_kd(i, l) * t_F(k, j)) + (t_kd(j, l) * t_F(k, i));
 
       auto TL =
           EigenMatrix::getDiffDiffMat(eig, eigen_vec, f, d_f, dd_f, T, nb_uniq);
 
       TL(i, j, k, l) *= 4;
       FTensor::Ddg<double, DIM, DIM> P_D_P_plus_TL;
       P_D_P_plus_TL(i, j, k, l) =
           TL(i, j, k, l) +
           (t_logC_dC(i, j, o, p) * t_D(o, p, m, n)) * t_logC_dC(m, n, k, l);
       P_D_P_plus_TL(i, j, k, l) *= 0.5;
       dP_dF(i, j, m, n) = t_kd(i, m) * (t_kd(k, n) * t_S(k, j));
       dP_dF(i, j, m, n) +=
           t_F(i, k) * (P_D_P_plus_TL(k, j, o, p) * dC_dF(o, p, m, n));
 
 #endif
 
       ++dP_dF;
 
       ++t_grad;
       ++t_eig_val;
       ++t_eig_vec;
       ++t_logC_dC;
       ++t_S;
       ++t_T;
       ++t_D;
     }
 
     MoFEMFunctionReturn(0);
   }
 
 private:
   boost::shared_ptr<CommonData> commonDataPtr;
   boost::shared_ptr<MatrixDouble> matDPtr;
 };
 
 template <typename DomainEleOp> struct HookeIntegrators {

 

   template <int DIM, IntegrationType I, int S>
   using OpCalculateCauchyStress =
       OpCalculateCauchyStressImpl<DIM, I, DomainEleOp, S>;
 
 
   template <int DIM, IntegrationType I, int S>
   using OpHookeTangent = OpHookeTangentImpl<DIM, GAUSS, DomainEleOp, S>;
 };
 
 template <int DIM>
 MoFEMErrorCode
 addMatBlockOps(MoFEM::Interface &m_field,
                boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
                std::string block_name,
                boost::shared_ptr<MatrixDouble> mat_D_Ptr, Sev sev,
                double scale = 1) {
   MoFEMFunctionBegin;
 
   struct OpMatBlocks : public DomainEleOp {
     OpMatBlocks(boost::shared_ptr<MatrixDouble> m, double bulk_modulus_K,
                 double shear_modulus_G, MoFEM::Interface &m_field, Sev sev,
                 std::vector<const CubitMeshSets *> meshset_vec_ptr,
                 double scale)
         : DomainEleOp(NOSPACE, DomainEleOp::OPSPACE), matDPtr(m),
           bulkModulusKDefault(bulk_modulus_K),
           shearModulusGDefault(shear_modulus_G), scaleYoungModulus(scale) {
       CHK_THROW_MESSAGE(extractBlockData(m_field, meshset_vec_ptr, sev),
                         "Can not get data from block");
     }
 
     MoFEMErrorCode doWork(int side, EntityType type,
                           EntitiesFieldData::EntData &data) {
       MoFEMFunctionBegin;
 
       for (auto &b : blockData) {
 
         if (b.blockEnts.find(getFEEntityHandle()) != b.blockEnts.end()) {
           CHKERR getMatDPtr(matDPtr, b.bulkModulusK * scaleYoungModulus,
                             b.shearModulusG * scaleYoungModulus);
           MoFEMFunctionReturnHot(0);
         }
       }
 
       CHKERR getMatDPtr(matDPtr, bulkModulusKDefault * scaleYoungModulus,
                         shearModulusGDefault * scaleYoungModulus);
       MoFEMFunctionReturn(0);
     }
 
   private:
     boost::shared_ptr<MatrixDouble> matDPtr;
     const double scaleYoungModulus;
 
     struct BlockData {
       double bulkModulusK;
       double shearModulusG;
       Range blockEnts;
     };
 
     double bulkModulusKDefault;
     double shearModulusGDefault;
     std::vector<BlockData> blockData;
 
     MoFEMErrorCode
     extractBlockData(MoFEM::Interface &m_field,
                      std::vector<const CubitMeshSets *> meshset_vec_ptr,
                      Sev sev) {
       MoFEMFunctionBegin;
 
       for (auto m : meshset_vec_ptr) {
         MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock") << *m;
         std::vector<double> block_data;
         CHKERR m->getAttributes(block_data);
         if (block_data.size() != 2) {
           SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                   "Expected that block has two attribute");
         }
         auto get_block_ents = [&]() {
           Range ents;
           CHKERR
           m_field.get_moab().get_entities_by_handle(m->meshset, ents, true);
           return ents;
         };
 
         double young_modulus = block_data[0];
         double poisson_ratio = block_data[1];
         double bulk_modulus_K = young_modulus / (3 * (1 - 2 * poisson_ratio));
         double shear_modulus_G = young_modulus / (2 * (1 + poisson_ratio));
 
         MOFEM_TAG_AND_LOG("WORLD", sev, "MatBlock")
             << "E = " << young_modulus << " nu = " << poisson_ratio;
 
         blockData.push_back(
             {bulk_modulus_K, shear_modulus_G, get_block_ents()});
       }
       MOFEM_LOG_CHANNEL("WORLD");
       MoFEMFunctionReturn(0);
     }
 
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
 
   double E = 1.0;
   double nu = 0.3;
 
   CHKERR PetscOptionsBegin(PETSC_COMM_WORLD, "", "", "none");
   CHKERR PetscOptionsScalar("-young_modulus", "Young modulus", "", E, &E,
                             PETSC_NULL);
   CHKERR PetscOptionsScalar("-poisson_ratio", "poisson ratio", "", nu, &nu,
                             PETSC_NULL);
   ierr = PetscOptionsEnd();
 
   double bulk_modulus_K = E / (3 * (1 - 2 * nu));
   double shear_modulus_G = E / (2 * (1 + nu));
   pip.push_back(new OpMatBlocks(
       mat_D_Ptr, bulk_modulus_K, shear_modulus_G, m_field, sev,
 
       // Get blockset using regular expression
       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(
 
           (boost::format("%s(.*)") % block_name).str()
 
               )),
       scale
 
       ));
 
   MoFEMFunctionReturn(0);
 }
 
 template <int DIM, IntegrationType I, typename DomainEleOp>
 auto commonDataFactory(
     MoFEM::Interface &m_field,
     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
     std::string field_name, std::string block_name, Sev sev, double scale = 1) {
 
   auto common_ptr = boost::make_shared<HenckyOps::CommonData>();
   common_ptr->matDPtr = boost::make_shared<MatrixDouble>();
   common_ptr->matGradPtr = boost::make_shared<MatrixDouble>();
 
   CHK_THROW_MESSAGE(addMatBlockOps<DIM>(m_field, pip, block_name,
                                         common_ptr->matDPtr, sev, scale),
                     "addMatBlockOps");
 
   using H = HenckyIntegrators<DomainEleOp>;
 
   pip.push_back(new OpCalculateVectorFieldGradient<DIM, DIM>(
       field_name, common_ptr->matGradPtr));
  

   // Assumes constant D matrix per entity
   pip.push_back(new typename H::template OpCalculateCauchyStress<DIM, I, 0>(
       field_name, common_ptr));
 
   return common_ptr;
 }
 
 template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
 MoFEMErrorCode opFactoryDomainRhs(
     MoFEM::Interface &m_field,
     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
     std::string field_name, boost::shared_ptr<HookeOps::CommonData> common_ptr,
     Sev sev) {
   MoFEMFunctionBegin;
 
   using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
       A>::template LinearForm<I>;
   using OpInternalForcePiola =
       typename B::template OpGradTimesTensor<1, DIM, DIM>;
   pip.push_back(
       new OpInternalForcePiola("U", common_ptr->getMatFirstPiolaStress()));
 
   MoFEMFunctionReturn(0);
 }
 
 template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
 MoFEMErrorCode opFactoryDomainRhs(
     MoFEM::Interface &m_field,
     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
     std::string field_name, std::string block_name, Sev sev, double scale = 1) {
   MoFEMFunctionBegin;
 
   auto common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
       m_field, pip, field_name, block_name, sev, scale);
   CHKERR opFactoryDomainRhs<DIM, A, I, DomainEleOp>(m_field, pip, field_name,
                                                     common_ptr, sev);
 
   MoFEMFunctionReturn(0);
 }
 
 template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
 MoFEMErrorCode opFactoryDomainLhs(
     MoFEM::Interface &m_field,
     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
     std::string field_name, boost::shared_ptr<HookeOps::CommonData> common_ptr,
     Sev sev) {
   MoFEMFunctionBegin;
 
   using B = typename FormsIntegrators<DomainEleOp>::template Assembly<
       A>::template BiLinearForm<I>;
   using OpKPiola = typename B::template OpGradTensorGrad<1, DIM, DIM, 1>;
 
   using H = HenckyIntegrators<DomainEleOp>;
   // Assumes constant D matrix per entity
   pip.push_back(
       new typename H::template OpHenckyTangent<DIM, I, 0>(field_name, common_ptr));
   pip.push_back(
       new OpKPiola(field_name, field_name, common_ptr->getMatTangent()));
 
   MoFEMFunctionReturn(0);
 }
 
 template <int DIM, AssemblyType A, IntegrationType I, typename DomainEleOp>
 MoFEMErrorCode opFactoryDomainLhs(
     MoFEM::Interface &m_field,
     boost::ptr_deque<ForcesAndSourcesCore::UserDataOperator> &pip,
     std::string field_name, std::string block_name, Sev sev, double scale = 1) {
   MoFEMFunctionBegin;
 
   auto common_ptr = commonDataFactory<DIM, I, DomainEleOp>(
       m_field, pip, field_name, block_name, sev, scale);
   CHKERR opFactoryDomainLhs<DIM, A, I, DomainEleOp>(m_field, pip, field_name,
                                                     common_ptr, sev);
 
   MoFEMFunctionReturn(0);
 }
 } // namespace HookeOps
 
 #endif // __HOOKE_OPS_HPP__