/*/**
 * @file DGProjection.hpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-08-23
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef DGPROJECTION_HPP
#define DGPROJECTION_HPP

namespace MoFEM {

/**
 * @brief Execute "this" element in the operator
 *
 * \code An example of usage in post-processing

using InteriorElem = VolumeElementForcesAndSourcesCore;
using PostProcFe = PostProcBrokenMeshInMoab<InteriorElem>;
auto fe_post_proc = boost::make_shared<PostProcFe>(mField);
auto op_this = new OpLoopThis<InteriorElem>(m_field, "DomainFE", Sev::noisy);
fe_post_proc->getOpPtrVector()->push_back(op_this);

auto fe_physics = op_this->getThisFEPtr();
fe_physics->getRuleHook = [](int, int, int o) {
  return 2 * o; };

int order = 1;

// entity data shared between physical and post proc elements
auto entity_data_l2 = boost::make_shared<EntitiesFieldData>(MBENTITYSET); 
// integrated mass matrix of post proc element
auto mass_ptr = boost::make_shared<MatrixDouble>(); 
// vector of coeff shared between physical and post proc elements
auto coeffs_ptr = boost::make_shared<MatrixDouble>();
// data stored at
// integration points of the physical element and evaluated at integration
// points of the post proc element
auto data_ptr = boost::make_shared<MatrixDouble>(); 
fe_physics->getOpPtrVector()->push_back(new
 OpDGProjectionMassMatrix(order, mass_ptr, entity_data_l2,
  AINSWORTH_LEGENDRE_BASE, L2));

// You need to call operator which will evaluate data_ptr

fe_physics->getOpPtrVector()->push_back(new
  OpCalculateVectorFieldValues("V", data_ptr);
fe_physics->getOpPtrVector()->push_back( new
  OpDGProjectionCoefficients(data_ptr, coeffs_ptr, mass_ptr, entity_data_l2,
  AINSWORTH_LEGENDRE_BASE, L2));

fe_post_proc->getOpPtrVector()->push_back(new
  OpDGProjectionEvaluation(data_ptr, coeffs_ptr, 
  entity_data_l2, AINSWORTH_LEGENDRE_BASE, L2));

 * \endcode
 *
 * @tparam E template for "this" element type
 *
 *
 *
 */
template <typename E>
struct OpLoopThis : public ForcesAndSourcesCore::UserDataOperator {

  using UserDataOperator = ForcesAndSourcesCore::UserDataOperator;

  /**
   * @brief Construct a new OpThis object
   *
   * @param m_field
   * @param fe_name name of "this" element (domain element)
   */
  OpLoopThis(MoFEM::Interface &m_field, const std::string fe_name,
         const LogManager::SeverityLevel sev = Sev::noisy)
      : UserDataOperator(NOSPACE, OPSPACE), thisFEPtr(new E(m_field)),
        feName(fe_name), sevLevel(sev) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {
    MoFEMFunctionBegin;
    CHKERR loopThis(feName, thisFEPtr.get(), VERBOSE, sevLevel);
    MoFEMFunctionReturn(0);
  };

  /**
   * Here user will push operator evaluating field, or data at "this" element
   * integration points.
   */
  boost::ptr_deque<UserDataOperator> &getOpPtrVector() {
    return thisFEPtr->getOpPtrVector();
  }

  boost::shared_ptr<E> &getThisFEPtr() { return thisFEPtr; }

protected:
  const std::string feName;
  boost::shared_ptr<E> thisFEPtr;
  const LogManager::SeverityLevel sevLevel;
};

struct OpDGProjectionMassMatrix : public OpBaseDerivativesBase {
  OpDGProjectionMassMatrix(int order, boost::shared_ptr<MatrixDouble> mass_ptr,
                           boost::shared_ptr<EntitiesFieldData> data_l2,
                           const FieldApproximationBase b, const FieldSpace s,
                           int verb = QUIET, Sev sev = Sev::verbose);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  int baseOrder;
};

struct OpDGProjectionCoefficients : public OpBaseDerivativesBase {

  OpDGProjectionCoefficients(boost::shared_ptr<MatrixDouble> data_ptr,
                             boost::shared_ptr<MatrixDouble> coeffs_ptr,
                             boost::shared_ptr<MatrixDouble> mass_ptr,
                             boost::shared_ptr<EntitiesFieldData> data_l2,
                             const FieldApproximationBase b, const FieldSpace s,
                             const LogManager::SeverityLevel sev = Sev::noisy);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

protected:
  boost::shared_ptr<MatrixDouble> dataPtr;
  boost::shared_ptr<MatrixDouble> coeffsPtr;
};

struct OpDGProjectionEvaluation : public OpDGProjectionCoefficients {

  OpDGProjectionEvaluation(boost::shared_ptr<MatrixDouble> data_ptr,
                           boost::shared_ptr<MatrixDouble> coeffs_ptr,
                           boost::shared_ptr<EntitiesFieldData> data_l2,
                           const FieldApproximationBase b, const FieldSpace s,
                           const LogManager::SeverityLevel sev = Sev::noisy);

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);
};
}; // namespace MoFEM

#endif // DGPROJECTION_HPP