/** \file CGGTonsorialBubbleBase.hpp

  \brief Implementation of tonsorial bubble base div(v) = 0.

  Implementation is based and motiveted by \cite cockburn2010new. This base
  is used to approximate stresses using Hdiv base with weakly enforced
  symmetry.

*/

namespace MoFEM {

/**
 * @brief Calculate CGGT tonsorial bubble base
 * 
 * See details in \cite cockburn2010new
 * 
 * @param p polynomial order
 * @param N shape functions
 * @param diffN direvatives of shape functions
 * @param l2_base base functions for l2 space
 * @param diff_l2_base direvatives base functions for l2 space
 * @param phi returned base functions
 * @param gdim number of integration points
 * @return MoFEMErrorCode 
 */
MoFEMErrorCode
CGG_BubbleBase_MBTET(const int p, const double *N, const double *diffN,
                     const double *l2_base, const double *diff_l2_base,
                     FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> &phi,
                     const int gdim);

}