/** \file TriPolynomialBase.hpp
\brief Implementation of  H1, Hcurl base on triangle

*/

#ifndef __H1TRIPOLYNOMIAL_HPP__
#define __H1TRIPOLYNOMIAL_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on triangle
 *
 * \ingroup mofem_base_functions
 */
struct TriPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  TriPolynomialBase() = default;
  virtual ~TriPolynomialBase() = default;

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

  /**
   * @brief Set map of dof to side number
   *
   * That is used for broken space to establish connection between dofs in the
   * interior of element/entity and side of element/entity to which that dof is
   * associated. That depends on implementation of the base for given space, and
   * has to be implemented while implementing base function for given space.
   *
   * @param space
   * @param continuity
   * @param base
   * @param DofsSideMap
   * @return MoFEMErrorCode
   */
  static MoFEMErrorCode setDofsSideMap(const FieldSpace space,
                                       const FieldContinuity continuity,
                                       const FieldApproximationBase base,
                                       DofsSideMap &);

private:
  EntPolynomialBaseCtx *cTx;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);
  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2(MatrixDouble &pts);
  MoFEMErrorCode getValueL2AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2BernsteinBezierBase(MatrixDouble &pts);

  ublas::matrix<MatrixDouble> N_face_edge;
  ublas::vector<MatrixDouble> N_face_bubble;
  ublas::matrix<MatrixDouble> diffN_face_edge;
  ublas::vector<MatrixDouble> diffN_face_bubble;

  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivAinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurlAinsworthBrokenBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBrokenBase(MatrixDouble &pts);

  static MoFEMErrorCode setDofsSideMapHcurl(const FieldSpace space,
                                            const FieldContinuity continuity,
                                            const FieldApproximationBase base,
                                            DofsSideMap &dofs_side_map);
};

} // namespace MoFEM

#endif //__H1TRIPOLYNOMIAL_HPP__
