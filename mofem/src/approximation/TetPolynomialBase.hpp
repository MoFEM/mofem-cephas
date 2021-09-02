/** \file TetPolynomialBase.hpp
\brief Implementation of Ainsworth-Coyle / Demkowicz or any other H1, Hcurl,
Hdiv and L2 base on tetrahedral

*/

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#ifndef __TETPOLYNOMIALBASE_HPP__
#define __TETPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on tetrahedral
 *
 * \ingroup mofem_base_functions
 */
struct TetPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  TetPolynomialBase() = default;
  ~TetPolynomialBase() = default;

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  /**
   * @brief Get base functions for H1 space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  /**
   * @brief Get base functions for L2 space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hdiv space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  /**
   * @brief Get base functions for Hcurl space
   *
   * @param pts matrix of intergation pts
   * @return MoFEMErrorCode
   *
   * \note matrix of integration points on rows has local coordinates of finite
   * element on columns are integration pts.
   */
  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  ublas::matrix<MatrixDouble> N_face_edge;
  ublas::vector<MatrixDouble> N_face_bubble;
  ublas::vector<MatrixDouble> N_volume_edge;
  ublas::vector<MatrixDouble> N_volume_face;
  MatrixDouble N_volume_bubble;

  ublas::matrix<MatrixDouble> diffN_face_edge;
  ublas::vector<MatrixDouble> diffN_face_bubble;
  ublas::vector<MatrixDouble> diffN_volume_edge;
  ublas::vector<MatrixDouble> diffN_volume_face;
  MatrixDouble diffN_volume_bubble;

private:
  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivAinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);

  MatrixInt senseFaceAlpha;
};

} // namespace MoFEM

#endif //__TETPOLYNOMIALBASE_HPP__
