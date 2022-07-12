/** \file TetPolynomialBase.hpp
\brief Implementation of Ainsworth-Coyle / Demkowicz or any other H1, Hcurl,
Hdiv and L2 base on tetrahedral

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

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
