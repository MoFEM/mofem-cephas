/** \file TriPolynomialBase.hpp
\brief Implementation of  H1, Hcurl base on triangle

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

#ifndef __H1TRIPOLYNOMIAL_HPP__
#define __H1TRIPOLYNOMIAL_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on triangle
 *
 * \ingroup mofem_base_functions
 */
struct TriPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index, UnknownInterface **iface) const;

  TriPolynomialBase();
  ~TriPolynomialBase();

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

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

};

} // namespace MoFEM

#endif //__H1TRIPOLYNOMIAL_HPP__
