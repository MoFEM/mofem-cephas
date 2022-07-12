/** \file QuadPolynomialBase.hpp
\brief Implementation of H1 base on a quad face

\todo Quad element can be integrated exploiting tonsorial product. Current
implementation do not take that opportunity. That can be viewed as a bug. 

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

#ifndef __H1QUADPOLYNOMIAL_HPP__
#define __H1QUADPOLYNOMIAL_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on triangle
 *
 * \ingroup mofem_base_functions
 */
struct QuadPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index, UnknownInterface **iface) const;

  QuadPolynomialBase() = default;
  ~QuadPolynomialBase() = default;

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);
  MoFEMErrorCode getValueL2(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);
  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);
  MoFEMErrorCode getValueH1DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueL2DemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);
  MoFEMErrorCode getValueHdivDemkowiczBase(MatrixDouble &pts);

  MatrixDouble faceFamily;
  MatrixDouble diffFaceFamily;

};

} // namespace MoFEM

#endif //__H1QUADPOLYNOMIAL_HPP__
