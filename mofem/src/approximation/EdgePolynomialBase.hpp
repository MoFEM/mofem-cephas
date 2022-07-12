/** \file EdgePolynomialBase.hpp
\brief Implementation of base on tetrahedral for H1 bases.

TODO:
\todo L2 base on edge

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

#ifndef __EDGEPOLYNOMIALBASE_HPP__
#define __EDGEPOLYNOMIALBASE_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on tetrahedral
 *
 * \ingroup mofem_base_functions
 */
struct EdgePolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  EdgePolynomialBase();
  ~EdgePolynomialBase();

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  EntPolynomialBaseCtx *cTx;

  VectorDouble L, diffL;

  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  MoFEMErrorCode getValueH1AinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueH1DemkowiczBase(MatrixDouble &pts);

  MoFEMErrorCode getValueH1BernsteinBezierBase(MatrixDouble &pts);

  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  MoFEMErrorCode getValueL2DemkowiczBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurlAinsworthBase(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurlDemkowiczBase(MatrixDouble &pts);
};

} // namespace MoFEM

#endif //__EDGEPOLYNOMIALBASE_HPP__
