/** \file QuadPolynomialBase.hpp
\brief Implementation of H1 base on a quad face

\todo Quad element can be integrated exploiting tonsorial product. Current
implementation do not take that opportunity. That can be viewed as a bug. 

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

#ifndef __H1QUADPOLYNOMIAL_HPP__
#define __H1QUADPOLYNOMIAL_HPP__

namespace MoFEM {

/**
 * \brief Calculate base functions on triangle
 *
 * \ingroup mofem_base_functions
 */
struct QuadPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(const MOFEMuuid &uuid,
                                 BaseFunctionUnknownInterface **iface) const;

  QuadPolynomialBase();
  ~QuadPolynomialBase();

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

  ublas::matrix<MatrixDouble> N_face_edge;
  ublas::vector<MatrixDouble> N_face_bubble;
  ublas::matrix<MatrixDouble> diffN_face_edge;
  ublas::vector<MatrixDouble> diffN_face_bubble;

};

} // namespace MoFEM

#endif //__H1QUADPOLYNOMIAL_HPP__
