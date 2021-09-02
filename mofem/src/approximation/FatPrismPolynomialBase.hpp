/** \file FatPrismPolynomialBase.hpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

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

#ifndef __FATPRISMPOLYNOMIALBASE_HPP__
#define __FATPRISMPOLYNOMIALBASE_HPP__

namespace MoFEM {

struct NumeredEntFiniteElement;

/**
 * \brief Class used to pass element data to calculate base functions on fat
 * prism
 *
 * \ingroup mofem_base_functions
 * FIXME: Need moab and mofem finite element structure to work (that not
 * perfect)
 */
struct FatPrismPolynomialBaseCtx : public EntPolynomialBaseCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  DataForcesAndSourcesCore &dataTrianglesOnly;
  DataForcesAndSourcesCore &dataTroughThickness;

  MatrixDouble &gaussPtsTrianglesOnly;
  MatrixDouble &gaussPtsThroughThickness;

  moab::Interface &mOab;
  const NumeredEntFiniteElement *fePtr;

  FatPrismPolynomialBaseCtx(
      DataForcesAndSourcesCore &data,
      DataForcesAndSourcesCore &data_triangles_only,
      DataForcesAndSourcesCore &data_trough_thickness,
      MatrixDouble &gauss_pts_triangles_only,
      MatrixDouble &gauss_pts_through_thickness, moab::Interface &moab,
      const NumeredEntFiniteElement *fe_ptr, const FieldSpace space,
      const FieldApproximationBase base,
      const FieldApproximationBase copy_node_base = LASTBASE);

  ~FatPrismPolynomialBaseCtx();
};

/**
 * \brief Calculate base functions on tetrahedral
 * \ingroup mofem_base_functions
 * FIXME: Need moab and mofem finite element structure to work (that not
 * perfect)
 */
struct FatPrismPolynomialBase : public BaseFunction {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  FatPrismPolynomialBase();
  ~FatPrismPolynomialBase();

  MoFEMErrorCode getValue(MatrixDouble &pts,
                          boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

private:
  FatPrismPolynomialBaseCtx *cTx;

  MoFEMErrorCode getValueH1TrianglesOnly();

  MoFEMErrorCode getValueH1ThroughThickness();

  MoFEMErrorCode getValueH1(MatrixDouble &pts);

  MoFEMErrorCode getValueL2(MatrixDouble &pts);

  MoFEMErrorCode getValueHdiv(MatrixDouble &pts);

  MoFEMErrorCode getValueHcurl(MatrixDouble &pts);

  // int faceNodes[2][3];
  MatrixDouble N;
  MatrixDouble diffN;
};

} // namespace MoFEM

#endif //__FATPRISMPOLYNOMIALBASE_HPP__
