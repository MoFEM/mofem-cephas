/** \file EntPolynomialBaseCtx.hpp
\brief Implementation of Ainsworth-Cole H1 base on tetrahedral

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

#ifndef __ENTPOLYNOMIALCTX_HPP__
#define __ENTPOLYNOMIALCTX_HPP__

namespace MoFEM {

struct EntitiesFieldData;
struct FEMethod;

/**
 * \brief Class used to pass element data to calculate base functions on
 * tet,triangle,edge
 *
 * \ingroup mofem_base_functions
 */
struct EntPolynomialBaseCtx : public BaseFunctionCtx {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  PetscErrorCode (*basePolynomialsType0)(int p, double s, double *diff_s,
                                         double *L, double *diffL,
                                         const int dim);

  PetscErrorCode (*basePolynomialsType1)(int p, double alpha, double x,
                                         double t, double *diff_x,
                                         double *diff_t, double *L,
                                         double *diffL, const int dim);

  EntitiesFieldData &dAta;
  const FieldSpace sPace;
  const FieldApproximationBase bAse;
  const std::string fieldName;
  const FieldApproximationBase copyNodeBase;

  EntPolynomialBaseCtx(EntitiesFieldData &data, const FieldSpace space,
                       const FieldApproximationBase base,
                       const FieldApproximationBase copy_node_base = LASTBASE);

  EntPolynomialBaseCtx(EntitiesFieldData &data,
                       const std::string field_name, const FieldSpace space,
                       const FieldApproximationBase base,
                       const FieldApproximationBase copy_node_base = LASTBASE);

protected:
  MoFEMErrorCode setBase();
};

} // namespace MoFEM

#endif //__ENTPOLYNOMIALCTX_HPP__
