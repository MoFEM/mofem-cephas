/** \file BaseFunction.hpp
\brief General implementation of base function

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

#ifndef __BASEFUNCTION_HPP__
#define __BASEFUNCTION_HPP__

namespace MoFEM {

struct BaseFunctionUnknownInterface : public UnknownInterface {

  using UnknownInterface::UnknownInterface;

  virtual MoFEMErrorCode
  query_interface(boost::typeindex::type_index type_index,
                  UnknownInterface **iface) const = 0;

  virtual ~BaseFunctionUnknownInterface() = default;
};

/**
 * \brief Base class used to exchange data between element data structures and
 * class calculating base functions \ingroup mofem_base_functions
 */
struct BaseFunctionCtx : public BaseFunctionUnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 UnknownInterface **iface) const;

  using BaseFunctionUnknownInterface::BaseFunctionUnknownInterface;
};

/**
 * \brief Base class if inherited used to calculate base functions
 * \ingroup mofem_base_functions
 */
struct BaseFunction : public BaseFunctionUnknownInterface {

  MoFEMErrorCode query_interface(boost::typeindex::type_index type_index,
                                 MoFEM::UnknownInterface **iface) const;

  using BaseFunctionUnknownInterface::BaseFunctionUnknownInterface;

  virtual MoFEMErrorCode getValue(MatrixDouble &pts,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr);

  virtual MoFEMErrorCode getValue(MatrixDouble &pts_x, MatrixDouble &pts_t,
                                  boost::shared_ptr<BaseFunctionCtx> ctx_ptr);
};

} // namespace MoFEM

#endif //__BASEFUNCTION_HPP__

/**
 * \defgroup mofem_base_functions Base functions
 *
 * \brief Calculation of base functions at integration points.
 *
 * \ingroup mofem
 **/
