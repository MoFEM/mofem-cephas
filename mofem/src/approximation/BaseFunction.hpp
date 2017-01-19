/** \file BaseFunction.hpp
\brief General implementation of base function

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

#ifndef __BASEFUNCTION_HPP__
#define __BASEFUNCTION_HPP__

namespace MoFEM {

  static const MOFEMuuid IDD_UNKNOWN_BASE_FUNCTION = MOFEMuuid(BitIntefaceId(UNKNOWN_BASE_FUNCTION_INTERFACE));

  /**
   * \brief Base class used to exchange data between element data structures and class calculating base functions
   * \ingroup mofem_base_functions
   */
  struct BaseFunctionCtx: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    BaseFunctionCtx() {}
    ~BaseFunctionCtx() {}

  };

  /**
   * \brief Base class if inherited used to calculate base functions
   * \ingroup mofem_base_functions
   */
  struct BaseFunction: public UnknownInterface {

    PetscErrorCode queryInterface(const MOFEMuuid& uuid,MoFEM::UnknownInterface** iface);

    BaseFunction() {}
    ~BaseFunction() {}

    virtual PetscErrorCode getValue(
      ublas::matrix<double> &pts,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );

    virtual PetscErrorCode getValue(
      ublas::matrix<double> &pts_x,
      ublas::matrix<double> &pts_t,
      boost::shared_ptr<BaseFunctionCtx> ctx_ptr
    );


  };

}

#endif //__BASEFUNCTION_HPP__

/***************************************************************************//**
 * \defgroup mofem_base_functions Base functions
 *
 * \brief Calculation of base functions at integration points.
 *
 * \ingroup mofem
 ******************************************************************************/
