/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __MOISTUREFE_RVELAGRANGE_DISP_HPP__
#define __MOISTUREFE_RVELAGRANGE_DISP_HPP__

#include "ElasticFE_RVELagrange_Disp.hpp"
#include <moab/ParallelComm.hpp>

namespace MoFEM {
  
  struct MoistureFE_RVELagrange_Disp: public ElasticFE_RVELagrange_Disp {
    MoistureFE_RVELagrange_Disp(FieldInterface& _mField,Mat &_A,Vec &_C, Vec& _F,ublas::vector<FieldData> _applied_strain);
    
//    //pre and post processors (commented because same as in MoistureFE_RVELagrange_Disp so no need to repeat)
//    PetscErrorCode preProcess();
//    PetscErrorCode postProcess();

    //Find indices and calculate shape function  and arrange in required form
    virtual PetscErrorCode GetN_and_Indices();
    
    //Calculate H matrix
    virtual PetscErrorCode Get_H_mat();

//    //Calculate and assemble NT x N matrix (commented because same as in MoistureFE_RVELagrange_Disp so no need to repeat)
//    virtual PetscErrorCode Stiffness();
//    virtual PetscErrorCode Lhs();
    
    //Calculate the right hand side vector, i.e. f=D_max * applied_strain and assemble it into the global force vector F
    virtual PetscErrorCode Rhs();
    
    //Loop over all the elements
    PetscErrorCode operator()();
    
  };
  
}

#endif  //__MOISTUREFE_RVELAGRANGE_DISP_HPP__
