/** \file moabField_Core.hpp
 * \brief Core moabField::FEMethod class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MOABFEMETHOD_COMPLEXFORLAZY_HPP__
#define __MOABFEMETHOD_COMPLEXFORLAZY_HPP__

#include "moabField.hpp"
#include "Core_dataStructures.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** \brief The user interface for NonLineae FE method (tangent is calulated
 * using complex direvatives)
 * 
 * This class give user some data structures and methods on those that
 * structures which could be useful. Method of calulating complex direvatives
 * is numerically inefficient, however it is simple to use.
 *
 * 
*/
struct FEMethod_ComplexForLazy: public FEMethod_UpLevelStudent {

  FEMethod_ComplexForLazy(Interface& _moab,int _verbose = 0): 
    FEMethod_UpLevelStudent(_moab,_verbose) {
    pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  }
    
  ErrorCode rval;  
  PetscErrorCode ierr;
  ParallelComm* pcomm;

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;

  vector<vector<DofIdx> > RowGlob;
  vector<vector<DofIdx> > ColGlob;
  vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
  vector<vector<ublas::matrix<FieldData> > > rowDiffNMatrices;
  vector<vector<ublas::matrix<FieldData> > > rowBMatrices;
  vector<vector<ublas::matrix<FieldData> > > colNMatrices;
  vector<vector<ublas::matrix<FieldData> > > colDiffNMatrices;
  vector<vector<ublas::matrix<FieldData> > > colBMatrices;

  PetscErrorCode GetIndices();
  PetscErrorCode GetTangent();
  PetscErrorCode GerResidual();


};

}

#endif //__MOABFEMETHOD_COMPLEXFORLAZY_HPP__
