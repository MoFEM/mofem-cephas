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

#ifndef __BCS_RVELAGRANGE_PERIIODIC_HPP
#define __BCS_RVELAGRANGE_PERIIODIC_HPP

namespace MoFEM {
  
  struct BCs_RVELagrange_Periodic: public BCs_RVELagrange_Trac {
      BCs_RVELagrange_Periodic(FieldInterface &m_field): BCs_RVELagrange_Trac(m_field){}
    
    
//     PetscErrorCode setRVEBCsOperators(string field_name,string lagrang_field_name,Mat _Aij, Vec _F1, Vec _F2, Vec _F3, Vec _F4, Vec _F5, Vec _F6, const string mesh_nodals_positions) {
//      PetscFunctionBegin;
//      
//      bool ho_geometry = false;
//      if(mField.check_field(mesh_nodals_positions)) {
//        ho_geometry = true;
//      }
////      cout<<"Hi 1 from setRVEBCsOperators "<<endl;
//      map<int,RVEBC_Data>::iterator sit = setOfRVEBC.begin();
//      for(;sit!=setOfRVEBC.end();sit++) {
////        cout<<"Hi from setRVEBCsOperators "<<endl;
//        //LHS
//        feRVEBCLhs.getRowColOpPtrVector().push_back(new OpRVEBCsLhs(field_name,lagrang_field_name, _Aij, sit->second, common_functions, ho_geometry));
//        
//        //RHS
//        //Caclculte D_mat
//        feRVEBCRhs.getColOpPtrVector().push_back(new OpRVEBCsRhs_Cal(field_name, sit->second, commonData, common_functions, ho_geometry));
//        //Caclculte f and assemplbe
//        feRVEBCRhs.getRowOpPtrVector().push_back(new OpRVEBCsRhs_Assemble(lagrang_field_name, _F1, _F2, _F3, _F4, _F5, _F6, sit->second, commonData, ho_geometry));
//      }
//      PetscFunctionReturn(0);
//    }
    
    
    
  };
  
}

#endif
