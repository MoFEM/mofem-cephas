/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#ifndef __MOABFEMETHOD_DIRIHLETBC_HPP__
#define __MOABFEMETHOD_DIRIHLETBC_HPP__

#include "CoreDataStructures.hpp"
#include "FieldInterface.hpp"
#include "FEMethod_LowLevelStudent.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** 
 * \brief The student user interface for Dirihlet boundary conditions
 * 
*/
struct BaseDirihletBC {

  BaseDirihletBC();


  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirihletBC,vector<DofIdx>& FaceNodeGlobalDofs,vector<vector<DofIdx> > &FaceEdgeGlobalDofs,vector<DofIdx> &FaceGlobalDofs);
  virtual PetscErrorCode SetDirihletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);
  virtual PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  virtual PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);

};

struct CubitDisplacementDirihletBC: public BaseDirihletBC {
  FieldInterface& mField;
  string problem_name;  
  string field_name;

  CubitDisplacementDirihletBC(FieldInterface& _mField,const string _problem_name,const string _field_name); 

  PetscErrorCode ierr;
  ErrorCode rval;

  map<int,Range> bc_map[3];
  map<int,double> bc_map_val[3];

  virtual PetscErrorCode Init();

  PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndicies(
    FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC);
  PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    vector<DofIdx>& DirihletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices);
  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);
  PetscErrorCode SetDirihletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);

};
    
    
//struct CubitTemperatureDirihletBC: public CubitDisplacementDirihletBC {
    //CubitTemperatureDirihletBC(FieldInterface& _mField,const string _problem_name,const string _field_name);
    
//    virtual PetscErrorCode Init() {
//        PetscFunctionBegin;
//            for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|TemperatureSet,it)) {
//
//              temperature_cubit_bc_data mydata;
//              ierr = it->get_cubit_bc_data_structure(mydata); CHKERRQ(ierr);
//              for(int dim = 0;dim<3;dim++) {
//            Range _ents;
//            ierr = it->get_Cubit_msId_entities_by_dimension(mField.get_moab(),dim,_ents,true); CHKERRQ(ierr);
//            if(dim>1) {
//              Range _edges;
//              ierr = mField.get_moab().get_adjacencies(_ents,1,false,_edges,Interface::UNION); CHKERRQ(ierr);
//              _ents.insert(_edges.begin(),_edges.end());
//            }
//            if(dim>0) {
//              Range _nodes;
//              rval = mField.get_moab().get_connectivity(_ents,_nodes,true); CHKERR_PETSC(rval);
//              _ents.insert(_nodes.begin(),_nodes.end());
//            }
//                  
//            if(dim>2) SETERRQ(PETSC_COMM_SELF,1,"not yet implemented");
//            
//              (bc_map[0])[it->get_msId()].insert(_ents.begin(),_ents.end());
//              (bc_map_val[0])[it->get_msId()] = mydata.data.value1;
//                  
//              }
//            }
//        }

//    };
}
#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
