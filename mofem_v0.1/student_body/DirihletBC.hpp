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

#include "FieldInterface.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct DisplacementBCFEMethodPreAndPostProc: public FieldInterface::FEMethod {

  FieldInterface& mField;
  const string fieldName;

  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,
    Mat &_Aij,Vec _X,Vec _F): mField(_mField),fieldName(_field_name) {
    snes_B = &_Aij;
    snes_x = _X;
    snes_f = _F;
    ts_B = &_Aij;
  };

  DisplacementBCFEMethodPreAndPostProc(FieldInterface& _mField,const string &_field_name): 
    mField(_mField),fieldName(_field_name) {}; 

  PetscErrorCode ierr;
  ErrorCode rval;

  map<DofIdx,FieldData> map_zero_rows;
  vector<int> dofsIndices;
  vector<double> dofsValues;

  virtual PetscErrorCode iNitalize();

  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

};


struct SpatialPositionsBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,Mat &_Aij,Vec _X,Vec _F): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F) {}

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name) {}

  vector<string> fixFields;

  ublas::vector<double> cOords;
  PetscErrorCode iNitalize();

};

struct FixMaterialPoints: public DisplacementBCFEMethodPreAndPostProc {

  Range &eNts;
  vector<string> fieldNames;
  FixMaterialPoints(
    FieldInterface& _mField,const string &_field_name,Mat &_Aij,Vec _X,Vec _F,Range &ents): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F),eNts(ents) {
    fieldNames.push_back(fieldName);
  }

  FixMaterialPoints(
    FieldInterface& _mField,const string &_field_name,Range &ents): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name),eNts(ents) {
    fieldNames.push_back(fieldName);
  }


  PetscErrorCode iNitalize();
  PetscErrorCode preProcess();

};

///******************************************

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
    FieldInterface::FEMethod *fe_method_ptr,
      vector<vector<DofIdx> > &RowGlobDofs,
      vector<vector<DofIdx> > &ColGlobDofs,
      vector<DofIdx>& DirihletBC);

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& RowGlobDofs,vector<DofIdx>& DirihletBC);
  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& ColGlobDofs,vector<DofIdx>& DirihletBC);

  virtual PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
    FieldInterface::FEMethod *fe_method_ptr,
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
    FieldInterface::FEMethod *fe_method_ptr,
    vector<DofIdx>& DirihletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices);
  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
  PetscErrorCode SetDirihletBC_to_RHS(FieldInterface::FEMethod *fe_method_ptr,Vec F);
  PetscErrorCode SetDirihletBC_to_FieldData(FieldInterface::FEMethod *fe_method_ptr,Vec D);

};

struct CubitTemperatureDirihletBC: public CubitDisplacementDirihletBC {

  CubitTemperatureDirihletBC(FieldInterface& _mField,const string _problem_name,const string _field_name); 
  PetscErrorCode Init();

};
    
struct CubitDisplacementDirihletBC_ZerosRowsColumns: public CubitDisplacementDirihletBC {
    CubitDisplacementDirihletBC_ZerosRowsColumns(FieldInterface& _mField,const string _problem_name,const string _field_name);
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesRow(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesCol(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndicies(
        FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlobDofs,vector<vector<DofIdx> > &ColGlobDofs,vector<DofIdx>& DirihletBC){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(
        FieldInterface::FEMethod *fe_method_ptr,
        vector<DofIdx>& DirihletBC,vector<DofIdx> &FaceNodeIndices, vector<vector<DofIdx> > &FaceEdgeIndices, vector<DofIdx> &FaceIndices){
        PetscFunctionBegin;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode SetDirihletBC_to_MatrixDiagonal(FieldInterface::FEMethod *fe_method_ptr,Mat Aij);
    
};
    

    
}
#endif //__MOABFEMETHOD_DIRIHLETBC_HPP__
