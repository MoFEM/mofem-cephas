/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * DirichletBCFromBlockSetFEMethodPreAndPostProc implmented by Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
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

#ifndef __MOABFEMETHOD_DIRICHLETBC_HPP__
#define __MOABFEMETHOD_DIRICHLETBC_HPP__

using namespace boost::numeric;

namespace MoFEM {

struct DisplacementBCFEMethodPreAndPostProc: public FEMethod {

  FieldInterface& mField;
  const string fieldName;

  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,
    Mat _Aij,Vec _X,Vec _F): mField(_mField),fieldName(_field_name) {
    snes_B = _Aij;
    snes_x = _X;
    snes_f = _F;
    ts_B = _Aij;
    ts_u = _X;
    ts_F = _F;
  };

  DisplacementBCFEMethodPreAndPostProc(FieldInterface& _mField,const string &_field_name): 
    mField(_mField),fieldName(_field_name) {
    snes_B = PETSC_NULL;
    snes_x = PETSC_NULL;
    snes_f = PETSC_NULL;
    ts_B = PETSC_NULL;
    ts_u = PETSC_NULL;
    ts_F = PETSC_NULL;
  }; 

  PetscErrorCode ierr;
  ErrorCode rval;

  map<DofIdx,FieldData> map_zero_rows;
  vector<int> dofsIndices;
  vector<double> dofsValues;
  vector<double> dofsXValues;
  virtual PetscErrorCode iNitalize();

  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

};


struct SpatialPositionsBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,Mat _Aij,Vec _X,Vec _F): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F) {}

  SpatialPositionsBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name) {}

  vector<string> fixFields;

  ublas::vector<double> cOords;
  PetscErrorCode iNitalize();

};

struct TemperatureBCFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {

  TemperatureBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,Mat _Aij,Vec _X,Vec _F): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F) {}

  TemperatureBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name) {}

  PetscErrorCode iNitalize();

};

  
struct FixBcAtEntities: public DisplacementBCFEMethodPreAndPostProc {

  Range &eNts;
  vector<string> fieldNames;
  FixBcAtEntities(
    FieldInterface& _mField,const string &_field_name,Mat _Aij,Vec _X,Vec _F,Range &ents): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F),eNts(ents) {
    fieldNames.push_back(fieldName);
  }

  FixBcAtEntities(
    FieldInterface& _mField,const string &_field_name,Range &ents): 
    DisplacementBCFEMethodPreAndPostProc(_mField,_field_name),eNts(ents) {
    fieldNames.push_back(fieldName);
  }


  PetscErrorCode iNitalize();
  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

};
  
  
/**
  * Implemntaiton of generalised Dirichlet Boundary Conditions from CUBIT Blockset
  * (or not using CUBIT buildin boundary conditons, e.g. Temprature or Displacements etc).
  * It can work for any Problem rank (1,2,3)
**/
struct DirichletBCFromBlockSetFEMethodPreAndPostProc: public DisplacementBCFEMethodPreAndPostProc {
  const string _blockset_name;
  const string blocksetName;
  DirichletBCFromBlockSetFEMethodPreAndPostProc(FieldInterface& _mField,const string &_field_name,const string &_blockset_name,Mat _Aij,Vec _X,Vec _F):
  DisplacementBCFEMethodPreAndPostProc(_mField,_field_name,_Aij,_X,_F),blocksetName(_blockset_name) {}
  
  DirichletBCFromBlockSetFEMethodPreAndPostProc(FieldInterface& _mField,const string &_field_name):
  DisplacementBCFEMethodPreAndPostProc(_mField,_field_name),blocksetName(_blockset_name) {}
  
  PetscErrorCode iNitalize();
  
};

    
}
#endif //__MOABFEMETHOD_DIRICHLETBC_HPP__
