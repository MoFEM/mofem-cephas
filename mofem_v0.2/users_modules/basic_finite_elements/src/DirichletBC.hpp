/* \file Dirichlet.hpp
 * \brief Implementation of Dirichlet boundary conditions
 *
 *
 * Structures and method in this file erase rows and column, set value on
 * matrix diagonal and on the right hand side vector to enforce boundary
 * condition.
 *
 * Current implementation is suboptimal, classes name too long. Need to
 * rethinking and improved, more elegant and more efficient implementation.
 *
 */

/* Notes:
 * DirichletBCFromBlockSetFEMethodPreAndPostProc implemented by Zahur Ullah (Zahur.Ullah@glasgow.ac.uk)
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

#ifndef __DIRICHLETBC_HPP__
#define __DIRICHLETBC_HPP__

using namespace boost::numeric;

/** \brief Set Dirichlet boundary conditions on displacements
  * \ingroup Dirichlet_bc
  */
struct DisplacementBCFEMethodPreAndPostProc: public FEMethod {

  FieldInterface& mField;
  const string fieldName;			///< field name to set Dirichlet BC
  double dIag;					///< diagonal value set on zeroed column and rows

  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name,
    Mat _Aij,Vec _X,Vec _F);
  DisplacementBCFEMethodPreAndPostProc(
    FieldInterface& _mField,const string &_field_name);

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

/** \brief Set Dirichlet boundary conditions on spatial displacements
  * \ingroup Dirichlet_bc
  */
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

/** \brief Fix dofs on entities
  * \ingroup Dirichlet_bc
  */
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


/** \brief blockset boundary conditions
  * \ingroup Dirichlet_bc
  *
  * Implementation of generalised Dirichlet Boundary Conditions from CUBIT Blockset
  * (or not using CUBIT building boundary conditions, e.g. Temperature or Displacements etc).
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

#endif //__DIRICHLETBC_HPP__

/***************************************************************************//**
 * \defgroup Dirichlet_bc Dirichlet boundary conditions
 * \ingroup user_modules
 ******************************************************************************/
