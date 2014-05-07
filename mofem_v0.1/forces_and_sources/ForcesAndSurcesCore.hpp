/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
* --------------------------------------------------------------
*
* DESCRIPTION: FIXME
*
* This is not exactly procedure for linear elatic dynamics, since jacobian is
* evaluated at every time step and snes procedure is involved. However it is
* implemented like that, to test methodology for general nonlinear problem.
*
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

#ifndef __CORE_FORCES_AND_SURCES_HPP
#define __CORE_FORCES_AND_SURCES_HPP

#include "FieldInterface.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace MoFEM {

struct dataForcesAndSurcesCore {

  struct entData {
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::matrix<FieldData> N,diffN;
    entData(): sEnse(0),oRder(0) {};
  };

  entData nOdes;
  vector<entData> eDges;
  vector<entData> fAces;
  entData vOlume;

  ublas::matrix<DofIdx> facesNodes;

  friend ostream& operator<<(ostream& os,const dataForcesAndSurcesCore &e);

};

ostream& operator<<(ostream& os,const dataForcesAndSurcesCore::entData &e);

struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};

  PetscErrorCode getSense(EntityType type,vector<dataForcesAndSurcesCore::entData> &data);
  PetscErrorCode getOrder(EntityType type,vector<dataForcesAndSurcesCore::entData> &data,const string &field_name);

  PetscErrorCode getEdgesSense(dataForcesAndSurcesCore &data);
  PetscErrorCode getFacesSense(dataForcesAndSurcesCore &data);
  PetscErrorCode getEdgesOrder(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesOrder(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getOrderVolume(dataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getNodesIndices(
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);

  PetscErrorCode getRowNodesIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name);

  PetscErrorCode getColNodesIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name);

  PetscErrorCode getTypeIndices(
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<int> &indices);
  PetscErrorCode getTypeIndices(
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,
    vector<dataForcesAndSurcesCore::entData> &data);

  PetscErrorCode getEdgeRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgeColIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesColIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetColIndices(dataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getNodesFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<FieldData> &nodes_indices);
  PetscErrorCode getTypeFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<FieldData> &ent_field_data);
  PetscErrorCode getTypeFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,
    vector<dataForcesAndSurcesCore::entData> &data);

  PetscErrorCode getNodesFieldData(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgeFieldData(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesFieldData(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetFieldData(dataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getFaceNodes(dataForcesAndSurcesCore &data);

  PetscErrorCode shapeTETFunctions_H1(
    dataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTETFunctions_L2(
    dataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTRIFunctions_H1(
    dataForcesAndSurcesCore &data,
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM);

};

struct dataOperator {

  virtual PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type,EntityType col_type,
    dataForcesAndSurcesCore::entData &row_data,
    dataForcesAndSurcesCore::entData &col_data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode opNH1NH1(dataForcesAndSurcesCore &row_data,dataForcesAndSurcesCore &col_data);

  virtual PetscErrorCode doWork(
    int side,
    EntityType type,
    dataForcesAndSurcesCore::entData &data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode opNH1(dataForcesAndSurcesCore &data);


};


}

#endif //__CORE_FORCES_AND_SURCES_HPP

