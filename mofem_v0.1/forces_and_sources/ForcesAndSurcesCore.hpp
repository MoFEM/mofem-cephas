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
  //sense
  ublas::vector<int> edgesSense;
  ublas::vector<int> facesSense;
  //order
  ublas::vector<ApproximationOrder> edgesOrder;
  ublas::vector<ApproximationOrder> facesOrder;
  DofIdx volumeOrder;
  //indices
  ublas::vector<DofIdx> nodesIndices;
  ublas::vector<ublas::vector<DofIdx> > edgesIndcies;
  ublas::vector<ublas::vector<DofIdx> > facesIndices;
  ublas::vector<DofIdx> volumeIndices;
  //face volume nodes
  ublas::matrix<DofIdx> facesNodes;
  //field data
  //h1
  ublas::matrix<FieldData> nodesNH1,diffNodesNH1;
  vector<ublas::matrix<FieldData> > edgesNH1;
  vector<ublas::matrix<FieldData> > facesNH1;
  ublas::matrix<FieldData> volumeNH1;
  //l2
  ublas::matrix<FieldData> nodesNL2,diffNodesNL2;
  ublas::matrix<FieldData> volumeNL2;
};

struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};

  PetscErrorCode getEdgesSense(dataForcesAndSurcesCore &data);
  PetscErrorCode getFacesSense(dataForcesAndSurcesCore &data);
  PetscErrorCode getEdgesOrder(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesOrder(dataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getOrderVolume(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getNodesIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);
  PetscErrorCode getRowNodesIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name);
  PetscErrorCode getColNodesIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name);
  PetscErrorCode getTypeIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<int> &indices);
  PetscErrorCode getTypeIndices(
    dataForcesAndSurcesCore &data,
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,ublas::vector<ublas::vector<int> > &indices);
  PetscErrorCode getEdgeRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgeColIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesColIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetRowIndices(dataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetColIndices(dataForcesAndSurcesCore &data,const string &field_name);

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
    ublas::vector<DofIdx> &row_indices,ublas::vector<DofIdx> &col_indices,
    ublas::matrix<FieldData> &rows_N,ublas::matrix<FieldData> &cols_N) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()(dataForcesAndSurcesCore &row_data,dataForcesAndSurcesCore &col_data);

  virtual PetscErrorCode doWork(
    int side,
    EntityType type,
    ublas::vector<DofIdx> &indices,
    ublas::matrix<FieldData> &N) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()(dataForcesAndSurcesCore &data);


};


}

#endif //__CORE_FORCES_AND_SURCES_HPP

