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
  ublas::vector<int> edgesSense;
  ublas::vector<int> facesSense;
  ublas::vector<int> edgesOrder;
  ublas::vector<int> facesOrder;
  ublas::vector<int> rowNodesIndices,colNodesIndices;
  ublas::vector<ublas::vector<int> > rowEdgesIndcies;
  ublas::vector<ublas::vector<int> > colEdgesIndcies;
  ublas::vector<ublas::vector<int> > rowFacesIndcies;
  ublas::vector<ublas::vector<int> > colFacesIndcies;
  ublas::vector<int> rowTetIndcies,colTetIndcies;
  ublas::matrix<int> facesNodes;
  ublas::matrix<FieldData> Nnodes_H1,diffNnodes_H1;
  vector<ublas::matrix<FieldData> > Nedges_H1;
  vector<ublas::matrix<FieldData> > Nfaces_H1;
  ublas::matrix<FieldData> Nvolume_H1;
  ublas::matrix<FieldData> Nnodes_L2,diffNnodes_L2;
  ublas::matrix<FieldData> Nvolume_L2;
  int volumeOrder;
};

struct ForcesAndSurcesCore: public FieldInterface::FEMethod {


  dataForcesAndSurcesCore &data;

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField,dataForcesAndSurcesCore &_data): 
    mField(_mField),data(_data) {};

  PetscErrorCode getEdgesSense();
  PetscErrorCode getFacesSense();
  PetscErrorCode getEdgesOrder(const string &field_name);
  PetscErrorCode getFacesOrder(const string &field_name);


  PetscErrorCode getOrderVolume(const string &field_name);
  PetscErrorCode getNodesIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);
  PetscErrorCode getRowNodesIndices(const string &field_name);
  PetscErrorCode getColNodesIndices(const string &field_name);
  PetscErrorCode getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices);
  PetscErrorCode getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,ublas::vector<ublas::vector<int> > &indices);
  PetscErrorCode getEdgeRowIndices(const string &field_name);
  PetscErrorCode getEdgeColIndices(const string &field_name);
  PetscErrorCode getFacesRowIndices(const string &field_name);
  PetscErrorCode getFacesColIndices(const string &field_name);
  PetscErrorCode getTetRowIndices(const string &field_name);
  PetscErrorCode getTetColIndices(const string &field_name);

  PetscErrorCode getFaceNodes();

  PetscErrorCode shapeTETFunctions_H1(
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTETFunctions_L2(const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTRIFunctions_H1(
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM);

  
};

struct mult_H1_H1 {


  dataForcesAndSurcesCore &data_row,&data_col;
  mult_H1_H1(
    dataForcesAndSurcesCore &_data_row,
    dataForcesAndSurcesCore &_data_col): 
    data_row(_data_row),data_col(data_col) {};

  PetscErrorCode calculate_H1_H1_nonsymetric(double *G_W);
  PetscErrorCode calculate_H1_H1_symmetric(double *G_W);
    
};

}

#endif //__CORE_FORCES_AND_SURCES_HPP

