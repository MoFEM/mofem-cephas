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


struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): mField(_mField) {};

  ublas::vector<int> edgesSense;
  PetscErrorCode getEdgesSense();

  ublas::vector<int> facesSense;
  PetscErrorCode getFacesSense();

  ublas::vector<int> edgesOrder;
  PetscErrorCode getEdgesOrder(const string &field_name);

  ublas::vector<int> facesOrder;
  PetscErrorCode getFacesOrder(const string &field_name);

  int volumeOrder;
  PetscErrorCode getOrderVolume(const string &field_name);

  PetscErrorCode getNodesIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);

  ublas::vector<int> rowNodesIndices,colNodesIndices;
  PetscErrorCode getRowNodesIndices(const string &field_name);
  PetscErrorCode getColNodesIndices(const string &field_name);

  PetscErrorCode getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,int side_number,ublas::vector<int> &indices);
  PetscErrorCode getTypeIndices(const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,vector<ublas::vector<int> > &indices);

  vector<ublas::vector<int> > rowEdgesIndcies,colEdgesIndcies;
  PetscErrorCode getEdgeRowIndices(const string &field_name);
  PetscErrorCode getEdgeColIndices(const string &field_name);

  vector<ublas::vector<int> > rowFacesIndcies,colFacesIndcies;
  PetscErrorCode getFacesRowIndices(const string &field_name);
  PetscErrorCode getFacesColIndices(const string &field_name);

  ublas::vector<int> rowTetIndcies,colTetIndcies;
  PetscErrorCode getTetRowIndices(const string &field_name);
  PetscErrorCode getTetColIndices(const string &field_name);

  ublas::matrix<int> facesNodes;
  PetscErrorCode getFaceNodes();

  ublas::matrix<FieldData> Nnodes_H1,diffNnodes_H1;
  vector<ublas::matrix<FieldData> > Nedges_H1;
  vector<ublas::matrix<FieldData> > Nfaces_H1;
  ublas::matrix<FieldData> Nvolume_H1;

  PetscErrorCode shapeTETFunctions_H1(
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);

  ublas::matrix<FieldData> Nnodes_L2,diffNnodes_L2;
  ublas::matrix<FieldData> Nvolume_L2;

  PetscErrorCode shapeTETFunctions_L2(const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);

  PetscErrorCode shapeTRIFunctions_H1(
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM);

};

struct mult_H1_H1 {

  ublas::matrix<FieldData> &B1_Nnodes,&B1_diffNnodes;
  vector<ublas::matrix<FieldData> > &B1_Nedges;
  vector<ublas::matrix<FieldData> > &B1_Nfaces;
  ublas::matrix<FieldData> &B1_Nvolume;

  ublas::matrix<FieldData> &B2_Nnodes,&B2_diffNnodes;
  vector<ublas::matrix<FieldData> > &B2_Nedges;
  vector<ublas::matrix<FieldData> > &B2_Nfaces;
  ublas::matrix<FieldData> &B2_Nvolume;

  mult_H1_H1(
    ublas::matrix<FieldData> &_B1_Nnodes,
    ublas::matrix<FieldData> &_B1_diffNnodes,
    vector<ublas::matrix<FieldData> > &_B1_Nedges,
    vector<ublas::matrix<FieldData> > &_B1_Nfaces,
    ublas::matrix<FieldData> &_B1_Nvolume,
    ublas::matrix<FieldData> &_B2_Nnodes,
    ublas::matrix<FieldData> &_B2_diffNnodes,
    vector<ublas::matrix<FieldData> > &_B2_Nedges,
    vector<ublas::matrix<FieldData> > &_B2_Nfaces,
    ublas::matrix<FieldData> &_B2_Nvolume):
      B1_Nnodes(_B1_Nnodes),B1_diffNnodes(_B1_diffNnodes),B1_Nedges(_B1_Nedges),B1_Nfaces(_B1_Nfaces),B1_Nvolume(_B1_Nvolume),
      B2_Nnodes(_B2_Nnodes),B2_diffNnodes(_B2_diffNnodes),B2_Nedges(_B2_Nedges),B2_Nfaces(_B2_Nfaces),B2_Nvolume(_B2_Nvolume) {}
  
  ublas::matrix<FieldData> H1_H1_nodes;

  PetscErrorCode calculate_H1_H1_nonsymetric(double *G_W);
  PetscErrorCode calculate_H1_H1_symmetric(double *G_W);
    
};

}

#endif //__CORE_FORCES_AND_SURCES_HPP

