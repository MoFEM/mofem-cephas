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

struct DataForcesAndSurcesCore {

  struct EntData {
    EntData(): sEnse(0),oRder(0) {};
    virtual int getSense() const { return sEnse; }
    virtual ApproximationOrder getOrder() const { return oRder; }
    virtual const ublas::vector<DofIdx>& getIndices() const { return iNdices; }
    virtual const ublas::vector<FieldData>& getFieldData() const { return fieldData; }
    virtual const ublas::matrix<FieldData>& getN() const { return N; }
    virtual const ublas::matrix<FieldData>& getDiffN() const { return diffN; }
    virtual int& getSense() { return sEnse; }
    virtual ApproximationOrder& getOrder() { return oRder; }
    virtual ublas::vector<DofIdx>& getIndices() { return iNdices; }
    virtual ublas::vector<FieldData>& getFieldData() { return fieldData; }
    virtual ublas::matrix<FieldData>& getN() { return N; }
    virtual ublas::matrix<FieldData>& getDiffN() { return diffN; }
    private:
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::matrix<FieldData> N;
    ublas::matrix<FieldData> diffN;
  };

  virtual const EntData& getNodes() const { return nOdes; }
  virtual const vector<EntData>& getEdges() const { return eDges; }
  virtual const EntData& getEdges(int ee) const { return eDges[ee]; }
  virtual const vector<EntData>& getFaces() const { return fAces; }
  virtual const EntData& getFaces(int ff) const { return fAces[ff]; }
  virtual const EntData& getVolume() const { return vOlume; }
  virtual EntData& getNodes() { return nOdes; }
  virtual vector<EntData>& getEdges() { return eDges; }
  virtual EntData& getEdges(int ee) { return eDges[ee]; }
  virtual vector<EntData>& getFaces() { return fAces; }
  virtual EntData& getFaces(int ff) { return fAces[ff]; }
  virtual EntData& getVolume() { return vOlume; }

  ublas::matrix<DofIdx> facesNodes;

  friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e);

  EntData nOdes;
  vector<EntData> eDges;
  vector<EntData> fAces;
  EntData vOlume;

};

/*struct DerivedDataForcesAndSurcesCore: public DataForcesAndSurcesCore {

 struct DerivedEntData: public EntData {
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::vector<FieldData>& getFieldData() { return iNdices; }
    ublas::vector<FieldData>& getFieldData() { return fieldData; }
    DerivedEntData(EntData &ent_data): EntData(ent_data) {}
  };

  DerivedEntData nOdes;
  vector<DerivedEntData> eDges;
  vector<DerivedEntData> fAces;
  DerivedEntData vOlume;
  virtual EntData& getNodes() { return nOdes; }
  virtual EntData& getEdges(int ee) { return eDges[ee]; }
  virtual EntData& getFaces(int ff) { return fAces[ff]; }
  virtual EntData& getVolume() { return vOlume; }

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

  friend ostream& operator<<(ostream& os,const DerivedDataForcesAndSurcesCore &e);

};*/

ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e);

struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};

  PetscErrorCode getSense(EntityType type,vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(EntityType type,vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getEdgesSense(DataForcesAndSurcesCore &data);
  PetscErrorCode getFacesSense(DataForcesAndSurcesCore &data);
  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data);
  PetscErrorCode getFacesOrder(DataForcesAndSurcesCore &data);
  PetscErrorCode getOrderVolume(DataForcesAndSurcesCore &data);

  PetscErrorCode getNodesIndices(
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,ublas::vector<int> &nodes_indices);

  PetscErrorCode getRowNodesIndices(
    DataForcesAndSurcesCore &data,
    const string &field_name);

  PetscErrorCode getColNodesIndices(
    DataForcesAndSurcesCore &data,
    const string &field_name);

  PetscErrorCode getTypeIndices(
    const string &field_name,
    FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<int> &indices);
  PetscErrorCode getTypeIndices(
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,EntityType type,
    vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getEdgeRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgeColIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesColIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getNodesFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<FieldData> &nodes_indices);
  PetscErrorCode getTypeFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<FieldData> &ent_field_data);
  PetscErrorCode getTypeFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,EntityType type,
    vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgeFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetFieldData(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getFaceNodes(DataForcesAndSurcesCore &data);

  PetscErrorCode shapeTETFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);
  PetscErrorCode shapeTRIFunctions_H1(
    DataForcesAndSurcesCore &data,
    ApproximationOrder order,const double *G_X,const double *G_Y,const int G_DIM);

};

struct DataOperator {

  virtual PetscErrorCode doWork(
    int row_side,int col_side,
    EntityType row_type,EntityType col_type,
    DataForcesAndSurcesCore::EntData &row_data,
    DataForcesAndSurcesCore::EntData &col_data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode opNH1NH1(DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data);

  virtual PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode opNH1(DataForcesAndSurcesCore &data);


};

struct OpSetJac: public DataOperator {

  ublas::matrix<double> &invJac;
  OpSetJac(ublas::matrix<double> &_invJac): invJac(_invJac) {}

  ublas::matrix<FieldData> diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

struct OpGetData: public DataOperator {

  ublas::matrix<FieldData> &data_at_GaussPt;
  ublas::matrix<FieldData> &dataGrad_at_GaussPt;
 
  const unsigned int dim;
  const ApproximationRank rank;

  OpGetData(
    ublas::matrix<FieldData> &_data_at_GaussPt,
    ublas::matrix<FieldData> &_dataGrad_at_GaussPt,
    ApproximationRank _rank,unsigned int _dim = 3): 
      data_at_GaussPt(_data_at_GaussPt),
      dataGrad_at_GaussPt(_dataGrad_at_GaussPt),
      rank(_rank),dim(_dim) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

};

struct VolumeH1H1ElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  VolumeH1H1ElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),opSetJac(invJac) {};

  DataForcesAndSurcesCore data;

  ErrorCode rval;
  PetscErrorCode ierr;
  ublas::matrix<double> invJac;
  ublas::matrix<double> dataFIELD1;
  ublas::matrix<double> dataDiffFIELD1;
  ublas::matrix<double> gaussPts;
  OpSetJac opSetJac;

  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name) {};
  };

  vector<UserDataOperator> vecUserOpNH1; 
  vector<UserDataOperator> vecUserOpNH1NH1;

  PetscErrorCode operator()();
  
};

}

#endif //__CORE_FORCES_AND_SURCES_HPP

