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
#include <boost/ptr_container/ptr_vector.hpp>

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

    friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e);

    private:
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::matrix<FieldData> N;
    ublas::matrix<FieldData> diffN;
    
  };

  ublas::matrix<DofIdx> facesNodes;

  boost::ptr_vector<EntData> nOdes;
  boost::ptr_vector<EntData> eDges;
  boost::ptr_vector<EntData> fAces;
  boost::ptr_vector<EntData> vOlumes;

  DataForcesAndSurcesCore(EntityType type);

  friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e);

  protected:
  DataForcesAndSurcesCore() {}

};

struct DerivedDataForcesAndSurcesCore: public DataForcesAndSurcesCore  {

  struct DerivedEntData: public DataForcesAndSurcesCore::EntData {
    DataForcesAndSurcesCore::EntData &entData;
    DerivedEntData(DataForcesAndSurcesCore::EntData &ent_data): entData(ent_data)  {}
    const ublas::vector<DofIdx>& getIndices() const { return iNdices; }
    ublas::vector<DofIdx>& getIndices() { return iNdices; }
    const ublas::vector<FieldData>& getFieldData() const { return fieldData; }
    ublas::vector<FieldData>& getFieldData() { return fieldData; }

    int getSense() const { return entData.getSense(); }
    ApproximationOrder getOrder() const { return entData.getOrder(); }
    const ublas::matrix<FieldData>& getN() const { return entData.getN(); }
    const ublas::matrix<FieldData>& getDiffN() const { return entData.getDiffN(); }
    int& getSense() { return entData.getSense(); }
    ApproximationOrder& getOrder() { return entData.getOrder(); }
    ublas::matrix<FieldData>& getN() { return entData.getN(); }
    ublas::matrix<FieldData>& getDiffN() { return entData.getDiffN(); }

    private:
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;

  };

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

};

struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};

  PetscErrorCode getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

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
    const string &field_name,FENumeredDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

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
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

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
    const double *G_X,const double *G_Y,const int G_DIM);

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
  PetscErrorCode opSymmetric(DataForcesAndSurcesCore &row_data,DataForcesAndSurcesCore &col_data);

  virtual PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    PetscFunctionReturn(0);
  }
  PetscErrorCode op(DataForcesAndSurcesCore &data);


};

struct OpSetInvJac: public DataOperator {

  ublas::matrix<double> &invJac;
  OpSetInvJac(ublas::matrix<double> &_invJac): invJac(_invJac) {}

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

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derived_data;

  VolumeH1H1ElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),opSetInvJac(invJac),
    data(MBTET),derived_data(data) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  double vOlume;
  ublas::vector<double> coords;
  ublas::matrix<double> invJac;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;
  OpSetInvJac opSetInvJac;

  virtual int getRule(int order) { return order; };

  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),ptrFE(NULL) {};
    inline double getVolume() { return ptrFE->vOlume; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fe_ptr; };
    PetscErrorCode setPtrFE(VolumeH1H1ElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    VolumeH1H1ElementForcesAndSurcesCore *ptrFE; 

  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpNH1NH1; }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  
};

struct OpGetNormals: public DataOperator {

  ublas::matrix<FieldData> &nOrmals_at_GaussPt;
  ublas::matrix<FieldData> &tAngent1_at_GaussPt;
  ublas::matrix<FieldData> &tAngent2_at_GaussPt;

  OpGetNormals(
    ublas::matrix<FieldData> &_nOrmals_at_GaussPt,
    ublas::matrix<FieldData> &_tAngent1_at_GaussPt,
    ublas::matrix<FieldData> &_tAngent2_at_GaussPt): 
    nOrmals_at_GaussPt(_nOrmals_at_GaussPt),
    tAngent1_at_GaussPt(_tAngent1_at_GaussPt),
    tAngent2_at_GaussPt(_tAngent2_at_GaussPt) {}

  ublas::matrix<FieldData> sPin;
  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

  PetscErrorCode calculateNormals();

};

struct TriangleH1H1ElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derived_data;
  string meshPositionsFieldName;

  ublas::matrix<FieldData> nOrmals_at_GaussPt;
  ublas::matrix<FieldData> tAngent1_at_GaussPt;
  ublas::matrix<FieldData> tAngent2_at_GaussPt;
  OpGetNormals opHONormals;

  TriangleH1H1ElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBTRI),derived_data(data),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHONormals(nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  double aRea;;
  ublas::vector<double> normal;
  ublas::vector<double> coords;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  virtual int getRule(int order) { return order; };

  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),ptrFE(NULL) {};
    inline double getArea() { return ptrFE->aRea; }
    inline ublas::vector<double>& getNormal() { return ptrFE->normal; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline ublas::matrix<FieldData>& getNormals_at_GaussPt() { return ptrFE->nOrmals_at_GaussPt; }
    inline ublas::matrix<FieldData>& getTangent1_at_GaussPt() { return ptrFE->tAngent1_at_GaussPt; }
    inline ublas::matrix<FieldData>& getTangent2_at_GaussPt() { return ptrFE->tAngent2_at_GaussPt; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fe_ptr; };
    PetscErrorCode setPtrFE(TriangleH1H1ElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    TriangleH1H1ElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Lhs() { return vecUserOpNH1NH1; }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()();
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

}

#endif //__CORE_FORCES_AND_SURCES_HPP

