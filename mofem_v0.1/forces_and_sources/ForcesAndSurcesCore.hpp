/** \file ForcesAndSurcesCore.hpp 
 * \brief Forces and sources data structures
 *
 * It is set of objects to implement finite elements, in particular it is used
 * to implement source and force therms on right hand side. 
 *
*/

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
#include <boost/ptr_container/ptr_map.hpp>

using namespace boost::numeric;

namespace MoFEM {

/** \brief data structure for finite element entity
  * \ingroup mofem_forces_and_sources
  *
  * It keeps that about indices of degrees of freedom, dofs data, shape
  * functions, entity side number, type of entities, approximation order, etc.
  *
  */
struct DataForcesAndSurcesCore {


  /** \brief data on single entity
    */
  struct EntData {

    EntData(): sEnse(0),oRder(0) {};
    virtual ~EntData() {}

    /// \brief get enetity sense, need to calulate shape functions with conforming approximation fields
    virtual int getSense() const { return sEnse; }

    /// \brief get approximation order
    virtual ApproximationOrder getOrder() const { return oRder; }

    /// \brief get gloabl inidces of dofs on entity
    virtual const ublas::vector<DofIdx>& getIndices() const { return iNdices; }

    /// \brief get dofs values 
    virtual const ublas::vector<FieldData>& getFieldData() const { return fieldData; }

    /// \brief get shape functions
    virtual const ublas::matrix<FieldData>& getN() const { return N; }

    /// \brief get direvatives of shape fucntiosn
    virtual const ublas::matrix<FieldData>& getDiffN() const { return diffN; }

    virtual int& getSense() { return sEnse; }
    virtual ApproximationOrder& getOrder() { return oRder; }
    virtual ublas::vector<DofIdx>& getIndices() { return iNdices; }
    virtual ublas::vector<FieldData>& getFieldData() { return fieldData; }
    virtual ublas::matrix<FieldData>& getN() { return N; }

    /** \brief get direvatives of shape functions
     *
     * Matrix at rows has nb. of Gauss pts, at columns it has direvative of
     * shape functions. Colummns are orgasised as follows, [ dN1/dx, dN1/dy,
     * dN1/dz, dN2/dx, dN2/dy, dN2/dz, ... ]
     *
     * Note that shape functions are calculated in file H1.c
     * Above description not apply for direvatives of nodal functions, since
     * direvative of nodal functions in case of simplexes, EDGES, TRIANGLES and
     * TETS are constant. So that matrix rows represents nb. of shape
     * functions, columns are direvatives. Nb. of columns depend on element
     * dimension, for EDGES is one, for TRIS is 2 and TETS is 3. 
     *
     * Note that for node element this function make no sense.
     *
     */
    virtual ublas::matrix<FieldData>& getDiffN() { return diffN; }

    // shallow adaptor classes
    typedef ublas::vector<FieldData,ublas::shallow_array_adaptor<FieldData> > vector_adaptor;
    typedef ublas::matrix<FieldData,ublas::row_major,ublas::shallow_array_adaptor<FieldData> > matrix_adaptor;

    /// \brief get shape functions at Gauss pts
    inline const vector_adaptor getN(int gg) {
      int size = getN().size2();
      FieldData *data = &getN()(gg,0);
      return vector_adaptor(size,ublas::shallow_array_adaptor<FieldData>(size,data));
    }

    /** \brief get direvative of shape functions at Gauss pts
      *
      * retruned matrx on rows has shape functions, in columen its direvatives.
      *
      * \param gg nb. of Gauss pts.
      *
      */
    inline const matrix_adaptor getDiffN(int gg) {
      if(getN().size1() == getDiffN().size1()) {
	int size = getN().size2();	
	int dim = getDiffN().size2()/size;
	FieldData *data = &getDiffN()(gg,0);
	return matrix_adaptor(getN().size2(),dim,ublas::shallow_array_adaptor<FieldData>(getDiffN().size2(),data));
      } else {
	// in some cases, f.e. for direvatives of nodal shape functions ony one
	// gauss point is needed
	return matrix_adaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<FieldData>(getDiffN().data().size(),&getDiffN().data()[0]));
      }
    }

    /** \brief get shape functions at Gauss pts
      *
      * Note that multi field element, two diffren field can have different
      * approximation orders. Since we use hierarhical approximation basis,
      * shape functions are calulared once for element, using maximal
      * apprimation order on given entity.
      *
      * Specifing addional parameters, only firsty nb_dofs are indicated as a
      * row of shape function matrix.
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const vector_adaptor getN(int gg,const int nb_dofs) {
      (void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
      FieldData *data = &getN()(gg,0);
      return vector_adaptor(nb_dofs,ublas::shallow_array_adaptor<FieldData>(nb_dofs,data));
    }

    /** \brief get derivatives of shape functions at Gauss pts
      *
      * Note that multi field element, two diffren field can have different
      * approximation orders. Since we use hierarhical approximation basis,
      * shape functions are calulared once for element, using maximal
      * apprimation order on given entity.
      *
      * Specifing addional parameters, only firsty nb_dofs are indicated as a
      * row of shape function derivative matrix.
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const matrix_adaptor getDiffN(int gg,const int nb_dofs) {
      if(getN().size1() == getDiffN().size1()) {
	(void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
	int dim = getDiffN().size2()/getN().size2();
	FieldData *data = &getDiffN()(gg,0);
	return matrix_adaptor(nb_dofs,dim,ublas::shallow_array_adaptor<FieldData>(dim*nb_dofs,data));
      } else {
	// in some cases, f.e. for direvatives of nodal shape functions ony one
	// gauss point is needed
	return matrix_adaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<FieldData>(getDiffN().data().size(),&getDiffN().data()[0]));

      }
    }

    friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e);

    private:
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::matrix<FieldData> N;
    ublas::matrix<FieldData> diffN;

  };

  ublas::matrix<DofIdx> facesNodes; ///< nodes on finite element faces

  boost::ptr_vector<EntData> nOdes; ///< data on nodes, shape function, dofs values, etc.
  boost::ptr_vector<EntData> eDges; ///< data on edges, shape function, dofs values, etc.
  boost::ptr_vector<EntData> fAces; ///< data on faces, shape function, dofs values, etc.
  boost::ptr_vector<EntData> vOlumes; ///< data on volume, shape function, dofs values, etc.

  DataForcesAndSurcesCore(EntityType type);
  virtual ~DataForcesAndSurcesCore() {}

  friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore &e);

  protected:
  DataForcesAndSurcesCore() {}

};

/** \brief this class derive data form other dats strucrure
  * \ingroup mofem_forces_and_sources
  *
  *
  * It behavies like normal data struture it is used to share infromation with
  * other data strutures abot shape functons. Dofs values, approx. order and
  * incices are not shared.
  *
  * shape functions, senses are shared with other data structure.
  *
  */
struct DerivedDataForcesAndSurcesCore: public DataForcesAndSurcesCore  {

  struct DerivedEntData: public DataForcesAndSurcesCore::EntData {
    DataForcesAndSurcesCore::EntData &entData;
    DerivedEntData(DataForcesAndSurcesCore::EntData &ent_data): 
      entData(ent_data),oRder(0) {}
    const ublas::vector<DofIdx>& getIndices() const { return iNdices; }
    ublas::vector<DofIdx>& getIndices() { return iNdices; }
    const ublas::vector<FieldData>& getFieldData() const { return fieldData; }
    ublas::vector<FieldData>& getFieldData() { return fieldData; }
    ApproximationOrder getOrder() const { return oRder; }
    ApproximationOrder& getOrder() { return oRder; }

    int getSense() const { return entData.getSense(); }
    const ublas::matrix<FieldData>& getN() const { return entData.getN(); }
    const ublas::matrix<FieldData>& getDiffN() const { return entData.getDiffN(); }
    int& getSense() { return entData.getSense(); }
    ublas::matrix<FieldData>& getN() { return entData.getN(); }
    ublas::matrix<FieldData>& getDiffN() { return entData.getDiffN(); }

    private:
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;

  };

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

};

/** \brief base clas to get information form mofem and but them into DataForcesAndSurcesCore
  * \ingroup mofem_forces_and_sources
  */
struct ForcesAndSurcesCore: public FieldInterface::FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};
  virtual ~ForcesAndSurcesCore() {}

  PetscErrorCode getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(const string &field_name,EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getEdgesSense(DataForcesAndSurcesCore &data);
  PetscErrorCode getFacesSense(DataForcesAndSurcesCore &data);

  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data);
  PetscErrorCode getFacesOrder(DataForcesAndSurcesCore &data);
  PetscErrorCode getVolumesOrder(DataForcesAndSurcesCore &data);
  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getFacesOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getVolumesOrder(DataForcesAndSurcesCore &data,const string &field_name);

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
  PetscErrorCode getVolumesFieldData(DataForcesAndSurcesCore &data,const string &field_name);

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
  PetscErrorCode shapeEDGEFunctions_H1(
    DataForcesAndSurcesCore &data,const double *G_X,const int G_DIM);

};

/** \brief base operator to do operations at Gauss Pt. leve
  * \ingroup mofem_forces_and_sources
  */
struct DataOperator {

  /** \brief operator for linear form, usaully to calulate values on right hand side
    */
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

  /** \brief operator for linear form, usaully to calulate values on left hand side
    */
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

/// \brief operator on Gauss pts level, calulates inverse of Jacobian
struct OpSetInvJac: public DataOperator {

  ublas::matrix<double> &invJac;
  OpSetInvJac(ublas::matrix<double> &_invJac): invJac(_invJac) {}

  ublas::matrix<FieldData> diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

};

/** \brief operator on Gauss pts level, calulates inverse of Jacobian for higher order geometry approximation
  */
struct OpSetHoInvJac: public DataOperator {

  ublas::matrix<double> &invHoJac;
  OpSetHoInvJac(ublas::matrix<double> &_invHoJac): invHoJac(_invHoJac) {}

  ublas::matrix<FieldData> diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data);
 
};

/** \brief operator to calculate function values and its gradients at Gauss points
  * \ingroup mofem_forces_and_sources
  */
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
      dim(_dim),rank(_rank) {}

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data);

};

/** \brief Tet finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from TetElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpNH1 and
 * vecUserOpNH1NH1. 
 *
 */
struct TetElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  OpSetInvJac opSetInvJac;

  string meshPositionsFieldName;
  ublas::matrix<FieldData> hoCoordsAtGaussPtsPts;
  ublas::matrix<FieldData> hoGaussPtsInvJac;
  ublas::vector<FieldData> hoGaussPtsDetJac;
  OpGetData opHOatGaussPoints; ///< higher order geometry data at Gauss pts
  OpSetHoInvJac opSetHoInvJac;

  TetElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBTET),
    derivedData(data),opSetInvJac(invJac),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHOatGaussPoints(hoCoordsAtGaussPtsPts,hoGaussPtsInvJac,3,3),
    opSetHoInvJac(hoGaussPtsInvJac) {};
    
  virtual ~TetElementForcesAndSurcesCore() {}

  ErrorCode rval;
  PetscErrorCode ierr;
  double vOlume;
  ublas::vector<double> coords;
  ublas::matrix<double> invJac;
  ublas::matrix<double> gaussPts;
  ublas::matrix<double> coordsAtGaussPts;

  /** \brief it is used to calulate nb. of Gauss integartion points
   *
   * for more details pleas look 
   *   Reference:
   *
   * Albert Nijenhuis, Herbert Wilf,
   * Combinatorial Algorithms for Computers and Calculators,
   * Second Edition,
   * Academic Press, 1978,
   * ISBN: 0-12-519260-6,
   * LC: QA164.N54.
   *
   * More details about algorithm 
   * http://people.sc.fsu.edu/~jburkardt/cpp_src/gm_rule/gm_rule.html
  **/
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
    virtual ~UserDataOperator() {}
    inline double getVolume() { return ptrFE->vOlume; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline ublas::matrix<double>& getHoCoordsAtGaussPtsPts() { return ptrFE->hoCoordsAtGaussPtsPts; }
    inline ublas::matrix<double>& getHoGaussPtsInvJac() { return ptrFE->hoGaussPtsInvJac; }
    inline ublas::vector<double>& getHoGaussPtsDetJac() { return ptrFE->hoGaussPtsDetJac; }
    inline const FieldInterface::FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(TetElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    TetElementForcesAndSurcesCore *ptrFE; 

  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calulate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calulate matrices or other quantities on mesh.
   */
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

/** \brief calulate normals at Gauss points of triangle element
  * \ingroup mofem_forces_and_sources
  */
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

/** \brief Tri finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from TriElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpNH1 and
 * vecUserOpNH1NH1. 
 *
 */
struct TriElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  string meshPositionsFieldName;

  ublas::matrix<FieldData> nOrmals_at_GaussPt;
  ublas::matrix<FieldData> tAngent1_at_GaussPt;
  ublas::matrix<FieldData> tAngent2_at_GaussPt;
  OpGetNormals opHONormals;

  TriElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBTRI),derivedData(data),
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
    virtual ~UserDataOperator() {}
    inline double getArea() { return ptrFE->aRea; }

    /** \bried get triangle normal
     */
    inline ublas::vector<double>& getNormal() { return ptrFE->normal; }

    /** \bried get triangle coords
     */
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }

    /** \bried get triangle Gauss pts.
     */
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }

    /** \bried get coordinates at Gauss pts.
     */
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \bried if higher order geometry return normals at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getNormals_at_GaussPt() { return ptrFE->nOrmals_at_GaussPt; }

    /** \bried if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<ublas::matrix<double> > getNormals_at_GaussPt(const int gg) { 
      return ublas::matrix_row<ublas::matrix<double> >(ptrFE->nOrmals_at_GaussPt,gg); 
    }

    /** \bried if higher order geometry return tangent vetor to triangle at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getTangent1_at_GaussPt() { return ptrFE->tAngent1_at_GaussPt; }

    /** \bried if higher order geometry return tangent vetor to triangle at Gauss pts.
     */
    inline ublas::matrix<FieldData>& getTangent2_at_GaussPt() { return ptrFE->tAngent2_at_GaussPt; }

    /** \bried return pointer to triangle finite element object 
     */
    inline const TriElementForcesAndSurcesCore* getTriElementForcesAndSurcesCore() { return ptrFE; }

    /** \bried return pointer to FEMthod object
     */
    inline const FieldInterface::FEMethod* getFEMethod() { return ptrFE; }

    /** \bried return pointer to NumeredMoFEMFiniteElement 
     */
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };

    PetscErrorCode setPtrFE(TriElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    TriElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calulate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calulate matrices or other quantities on mesh.
   */
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

/** \brief Edge finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from EdgeElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpNH1 and
 * vecUserOpNH1NH1. 
 *
 */
struct EdgeElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;

  EdgeElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBEDGE),derivedData(data) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  double lEngth;;
  ublas::vector<double> dIrection;
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
    virtual ~UserDataOperator() {}
    inline double getLength() { return ptrFE->lEngth; }
    inline ublas::vector<double>& getDirection() { return ptrFE->dIrection; }
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline ublas::matrix<double>& getGaussPts() { return ptrFE->gaussPts; }
    inline ublas::matrix<double>& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }
    inline const FieldInterface::FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(EdgeElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    EdgeElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calulate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calulate matrices or other quantities on mesh.
   */
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

/** \brief Vertex finite element  
 * \ingroup mofem_forces_and_sources
 *
 * User is implementing own operator at Guass piint level, by own object
 * derived from VertexElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 * number of operator added pushing objects to vecUserOpNH1 and
 * vecUserOpNH1NH1. 
 *
 */
struct VertexElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  DataForcesAndSurcesCore data;
  DerivedDataForcesAndSurcesCore derivedData;
  string meshPositionsFieldName;

  VertexElementForcesAndSurcesCore(FieldInterface &_mField):
    ForcesAndSurcesCore(_mField),data(MBVERTEX),derivedData(data) {};

  ErrorCode rval;
  PetscErrorCode ierr;
  ublas::vector<double> coords;

  struct UserDataOperator: public DataOperator {
    string row_field_name;
    string col_field_name;
    UserDataOperator(
      const string &_field_name):
	row_field_name(_field_name),col_field_name(_field_name),ptrFE(NULL) {};
    UserDataOperator(
      const string &_row_field_name,const string &_col_field_name):
	row_field_name(_row_field_name),col_field_name(_col_field_name),ptrFE(NULL) {};
    virtual ~UserDataOperator() {}
    inline ublas::vector<double>& getCoords() { return ptrFE->coords; }
    inline const FieldInterface::FEMethod* getFEMethod() { return ptrFE; }
    inline const NumeredMoFEMFiniteElement* getMoFEMFEPtr() { return ptrFE->fePtr; };
    PetscErrorCode setPtrFE(VertexElementForcesAndSurcesCore *ptr) { 
      PetscFunctionBegin;
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }
    private:
    VertexElementForcesAndSurcesCore *ptrFE; 
  };

  boost::ptr_vector<UserDataOperator> vecUserOpNH1; 
  boost::ptr_vector<UserDataOperator> vecUserOpNH1NH1;

  /** \brief Use to push back operator for right hand side
   * It can be ussed to calulate nodal forces or other quantities on the mesh.
   */
  boost::ptr_vector<UserDataOperator>& get_op_to_do_Rhs() { return vecUserOpNH1; }

  /** \brief Use to push back operator for left hand side
   * It can be ussed to calulate matrices or other quantities on mesh.
   */
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

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/




