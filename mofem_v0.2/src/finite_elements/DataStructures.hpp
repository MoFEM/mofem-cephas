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


#ifndef __DATASTRUCTURES_HPP
#define __DATASTRUCTURES_HPP

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

  // shallow adaptor classes
  typedef ublas::vector<FieldData,ublas::shallow_array_adaptor<FieldData> > VectorAdaptor;
  typedef ublas::matrix<FieldData,ublas::row_major,ublas::shallow_array_adaptor<FieldData> > MatrixAdaptor;

  /** \brief data on single entity
    */
  struct EntData {

    EntData(): sEnse(0),oRder(0) {};
    virtual ~EntData() {}

    /// \brief get enetity sense, need to calculate shape functions with conforming approximation fields
    virtual int getSense() const { return sEnse; }

    /// \brief get approximation order
    virtual ApproximationOrder getOrder() const { return oRder; }

    /// \brief get gloabl inidces of dofs on entity
    virtual const ublas::vector<DofIdx>& getIndices() const { return iNdices; }

    /// \brief get dofs values 
    virtual const ublas::vector<FieldData>& getFieldData() const { return fieldData; }

    /// \brief get dofs data strature FEDofMoFEMEntity
    virtual const ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() const { return dOfs; }

    /** \brief get shape functions
      * this return matrix (nb. of rows is equal to nb. of Gauss pts, nb. of
      * columns is equalt to number of shape functions on this entity 
      */
    virtual const ublas::matrix<FieldData>& getN() const { return N; }

    /** \brief get direvatives of shape fucntiosn
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
    virtual const ublas::matrix<FieldData>& getDiffN() const { return diffN; }

    virtual int& getSense() { return sEnse; }
    virtual ApproximationOrder& getOrder() { return oRder; }
    virtual ublas::vector<DofIdx>& getIndices() { return iNdices; }
    virtual ublas::vector<FieldData>& getFieldData() { return fieldData; }
    virtual ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() { return dOfs; }
    virtual ublas::matrix<FieldData>& getN() { return N; }

    /** \brief get direvatives of shape functions
     *
     */
    virtual ublas::matrix<FieldData>& getDiffN() { return diffN; }

    /// \brief get shape functions at Gauss pts
    inline const VectorAdaptor getN(int gg) {
      int size = getN().size2();
      FieldData *data = &getN()(gg,0);
      return VectorAdaptor(size,ublas::shallow_array_adaptor<FieldData>(size,data));
    }

    /** \brief get direvative of shape functions at Gauss pts
      *
      * retruned matrx on rows has shape functions, in columen its direvatives.
      *
      * \param gg nb. of Gauss pts.
      *
      */
    inline const MatrixAdaptor getDiffN(int gg) {
      if(getN().size1() == getDiffN().size1()) {
	int size = getN().size2();	
	int dim = getDiffN().size2()/size;
	FieldData *data = &getDiffN()(gg,0);
	return MatrixAdaptor(getN().size2(),dim,ublas::shallow_array_adaptor<FieldData>(getDiffN().size2(),data));
      } else {
	// in some cases, f.e. for direvatives of nodal shape functions ony one
	// gauss point is needed
	return MatrixAdaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<FieldData>(getDiffN().data().size(),&getDiffN().data()[0]));
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
    inline const VectorAdaptor getN(int gg,const int nb_dofs) {
      (void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
      FieldData *data = &getN()(gg,0);
      return VectorAdaptor(nb_dofs,ublas::shallow_array_adaptor<FieldData>(nb_dofs,data));
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
    inline const MatrixAdaptor getDiffN(int gg,const int nb_dofs) {
      if(getN().size1() == getDiffN().size1()) {
	(void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
	int dim = getDiffN().size2()/getN().size2();
	FieldData *data = &getDiffN()(gg,0);
	return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<FieldData>(dim*nb_dofs,data));
      } else {
	// in some cases, f.e. for direvatives of nodal shape functions ony one
	// gauss point is needed
	return MatrixAdaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<FieldData>(getDiffN().data().size(),&getDiffN().data()[0]));
      }
    }

    /** \brief get shape functions for Hdiv space 
      */
    inline const ublas::matrix<FieldData>&  getHdivN() const { return getN(); };

    /** \brief get direvatives of shape functions for Hdiv space 
      */
    inline const ublas::matrix<FieldData>&  getDiffHdivN() const { return getDiffN(); };

    /** \brief get shape functions for Hdiv space 
      */
    inline ublas::matrix<FieldData>&  getHdivN() { return getN(); };

    /** \brief get direvatives of shape functions for Hdiv space 
      */
    inline ublas::matrix<FieldData>&  getDiffHdivN() { return getDiffN(); };

    /** \brief get Hdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getHdivN(int gg) {
      int dim = 3;
      int nb_dofs = getHdivN().size2()/dim;
      FieldData *data = &getHdivN()(gg,0);
      return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<FieldData>(dim*nb_dofs,data));
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getDiffHdivN(int gg) {
      int nb_dofs = getDiffHdivN().size2()/9;
      FieldData *data = &getDiffHdivN()(gg,0);
      return MatrixAdaptor(nb_dofs,9,ublas::shallow_array_adaptor<FieldData>(9*nb_dofs,data));
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getDiffHdivN(int dof,int gg) {
      FieldData *data = &getDiffHdivN()(gg,9*dof);
      return MatrixAdaptor(3,3,ublas::shallow_array_adaptor<FieldData>(9,data));
    }

    friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e);

    private:
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<const FEDofMoFEMEntity*> dOfs;
    ublas::vector<FieldData> fieldData;
    ublas::matrix<FieldData> N;
    ublas::matrix<FieldData> diffN;

  };

  ublas::matrix<DofIdx> facesNodes; 			///< nodes on finite element faces
  bitset<LASTSPACE> spacesOnEntities[MBMAXTYPE]; 	///< spaces on entity types
  boost::ptr_vector<EntData> dataOnEntities[MBMAXTYPE]; ///< data on nodes, shape function, dofs values, etc.

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
    const ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() const { return dOfs; }
    ublas::vector<FieldData>& getFieldData() { return fieldData; }
    ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() { return dOfs; }
    ApproximationOrder getOrder() const { return oRder; }
    ApproximationOrder& getOrder() { return oRder; }

    inline int getSense() const { return entData.getSense(); }
    inline const ublas::matrix<FieldData>& getN() const { return entData.getN(); }
    inline const ublas::matrix<FieldData>& getDiffN() const { return entData.getDiffN(); }
    inline ublas::matrix<FieldData>& getN() { return entData.getN(); }
    inline ublas::matrix<FieldData>& getDiffN() { return entData.getDiffN(); }
    inline const ublas::matrix<FieldData>&  getHdivN() const { return entData.getHdivN(); };
    inline ublas::matrix<FieldData>&  getHdivN() { return entData.getHdivN(); };

    private:
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<FieldData> fieldData;
    ublas::vector<const FEDofMoFEMEntity*> dOfs;

  };

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

};

/** \brief base clas to get information form mofem and but them into DataForcesAndSurcesCore
  * \ingroup mofem_forces_and_sources
  */
struct ForcesAndSurcesCore: public FEMethod {

  FieldInterface& mField;
  ForcesAndSurcesCore(FieldInterface& _mField): 
    mField(_mField) {};
  virtual ~ForcesAndSurcesCore() {}

  PetscErrorCode getSense(EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(const EntityType type,const FieldSpace space,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);
  PetscErrorCode getOrder(const string &field_name,const EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getEdgesSense(DataForcesAndSurcesCore &data);
  PetscErrorCode getTrisSense(DataForcesAndSurcesCore &data);

  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getTrisOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getTetsOrder(DataForcesAndSurcesCore &data,const FieldSpace space);
  PetscErrorCode getEdgesOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisOrder(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsOrder(DataForcesAndSurcesCore &data,const string &field_name);

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

  PetscErrorCode getEdgesRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesColIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisColIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsRowIndices(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsColIndices(DataForcesAndSurcesCore &data,const string &field_name);

  //data
  PetscErrorCode getNodesFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<FieldData> &nodes_data);
  PetscErrorCode getTypeFieldData(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<FieldData> &ent_field_data);
  PetscErrorCode getTypeFieldData(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getNodesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldData(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldData(DataForcesAndSurcesCore &data,const string &field_name);

  //dofs
  PetscErrorCode getNodesFieldDofs(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,ublas::vector<const FEDofMoFEMEntity*> &nodes_dofs);
  PetscErrorCode getTypeFieldDofs(
    const string &field_name,
    FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,int side_number,ublas::vector<const FEDofMoFEMEntity*> &ent_field_dofs);
  PetscErrorCode getTypeFieldDofs(
    const string &field_name,FEDofMoFEMEntity_multiIndex &dofs,
    EntityType type,boost::ptr_vector<DataForcesAndSurcesCore::EntData> &data);

  PetscErrorCode getNodesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getEdgesFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTrisFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);
  PetscErrorCode getTetsFieldDofs(DataForcesAndSurcesCore &data,const string &field_name);

  PetscErrorCode getFaceNodes(DataForcesAndSurcesCore &data);
  PetscErrorCode getSpacesOnEntities(DataForcesAndSurcesCore &data);

  /** \brief computes approximation functions for tetrahedral and H1 space
    */
  PetscErrorCode shapeTETFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);

  /** \brief computes approximation functions for tetrahedral and L2 space
    */
  PetscErrorCode shapeTETFunctions_L2(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);


  ublas::matrix<ublas::matrix<FieldData> > N_face_edge;
  ublas::vector<ublas::matrix<FieldData> > N_face_bubble;
  ublas::vector<ublas::matrix<FieldData> > N_volume_edge;
  ublas::vector<ublas::matrix<FieldData> > N_volume_face;
  ublas::matrix<FieldData> N_volume_bubble;

  ublas::matrix<ublas::matrix<FieldData> > diffN_face_edge;
  ublas::vector<ublas::matrix<FieldData> > diffN_face_bubble;
  ublas::vector<ublas::matrix<FieldData> > diffN_volume_edge;
  ublas::vector<ublas::matrix<FieldData> > diffN_volume_face;
  ublas::matrix<FieldData> diffN_volume_bubble;


  /** \brief computes approximation functions for tetrahedral and H1 space
    */
  PetscErrorCode shapeTETFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const double *G_Z,const int G_DIM);


  /** \brief computes approximation functions for triangle and H1 space
    */
  PetscErrorCode shapeTRIFunctions_H1(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);


  /** \brief computes approximation functions for triangle and H1 space
    */
  PetscErrorCode shapeTRIFunctions_Hdiv(
    DataForcesAndSurcesCore &data,
    const double *G_X,const double *G_Y,const int G_DIM);

  /** \brief computes approximation functions for edge and H1 space
    */
  PetscErrorCode shapeEDGEFunctions_H1(
    DataForcesAndSurcesCore &data,const double *G_X,const int G_DIM);

  /** \brief it is used to calculate nb. of Gauss integartion points
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

  virtual PetscErrorCode setGaussPts(int order) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"sorry, not implemented");
    PetscFunctionReturn(0);
  }

};

}

#endif //__DATASTRUCTURES_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/


