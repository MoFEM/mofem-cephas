/** \file DataStructures.hpp

  \brief Data structures for accessing information about finite element and its
  degrees of freedom.

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
  typedef ublas::vector<double,ublas::shallow_array_adaptor<double> > VectorAdaptor;
  typedef ublas::matrix<double,ublas::row_major,ublas::shallow_array_adaptor<double> > MatrixAdaptor;

  /** \brief data on single entity
    */
  struct EntData {

    EntData(): sEnse(0),oRder(0) {};
    virtual ~EntData() {}

    /// \brief get entity sense, need to calculate shape functions with conforming approximation fields
    virtual int getSense() const { return sEnse; }

    /// \brief get approximation order
    virtual ApproximationOrder getOrder() const { return oRder; }

    /// \brief get global indices of dofs on entity
    virtual const ublas::vector<DofIdx>& getIndices() const { return iNdices; }

    /// \brief get dofs values
    virtual const ublas::vector<double>& getFieldData() const { return fieldData; }

    /// \brief get dofs data stature FEDofMoFEMEntity
    virtual const ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() const { return dOfs; }

    /** \brief get shape functions
      * this return matrix (nb. of rows is equal to nb. of Gauss pts, nb. of
      * columns is equalt to number of shape functions on this entity
      */
    virtual const ublas::matrix<double>& getN() const { return N; }

    /** \brief get derivatives of shape functions
     *
     * Matrix at rows has nb. of Gauss pts, at columns it has derivative of
     * shape functions. Columns are structured as follows, [ dN1/dx, dN1/dy,
     * dN1/dz, dN2/dx, dN2/dy, dN2/dz, ... ]
     *
     * Note that shape functions are calculated in file H1.c
     * Above description not apply for derivatives of nodal functions, since
     * derivative of nodal functions in case of simplexes, EDGES, TRIANGLES and
     * TETS are constant. So that matrix rows represents nb. of shape
     * functions, columns are derivatives. Nb. of columns depend on element
     * dimension, for EDGES is one, for TRIS is 2 and TETS is 3.
     *
     * Note that for node element this function make no sense.
     *
     */
    virtual const ublas::matrix<double>& getDiffN() const { return diffN; }

    virtual int& getSense() { return sEnse; }
    virtual ApproximationOrder& getOrder() { return oRder; }
    virtual ublas::vector<DofIdx>& getIndices() { return iNdices; }
    virtual ublas::vector<double>& getFieldData() { return fieldData; }
    virtual ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() { return dOfs; }
    virtual ublas::matrix<double>& getN() { return N; }

    /** \brief get derivatives of shape functions
     *
     */
    virtual ublas::matrix<double>& getDiffN() { return diffN; }

    /// \brief get shape functions at Gauss pts
    inline const VectorAdaptor getN(int gg) {
      int size = getN().size2();
      double *data = &getN()(gg,0);
      return VectorAdaptor(size,ublas::shallow_array_adaptor<double>(size,data));
    }

    /** \brief get derivative of shape functions at Gauss pts
      *
      * returned matrix on rows has shape functions, in column its derivatives.
      *
      * \param gg nb. of Gauss pts.
      *
      */
      inline const MatrixAdaptor getDiffN(int gg) {
        if(getN().size1() == getDiffN().size1()) {
          int size = getN().size2();
          int dim = getDiffN().size2()/size;
          double *data = &getDiffN()(gg,0);
          return MatrixAdaptor(getN().size2(),dim,ublas::shallow_array_adaptor<double>(getDiffN().size2(),data));
        } else {
          // in some cases, f.e. for derivatives of nodal shape functions at only one
          // gauss point is needed
          return MatrixAdaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<double>(getDiffN().data().size(),&getDiffN().data()[0]));
        }
      }

    /** \brief get shape functions at Gauss pts
      *
      * Note that multi field element, two different field can have different
      * approximation orders. Since we use hierarchical approximation basis,
      * shape functions are calculated once for element, using maximal
      * approximation order on given entity.
      *
      * Specifying add ional parameters, only first nb_dofs are indicated as a
      * row of shape function matrix.
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const VectorAdaptor getN(int gg,const int nb_dofs) {
      (void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
      double *data = &getN()(gg,0);
      return VectorAdaptor(nb_dofs,ublas::shallow_array_adaptor<double>(nb_dofs,data));
    }

    /** \brief get derivatives of shape functions at Gauss pts
      *
      * Note that multi field element, two different field can have different
      * approximation orders. Since we use hierarchical approximation basis,
      * shape functions are calculated once for element, using maximal
      * approximation order on given entity.
      *
      * Specifying add ional parameters, only first nb_dofs are indicated as a
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
	double *data = &getDiffN()(gg,0);
	return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<double>(dim*nb_dofs,data));
      } else {
	// in some cases, f.e. for derivatives of nodal shape functions ony one
	// gauss point is needed
	return MatrixAdaptor(getN().size1(),getN().size2(),ublas::shallow_array_adaptor<double>(getDiffN().data().size(),&getDiffN().data()[0]));
      }
    }

    /** \brief get shape functions for Hdiv space
      */
    inline const ublas::matrix<double>&  getHdivN() const { return getN(); };

    /** \brief get derivatives of shape functions for Hdiv space
      */
    inline const ublas::matrix<double>&  getDiffHdivN() const { return getDiffN(); };

    /** \brief get shape functions for Hdiv space
      */
    inline ublas::matrix<double>&  getHdivN() { return getN(); };

    /** \brief get derivatives of shape functions for Hdiv space
      */
    inline ublas::matrix<double>&  getDiffHdivN() { return getDiffN(); };

    /** \brief get Hdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getHdivN(int gg) {
      int dim = 3;
      int nb_dofs = getHdivN().size2()/dim;
      double *data = &getHdivN()(gg,0);
      return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<double>(dim*nb_dofs,data));
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getDiffHdivN(int gg) {
      int nb_dofs = getDiffHdivN().size2()/9;
      double *data = &getDiffHdivN()(gg,0);
      return MatrixAdaptor(nb_dofs,9,ublas::shallow_array_adaptor<double>(9*nb_dofs,data));
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
      *
      * \param gg nb. of Gauss point
      * \param number of of shape functions
      *
      */
    inline const MatrixAdaptor getDiffHdivN(int dof,int gg) {
      double *data = &getDiffHdivN()(gg,9*dof);
      return MatrixAdaptor(3,3,ublas::shallow_array_adaptor<double>(9,data));
    }

    friend ostream& operator<<(ostream& os,const DataForcesAndSurcesCore::EntData &e);

    private:
    int sEnse;
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<const FEDofMoFEMEntity*> dOfs;
    ublas::vector<double> fieldData;
    ublas::matrix<double> N;
    ublas::matrix<double> diffN;

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

/** \brief this class derive data form other data structure
  * \ingroup mofem_forces_and_sources
  *
  *
  * It behaves like normal data structure it is used to share information with
  * other data structures about shape functions. Dofs values, approx. order and
  * indices are not shared.
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
    const ublas::vector<double>& getFieldData() const { return fieldData; }
    const ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() const { return dOfs; }
    ublas::vector<double>& getFieldData() { return fieldData; }
    ublas::vector<const FEDofMoFEMEntity*>& getFieldDofs() { return dOfs; }
    ApproximationOrder getOrder() const { return oRder; }
    ApproximationOrder& getOrder() { return oRder; }

    inline int getSense() const { return entData.getSense(); }
    inline const ublas::matrix<double>& getN() const { return entData.getN(); }
    inline const ublas::matrix<double>& getDiffN() const { return entData.getDiffN(); }
    inline ublas::matrix<double>& getN() { return entData.getN(); }
    inline ublas::matrix<double>& getDiffN() { return entData.getDiffN(); }
    inline const ublas::matrix<double>&  getHdivN() const { return entData.getHdivN(); };
    inline ublas::matrix<double>&  getHdivN() { return entData.getHdivN(); };

    private:
    ApproximationOrder oRder;
    ublas::vector<DofIdx> iNdices;
    ublas::vector<double> fieldData;
    ublas::vector<const FEDofMoFEMEntity*> dOfs;

  };

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

};


}

#endif //__DATASTRUCTURES_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/
