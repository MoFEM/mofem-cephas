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

typedef ublas::unbounded_array<const FEDofEntity*,std::allocator<const FEDofEntity*> > DofsAllacator;
typedef ublas::vector<const FEDofEntity*,DofsAllacator > VectorDofs;

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

    EntData();
    virtual ~EntData();

    /// \brief get entity sense, need to calculate shape functions with conforming approximation fields
    virtual int getSense() const { return sEnse; }

    /// \brief get approximation order
    inline ApproximationOrder getOrder() const { return oRder; }

    /// \brief get global indices of dofs on entity
    inline const VectorInt& getIndices() const { return iNdices; }

    /// \brief get global indices of dofs on entity up to given order
    inline const VectorIntAdaptor getIndicesUpToOrder(int order) {
      unsigned int size = 0;
      if(iNdices.size()) {
        size = dOfs[0]->get_order_nb_dofs(order)*dOfs[0]->get_nb_of_coeffs();
        size = size < iNdices.size() ? size : iNdices.size();
      }
      int *data = &*iNdices.data().begin();
      return VectorIntAdaptor(size,ublas::shallow_array_adaptor<int>(size,data));
    }

    /// \brief get local indices of dofs on entity
    inline const VectorInt& getLocalIndices() const { return localIndices; }

    /// \brief get local indices of dofs on entity up to given order
    inline const VectorIntAdaptor getLocalIndicesUpToOrder(int order) {
      unsigned int size = 0;
      if(localIndices.size()) {
        size = dOfs[0]->get_order_nb_dofs(order)*dOfs[0]->get_nb_of_coeffs();
        size = size < localIndices.size() ? size : localIndices.size();
      }
      int *data = &*localIndices.data().begin();
      return VectorIntAdaptor(size,ublas::shallow_array_adaptor<int>(size,data));
    }

    /// \brief get dofs values
    inline const VectorDouble& getFieldData() const { return fieldData; }

    /// \brief get dofs values up to given order
    inline const VectorAdaptor getFieldDataUpToOrder(int order) {
      unsigned int size = 0;
      if(fieldData.size()) {
        size = dOfs[0]->get_order_nb_dofs(order)*dOfs[0]->get_nb_of_coeffs();
        size = size < fieldData.size() ? size : fieldData.size();
      }
      double *data = &*fieldData.data().begin();
      return VectorAdaptor(size,ublas::shallow_array_adaptor<double>(size,data));
    }

    /// \brief get dofs data stature FEDofEntity
    inline const VectorDofs& getFieldDofs() const { return dOfs; }

    inline int& getSense() { return sEnse; }
    inline ApproximationOrder& getDataOrder() { return oRder; }
    inline VectorInt& getIndices() { return iNdices; }
    inline VectorInt& getLocalIndices() { return localIndices; }
    inline VectorDouble& getFieldData() { return fieldData; }
    inline VectorDofs& getFieldDofs() { return dOfs; }

    /**
     * \brief Get approximation base
     * @return Approximation base
     */
    inline FieldApproximationBase& getBase() { return bAse; }

    /**
     * \brief Get field space
     * @return Field space
     */
    inline FieldSpace& getSpace() { return sPace; }

    /**
     * Get shared pointer to base shape functions
     */
    virtual boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) {
      return N[base];
    }

    /**
     * Get shared pointer to base shape functions
     */
    virtual const boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) const {
      return N[base];
    }

    /**
    * Get shared pointer to derivatives of base shape functions
    */
    virtual boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) {
      return diffN[base];
    }

    /**
    * Get shared pointer to derivatives of base shape functions
    */
    virtual const boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) const {
      return diffN[base];
    }

    /** \brief get shape functions
    * this return matrix (nb. of rows is equal to nb. of Gauss pts, nb. of
    * columns is equal to number of shape functions on this entity
    */
    virtual const MatrixDouble& getN(const FieldApproximationBase base) const {
      return *(getNSharedPtr(base));
    }

    inline const MatrixDouble& getN() const { return getN(bAse); }

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
    virtual const MatrixDouble& getDiffN(const FieldApproximationBase base) const {
      return *(getDiffNSharedPtr(base));
    }

    inline const MatrixDouble& getDiffN() const { return getDiffN(bAse); }

    /**
     * \brief Get shape functions
     * @param  base Approximation base
     * @return      Error code
     */
    inline MatrixDouble& getN(const FieldApproximationBase base) { return *(getNSharedPtr(base)); }

    /**
     * \brief Get shape functions
     *
     * It assumed that approximation base for given field is known and stored in this data structure
     *
     * @return Error code
     */
     inline MatrixDouble& getN() { return getN(bAse); }

    /**
     * \brief Get derivatives of shape functions
     * @param  base Approximation base
     * @return      Error code
     */
     inline MatrixDouble& getDiffN(const FieldApproximationBase base) { return *(getDiffNSharedPtr(base)); }

    /**
     * \brief Get derivatives of shape functions
     *
     * It assumed that approximation base for given field is known and stored in this data structure
     *
     * @return Error code
     */
     inline MatrixDouble& getDiffN() { return getDiffN(bAse); }

    /// \brief get shape functions at Gauss pts
    inline const VectorAdaptor getN(const FieldApproximationBase base,const int gg) {
      int size = getN(base).size2();
      double *data = &getN(base)(gg,0);
      return VectorAdaptor(size,ublas::shallow_array_adaptor<double>(size,data));
    }

    /// \brief get shape functions at Gauss pts
    inline const VectorAdaptor getN(const int gg) { return getN(bAse,gg); }

    /** \brief get derivative of shape functions at Gauss pts

    * returned matrix on rows has shape functions, in column its derivatives.
    *
    * \param base Approximation base
    * \param gg Nb. of Gauss pts.
    *
    */
    inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,const int gg) {
      if(getN(base).size1() == getDiffN(base).size1()) {
        int size = getN(base).size2();
        int dim = getDiffN(base).size2()/size;
        double *data = &getDiffN(base)(gg,0);
        return MatrixAdaptor(
          getN(base).size2(),dim,ublas::shallow_array_adaptor<double>(getDiffN(base).size2(),data)
        );
      } else {
        // in some cases, f.e. for derivatives of nodal shape functions at only one
        // gauss point is needed
        return MatrixAdaptor(
          getN(base).size1(),
          getN(base).size2(),
          ublas::shallow_array_adaptor<double>(
            getDiffN(base).data().size(),&getDiffN(base).data()[0]
          )
        );
      }
    }

    /** \brief get derivative of shape functions at Gauss pts

    * returned matrix on rows has shape functions, in column its derivatives.
    *
    * \param gg nb. of Gauss pts.
    *
    */
    inline const MatrixAdaptor getDiffN(const int gg) { return getDiffN(bAse,gg); }

    /** \brief get shape functions at Gauss pts

    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * shape functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param base Approximation base
    * \param gg number of Gauss point
    * \param nb_dofs number of of shape functions returned

    */
    inline const VectorAdaptor getN(const FieldApproximationBase base,const int gg,const int nb_dofs) {
      (void)getN()(gg,nb_dofs-1); // throw error if nb_dofs is to big
      double *data = &getN(base)(gg,0);
      return VectorAdaptor(nb_dofs,ublas::shallow_array_adaptor<double>(nb_dofs,data));
    }

    /** \brief get shape functions at Gauss pts

    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * shape functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param gg number of Gauss point
    * \param nb_dofs number of of shape functions returned

    */
    inline const VectorAdaptor getN(const int gg,const int nb_dofs) {
      return getN(bAse,gg,nb_dofs);
    }

    /** \brief get derivatives of shape functions at Gauss pts
    *
    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * shape functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param nb_dofs number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,const int gg,const int nb_dofs) {
      if(getN(base).size1() == getDiffN(base).size1()) {
        (void)getN(base)(gg,nb_dofs-1); // throw error if nb_dofs is to big
        int dim = getDiffN(base).size2()/getN(base).size2();
        double *data = &getDiffN(base)(gg,0);
        return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<double>(dim*nb_dofs,data));
      } else {
        // in some cases, f.e. for derivatives of nodal shape functions only one
        // gauss point is needed
        return MatrixAdaptor(
          getN(base).size1(),
          getN(base).size2(),
          ublas::shallow_array_adaptor<double>(
            getDiffN(base).data().size(),
            &getDiffN(base).data()[0]
          )
        );
      }
    }

    /** \brief get derivatives of shape functions at Gauss pts
    *
    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * shape functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param gg nb. of Gauss point
    * \param nb_dofs number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffN(const int gg,const int nb_dofs) {
      return getDiffN(bAse,gg,nb_dofs);
    }


    inline const MatrixDouble& getHdivN(const FieldApproximationBase base) const { return getN(base); };
    inline const MatrixDouble& getDiffHdivN(const FieldApproximationBase base) const { return getDiffN(base); };
    inline MatrixDouble& getHdivN(const FieldApproximationBase base) { return getN(base); };
    inline MatrixDouble& getDiffHdivN(const FieldApproximationBase base) { return getDiffN(base); };

    /** \brief get shape functions for Hdiv space
    */
    inline const MatrixDouble& getHdivN() const { return getN(bAse); };

    /** \brief get derivatives of shape functions for Hdiv space
    */
    inline const MatrixDouble& getDiffHdivN() const { return getDiffN(bAse); };

    /** \brief get shape functions for Hdiv space
    */
    inline MatrixDouble& getHdivN() { return getN(bAse); };

    /** \brief get derivatives of shape functions for Hdiv space
    */
    inline MatrixDouble& getDiffHdivN() { return getDiffN(bAse); };

    /** \brief get Hdiv of shape functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getHdivN(const FieldApproximationBase base,const int gg) {
      int dim = 3;
      int nb_dofs = getHdivN(base).size2()/dim;
      double *data = &getHdivN(base)(gg,0);
      return MatrixAdaptor(nb_dofs,dim,ublas::shallow_array_adaptor<double>(dim*nb_dofs,data));
    }

    /** \brief get Hdiv of shape functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getHdivN(const int gg) {
      return getHdivN(bAse,gg);
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(FieldApproximationBase base,const int gg) {
      int nb_dofs = getDiffHdivN(base).size2()/9;
      double *data = &getDiffHdivN(base)(gg,0);
      return MatrixAdaptor(nb_dofs,9,ublas::shallow_array_adaptor<double>(9*nb_dofs,data));
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(const int gg) {
      return getDiffHdivN(bAse,gg);
    }

    /** \brief get DiffHdiv of shape functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(
      const FieldApproximationBase base,const int dof,const int gg
    ) {
      double *data = &getDiffHdivN(base)(gg,9*dof);
      return MatrixAdaptor(3,3,ublas::shallow_array_adaptor<double>(9,data));
    }


    /** \brief get DiffHdiv of shape functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of shape functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(const int dof,const int gg) {
      return getDiffHdivN(bAse,dof,gg);
    }

    inline FTensor::Tensor0<double*> getFTensorN(const FieldApproximationBase base) {
      double *ptr = &*getN(base).data().begin();
      return FTensor::Tensor0<double*>(ptr);
    }


    friend std::ostream& operator<<(std::ostream& os,const DataForcesAndSurcesCore::EntData &e);

  protected:
    int sEnse;                    ///< Entity sense (orientation)
    ApproximationOrder oRder;     ///< Entity order
    FieldSpace sPace;             ///< Entity space
    FieldApproximationBase bAse;  ///< Field approximation base
    VectorInt iNdices;            ///< Global indices on entity
    VectorInt localIndices;       ///< Local indices on entity
    VectorDofs dOfs;              ///< DoFs on entity
    VectorDouble fieldData;       ///< Field data on entity
    ShapeFunctionBasesVector N;     ///< Base functions
    ShapeFunctionBasesVector diffN; ///< Derivatives of base functions
  };

  std::bitset<LASTSPACE> sPace;   ///< spaces on element
  std::bitset<LASTBASE> bAse;    ///< bases on element
  ublas::matrix<int> facesNodes; 			                  ///< nodes on finite element faces
  std::bitset<LASTSPACE> spacesOnEntities[MBMAXTYPE]; 	      ///< spaces on entity types
  std::bitset<LASTBASE> basesOnEntities[MBMAXTYPE]; 	        ///< bases on entity types
  boost::ptr_vector<EntData> dataOnEntities[MBMAXTYPE]; ///< data on nodes, shape function, dofs values, etc.

  DataForcesAndSurcesCore(EntityType type);
  virtual ~DataForcesAndSurcesCore() {}

  friend std::ostream& operator<<(std::ostream& os,const DataForcesAndSurcesCore &e);

  protected:
  DataForcesAndSurcesCore() {}

};

/** \brief this class derive data form other data structure
  * \ingroup mofem_forces_and_sources
  *
  *
  * It behaves like normal data structure it is used to share base functions with
  * other data structures. Dofs values, approx. order and
  * indices are not shared.
  *
  * Shape functions, senses are shared with other data structure.
  *
  */
struct DerivedDataForcesAndSurcesCore: public DataForcesAndSurcesCore  {

  struct DerivedEntData: public DataForcesAndSurcesCore::EntData {

    DataForcesAndSurcesCore::EntData &entData;
    DerivedEntData(DataForcesAndSurcesCore::EntData &ent_data):
    entData(ent_data) {
    }

    int getSense() const { return entData.getSense(); }

    boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) {
      return entData.getNSharedPtr(base);
    }
    boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) {
      return entData.getDiffNSharedPtr(base);
    }
    const boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) const {
      return entData.getNSharedPtr(base);
    }
    const boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) const {
      return entData.getDiffNSharedPtr(base);
    }


  };

  DerivedDataForcesAndSurcesCore(DataForcesAndSurcesCore &data);

};

}

#endif //__DATASTRUCTURES_HPP

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources Forces and sources
 ******************************************************************************/
