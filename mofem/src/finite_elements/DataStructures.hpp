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

/**
* \brief Get tensor rank 0 (scalar) form data vector
* \ingroup mofem_forces_and_sources_user_data_operators

Example how to use it.
\code
VectorDouble vec;
vec.resize(nb_gauss_pts,false);
vec.clear();
FTensor::Tensor0<double*> t0 = getTensor0FormData(data);
for(int gg = 0;gg!=nb_gauss_pts;gg++) {

  ++t0;
}
\endcode

*/
template<class T, class A>
FTensor::Tensor0<T*> getTensor0FormData(
  ublas::vector<T,A> &data
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor0<T*>();
}

template<>
FTensor::Tensor0<double*> getTensor0FormData<double,ublas::unbounded_array<double> >(
  ublas::vector<double,ublas::unbounded_array<double> > &data
);

/**
 * \brief Get tensor rank 1 (vector) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim, class T, class L, class A>
FTensor::Tensor1<T*,Tensor_Dim> getTensor1FormData(
  ublas::matrix<T,L,A> &data
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  s << " and dim = " << Tensor_Dim;
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor1<T*,Tensor_Dim>();
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim>
FTensor::Tensor1<double*,Tensor_Dim> getTensor1FormData(
  MatrixDouble &data
) {
  return getTensor1FormData<
  Tensor_Dim,double,ublas::row_major,ublas::unbounded_array<double>
  >(data);
}


template<>
FTensor::Tensor1<double*,3> getTensor1FormData<3,double,ublas::row_major,ublas::unbounded_array<double> >(
  MatrixDouble &data
);

template<>
FTensor::Tensor1<double*,2> getTensor1FormData<2,double,ublas::row_major,ublas::unbounded_array<double> >(
  MatrixDouble &data
);

/**
 * \brief Get tensor rank 2 (matrix) form data matrix
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim0, int Tensor_Dim1, class T, class L, class A>
FTensor::Tensor2<T*,Tensor_Dim0,Tensor_Dim1> getTensor2FormData(
  ublas::matrix<T,L,A> &data
) {
  std::stringstream s;
  s << "Not implemented for T = " << typeid(T).name();
  s << " and dim0 = " << Tensor_Dim0;
  s << " dim1 = " << Tensor_Dim1;
  THROW_MESSAGE(s.str());
  // return FTensor::Tensor1<T*,Tensor_Dim>();
}

template<>
FTensor::Tensor2<double*,3,3> getTensor2FormData(
  MatrixDouble &data
);

/**
 * \brief Get tensor rank 2 (matrix) form data matrix (specialization)
 * \ingroup mofem_forces_and_sources_user_data_operators
 */
template<int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getTensor2FormData(
  MatrixDouble &data
) {
  return getTensor2FormData<
  Tensor_Dim0,Tensor_Dim1,double,ublas::row_major,ublas::unbounded_array<double>
  >(data);
}

/** \brief data structure for finite element entity
  * \ingroup mofem_forces_and_sources
  *
  * It keeps that about indices of degrees of freedom, dofs data, base functions
  * functions, entity side number, type of entities, approximation order, etc.
  *
  */
struct DataForcesAndSurcesCore {

  /** \brief Data on single entity (This is passed as argument to DataOperator::doWork)
    * \ingroup mofem_forces_and_sources_user_data_operators
    */
  struct EntData {

    EntData();
    virtual ~EntData();

    /// \brief get entity sense, need to calculate base functions with conforming approximation fields
    virtual int getSense() const { return sEnse; }

    /// \brief get approximation order
    inline ApproximationOrder getOrder() const { return oRder; }

    /// \brief Get global indices of dofs on entity
    inline const VectorInt& getIndices() const { return iNdices; }

    /// \brief get global indices of dofs on entity up to given order
    inline const VectorIntAdaptor getIndicesUpToOrder(int order) {
      unsigned int size = 0;
      if(iNdices.size()) {
        size = dOfs[0]->getOrderNbDofs(order)*dOfs[0]->getNbOfCoeffs();
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
        size = dOfs[0]->getOrderNbDofs(order)*dOfs[0]->getNbOfCoeffs();
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
        size = dOfs[0]->getOrderNbDofs(order)*dOfs[0]->getNbOfCoeffs();
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

    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1FieldData()  {
      std::stringstream s;
      s << "Not implemented for this dimension dim = " << Tensor_Dim;
      THROW_MESSAGE(s.str());
    }

    FTensor::Tensor0<double*> getFTensor0FieldData();

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
     * Get shared pointer to base base functions
     */
    virtual boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) {
      return N[base];
    }

    /**
     * Get shared pointer to base base functions
     */
    virtual const boost::shared_ptr<MatrixDouble>& getNSharedPtr(const FieldApproximationBase base) const {
      return N[base];
    }

    /**
    * Get shared pointer to derivatives of base base functions
    */
    virtual boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) {
      return diffN[base];
    }

    /**
    * Get shared pointer to derivatives of base base functions
    */
    virtual const boost::shared_ptr<MatrixDouble>& getDiffNSharedPtr(const FieldApproximationBase base) const {
      return diffN[base];
    }

    // ************ H1/L2 ************

    /** \brief get base functions
    * this return matrix (nb. of rows is equal to nb. of Gauss pts, nb. of
    * columns is equal to number of base functions on this entity
    */
    virtual const MatrixDouble& getN(const FieldApproximationBase base) const {
      return *(getNSharedPtr(base));
    }

    inline const MatrixDouble& getN() const { return getN(bAse); }

    /** \brief get derivatives of base functions
     *
     * Matrix at rows has nb. of Gauss pts, at columns it has derivative of
     * base functions. Columns are structured as follows, [ dN1/dx, dN1/dy,
     * dN1/dz, dN2/dx, dN2/dy, dN2/dz, ... ]
     *
     * Note that base functions are calculated in file H1.c
     * Above description not apply for derivatives of nodal functions, since
     * derivative of nodal functions in case of simplexes, EDGES, TRIANGLES and
     * TETS are constant. So that matrix rows represents nb. of base
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
     * \brief Get base functions
     * @param  base Approximation base
     * @return      Error code
     */
    inline MatrixDouble& getN(const FieldApproximationBase base) { return *(getNSharedPtr(base)); }

    /**
     * \brief Get base functions
     *
     * It assumed that approximation base for given field is known and stored in this data structure
     *
     * @return Error code
     */
     inline MatrixDouble& getN() { return getN(bAse); }

    /**
     * \brief Get derivatives of base functions
     * @param  base Approximation base
     * @return      Error code
     */
     inline MatrixDouble& getDiffN(const FieldApproximationBase base) { return *(getDiffNSharedPtr(base)); }

    /**
     * \brief Get derivatives of base functions
     *
     * It assumed that approximation base for given field is known and stored in this data structure
     *
     * @return Error code
     */
     inline MatrixDouble& getDiffN() { return getDiffN(bAse); }

    /// \brief get base functions at Gauss pts
    inline const VectorAdaptor getN(const FieldApproximationBase base,const int gg) {
      int size = getN(base).size2();
      double *data = &getN(base)(gg,0);
      return VectorAdaptor(size,ublas::shallow_array_adaptor<double>(size,data));
    }

    /// \brief get base functions at Gauss pts
    inline const VectorAdaptor getN(const int gg) { return getN(bAse,gg); }

    /** \brief get derivative of base functions at Gauss pts

    * returned matrix on rows has base functions, in column its derivatives.
    *
    * \param base Approximation base
    * \param gg Nb. of Gauss pts.
    *
    */
    inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,const int gg) {
      //FIXME: That is bug, it will not work if number of integration pts is equal to number of nodes on entity.
      //User who not implementing low level DataOperator will not experience this.
      if(getN(base).size1() == getDiffN(base).size1()) {
        int size = getN(base).size2();
        int dim = getDiffN(base).size2()/size;
        double *data = &getDiffN(base)(gg,0);
        return MatrixAdaptor(
          getN(base).size2(),dim,ublas::shallow_array_adaptor<double>(getDiffN(base).size2(),data)
        );
      } else {
        // in some cases, f.e. for derivatives of nodal base functions at only one
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

    /** \brief get derivative of base functions at Gauss pts

    * returned matrix on rows has base functions, in column its derivatives.
    *
    * \param gg nb. of Gauss pts.
    *
    */
    inline const MatrixAdaptor getDiffN(const int gg) { return getDiffN(bAse,gg); }

    /** \brief get base functions at Gauss pts

    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * base functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param base Approximation base
    * \param gg number of Gauss point
    * \param nb_base_functions number of of base functions returned

    */
    inline const VectorAdaptor getN(const FieldApproximationBase base,const int gg,const int nb_base_functions) {
      (void)getN()(gg,nb_base_functions-1); // throw error if nb_base_functions is to big
      double *data = &getN(base)(gg,0);
      return VectorAdaptor(nb_base_functions,ublas::shallow_array_adaptor<double>(nb_base_functions,data));
    }

    /** \brief get base functions at Gauss pts

    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * base functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param gg number of Gauss point
    * \param nb_base_functions number of of base functions returned

    */
    inline const VectorAdaptor getN(const int gg,const int nb_base_functions) {
      return getN(bAse,gg,nb_base_functions);
    }

    /** \brief get derivatives of base functions at Gauss pts
    *
    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * base functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param nb_base_functions number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,const int gg,const int nb_base_functions) {
      //FIXME: That is bug, it will not work if number of integration pts is equal to number of nodes on entity.
      //User who not implementing low level DataOperator will not experience this.
      if(getN(base).size1() == getDiffN(base).size1()) {
        (void)getN(base)(gg,nb_base_functions-1); // throw error if nb_base_functions is to big
        int dim = getDiffN(base).size2()/getN(base).size2();
        double *data = &getDiffN(base)(gg,0);
        return MatrixAdaptor(
          nb_base_functions,dim,ublas::shallow_array_adaptor<double>(
            dim*nb_base_functions,data
          )
        );
      } else {
        // in some cases, f.e. for derivatives of nodal base functions only one
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

    /** \brief get derivatives of base functions at Gauss pts
    *
    * Note that multi field element, two different field can have different
    * approximation orders. Since we use hierarchical approximation basis,
    * base functions are calculated once for element, using maximal
    * approximation order on given entity.
    *
    * \param gg nb. of Gauss point
    * \param nb_base_functions number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffN(const int gg,const int nb_base_functions) {
      return getDiffN(bAse,gg,nb_base_functions);
    }

    // ************ HDIV ************

    inline const MatrixDouble& getHdivN(const FieldApproximationBase base) const { return getN(base); };
    inline const MatrixDouble& getDiffHdivN(const FieldApproximationBase base) const { return getDiffN(base); };
    inline MatrixDouble& getHdivN(const FieldApproximationBase base) { return getN(base); };
    inline MatrixDouble& getDiffHdivN(const FieldApproximationBase base) { return getDiffN(base); };

    /** \brief get base functions for Hdiv space
    */
    inline const MatrixDouble& getHdivN() const { return getN(bAse); };

    /** \brief get derivatives of base functions for Hdiv space
    *
    * Note: In rows ale integration pts, columns are formatted that that
    * components of vectors and then derivatives, for example row for given
    * integration points is formatted in array
    * \f[
    * t_{0,0}, t_{1,0}, t_{1,0}, t_{0,1}, t_{1,1}, t_{1,1}, t_{0,2}, t_{1,2}, t_{1,2}
    * \f]
    * where comma express derivative, i.e. \f$t_{2,1} = \frac{\partial t_2}{\partial \xi_1}\f$
    *
    */
    inline const MatrixDouble& getDiffHdivN() const { return getDiffN(bAse); };

    /** \brief get base functions for Hdiv space
    */
    inline MatrixDouble& getHdivN() { return getN(bAse); };

    /** \brief Get derivatives of base functions for Hdiv space
    *
    */
    inline MatrixDouble& getDiffHdivN() { return getDiffN(bAse); };

    /** \brief get Hdiv of base functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    *
    */
    inline const MatrixAdaptor getHdivN(const FieldApproximationBase base,const int gg) {
      const int dim = 3;
      int nb_base_functions = getHdivN(base).size2()/dim;
      double *data = &getHdivN(base)(gg,0);
      return MatrixAdaptor(nb_base_functions,dim,ublas::shallow_array_adaptor<double>(dim*nb_base_functions,data));
    }

    /** \brief get Hdiv of base functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getHdivN(const int gg) {
      return getHdivN(bAse,gg);
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(FieldApproximationBase base,const int gg) {
      int nb_base_functions = getDiffHdivN(base).size2()/9;
      double *data = &getDiffHdivN(base)(gg,0);
      return MatrixAdaptor(
        nb_base_functions,9,ublas::shallow_array_adaptor<double>(9*nb_base_functions,data)
      );
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(const int gg) {
      return getDiffHdivN(bAse,gg);
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(
      const FieldApproximationBase base,const int dof,const int gg
    ) {
      double *data = &getDiffHdivN(base)(gg,9*dof);
      return MatrixAdaptor(3,3,ublas::shallow_array_adaptor<double>(9,data));
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHdivN(const int dof,const int gg) {
      return getDiffHdivN(bAse,dof,gg);
    }

    // ************ HCURL ************

    inline const MatrixDouble& getHcurlN(const FieldApproximationBase base) const { return getN(base); };
    inline const MatrixDouble& getDiffHcurlN(const FieldApproximationBase base) const { return getDiffN(base); };
    inline MatrixDouble& getHcurlN(const FieldApproximationBase base) { return getN(base); };
    inline MatrixDouble& getDiffHcurlN(const FieldApproximationBase base) { return getDiffN(base); };

    /** \brief get base functions for Hdiv space
    */
    inline const MatrixDouble& getHcurlN() const { return getN(bAse); };

    /** \brief get derivatives of base functions for Hdiv space
    *
    * Note: In rows ale integration pts, columns are formatted that that
    * components of vectors and then derivatives, for example row for given
    * integration points is formatted in array
    * \f[
    * t_{0,0}, t_{1,0}, t_{1,0}, t_{0,1}, t_{1,1}, t_{1,1}, t_{0,2}, t_{1,2}, t_{1,2}
    * \f]
    * where comma express derivative, i.e. \f$t_{2,1} = \frac{\partial t_2}{\partial \xi_1}\f$
    *
    */
    inline const MatrixDouble& getDiffHcurlN() const { return getDiffN(bAse); };

    /** \brief get base functions for Hdiv space
    */
    inline MatrixDouble& getHcurlN() { return getN(bAse); };

    /** \brief Get derivatives of base functions for Hdiv space
    *
    */
    inline MatrixDouble& getDiffHcurlN() { return getDiffN(bAse); };


    /** \brief get Hcurl of base functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getHcurlN(const FieldApproximationBase base,const int gg) {
      const int dim = 3;
      int nb_base_functions = getHcurlN(base).size2()/dim;
      double *data = &getHcurlN(base)(gg,0);
      return MatrixAdaptor(nb_base_functions,dim,ublas::shallow_array_adaptor<double>(dim*nb_base_functions,data));
    }

    /** \brief get Hcurl of base functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getHcurlN(const int gg) {
      return getHcurlN(bAse,gg);
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param base Approximation base
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHcurlN(FieldApproximationBase base,const int gg) {
      int nb_base_functions = getDiffHcurlN(base).size2()/9;
      double *data = &getDiffHcurlN(base)(gg,0);
      return MatrixAdaptor(
        nb_base_functions,9,ublas::shallow_array_adaptor<double>(9*nb_base_functions,data)
      );
    }

    /** \brief get DiffHdiv of base functions at Gauss pts
    *
    * \param gg nb. of Gauss point
    * \param number of of base functions
    *
    */
    inline const MatrixAdaptor getDiffHcurlN(const int gg) {
      return getDiffHcurlN(bAse,gg);
    }

    // ********* Tensors *******

    /**
     * \brief Get base function as Tensor0
     *
     * \param base
     * \return Tensor0
     *
     * operator++() is moving pointer to next base function.
     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N(const FieldApproximationBase base) {
      double *ptr = &*getN(base).data().begin();
      return FTensor::Tensor0<double*>(ptr);
    };

    /**
     * \brief Get base function as Tensor0
     *
     * Return base functions for field base
     *
     * \return Tensor0
     *
     * operator++() is moving pointer to next base function.
     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N() {
      return getFTensor0N(bAse);
    };

    /**
     * \brief Get base function as Tensor0 (Loop by integration points)
     *
     * \param base
     * \param bb base function
     * \return Tensor0

     Note that:
     \code
     t0 = data.getFTensor0N(base,bb);
     ++t0
     \endcode
     Increment in above code will move pointer to base function in next integration
     point.

     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N(const FieldApproximationBase base,const int bb) {
      double *ptr = &getN(base)(0,bb);
      return FTensor::Tensor0<double*>(ptr,getN(base).size2());
    };

    /**
     * \brief Get base function as Tensor0 (Loop by integration points)
     *
     * Return base functions for field base
     *
     * \param bb base function
     * \return Tensor0
     *
     * operator++() is moving pointer to next integration point.
     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N(const int bb) {
      return getFTensor0N(bAse,bb);
    };

    /**
     * \brief Get base function as Tensor0 (Loop by integration points)
     *
     * \param base
     * \param gg integration points
     * \param bb base function
     * \return Tensor0

     Note that:
     \code
     t0 = data.getFTensor0N(base,bb);
     ++t0
     \endcode
     Increment in above code will move pointer to base function in next integration
     point.

     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N(const FieldApproximationBase base,const int gg,const int bb) {
      double *ptr = &getN(base)(gg,bb);
      return FTensor::Tensor0<double*>(ptr);
    };

    /**
     * \brief Get base function as Tensor0 (Loop by integration points)
     *
     * Return base functions for field base
     *
     * \param bb base function
     * \return Tensor0
     *
     * operator++() is moving pointer to next base function.
     *
     */
    inline FTensor::Tensor0<double*> getFTensor0N(const int gg,const int bb) {
      return getFTensor0N(bAse,gg,bb);
    };

    /**
     * \brief Get derivatives of base functions
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>();
     * \endcode
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>();
     * \endcode
     *
     * \param base functions
     * \return Tensor rank 1 (vector)
     *
     * operator++() is moving pointer to next base function.
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN(const FieldApproximationBase base);

    /**
     * \brief Get derivatives of base functions
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>();
     * \endcode
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>();
     * \endcode
     *
     * \return Tensor rank 1 (vector)
     *
     * operator++() is moving pointer to next base function.
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN();

    /**
     * \brief Get derivatives of base functions (Loop by integration points)
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>(base,bb);
     * \endcode
     * where bb is base function. Operator ++diff_base will move tensor pointer
     * to next integration point.
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>(base,bb);
     * \endcode
     *
     * \param base functions
     * \return Tensor rank 1 (vector)
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN(
      const FieldApproximationBase base,const int bb
    );

    /**
     * \brief Get derivatives of base functions (Loop by integration points)
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>(bb);
     * \endcode
     * where bb is base function. Operator ++diff_base will move tensor pointer
     * to next integration point.
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>(bb);
     * \endcode
     *
     * \return Tensor rank 1 (vector)
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN(const int bb);

    /**
     * \brief Get derivatives of base functions (Loop by integration points)
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>(base,gg,bb);
     * \endcode
     * where bb is base function and gg is integration pt. Operator ++diff_base
     * will move tensor pointer to next integration point.
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>(base,gg,bb);
     * \endcode
     *
     * \return Tensor rank 1 (vector)
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN(
      const FieldApproximationBase base,const int gg,const int bb
    );

    /**
     * \brief Get derivatives of base functions (Loop by integration points)
     *
     * For volume element like tetrahedral or prism,
     * \code
     * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>(gg,bb);
     * \endcode
     * where bb is base function and gg is integration pt. Operator ++diff_base
     * will move tensor pointer to next integration point.
     *
     * For face element like triangle or quad
     * \code
     * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>(gg,bb);
     * \endcode
     *
     * \return Tensor rank 1 (vector)
     *
     */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1DiffN(const int gg,const int bb);

    /** \brief Get base functions for Hdiv space

    \note You probably like to use getFTensor1HdivN(), in typical use base is set
    automatically based on base set to field.

    * @param  base Approximation base

    Example:
    \code
    FTensor::Index<'i',3> i;
    int nb_dofs = data.getFieldData().size();
    FTensor::Tensor1<double*,3> t_n_hdiv = data.getFTensor1HdivN<3>();
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      int ll = 0;
      for(;ll!=nb_dofs;ll++) {
        double dot_product = t_n_hdiv(i)*t_n_hdiv(i);
        ++t_n_hdiv;
      }
      for(;ll!=data.getHdivN().size2()/3;ll++) {
        ++t_n_hdiv;
      }
    }
    \endcode

    */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1HdivN(FieldApproximationBase base);

    /** \brief Get base functions for Hdiv space

    Example:
    \code
    FTensor::Index<'i',3> i;
    int nb_dofs = data.getFieldData().size();
    FTensor::Tensor1<double*,3> t_n_hdiv = data.getFTensor1HdivN<3>();
    for(int gg = 0;gg!=nb_gauss_pts;gg++) {
      int ll = 0;
      for(;ll!=nb_dofs;ll++) {
        double dot_product = t_n_hdiv(i)*t_n_hdiv(i);
        ++t_n_hdiv;
      }
      for(;ll!=data.getHdivN().size2()/3;ll++) {
        ++t_n_hdiv;
      }
    }
    \endcode

    */
    template<int Tensor_Dim>
    FTensor::Tensor1<double*,Tensor_Dim> getFTensor1HdivN();

    /** \brief Get derivatives of base functions for Hdiv space
    */
    template<int Tensor_Dim0,int Tensor_Dim1>
    FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getFTensor2DiffHdivN(FieldApproximationBase base);

    /** \brief Get derivatives of base functions for Hdiv space
    */
    template<int Tensor_Dim0,int Tensor_Dim1>
    FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getFTensor2DiffHdivN();

    /** \brief Get base functions for Hcurl space
    */
    template<int Tensor_Dim>
    inline FTensor::Tensor1<double*,Tensor_Dim> getFTensor1HcurlN(FieldApproximationBase base) {
      return getFTensor1HdivN<Tensor_Dim>(base);
    }

    /** \brief Get base functions for Hcurl space
    */
    template<int Tensor_Dim>
    inline FTensor::Tensor1<double*,Tensor_Dim> getFTensor1HcurlN() {
      return getFTensor1HcurlN<Tensor_Dim>(bAse);
    }

    /** \brief Get derivatives of base functions for Hcurl space
    */
    template<int Tensor_Dim0,int Tensor_Dim1>
    inline FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getFTensor2DiffHcurlN(FieldApproximationBase base) {
      return getFTensor2DiffHdivN<Tensor_Dim0,Tensor_Dim1>(base);
    }

    /** \brief Get derivatives of base functions for Hcurl space
    */
    template<int Tensor_Dim0,int Tensor_Dim1>
    FTensor::Tensor2<double*,Tensor_Dim0,Tensor_Dim1> getFTensor2DiffHcurlN() {
      return getFTensor2DiffHdivN<Tensor_Dim0,Tensor_Dim1>(bAse);
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
  boost::ptr_vector<EntData> dataOnEntities[MBMAXTYPE]; ///< data on nodes, base function, dofs values, etc.

  DataForcesAndSurcesCore(EntityType type);
  virtual ~DataForcesAndSurcesCore() {}

  friend std::ostream& operator<<(std::ostream& os,const DataForcesAndSurcesCore &e);

  protected:
  DataForcesAndSurcesCore() {}

};

template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1FieldData<3>();
template<>
FTensor::Tensor1<double*,2> DataForcesAndSurcesCore::EntData::getFTensor1FieldData<2>();

template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<3>(
  const FieldApproximationBase base
);
template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<3>();

template<>
FTensor::Tensor1<double*,2> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<2>(
  const FieldApproximationBase base
);
template<>
FTensor::Tensor1<double*,2> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<2>();
template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<3>(
  const FieldApproximationBase base,const int bb
);
template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<3>(
  const int bb
);
template<>
FTensor::Tensor1<double*,2> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<2>(
  const FieldApproximationBase base,const int bb
);
template<>
FTensor::Tensor1<double*,2> DataForcesAndSurcesCore::EntData::getFTensor1DiffN<2>(
  const int bb
);

// Specifications for HDiv/HCrul

template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1HdivN<3>(
  FieldApproximationBase base
);
template<>
FTensor::Tensor1<double*,3> DataForcesAndSurcesCore::EntData::getFTensor1HdivN<3>();
template<>
FTensor::Tensor2<double*,3,3> DataForcesAndSurcesCore::EntData::getFTensor2DiffHdivN<3,3>(
  FieldApproximationBase base
);
template<>
FTensor::Tensor2<double*,3,3> DataForcesAndSurcesCore::EntData::getFTensor2DiffHdivN<3,3>();


// Specifications 2d 

template<>
FTensor::Tensor2<double*,3,2> DataForcesAndSurcesCore::EntData::getFTensor2DiffHdivN<3,2>(
  FieldApproximationBase base
);
template<>
FTensor::Tensor2<double*,3,2> DataForcesAndSurcesCore::EntData::getFTensor2DiffHdivN<3,2>();


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
 *
 * \brief Manages complexities related to assembly of vector and matrices at single finite element level.
 *
 ******************************************************************************/
