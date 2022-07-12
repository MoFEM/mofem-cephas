/** \file EntitiesFieldData.hpp

\brief Data structures for accessing information about finite element and its
degrees of freedom.

*/

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __ENTITIES_FIELD_DATA_HPP__
#define __ENTITIES_FIELD_DATA_HPP__

using namespace boost::numeric;

namespace MoFEM {

using DofsAllocator = ublas::unbounded_array<

    FEDofEntity *, std::allocator<FEDofEntity *>

    >;

using VectorDofs = ublas::vector<FEDofEntity *, DofsAllocator>;

using FieldEntAllocator = ublas::unbounded_array<

    FieldEntity *, std::allocator<FieldEntity *>

    >;

using VectorFieldEntities = ublas::vector<FieldEntity *, FieldEntAllocator>;

/** \brief data structure for finite element entity
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 * It keeps that about indices of degrees of freedom, dofs data, base functions
 * functions, entity side number, type of entities, approximation order, etc.
 *
 */
struct EntitiesFieldData {

  struct EntData;

  std::bitset<LASTSPACE> sPace; ///< spaces on element
  std::bitset<LASTBASE> bAse;   ///< bases on element
  MatrixInt facesNodes;         ///< nodes on finite element faces
  MatrixInt facesNodesOrder;    ///< order of face nodes on element

  std::array<std::bitset<LASTSPACE>, MBMAXTYPE>
      spacesOnEntities; ///< spaces on entity types
  std::array<std::bitset<LASTBASE>, MBMAXTYPE>
      basesOnEntities; ///< bases on entity types
  std::array<std::bitset<LASTBASE>, LASTSPACE>
      basesOnSpaces; ///< base on spaces
  std::array<boost::ptr_vector<EntData>, MBMAXTYPE>
      dataOnEntities; ///< data on nodes, base
                      ///< function, dofs
                      ///< values, etc.

  EntitiesFieldData(const EntityType type);
  virtual ~EntitiesFieldData() = default;

  virtual MoFEMErrorCode setElementType(const EntityType type);

  /**
   * Reset data associated with particular field name
   * @return error code
   */
  inline MoFEMErrorCode resetFieldDependentData();

  /**
   * @brief Swap approximation base
   *
   * Bernstein-Bezier (BB) base is not hierarchical, and is calculated for
   * particular field, since it all shape functions change with the order. BB
   * base is precalculated for every field, and when user push operator with
   * paricular field using BB base, pointers to shape funtions and
   * BaseDerivatives of shape functions are set to particular location, once
   * operator is executed, pointers are switch back to its oroginal position.
   *
   * getNSharedPtr(base) <=== getBBNSharedPtr(field_name);
   * // DO OPERATOR WORK
   * getNSharedPtr(base) ==> getBBNSharedPtr(field_name);
   *
   * @param field_name
   * @param base
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode baseSwap(const std::string &field_name,
                                  const FieldApproximationBase base);

  friend std::ostream &operator<<(std::ostream &os, const EntitiesFieldData &e);

protected:
  EntitiesFieldData() {}
};

/** \brief this class derive data form other data structure
 * \ingroup mofem_forces_and_sources_user_data_operators
 *
 *
 * It behaves like normal data structure it is used to share base functions with
 * other data structures. Dofs values, approx. order and
 * indices are not shared.
 *
 * Shape functions, senses are shared with other data structure.
 *
 */
struct DerivedEntitiesFieldData : public EntitiesFieldData {

  struct DerivedEntData;

  DerivedEntitiesFieldData(
      const boost::shared_ptr<EntitiesFieldData> &data_ptr);
  MoFEMErrorCode setElementType(const EntityType type);

private:
  const boost::shared_ptr<EntitiesFieldData> dataPtr;
};

/** \brief Data on single entity (This is passed as argument to
 * DataOperator::doWork) \ingroup mofem_forces_and_sources_user_data_operators
 * \nosubgrouping
 *
 * \todo Hdiv and Hcurl functions should be accessed through common interface.
 */
struct EntitiesFieldData::EntData {

  enum BaseDerivatives {
    ZeroDerivative = 0,
    FirstDerivative,
    SecondDerivative,
    ThirdDerivative,
    ForthDerivative,
    LastDerivative
  };

  /** \name Constructor and destructor */

  /**@{*/

  EntData(const bool allocate_base_matrices = true);
  virtual ~EntData() = default;

  /**@}*/

  /** \name Sense, order and indices */

  /**@{*/

  /// \brief get entity sense, need to calculate base functions with
  /// conforming approximation fields
  virtual int getSense() const;

  /// \brief get approximation order
  inline ApproximationOrder getOrder() const;

  /// \brief Get global indices of dofs on entity
  inline const VectorInt &getIndices() const;

  /// \brief get global indices of dofs on entity up to given order
  inline const VectorIntAdaptor getIndicesUpToOrder(int order);

  /// \brief get local indices of dofs on entity
  inline const VectorInt &getLocalIndices() const;

  /// \brief get local indices of dofs on entity up to given order
  inline const VectorIntAdaptor getLocalIndicesUpToOrder(int order);

  inline int &getSense();

  inline ApproximationOrder &getOrder();

  inline VectorInt &getIndices();

  inline VectorInt &getLocalIndices();

  /**@}*/

  /** \name Data on entity */

  /**@{*/

  /// \brief get dofs values
  inline const VectorDouble &getFieldData() const;

  /// \brief get dofs values up to given order
  inline const VectorAdaptor getFieldDataUpToOrder(int order);

  /// \brief get dofs data stature FEDofEntity
  inline const VectorDofs &getFieldDofs() const;

  /// \brief get dofs data stature FEDofEntity
  inline VectorDouble &getFieldData();

  /// \brief get field entities
  inline const VectorFieldEntities &getFieldEntities() const;

  /// \brief get field entities
  inline VectorFieldEntities &getFieldEntities();

  //// \brief get entity bit ref level
  virtual BitRefLevel &getEntDataBitRefLevel();

  /**
   * @brief Return FTensor of rank 1, i.e. vector from filed data coeffinects
   *
   * \code
   * auto t_vec = data.getFTensor1FieldData<3>();
   * \endcode
   *
   * @tparam Tensor_Dim size of vector
   * @return FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>,
   * Tensor_Dim>
   */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1FieldData();

  /**
   * @brief  Return FTensor rank 2, i.e. matrix from filed data coeffinects
   *
   * \code
   * auto t_mat = data.getFTensor2FieldData<3,3>();
   * \endcode
   *
   * @tparam Tensor_Dim0
   * @tparam Tensor_Dim1
   * @return FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 *
   * Tensor_Dim1>, Tensor_Dim0, Tensor_Dim1>
   */
  template <int Tensor_Dim0, int Tensor_Dim1>
  FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                   Tensor_Dim0, Tensor_Dim1>
  getFTensor2FieldData();

  /**
   * @brief  Return symmetric FTensor rank 2, i.e. matrix from filed data
   * coeffinects
   *
   * \code
   * auto t_mat = data.getFTensor2SymmetricFieldData<3>();
   * \endcode
   *
   * @tparam Tensor_Dim dimension of the tensor
   * @return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, (Tensor_Dim
   * * (Tensor_Dim + 1)) / 2>, Tensor_Dim>
   */
  template <int Tensor_Dim>
  FTensor::Tensor2_symmetric<
      FTensor::PackPtr<double *, (Tensor_Dim * (Tensor_Dim + 1)) / 2>,
      Tensor_Dim>
  getFTensor2SymmetricFieldData();

  /**
   * @brief Resturn scalar files as a FTensor of rank 0
   *
   * @return FTensor::Tensor0<FTensor::PackPtr<double *,1> >
   */
  FTensor::Tensor0<FTensor::PackPtr<double *, 1>> getFTensor0FieldData();

  inline VectorDofs &getFieldDofs();

  /**@}*/

  /** \name Base and space */

  /**@{*/

  /**
   * \brief Get approximation base
   * @return Approximation base
   */
  inline FieldApproximationBase &getBase();

  /**
   * \brief Get field space
   * @return Field space
   */
  inline FieldSpace &getSpace();

  /**
   * Get shared pointer to base base functions
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getNSharedPtr(const FieldApproximationBase base,
                const BaseDerivatives direvatie);

  /**
   * Get shared pointer to base base functions
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getNSharedPtr(const FieldApproximationBase base);

  /**
   * Get shared pointer to derivatives of base base functions
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getDiffNSharedPtr(const FieldApproximationBase base);

  /**@}*/

  /** \name Get base functions for H1/L2 */

  /**@{*/

  /** \brief get base functions
   * this return matrix (nb. of rows is equal to nb. of Gauss pts, nb. of
   * columns is equal to number of base functions on this entity.
   *
   * \note Note that for vectorial base, like Hdiv or Hcurl, in columns are
   * vectorial base functions. For tonsorial would be tonsorial base
   * functions. Interpretation depends on type of base, scalar, vectorial or
   * tonsorial and dimension fo problem.
   *
   */
  inline MatrixDouble &getN(const FieldApproximationBase base);

  /**
   * @copydoc MoFEM::EntitiesFieldData::EntData::getN
   */
  inline MatrixDouble &getN(const std::string &field_name);

  /**
   * @copydoc MoFEM::EntitiesFieldData::EntData::getN
   */
  inline MatrixDouble &getN();

  /** \brief get derivatives of base functions
   *
   * Matrix at rows has nb. of Gauss pts, at columns it has derivative of
   * base functions. Columns are structured as follows, [ dN1/dx, dN1/dy,
   * dN1/dz, dN2/dx, dN2/dy, dN2/dz, ... ]
   *
   * Scalar base functions:
   * Note that base functions are calculated in file H1.c
   * Above description not apply for derivatives of nodal functions, since
   * derivative of nodal functions in case of simplexes, EDGES, TRIANGLES and
   * TETS are constant. So that matrix rows represents nb. of base
   * functions, columns are derivatives. Nb. of columns depend on element
   * dimension, for EDGES is one, for TRIS is 2 and TETS is 3.
   *
   * \note Note that for node element this function make no sense.
   *
   * Tonsorial base functions:
   * \note Note: In rows ale integration pts, columns are formatted that that
   * components of vectors and then derivatives, for example row for given
   * integration points is formatted in array
   * \f[
   * t_{0,0}, t_{1,0}, t_{1,0}, t_{0,1}, t_{1,1}, t_{1,1}, t_{0,2}, t_{1,2},
   * t_{1,2} \f] where comma express derivative, i.e. \f$t_{2,1} =
   * \frac{\partial t_2}{\partial \xi_1}\f$
   *
   */
  inline MatrixDouble &getDiffN(const FieldApproximationBase base);

  /**
   * @brief Get base function derivative
   * 
   * @param base  base  
   * @param derivative derivative
   * @return MatrixDouble& 
   */
  inline MatrixDouble &getN(const FieldApproximationBase base,
                                const BaseDerivatives derivative);

  /**
   * @copydoc MoFEM::EntitiesFieldData::EntData::getDiffN
   */
  inline MatrixDouble &getDiffN(const std::string &field_name);

  /**
   * @copydoc MoFEM::EntitiesFieldData::EntData::getDiffN
   */
  inline MatrixDouble &getDiffN();

  /**
   * @brief Get base function derivative
   * 
   * @param derivative 
   * @return MatrixDouble& 
   */
  inline MatrixDouble &getN(const BaseDerivatives derivative);

  /// \brief get base functions at Gauss pts
  inline const VectorAdaptor getN(const FieldApproximationBase base,
                                  const int gg);

  /// \brief get base functions at Gauss pts
  inline const VectorAdaptor getN(const int gg);

  /** \brief get derivative of base functions at Gauss pts

  * returned matrix on rows has base functions, in column its derivatives.
  *
  * \param base Approximation base
  * \param gg Nb. of Gauss pts.
  *
  */
  inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,
                                      const int gg);

  /** \brief get derivative of base functions at Gauss pts

  * returned matrix on rows has base functions, in column its derivatives.
  *
  * \param gg nb. of Gauss pts.
  *
  */
  inline const MatrixAdaptor getDiffN(const int gg);

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
  inline const VectorAdaptor getN(const FieldApproximationBase base,
                                  const int gg, const int nb_base_functions);

  /** \brief get base functions at Gauss pts

  * Note that multi field element, two different field can have different
  * approximation orders. Since we use hierarchical approximation basis,
  * base functions are calculated once for element, using maximal
  * approximation order on given entity.
  *
  * \param gg number of Gauss point
  * \param nb_base_functions number of of base functions returned

  */
  inline const VectorAdaptor getN(const int gg, const int nb_base_functions);

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
  inline const MatrixAdaptor getDiffN(const FieldApproximationBase base,
                                      const int gg,
                                      const int nb_base_functions);

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
  inline const MatrixAdaptor getDiffN(const int gg,
                                      const int nb_base_functions);

  /**@}*/

  /** \name Get base functions for vectorial approximation basese, i.e.
   * Hdiv/Hcurl */

  /**@{*/

  /** \brief get Hdiv of base functions at Gauss pts
   *
   * \param base Approximation base
   * \param gg nb. of Gauss point
   *
   */
  template <int DIM>
  inline const MatrixAdaptor getVectorN(const FieldApproximationBase base,
                                        const int gg);

  /** \brief get Hdiv of base functions at Gauss pts
   *
   * \param gg nb. of Gauss point
   * \param number of of base functions
   *
   */
  template <int DIM> inline const MatrixAdaptor getVectorN(const int gg);

  /** \brief get DiffHdiv of base functions at Gauss pts
   *
   * \param base Approximation base
   * \param gg nb. of Gauss point
   * \param number of of base functions
   *
   */
  template <int DIM0, int DIM1>
  inline const MatrixAdaptor getVectorDiffN(FieldApproximationBase base,
                                            const int gg);

  /** \brief get DiffHdiv of base functions at Gauss pts
   *
   * \param gg nb. of Gauss point
   * \param number of of base functions
   *
   */
  template <int DIM0, int DIM1>
  inline const MatrixAdaptor getVectorDiffN(const int gg);

  /** \brief get DiffHdiv of base functions at Gauss pts
   *
   * \param base Approximation base
   * \param gg nb. of Gauss point
   * \param number of of base functions
   *
   */
  template <int DIM0, int DIM1>
  inline const MatrixAdaptor getVectorDiffN(const FieldApproximationBase base,
                                            const int dof, const int gg);

  /** \brief get DiffHdiv of base functions at Gauss pts
   *
   * \param gg nb. of Gauss point
   * \param number of of base functions
   *
   */
  template <int DIM0, int DIM1>
  inline const MatrixAdaptor getVectorDiffN(const int dof, const int gg);

  /**@}*/

  /** \name Get base functions with FTensor */

  /**@{*/

  /**
   * \brief Get base function as Tensor0
   *
   * \param base
   * \return Tensor0
   *
   */
  inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
  getFTensor0N(const FieldApproximationBase base);

  /**
   * \brief Get base function as Tensor0
   *
   * Return base functions for field base
   *
   * \return Tensor0
   *
   */
  inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>> getFTensor0N();

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
   Increment in above code will move pointer to base function in next
   integration point.

   *
   */
  inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
  getFTensor0N(const FieldApproximationBase base, const int gg, const int bb);

  /**
   * \brief Get base function as Tensor0 (Loop by integration points)
   *
   * Return base functions for field base
   *
   * \param bb base function
   * \return Tensor0
   *
   */
  inline FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
  getFTensor0N(const int gg, const int bb);

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
   */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1DiffN(const FieldApproximationBase base);

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
   */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1DiffN();

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
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1DiffN(const FieldApproximationBase base, const int gg,
                   const int bb);

  /**
   * \brief Get derivatives of base functions (Loop by integration points)
   *
   * For volume element like tetrahedral or prism,
   * \code
   * Tensor1<double*,3> diff_base = data.getFTensor1DiffN<3>(gg,bb);
   * \endcode
   * where bb is base function and gg is integration pt. Operator ++diff_base
   * will move tensor pointer to next base function.
   *
   * For face element like triangle or quad
   * \code
   * Tensor1<double*,2> diff_base = data.getFTensor1DiffN<2>(gg,bb);
   * \endcode
   *
   * \return Tensor rank 1 (vector)
   *
   */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1DiffN(const int gg, const int bb);

  /** \brief Get base functions for Hdiv/Hcurl spaces

  \note You probably like to use getFTensor1N(), in typical use base is
  set automatically based on base set to field.

  * @param  base Approximation base

  Example:
  \code
  FTensor::Index<'i',3> i;
  int nb_dofs = data.getFieldData().size();
  auto t_n_hdiv = data.getFTensor1N<3>();
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    int ll = 0;
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i)*t_n_hdiv(i);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1N(FieldApproximationBase base);

  /** \brief Get base functions for Hdiv space

  Example:
  \code
  FTensor::Index<'i',3> i;
  int nb_dofs = data.getFieldData().size();
  auto t_n_hdiv = data.getFTensor1N<3>();
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    int ll = 0;
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i)*t_n_hdiv(i);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim> auto getFTensor1N();

  /** \brief Get derivatives of base functions for Hdiv space
   */
  template <int Tensor_Dim0, int Tensor_Dim1>
  FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                   Tensor_Dim0, Tensor_Dim1>
  getFTensor2DiffN(FieldApproximationBase base);

  /** \brief Get derivatives of base functions for Hdiv space at integration
   * pts
   */
  template <int Tensor_Dim0, int Tensor_Dim1>
  FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                   Tensor_Dim0, Tensor_Dim1>
  getFTensor2DiffN(FieldApproximationBase base, const int gg, const int bb);

  /** \brief Get derivatives of base functions for Hdiv space
   */
  template <int Tensor_Dim0, int Tensor_Dim1>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                          Tensor_Dim0, Tensor_Dim1>
  getFTensor2DiffN() {
    return getFTensor2DiffN<Tensor_Dim0, Tensor_Dim1>(bAse);
  }

  /** \brief Get derivatives of base functions for Hdiv space at integration
   * pts
   */
  template <int Tensor_Dim0, int Tensor_Dim1>
  inline FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                          Tensor_Dim0, Tensor_Dim1>
  getFTensor2DiffN(const int gg, const int bb) {
    return getFTensor2DiffN<Tensor_Dim0, Tensor_Dim1>(bAse, gg, bb);
  }

  /** \brief Get second derivatives of base functions for Hvec space
   */
  template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  FTensor::Tensor3<
      FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1 * Tensor_Dim2>,
      Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>
  getFTensor3Diff2N(FieldApproximationBase base);

  /** \brief Get second derivatives of base functions for Hvec space
   */
  template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2>
  inline FTensor::Tensor3<
      FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1 * Tensor_Dim2>,
      Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>
  getFTensor3Diff2N() {
    return getFTensor3Diff2N<Tensor_Dim0, Tensor_Dim1, Tensor_Dim2>(bAse);
  }

  /**
   * \brief Get Hdiv base functions at integration point

   \code
   FTensor::Index<'i',3> i;
   for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    auto t_base = data.getFTensor1N(base,gg,bb);
    for(int bb = 0;bb!=nb_base_functions;bb++) {
      auto dot = t_base(i)*t_base(i);
    }
   }
   \endcode

   */
  template <int Tensor_Dim>
  FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
  getFTensor1N(FieldApproximationBase base, const int gg, const int bb);

  /**
   * \brief Get Hdiv base functions at integration point

   \code
   FTensor::Index<'i',3> i;
   for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    auto t_base = data.getFTensor1N(gg,0);
    for(int bb = 0;bb!=nb_base_functions;bb++) {
      double dot = t_base(i)*t_base(i);
    }
   }
   \endcode

   */
  template <int Tensor_Dim>
  inline auto getFTensor1N(const int gg, const int bb);

  /** \brief Get base functions for Hdiv/Hcurl spaces

  \note You probably like to use getFTensor1N(), in typical use base is
  set automatically based on base set to field.

  * @param  base Approximation base

  Example:
  \code
  FTensor::Index<'i',3> i;
  FTensor::Index<'i',3> j;
  int nb_dofs = data.getFieldData().size();
  auto t_n_hdiv = data.getFTensor2N<3,3>();
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    int ll = 0;
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i,j)*t_n_hdiv(i,j);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim0, int Tensor_Dim1>
  FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                   Tensor_Dim0, Tensor_Dim1>
  getFTensor2N(FieldApproximationBase base);

  /** \brief Get base functions for Hdiv space

  Example:
  \code
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;

  int nb_dofs = data.getFieldData().size();
  auto t_n_hdiv = data.getFTensor2N<3,3>();
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    int ll = 0;
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i,j)*t_n_hdiv(i,j);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim0, int Tensor_Dim1> auto getFTensor2N();

  /** \brief Get base functions for tensor Hdiv/Hcurl spaces

  \note You probably like to use getFTensor2N(), in typical use base is
  set automatically based on base set to field.

  @param  base Approximation base

  Example:
  \code
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> i;
  int nb_dofs = data.getFieldData().size();
  for(int gg = 0;gg!=nb_gauss_pts;gg++) {
    auto t_n_hdiv = data.getFTensor2N<3>(base,gg,bb);
    int ll = 0;
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i,j)*t_n_hdiv(i,j);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim0, int Tensor_Dim1>
  FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                   Tensor_Dim0, Tensor_Dim1>
  getFTensor2N(FieldApproximationBase base, const int gg, const int bb);

  /** \brief Get base functions for Hdiv space

  Example:
  \code
  FTensor::Index<'i',3> i;
  FTensor::Index<'j',3> j;
  int nb_dofs = data.getFieldData().size();
  for(int gg = 0;gg!=nb_gauss_pts;++gg) {
    int ll = 0;
    auto t_n_hdiv = data.getFTensor2N<3,3>(gg,0);
    for(;ll!=nb_dofs;ll++) {
      double dot_product = t_n_hdiv(i)*t_n_hdiv(i);
      ++t_n_hdiv;
    }
    for(;ll!=data.getVectorN().size2()/3;ll++) {
      ++t_n_hdiv;
    }
  }
  \endcode

  */
  template <int Tensor_Dim0, int Tensor_Dim1>
  auto getFTensor2N(const int gg, const int bb);

  /**@}*/

  /** \name Auxiliary functions */

  /**@{*/

  friend std::ostream &operator<<(std::ostream &os,
                                  const EntitiesFieldData::EntData &e);

  /**
   * Reset data associated with particular field name
   * @return error code
   */
  inline MoFEMErrorCode resetFieldDependentData();

  /**@}*/

  /** \name Bernstein-Bezier base only functions */

  /**@{*/

  /**
   * @brief Get orders at the nodes
   *
   * @return VectorInt&
   */
  inline VectorInt &getBBNodeOrder();

  /**
   * @brief Get file BB indices
   *
   * @return MatrixInt&
   */
  inline MatrixInt &getBBAlphaIndices();

  virtual boost::shared_ptr<MatrixInt> &
  getBBAlphaIndicesSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to BB base base functions
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getBBNSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to BB base base functions
   */
  virtual const boost::shared_ptr<MatrixDouble> &
  getBBNSharedPtr(const std::string &field_name) const;

  /**
   * Get shared pointer to BB derivatives of base base functions
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getBBDiffNSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to derivatives of BB base base functions
   */
  virtual const boost::shared_ptr<MatrixDouble> &
  getBBDiffNSharedPtr(const std::string &field_name) const;

  virtual std::map<std::string, boost::shared_ptr<MatrixInt>> &
  getBBAlphaIndicesMap();

  /**
   * @brief get hash map of base function for BB base, key is a field name
   * 
   * @return std::map<std::string, boost::shared_ptr<MatrixDouble>>& 
   */
  virtual std::map<std::string, boost::shared_ptr<MatrixDouble>> &getBBNMap();

  /**
   * @brief get hash map of direvarives base function for BB base, key is a
   * field name
   *
   * @return std::map<std::string, boost::shared_ptr<MatrixDouble>>&
   */
  virtual std::map<std::string, boost::shared_ptr<MatrixDouble>> &
  getBBDiffNMap();

  /**
   * @brief get ALpha indices for BB base by order
   * 
   * @param o approximation order
   * @return boost::shared_ptr<MatrixInt>& 
   */
  virtual boost::shared_ptr<MatrixInt> &
  getBBAlphaIndicesByOrderSharedPtr(const size_t o);

  /**
   * @brief get BB base by order
   * 
   * @param o 
   * @return boost::shared_ptr<MatrixDouble>& 
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getBBNByOrderSharedPtr(const size_t o);

  /**
   * @brief get BB base derivative by order
   * 
   * @param o 
   * @return boost::shared_ptr<MatrixDouble>& 
   */
  virtual boost::shared_ptr<MatrixDouble> &
  getBBDiffNByOrderSharedPtr(const size_t o);

  static constexpr size_t MaxBernsteinBezierOrder = BITFEID_SIZE;

  virtual std::array<boost::shared_ptr<MatrixInt>, MaxBernsteinBezierOrder> &
  getBBAlphaIndicesByOrderArray();

  virtual std::array<boost::shared_ptr<MatrixDouble>, MaxBernsteinBezierOrder> &
  getBBNByOrderArray();

  virtual std::array<boost::shared_ptr<MatrixDouble>, MaxBernsteinBezierOrder> &
  getBBDiffNByOrderArray();

  /**
   * @brief Swap bases functions
   *
   * Some base are not hierarchical and depene on approximation order. Such case
   * demand special handling, that appropiate base order is set depending on
   * field, such that is accessible in operator.
   *
   * @note Base is not swap on meshsets
   *
   * @param field_name
   * @param base
   * @return MoFEMErrorCode
   */
  virtual MoFEMErrorCode baseSwap(const std::string &field_name,
                                  const FieldApproximationBase base);

  /**@}*/

protected:
  int sEnse;                         ///< Entity sense (orientation)
  ApproximationOrder oRder;          ///< Entity order
  FieldSpace sPace;                  ///< Entity space
  FieldApproximationBase bAse;       ///< Field approximation base
  BitRefLevel entDataBitRefLevel;    ///< Bit ref level in entity
  VectorInt iNdices;                 ///< Global indices on entity
  VectorInt localIndices;            ///< Local indices on entity
  VectorDofs dOfs;                   ///< DoFs on entity
  VectorFieldEntities fieldEntities; ///< Field entities
  VectorDouble fieldData;            ///< Field data on entity

  std::array<std::array<boost::shared_ptr<MatrixDouble>, LASTBASE>,
             LastDerivative>
      baseFunctionsAndBaseDerivatives;

  std::array<boost::shared_ptr<MatrixDouble>, LASTBASE> &N; ///< Base functions
  std::array<boost::shared_ptr<MatrixDouble>, LASTBASE>
      &diffN; ///< Derivatives of base functions

  std::string bbFieldName; ///< field name
  VectorInt bbNodeOrder;   ///< order of nodes
  std::map<std::string, boost::shared_ptr<MatrixDouble>> bbN;
  std::map<std::string, boost::shared_ptr<MatrixDouble>> bbDiffN;
  std::map<std::string, boost::shared_ptr<MatrixInt>>
      bbAlphaIndices; ///< Indices for Bernstein-Bezier (BB) base

  std::array<boost::shared_ptr<MatrixDouble>, MaxBernsteinBezierOrder>
      bbNByOrder; ///< BB base functions by order
  std::array<boost::shared_ptr<MatrixDouble>, MaxBernsteinBezierOrder>
      bbDiffNByOrder; ///< BB base functions direvatives by order
  std::array<boost::shared_ptr<MatrixInt>, MaxBernsteinBezierOrder>
      bbAlphaIndicesByOrder; ///< BB alpha indices by order

protected:
  /**
   * @brief Used by Bernstein base to keep temporally pointer
   *
   * @copydoc MoFEM::EntitiesFieldData::baseSwap
   */
  boost::shared_ptr<MatrixDouble> swapBaseNPtr;

  /**
   * @brief Used by Bernstein base to keep temporally pointer
   *
   * @copydoc MoFEM::EntitiesFieldData::baseSwap
   */
  boost::shared_ptr<MatrixDouble> swapBaseDiffNPtr;

  friend struct OpAddParentEntData;
};

using BaseDerivatives = EntitiesFieldData::EntData::BaseDerivatives;

/** \brief Derived ata on single entity (This is passed as argument to
 * DataOperator::doWork) \ingroup mofem_forces_and_sources_user_data_operators
 * \nosubgrouping
 *
 * DerivedEntData share part information with EntData except infomation about
 * base functions.
 *
 */
struct DerivedEntitiesFieldData::DerivedEntData
    : public EntitiesFieldData::EntData {

  DerivedEntData(
      const boost::shared_ptr<EntitiesFieldData::EntData> &ent_data_ptr);

  int getSense() const;

  //// \brief get entity bit ref level
  BitRefLevel &getEntDataBitRefLevel();

  boost::shared_ptr<MatrixDouble> &
  getNSharedPtr(const FieldApproximationBase base,
                const BaseDerivatives derivative);

  boost::shared_ptr<MatrixDouble> &
  getNSharedPtr(const FieldApproximationBase base);

  boost::shared_ptr<MatrixDouble> &
  getDiffNSharedPtr(const FieldApproximationBase base);

  const boost::shared_ptr<MatrixDouble> &
  getNSharedPtr(const FieldApproximationBase base) const;

  const boost::shared_ptr<MatrixDouble> &
  getDiffNSharedPtr(const FieldApproximationBase base) const;

  inline boost::shared_ptr<MatrixDouble> &
  getDerivedNSharedPtr(const FieldApproximationBase base);

  inline boost::shared_ptr<MatrixDouble> &
  getDerivedDiffNSharedPtr(const FieldApproximationBase base);

  boost::shared_ptr<MatrixInt> &
  getBBAlphaIndicesSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to BB base base functions
   */
  boost::shared_ptr<MatrixDouble> &
  getBBNSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to BB base base functions
   */
  const boost::shared_ptr<MatrixDouble> &
  getBBNSharedPtr(const std::string &field_name) const;

  /**
   * Get shared pointer to BB derivatives of base base functions
   */
  boost::shared_ptr<MatrixDouble> &
  getBBDiffNSharedPtr(const std::string &field_name);

  /**
   * Get shared pointer to derivatives of BB base base functions
   */
  const boost::shared_ptr<MatrixDouble> &
  getBBDiffNSharedPtr(const std::string &field_name) const;

  /**
   * @copydoc MoFEM::EntitiesFieldData::EntData::swapBaseNPtr
   */
  MoFEMErrorCode baseSwap(const std::string &field_name,
                          const FieldApproximationBase base);

protected:
  const boost::shared_ptr<EntitiesFieldData::EntData> entDataPtr;
};

ApproximationOrder EntitiesFieldData::EntData::getOrder() const {
  return oRder;
}

const VectorInt &EntitiesFieldData::EntData::getIndices() const {
  return iNdices;
}

const VectorIntAdaptor
EntitiesFieldData::EntData::getIndicesUpToOrder(int order) {
  unsigned int size = 0;
  if (auto dof = dOfs[0]) {
    size = dof->getOrderNbDofs(order) * dof->getNbOfCoeffs();
    size = size < iNdices.size() ? size : iNdices.size();
  }
  int *data = &*iNdices.data().begin();
  return VectorIntAdaptor(size, ublas::shallow_array_adaptor<int>(size, data));
}

const VectorInt &EntitiesFieldData::EntData::getLocalIndices() const {
  return localIndices;
}

const VectorIntAdaptor
EntitiesFieldData::EntData::getLocalIndicesUpToOrder(int order) {
  unsigned int size = 0;
  if (auto dof = dOfs[0]) {
    size = dof->getOrderNbDofs(order) * dof->getNbOfCoeffs();
    size = size < localIndices.size() ? size : localIndices.size();
  }
  int *data = &*localIndices.data().begin();
  return VectorIntAdaptor(size, ublas::shallow_array_adaptor<int>(size, data));
}

int &EntitiesFieldData::EntData::getSense() { return sEnse; }

ApproximationOrder &EntitiesFieldData::EntData::getOrder() { return oRder; }

VectorInt &EntitiesFieldData::EntData::getIndices() { return iNdices; }

VectorInt &EntitiesFieldData::EntData::getLocalIndices() {
  return localIndices;
}

const VectorDouble &EntitiesFieldData::EntData::getFieldData() const {
  return fieldData;
}

const VectorAdaptor
EntitiesFieldData::EntData::getFieldDataUpToOrder(int order) {
  unsigned int size = 0;
  if (auto dof = dOfs[0]) {
    size = dof->getOrderNbDofs(order) * dof->getNbOfCoeffs();
    size = size < fieldData.size() ? size : fieldData.size();
  }
  double *data = &*fieldData.data().begin();
  return getVectorAdaptor(data, size);
}

const VectorDofs &EntitiesFieldData::EntData::getFieldDofs() const {
  return dOfs;
}

VectorDofs &EntitiesFieldData::EntData::getFieldDofs() { return dOfs; }

VectorDouble &EntitiesFieldData::EntData::getFieldData() { return fieldData; }

VectorFieldEntities &EntitiesFieldData::EntData::getFieldEntities() {
  return fieldEntities;
}

const VectorFieldEntities &
EntitiesFieldData::EntData::getFieldEntities() const {
  return fieldEntities;
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1FieldData() {
  std::stringstream s;
  s << "Not implemented for this dimension dim = " << Tensor_Dim;
  THROW_MESSAGE(s.str());
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
EntitiesFieldData::EntData::getFTensor2FieldData() {
  std::stringstream s;
  s << "Not implemented for this dimension dim0 = " << Tensor_Dim0;
  s << " and dim1 " << Tensor_Dim1;
  THROW_MESSAGE(s.str());
}

template <int Tensor_Dim>
FTensor::Tensor2_symmetric<
    FTensor::PackPtr<double *, (Tensor_Dim * (Tensor_Dim + 1)) / 2>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor2SymmetricFieldData() {
  std::stringstream s;
  s << "Not implemented for this dimension dim = " << Tensor_Dim;
  THROW_MESSAGE(s.str());
}

FieldApproximationBase &EntitiesFieldData::EntData::getBase() { return bAse; }

FieldSpace &EntitiesFieldData::EntData::getSpace() { return sPace; }

MatrixDouble &
EntitiesFieldData::EntData::getN(const FieldApproximationBase base) {
  return *(getNSharedPtr(base));
}

MatrixDouble &EntitiesFieldData::EntData::getN(const std::string &field_name) {
  return *(getBBNSharedPtr(field_name));
}

MatrixDouble &EntitiesFieldData::EntData::getN() { return getN(bAse); }

MatrixDouble &
EntitiesFieldData::EntData::getDiffN(const FieldApproximationBase base) {
  return *(getDiffNSharedPtr(base));
}

MatrixDouble &
EntitiesFieldData::EntData::getN(const FieldApproximationBase base,
                                 const BaseDerivatives derivative) {
#ifndef NDEBUG
  if (!getNSharedPtr(base, derivative)) {
    MOFEM_LOG_C("SELF", Sev::error,
                "Ptr to base %s functions derivative %d is null",
                ApproximationBaseNames[base], derivative);
    THROW_MESSAGE("Null pointer");
  }
#endif
  return *(getNSharedPtr(base, derivative));
}

MatrixDouble &
EntitiesFieldData::EntData::getDiffN(const std::string &field_name) {
  return *(getBBDiffNSharedPtr(field_name));
}

MatrixDouble &EntitiesFieldData::EntData::getDiffN() { return getDiffN(bAse); }

MatrixDouble &
EntitiesFieldData::EntData::getN(const BaseDerivatives derivative) {
  return getN(bAse, derivative);
}

const VectorAdaptor
EntitiesFieldData::EntData::getN(const FieldApproximationBase base,
                                 const int gg) {
  int size = getN(base).size2();
  double *data = &getN(base)(gg, 0);
  return VectorAdaptor(size, ublas::shallow_array_adaptor<double>(size, data));
}

const VectorAdaptor EntitiesFieldData::EntData::getN(const int gg) {
  return getN(bAse, gg);
}

const MatrixAdaptor
EntitiesFieldData::EntData::getDiffN(const FieldApproximationBase base,
                                     const int gg) {
  // FIXME: That is bug, it will not work if number of integration pts is
  // equal to number of nodes on entity.  User who not implementing low
  // level DataOperator will not experience this.
  if (getN(base).size1() == getDiffN(base).size1()) {
    int size = getN(base).size2();
    int dim = getDiffN(base).size2() / size;
    double *data = &getDiffN(base)(gg, 0);
    return MatrixAdaptor(
        getN(base).size2(), dim,
        ublas::shallow_array_adaptor<double>(getDiffN(base).size2(), data));
  } else {
    // in some cases, f.e. for derivatives of nodal base functions at only
    // one gauss point is needed
    return MatrixAdaptor(
        getN(base).size1(), getN(base).size2(),
        ublas::shallow_array_adaptor<double>(getDiffN(base).data().size(),
                                             &getDiffN(base).data()[0]));
  }
}

const MatrixAdaptor EntitiesFieldData::EntData::getDiffN(const int gg) {
  return getDiffN(bAse, gg);
}

const VectorAdaptor
EntitiesFieldData::EntData::getN(const FieldApproximationBase base,
                                 const int gg, const int nb_base_functions) {
  (void)getN()(gg, nb_base_functions -
                       1); // throw error if nb_base_functions is to big
  double *data = &getN(base)(gg, 0);
  return VectorAdaptor(nb_base_functions, ublas::shallow_array_adaptor<double>(
                                              nb_base_functions, data));
}

const VectorAdaptor
EntitiesFieldData::EntData::getN(const int gg, const int nb_base_functions) {
  return getN(bAse, gg, nb_base_functions);
}

const MatrixAdaptor
EntitiesFieldData::EntData::getDiffN(const FieldApproximationBase base,
                                     const int gg,
                                     const int nb_base_functions) {
  // FIXME: That is bug, it will not work if number of integration pts is
  // equal to number of nodes on entity.  User who not implementing low
  // level DataOperator will not experience this.
  if (getN(base).size1() == getDiffN(base).size1()) {
    (void)getN(base)(gg,
                     nb_base_functions -
                         1); // throw error if nb_base_functions is to big
    int dim = getDiffN(base).size2() / getN(base).size2();
    double *data = &getDiffN(base)(gg, 0);
    return MatrixAdaptor(
        nb_base_functions, dim,
        ublas::shallow_array_adaptor<double>(dim * nb_base_functions, data));
  } else {
    // in some cases, f.e. for derivatives of nodal base functions only one
    // gauss point is needed
    return MatrixAdaptor(
        getN(base).size1(), getN(base).size2(),
        ublas::shallow_array_adaptor<double>(getDiffN(base).data().size(),
                                             &getDiffN(base).data()[0]));
  }
}

const MatrixAdaptor
EntitiesFieldData::EntData::getDiffN(const int gg,
                                     const int nb_base_functions) {
  return getDiffN(bAse, gg, nb_base_functions);
}

template <int DIM>
const MatrixAdaptor
EntitiesFieldData::EntData::getVectorN(const FieldApproximationBase base,
                                       const int gg) {
  if (PetscUnlikely(getN(base).size2() % DIM)) {
    THROW_MESSAGE("Wrong dimension");
  }

  const int nb_base_functions = getN(base).size2() / DIM;
  double *data = &getN(base)(gg, 0);
  return MatrixAdaptor(
      nb_base_functions, DIM,
      ublas::shallow_array_adaptor<double>(DIM * nb_base_functions, data));
}

template <int DIM>
const MatrixAdaptor EntitiesFieldData::EntData::getVectorN(const int gg) {
  return getVectorN<DIM>(bAse, gg);
}

template <int DIM0, int DIM1>
const MatrixAdaptor
EntitiesFieldData::EntData::getVectorDiffN(FieldApproximationBase base,
                                           const int gg) {
  if (PetscUnlikely(getDiffN(base).size2() % (DIM0 * DIM1))) {
    THROW_MESSAGE("Wrong dimension");
  }

  const int nb_base_functions = getN(base).size2() / (DIM0 * DIM1);
  double *data = &getN(base)(gg, 0);
  return MatrixAdaptor(nb_base_functions, DIM0 * DIM1,
                       ublas::shallow_array_adaptor<double>(
                           DIM0 * DIM1 * nb_base_functions, data));
}

template <int DIM0, int DIM1>
const MatrixAdaptor EntitiesFieldData::EntData::getVectorDiffN(const int gg) {
  return getVectorDiffN<DIM0, DIM1>(bAse, gg);
}

template <int DIM0, int DIM1>
const MatrixAdaptor
EntitiesFieldData::EntData::getVectorDiffN(const FieldApproximationBase base,
                                           const int dof, const int gg) {
  double *data =
      &EntitiesFieldData::EntData::getDiffN(base)(gg, DIM0 * DIM1 * dof);
  return MatrixAdaptor(DIM0, DIM1,
                       ublas::shallow_array_adaptor<double>(DIM0 * DIM1, data));
}

template <int DIM0, int DIM1>
const MatrixAdaptor EntitiesFieldData::EntData::getVectorDiffN(const int dof,
                                                               const int gg) {
  return getVectorDiffN<DIM0, DIM1>(bAse, dof, gg);
}

FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
EntitiesFieldData::EntData::getFTensor0N(const FieldApproximationBase base) {
  double *ptr = &*getN(base).data().begin();
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(ptr);
};

FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
EntitiesFieldData::EntData::getFTensor0N() {
  return getFTensor0N(bAse);
};

FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
EntitiesFieldData::EntData::getFTensor0N(const FieldApproximationBase base,
                                         const int gg, const int bb) {
  double *ptr = &getN(base)(gg, bb);
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(ptr);
};

FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
EntitiesFieldData::EntData::getFTensor0N(const int gg, const int bb) {
  return getFTensor0N(bAse, gg, bb);
};

template <int Tensor_Dim> auto EntitiesFieldData::EntData::getFTensor1N() {
  return getFTensor1N<Tensor_Dim>(bAse);
}

template <int Tensor_Dim>
auto EntitiesFieldData::EntData::getFTensor1N(const int gg, const int bb) {
  return getFTensor1N<Tensor_Dim>(bAse, gg, bb);
}

template <int Tensor_Dim0, int Tensor_Dim1>
auto EntitiesFieldData::EntData::getFTensor2N() {
  return getFTensor2N<Tensor_Dim0, Tensor_Dim1>(bAse);
}

template <int Tensor_Dim0, int Tensor_Dim1>
auto EntitiesFieldData::EntData::getFTensor2N(const int gg, const int bb) {
  return getFTensor2N<Tensor_Dim0, Tensor_Dim1>(bAse, gg, bb);
}

/** \name Bernstein-Bezier base only functions */

/**@{*/

VectorInt &EntitiesFieldData::EntData::getBBNodeOrder() { return bbNodeOrder; }

MatrixInt &EntitiesFieldData::EntData::getBBAlphaIndices() {
  return *getBBAlphaIndicesSharedPtr(bbFieldName);
}

/**@}*/

/** \name DerivedEntData */

/**@{*/

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getDerivedNSharedPtr(
    const FieldApproximationBase base) {
  return N[base];
}

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getDerivedDiffNSharedPtr(
    const FieldApproximationBase base) {
  return diffN[base];
}

/**@}*/

/**
 * @brief Assemble PETSc vector
 *
 * Function extract indices from entity data and assemble vector
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Vec/VecSetValues.html>See
 * PETSc documentation</a>
 *
 * @param V
 * @param data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <typename T = EntityStorage>
inline MoFEMErrorCode VecSetValues(Vec V,
                                   const EntitiesFieldData::EntData &data,
                                   const double *ptr, InsertMode iora) {
  static_assert(!std::is_same<T, T>::value,
                "VecSetValues value for this data storage is not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template <>
inline MoFEMErrorCode
VecSetValues<EntityStorage>(Vec V, const EntitiesFieldData::EntData &data,
                            const double *ptr, InsertMode iora) {
  return VecSetValues(V, data.getIndices().size(), &*data.getIndices().begin(),
                      ptr, iora);
}

/**
 * @brief Assemble PETSc matrix
 *
 * Function extract indices from entity data and assemble vector
 *
 * <a
 * href=https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatSetValues.html>See
 * PETSc documentation</a>
 *
 * @param M
 * @param row_data
 * @param col_data
 * @param ptr
 * @param iora
 * @return MoFEMErrorCode
 */
template <typename T = EntityStorage>
inline MoFEMErrorCode MatSetValues(Mat M,
                                   const EntitiesFieldData::EntData &row_data,
                                   const EntitiesFieldData::EntData &col_data,
                                   const double *ptr, InsertMode iora) {
  static_assert(!std::is_same<T, T>::value,
                "MatSetValues value for this data storage is not implemented");
  return MOFEM_NOT_IMPLEMENTED;
}

template <>
inline MoFEMErrorCode
MatSetValues<EntityStorage>(Mat M, const EntitiesFieldData::EntData &row_data,
                            const EntitiesFieldData::EntData &col_data,
                            const double *ptr, InsertMode iora) {
  return MatSetValues(
      M, row_data.getIndices().size(), &*row_data.getIndices().begin(),
      col_data.getIndices().size(), &*col_data.getIndices().begin(), ptr, iora);
}

/** \name Specializations for tensor base function */

/**@{*/

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1N<3>(FieldApproximationBase base);

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1N<3>(FieldApproximationBase base,
                                            const int gg, const int bb);

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2N<3, 3>(FieldApproximationBase base);

/**@}*/

/** \name Specializations for direcatives of base functions */

/**@{*/

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base);
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>();

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base);
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>();

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 2>(FieldApproximationBase base);
template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 2>(FieldApproximationBase base,
                                                   const int gg, const int bb);

template <>
FTensor::Tensor3<FTensor::PackPtr<double *, 12>, 3, 2, 2>
EntitiesFieldData::EntData::getFTensor3Diff2N(FieldApproximationBase base);

/**@}*/

/** \name Specializations for field data */

/**@{*/

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1FieldData<3>();

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1FieldData<2>();

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 1>
EntitiesFieldData::EntData::getFTensor1FieldData<1>();

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2FieldData<3, 3>();

template <>
FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>
EntitiesFieldData::EntData::getFTensor2SymmetricFieldData<3>();

template <>
FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>
EntitiesFieldData::EntData::getFTensor2SymmetricFieldData<2>();

/**@}*/

/**
 * @deprecated Use EntitiesFieldData
 */
DEPRECATED typedef EntitiesFieldData DataForcesAndSourcesCore;

/**
 * @deprecated Use DerivedEntitiesFieldData
 */
DEPRECATED typedef DerivedEntitiesFieldData DerivedDataForcesAndSourcesCore;

} // namespace MoFEM

#endif //__ENTITIES_FIELD_DATA_HPP__

/**
 * \defgroup mofem_forces_and_sources_user_data_operators User data operator
 * data structures \ingroup
 *
 * \brief Users data structures and operator
 *
 * Data structures passed by argument to MoFEM::DataOperator::doWork and generic
 * user operators operating on those structures.
 *
 */
