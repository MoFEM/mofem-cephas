/** \file ForcesAndSourcesCore.hpp

  \brief Implementation of elements on entities.

  Those element are inherited by user to implement specific implementation of
  particular problem.

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

#ifndef __FATPRISMELEMENTFORCESANDSURCESCORE_HPP__
#define __FATPRISMELEMENTFORCESANDSURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief FatPrism finite element
 \ingroup mofem_forces_and_sources_prism_element

 User is implementing own operator at Gauss points level, by own object
 derived from FatPrismElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 \todo Need to implement operators that will make this element work as Volume
 element

 */
struct FatPrismElementForcesAndSourcesCore
    : public VolumeElementForcesAndSourcesCore {

  FatPrismElementForcesAndSourcesCore(Interface &m_field);

  virtual int getRuleTrianglesOnly(int order) { return 2 * order; };
  virtual int getRuleThroughThickness(int order) { return 2 * order; };

  virtual MoFEMErrorCode setGaussPtsTrianglesOnly(int order_triangles_only) {
    MoFEMFunctionBeginHot;
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    MoFEMFunctionReturnHot(0);
  }

  virtual MoFEMErrorCode setGaussPtsThroughThickness(int order_thickness) {
    MoFEMFunctionBeginHot;
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
    MoFEMFunctionReturnHot(0);
  }

  /** \brief default operator for Flat Prism element
   * \ingroup mofem_forces_and_sources_prism_element
   */
  struct UserDataOperator
      : public VolumeElementForcesAndSourcesCore::UserDataOperator {

    using VolumeElementForcesAndSourcesCore::UserDataOperator::UserDataOperator;

    /** \brief get face aRea
    \param dd if dd == 0 it is for face F3 if dd == 1 is for face F4
    */
    inline double getArea(const int dd);

    inline double getAreaF3();
    inline double getAreaF4();

    /** \brief get triangle normal

    Normal has 6 elements, first 3 are for face F3 another three for face F4

     */
    inline VectorDouble &getNormal();

    inline VectorAdaptor getNormalF3();

    inline VectorAdaptor getNormalF4();
    
    /** \brief get Gauss pts. in the prism
     */
    inline MatrixDouble &getGaussPts();

    /** \brief get Gauss pts. on triangles
     */
    inline MatrixDouble &getGaussPtsTrianglesOnly();

    /** \brief get Gauss pts. through thickness
     */
    inline MatrixDouble &getGaussPtsThroughThickness();

    /** \brief get coordinates at Gauss pts.

      Matrix has size (nb integration points)x(coordinates on F3 and F4 = 6),
      i.e. coordinates on face F3 and F4

     */
    inline MatrixDouble &getCoordsAtGaussPts();

    /** \brief get coordinates at Gauss pts.

      Matrix has size (nb integration points)x(coordinates on F3 and F4 = 6),
      i.e. coordinates on face F3 and F4

     */
    inline MatrixDouble &getCoordsAtGaussPtsTrianglesOnly();

    /** \brief coordinate at Gauss points on face 3 (if hierarchical
     * approximation of element geometry)
     */
    inline MatrixDouble &getHOCoordsAtGaussPtsF3();

    /** \brief coordinate at Gauss points on face 4 (if hierarchical
     * approximation of element geometry)
     */
    inline MatrixDouble &getHOCoordsAtGaussPtsF4();

    /** \brief if higher order geometry return normals at face F3 at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see \cite
     * tautges2010canonical
     *
     */
    inline MatrixDouble &getNormalsAtGaussPtF3();

    /** \brief if higher order geometry return normals at face F4 at Gauss pts.
     *
     * Face 4 is top face in canonical triangle numeration, see \cite
     * tautges2010canonical
     *
     */
    inline MatrixDouble &getNormalsAtGaussPtF4();

    /** \brief if higher order geometry return normals at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see \cite
     * tautges2010canonical
     *
     * \param gg gauss point number
     */
    inline ublas::matrix_row<MatrixDouble> getNormalsAtGaussPtF3(const int gg);

    /** \brief if higher order geometry return normals at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see \cite
     * tautges2010canonical
     *
     * \param gg gauss point number
     */
    inline ublas::matrix_row<MatrixDouble> getNormalsAtGaussPtF4(const int gg);

    /** \brief if higher order geometry return tangent vector to triangle at
     * Gauss pts.
     */
    inline MatrixDouble &getTangent1AtGaussPtF3();

    /** \brief if higher order geometry return tangent vector to triangle at
     * Gauss pts.
     */
    inline MatrixDouble &getTangent2AtGaussPtF3();

    /** \brief if higher order geometry return tangent vector to triangle at
     * Gauss pts.
     */
    inline MatrixDouble &getTangent1AtGaussPtF4();

    /** \brief if higher order geometry return tangent vector to triangle at
     * Gauss pts.
     */
    inline MatrixDouble &getTangent2AtGaussPtF4();

    inline DataForcesAndSourcesCore &getTrianglesOnlyDataStructure();

    inline DataForcesAndSourcesCore &getTroughThicknessDataStructure();

    /** \brief return pointer to fat prism finite element
     */
    inline const FatPrismElementForcesAndSourcesCore *getPrismFE();

  protected:
    MoFEMErrorCode setPtrFE(ForcesAndSourcesCore *ptr);

  };

  MoFEMErrorCode operator()();

protected:
  double aRea[2];
  VectorDouble normal;

  MatrixDouble gaussPtsTrianglesOnly;
  MatrixDouble coordsAtGaussPtsTrianglesOnly;
  MatrixDouble gaussPtsThroughThickness;

  DataForcesAndSourcesCore dataH1TrianglesOnly;
  DataForcesAndSourcesCore dataH1TroughThickness;

  MatrixDouble hoCoordsAtGaussPtsF3;
  MatrixDouble nOrmals_at_GaussPtF3;
  MatrixDouble tAngent1_at_GaussPtF3;
  MatrixDouble tAngent2_at_GaussPtF3;
  MatrixDouble hoCoordsAtGaussPtsF4;
  MatrixDouble nOrmals_at_GaussPtF4;
  MatrixDouble tAngent1_at_GaussPtF4;
  MatrixDouble tAngent2_at_GaussPtF4;
  OpGetCoordsAndNormalsOnPrism opHOCoordsAndNormals;

  friend class UserDataOperator;
};

inline double
FatPrismElementForcesAndSourcesCore::UserDataOperator::getArea(const int dd) {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->aRea[0];
}

inline double
FatPrismElementForcesAndSourcesCore::UserDataOperator::getAreaF3() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->aRea[0];
}
inline double
FatPrismElementForcesAndSourcesCore::UserDataOperator::getAreaF4() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->aRea[1];
}

inline VectorDouble &
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormal() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->normal;
}

inline VectorAdaptor
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalF3() {
  double *data =
      &(static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->normal[0]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline VectorAdaptor
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalF4() {
  double *data =
      &(static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->normal[3]);
  return VectorAdaptor(3, ublas::shallow_array_adaptor<double>(3, data));
}

inline MatrixDouble &
FatPrismElementForcesAndSourcesCore::UserDataOperator::getGaussPts() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)->gaussPts;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getGaussPtsTrianglesOnly() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->gaussPtsTrianglesOnly;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getGaussPtsThroughThickness() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->gaussPtsThroughThickness;
}

inline MatrixDouble &
FatPrismElementForcesAndSourcesCore::UserDataOperator::getCoordsAtGaussPts() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->coordsAtGaussPts;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getCoordsAtGaussPtsTrianglesOnly() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->coordsAtGaussPtsTrianglesOnly;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getHOCoordsAtGaussPtsF3() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->hoCoordsAtGaussPtsF3;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getHOCoordsAtGaussPtsF4() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->hoCoordsAtGaussPtsF4;
}

inline MatrixDouble &
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPtF3() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->nOrmals_at_GaussPtF3;
}

inline MatrixDouble &
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPtF4() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->nOrmals_at_GaussPtF4;
}

inline ublas::matrix_row<MatrixDouble>
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPtF3(
    const int gg) {
  return ublas::matrix_row<MatrixDouble>(
      static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
          ->nOrmals_at_GaussPtF3,
      gg);
}

inline ublas::matrix_row<MatrixDouble>
FatPrismElementForcesAndSourcesCore::UserDataOperator::getNormalsAtGaussPtF4(
    const int gg) {
  return ublas::matrix_row<MatrixDouble>(
      static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
          ->nOrmals_at_GaussPtF4,
      gg);
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangent1AtGaussPtF3() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->tAngent1_at_GaussPtF3;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangent2AtGaussPtF3() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->tAngent2_at_GaussPtF3;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangent1AtGaussPtF4() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->tAngent1_at_GaussPtF4;
}

inline MatrixDouble &FatPrismElementForcesAndSourcesCore::UserDataOperator::
    getTangent2AtGaussPtF4() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->tAngent2_at_GaussPtF4;
}

inline DataForcesAndSourcesCore &FatPrismElementForcesAndSourcesCore::
    UserDataOperator::getTrianglesOnlyDataStructure() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->dataH1TrianglesOnly;
}

inline DataForcesAndSourcesCore &FatPrismElementForcesAndSourcesCore::
    UserDataOperator::getTroughThicknessDataStructure() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE)
      ->dataH1TroughThickness;
}

inline const FatPrismElementForcesAndSourcesCore *
FatPrismElementForcesAndSourcesCore::UserDataOperator::getPrismFE() {
  return static_cast<FatPrismElementForcesAndSourcesCore *>(ptrFE);
}

} // namespace MoFEM

#endif //__FATPRISMELEMENTFORCESANDSURCESCORE_HPP__
