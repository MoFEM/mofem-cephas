/** \file ElementsOnEntities.hpp

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

using namespace boost::numeric;

#ifndef __FLATPRISMELEMENTFORCESANDSURCESCORE_HPP__
#define __FLATPRISMELEMENTFORCESANDSURCESCORE_HPP__

namespace MoFEM {

/** \brief FlatPrism finite element
 \ingroup mofem_forces_and_sources_prism_element

 User is implementing own operator at Gauss points level, by own object
 derived from FlatPrismElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct FlatPrismElementForcesAndSurcesCore: public ForcesAndSurcesCore {

  MoABErrorCode rval;
  double aRea[2];
  VectorDouble normal;
  VectorDouble coords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;

  std::string meshPositionsFieldName;

  MatrixDouble hoCoordsAtGaussPtsF3;
  MatrixDouble nOrmals_at_GaussPtF3;
  MatrixDouble tAngent1_at_GaussPtF3;
  MatrixDouble tAngent2_at_GaussPtF3;
  MatrixDouble hoCoordsAtGaussPtsF4;
  MatrixDouble nOrmals_at_GaussPtF4;
  MatrixDouble tAngent1_at_GaussPtF4;
  MatrixDouble tAngent2_at_GaussPtF4;
  OpGetCoordsAndNormalsOnPrism opHOCoordsAndNormals;

  FlatPrismElementForcesAndSurcesCore(FieldInterface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBPRISM),derivedDataH1(dataH1),
    dataHdiv(MBPRISM),derivedDataHdiv(dataHdiv),
    dataNoField(MBPRISM),dataNoFieldCol(MBPRISM),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHOCoordsAndNormals(
      hoCoordsAtGaussPtsF3,nOrmals_at_GaussPtF3,tAngent1_at_GaussPtF3,tAngent2_at_GaussPtF3,
      hoCoordsAtGaussPtsF4,nOrmals_at_GaussPtF4,tAngent1_at_GaussPtF4,tAngent2_at_GaussPtF4
    ) {
    }

  /** \brief default operator for Flat Prism element
    * \ingroup mofem_forces_and_sources_prism_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(const std::string &field_name,const char type):
    ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type
    ):
    ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {
    }

    /** \brief get face aRea
    \param dd if dd == 0 it is for face F3 if dd == 1 is for face F4
    */
    inline double getArea(const int dd) {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->aRea[0];
    }

    inline double getAreaF3() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->aRea[0];
    }
    inline double getAreaF4() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->aRea[1];
    }

    /** \brief get triangle normal

    Normal has 6 elements, first 3 are for face F3 another three for face F4

     */
    inline VectorDouble& getNormal() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->normal;
    }

    inline VectorAdaptor getNormalF3() {
      double *data  = &(static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->normal[0]);
      return VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,data));
    }

    inline VectorAdaptor getNormalF4() {
      double *data  = &(static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->normal[3]);
      return VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,data));
    }

    /** \brief get triangle coordinates

      Vector has 6 elements, i.e. coordinates on face F3 and F4

     */
    inline VectorDouble& getCoords() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->coords;
    }

    /** \brief get triangle Gauss pts.
     */
    inline MatrixDouble& getGaussPts() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->gaussPts;
    }

    /** \brief get coordinates at Gauss pts.

      Matrix has size (nb integration points)x(coordinates on F3 and F4 = 6), i.e. coordinates on face F3 and F4

     */
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->coordsAtGaussPts;
    }

    /** \brief coordinate at Gauss points on face 3 (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPtsF3() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->hoCoordsAtGaussPtsF3;
    }

    /** \brief coordinate at Gauss points on face 4 (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPtsF4() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->hoCoordsAtGaussPtsF4;
    }

    /** \brief if higher order geometry return normals at face F3 at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
     *
     */
    inline MatrixDouble& getNormals_at_GaussPtF3() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->nOrmals_at_GaussPtF3;
    }

    /** \brief if higher order geometry return normals at face F4 at Gauss pts.
     *
     * Face 4 is top face in canonical triangle numeration, see \cite tautges2010canonical
     *
     */
    inline MatrixDouble& getNormals_at_GaussPtF4() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->nOrmals_at_GaussPtF4;
    }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF3(const int gg) {
      return ublas::matrix_row<MatrixDouble >(static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->nOrmals_at_GaussPtF3,gg);
    }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF4(const int gg) {
      return ublas::matrix_row<MatrixDouble >(static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->nOrmals_at_GaussPtF4,gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF3() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->tAngent1_at_GaussPtF3;
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF3() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->tAngent2_at_GaussPtF3;
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF4() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->tAngent1_at_GaussPtF4;
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF4() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE)->tAngent2_at_GaussPtF4;
    }

    /** \brief return pointer to triangle finite element object
     */
    inline const FlatPrismElementForcesAndSurcesCore* getFlatPrismElementForcesAndSurcesCore() {
      return static_cast<FlatPrismElementForcesAndSurcesCore*>(ptrFE);
    }

  };

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

#endif //__FLATPRISMELEMENTFORCESANDSURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_prism_element Prism Element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
