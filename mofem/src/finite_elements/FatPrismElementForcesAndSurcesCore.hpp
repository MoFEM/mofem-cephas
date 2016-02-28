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

#ifndef __FATPRISMELEMENTFORCESANDSURCESCORE_HPP__
#define __FATPRISMELEMENTFORCESANDSURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief FatPrism finite element
 \ingroup mofem_forces_and_sources_prism_element

 User is implementing own operator at Gauss points level, by own object
 derived from FatPrismElementForcesAndSurcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 \todo Need to implement operators that will make this element work as Volume element

 */
struct FatPrismElementForcesAndSurcesCore: public VolumeElementForcesAndSourcesCore {

  MoABErrorCode rval;
  double aRea[2];
  VectorDouble normal;

  MatrixDouble gaussPtsTrianglesOnly;
  MatrixDouble coordsAtGaussPtsTrianglesOnly;
  MatrixDouble gaussPtsThroughThickness;

  DataForcesAndSurcesCore dataH1TrianglesOnly;
  DataForcesAndSurcesCore dataH1TroughThickness;

  MatrixDouble hoCoordsAtGaussPtsF3;
  MatrixDouble nOrmals_at_GaussPtF3;
  MatrixDouble tAngent1_at_GaussPtF3;
  MatrixDouble tAngent2_at_GaussPtF3;
  MatrixDouble hoCoordsAtGaussPtsF4;
  MatrixDouble nOrmals_at_GaussPtF4;
  MatrixDouble tAngent1_at_GaussPtF4;
  MatrixDouble tAngent2_at_GaussPtF4;
  OpGetCoordsAndNormalsOnPrism opHOCoordsAndNormals;

  FatPrismElementForcesAndSurcesCore(FieldInterface &m_field):
  VolumeElementForcesAndSourcesCore(m_field,MBPRISM),
  dataH1TrianglesOnly(MBPRISM),
  dataH1TroughThickness(MBPRISM),
  opHOCoordsAndNormals(
    hoCoordsAtGaussPtsF3,nOrmals_at_GaussPtF3,tAngent1_at_GaussPtF3,tAngent2_at_GaussPtF3,
    hoCoordsAtGaussPtsF4,nOrmals_at_GaussPtF4,tAngent1_at_GaussPtF4,tAngent2_at_GaussPtF4
  ) {
  }

  virtual int getRuleTrianglesOnly(int order) { return 2*order; };
  virtual int getRuleThroughThickness(int order) { return 2*order; };

  virtual PetscErrorCode setGaussPtsTrianglesOnly(int order_triangles_only)  {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode setGaussPtsThroughThickness(int order_thickness) {
    PetscFunctionBegin;
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"not implemented");
    PetscFunctionReturn(0);
  }

  /** \brief default operator for Flat Prism element
    * \ingroup mofem_forces_and_sources_prism_element
    */
  struct UserDataOperator: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    UserDataOperator(const string &field_name,const char type):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const string &row_field_name,const string &col_field_name,const char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(row_field_name,col_field_name,type) {
    }

    /** \brief get face aRea
    \param dd if dd == 0 it is for face F3 if dd == 1 is for face F4
    */
    inline double getArea(const int dd) { return ptrFE->aRea[0]; }

    inline double getAreaF3() { return ptrFE->aRea[0]; }
    inline double getAreaF4() { return ptrFE->aRea[1]; }

    /** \brief get triangle normal

    Normal has 6 elements, first 3 are for face F3 another three for face F4

     */
    inline VectorDouble& getNormal() { return ptrFE->normal; }

    inline VectorAdaptor getNormalF3() {
      double *data  = &(ptrFE->normal[0]);
      return VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,data));
    }

    inline VectorAdaptor getNormalF4() {
      double *data  = &(ptrFE->normal[3]);
      return VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,data));
    }

    /** \brief get Gauss pts. in the prism
     */
    inline MatrixDouble& getGaussPts() { return ptrFE->gaussPts; }

    /** \brief get Gauss pts. on triangles
    */
    inline MatrixDouble& getGaussPtsTrianglesOnly() { return ptrFE->gaussPtsTrianglesOnly; }

    /** \brief get Gauss pts. through thickness
    */
    inline MatrixDouble& getGaussPtsThroughThickness() { return ptrFE->gaussPtsThroughThickness; }

    /** \brief get coordinates at Gauss pts.

      Matrix has size (nb integration points)x(coordinates on F3 and F4 = 6), i.e. coordinates on face F3 and F4

     */
    inline MatrixDouble& getCoordsAtGaussPts() { return ptrFE->coordsAtGaussPts; }

    /** \brief get coordinates at Gauss pts.

      Matrix has size (nb integration points)x(coordinates on F3 and F4 = 6), i.e. coordinates on face F3 and F4

     */
    inline MatrixDouble& getCoordsAtGaussPtsTrianglesOnly() { return ptrFE->coordsAtGaussPtsTrianglesOnly; }


    /** \brief coordinate at Gauss points on face 3 (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPtsF3() { return ptrFE->hoCoordsAtGaussPtsF3; }

    /** \brief coordinate at Gauss points on face 4 (if hierarchical approximation of element geometry)
      */
    inline MatrixDouble& getHoCoordsAtGaussPtsF4() { return ptrFE->hoCoordsAtGaussPtsF4; }

    /** \brief if higher order geometry return normals at face F3 at Gauss pts.
     *
     * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
     *
     */
    inline MatrixDouble& getNormals_at_GaussPtF3() { return ptrFE->nOrmals_at_GaussPtF3; }

    /** \brief if higher order geometry return normals at face F4 at Gauss pts.
     *
     * Face 4 is top face in canonical triangle numeration, see \cite tautges2010canonical
     *
     */
    inline MatrixDouble& getNormals_at_GaussPtF4() { return ptrFE->nOrmals_at_GaussPtF4; }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF3(const int gg) {
      return ublas::matrix_row<MatrixDouble >(ptrFE->nOrmals_at_GaussPtF3,gg);
    }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * Face 3 is top face in canonical triangle numeration, see \cite tautges2010canonical
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPtF4(const int gg) {
      return ublas::matrix_row<MatrixDouble >(ptrFE->nOrmals_at_GaussPtF4,gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF3() { return ptrFE->tAngent1_at_GaussPtF3; }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF3() { return ptrFE->tAngent2_at_GaussPtF3; }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent1_at_GaussPtF4() { return ptrFE->tAngent1_at_GaussPtF4; }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.
     */
    inline MatrixDouble& getTangent2_at_GaussPtF4() { return ptrFE->tAngent2_at_GaussPtF4; }

    inline DataForcesAndSurcesCore& getTrianglesOnlyDataStructure() { return ptrFE->dataH1TrianglesOnly; }

    inline DataForcesAndSurcesCore& getTroughThicknessDataStructure() { return ptrFE->dataH1TroughThickness; }

    // /** \brief return pointer to triangle finite element object
    //  */
    // inline const FatPrismElementForcesAndSurcesCore* getFlatPrismElementForcesAndSurcesCore() {
    //   return ptrFE;
    // }

    /** \brief return pointer to fat prism finite element
     */
    inline const FatPrismElementForcesAndSurcesCore* getPrismFE() { return ptrFE; }

    PetscErrorCode setPtrFE(FatPrismElementForcesAndSurcesCore *ptr) {
      PetscFunctionBegin;
      VolumeElementForcesAndSourcesCore::UserDataOperator::setPtrFE(ptr);
      ptrFE = ptr;
      PetscFunctionReturn(0);
    }

    protected:
    FatPrismElementForcesAndSurcesCore *ptrFE;

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

#endif //__FATPRISMELEMENTFORCESANDSURCESCORE_HPP__
