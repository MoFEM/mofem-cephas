/** \file FaceElementForcesAndSourcesCore.hpp

  \brief Implementation of face element.

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

#ifndef __FACEELEMENTFORCESANDSOURCESCORE_HPP__
#define __FACEELEMENTFORCESANDSOURCESCORE_HPP__

using namespace boost::numeric;

namespace MoFEM {

/** \brief Face finite element
 \ingroup mofem_forces_and_sources_tri_element

 User is implementing own operator at Gauss point level, by own object
 derived from FaceElementForcesAndSourcesCoreL::UserDataOperator.  Arbitrary
 number of operator added pushing objects to rowOpPtrVector and
 rowColOpPtrVector.

 */
struct FaceElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  MoABErrorCode rval;
  double aRea;;
  int num_nodes;
  const EntityHandle* conn;
  VectorDouble normal;
  VectorDouble coords;
  MatrixDouble gaussPts;
  MatrixDouble coordsAtGaussPts;

  DataForcesAndSurcesCore dataH1;
  DerivedDataForcesAndSurcesCore derivedDataH1;
  DataForcesAndSurcesCore dataHdiv;
  DerivedDataForcesAndSurcesCore derivedDataHdiv;
  DataForcesAndSurcesCore dataHcurl;
  DerivedDataForcesAndSurcesCore derivedDataHcurl;
  DataForcesAndSurcesCore dataL2;
  DerivedDataForcesAndSurcesCore derivedDataL2;
  DataForcesAndSurcesCore dataNoField,dataNoFieldCol;

  std::string meshPositionsFieldName;

  MatrixDouble hoCoordsAtGaussPts;
  MatrixDouble nOrmals_at_GaussPt;
  MatrixDouble tAngent1_at_GaussPt;
  MatrixDouble tAngent2_at_GaussPt;
  OpGetCoordsAndNormalsOnFace opHOCoordsAndNormals;
  OpSetPiolaTransoformOnTriangle opSetPiolaTransoformOnTriangle;

  FaceElementForcesAndSourcesCore(Interface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    dataHcurl(MBTRI),derivedDataHcurl(dataHdiv),
    dataL2(MBTRI),derivedDataL2(dataHdiv),
    dataNoField(MBTRI),dataNoFieldCol(MBTRI),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHOCoordsAndNormals(
      hoCoordsAtGaussPts,nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt
    ),
    opSetPiolaTransoformOnTriangle(normal,nOrmals_at_GaussPt) {};

  /** \brief default operator for TRI element
    * \ingroup mofem_forces_and_sources_tri_element
    */
  struct UserDataOperator: public ForcesAndSurcesCore::UserDataOperator {

    UserDataOperator(const std::string &field_name,const char type):
    ForcesAndSurcesCore::UserDataOperator(field_name,type) {}

    UserDataOperator(
      const std::string &row_field_name,const std::string &col_field_name,const char type):
      ForcesAndSurcesCore::UserDataOperator(row_field_name,col_field_name,type) {};

    inline double getArea() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->aRea;
    }

    /** \brief get triangle normal
     */
    inline VectorDouble& getNormal() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->normal;
    }

    /** \brief get element number of nodes
    */
    inline int getNumNodes() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->num_nodes;
    }

    /** \brief get element connectivity
     */
    inline const EntityHandle* getConn() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->conn;
    }

    /** \brief get triangle coordinates
     */
    inline VectorDouble& getCoords() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->coords;
    }

    /** \brief get triangle Gauss pts.
     */
    inline MatrixDouble& getGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->gaussPts;
    }

    /** \brief get coordinates at Gauss pts.
     */
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->coordsAtGaussPts;
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of element geometry)

    Note: Returenet matrix has size 0 in rows and columns if no HO approxmiation
    of geometry is avaliable.

      */
    inline MatrixDouble& getHoCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->hoCoordsAtGaussPts;
    }

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: Returenet matrix has size 0 in rows and columns if no HO approxmiation
    of geometry is avaliable.

     */
    inline MatrixDouble& getNormals_at_GaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->nOrmals_at_GaussPt;
    }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormals_at_GaussPt(const int gg) {
      return ublas::matrix_row<MatrixDouble >(static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->nOrmals_at_GaussPt,gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.

    Note: Returenet matrix has size 0 in rows and columns if no HO approxmiation
    of geometry is avaliable.

     */
    inline MatrixDouble& getTangent1_at_GaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tAngent1_at_GaussPt;
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.

    Note: Returenet matrix has size 0 in rows and columns if no HO approxmiation
    of geometry is avaliable.

     */
    inline MatrixDouble& getTangent2_at_GaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tAngent2_at_GaussPt;
    }

    /** \brief return pointer to triangle finite element object
     */
    inline const FaceElementForcesAndSourcesCore* getFaceElementForcesAndSourcesCore() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE);
    }

    /** \brief return pointer to Generic Triangle Finite Element object
     */
    inline const FaceElementForcesAndSourcesCore* getTriFE() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE);
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

/** \brief Calculate inverse of jacobian for face element

  It is assumed that face element is XY plane. Applied
  only for 2d problems.

  \todo Generalize function for arbitrary face orientation in 3d space

*/
struct OpCalculateInvJacForFace: public FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &invJac;
  OpCalculateInvJacForFace(const std::string &field_name,MatrixDouble &inv_jac):
  FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  invJac(inv_jac) {}
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );
};

/** \brief Transform local reference derivatives of shape functions to global derivatives

It is used for 2d problems.

*/
struct OpSetInvJacH1ForFace: public FaceElementForcesAndSourcesCore::UserDataOperator {
  MatrixDouble &invJac;
  OpSetInvJacH1ForFace(const std::string &field_name,MatrixDouble &inv_jac):
  FaceElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
  invJac(inv_jac) {}
  MatrixDouble diffNinvJac;
  PetscErrorCode doWork(
    int side,EntityType type,DataForcesAndSurcesCore::EntData &data
  );
};



}

#endif //__FACEELEMENTFORCESANDSOURCESCORE_HPP__

/***************************************************************************//**
 * \defgroup mofem_forces_and_sources_tri_element Face Element
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
