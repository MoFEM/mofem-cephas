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
 number of operator added pushing objects to OpPtrVector

 */
struct FaceElementForcesAndSourcesCore: public ForcesAndSurcesCore {

  MoABErrorCode rval;
  double aRea;;
  int num_nodes;
  const EntityHandle* conn;
  VectorDouble normal,tangent1,tangent2;
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
  OpSetContravariantPiolaTransoformOnTriangle opContravariantTransoform;
  OpSetCovariantPiolaTransoformOnTriangle opCovariantTransoform;

  FaceElementForcesAndSourcesCore(Interface &m_field):
    ForcesAndSurcesCore(m_field),
    dataH1(MBTRI),derivedDataH1(dataH1),
    dataHdiv(MBTRI),derivedDataHdiv(dataHdiv),
    dataHcurl(MBTRI),derivedDataHcurl(dataHcurl),
    dataL2(MBTRI),derivedDataL2(dataHdiv),
    dataNoField(MBTRI),dataNoFieldCol(MBTRI),
    meshPositionsFieldName("MESH_NODE_POSITIONS"),
    opHOCoordsAndNormals(
      hoCoordsAtGaussPts,nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt
    ),
    opContravariantTransoform(normal,nOrmals_at_GaussPt),
    opCovariantTransoform(
      normal,nOrmals_at_GaussPt,
      tangent1,tAngent1_at_GaussPt,
      tangent2,tAngent2_at_GaussPt
    ) {
    }

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

    /** \brief get triangle tangent 1
     */
    inline VectorDouble& getTangent1() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tangent1;
    }

    /** \brief get triangle tangent 2
     */
    inline VectorDouble& getTangent2() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tangent2;
    }

    /** \brief get normal as tensor
    */
    inline FTensor::Tensor1<double*,3> getTensor1Normal() {
      double *ptr = &*getNormal().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2]);
    }

    /** \brief get tangent1 as tensor
    */
    inline FTensor::Tensor1<double*,3> getTensor1Tangent1() {
      double *ptr = &*getTangent1().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2]);
    }

    /** \brief get tangent2 as tensor
    */
    inline FTensor::Tensor1<double*,3> getTensor2Tangent1() {
      double *ptr = &*getTangent2().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2]);
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

    /**
     * \brief get get coords at gauss points

     \code
     FTensor::Index<'i',3> i;
     FTensor::Tensor1<double,3> t_center;
     FTensor::Tensor1<double*,3> t_coords = getTensor1Coords();
     t_center(i) = 0;
     for(int nn = 0;nn!=3;nn++) {
        t_center(i) += t_coords(i);
        ++t_coords;
      }
      t_center(i) /= 3;
    \endcode

     */
    inline FTensor::Tensor1<double*,3> getTensor1Coords() {
      double *ptr = getCoords().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief get matrix of integration (Gauss) points on Face Element
    *  where columns 0,1 are x,y coordinates respectively and column 2 is a value of weight
    * for example getGaussPts()(1,13) returns y coordinate of 13th Gauss point on particular face element
    */
    inline MatrixDouble& getGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->gaussPts;
    }

    /** \brief get coordinates at Gauss pts.
     */
    inline MatrixDouble& getCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->coordsAtGaussPts;
    }

    /** \brief get coordinates at Gauss pts.
     */
    inline FTensor::Tensor1<double*,3> getTensor1CoordsAtGaussPts() {
      double *ptr = getCoordsAtGaussPts().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief coordinate at Gauss points (if hierarchical approximation of element geometry)

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

      */
    inline MatrixDouble& getHoCoordsAtGaussPts() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->hoCoordsAtGaussPts;
    }

    /** \brief get coordinates at Gauss pts (takes in account ho approx. of geometry)
     */
    inline FTensor::Tensor1<double*,3> getTensor1HoCoordsAtGaussPts() {
      if(
        getHoCoordsAtGaussPts().size1()==0 &&
        getHoCoordsAtGaussPts().size2()!=3
      ) {
        return getTensor1Coords();
      }
      double *ptr = getHoCoordsAtGaussPts().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief if higher order geometry return normals at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is available.

     */
    inline MatrixDouble& getNormalsAtGaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->nOrmals_at_GaussPt;
    }

    /** \brief if higher order geometry return normals at Gauss pts.
      *
      * \param gg gauss point number
      */
    inline ublas::matrix_row<MatrixDouble > getNormalsAtGaussPt(const int gg) {
      return ublas::matrix_row<MatrixDouble >(static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->nOrmals_at_GaussPt,gg);
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble& getTangent1AtGaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tAngent1_at_GaussPt;
    }

    /** \brief if higher order geometry return tangent vector to triangle at Gauss pts.

    Note: returned matrix has size 0 in rows and columns if no HO approximation
    of geometry is avaliable.

     */
    inline MatrixDouble& getTangent2AtGaussPt() {
      return static_cast<FaceElementForcesAndSourcesCore*>(ptrFE)->tAngent2_at_GaussPt;
    }

    /** \brief get normal at integration points

      Example:
      \code
      double nrm2;
      FTensor::Index<'i',3> i;
      FTensor::Tensor1<double*,3> t_normal = getTensor1NormalsAtGaussPt();
      for(int gg = gg!=data.getN().size1();gg++) {
        nrm2 = sqrt(t_normal(i)*t_normal(i));
        ++t_normal;
      }
      \endcode

    */
    inline FTensor::Tensor1<double*,3> getTensor1NormalsAtGaussPt() {
      double *ptr = &*getNormalsAtGaussPt().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief get tangent 1 at integration points

    */
    inline FTensor::Tensor1<double*,3> getTensor1Tangent1AtGaussPt() {
      double *ptr = &*getTangent1AtGaussPt().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
    }

    /** \brief get tangent 2 at integration points

    */
    inline FTensor::Tensor1<double*,3> getTensor1Tangent2AtGaussPt() {
      double *ptr = &*getTangent2AtGaussPt().data().begin();
      return FTensor::Tensor1<double*,3>(ptr,&ptr[1],&ptr[2],3);
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

    PetscErrorCode loopSideVolumes(
      const string &fe_name,VolumeElementForcesAndSourcesCoreOnSide &method
    );

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

  \ingroup mofem_forces_and_sources_tri_element

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

\ingroup mofem_forces_and_sources_tri_element

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
