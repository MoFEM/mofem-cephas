/** \file NitscheMethod.hpp
 * \ingroup nitsche_method
 * \brief Basic implementation of Nitsche method
 */

/*
 * This file is part of MoFEM.
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

#ifndef __NITCHE_BOUNDARY_CONDITIONS_HPP__
#define __NITCHE_BOUNDARY_CONDITIONS_HPP__

/** \brief Basic implementation of Nitsche's method
 * \ingroup nitsche_method

  This implementation have been used to enforce periodic constrains, see \ref
  nitsche_periodic. For theoretical founding of the method explanations see
  \cite nitsche_method_hal and \cite juntunen2009nitsche. For scratch book derivations
  see <a href="nitsche_bc_for_ch.pdf" target="_blank"><b>link</b></a>.

*/
struct NitscheMethod {

  struct MyFace: FaceElementForcesAndSourcesCore {
    int addToRule;
    MyFace(FieldInterface &m_field):
    FaceElementForcesAndSourcesCore(m_field),
    addToRule(0) {}
    int getRule(int order) { return order+addToRule; }
    /*int getRule(int order) { return -1; }
    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;
      int rule = order+addToRule+1;
      int nb_gauss_pts = triangle_ncc_order_num(rule);
      gaussPts.resize(3,nb_gauss_pts,false);
      triangle_ncc_rule(rule,nb_gauss_pts,&gaussPts(0,0),&gaussPts(2,0));
      PetscFunctionReturn(0);
    }*/

  };

  /** \brief Block data for Nitsche method
  * \ingroup nitsche_method
  */
  struct BlockData {
    double gamma;         ///< Penalty term, see \cite nitsche_method_hal
    double phi;           ///< Nitsche method parameter, see \cite nitsche_method_hal
    string faceElemName;  ///< name of element face
    Range fAces;          ///< faces on which constrain is applied
  };

  /** \brief Common data shared between finite element operators
  */
  struct CommonData {
    int nbActiveFaces;
    vector<EntityHandle> fAces;
    vector<const NumeredMoFEMFiniteElement *> facesFePtr;
    vector<VectorDouble> cOords;
    vector<MatrixDouble> faceNormals;
    vector<MatrixDouble> faceVertexShapeFunctions;
    vector<MatrixDouble> faceGaussPts;
    vector<MatrixDouble> coordsAtGaussPts;
    vector<MatrixDouble> hoCoordsAtGaussPts;
    vector<VectorDouble> rAy;

    /** \brief projection matrix

      First index is face index, then integration point and
      projection matrix 3x3

    */
    vector<vector<MatrixDouble > > P;

    /** \brief derivative of projection matrix in respect DoFs
     This is EntityType, face, gauss point at face.
     Matrix has 3 rows (components of displacements)
     Matrix has columns equal to number of DoFs on entity
    */
    map<EntityType,vector<vector<MatrixDouble> > > dP;

    CommonData() {
    }

  };

  /** \brief Get integration pts data on face
  */
  struct OpGetFaceData: FaceElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;

    OpGetFaceData(CommonData &common_data):
    FaceElementForcesAndSourcesCore::UserDataOperator("DISPLACEMENT",OPROW),
    commonData(common_data) {
    }

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      try {
        int faceInRespectToTet = getFEMethod()->nInTheLoop;
        if(type == MBVERTEX) {
          commonData.faceGaussPts.resize(4);
          commonData.faceGaussPts[faceInRespectToTet] = getGaussPts();
          commonData.faceVertexShapeFunctions.resize(4);
          commonData.faceVertexShapeFunctions[faceInRespectToTet] = data.getN();
          commonData.faceNormals.resize(4);
          commonData.faceNormals[faceInRespectToTet] = 0.5*getNormals_at_GaussPt();
          commonData.hoCoordsAtGaussPts.resize(4);
          commonData.hoCoordsAtGaussPts[faceInRespectToTet] = getHoCoordsAtGaussPts();
          commonData.cOords.resize(4);
          commonData.cOords[faceInRespectToTet] = getCoords();
          commonData.coordsAtGaussPts.resize(4);
          commonData.coordsAtGaussPts[faceInRespectToTet] = getCoordsAtGaussPts();
          commonData.rAy.resize(4);
          commonData.rAy[faceInRespectToTet] = -getNormal();
          commonData.rAy[faceInRespectToTet] /= norm_2(getNormal());
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  /// \brief Definition of volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {

    BlockData &blockData;
    CommonData &commonData;
    MyFace faceFE;
    int addToRule;

    virtual PetscErrorCode doAdditionalJobWhenGuassPtsAreCalulated() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    MyVolumeFE(
      FieldInterface &m_field,
      BlockData &block_data,
      CommonData &common_data
    ):
    VolumeElementForcesAndSourcesCore(m_field),
    blockData(block_data),
    commonData(common_data),
    faceFE(m_field),
    addToRule(0) {
      faceFE.getOpPtrVector().push_back(new OpGetFaceData(commonData));
    }

    int getRule(int order) { return -1; };

    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;

      gaussPts.resize(4,0,false);

      try {
        commonData.nbActiveFaces = 0;
        commonData.fAces.resize(4);
        EntityHandle tet = fePtr->get_ent();
        for(int ff = 0;ff<4;ff++) {
          EntityHandle face;
          rval = mField.get_moab().side_element(tet,2,ff,face); CHKERR_PETSC(rval);
          if(blockData.fAces.find(face)!=blockData.fAces.end()) {
            commonData.fAces[ff] = face;
            commonData.nbActiveFaces++;
          } else {
            commonData.fAces[ff] = 0;
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      try {
        commonData.facesFePtr.resize(4);
        for(int ff = 0;ff<4;ff++) {
          if(commonData.fAces[ff] != 0) {
            NumeredMoFEMFiniteElement_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator it,hi_it;
            it = problemPtr->numeredFiniteElements.get<Composite_Name_And_Ent_mi_tag>().
              lower_bound(boost::make_tuple(blockData.faceElemName,commonData.fAces[ff]));
            hi_it = problemPtr->numeredFiniteElements.get<Composite_Name_And_Ent_mi_tag>().
              upper_bound(boost::make_tuple(blockData.faceElemName,commonData.fAces[ff]));
            if(it == problemPtr->numeredFiniteElements.get<Composite_Name_And_Ent_mi_tag>().end()) {
              SETERRQ1(
                PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"No finite element found < %s >",
                blockData.faceElemName.c_str()
              );
            }
            commonData.facesFePtr[ff] = &*it;
          } else {
            commonData.facesFePtr[ff] = NULL;
          }
        }
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]!=NULL) {
            const NumeredMoFEMFiniteElement *faceFEPtr = commonData.facesFePtr[ff];
            faceFE.copy_basic_method(*this);
            faceFE.feName = blockData.faceElemName;
            faceFE.nInTheLoop = ff;
            faceFE.fePtr = faceFEPtr;
            faceFE.dataPtr = const_cast<FEDofMoFEMEntity_multiIndex*>(&faceFEPtr->fe_ptr->data_dofs);
            faceFE.rowPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&faceFEPtr->rows_dofs);
            faceFE.colPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&faceFEPtr->cols_dofs);
            faceFE.addToRule = addToRule;
            ierr = faceFE(); CHKERRQ(ierr);
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      try {
        int nb_gauss_pts = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          nb_gauss_pts += commonData.faceGaussPts[ff].size2();
        }
        gaussPts.resize(4,nb_gauss_pts,false);
        const double coords_tet[12] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };
        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          int nb_gauss_face_pts = commonData.faceGaussPts[ff].size2();
          for(int fgg = 0;fgg<nb_gauss_face_pts;fgg++,gg++) {
            for(int dd = 0;dd<3;dd++) {
              gaussPts(dd,gg) =
              commonData.faceVertexShapeFunctions[ff](fgg,0)*coords_tet[3*dataH1.facesNodes(ff,0)+dd] +
              commonData.faceVertexShapeFunctions[ff](fgg,1)*coords_tet[3*dataH1.facesNodes(ff,1)+dd] +
              commonData.faceVertexShapeFunctions[ff](fgg,2)*coords_tet[3*dataH1.facesNodes(ff,2)+dd];
            }
            gaussPts(3,gg) = commonData.faceGaussPts[ff](2,fgg);
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      ierr = doAdditionalJobWhenGuassPtsAreCalulated(); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  /** \brief Basic operated shared between all Nitsche operators
  */
  struct OpBasicCommon: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &nitscheBlockData;
    CommonData &nitscheCommonData;
    bool fieldDisp;

    OpBasicCommon(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      bool field_disp,
      const char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name,type
    ),
    nitscheBlockData(nitsche_block_data),
    nitscheCommonData(nitsche_common_data),
    fieldDisp(field_disp) {
      sYmm = false;
    }

    double faceRadius;
    virtual PetscErrorCode getFaceRadius(int ff) {
      PetscFunctionBegin;
      VectorDouble &coords = nitscheCommonData.cOords[ff];
      double center[3];
      tricircumcenter3d_tp(&coords[0],&coords[3],&coords[6],center,NULL,NULL);
      cblas_daxpy(3,-1,&coords[0],1,center,1);
      faceRadius = cblas_dnrm2(3,center,1);
      PetscFunctionReturn(0);
    }
    double gammaH;
    virtual PetscErrorCode getGammaH(double gamma,int gg) {
      PetscFunctionBegin;
      gammaH = gamma;
      //gammaH*= faceRadius;
      PetscFunctionReturn(0);
    }

  };

  /** \brief Calculate jacobian and variation of tractions
  */
  struct OpCommon: public OpBasicCommon {

    NonlinearElasticElement::BlockData &dAta;
    NonlinearElasticElement::CommonData &commonData;

    OpCommon(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      bool field_disp,
      const char type
    ):
    OpBasicCommon(
      field_name,nitsche_block_data,nitsche_common_data,field_disp,type
    ),
    dAta(data),
    commonData(common_data) {
    }

    VectorDouble dIsp;
    VectorDouble tRaction;
    MatrixDouble jAc_row;
    MatrixDouble jAc_col;
    MatrixDouble tRac_v;
    MatrixDouble tRac_u;

    PetscErrorCode getJac(
      DataForcesAndSurcesCore::EntData &data,int gg,MatrixDouble &jac
    ) {
      PetscFunctionBegin;
      try {
        int nb = data.getFieldData().size();
        jac.resize(9,nb,false);
        jac.clear();
        const MatrixAdaptor diffN = data.getDiffN(gg,nb/3);
        ublas::matrix<double> &jac_stress = commonData.jacStress[gg];
        for(int dd = 0;dd<nb/3;dd++) {
          for(int rr = 0;rr<3;rr++) {
            for(int ii = 0;ii<9;ii++) {
              for(int jj = 0;jj<3;jj++) {
                jac(ii,3*dd+rr) += jac_stress(ii,3*rr+jj)*diffN(dd,jj);
              }
            }
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode getTractionVariance(
      int gg,int fgg,int ff,MatrixDouble &jac,MatrixDouble &trac
    ) {
      PetscFunctionBegin;
      try {
        VectorAdaptor normal = VectorAdaptor(
          3,ublas::shallow_array_adaptor<double>(
            3,&nitscheCommonData.faceNormals[ff](fgg,0)
          )
        );
        trac.resize(3,jac.size2());
        trac.clear();
        for(unsigned int dd2 = 0;dd2<jac.size2();dd2++) {
          for(int nn = 0;nn<3;nn++) {
            trac(nn,dd2) = cblas_ddot(
              3,&jac(3*nn,dd2),jac.size2(),&normal[0],1
            );
          }
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      PetscFunctionReturn(0);
    }

  };

  /** \brief Calculate Nitsche method terms on left hand side
  * \ingroup nitsche_method
  */
  struct OpLhsNormal: public OpCommon {

    OpLhsNormal(
        const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      bool field_disp
    ):
    OpCommon(
      field_name,
      nitsche_block_data,
      nitsche_common_data,
      data,
      common_data,
      field_disp,
      UserDataOperator::OPROWCOL
    ) {
    }

    virtual PetscErrorCode calculateP(int gg,int fgg,int ff) {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    MatrixDouble kMatrix,kMatrix0,kMatrix1;

    PetscErrorCode doWork(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      if(col_data.getIndices().size()==0) PetscFunctionReturn(0);
      int nb_dofs_row = row_data.getIndices().size();
      int nb_dofs_col = col_data.getIndices().size();
      double gamma = nitscheBlockData.gamma;
      double phi = nitscheBlockData.phi;

      kMatrix0.resize(nb_dofs_row,nb_dofs_col,false);
      kMatrix1.resize(nb_dofs_row,nb_dofs_col,false);
      kMatrix.resize(nb_dofs_row,nb_dofs_col,false);
      kMatrix.clear();

      try {

        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(nitscheCommonData.facesFePtr[ff]==NULL) continue;
          int nb_face_gauss_pts = nitscheCommonData.faceGaussPts[ff].size2();
          ierr = getFaceRadius(ff); CHKERRQ(ierr);
          kMatrix0.clear();
          kMatrix1.clear();
          for(int fgg = 0;fgg<nb_face_gauss_pts;fgg++,gg++) {
            ierr = getGammaH(gamma,gg); CHKERRQ(ierr);
            double val = getGaussPts()(3,gg);
            ierr = getJac(row_data,gg,jAc_row); CHKERRQ(ierr);
            ierr = getTractionVariance(gg,fgg,ff,jAc_row,tRac_v); CHKERRQ(ierr);
            ierr = getJac(col_data,gg,jAc_col); CHKERRQ(ierr);
            ierr = getTractionVariance(gg,fgg,ff,jAc_col,tRac_u); CHKERRQ(ierr);
            VectorAdaptor normal = VectorAdaptor(
              3,ublas::shallow_array_adaptor<double>(
                3,&nitscheCommonData.faceNormals[ff](fgg,0)
              )
            );
            double area = cblas_dnrm2(3,&normal[0],1);

            ierr = calculateP(gg,fgg,ff); CHKERRQ(ierr);
            MatrixDouble &P = nitscheCommonData.P[ff][fgg];

            //P
            for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
              double n_row = row_data.getN()(gg,dd1);
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  for(int dd4 = 0;dd4<3;dd4++) {
                    if(!P(dd3,dd4)) continue;
                    kMatrix0(3*dd1+dd3,3*dd2+dd4) += val*n_row*P(dd3,dd4)*n_col*area;
                  }
                }
              }
            }
            for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
              double n_row = row_data.getN()(gg,dd1);
              for(int dd2 = 0;dd2<nb_dofs_col;dd2++) {
                for(int dd3 = 0;dd3<3;dd3++) {
                  double t = cblas_ddot(
                    3,&P(dd3,0),1,&tRac_u(0,dd2),tRac_u.size2()
                  );
                  kMatrix1(3*dd1+dd3,dd2) -= val*n_row*t;
                }
              }
            }
            for(int dd1 = 0;dd1<nb_dofs_row;dd1++) {
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  double t = cblas_ddot(
                    3,&P(0,dd3),3,&tRac_v(0,dd1),tRac_v.size2()
                  );
                  kMatrix1(dd1,3*dd2+dd3) -= phi*val*t*n_col;
                }
              }
            }
            //dP
            if(nitscheCommonData.dP[col_type].empty()!=0) {
              if(nitscheCommonData.dP.size()==4) {
                if(nitscheCommonData.dP[col_type][ff].size()==nb_face_gauss_pts) {
                  MatrixDouble &dP = nitscheCommonData.dP[col_type][ff][fgg];
                  if(dP.size1()==3 && dP.size2()==nb_dofs_col) {
                    VectorDouble &u = (commonData.dataAtGaussPts[commonData.spatialPositions])[gg];
                    for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
                      double n_row = row_data.getN()(gg,dd1);
                      for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                        for(int dd3 = 0;dd3<3;dd3++) {
                          for(int dd4 = 0;dd4<3;dd4++) {
                            kMatrix0(3*dd1+dd3,3*dd2+dd4) += val*n_row*u[dd4]*dP(dd4,3*dd2+dd4);
                          }
                        }
                      }
                    }
                    VectorDouble t = prod(trans(commonData.sTress[gg]),normal); // this is Piola-Kirchhoff I stress
                    for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
                      double n_row = row_data.getN()(gg,dd1);
                      for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                        for(int dd3 = 0;dd3<3;dd3++) {
                          for(int dd4 = 0;dd4<3;dd4++) {
                            kMatrix1(3*dd1+dd3,3*dd2+dd4) += val*n_row*t[dd4]*dP(dd4,3*dd2+dd4);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          kMatrix0 /= gammaH;
          noalias(kMatrix) += kMatrix0+kMatrix1;

        }


        if(gg != (int)getGaussPts().size2()) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"wrong number of gauss pts");
        }

        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_dofs_row,
          &row_data.getIndices()[0],
          nb_dofs_col,
          &col_data.getIndices()[0],
          &kMatrix(0,0),
          ADD_VALUES
        ); CHKERRQ(ierr);

      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

};

#endif // __NITCHE_BOUNDARY_CONDITIONS_HPP__

/***************************************************************************//**
 * \defgroup nitsche_method Nitsche Method
 * \ingroup user_modules
 ******************************************************************************/
