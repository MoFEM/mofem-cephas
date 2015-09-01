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

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

struct NitscheMethod   {

  struct MyFace: FaceElementForcesAndSourcesCore {
    MyFace(FieldInterface m_field): FaceElementForcesAndSourcesCore() {}
     int getRule (int order) { return order}
  };

  struct BlockData {
    double gamma;
    double phi;
    string faceElemName;
    Range fAces;
  };

  struct CommonData {

    MyFace faceFE;
    CommonData(FieldInterface m_field):
    faceFE(m_field) {
    }

    int nbActiveFaces;
    vector<EntityHandle> fAces;
    vector<const NumeredMoFEMFiniteElement *> facesFePtr;
    vector<ublas::matrix<double> > faceNormals;
    vector<ublas::matrix<double> > faceVertexShapeFunctions;
    vector<ublas::matrix<double> > faceGuassPts;

  };

  struct OpGetFaceData: FaceElementForcesAndSourcesCore::UserDataOperator {

    CommonData &commonData;

    OpGetFaceData(commonData &common_data):
    FaceElementForcesAndSourcesCore::UserDataOperator("MESH_NODE_POSITIONS",OPROW),
    commonData(common_data) {
    }

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      int faceInRespectToTet = getFEMethod()->nInTheLoop;

      if(type == MBVERTEX) {
        commonData.faceNormals.resize(4);
        commonData.faceNormals[faceInRespectToTet] = getNormals_at_GaussPt();
        commonData.faceVertexShapeFunctions[faceInRespectToTet] = data.getN();
        commonData.faceGuassPts[faceInRespectToTet] = getGaussPts();
      }

      PetscFunctionReturn(0);
    }

  };

  /// \brief  definition of volume element
  struct MyVolumeFE: public VolumeElementForcesAndSourcesCore {

    BlockData &blockData;
    CommonData &commonData;

    int addToRule;
    MyVolumeFE(
      FieldInterface &m_field,
      BlockData &bloack_data,
      CommonData &common_data
    ):
    VolumeElementForcesAndSourcesCore(m_field),
    blockData(block_data),
    commonData(common_data)
    addToRule(1) {
      commonData.faceFE.getOpPtrVector().push_back(OpGetFaceData(commonData));
    }

    int getRule(int order) { return -1; };

    PetscErrorCode setGaussPts(int order) {
      PetscFunctionBegin;

      int commonData.nbActiveFaces = 0;
      commonData.fAces.resize(4);
      EntityHandle tet = fePtr->get_ent();
      for(int ff = 0;ff<4;ff++) {
        EntityHandle face;
        rval = mField->side_element(tet,2,ff,face); CHKERR_PETSC(rval);
        if(blockData.fAces.find(face)!=blockData.fAces.end()) {
          commonData.fAces[ff] = face;
          commonData.nbActiveFaces++;
        } else {
          commonData.fAces[ff] = 0;
        }
      }

      {
        commonData.facesFePtr.resize(4);
        for(int ff = 0;ff<4;ff++) {
          if(commonData.fAces[ff] != 0) {
            commonData.facesFePtr[ff] =
            &*problemPtr->numeredFiniteElements.get<Composite_Name_And_Ent_mi_tag>().
            find(boost::make_tuple(blockData.faceElemName,commonData.fAces[ff]));
          } else {
            commonData.facesFePtr[ff] = NULL;
          }
        }
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]!=NULL) {
            commonData.faceFE.copy_basic_method(*this);
            commonData.faceFE.feName = blockData.faceElemName;
            commonData.faceFE.nInTheLoop = ff;
            const NumeredMoFEMFiniteElement *faceFEPtr = commonData.facesFePtr[ff];
            commonData.faceFE.fePtr = faceFEPtr;
            commonData.faceFE.dataPtr = const_cast<FEDofMoFEMEntity_multiIndex*>(&faceFEPtr->fe_ptr->data_dofs);
            commonData.faceFE.rowPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&faceFEPtr->rows_dofs);
            commonData.faceFE.colPtr = const_cast<FENumeredDofMoFEMEntity_multiIndex*>(&faceFEPtr->cols_dofs);
            ierr = commonData.faceFE.operator()(); CHKERRQ(ierr);
          }
        }
      }

      {
        int nb_gauss_pts = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          nb_gauss_pts += commonData.faceGaussPts[ff](2,fgg);
        }
        gaussPts.resize(4,nb_gauss_pts,false);
        const double coords_tet[12] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };
        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          for(int fgg = 0;fgg<commonData.faceGaussPts[ff].size2();fgg++) {
            for(int dd = 0;dd<3;dd++) {
              gaussPts(dd,gg) =
              commonData.faceVertexShapeFunctions(fgg,0)*coords_tet[3*data.facesNodes(ff,0)+dd] +
              commonData.faceVertexShapeFunctions(fgg,1)*coords_tet[3*data.facesNodes(ff,1)+dd] +
              commonData.faceVertexShapeFunctions(fgg,2)*coords_tet[3*data.facesNodes(ff,2)+dd];
            }
            gaussPts(3,gg) = commonData.faceGaussPts[ff](2,fgg);
            gg++;
          }
        }
      }

      PetscFunctionReturn(0;)
    }

  };

  struct OpCommon: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    BlockData &nitscheBlockData;
    CommonData &nitscheCommonData;
    NonLinearElasticElement::BlockData &dAta;
    NonLinearElasticElement::CommonData &commonData;
    bool fieldDisp;

    OpCommon(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonLinearElasticElement::BlockData &data,
      NonLinearElasticElement::CommonData &common_data
      bool field_disp,
      const char type
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(
      field_name,type
    ),
    nitscheBlockData(nitsche_block_data),
    nitscheCommonData(nitsche_common_data),
    dAta(data),
    commonData(common_data),
    fieldDisp(field_disp) {
    }

    VcetorDouble dIsp;
    VectorDouble tRaction;
    MatrixDouble jAc_row;
    MatrixDouble jAc_col;
    MatrixDouble tRac_v;
    MatrixDouble tRac_u;

    PetscErrorCode getJac(
      DataForcesAndSurcesCore::EntData &data,int gg,MatrixDouble &jac
    ) {
      PetscFunctionBegin;
      int nb = data.getFieldData().size();
      jac.resize(9,nb,false);
      jac.clear();
      const MatrixAdaptor diffN = data.getDiffN(gg,nb/3);
      ublas::matrix<double> &jac_stress = commonData.jacStress[gg];
      for(int dd = 0;dd<nb/3;dd++) {
        for(int rr = 0;rr<3;rr++) {
          for(int ii = 0;ii<9;ii++) {
            for(int jj = 0;jj<3;jj++) {
              jAc(ii,3*dd+rr) += jac_stress(ii,3*rr+jj)*diffN(dd,jj);
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode getTractionVariance(
      DataForcesAndSurcesCore::EntData &data,int gg,int ff,
      MatrixDouble &jac,VectorDouble &trac
    ) {
      PetscFunctionBegin;
      VectorDouble &normal = nitscheCommonData.faceNormals[ff];
      trac.resize(3,jac.size2());
      trac.clear();
      for(unsigned int dd2 = 0;dd2<jAc.size2();dd2++) {
        for(unsigned int dd1 = 0;dd1<9;dd1++) {
          trac(0,dd2) = 0.5*jac(0,dd2)*normal[0]+jac(1,dd2)*normal[1]+jac(2,dd2)*normal[2];
          trac(1,dd2) = 0.5*jac(3,dd2)*normal[0]+jac(4,dd2)*normal[1]+jac(5,dd2)*normal[2];
          trac(2,dd2) = 0.5*jac(6,dd2)*normal[0]+jac(7,dd2)*normal[1]+jac(8,dd2)*normal[2];
        }
      }
      PetscFunctionReturn(0);
    }

  };

  struct OpRhsNormal: public OpCommon {

    OpRhsNormalPiolaKirchhoff(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonLinearElasticElement::BlockData &data,
      NonLinearElasticElement::CommonData &common_data
      bool field_disp
    ):
    OpCommon(
      field_name,
      nitsche_block_data,
      nitsche_common_data,
      data,
      common_data,
      bool field_disp,
      UserDataOperator::OPROW
    ) {
    }

    VectorDouble nF;

    PetscErrorCode doWork(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
        PetscFunctionReturn(0);
      }
      if(row_data.getIndices().size()==0) PetscFunctionReturn(0);
      int nb_dofs = row_data.getIndices().size();
      double gamma = nitscheBlockData.gamma;
      double phi = nitscheBlockData.phi;

      try {

        nF.resize(nb_dofs);
        nF.clear();

        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          const MatrixDouble& stress = commonData.sTress[gg];
          int nb_face_gauss_pts = nitscheCommonData.faceGuassPts[ff].size2();
          for(int fgg = 0;fgg<nb_face_gauss_pts;fgg++,gg++) {
            double val = getGaussPts()(3,gg);
            
            ierr = getJac(row_data,gg,jAc_row); CHKERRQ(ierr);
            ierr = getTractionVariance(row_data,gg,ff,jAc_row,tRac_v); CHKERRQ(ierr);
            VectorDouble &normal = nitscheCommonData.faceNormals[ff];
            double area = cblas_dnrm2(3,&normal[0],1)*0.5;
            tRaction.resize(3,false);
            noalias(tRaction) = prod(stress,normal);
            dIsp.resize(3,false);
            noalias(dIsp) = commonData.dataAtGaussPts[rowFieldName][gg];
            if(fieldDisp) {
              dIsp -= commonData.dataAtGaussPts["MESH_NODE_POSITIONS"][gg];
            }
            for(int dd = 0;dd<3;dd++) {
              nF += (1./gamma)*dIsp[dd]*data.getN(gg)*val*area;
              nF += tRaction[dd]*data.getN(gg)*val;
            }
            for(int dd = 0;dd<nb_dofs;dd++) {
              nF[dd] += val*phi*(dIsp[0]*tRac_v(0,dd)+dIsp[1]*tRac_v(1,dd)+dIsp[2]*tRac_v(2,dd));
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

  };

  struct OpLhsNormal: public OpCommon {

    OpLhsNormalPiolaKirchhoff(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonLinearElasticElement::BlockData &data,
      NonLinearElasticElement::CommonData &common_data
      bool field_disp
    ):
    OpCommon(
      field_name,
      nitsche_block_data,
      nitsche_common_data,
      data,
      common_data,
      bool field_disp,
      UserDataOperator::OPROWCOL
    ) {
    }

    MatrixDouble kMatrix;

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

      kMatrix.resize(nb_dofs_row,nb_dofs_col,false);
      kMatrix.clear();

      try {

        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          const MatrixDouble& stress = commonData.sTress[gg];
          int nb_face_gauss_pts = nitscheCommonData.faceGuassPts[ff].size2();
          for(int fgg = 0;fgg<nb_face_gauss_pts;fgg++,gg++) {
            double val = getGaussPts()(3,gg);
            ierr = getJac(row_data,gg,jAc_row); CHKERRQ(ierr);
            ierr = getTractionVariance(row_data,gg,ff,jAc_row,tRac_v); CHKERRQ(ierr);
            ierr = getJac(col_data,gg,jAc_col); CHKERRQ(ierr);
            ierr = getTractionVariance(col_data,gg,ff,jAc_col,tRac_u); CHKERRQ(ierr);
            VectorDouble &normal = nitscheCommonData.faceNormals[ff];
            double area = cblas_dnrm2(3,&normal[0],1)*0.5;

            for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
              double n_row = row_data.getN()(gg,dd1);
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  for(int dd4 = 0;dd4<3;dd4++) {
                    kMatrix(3*dd1+dd3,3*dd2+dd4) += (1./gamma)*val*n_row*n_col;
                  }
                }
              }
            }

            for(int dd1 = 0;dd1<nb_dofs_row;dd1++) {
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  kMatrix(dd1,3*dd2+dd3) += phi*val*tRac_v(dd1,dd3)*n_col;
                }
              }
            }

            for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
              double n_row = row_data.getN()(gg,dd1);
              for(int dd2 = 0;dd2<3;dd2++) {
                for(int dd3 = 0;dd3<nb_dofs_col;dd3++) {
                  kMatrix(3*dd1+dd2,dd3) += val*n_row*tRac_u(dd2,dd3);
                }
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

  };

  struct OpRhsGMatrixNormal: public OpCommon {

    OpRhsNormalPiolaKirchhoff(
      const string field_name,
      BlockData &nitsche_block_data,
      CommonData &nitsche_common_data,
      NonLinearElasticElement::BlockData &data,
      NonLinearElasticElement::CommonData &common_data
      bool field_disp
    ):
    OpCommon(
      field_name,
      nitsche_block_data,
      nitsche_common_data,
      data,
      common_data,
      bool field_disp,
      UserDataOperator::OPROWCOL
    ) {
    }

    MatrixDouble gMatrix;

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

      gMatrix.resize(nb_dofs_row,nb_dofs_col,false);
      gMatrix.clear();

      try {

        int gg = 0;
        for(int ff = 0;ff<4;ff++) {
          if(commonData.facesFePtr[ff]==NULL) continue;
          const MatrixDouble& stress = commonData.sTress[gg];
          int nb_face_gauss_pts = nitscheCommonData.faceGuassPts[ff].size2();
          for(int fgg = 0;fgg<nb_face_gauss_pts;fgg++,gg++) {
            double val = getGaussPts()(3,gg);
            ierr = getJac(row_data,gg,jAc_row); CHKERRQ(ierr);
            ierr = getTractionVariance(row_data,gg,ff,jAc_row,tRac_v); CHKERRQ(ierr);
            VectorDouble &normal = nitscheCommonData.faceNormals[ff];
            double area = cblas_dnrm2(3,&normal[0],1);

            for(int dd1 = 0;dd1<nb_dofs_row/3;dd1++) {
              double n_row = row_data.getN()(gg,dd1);
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  for(int dd4 = 0;dd4<3;dd4++) {
                    gMatrix(3*dd1+dd3,3*dd2+dd4) += (1./gamma)*val*n_row*n_col;
                  }
                }
              }
            }

            for(int dd1 = 0;dd1<nb_dofs_row;dd1++) {
              for(int dd2 = 0;dd2<nb_dofs_col/3;dd2++) {
                double n_col = col_data.getN()(gg,dd2);
                for(int dd3 = 0;dd3<3;dd3++) {
                  gMatrix(dd1,3*dd2+dd3) += phi*val*tRac_v(dd1,dd3)*n_col;
                }
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

  };

}

#endif // __NITCHE_BOUNDARY_CONDITIONS_HPP__
