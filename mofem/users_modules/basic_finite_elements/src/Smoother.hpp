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

#ifndef __SMOOTHER_HPP__
#define __SMOOTHER_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

struct Smoother {

  struct SmootherBlockData {

    bool sTabilised;
    Vec frontF;
    Vec tangentFrontF;
    bool ownVectors;

    SmootherBlockData():
    sTabilised(false),
    frontF(PETSC_NULL),
    tangentFrontF(PETSC_NULL),
    ownVectors(false) {
    }

    ~SmootherBlockData() {
      if(ownVectors) {
        ierr = VecDestroy(&frontF); CHKERRABORT(PETSC_COMM_WORLD,ierr);
        ierr = VecDestroy(&tangentFrontF); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      }
    }

  };
  SmootherBlockData smootherData;

  std::map<int,NonlinearElasticElement::BlockData> setOfBlocks;
  NonlinearElasticElement::CommonData commonData;

  struct MyVolumeFE: public NonlinearElasticElement::MyVolumeFE {

    SmootherBlockData &smootherData;

    MyVolumeFE(MoFEM::Interface &m_field,SmootherBlockData &smoother_data):
    NonlinearElasticElement::MyVolumeFE(m_field),
    smootherData(smoother_data) {}

    PetscErrorCode preProcess() {
      PetscFunctionBegin;


      ierr = VolumeElementForcesAndSourcesCore::preProcess(); CHKERRQ(ierr);

      if(A != PETSC_NULL) {
        snes_B = A;
      }

      if(F != PETSC_NULL) {
        snes_f = F;
      }

      switch (ts_ctx) {
        case CTX_TSSETIFUNCTION: {
          if(!F) {
            snes_ctx = CTX_SNESSETFUNCTION;
            snes_f = ts_F;
          }
          break;
        }
        case CTX_TSSETIJACOBIAN: {
          if(!A) {
            snes_ctx = CTX_SNESSETJACOBIAN;
            snes_B = ts_B;
          }
          break;
        }
        default:
        break;
      }

      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
          if(smootherData.frontF) {
            ierr = VecZeroEntries(smootherData.frontF); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(smootherData.frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(smootherData.frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          }
        }
        break;
        default:
        break;
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;


      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
          if(smootherData.frontF) {
            ierr = VecAssemblyBegin(smootherData.frontF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(smootherData.frontF); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(smootherData.frontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(smootherData.frontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(smootherData.frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(smootherData.frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          }
          break;
          default:
          break;
        }
      }

      ierr = VolumeElementForcesAndSourcesCore::postProcess(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }


  };

  boost::shared_ptr<MyVolumeFE> feRhsPtr;
  boost::shared_ptr<MyVolumeFE> feLhsPtr;

  MyVolumeFE& feRhs; ///< calculate right hand side for tetrahedral elements
  MyVolumeFE& getLoopFeRhs() { return feRhs; } ///< get rhs volume element
  MyVolumeFE& feLhs; //< calculate left hand side for tetrahedral elements
  MyVolumeFE& getLoopFeLhs() { return feLhs; } ///< get lhs volume element

  Smoother(MoFEM::Interface &m_field):
  feRhsPtr(new MyVolumeFE(m_field,smootherData)),
  feLhsPtr(new MyVolumeFE(m_field,smootherData)),
  feRhs(*feRhsPtr),
  feLhs(*feLhsPtr) {
  }

  struct OpJacobianSmoother: public NonlinearElasticElement::OpJacobianPiolaKirchhoffStress {

    OpJacobianSmoother(
      const std::string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      int tag,
      bool jacobian
    ):
    NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
      field_name,data,common_data,tag,jacobian,false,false
    ) {
    }

    PetscErrorCode calculateStress(const int gg) {
      PetscFunctionBegin;

      try {


        ierr = dAta.materialAdoublePtr->calculateP_PiolaKirchhoffI(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);

        commonData.sTress[gg].resize(3,3,false);
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            dAta.materialAdoublePtr->P(dd1,dd2) >>= (commonData.sTress[gg])(dd1,dd2);
          }
        }

      } catch (const std::exception& ex) {
        std::ostringstream ss;
        ss << "throw in method: " << ex.what() << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpRhsSmoother: public NonlinearElasticElement::OpRhsPiolaKirchhoff {

    SmootherBlockData &smootherData;

    OpRhsSmoother(
      const std::string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      SmootherBlockData &smoother_data
    ):
    NonlinearElasticElement::OpRhsPiolaKirchhoff(field_name,data,common_data),
    smootherData(smoother_data)
    {}


    ublas::vector<int> frontIndices;

    PetscErrorCode aSemble(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;


      int nb_dofs = row_data.getIndices().size();

      int *indices_ptr = &row_data.getIndices()[0];
      if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
        iNdices.resize(nb_dofs,false);
        noalias(iNdices) = row_data.getIndices();
        if(!smootherData.sTabilised) {
          indices_ptr = &iNdices[0];
        }
        frontIndices.resize(nb_dofs,false);
        noalias(frontIndices) = row_data.getIndices();
        VectorDofs& dofs = row_data.getFieldDofs();
        VectorDofs::iterator dit = dofs.begin();
        for(int ii = 0;dit!=dofs.end();dit++,ii++) {
          if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->getEnt())!=dAta.forcesOnlyOnEntitiesRow.end()) {
            iNdices[ii] = -1;
          } else {
            frontIndices[ii] = -1;
          }
        }
        if(smootherData.frontF) {
          ierr = VecSetValues(
            smootherData.frontF,
            nb_dofs,
            &frontIndices[0],&nf[0],
            ADD_VALUES
          ); CHKERRQ(ierr);
        }
      }

      ierr = VecSetOption(
        getFEMethod()->snes_f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE
      );  CHKERRQ(ierr);

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        indices_ptr,&nf[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpLhsSmoother: public NonlinearElasticElement::OpLhsPiolaKirchhoff_dx {

    SmootherBlockData &smootherData;
    const std::string fieldCrackAreaTangentConstrains;

    OpLhsSmoother(
      const std::string vel_field,
      const std::string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      SmootherBlockData &smoother_data,
      const std::string crack_area_tangent_constrains// = "LAMBDA_CRACK_TANGENT_CONSTRAIN"
    ):
    NonlinearElasticElement::OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data),
    smootherData(smoother_data),
    fieldCrackAreaTangentConstrains(crack_area_tangent_constrains) {
    }

    ublas::vector<int> rowFrontIndices;

    PetscErrorCode aSemble(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;


      int nb_row = row_data.getIndices().size();
      int nb_col = col_data.getIndices().size();
      int *row_indices_ptr = &row_data.getIndices()[0];
      int *col_indices_ptr = &col_data.getIndices()[0];

      if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
        rowIndices.resize(nb_row,false);
        noalias(rowIndices) = row_data.getIndices();
        if(!smootherData.sTabilised) {
          row_indices_ptr = &rowIndices[0];
        }
        rowFrontIndices.resize(nb_row,false);
        noalias(rowFrontIndices) = row_data.getIndices();
        VectorDofs& dofs = row_data.getFieldDofs();
        VectorDofs::iterator dit = dofs.begin();
        for(int ii = 0;dit!=dofs.end();dit++,ii++) {
          if(dAta.forcesOnlyOnEntitiesRow.find(
            (*dit)->getEnt())!=dAta.forcesOnlyOnEntitiesRow.end()
          ) {
            rowIndices[ii] = -1;
          } else {
            rowFrontIndices[ii] = -1;
          }
        }
      }

      ierr = MatSetValues(
        getFEMethod()->snes_B,
        nb_row,row_indices_ptr,
        nb_col,col_indices_ptr,
        &k(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      if(smootherData.tangentFrontF) {

        // get tangent vector array
        double *f_tangent_front_mesh_array;
        ierr = VecGetArray(smootherData.tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);
        // iterate nodes on tet
        for(int nn = 0;nn<4;nn++) {

          // get indices with Lagrange multiplier at node nn
          FENumeredDofEntityByNameAndEnt::iterator dit,hi_dit;
          dit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
          lower_bound(boost::make_tuple(fieldCrackAreaTangentConstrains,getConn()[nn]));
          hi_dit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
          upper_bound(boost::make_tuple(fieldCrackAreaTangentConstrains,getConn()[nn]));

          // continue if Lagrange are on element
          if(distance(dit,hi_dit)>0) {

            FENumeredDofEntityByNameAndEnt::iterator diit,hi_diit;

            // get mesh node positions at node nn
            diit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
            lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",getConn()[nn]));
            hi_diit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
            upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",getConn()[nn]));

            // iterate over dofs on node nn
            for(;diit!=hi_diit;diit++) {
              // iterate overt dofs in element column
              for(int ddd = 0;ddd<nb_col;ddd++) {
                // check consistency, node has to be at crack front
                if(rowFrontIndices[3*nn+diit->get()->getDofCoeffIdx()]!=diit->get()->getPetscGlobalDofIdx()) {
                  SETERRQ2(
                    PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency %d != %d",
                    3*nn+diit->get()->getDofCoeffIdx(),diit->get()->getPetscGlobalDofIdx()
                  );
                }
                // dof is not on this partition
                if(diit->get()->getPetscLocalDofIdx()==-1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
                double g =
                f_tangent_front_mesh_array[diit->get()->getPetscLocalDofIdx()]*
                k(3*nn+diit->get()->getDofCoeffIdx(),ddd);
                int lambda_idx = dit->get()->getPetscGlobalDofIdx();
                ierr = MatSetValues(
                  getFEMethod()->snes_B,1,&lambda_idx,1,&col_indices_ptr[ddd],&g,ADD_VALUES
                ); CHKERRQ(ierr);
              }
            }
          }

        }
        ierr = VecRestoreArray(smootherData.tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);

      }

      PetscFunctionReturn(0);
    }

  };

};

#endif //__SMOOTHER_HPP__
