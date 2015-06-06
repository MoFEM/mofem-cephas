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

  struct OpJacobianSmoother: public NonlinearElasticElement::OpJacobianPiolaKirchhoffStress {

    OpJacobianSmoother(
      const string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data,
      int tag,
      bool jacobian
    ):
    NonlinearElasticElement::OpJacobianPiolaKirchhoffStress(
      field_name,data,common_data,tag,jacobian,false,false
    ) {}

    PetscErrorCode calculateStress() {
      PetscFunctionBegin;

      try {

        PetscErrorCode ierr;
        ierr = dAta.materialAdoublePtr->calculateP_PiolaKirchhoffI(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);

        commonData.sTress[0].resize(3,3,false);
        for(int dd1 = 0;dd1<3;dd1++) {
          for(int dd2 = 0;dd2<3;dd2++) {
            dAta.materialAdoublePtr->SiGma(dd1,dd2) >>= (commonData.sTress[0])(dd1,dd2);
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

  struct OpRhsSmoother: public NonlinearElasticElement::OpRhsPiolaKirchhoff {

    bool sTabilised;
    Vec frontF;

    OpRhsSmoother(
      const string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data
    ):
    NonlinearElasticElement::OpRhsPiolaKirchhoff(field_name,data,common_data),
    sTabilised(false),
    frontF(PETSC_NULL)
    {}


    ublas::vector<int> frontIndices;

    PetscErrorCode aSemble(
      int row_side,EntityType row_type,DataForcesAndSurcesCore::EntData &row_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      int nb_dofs = row_data.getIndices().size();

      int *indices_ptr = &row_data.getIndices()[0];
      if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
        iNdices.resize(nb_dofs,false);
        noalias(iNdices) = row_data.getIndices();
        if(!sTabilised) {
          indices_ptr = &iNdices[0];
        }
        frontIndices.resize(nb_dofs,false);
        noalias(frontIndices) = row_data.getIndices();
        ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
        ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
        for(int ii = 0;dit!=dofs.end();dit++,ii++) {
          if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())!=dAta.forcesOnlyOnEntitiesRow.end()) {
            iNdices[ii] = -1;
          } else {
            frontIndices[ii] = -1;
          }
        }
      }

      ierr = VecSetValues(
        getFEMethod()->snes_f,
        nb_dofs,
        indices_ptr,&nf[0],
        ADD_VALUES
      ); CHKERRQ(ierr);

      if(frontF) {
        ierr = VecSetValues(
          frontF,
          nb_dofs,
          &frontIndices[0],&nf[0],
          ADD_VALUES
        ); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

  };

  struct OpLhsSmoother: public NonlinearElasticElement::OpLhsPiolaKirchhoff_dx {

    bool sTabilised;
    Vec tangentFrontF;

    OpLhsSmoother(
      const string vel_field,
      const string field_name,
      NonlinearElasticElement::BlockData &data,
      NonlinearElasticElement::CommonData &common_data
    ):
    NonlinearElasticElement::OpLhsPiolaKirchhoff_dx(vel_field,field_name,data,common_data),
    sTabilised(false),
    tangentFrontF(PETSC_NULL)
    {}

    ublas::vector<int> rowFrontIndices;
    ublas::vector<int> colFrontIndices;

    PetscErrorCode aSemble(
      int row_side,int col_side,
      EntityType row_type,EntityType col_type,
      DataForcesAndSurcesCore::EntData &row_data,
      DataForcesAndSurcesCore::EntData &col_data
    ) {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      int nb_row = row_data.getIndices().size();
      int nb_col = col_data.getIndices().size();

      int *row_indices_ptr = &row_data.getIndices()[0];
      int *col_indices_ptr = &col_data.getIndices()[0];

      if(!sTabilised) {

        if(!dAta.forcesOnlyOnEntitiesRow.empty()) {
          rowIndices.resize(nb_row,false);
          noalias(rowIndices) = row_data.getIndices();
          if(!sTabilised) {
            row_indices_ptr = &rowIndices[0];
          }
          rowFrontIndices.resize(nb_row,false);
          noalias(rowFrontIndices) = row_data.getIndices();
          ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
          ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
          for(int ii = 0;dit!=dofs.end();dit++,ii++) {
            if(dAta.forcesOnlyOnEntitiesRow.find((*dit)->get_ent())!=dAta.forcesOnlyOnEntitiesRow.end()) {
              rowIndices[ii] = -1;
            } else {
              colFrontIndices[ii] = -1;
            }
          }
        }

      }

      ierr = MatSetValues(
        getFEMethod()->snes_B,
        nb_row,row_indices_ptr,
        nb_col,col_indices_ptr,
        &k(0,0),ADD_VALUES
      ); CHKERRQ(ierr);

      //is symmetric
      if(row_side != col_side || row_type != col_type) {

        row_indices_ptr = &row_data.getIndices()[0];
        col_indices_ptr = &col_data.getIndices()[0];

        if(sTabilised) {

          if(!dAta.forcesOnlyOnEntitiesCol.empty()) {
            rowIndices.resize(nb_row,false);
            noalias(rowIndices) = row_data.getIndices();
            if(!sTabilised) {
              row_indices_ptr = &rowIndices[0];
            }
            rowFrontIndices.resize(nb_row,false);
            noalias(rowFrontIndices) = row_data.getIndices();
            ublas::vector<const FEDofMoFEMEntity*>& dofs = row_data.getFieldDofs();
            ublas::vector<const FEDofMoFEMEntity*>::iterator dit = dofs.begin();
            for(int ii = 0;dit!=dofs.end();dit++,ii++) {
              if(dAta.forcesOnlyOnEntitiesCol.find((*dit)->get_ent())!=dAta.forcesOnlyOnEntitiesCol.end()) {
                rowIndices[ii] = -1;
              } else {
                rowFrontIndices[ii] = -1;
              }
            }
          }

        }

        trans_k.resize(nb_col,nb_row,false);
        noalias(trans_k) = trans(k);
        ierr = MatSetValues(
          getFEMethod()->snes_B,
          nb_col,col_indices_ptr,
          nb_row,row_indices_ptr,
          &trans_k(0,0),ADD_VALUES
        ); CHKERRQ(ierr);

      }

      if(tangentFrontF && row_type == MBVERTEX) {

        double *f_tangent_front_mesh_array;
        if(tangentFrontF==PETSC_NULL) SETERRQ(PETSC_COMM_SELF,1,"vector for crack front not created");
        ierr = VecGetArray(tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);
        for(int nn = 0;nn<4;nn++) {
          FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;

          dit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
          lower_bound(boost::make_tuple("LAMBDA_CRACK_TANGENT_CONSTRAIN",getConn()[nn]));
          hi_dit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
          upper_bound(boost::make_tuple("LAMBDA_CRACK_TANGENT_CONSTRAIN",getConn()[nn]));

          if(distance(dit,hi_dit)>0) {

            FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator diit,hi_diit;

            diit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
            lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",getConn()[nn]));
            hi_diit = getFEMethod()->rowPtr->get<Composite_Name_And_Ent_mi_tag>().
            upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",getConn()[nn]));

            for(;diit!=hi_diit;diit++) {
              for(unsigned int ddd = 0;ddd<nb_col;ddd++) {
              if(diit->get_petsc_local_dof_idx()==-1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
              double g = f_tangent_front_mesh_array[diit->get_petsc_local_dof_idx()]*k(3*nn+diit->get_dof_rank(),ddd);
              int lambda_idx = dit->get_petsc_gloabl_dof_idx();
              ierr = MatSetValues(
                getFEMethod()->snes_B,1,&lambda_idx,1,&col_indices_ptr[ddd],&g,ADD_VALUES
              ); CHKERRQ(ierr);
            }
          }
        }

      }
      ierr = VecRestoreArray(tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);

    }

      PetscFunctionReturn(0);
    }

  };

};

#endif //__SMOOTHER_HPP__
