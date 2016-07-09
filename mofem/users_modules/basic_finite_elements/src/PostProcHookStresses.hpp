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


struct PostPorcStress: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

  Interface& mField;
  Interface &postProcMesh;
  std::vector<EntityHandle> &mapGaussPts;

  PostProcVolumeOnRefinedMesh::CommonData &commonData;

  PostPorcStress(
    MoFEM::Interface& m_field,
    Interface &post_proc_mesh,
    std::vector<EntityHandle> &map_gauss_pts,
    const std::string field_name,
    PostProcVolumeOnRefinedMesh::CommonData &common_data):
    MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
    mField(m_field),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts),
    commonData(common_data) {}

  PetscErrorCode getMatParameters(double *_lambda,double *_mu,int *_block_id) {
    PetscFunctionBegin;

    PetscErrorCode ierr;
    ErrorCode rval;

    *_lambda = 1;
    *_mu = 1;

    EntityHandle ent = getNumeredEntFiniteElementPtr()->getEnt();
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);

      Range meshsets;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERRQ_MOAB(rval);
      meshsets.insert(it->meshset);
      for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
        if( mField.get_moab().contains_entities(*mit,&ent,1) ) {
          *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
          *_mu = MU(mydata.data.Young,mydata.data.Poisson);
          *_block_id = it->get_msId();
          PetscFunctionReturn(0);
        }
      }
    }

    SETERRQ(PETSC_COMM_SELF,1,
      "Element is not in elastic block, however you run linear elastic analysis with that element\n"
      "top tip: check if you update block sets after mesh refinments or interface insertion");

    PetscFunctionReturn(0);
  }


  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;

    if(type != MBVERTEX) PetscFunctionReturn(0);
    if(data.getFieldData().size()==0) PetscFunctionReturn(0);

    ErrorCode rval;
    PetscErrorCode ierr;

    //const MoFEM::FEDofEntity *dof_ptr = data.getFieldDofs()[0];

    int id;
    double lambda,mu;
    ierr = getMatParameters(&lambda,&mu,&id); CHKERRQ(ierr);

    ublas::matrix<FieldData> D_lambda,D_mu,D;
    D_lambda.resize(6,6);
    D_lambda.clear();
    for(int rr = 0;rr<3;rr++) {
      for(int cc = 0;cc<3;cc++) {
        D_lambda(rr,cc) = 1;
      }
    }
    D_mu.resize(6,6);
    D_mu.clear();
    for(int rr = 0;rr<6;rr++) {
      D_mu(rr,rr) = rr<3 ? 2 : 1;
    }
    D = lambda*D_lambda + mu*D_mu;

    int tag_length = 9;
    double def_VAL[tag_length];
    bzero(def_VAL,tag_length*sizeof(double));
    Tag th_stress;
    rval = postProcMesh.tag_get_handle(
      "STRESS",9,MB_TYPE_DOUBLE,th_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);

    Tag th_prin_stress_vect1,th_prin_stress_vect2,th_prin_stress_vect3;
    rval = postProcMesh.tag_get_handle(
      "S1",3,MB_TYPE_DOUBLE,th_prin_stress_vect1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);
    rval = postProcMesh.tag_get_handle(
      "S2",3,MB_TYPE_DOUBLE,th_prin_stress_vect2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);
    rval = postProcMesh.tag_get_handle(
      "S3",3,MB_TYPE_DOUBLE,th_prin_stress_vect3,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);

    Tag th_prin_stress_vals;
    rval = postProcMesh.tag_get_handle(
      "PRINCIPAL_STRESS",3,MB_TYPE_DOUBLE,th_prin_stress_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);

    Tag th_id;
    int def_block_id = -1;
    rval = postProcMesh.tag_get_handle(
      "BLOCK_ID",1,MB_TYPE_INTEGER,th_id,MB_TAG_CREAT|MB_TAG_SPARSE,&def_block_id
    ); CHKERRQ_MOAB(rval);
    Range::iterator tit = commonData.tEts.begin();
    for(;tit!=commonData.tEts.end();tit++) {
      rval = postProcMesh.tag_set_data(th_id,&*tit,1,&id);  CHKERRQ_MOAB(rval);
    }

    ublas::vector<double> strain;
    ublas::vector<double> stress;
    ublas::matrix<double> Stress;

    //Combine eigenvalues and vectors to create principal stress vector
    ublas::vector<double> prin_stress_vect1(3);
    ublas::vector<double> prin_stress_vect2(3);
    ublas::vector<double> prin_stress_vect3(3);
    ublas::vector<double> prin_vals_vect(3);

    int nb_gauss_pts = data.getN().size1();
    if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    for(int gg = 0;gg<nb_gauss_pts;gg++) {

      strain.resize(6);
      strain[0] = (commonData.gradMap[rowFieldName][gg])(0,0);
      strain[1] = (commonData.gradMap[rowFieldName][gg])(1,1);
      strain[2] = (commonData.gradMap[rowFieldName][gg])(2,2);
      strain[3] = (commonData.gradMap[rowFieldName][gg])(0,1)+(commonData.gradMap[rowFieldName][gg])(1,0);
      strain[4] = (commonData.gradMap[rowFieldName][gg])(1,2)+(commonData.gradMap[rowFieldName][gg])(2,1);
      strain[5] = (commonData.gradMap[rowFieldName][gg])(0,2)+(commonData.gradMap[rowFieldName][gg])(2,0);

      stress.resize(6);
      noalias(stress) = prod(D,strain);

      Stress.resize(3,3);
      Stress(0,0) = stress[0];
      Stress(1,1) = stress[1];
      Stress(2,2) = stress[2];
      Stress(0,1) = Stress(1,0) = stress[3];
      Stress(1,2) = Stress(2,1) = stress[4];
      Stress(2,0) = Stress(0,2) = stress[5];

      rval = postProcMesh.tag_set_data(th_stress,&mapGaussPts[gg],1,&Stress(0,0)); CHKERRQ_MOAB(rval);

      ublas::matrix< FieldData > eigen_vectors = Stress;
      ublas::vector<double> eigen_values(3);

      //LAPACK - eigenvalues and vectors. Applied twice for initial creates memory space
      int n = 3, lda = 3, info, lwork = -1;
      double wkopt;
      info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),&wkopt,lwork);
      if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
      lwork = (int)wkopt;
      double work[lwork];
      info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),work,lwork);
      if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);

      //eigen_vectors = trans(eigen_vectors);
      for (int ii=0; ii < 3; ii++) {
        prin_vals_vect[0] = eigen_values[0];
        prin_vals_vect[1] = eigen_values[1];
        prin_vals_vect[2] = eigen_values[2];
        prin_stress_vect1[ii] = eigen_vectors.data()[ii+3*0];
        prin_stress_vect2[ii] = eigen_vectors.data()[ii+3*1];
        prin_stress_vect3[ii] = eigen_vectors.data()[ii+3*2];
      }

      //Tag principle stress vectors 1, 2, 3
      rval = postProcMesh.tag_set_data(th_prin_stress_vect1,&mapGaussPts[gg],1,&prin_stress_vect1[0]); CHKERRQ_MOAB(rval);
      rval = postProcMesh.tag_set_data(th_prin_stress_vect2,&mapGaussPts[gg],1,&prin_stress_vect2[0]); CHKERRQ_MOAB(rval);
      rval = postProcMesh.tag_set_data(th_prin_stress_vect3,&mapGaussPts[gg],1,&prin_stress_vect3[0]); CHKERRQ_MOAB(rval);
      rval = postProcMesh.tag_set_data(th_prin_stress_vals,&mapGaussPts[gg],1,&prin_vals_vect[0]); CHKERRQ_MOAB(rval);

    }

    PetscFunctionReturn(0);
  }

};
