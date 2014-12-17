/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __POSTPROCDISPLACEMENTANDSTRAINONREFINDEDMESH_HPP__
#define __POSTPROCDISPLACEMENTANDSTRAINONREFINDEDMESH_HPP__

namespace ObosleteUsersModules {

struct PostProcOnRefMesh_Base {

    //this is moab mesh of all refined elements
    Interface& moab_post_proc;
    moab::Core mb_instance_post_proc;

    //this is moab mesh for reference element
    Interface& moab_ref;
    moab::Core mb_instance_ref;

    int max_level;
    vector<EntityHandle> meshset_level;
    bool init_ref;
    bool do_broadcast;

    PostProcOnRefMesh_Base(): 
      moab_post_proc(mb_instance_post_proc),moab_ref(mb_instance_ref),
      max_level(0),init_ref(false),do_broadcast(true) {
      PetscBool flg = PETSC_TRUE;
      PetscOptionsGetInt(PETSC_NULL,"-my_max_post_proc_ref_level",&max_level,&flg);
      meshset_level.resize(max_level+1);
    }

    vector<double> g_NTET;
    PetscErrorCode do_preprocess() {
      PetscFunctionBegin;

      ErrorCode rval;
      PetscErrorCode ierr;

      double base_coords[] = {
	0,0,0,
	1,0,0,
	0,1,0,
	0,0,1 };
      
      //
      EntityHandle nodes[4];
      for(int nn = 0;nn<4;nn++) {
	rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
      }
      EntityHandle tet;
      rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERR_PETSC(rval);

      //
      MoFEM::Core core_ref(moab_ref,-1);
      FieldInterface& m_field_ref = core_ref;
      MeshRefinment& m_ref_ref = core_ref; 
      ierr = m_field_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

      for(int ll = 0;ll<max_level;ll++) {
	PetscPrintf(PETSC_COMM_WORLD,"Refine Level %d\n",ll);
	rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = m_field_ref.get_entities_by_ref_level(BitRefLevel().set(ll),BitRefLevel().set(),meshset_level[ll]); CHKERRQ(ierr);
	ierr = m_ref_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = m_ref_ref.refine_TET(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = m_field_ref.get_entities_by_ref_level(BitRefLevel().set(max_level),BitRefLevel().set(),meshset_level[max_level]); CHKERRQ(ierr);

      std::vector<double> ref_coords;
      rval = moab_ref.get_vertex_coordinates(ref_coords); CHKERR_PETSC(rval);
      g_NTET.resize(4*ref_coords.size()/3);
      ShapeMBTET(&g_NTET[0],&ref_coords[0],&ref_coords[ref_coords.size()/3],&ref_coords[2*ref_coords.size()/3],ref_coords.size()/3);


      PetscFunctionReturn(0);
    }

    PetscErrorCode do_postproc() {
      PetscFunctionBegin;
      ParallelComm* pcomm_post_proc = ParallelComm::get_pcomm(&moab_post_proc,MYPCOMM_INDEX);
      if(pcomm_post_proc == NULL) pcomm_post_proc =  new ParallelComm(&moab_post_proc,PETSC_COMM_WORLD);
      if(do_broadcast) {
	ErrorCode rval;
	Range tets;
	rval = moab_post_proc.get_entities_by_type(0,MBTET,tets); CHKERR_PETSC(rval);
	for(unsigned int rr = 0; rr<pcomm_post_proc->size();rr++) {
	  rval = pcomm_post_proc->broadcast_entities(rr,tets); CHKERR(rval);
	}
      }
      PetscFunctionReturn(0);
    }

};

struct PostProcDisplacementsOnRefMesh: public FEMethod_UpLevelStudent,PostProcOnRefMesh_Base {
    ParallelComm* pcomm;

    Interface& mOab;
    string field_name;

    PostProcDisplacementsOnRefMesh(Interface& _moab,string _field_name = "DISPLACEMENT"): 
      FEMethod_UpLevelStudent(_moab),PostProcOnRefMesh_Base(),mOab(_moab),field_name(_field_name) {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    }

    Tag th_disp,th_positions;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      if(init_ref) PetscFunctionReturn(0);
      
      ierr = do_preprocess(); CHKERRQ(ierr);

      double def_VAL[3] = {0,0,0};
      // create TAG
      string tag_name = field_name+"_VAL";
      rval = moab_post_proc.tag_get_handle(tag_name.c_str(),3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("MESH_NODAL_POSITIONS_VAL",3,MB_TYPE_DOUBLE,th_positions,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

      init_ref = true;

      PetscFunctionReturn(0);
    }


    map<EntityHandle,EntityHandle> node_map;

    PetscErrorCode do_operator() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Tag th_block_id;
      int def_marker = 0;
      rval = mOab.tag_get_handle("BLOCKID",1,MB_TYPE_INTEGER,th_block_id,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_PETSC(rval); 
      Tag th_block_id_pp;
      rval = moab_post_proc.tag_get_handle("BLOCKID",1,MB_TYPE_INTEGER,th_block_id_pp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_PETSC(rval); 

      EntityHandle fe_ent = fePtr->get_ent();
      int block_id;
      rval = mOab.tag_get_data(th_block_id,&fe_ent,1,&block_id); CHKERR_PETSC(rval);

      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      Range::iterator nit = ref_nodes.begin();
      node_map.clear();
      for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      Range ref_tets;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
      Range::iterator tit = ref_tets.begin();
      for(;tit!=ref_tets.end();tit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_tet;
	rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_block_id_pp,&ref_tet,1,&block_id); CHKERR_PETSC(rval);
      }

      //Get displacements at Gauss points
      H1L2_Data_at_Gauss_pt::iterator diit = h1l2_data_at_gauss_pt.find(field_name);
      if(diit==h1l2_data_at_gauss_pt.end()) SETERRQ1(PETSC_COMM_SELF,1,"no field_name %s !!!",field_name.c_str());
      vector< ublas::vector<FieldData> > &data = diit->second;
      vector< ublas::vector<FieldData> >::iterator vit = data.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;vit!=data.end();vit++,mit++) {
	rval = moab_post_proc.tag_set_data(th_disp,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
      }


      //Get mesh nodal positions at Gauss points
      diit = h1l2_data_at_gauss_pt.find("MESH_NODE_POSITIONS");
      if(diit!=h1l2_data_at_gauss_pt.end()) {
	  vector< ublas::vector<FieldData> > &data = diit->second;
	  vector< ublas::vector<FieldData> >::iterator vit = data.begin();
	  map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
	  for(;vit!=data.end();vit++,mit++) {
	    rval = moab_post_proc.tag_set_data(th_positions,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
	}
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = do_operator(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }


    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = do_postproc(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

};

struct PostProcDisplacemenysAndStarinOnRefMesh: public PostProcDisplacementsOnRefMesh {


    Tag th_strain;
    PostProcDisplacemenysAndStarinOnRefMesh(Interface& _moab,string _field_name): PostProcDisplacementsOnRefMesh(_moab,_field_name) {
      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0};
      rval = moab_post_proc.tag_get_handle("STRAIN_VAL",9,MB_TYPE_DOUBLE,th_strain,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);

      //Strains to Nodes in PostProc Mesh
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;viit!=GradU_at_GaussPt.end();viit++,mit++) {
	ublas::matrix< FieldData > GradU = *viit;
	ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
	rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
      }

      PetscFunctionReturn(0);
    }


};


struct PostProcFieldsAndGradientOnRefMesh: public PostProcDisplacementsOnRefMesh {

    Tag th_strain;
    PostProcFieldsAndGradientOnRefMesh(Interface& _moab): PostProcDisplacementsOnRefMesh(_moab,"SPATIAL_POSITION") {
      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0};
      rval = moab_post_proc.tag_get_handle("F_VAL",9,MB_TYPE_DOUBLE,th_strain,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);

      //Strains to Nodes in PostProc Mesh
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;viit!=GradU_at_GaussPt.end();viit++,mit++) {
	ublas::matrix< FieldData > F = *viit;
	rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(F.data()[0])); CHKERR_PETSC(rval);
      }

      PetscFunctionReturn(0);
    }

};

struct PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh: public PostProcDisplacemenysAndStarinOnRefMesh {

  FieldInterface& mField;
  double lambda,mu;

  ublas::matrix<double> D,D_lambda,D_mu;

  Tag th_stress,th_prin_stress_vect1,th_prin_stress_vect2,th_prin_stress_vect3,th_prin_stress_vals;
  bool propeties_from_BLOCKSET_MAT_ELASTICSET;
	
  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh(
    FieldInterface& _mField, string _field_name,double _lambda,double _mu): 
      PostProcDisplacemenysAndStarinOnRefMesh(_mField.get_moab(),_field_name),mField(_mField),lambda(_lambda),mu(_mu) {
    double def_VAL2[3] = { 0.0, 0.0, 0.0 };
    rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT1",3,MB_TYPE_DOUBLE,th_prin_stress_vect1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT2",3,MB_TYPE_DOUBLE,th_prin_stress_vect2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT3",3,MB_TYPE_DOUBLE,th_prin_stress_vect3,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VALUES",3,MB_TYPE_DOUBLE,th_prin_stress_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
    rval = moab_post_proc.tag_get_handle("STRESS_VAL",9,MB_TYPE_DOUBLE,th_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);


    // See FEAP - - A Finite Element Analysis Program
    D_lambda = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<3;rr++) {
      for(int cc = 0;cc<3;cc++) {
	D_lambda(rr,cc) = 1;
      }
    }
    D_mu = ublas::zero_matrix<FieldData>(6,6);
    for(int rr = 0;rr<6;rr++) {
      D_mu(rr,rr) = rr<3 ? 2 : 1;
    }
    //    D = lambda*D_lambda + mu*D_mu;
    
    propeties_from_BLOCKSET_MAT_ELASTICSET = false;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
      propeties_from_BLOCKSET_MAT_ELASTICSET = true;
    }


  }

  virtual PetscErrorCode calculateD(double _lambda,double _mu) {
    PetscFunctionBegin;
    
    D = _lambda*D_lambda + _mu*D_mu;
    //cerr << D_lambda << endl;
    //cerr << D_mu << endl;
    //cerr << D << endl;
    
    PetscFunctionReturn(0);
  }
  
  PetscErrorCode GetMatParameters(double *_lambda,double *_mu) {
    PetscFunctionBegin;
    
    *_lambda = lambda;
    *_mu = mu;
    
    
    if(propeties_from_BLOCKSET_MAT_ELASTICSET) {
      EntityHandle ent = fePtr->get_ent();
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
        
        Mat_Elastic mydata;
        ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
        
        Range meshsets;
        rval = mOab.get_entities_by_type(it->meshset,MBENTITYSET,meshsets,true); CHKERR_PETSC(rval);
        meshsets.insert(it->meshset);
        for(Range::iterator mit = meshsets.begin();mit != meshsets.end(); mit++) {
          if( mOab.contains_entities(*mit,&ent,1) ) {
            *_lambda = LAMBDA(mydata.data.Young,mydata.data.Poisson);
            *_mu = MU(mydata.data.Young,mydata.data.Poisson);
            PetscFunctionReturn(0);
          }
        }
        
    }
      
    SETERRQ(PETSC_COMM_SELF,1,"Element is not in elestic block, however you run linear elastic analysis with that element\n"
              "top tip: check if you update block sets after mesh refinments or interface insertion");
      
    }
    
    PetscFunctionReturn(0);
  }

  vector< ublas::matrix<FieldData> > invH;
  vector< FieldData > detH;

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    try {

    //Loop over elements
    ierr = do_operator(); CHKERRQ(ierr);
    
    //Calculated D Matrix
    double _lambda,_mu;
    ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
    ierr = calculateD(_lambda,_mu); CHKERRQ(ierr);

    //Higher order approximation of geometry
    ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

    //Strains to Nodes in PostProc Mesh: create vector containing matrices
    vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
    //du/dx
    ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
    vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
    map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();

    int gg = 0;
    for(;viit!=GradU_at_GaussPt.end();viit++,mit++,gg++) {


      ublas::matrix< FieldData > GradU = *viit;
      if(!invH.empty()) {
        //GradU = 
        //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
        //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
        //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
        //H = 
        //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
         //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
         //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]    
        //invH = 
        //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
        //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
        //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
        //GradU = 
        //[ dU/dX1 dU/dX2 dU/dX3 ]
        //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
        //[ dW/dX1 dW/dX2 dW/dX3 ] 
        GradU = prod( GradU, invH[gg] ); 
      }
      ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
      rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
      //caluate stress and save it into tag
      ublas::vector<FieldData> Strain_VectorNotation(6);
      Strain_VectorNotation[0] = Strain(0,0);
      Strain_VectorNotation[1] = Strain(1,1);
      Strain_VectorNotation[2] = Strain(2,2);
      Strain_VectorNotation[3] = 2*Strain(0,1);
      Strain_VectorNotation[4] = 2*Strain(1,2);
      Strain_VectorNotation[5] = 2*Strain(2,0);
      ublas::vector< FieldData > Stress_VectorNotation = prod( D, Strain_VectorNotation );
      ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);
      Stress(0,0) = Stress_VectorNotation[0];
      Stress(1,1) = Stress_VectorNotation[1];
      Stress(2,2) = Stress_VectorNotation[2];
      Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
      Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
      Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];
      
      rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
  
      //Principal stress vectors
      //Calculated from finding the eigenvalues and eigenvectors
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
      
      //Combine eigenvalues and vectors to create principal stress vector
      ublas::vector<double> prin_stress_vect1(3);
      ublas::vector<double> prin_stress_vect2(3);
      ublas::vector<double> prin_stress_vect3(3);
      ublas::vector<double> prin_vals_vect(3);
            
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
      rval = moab_post_proc.tag_set_data(th_prin_stress_vect1,&mit->second,1,&prin_stress_vect1[0]); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_prin_stress_vect2,&mit->second,1,&prin_stress_vect2[0]); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_prin_stress_vect3,&mit->second,1,&prin_stress_vect3[0]); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_prin_stress_vals,&mit->second,1,&prin_vals_vect[0]); CHKERR_PETSC(rval);
    
    }

    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }

};

struct PostProcScalarFieldsAndGradientOnRefMesh: public FEMethod_UpLevelStudent,PostProcOnRefMesh_Base {
  ParallelComm* pcomm;
  string field_name;
  PostProcScalarFieldsAndGradientOnRefMesh(Interface& _moab,string _field_name): 
      FEMethod_UpLevelStudent(_moab),PostProcOnRefMesh_Base(), field_name(_field_name) {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    }

  Tag th_scalar,th_grad_scalar,th_positions;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    if(init_ref) PetscFunctionReturn(0);
    ierr = do_preprocess(); CHKERRQ(ierr);
    double def_VAL[3] = {0,0,0};
    // create TAG
    string tag_name = field_name+"_VAL";
    string tag_name_gradient = field_name+"_GRADIENT_VAL";
    rval = moab_post_proc.tag_get_handle(tag_name.c_str(),1,MB_TYPE_DOUBLE,th_scalar,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle(tag_name_gradient.c_str(),3,MB_TYPE_DOUBLE,th_grad_scalar,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("MESH_NODAL_POSITIONS_VAL",3,MB_TYPE_DOUBLE,th_positions,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    init_ref = true;
    PetscFunctionReturn(0);
  }

  map<EntityHandle,EntityHandle> node_map;
  vector< ublas::matrix<FieldData> > invH;
  vector< FieldData > detH;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

    Range ref_nodes;
    rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
    if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    Range::iterator nit = ref_nodes.begin();
    node_map.clear();
    for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
      EntityHandle &node = node_map[*nit];
      rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
    }
    Range ref_tets;
    rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
    Range::iterator tit = ref_tets.begin();
    for(;tit!=ref_tets.end();tit++) {
      const EntityHandle *conn_ref;
      int num_nodes;
      rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
      EntityHandle conn_post_proc[num_nodes];
      for(int nn = 0;nn<num_nodes;nn++) {
	conn_post_proc[nn] = node_map[conn_ref[nn]];
      }
      EntityHandle ref_tet;
      rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
    }

    //Get scalar at Gauss points
    H1L2_Data_at_Gauss_pt::iterator diit = h1l2_data_at_gauss_pt.find(field_name);
    if(diit==h1l2_data_at_gauss_pt.end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"no field_name %s !!!",field_name.c_str());
    }
    vector< ublas::vector<FieldData> > &data = diit->second;
    vector< ublas::vector<FieldData> >::iterator vit = data.begin();
    map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
    for(;vit!=data.end();vit++,mit++) {
      rval = moab_post_proc.tag_set_data(th_scalar,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
    }

    //Higher order approximation of geometry
    ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

    //Gradient to Nodes in PostProc Mesh
    vector< ublas::matrix< FieldData > > Grad_at_GaussPt;
    ierr = GetGaussDiffDataVector(field_name,Grad_at_GaussPt); CHKERRQ(ierr);
    vector< ublas::matrix< FieldData > >::iterator viit = Grad_at_GaussPt.begin();
    mit = node_map.begin();
    int gg = 0;
    for(;viit!=Grad_at_GaussPt.end();viit++,mit++,gg++) {
      ublas::matrix< FieldData > grad = *viit;
      if(!invH.empty()) {
	grad = prod( trans( invH[gg] ), trans(grad) ); 
      }
      rval = moab_post_proc.tag_set_data(th_grad_scalar,&mit->second,1,&(grad.data()[0])); CHKERR_PETSC(rval);
    }

    //Get mesh nodal positions at Gauss points
    diit = h1l2_data_at_gauss_pt.find("MESH_NODE_POSITIONS");
    if(diit!=h1l2_data_at_gauss_pt.end()) {
      vector< ublas::vector<FieldData> > &data = diit->second;
      vector< ublas::vector<FieldData> >::iterator vit = data.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;vit!=data.end();vit++,mit++) {
	rval = moab_post_proc.tag_set_data(th_positions,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    ierr = do_postproc(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct PostProcL2VelocitiesFieldsAndGradientOnRefMesh: public PostProcDisplacementsOnRefMesh {

    PostProcL2VelocitiesFieldsAndGradientOnRefMesh(Interface& _moab): PostProcDisplacementsOnRefMesh(_moab,"VELOCITIES") {}

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      Range::iterator nit = ref_nodes.begin();
      node_map.clear();
      for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      Range ref_tets;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
      Range::iterator tit = ref_tets.begin();
      for(;tit!=ref_tets.end();tit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_tet;
	rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
      }

      //Get displacements at Gauss points
      H1L2_Data_at_Gauss_pt::iterator diit = h1l2_data_at_gauss_pt.find(field_name);
      if(diit==h1l2_data_at_gauss_pt.end()) SETERRQ1(PETSC_COMM_SELF,1,"no field_name %s !!!",field_name.c_str());
      vector< ublas::vector<FieldData> > &data = diit->second;
      vector< ublas::vector<FieldData> >::iterator vit = data.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;vit!=data.end();vit++,mit++) {
	rval = moab_post_proc.tag_set_data(th_disp,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
      }


      PetscFunctionReturn(0);
    }

};

}

#endif //__POSTPROCDISPLACEMENTANDSTRAINONREFINDEDMESH_HPP__

