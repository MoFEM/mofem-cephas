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

#ifndef __POSTPROCNONLINEARELASTICITYSTRESSEONREFINDEDMESH_HPP__
#define __POSTPROCNONLINEARELASTICITYSTRESSEONREFINDEDMESH_HPP__

namespace ObosleteUsersModules {

struct PostProcStressNonLinearElasticity: public PostProcDisplacementsOnRefMesh {

  FEMethod_ComplexForLazy &fe_method;

  Tag th_F,th_cauchy_stress,th_piola_stress,th_eshelby_stress,th_psi,th_j,th_themp,th_positions,th_fibreDirection1,th_fibreDirection2,th_stretch1,th_stretch2,th_fibreEnergy1,th_fibreEnergy2;
  PostProcStressNonLinearElasticity(Interface& _moab,FEMethod_ComplexForLazy &_fe_method):
    PostProcDisplacementsOnRefMesh(_moab,_fe_method.spatial_field_name),fe_method(_fe_method) {

    double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
    rval = moab_post_proc.tag_get_handle("F_VAL",9,MB_TYPE_DOUBLE,th_F,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("CAUCHY_STRESS_VAL",9,MB_TYPE_DOUBLE,th_cauchy_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PIOLA1_STRESS_VAL",9,MB_TYPE_DOUBLE,th_piola_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("ESHELBY_STRESS_VAL",9,MB_TYPE_DOUBLE,th_eshelby_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PSI_VAL",1,MB_TYPE_DOUBLE,th_psi,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("J_VAL",1,MB_TYPE_DOUBLE,th_j,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("TEMPERATURE_VAL",1,MB_TYPE_DOUBLE,th_themp,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    rval = moab_post_proc.tag_get_handle("MESH_NODAL_POSITIONS_VAL",3,MB_TYPE_DOUBLE,th_positions,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION_1",3,MB_TYPE_DOUBLE,th_fibreDirection1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("FIBRE_DIRECTION_2",3,MB_TYPE_DOUBLE,th_fibreDirection2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    rval = moab_post_proc.tag_get_handle("FIBRE_STRETCH_1",1,MB_TYPE_DOUBLE,th_stretch1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("FIBRE_STRETCH_2",1,MB_TYPE_DOUBLE,th_stretch2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    rval = moab_post_proc.tag_get_handle("FIBRE_ENERGY_1",1,MB_TYPE_DOUBLE,th_fibreEnergy1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("FIBRE_ENERGY_2",1,MB_TYPE_DOUBLE,th_fibreEnergy2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
  }

  vector<double> remeber_g_NTET;

  PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);

      ierr = fe_method.setProblem(problemPtr); CHKERRQ(ierr);
      ierr = fe_method.setFields(fieldsPtr); CHKERRQ(ierr);
      ierr = fe_method.setEnts(refinedEntitiesPtr,entitiesPtr); CHKERRQ(ierr);
      ierr = fe_method.setDofs(dofsPtr); CHKERRQ(ierr);
      ierr = fe_method.setFiniteElements(refinedFiniteElementsPtr,finiteElementsPtr); CHKERRQ(ierr);
      ierr = fe_method.setFiniteElementsEntities(finiteElementsEntitiesPtr); CHKERRQ(ierr);
      ierr = fe_method.setAdjacencies(adjacenciesPtr); CHKERRQ(ierr);

      ierr = fe_method.setFE(fePtr); CHKERRQ(ierr);
      ierr = fe_method.setData(dataPtr); CHKERRQ(ierr);
      ierr = fe_method.setRowData(rowPtr); CHKERRQ(ierr);
      ierr = fe_method.setColData(colPtr); CHKERRQ(ierr);

      remeber_g_NTET = fe_method.g_NTET;
      fe_method.g_NTET = g_NTET;

      ierr = fe_method.OpComplexForLazyStart(); CHKERRQ(ierr);
      ierr = fe_method.GetData(
	fe_method.order_x_edges,fe_method.order_x_faces,fe_method.order_x_volume,
	fe_method.dofs_x_edge_data,fe_method.dofs_x_edge,
	fe_method.dofs_x_face_data,fe_method.dofs_x_face,
	fe_method.dofs_x_volume,fe_method.dofs_x,
	fe_method.spatial_field_name); CHKERRQ(ierr);
      //for(int ee = 0;ee<6;ee++) {
	//cout << "ee " << ee << " spatial " << fe_method.dofs_x_edge_data[ee] << endl;
      //}
      ierr = fe_method.GetData(
	fe_method.order_X_edges,fe_method.order_X_faces,fe_method.order_X_volume,
	fe_method.dofs_X_edge_data,fe_method.dofs_X_edge,
	fe_method.dofs_X_face_data,fe_method.dofs_X_face,
	fe_method.dofs_X_volume,fe_method.dofs_X,
	fe_method.material_field_name); CHKERRQ(ierr);
      //for(int ee = 0;ee<6;ee++) {
	//cout << "ee " << ee << " material " << fe_method.dofs_X_edge_data[ee] << endl;
      //}

      int ee = 0;
      for(;ee<6;ee++) {
	if(fe_method.diffH1edgeNinvJac[ee].size()==0) {
	  fe_method.diff_edgeNinvJac[ee] = NULL;
	} else {
	 fe_method.diff_edgeNinvJac[ee] = &(fe_method.diffH1edgeNinvJac[ee])[0];
	}
	if(fe_method.H1edgeN[ee].size()==0) {
	  fe_method.edgeN[ee] = NULL;
	} else {
	  fe_method.edgeN[ee] = &(fe_method.H1edgeN[ee])[0];
	}
      }
      int ff = 0;
      if(fe_method.H1faceN.size() != 4) {
	SETERRQ1(PETSC_COMM_SELF,1,"size of should be 4 but is %u",fe_method.H1faceN.size());
      }
      for(;ff<4;ff++) {
	if(fe_method.diffH1faceNinvJac[ff].size() == 0) {
	  fe_method.diff_faceNinvJac[ff] = NULL;
	} else {
	  fe_method.diff_faceNinvJac[ff] = &(fe_method.diffH1faceNinvJac[ff])[0];
	}
	if(fe_method.H1faceN[ff].size()==0) {
	  fe_method.faceN[ff] = NULL;
	} else {
	  fe_method.faceN[ff] = &(fe_method.H1faceN[ff])[0];
	}
      }
      fe_method.diff_volumeNinvJac = &fe_method.diffH1elemNinvJac[0];
      fe_method.volumeN = &fe_method.H1elemN[0];

      ierr = fe_method.GetDofs_Termal_FromElementData(); CHKERRQ(ierr);

      double _lambda,_mu,_thermal_expansion;
      ierr = fe_method.GetMatParameters(&_lambda,&_mu,&_thermal_expansion,&fe_method.ptr_matctx); CHKERRQ(ierr);

      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      int gg =0;
      for(;mit!=node_map.end();mit++) {

	ublas::matrix< double > F(3,3);
	ublas::matrix< double > Piola1Stress(3,3);
	ublas::matrix< double > CauhyStress(3,3);
	ublas::matrix< double > EshelbyStress(3,3);
	double Psi,J,themp;
	int order_T_volume = 0;

	ierr = Calulate_Stresses_at_GaussPoint(
	      &fe_method.maxOrderEdgeH1[0],&fe_method.maxOrderFaceH1[0],fe_method.maxOrderElemH1,
	      &fe_method.order_X_edges[0],&fe_method.order_X_faces[0],fe_method.order_X_volume,
	      &fe_method.order_x_edges[0],&fe_method.order_x_faces[0],fe_method.order_x_volume,
	      fe_method.V,_lambda,_mu,fe_method.ptr_matctx,
	      &fe_method.diffNTETinvJac[0],&fe_method.diff_edgeNinvJac[0],&fe_method.diff_faceNinvJac[0],fe_method.diff_volumeNinvJac,
	      &fe_method.dofs_X.data()[0],&fe_method.dofs_X_edge[0],&fe_method.dofs_X_face[0],&*fe_method.dofs_X_volume.data().begin(),
	      &*fe_method.dofs_x.data().begin(),&fe_method.dofs_x_edge[0],&fe_method.dofs_x_face[0],&*fe_method.dofs_x_volume.data().begin(),
	      //temperature
	      _thermal_expansion,1,
	      &g_NTET[0],&fe_method.edgeN[0],&fe_method.faceN[0],fe_method.volumeN,
	      NULL,NULL,order_T_volume, &fe_method.dofs_temp.data()[0],NULL,NULL,NULL,
	      &*F.data().begin(),
	      &*Piola1Stress.data().begin(),
	      &*CauhyStress.data().begin(),
	      &*EshelbyStress.data().begin(),
	      &Psi,&J,&themp,gg); CHKERRQ(ierr);
	gg++;

	rval = moab_post_proc.tag_set_data(th_F,&mit->second,1,&(F.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_cauchy_stress,&mit->second,1,&(CauhyStress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_piola_stress,&mit->second,1,&(Piola1Stress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_eshelby_stress,&mit->second,1,&(EshelbyStress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_psi,&mit->second,1,&Psi); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_j,&mit->second,1,&J); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_themp,&mit->second,1,&themp); CHKERR_PETSC(rval);

	//Get mesh nodal positions at Gauss points
	H1L2_Data_at_Gauss_pt::iterator diit = h1l2_data_at_gauss_pt.find(fe_method.material_field_name.c_str());
	if(diit!=h1l2_data_at_gauss_pt.end()) {
	  vector< ublas::vector<FieldData> > &data = diit->second;
	  vector< ublas::vector<FieldData> >::iterator vit = data.begin();
	  map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
	  for(;vit!=data.end();vit++,mit++) {
	    rval = moab_post_proc.tag_set_data(th_positions,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
	  }
	}

	if(get_PhysicalEquationNumber()==eberleinholzapfel1) {
      ctx_EberleinHolzapfel1 *material_data = (ctx_EberleinHolzapfel1*)fe_method.ptr_matctx;

      double k1;
      double k2;
      k1 = material_data->k1;
      k2 = material_data->k2;

      //Fibre 1
      ublas::vector<double> A1;
      ublas::vector<double> a1;
      A1.resize(3);
      a1.resize(3);
      cblas_dcopy(3,&material_data->fibre_vector_a1[0],1,&*A1.data().begin(),1);
      if (norm_2(A1)>0){
      //a1 current configuration
      a1=prod(F,A1);
      rval = moab_post_proc.tag_set_data(th_fibreDirection1,&mit->second,1,&*a1.data().begin()); CHKERR_PETSC(rval);

      //Fibre Stretch
      double stretch1 = norm_2(a1)/norm_2(A1);
      rval = moab_post_proc.tag_set_data(th_stretch1,&mit->second,1,&stretch1); CHKERR_PETSC(rval);

      //Fibre Energy
      double psi_f1;
      psi_f1 = (k1/(2*k2))*(exp(k2*pow(pow(stretch1,2)-1,2))-1);
      rval = moab_post_proc.tag_set_data(th_fibreEnergy1,&mit->second,1,&psi_f1); CHKERR_PETSC(rval);
      }

      //Fibre 2
      ublas::vector<double> A2;
      ublas::vector<double> a2;
      A2.resize(3);
      a2.resize(3);
      cblas_dcopy(3,&material_data->fibre_vector_a2[0],1,&*A2.data().begin(),1);
      if (norm_2(A2)>0){
      //a2 current configuration
      a2=prod(F,A2);
      rval = moab_post_proc.tag_set_data(th_fibreDirection2,&mit->second,1,&*a2.data().begin()); CHKERR_PETSC(rval);

      //Fibre Stretch
      double stretch2 = norm_2(a2)/norm_2(A2);
      rval = moab_post_proc.tag_set_data(th_stretch2,&mit->second,1,&stretch2); CHKERR_PETSC(rval);

      //Fibre Energy
      double psi_f2;
      psi_f2 = (k1/(2*k2))*(exp(k2*pow(pow(stretch2,2)-1,2))-1);
      rval = moab_post_proc.tag_set_data(th_fibreEnergy2,&mit->second,1,&psi_f2); CHKERR_PETSC(rval);
      }
  }
      }

      fe_method.g_NTET = remeber_g_NTET;

      PetscFunctionReturn(0);
  }
};

}

#endif //__POSTPROCNONLINEARELASTICITYSTRESSEONREFINDEDMESH_HPP__

