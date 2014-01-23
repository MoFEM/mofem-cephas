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


#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "FEMethod_ComplexForLazy.hpp"
#include "complex_for_lazy.h"

using namespace MoFEM;

struct PostProcStressNonLinearElasticity: public PostProcDisplacementsOnRefMesh {

  FEMethod_ComplexForLazy &fe_method;

  Tag th_cauchy_stress,th_piola_stress,th_eshelby_stress,th_psi,th_j;
  PostProcStressNonLinearElasticity(Interface& _moab,FEMethod_ComplexForLazy &_fe_method): 
    PostProcDisplacementsOnRefMesh(_moab,"SPATIAL_POSITION"),fe_method(_fe_method) {

    double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
    rval = moab_post_proc.tag_get_handle("CAUCHY_STRESS_VAL",9,MB_TYPE_DOUBLE,th_cauchy_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PIOLA1_STRESS_VAL",9,MB_TYPE_DOUBLE,th_piola_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("ESHELBY_STRESS_VAL",9,MB_TYPE_DOUBLE,th_eshelby_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("PSI_VAL",1,MB_TYPE_DOUBLE,th_psi,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("J_VAL",1,MB_TYPE_DOUBLE,th_j,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

  }

  vector<double> remeber_g_NTET;

  PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);

      ierr = fe_method.set_problem(problem_ptr); CHKERRQ(ierr);
      ierr = fe_method.set_fields(moabfields); CHKERRQ(ierr);
      ierr = fe_method.set_ents_multiIndex(ents_moabfield); CHKERRQ(ierr);
      ierr = fe_method.set_dofs_multiIndex(dofs_moabfield); CHKERRQ(ierr);
      ierr = fe_method.set_fes_multiIndex(finite_elements); CHKERRQ(ierr);
      ierr = fe_method.set_fes_data_multiIndex(finite_elements_moabents); CHKERRQ(ierr);
      ierr = fe_method.set_adjacencies(fem_adjacencies); CHKERRQ(ierr);

      ierr = fe_method.set_fe(fe_ptr); CHKERRQ(ierr);
      ierr = fe_method.set_data_multIndex(data_multiIndex); CHKERRQ(ierr);
      ierr = fe_method.set_row_multIndex(row_multiIndex); CHKERRQ(ierr);
      ierr = fe_method.set_col_multIndex(col_multiIndex); CHKERRQ(ierr);

      remeber_g_NTET = fe_method.g_NTET;
      fe_method.g_NTET = g_NTET;

      ierr = fe_method.OpComplexForLazyStart(); CHKERRQ(ierr);
      ierr = fe_method.GetData(fe_method.dofs_x_edge_data,fe_method.dofs_x_edge,
	fe_method.dofs_x_face_data,fe_method.dofs_x_face,
	fe_method.dofs_x_volume,fe_method.dofs_x,
	fe_method.spatial_field_name); CHKERRQ(ierr);

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

      ierr = fe_method.GetDofs_X_FromElementData(); CHKERRQ(ierr);
      ierr = fe_method.GetDofs_Termal_FromElementData(); CHKERRQ(ierr);

      double _lambda,_mu,_termal_expansion;
      ierr = fe_method.GetMatParameters(&_lambda,&_mu,&_termal_expansion,fe_method.ptr_matctx); CHKERRQ(ierr);

      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      int gg =0;
      for(;mit!=node_map.end();mit++) {

	ublas::matrix< double > Piola1Stress(3,3);
	ublas::matrix< double > CauhyStress(3,3);
	ublas::matrix< double > EshelbyStress(3,3);
	double Psi,J;

	ierr = Calulate_Stresses_at_GaussPoint(
	      &fe_method.order_edges[0],&fe_method.order_faces[0],fe_method.order_volume,fe_method.V,_lambda,_mu,fe_method.ptr_matctx, 
	      &fe_method.diffNTETinvJac[0],&fe_method.diff_edgeNinvJac[0],&fe_method.diff_faceNinvJac[0],fe_method.diff_volumeNinvJac, 
	      &fe_method.dofs_X.data()[0],&*fe_method.dofs_x.data().begin(),
	      &fe_method.dofs_x_edge[0],&fe_method.dofs_x_face[0],&*fe_method.dofs_x_volume.data().begin(), 
	      //temperature
	      _termal_expansion,
	      &g_NTET[0],&fe_method.edgeN[0],&fe_method.faceN[0],fe_method.volumeN,
	      NULL,NULL,NULL, &fe_method.dofs_temp.data()[0],NULL,NULL,NULL,
	      &*Piola1Stress.data().begin(),&*CauhyStress.data().begin(),&*EshelbyStress.data().begin(),&Psi,&J,gg); CHKERRQ(ierr);
	gg++;

	rval = moab_post_proc.tag_set_data(th_cauchy_stress,&mit->second,1,&(CauhyStress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_piola_stress,&mit->second,1,&(Piola1Stress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_eshelby_stress,&mit->second,1,&(EshelbyStress.data()[0])); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_psi,&mit->second,1,&Psi); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_j,&mit->second,1,&J); CHKERR_PETSC(rval);


      }

      fe_method.g_NTET = remeber_g_NTET;

      PetscFunctionReturn(0);
  }

};

#endif //__POSTPROCNONLINEARELASTICITYSTRESSEONREFINDEDMESH_HPP__

