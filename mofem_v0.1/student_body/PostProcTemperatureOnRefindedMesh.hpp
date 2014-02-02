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

#ifndef __POSTPROCTERMOERATUREONREFINDEDMESH_HPP__
#define __POSTPROCTERMOERATUREONREFINDEDMESH_HPP__

#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

struct PostProcTemperatureOnRefMesh: public PostProcOnRefMesh_Base,FEMethod_UpLevelStudent {

    Tag th_themp,th_grad_themp;
    PostProcTemperatureOnRefMesh(Interface& _moab): 
      PostProcOnRefMesh_Base(),FEMethod_UpLevelStudent(_moab)
      {

      double def_VAL[3] = {0,0,0};
      rval = moab_post_proc.tag_get_handle("TEMPERATURE_VAL",1,
	MB_TYPE_DOUBLE,th_themp,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("TEMPERATURE_GARDIENT_VAL",3,
	MB_TYPE_DOUBLE,th_grad_themp,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    }

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Temperature PostProc on refined mesh\n");
      if(init_ref) PetscFunctionReturn(0);
      
      ierr = do_preprocess(); CHKERRQ(ierr);

      init_ref = true;

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = do_postproc(); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"End Tempereatutre PostProc on refined mesh\n");
      PetscFunctionReturn(0);
    }



    map<EntityHandle,EntityHandle> node_map;

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
      {

	H1L2_Data_at_Gauss_pt::iterator diit = h1l2_data_at_gauss_pt.find("TEMPERATURE");
	if(diit==h1l2_data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no field_name TEMPERATURE !!!");
	vector< ublas::vector<FieldData> > &data = diit->second;
	vector< ublas::vector<FieldData> >::iterator vit = data.begin();
	map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
	for(;vit!=data.end();vit++,mit++) {
	  rval = moab_post_proc.tag_set_data(th_themp,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
	}

      }

      {
  
	H1_DiffData_at_Gauss_pt::iterator  giit = h1_diff_data_at_gauss_pt.find("TEMPERATURE");
	if(giit == h1_diff_data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no field_name TEMPERATURE !!!");
	vector< ublas::matrix<FieldData> > &data = giit->second;
	vector< ublas::matrix<FieldData> >::iterator m_dit = data.begin();
	map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
	for(;m_dit!=data.end();m_dit++,mit++) {
	  rval = moab_post_proc.tag_set_data(th_grad_themp,&mit->second,1,&m_dit->data()[0]); CHKERR_PETSC(rval);
	}

      }

      PetscFunctionReturn(0);
    }

};


#endif //__POSTPROCTERMOERATUREONREFINDEDMESH_HPP__


