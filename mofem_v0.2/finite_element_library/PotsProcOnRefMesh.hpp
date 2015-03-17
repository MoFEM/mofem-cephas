/** \file PostProcOnRefMesh.hpp 
 * \brief Postprocess fields on refined mesh made for 10 Node tets
 *
 * Create refined mesh, without enforcing continuity between element. Calculate
 * field values on nodes of that mesh.
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
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

#ifndef __POSTPROC_ON_REF_MESH_HPP
#define __POSTPROC_ON_REF_MESH_HPP

using namespace boost::numeric;

namespace MoFEM {

/** \brief Post processing 
  * \ingroup mofem_fs_post_proc  
  */
struct PostPocOnRefinedMesh: public TetElementForcesAndSourcesCore {

  moab::Core coreMesh;
  Interface &postProcMesh;

  bool tenNodesPostProcTets;
  int nbOfRefLevels;

  PostPocOnRefinedMesh(FieldInterface &m_field,bool ten_nodes_post_proc_tets = true,int nb_ref_levels = -1):
    TetElementForcesAndSourcesCore(m_field),postProcMesh(coreMesh),
    tenNodesPostProcTets(ten_nodes_post_proc_tets),nbOfRefLevels(nb_ref_levels) {

  }

  ~PostPocOnRefinedMesh() {
    ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
    if(pcomm_post_proc_mesh != NULL) {
      delete pcomm_post_proc_mesh;
    }
  }
  
  ublas::matrix<int> refTets;
  ublas::matrix<double> gaussPts_FirstOrder;
  vector<EntityHandle> mapGaussPts;

  // Gauss pts set on refined mesh
  int getRule(int order) { return -1; };

  struct CommonData {
    Range tEts;
    map<string,vector<ublas::vector<double> > > fieldMap;
    map<string,vector<ublas::matrix<double> > > gradMap;
  };
  CommonData commonData;

  PetscErrorCode generateRefereneElemenMesh() {
    PetscFunctionBegin;

    ErrorCode rval;
    PetscErrorCode ierr;

    int max_level = 0;
    if(nbOfRefLevels == -1) {
      PetscBool flg = PETSC_TRUE;
      PetscOptionsGetInt(PETSC_NULL,"-my_max_post_proc_ref_level",&max_level,&flg);
    } else {
      max_level = nbOfRefLevels;
    }

    double base_coords[] = {
      0,0,0,
      1,0,0,
      0,1,0,
      0,0,1 
    };

    moab::Core core_ref;
    Interface& moab_ref = core_ref;
      
    EntityHandle nodes[4];
    for(int nn = 0;nn<4;nn++) {
      rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
    }
    EntityHandle tet;
    rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERR_PETSC(rval);

    Core m_core_ref(moab_ref,PETSC_COMM_SELF,-2);
    FieldInterface& m_field_ref = m_core_ref;

    ierr = m_field_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

    for(int ll = 0;ll<max_level;ll++) {
      PetscPrintf(mField.get_comm(),"Refine Level %d\n",ll);
      Range edges;
      ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(ll),BitRefLevel().set(),MBEDGE,edges); CHKERRQ(ierr);
      Range tets;
      ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(ll),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);
      //refine mesh
      MeshRefinment& m_ref = m_core_ref;
      ierr = m_ref.add_verices_in_the_middel_of_edges(edges,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      ierr = m_ref.refine_TET(tets,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
    }
    
    Range elem_nodes;
    ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(max_level),BitRefLevel().set(),MBVERTEX,elem_nodes); CHKERRQ(ierr);

    map<EntityHandle,int> little_map;
    gaussPts_FirstOrder.resize(elem_nodes.size(),4,0);
    Range::iterator nit = elem_nodes.begin();
    for(int gg = 0;nit!=elem_nodes.end();nit++,gg++) {
      rval = moab_ref.get_coords(&*nit,1,&gaussPts_FirstOrder(gg,0)); CHKERR_PETSC(rval);
      little_map[*nit] = gg;
    }
    gaussPts_FirstOrder = trans(gaussPts_FirstOrder);

    Range tets;
    ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(max_level),BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr);

    refTets.resize(tets.size(),4);
    Range::iterator tit = tets.begin();
    for(int tt = 0;tit!=tets.end();tit++,tt++) {
      const EntityHandle *conn;
      int num_nodes;
      rval = moab_ref.get_connectivity(*tit,conn,num_nodes,false); CHKERR_PETSC(rval);
      for(int nn = 0;nn<num_nodes;nn++) {
	refTets(tt,nn) = little_map[conn[nn]];
      }
    }

    //moab_ref.list_entities(tets);

    PetscFunctionReturn(0);
  }

  PetscErrorCode setGaussPts(int order) {
    PetscFunctionBegin;

    try {

    PetscErrorCode ierr;
    ErrorCode rval;
 
    gaussPts_FirstOrder = trans(gaussPts_FirstOrder);
    mapGaussPts.resize(gaussPts_FirstOrder.size1());
    for(unsigned int gg = 0;gg<gaussPts_FirstOrder.size1();gg++) {
      rval = postProcMesh.create_vertex(&gaussPts_FirstOrder(gg,0),mapGaussPts[gg]); CHKERR_PETSC(rval);
    }
    gaussPts_FirstOrder = trans(gaussPts_FirstOrder);

    commonData.tEts.clear();
    for(unsigned int tt = 0;tt<refTets.size1();tt++) {
      EntityHandle conn[] = { 
	mapGaussPts[refTets(tt,0)], mapGaussPts[refTets(tt,1)],
	mapGaussPts[refTets(tt,2)], mapGaussPts[refTets(tt,3)] };
      EntityHandle tet;
      rval = postProcMesh.create_element(MBTET,conn,4,tet); CHKERR_PETSC(rval);
      commonData.tEts.insert(tet);
    }
 
    //cerr << commonData.tEts.size() << endl;

    EntityHandle meshset;
    rval = postProcMesh.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
    rval = postProcMesh.add_entities(meshset,commonData.tEts); CHKERR_PETSC(rval);
    //create higher order entities 
    if(tenNodesPostProcTets) {
      rval = postProcMesh.convert_entities(meshset,true,false,false); CHKERR_PETSC(rval);
    }

    commonData.tEts.clear();
    rval = postProcMesh.get_entities_by_type(meshset,MBTET,commonData.tEts,true); CHKERR_PETSC(rval);

    //cerr << "<-- " << commonData.tEts.size() << endl;
    Range nodes;
    rval = postProcMesh.get_connectivity(commonData.tEts,nodes,false); CHKERR_PETSC(rval);

    gaussPts.resize(nodes.size(),4);
    Range::iterator nit = nodes.begin();
    for(int gg = 0;nit!=nodes.end();nit++,gg++) {
      rval = postProcMesh.get_coords(&*nit,1,&gaussPts(gg,0)); CHKERR_PETSC(rval);
      gaussPts(gg,3) = 0;
    }
    gaussPts = trans(gaussPts);

    //cerr << gaussPts << endl;

    ublas::matrix<FieldData> N;
    N.resize(nodes.size(),4);
    ierr = ShapeMBTET(&*N.data().begin(),&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nodes.size()); CHKERRQ(ierr);
    //cerr << N << endl;

    ublas::matrix<double> coords_at_gauss_pts;
    coords_at_gauss_pts.resize(nodes.size(),3);

    EntityHandle fe_ent = fePtr->get_ent();

    ublas::vector<double> coords(12);
    {
      const EntityHandle *conn;
      int num_nodes;
      mField.get_moab().get_connectivity(fe_ent,conn,num_nodes,false);
      coords.resize(3*num_nodes);
      rval = mField.get_moab().get_coords(conn,num_nodes,&coords[0]); CHKERR_PETSC(rval);
      //cerr << coords << endl;
    }

    for(unsigned int gg = 0;gg<nodes.size();gg++) {
      for(int dd = 0;dd<3;dd++) {
	coords_at_gauss_pts(gg,dd) = cblas_ddot(4,&N(gg,0),1,&coords[dd],3);
      }
    }

    //cerr << coords_at_gauss_pts << endl;

    mapGaussPts.resize(nodes.size());
    nit = nodes.begin();
    for(int gg = 0;nit!=nodes.end();nit++,gg++) {
      rval = postProcMesh.set_coords(&*nit,1,&coords_at_gauss_pts(gg,0)); CHKERR_PETSC(rval);
      mapGaussPts[gg] = *nit;
    }

    //tEts.clear();
    //rval = postProcMesh.get_entities_by_type(0,MBTET,tEts,true); CHKERR_PETSC(rval);
    //cerr << "<--- <--- " << tEts.size() << endl;

    } catch (exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }


    PetscFunctionReturn(0);
  }


  struct OpHdivFunctions: public TetElementForcesAndSourcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;

    OpHdivFunctions(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name): 
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      postProcMesh(post_proc_mesh),mapGaussPts(map_gauss_pts) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      ErrorCode rval;
      PetscErrorCode ierr;

      vector<Tag> th;
      th.resize(data.getFieldData().size());

      double def_VAL[9] = { 0,0,0 };

      switch(type) {
	case MBTRI:
	  for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
	    ostringstream ss;
	    ss << "HDIV_FACE_" << side << "_" << dd;
	    rval = postProcMesh.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,th[dd],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);
	  }
	  break;
	case MBTET:
	  for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
	    ostringstream ss;
	    ss << "HDIV_TET_" << dd;
	    rval = postProcMesh.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,th[dd],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);
	  }
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }

      for(unsigned int gg = 0;gg<data.getHdivN().size1();gg++) {
	for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
	  ierr = postProcMesh.tag_set_data(th[dd],&mapGaussPts[gg],1,&data.getHdivN(gg)(dd,0)); CHKERRQ(ierr);
	}
      }
  
      PetscFunctionReturn(0);
    }

  };

  struct OpGetFieldValues: public TetElementForcesAndSourcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;
    CommonData &commonData;

    OpGetFieldValues(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name,CommonData &common_data): 
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      postProcMesh(post_proc_mesh),mapGaussPts(map_gauss_pts),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      ErrorCode rval;
      //PetscErrorCode ierr;
      
      const MoFEM::FEDofMoFEMEntity *dof_ptr = data.getFieldDofs()[0];
      string tag_name = dof_ptr->get_name();
      int rank = dof_ptr->get_max_rank();    

      int tag_length = rank;
      FieldSpace space = dof_ptr->get_space();
      switch(space) {
	case L2:
	case H1:
	  break;
	case HCURL:
	case HDIV:
	  tag_length *= 3;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }

      double def_VAL[tag_length];
      bzero(def_VAL,tag_length*sizeof(double));
      Tag th;
      rval = postProcMesh.tag_get_handle(tag_name.c_str(),tag_length,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

      // zero tags, this for Vertex if H1 and TRI if Hdiv, EDGE for Hcurl
      // no need for L2
      const void* tags_ptr[mapGaussPts.size()];
      int nb_gauss_pts = data.getN().size1();
      if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }

      switch(space) {
	case H1:
	  commonData.fieldMap[row_field_name].resize(nb_gauss_pts);
	  if(type == MBVERTEX) {
	    for(int gg = 0;gg<nb_gauss_pts;gg++) {
	      rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERR_PETSC(rval);
	      (commonData.fieldMap[row_field_name])[gg].resize(rank);
	      (commonData.fieldMap[row_field_name])[gg].clear();
	    }
	  }
	  rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERR_PETSC(rval);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    for(int rr = 0;rr<rank;rr++) {
	      ((double*)tags_ptr[gg])[rr] += cblas_ddot(
		(data.getFieldData().size()/rank),&(data.getN(gg)[0]),1,&(data.getFieldData()[rr]),rank);
	      (commonData.fieldMap[row_field_name])[gg][rr] += data.getFieldData()[rr];
	    }
	  }
	  break;  
	case L2:
	  rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERR_PETSC(rval);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    for(int rr = 0;rr<rank;rr++) {
	      ((double*)tags_ptr[gg])[rr] = cblas_ddot(
		(data.getFieldData().size()/rank),&(data.getN(gg)[0]),1,&(data.getFieldData()[rr]),rank);
	    }
	  }
	  break;
	case HDIV:
	  if(type == MBTRI && side == 0) {
	    for(int gg = 0;gg<nb_gauss_pts;gg++) {
	      rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERR_PETSC(rval);
	    }
	  }
	  rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERR_PETSC(rval);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    for(int rr = 0;rr<rank;rr++) {
	      for(int dd = 0;dd<3;dd++) {
		((double*)tags_ptr[gg])[3*rr+dd] += cblas_ddot(
		  (data.getFieldData().size()/rank),&(data.getHdivN(gg)(0,dd)),3,&(data.getFieldData()[rr]),rank);
		
	      }
	    }
	  }
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }

      PetscFunctionReturn(0);

    }

  };

  struct OpGetFieldGradientValues: public TetElementForcesAndSourcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;
    CommonData &commonData;

    OpGetFieldGradientValues(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name,CommonData &common_data): 
      TetElementForcesAndSourcesCore::UserDataOperator(field_name),
      postProcMesh(post_proc_mesh),mapGaussPts(map_gauss_pts),
      commonData(common_data) {}

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getFieldData().size()==0) PetscFunctionReturn(0);

      ErrorCode rval;
      //PetscErrorCode ierr;
       
      const MoFEM::FEDofMoFEMEntity *dof_ptr = data.getFieldDofs()[0];
      string tag_name = dof_ptr->get_name()+"_GRAD";
      int rank = dof_ptr->get_max_rank();    

      int tag_length = rank*3;
      FieldSpace space = dof_ptr->get_space();
      switch(space) {
	case L2:
	case H1:
	  break;
	case HCURL:
	case HDIV:
	  tag_length *= 3;
	  break;
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }

      double def_VAL[tag_length];
      bzero(def_VAL,tag_length*sizeof(double));
      Tag th;
      rval = postProcMesh.tag_get_handle(tag_name.c_str(),tag_length,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

      // zero tags, this for Vertex if H1 and TRI if Hdiv, EDGE for Hcurl
      // no need for L2
      const void* tags_ptr[mapGaussPts.size()];
      int nb_gauss_pts = data.getN().size1();
      if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }

      try {

      switch(space) {
	case H1:
	  commonData.gradMap[row_field_name].resize(nb_gauss_pts);
	  if(type == MBVERTEX) {
	    for(int gg = 0;gg<nb_gauss_pts;gg++) {
	      rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERR_PETSC(rval);
	      (commonData.gradMap[row_field_name])[gg].resize(rank,3);
	      (commonData.gradMap[row_field_name])[gg].clear();
	    }
	  }
	  rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERR_PETSC(rval);
	  for(int gg = 0;gg<nb_gauss_pts;gg++) {
	    for(int rr = 0;rr<rank;rr++) {
	      for(int dd = 0;dd<3;dd++) {
		for(unsigned int dof = 0;dof<(data.getFieldData().size()/rank);dof++) {
		  ((double*)tags_ptr[gg])[rank*rr+dd] += data.getDiffN(gg)(dof,dd)*data.getFieldData()[rank*dof+rr];
		  (commonData.gradMap[row_field_name])[gg](rr,dd) += data.getDiffN(gg)(dof,dd)*data.getFieldData()[rank*dof+rr];
		}
	      }
	    }
	  }
	  break;  
	default:
	  SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
      }

      } catch (exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
	SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }


      PetscFunctionReturn(0);

    }

  };

  PetscErrorCode addHdivFunctionsPostProc(const string field_name) {
    PetscFunctionBegin;
    get_op_to_do_Rhs().push_back(new OpHdivFunctions(postProcMesh,mapGaussPts,field_name));
    PetscFunctionReturn(0);
  }

  PetscErrorCode addFieldValuesPostProc(const string field_name) {
    PetscFunctionBegin;
    get_op_to_do_Rhs().push_back(new OpGetFieldValues(postProcMesh,mapGaussPts,field_name,commonData));
    PetscFunctionReturn(0);
  }

  PetscErrorCode addFieldValuesGradientPostProc(const string field_name) {
    PetscFunctionBegin;
    get_op_to_do_Rhs().push_back(new OpGetFieldGradientValues(postProcMesh,mapGaussPts,field_name,commonData));
    PetscFunctionReturn(0);
  }

  PetscErrorCode clearOperators() {
    PetscFunctionBegin;
    get_op_to_do_Rhs().clear();
    PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    ErrorCode rval;
    rval = postProcMesh.delete_mesh(); CHKERR_PETSC(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;

    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
    if(pcomm_post_proc_mesh == NULL) {
      pcomm_post_proc_mesh = new ParallelComm(&postProcMesh,mField.get_comm());
    }

    Range edges;
    rval = postProcMesh.get_entities_by_type(0,MBEDGE,edges,false);  CHKERR_PETSC(rval);
    rval = postProcMesh.delete_entities(edges); CHKERR_PETSC(rval);
    Range tris;
    rval = postProcMesh.get_entities_by_type(0,MBTRI,tris,false);  CHKERR_PETSC(rval);
    rval = postProcMesh.delete_entities(tris); CHKERR_PETSC(rval);

    Range tets;
    rval = postProcMesh.get_entities_by_type(0,MBTET,tets,false);  CHKERR_PETSC(rval);

    //cerr << "total tets size " << tets.size() << endl;

    int rank = pcomm->rank();
    Range::iterator tit = tets.begin();
    for(;tit!=tets.end();tit++) {
      rval = postProcMesh.tag_set_data(pcomm_post_proc_mesh->part_tag(),&*tit,1,&rank); CHKERR_PETSC(rval);
    }

    rval = pcomm->resolve_shared_ents(0); CHKERR_PETSC(rval);

    PetscFunctionReturn(0);
  }


};

}

#endif //__POSTPROC_ON_REF_MESH_HPP

/***************************************************************************//**
 * \defgroup mofem_fs_post_proc Post Process
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/


