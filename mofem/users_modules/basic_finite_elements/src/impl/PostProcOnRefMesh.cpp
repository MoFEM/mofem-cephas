/** \file PostProcOnRefMesh.cpp
 * \brief Postprocess fields on refined mesh made for 10 Node tets
 *
 * Create refined mesh, without enforcing continuity between element. Calculate
 * field values on nodes of that mesh.
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

#include <MoFEM.hpp>
using namespace MoFEM;
using namespace boost::numeric;
#include <PostProcOnRefMesh.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif


PetscErrorCode PostProcCommonOnRefMesh::OpGetFieldValues::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  if(data.getFieldData().size()==0) PetscFunctionReturn(0);
  if(V) {
    vAlues.resize(data.getFieldData().size());
    double *a;
    ierr = VecGetArray(V,&a); CHKERRQ(ierr);
    VectorDofs::iterator it,hi_it;
    it = data.getFieldDofs().begin();
    hi_it = data.getFieldDofs().end();
    for(int ii = 0;it!=hi_it;it++,ii++) {
      int local_idx = getFEMethod()->rowPtr->find((*it)->getGlobalUniqueId())->get()->getPetscLocalDofIdx();
      vAlues[ii] = a[local_idx];
    }
    ierr = VecRestoreArray(V,&a); CHKERRQ(ierr);
    vAluesPtr = &vAlues;
  } else {
    vAluesPtr = &data.getFieldData();
  }

  const MoFEM::FEDofEntity *dof_ptr = data.getFieldDofs()[0];
  int rank = dof_ptr->getNbOfCoeffs();

  int tag_length = rank;
  FieldSpace space = dof_ptr->getSpace();
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

  if(tag_length>1 && tag_length < 3) {
    tag_length = 3;
  } else if(tag_length > 3 && tag_length < 9) {
    tag_length = 9;
  }

  double def_VAL[tag_length];
  bzero(def_VAL,tag_length*sizeof(double));
  Tag th;
  rval = postProcMesh.tag_get_handle(
    tagName.c_str(),tag_length,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL
  ); CHKERRQ_MOAB(rval);

  // zero tags, this for Vertex if H1 and TRI if Hdiv, EDGE for Hcurl
  // no need for L2
  const void* tags_ptr[mapGaussPts.size()];
  int nb_gauss_pts = data.getN().size1();
  if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
    SETERRQ2(
      PETSC_COMM_SELF,
      MOFEM_DATA_INCONSISTENCY,
      "data inconsistency %d!=%d",
      mapGaussPts.size(),nb_gauss_pts
    );
  }

  switch(space) {
    case H1:
    commonData.fieldMap[rowFieldName].resize(nb_gauss_pts);
    if(type == MBVERTEX) {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERRQ_MOAB(rval);
        (commonData.fieldMap[rowFieldName])[gg].resize(rank);
        (commonData.fieldMap[rowFieldName])[gg].clear();
      }
    }
    rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERRQ_MOAB(rval);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int rr = 0;rr<rank;rr++) {
        ((double*)tags_ptr[gg])[rr] += cblas_ddot(
          (vAluesPtr->size()/rank),&(data.getN(gg)[0]),1,&((*vAluesPtr)[rr]),rank
        );
        (commonData.fieldMap[rowFieldName])[gg][rr] += (*vAluesPtr)[rr];
      }
    }
    break;
    case L2:
    rval = postProcMesh.tag_get_by_ptr(
      th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr
    ); CHKERRQ_MOAB(rval);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      bzero((double*)tags_ptr[gg],sizeof(double)*tag_length);
      for(int rr = 0;rr<rank;rr++) {
        ((double*)tags_ptr[gg])[rr] = cblas_ddot(
          (vAluesPtr->size()/rank),&(data.getN(gg)[0]),1,&((*vAluesPtr)[rr]),rank
        );
      }
    }
    break;
    case HDIV:
    if(type == MBTRI && side == 0) {
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERRQ_MOAB(rval);
      }
    }
    rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERRQ_MOAB(rval);
    for(int gg = 0;gg<nb_gauss_pts;gg++) {
      for(int rr = 0;rr<rank;rr++) {
        for(int dd = 0;dd<3;dd++) {
          ((double*)tags_ptr[gg])[3*rr+dd] += cblas_ddot(
            (vAluesPtr->size()/rank),&(data.getHdivN(gg)(0,dd)),3,&((*vAluesPtr)[rr]),rank
          );

        }
      }
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcCommonOnRefMesh::OpGetFieldGradientValues::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  if(data.getFieldData().size()==0) PetscFunctionReturn(0);
  if(V) {
    vAlues.resize(data.getFieldData().size());
    double *a;
    ierr = VecGetArray(V,&a); CHKERRQ(ierr);
    VectorDofs::iterator it,hi_it;
    it = data.getFieldDofs().begin();
    hi_it = data.getFieldDofs().end();
    for(int ii = 0;it!=hi_it;it++,ii++) {
      int local_idx = getFEMethod()->rowPtr->find((*it)->getGlobalUniqueId())->get()->getPetscLocalDofIdx();
      vAlues[ii] = a[local_idx];
    }
    ierr = VecRestoreArray(V,&a); CHKERRQ(ierr);
    vAluesPtr = &vAlues;
  } else {
    vAluesPtr = &data.getFieldData();
  }

  const MoFEM::FEDofEntity *dof_ptr = data.getFieldDofs()[0];
  int rank = dof_ptr->getNbOfCoeffs();

  int tag_length = rank*3;
  FieldSpace space = dof_ptr->getSpace();
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
  rval = postProcMesh.tag_get_handle(tagName.c_str(),tag_length,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);

  // zero tags, this for Vertex if H1 and TRI if Hdiv, EDGE for Hcurl
  // no need for L2
  const void* tags_ptr[mapGaussPts.size()];
  int nb_gauss_pts = data.getN().size1();
  if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }

  try {

    switch(space) {
      case H1:
      commonData.gradMap[rowFieldName].resize(nb_gauss_pts);
      if(type == MBVERTEX) {
        for(int gg = 0;gg<nb_gauss_pts;gg++) {
          rval = postProcMesh.tag_set_data(th,&mapGaussPts[gg],1,def_VAL); CHKERRQ_MOAB(rval);
          (commonData.gradMap[rowFieldName])[gg].resize(rank,3);
          (commonData.gradMap[rowFieldName])[gg].clear();
        }
      }
      rval = postProcMesh.tag_get_by_ptr(th,&mapGaussPts[0],mapGaussPts.size(),tags_ptr); CHKERRQ_MOAB(rval);
      for(int gg = 0;gg<nb_gauss_pts;gg++) {
        for(int rr = 0;rr<rank;rr++) {
          for(int dd = 0;dd<3;dd++) {
            for(unsigned int dof = 0;dof<(vAluesPtr->size()/rank);dof++) {
              ((double*)tags_ptr[gg])[rank*rr+dd] += data.getDiffN(gg)(dof,dd)*(*vAluesPtr)[rank*dof+rr];
              (commonData.gradMap[rowFieldName])[gg](rr,dd) += data.getDiffN(gg)(dof,dd)*(*vAluesPtr)[rank*dof+rr];
            }
          }
        }
      }
      break;
      default:
      SETERRQ(PETSC_COMM_SELF,MOFEM_NOT_IMPLEMENTED,"field with that space is not implemented");
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::generateReferenceElementMesh() {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  int max_level = 0;
  if(nbOfRefLevels == -1) {
    PetscBool flg = PETSC_TRUE;
    PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_max_post_proc_ref_level",&max_level,&flg);
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
  moab::Interface& moab_ref = core_ref;

  EntityHandle nodes[4];
  for(int nn = 0;nn<4;nn++) {
    rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERRQ_MOAB(rval);
  }
  EntityHandle tet;
  rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERRQ_MOAB(rval);

  MoFEM::Core m_core_ref(moab_ref,PETSC_COMM_SELF,MB_TAG_DENSE,-2);
  MoFEM::Interface& m_field_ref = m_core_ref;

  ierr = m_field_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

  for(int ll = 0;ll<max_level;ll++) {
    PetscPrintf(mField.get_comm(),"Refine Level %d\n",ll);
    Range edges;
    ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(ll),BitRefLevel().set(),MBEDGE,edges); CHKERRQ(ierr);
    Range tets;
    ierr = m_field_ref.get_entities_by_type_and_ref_level(BitRefLevel().set(ll),BitRefLevel(ll).set(),MBTET,tets); CHKERRQ(ierr);
    //refine mesh
    MeshRefinment& m_ref = m_core_ref;
    ierr = m_ref.add_verices_in_the_middel_of_edges(edges,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
    ierr = m_ref.refine_TET(tets,BitRefLevel().set(ll+1)); CHKERRQ(ierr);
  }

  Range tets;
  ierr = m_field_ref.get_entities_by_type_and_ref_level(
    BitRefLevel().set(max_level),BitRefLevel().set(max_level),MBTET,tets
  ); CHKERRQ(ierr);

  if(tenNodesPostProcTets) {
    // Range edges;
    // rval = moab_ref.get_adjacencies(tets,1,true,edges); CHKERRQ_MOAB(rval);
    EntityHandle meshset;
    rval = moab_ref.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
    rval = moab_ref.add_entities(meshset,tets); CHKERRQ_MOAB(rval);
    rval = moab_ref.convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
    rval = moab_ref.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
  }

  Range elem_nodes;
  rval = moab_ref.get_connectivity(tets,elem_nodes,false); CHKERRQ_MOAB(rval);
  // ierr = m_field_ref.get_entities_by_type_and_ref_level(
  //   BitRefLevel().set(max_level),BitRefLevel().set(),MBVERTEX,elem_nodes
  // ); CHKERRQ(ierr);

  std::map<EntityHandle,int> little_map;
  gaussPts_FirstOrder.resize(elem_nodes.size(),4,0);
  Range::iterator nit = elem_nodes.begin();
  for(int gg = 0;nit!=elem_nodes.end();nit++,gg++) {
    rval = moab_ref.get_coords(&*nit,1,&gaussPts_FirstOrder(gg,0)); CHKERRQ_MOAB(rval);
    little_map[*nit] = gg;
  }
  gaussPts = gaussPts_FirstOrder;
  gaussPts = trans(gaussPts);

  Range::iterator tit = tets.begin();
  for(int tt = 0;tit!=tets.end();tit++,tt++) {
    const EntityHandle *conn;
    int num_nodes;
    rval = moab_ref.get_connectivity(*tit,conn,num_nodes,false); CHKERRQ_MOAB(rval);
    if(tt == 0) {
      refTets.resize(tets.size(),num_nodes);
    }
    for(int nn = 0;nn<num_nodes;nn++) {
      refTets(tt,nn) = little_map[conn[nn]];
    }
  }

  ublas::matrix<double> N;
  shapeFunctions.resize(elem_nodes.size(),4);
  ierr = ShapeMBTET(
    &*shapeFunctions.data().begin(),
    &gaussPts(0,0),
    &gaussPts(1,0),
    &gaussPts(2,0),
    elem_nodes.size()
  ); CHKERRQ(ierr);

  // EntityHandle meshset;
  // rval = moab_ref.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
  // rval = moab_ref.add_entities(meshset,tets); CHKERRQ_MOAB(rval);
  // rval = moab_ref.write_file("test_reference_mesh.vtk","VTK","",&meshset,1); CHKERRQ_MOAB(rval);
  //moab_ref.list_entities(tets);

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::setGaussPts(int order) {
  PetscFunctionBegin;

  try {

    ErrorCode rval;

    mapGaussPts.resize(gaussPts_FirstOrder.size1());
    for(unsigned int gg = 0;gg<gaussPts_FirstOrder.size1();gg++) {
      rval = postProcMesh.create_vertex(
        &gaussPts_FirstOrder(gg,0),mapGaussPts[gg]
      ); CHKERRQ_MOAB(rval);
    }

    commonData.tEts.clear();
    for(unsigned int tt = 0;tt<refTets.size1();tt++) {
      int num_nodes = refTets.size2();
      EntityHandle conn[num_nodes];
      for(int nn = 0;nn!=num_nodes;nn++) {
        conn[nn] = mapGaussPts[refTets(tt,nn)];
      }
      EntityHandle tet;
      rval = postProcMesh.create_element(MBTET,conn,num_nodes,tet); CHKERRQ_MOAB(rval);
      commonData.tEts.insert(tet);
    }

    EntityHandle fe_ent = numeredEntFiniteElementPtr->getEnt();
    coords.resize(12,false);
    {
      const EntityHandle *conn;
      int num_nodes;
      mField.get_moab().get_connectivity(fe_ent,conn,num_nodes,true);
      // coords.resize(3*num_nodes,false);
      rval = mField.get_moab().get_coords(conn,num_nodes,&coords[0]); CHKERRQ_MOAB(rval);
    }

    Range nodes;
    rval = postProcMesh.get_connectivity(commonData.tEts,nodes,false); CHKERRQ_MOAB(rval);

    coordsAtGaussPts.resize(nodes.size(),3,false);
    for(unsigned int gg = 0;gg<nodes.size();gg++) {
      for(int dd = 0;dd<3;dd++) {
        coordsAtGaussPts(gg,dd) = cblas_ddot(4,&shapeFunctions(gg,0),1,&coords[dd],3);
      }
    }

    mapGaussPts.resize(nodes.size());
    Range::iterator nit = nodes.begin();
    for(int gg = 0;nit!=nodes.end();nit++,gg++) {
      rval = postProcMesh.set_coords(&*nit,1,&coordsAtGaussPts(gg,0)); CHKERRQ_MOAB(rval);
      mapGaussPts[gg] = *nit;
    }

  } catch (std::exception& ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::clearOperators() {
  PetscFunctionBegin;
  getOpPtrVector().clear();
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::preProcess() {
  PetscFunctionBegin;
  ErrorCode rval;
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh != NULL) {
    delete pcomm_post_proc_mesh;
  }
  rval = postProcMesh.delete_mesh(); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::postProcess() {
  PetscFunctionBegin;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh == NULL) {
    pcomm_post_proc_mesh = new ParallelComm(&postProcMesh,mField.get_comm());
  }

  Range edges;
  rval = postProcMesh.get_entities_by_type(0,MBEDGE,edges,false);  CHKERRQ_MOAB(rval);
  rval = postProcMesh.delete_entities(edges); CHKERRQ_MOAB(rval);
  Range tris;
  rval = postProcMesh.get_entities_by_type(0,MBTRI,tris,false);  CHKERRQ_MOAB(rval);
  rval = postProcMesh.delete_entities(tris); CHKERRQ_MOAB(rval);

  Range tets;
  rval = postProcMesh.get_entities_by_type(0,MBTET,tets,false);  CHKERRQ_MOAB(rval);

  //std::cerr << "total tets size " << tets.size() << std::endl;

  int rank = pcomm->rank();
  Range::iterator tit = tets.begin();
  for(;tit!=tets.end();tit++) {
    rval = postProcMesh.tag_set_data(pcomm_post_proc_mesh->part_tag(),&*tit,1,&rank); CHKERRQ_MOAB(rval);
  }

  rval = pcomm_post_proc_mesh->resolve_shared_ents(0); CHKERRQ_MOAB(rval);

  // #ifndef MOAB_HDF5_PARALLEL
  // #warning "No parallel HDF5, not most efficient way of writing files"
  // for(int r = 0;r<pcomm_post_proc_mesh->size();r++) {
  //   // FIXME make better communication send only to proc 0
  //   rval = pcomm_post_proc_mesh->broadcast_entities(r,tets,false,true); CHKERRQ_MOAB(rval);
  // }
  // #endif //

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::writeFile(const std::string file_name) {
 PetscFunctionBegin;
 MoABErrorCode rval;
 // #ifdef MOAB_HDF5_PARALLEL
  rval = postProcMesh.write_file(file_name.c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERRQ_MOAB(rval);
 // #else
 //  #warning "No parallel HDF5, not most efficient way of writing files"
 //  if(mField.getCommRank()==0) {
 //    rval = postProcMesh.write_file(file_name.c_str(),"MOAB",""); CHKERRQ_MOAB(rval);
 //  }
 // #endif
 PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::OpHdivFunctions::doWork(
  int side,
  EntityType type,
  DataForcesAndSurcesCore::EntData &data
) {
  PetscFunctionBegin;

  if(data.getIndices().size()==0) PetscFunctionReturn(0);

  ErrorCode rval;
  PetscErrorCode ierr;

  std::vector<Tag> th;
  th.resize(data.getFieldData().size());

  double def_VAL[9] = { 0,0,0 };

  switch(type) {
    case MBTRI:
    for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
      std::ostringstream ss;
      ss << "HDIV_FACE_" << side << "_" << dd;
      rval = postProcMesh.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,th[dd],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);
    }
    break;
    case MBTET:
    for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
      std::ostringstream ss;
      ss << "HDIV_TET_" << dd;
      rval = postProcMesh.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,th[dd],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);
    }
    break;
    default:
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
  }

  for(unsigned int gg = 0;gg<data.getHdivN().size1();gg++) {
    for(unsigned int dd = 0;dd<data.getHdivN().size2()/3;dd++) {
      ierr = postProcMesh.tag_set_data(th[dd],&mapGaussPts[gg],1,&data.getHdivN(gg)(dd,0)); CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcVolumeOnRefinedMesh::addHdivFunctionsPostProc(const std::string field_name) {
  PetscFunctionBegin;
  getOpPtrVector().push_back(new OpHdivFunctions(postProcMesh,mapGaussPts,field_name));
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFatPrismOnRefinedMesh::generateReferenceElementMesh() {
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

struct PointsMap3D {
  const int kSi;
  const int eTa;
  const int zEta;
  int nN;
  PointsMap3D(
    const int ksi,
    const int eta,
    const int zeta,
    const int nn
  ):
  kSi(ksi),
  eTa(eta),
  zEta(zeta),
  nN(nn) {
  }
};
typedef multi_index_container<
  PointsMap3D,
  indexed_by<
    ordered_unique<
      composite_key<
      PointsMap3D,
      member<PointsMap3D,const int,&PointsMap3D::kSi>,
      member<PointsMap3D,const int,&PointsMap3D::eTa>,
      member<PointsMap3D,const int,&PointsMap3D::zEta>
    >
  >
> > PointsMap3D_multiIndex;

PetscErrorCode PostProcFatPrismOnRefinedMesh::setGaussPtsTrianglesOnly(int order_triangles_only) {
  PetscFunctionBegin;
  // if(gaussPtsTrianglesOnly.size1()==0 || gaussPtsTrianglesOnly.size2()==0) {
  //   SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"post-process mesh not generated");
  // }

  //FIXME: Refinement not implement
  //FIXME: This is inefficient implementation

  // PetscErrorCode ierr;
  MoABErrorCode rval;
  const EntityHandle *conn;
  int num_nodes;
  EntityHandle prism;

  if(elementsMap.find(numeredEntFiniteElementPtr->getEnt())!=elementsMap.end()) {
    prism = elementsMap[numeredEntFiniteElementPtr->getEnt()];
  } else {

    {
      gaussPtsTrianglesOnly.resize(3,3,false);
      gaussPtsTrianglesOnly.clear();
      gaussPtsTrianglesOnly(0,0) = 0;
      gaussPtsTrianglesOnly(1,0) = 0;
      gaussPtsTrianglesOnly(0,1) = 1;
      gaussPtsTrianglesOnly(1,1) = 0;
      gaussPtsTrianglesOnly(0,2) = 0;
      gaussPtsTrianglesOnly(1,2) = 1;
      gaussPtsThroughThickness.resize(2,2,false);
      gaussPtsThroughThickness(0,0) = 0;
      gaussPtsThroughThickness(0,1) = 1;
      mapGaussPts.clear();
    }

    ublas::vector<EntityHandle> prism_conn(6);
    ublas::vector<double> coords(3);
    for(int ggf = 0;ggf!=3;ggf++) {
      double ksi = gaussPtsTrianglesOnly(0,ggf);
      double eta = gaussPtsTrianglesOnly(1,ggf);
      coords[0] = ksi;
      coords[1] = eta;
      for(int ggt = 0;ggt!=2;ggt++) {
        double zeta = gaussPtsThroughThickness(0,ggt);
        coords[2] = zeta;
        int side = ggt*3+ggf;
        rval = postProcMesh.create_vertex(
          &coords[0],prism_conn[side]
        ); CHKERRQ_MOAB(rval);
      }
    }
    rval = postProcMesh.create_element(MBPRISM,&prism_conn[0],6,prism); CHKERRQ_MOAB(rval);

    elementsMap[numeredEntFiniteElementPtr->getEnt()] = prism;
    // Range faces;
    // rval = postProcMesh.get_adjacencies(&prism,1,2,true,faces); CHKERRQ_MOAB(rval);
    Range edges;
    rval = postProcMesh.get_adjacencies(&prism,1,1,true,edges); CHKERRQ_MOAB(rval);
    EntityHandle meshset;
    rval = postProcMesh.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
    rval = postProcMesh.add_entities(meshset,&prism,1); CHKERRQ_MOAB(rval);
    // rval = postProcMesh.add_entities(meshset,faces); CHKERRQ_MOAB(rval);
    rval = postProcMesh.add_entities(meshset,edges); CHKERRQ_MOAB(rval);
    if(tenNodesPostProcTets) {
      rval = postProcMesh.convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
    }
    rval = postProcMesh.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    rval = postProcMesh.delete_entities(edges); CHKERRQ_MOAB(rval);
    // rval = postProcMesh.delete_entities(faces); CHKERRQ_MOAB(rval);

    rval = mField.get_moab().get_connectivity(
      numeredEntFiniteElementPtr->getEnt(),conn,num_nodes,true
    ); CHKERRQ_MOAB(rval);
    ublas::matrix<double> coords_prism(num_nodes,3);
    rval = mField.get_moab().get_coords(conn,num_nodes,&coords_prism(0,0));

    rval = postProcMesh.get_connectivity(prism,conn,num_nodes,false); CHKERRQ_MOAB(rval);
    ublas::matrix<double> coords_prism_local(num_nodes,3);
    rval = postProcMesh.get_coords(conn,num_nodes,&coords_prism_local(0,0));
    ublas::matrix<double> coords_prism_global(num_nodes,3);

    const double eps = 1e-6;

    int nb_on_triangle = 0;
    int nb_through_thickness = 0;
    for(int nn = 0;nn!=num_nodes;nn++) {
      if(fabs(coords_prism_local(nn,2))<eps) {
        nb_on_triangle++;
      }
      if(
        fabs(coords_prism_local(nn,0))<eps&&
        fabs(coords_prism_local(nn,1))<eps
      ) {
        nb_through_thickness++;
      }
    }

    PointsMap3D_multiIndex points_map;

    gaussPtsTrianglesOnly.resize(3,nb_on_triangle,false);
    gaussPtsTrianglesOnly.clear();
    gaussPtsThroughThickness.resize(2,nb_through_thickness,false);
    gaussPtsThroughThickness.clear();
    {
      int ggf = 0,ggt = 0;
      for(int nn = 0;nn!=num_nodes;nn++) {
        double ksi = coords_prism_local(nn,0);
        double eta = coords_prism_local(nn,1);
        double zeta = coords_prism_local(nn,2);
        points_map.insert(PointsMap3D(100*ksi,100*eta,100*zeta,nn));
        if(fabs(zeta)< eps) {
          gaussPtsTrianglesOnly(0,ggf) = ksi;
          gaussPtsTrianglesOnly(1,ggf) = eta;
          ggf++;
        }
        if(fabs(ksi)<eps&&fabs(eta)<eps) {
          gaussPtsThroughThickness(0,ggt) = zeta;
          ggt++;
        }
        double n0 = N_MBTRI0(ksi,eta);
        double n1 = N_MBTRI1(ksi,eta);
        double n2 = N_MBTRI2(ksi,eta);
        double e0 = N_MBEDGE0(zeta);
        double e1 = N_MBEDGE1(zeta);
        for(int dd = 0;dd!=3;dd++) {
          coords_prism_global(nn,dd) =
          n0*e0*coords_prism(0,dd)+
          n1*e0*coords_prism(1,dd)+
          n2*e0*coords_prism(2,dd)+
          n0*e1*coords_prism(3,dd)+
          n1*e1*coords_prism(4,dd)+
          n2*e1*coords_prism(5,dd);
        }
        rval = postProcMesh.set_coords(
          &conn[nn],1,&coords_prism_global(nn,0)
        ); CHKERRQ_MOAB(rval);
      }
    }

    mapGaussPts.resize(nb_through_thickness*nb_on_triangle);
    std::fill(mapGaussPts.begin(),mapGaussPts.end(),0);
    {
      int gg = 0;
      for(unsigned int ggf = 0;ggf!=gaussPtsTrianglesOnly.size2();ggf++) {
        const int ksi = 100*gaussPtsTrianglesOnly(0,ggf);
        const int eta = 100*gaussPtsTrianglesOnly(1,ggf);
        for(unsigned int ggt = 0;ggt!=gaussPtsThroughThickness.size2();ggt++,gg++) {
          const int zeta = 100*gaussPtsThroughThickness(0,ggt);
          PointsMap3D_multiIndex::iterator it;
          it = points_map.find(boost::make_tuple(ksi,eta,zeta));
          if(it != points_map.end()) {
            mapGaussPts[gg] = conn[it->nN];
          };
        }
      }
    }

  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFatPrismOnRefinedMesh::setGaussPtsThroughThickness(int order_thickness) {
  PetscFunctionBegin;
  if(gaussPtsThroughThickness.size1()==0 || gaussPtsThroughThickness.size2()==0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"post-process mesh not generated");
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFatPrismOnRefinedMesh::preProcess() {
  PetscFunctionBegin;
  // MoABErrorCode rval;
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh != NULL) {
    delete pcomm_post_proc_mesh;
  }
  // rval = postProcMesh.delete_mesh(); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFatPrismOnRefinedMesh::postProcess() {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh == NULL) {
    pcomm_post_proc_mesh = new ParallelComm(&postProcMesh,mField.get_comm());
  }
  Range prims;
  rval = postProcMesh.get_entities_by_type(0,MBPRISM,prims,false);  CHKERRQ_MOAB(rval);
  //std::cerr << "total prims size " << prims.size() << std::endl;
  int rank = pcomm->rank();
  Range::iterator pit = prims.begin();
  for(;pit!=prims.end();pit++) {
    rval = postProcMesh.tag_set_data(
      pcomm_post_proc_mesh->part_tag(),&*pit,1,&rank
    ); CHKERRQ_MOAB(rval);
  }
  rval = pcomm_post_proc_mesh->resolve_shared_ents(0); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}


PetscErrorCode PostProcFaceOnRefinedMesh::generateReferenceElementMesh() {
  PetscFunctionBegin;

  gaussPts.resize(3,3,false);
  gaussPts.clear();
  gaussPts(0,0) = 0;
  gaussPts(1,0) = 0;
  gaussPts(0,1) = 1;
  gaussPts(1,1) = 0;
  gaussPts(0,2) = 0;
  gaussPts(1,2) = 1;
  mapGaussPts.resize(gaussPts.size2());

  moab::Core core_ref;
  moab::Interface& moab_ref = core_ref;
  const EntityHandle *conn;
  int num_nodes;
  EntityHandle tri_conn[3];
  MatrixDouble coords(6,3);
  for(int gg = 0;gg!=3;gg++) {
    coords(gg,0) = gaussPts(0,gg);
    coords(gg,1) = gaussPts(1,gg);
    coords(gg,2) = 0;
    rval = moab_ref.create_vertex(&coords(gg,0),tri_conn[gg]); CHKERRQ_MOAB(rval);
  }

  EntityHandle tri;
  rval = moab_ref.create_element(MBTRI,tri_conn,3,tri); CHKERRQ_MOAB(rval);
  Range edges;
  rval = moab_ref.get_adjacencies(&tri,1,1,true,edges); CHKERRQ_MOAB(rval);
  EntityHandle meshset;
  rval = moab_ref.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  rval = moab_ref.add_entities(meshset,&tri,1); CHKERRQ_MOAB(rval);
  rval = moab_ref.add_entities(meshset,edges); CHKERRQ_MOAB(rval);
  if(sixNodePostProcTris) {
    rval = moab_ref.convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
  }
  rval = moab_ref.get_connectivity(tri,conn,num_nodes,false); CHKERRQ_MOAB(rval);
  rval = moab_ref.get_coords(conn,num_nodes,&coords(0,0));

  gaussPts.resize(3,num_nodes,false);
  gaussPts.clear();
  for(int nn = 0;nn<3;nn++) {
    gaussPts(0,nn) = coords(nn,0);
    gaussPts(1,nn) = coords(nn,1);
    gaussPts(0,3+nn) = coords(3+nn,0);
    gaussPts(1,3+nn) = coords(3+nn,1);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFaceOnRefinedMesh::setGaussPts(int order) {
  PetscFunctionBegin;
  if(gaussPts.size1()==0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"post-process mesh not generated");
  }

  // PetscErrorCode ierr;
  MoABErrorCode rval;
  const EntityHandle *conn;
  int num_nodes;
  EntityHandle tri;

  if(elementsMap.find(numeredEntFiniteElementPtr->getEnt())!=elementsMap.end()) {
    tri = elementsMap[numeredEntFiniteElementPtr->getEnt()];
  } else {
    ublas::vector<EntityHandle> tri_conn(3);
    ublas::matrix<double> coords_tri(3,3);
    ublas::vector<double> coords(3);
    rval = mField.get_moab().get_connectivity(numeredEntFiniteElementPtr->getEnt(),conn,num_nodes,true); CHKERRQ_MOAB(rval);
    rval = mField.get_moab().get_coords(conn,num_nodes,&coords_tri(0,0));
    for(int gg = 0;gg!=3;gg++) {
      double ksi = gaussPts(0,gg);
      double eta = gaussPts(1,gg);
      double n0 = N_MBTRI0(ksi,eta);
      double n1 = N_MBTRI1(ksi,eta);
      double n2 = N_MBTRI2(ksi,eta);
      double x = n0*coords_tri(0,0)+n1*coords_tri(1,0)+n2*coords_tri(2,0);
      double y = n0*coords_tri(0,1)+n1*coords_tri(1,1)+n2*coords_tri(2,1);
      coords[0] = x;
      coords[1] = y;
      coords[2] = 0;
      rval = postProcMesh.create_vertex(
        &coords[0],tri_conn[gg]
      ); CHKERRQ_MOAB(rval);
    }
    rval = postProcMesh.create_element(MBTRI,&tri_conn[0],3,tri); CHKERRQ_MOAB(rval);
    elementsMap[numeredEntFiniteElementPtr->getEnt()] = tri;
    Range edges;
    rval = postProcMesh.get_adjacencies(&tri,1,1,true,edges); CHKERRQ_MOAB(rval);
    EntityHandle meshset;
    rval = postProcMesh.create_meshset(MESHSET_SET,meshset); CHKERRQ_MOAB(rval);
    rval = postProcMesh.add_entities(meshset,&tri,1); CHKERRQ_MOAB(rval);
    rval = postProcMesh.add_entities(meshset,edges); CHKERRQ_MOAB(rval);
    if(sixNodePostProcTris) {
      rval = postProcMesh.convert_entities(meshset,true,false,false); CHKERRQ_MOAB(rval);
    }
    rval = postProcMesh.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    rval = postProcMesh.delete_entities(edges); CHKERRQ_MOAB(rval);
  }

  // Set values which map nodes with integration points on the prism
  {
    rval = postProcMesh.get_connectivity(tri,conn,num_nodes,false); CHKERRQ_MOAB(rval);
    mapGaussPts.resize(num_nodes);
    for(int nn = 0;nn<num_nodes;nn++) {
      mapGaussPts[nn] = conn[nn];
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFaceOnRefinedMesh::preProcess() {
  PetscFunctionBegin;
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh != NULL) {
    delete pcomm_post_proc_mesh;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PostProcFaceOnRefinedMesh::postProcess() {
  PetscFunctionBegin;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
  if(pcomm_post_proc_mesh == NULL) {
    pcomm_post_proc_mesh = new ParallelComm(&postProcMesh,mField.get_comm());
  }
  Range tris;
  rval = postProcMesh.get_entities_by_type(0,MBTRI,tris,false);  CHKERRQ_MOAB(rval);
  int rank = pcomm->rank();
  Range::iterator pit = tris.begin();
  for(;pit!=tris.end();pit++) {
    rval = postProcMesh.tag_set_data(
      pcomm_post_proc_mesh->part_tag(),&*pit,1,&rank
    ); CHKERRQ_MOAB(rval);
  }
  rval = pcomm_post_proc_mesh->resolve_shared_ents(0); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}
