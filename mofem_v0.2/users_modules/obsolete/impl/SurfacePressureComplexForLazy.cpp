/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
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


#include <MoFEM.hpp>

using namespace MoFEM;
#include <SurfacePressureComplexForLazy.hpp>

extern "C" {
  #include "complex_for_lazy.h"
  void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
    double circumcenter[3],double *xi,double *eta,double *zeta);
  void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
    double circumcenter[3],double *xi,double *eta);
}

namespace ObosleteUsersModules {

PetscErrorCode NeummanForcesSurfaceComplexForLazy::
  AuxMethodSpatial::doWork(int side, EntityType type, DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  try {

  switch (type) {
    case MBVERTEX: {
      if(data.getFieldData().size()!=9) {
	SETERRQ2(PETSC_COMM_SELF,1,"it should be 9 dofs on vertices but is %d of field < %s >",
	  data.getFieldData().size(),row_field_name.c_str());
      }
      myPtr->N = &*data.getN().data().begin();
      myPtr->diffN = &*data.getDiffN().data().begin();
      myPtr->dOfs_x.resize(data.getFieldData().size());
      ublas::noalias(myPtr->dOfs_x) = data.getFieldData();
      myPtr->dofs_x = &*myPtr->dOfs_x.data().begin();
      myPtr->dOfs_x_indices.resize(data.getIndices().size());
      ublas::noalias(myPtr->dOfs_x_indices) = data.getIndices();
      myPtr->dofs_x_indices = &*myPtr->dOfs_x_indices.data().begin();
    }
    break;
    case MBEDGE: {
      myPtr->order_edge[side] = data.getOrder();
      myPtr->N_edge[side] = &*data.getN().data().begin();
      myPtr->diffN_edge[side] = &*data.getDiffN().data().begin();
      myPtr->dOfs_x_edge.resize(3);
      myPtr->dOfs_x_edge[side].resize(data.getFieldData().size());
      myPtr->dofs_x_edge[side] = &*myPtr->dOfs_x_edge[side].data().begin();
      myPtr->dOfs_x_edge_indices.resize(3);
      myPtr->dOfs_x_edge_indices[side].resize(data.getIndices().size());
      ublas::noalias(myPtr->dOfs_x_edge_indices[side]) = data.getIndices();
      myPtr->dofs_x_edge_indices[side] = &*myPtr->dOfs_x_edge_indices[side].data().begin();
    }
    break;
    case MBTRI: {
      myPtr->order_face = data.getOrder();
      myPtr->N_face = &*data.getN().data().begin();
      myPtr->diffN_face = &*data.getDiffN().data().begin();
      myPtr->dOfs_x_face.resize(data.getFieldData().size());
      ublas::noalias(myPtr->dOfs_x_face) = data.getFieldData();
      myPtr->dofs_x_face = &*myPtr->dOfs_x_face.data().begin();
      myPtr->dOfs_x_face_indices.resize(data.getIndices().size());
      ublas::noalias(myPtr->dOfs_x_face_indices) = data.getIndices();
      myPtr->dofs_x_face_indices = &*myPtr->dOfs_x_face_indices.data().begin();
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"unknown entity type");
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "side: " << side << " type: " << type << endl;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::
  AuxMethodMaterial::doWork(int side, EntityType type, DataForcesAndSurcesCore::EntData &data) {
  PetscFunctionBegin;

  //cerr << "AuxMethodMaterial\n";

  try {

  switch (type) {
    case MBVERTEX: {
      if(data.getFieldData().size()!=9) {
	SETERRQ(PETSC_COMM_SELF,1,"it should be 9 dofs on vertices");
      }
      if(data.getN().size2()!=3) {
	SETERRQ(PETSC_COMM_SELF,1,"it shoule 3 shape functiond for 3 nodes");
      }
      myPtr->N = &*data.getN().data().begin();
      myPtr->diffN = &*data.getDiffN().data().begin();
      myPtr->dOfs_X_indices.resize(data.getIndices().size());
      ublas::noalias(myPtr->dOfs_X_indices) = data.getIndices();
      myPtr->dofs_X_indices = &*myPtr->dOfs_X_indices.data().begin();
      myPtr->dOfs_X.resize(data.getFieldData().size());
      ublas::noalias(myPtr->dOfs_X) = data.getFieldData();
      myPtr->dofs_X = &*myPtr->dOfs_x.data().begin();
    }
    break;
    case MBEDGE: {	
      myPtr->order_edge_material[side] = data.getOrder();
      myPtr->dOfs_X_edge.resize(3);
      myPtr->dOfs_X_edge[side].resize(data.getFieldData().size());
      ublas::noalias(myPtr->dOfs_X_edge[side]) = data.getFieldData();
      myPtr->dofs_X_edge[side] = &*myPtr->dOfs_X_edge[side].data().begin();
    }
    break;
    case MBTRI: {
      myPtr->order_face_material = data.getOrder();
      myPtr->dOfs_X_face.resize(data.getFieldData().size());
      ublas::noalias(myPtr->dOfs_X_face) = data.getFieldData();
      myPtr->dofs_X_face = &*myPtr->dOfs_X_face.data().begin();
    }
    break;
    default:
      SETERRQ(PETSC_COMM_SELF,1,"unknown entity type");
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::rHs() {
  PetscFunctionBegin;

  //cerr << "MyTriangleSpatialFE::rHs\n";

  try {

  fExtNode.resize(9);	
  fExtFace.resize(dataH1.dataOnEntities[MBTRI][0].getFieldData().size());
  fExtEdge.resize(3);
  for(int ee = 0;ee<3;ee++) {
    int nb_edge_dofs = dOfs_x_edge_indices[ee].size();
    if(nb_edge_dofs > 0) {
      fExtEdge[ee].resize(nb_edge_dofs);
      Fext_edge[ee] = &*fExtEdge[ee].data().begin();
    } else {
      Fext_edge[ee] = NULL;
    }
  }
    
  //cerr << "dOfs_x: " << dOfs_x << endl;
  //for(int ee = 0;ee<3;ee++) {
    //cerr << dOfs_x_edge[ee] << endl;
  //}
    
  switch(typeOfForces) {
    case CONSERVATIVE:
      ierr = Fext_h_hierarchical(
	order_face,order_edge,//2
	N,N_face,N_edge,diffN,diffN_face,diffN_edge,//8
	t_loc,NULL,NULL,//11
	dofs_x,dofs_x_edge,dofs_x_face,//14
	NULL,NULL,NULL,//17
	&*fExtNode.data().begin(),Fext_edge,&*fExtFace.data().begin(),//20
	NULL,NULL,NULL,//23
	gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
      break;
    case NONCONSERVATIVE:
      ierr = Fext_h_hierarchical(
	order_face,order_edge,//2
	N,N_face,N_edge,diffN,diffN_face,diffN_edge,//8
	t_loc,NULL,NULL,//11
	dofs_X,dofs_X_edge,dofs_X_face,//14
	NULL,NULL,NULL,//17
	&*fExtNode.data().begin(),Fext_edge,&*fExtFace.data().begin(),//20
	NULL,NULL,NULL,//23
	gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
      break;
  }

  //cerr << "fExtNode: " << fExtNode << endl;
  //cerr << "fExtFace: " << fExtFace << endl;
  //for(int ee = 0;ee<3;ee++) {
    //cerr << "fExtEdge " << ee << " " << fExtEdge[ee] << endl;
  //}

 
  Vec f = snes_f;
  if(uSeF) f = F; 

  ierr = VecSetValues(f,
    9,dofs_x_indices,
    &*fExtNode.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  if(dOfs_x_face_indices.size()>0) {
    ierr = VecSetValues(f,
      dOfs_x_face_indices.size(),dofs_x_face_indices,
      &*fExtFace.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  }
  for(int ee = 0;ee<3;ee++) {
    if(dOfs_x_edge_indices[ee].size()>0) {
      ierr = VecSetValues(f,
	dOfs_x_edge_indices[ee].size(),dofs_x_edge_indices[ee],
	Fext_edge[ee],ADD_VALUES); CHKERRQ(ierr);
    }
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::lHs() {
  PetscFunctionBegin;

  try {

  double center[3];
  tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
  cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
  double r = cblas_dnrm2(3,center,1);

  kExtNodeNode.resize(9,9);
  kExtEdgeNode.resize(3);
  for(int ee = 0;ee<3;ee++) {
    kExtEdgeNode[ee].resize(dOfs_x_edge_indices[ee].size(),9);
    Kext_edge_node[ee] = &*kExtEdgeNode[ee].data().begin();
  }
  kExtFaceNode.resize(dOfs_x_face_indices.size(),9);
  ierr = KExt_hh_hierarchical(
    r*eps,order_face,order_edge,
    N,N_face,N_edge,diffN,diffN_face,diffN_edge,
    t_loc,NULL,NULL,
    dofs_x,dofs_x_edge,dofs_x_face,
    &*kExtNodeNode.data().begin(),Kext_edge_node,&*kExtFaceNode.data().begin(),
    gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
  //cerr << kExtNodeNode << endl;
  ierr = MatSetValues(snes_B,
    9,dofs_x_indices,
    9,dofs_x_indices,
    &*kExtNodeNode.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  ierr = MatSetValues(snes_B,
    kExtFaceNode.size1(),dofs_x_face_indices,
    9,dofs_x_indices,
    &*kExtFaceNode.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  //cerr << kExtFaceNode << endl;
  for(int ee = 0;ee<3;ee++) {
    //cerr << kExtEdgeNode[ee] << endl;
    ierr = MatSetValues(snes_B,
      kExtEdgeNode[ee].size1(),dofs_x_edge_indices[ee],
      9,dofs_x_indices,
      Kext_edge_node[ee],ADD_VALUES); CHKERRQ(ierr);
  }

  kExtNodeFace.resize(9,dOfs_x_face_indices.size());
  kExtEdgeFace.resize(3);
  for(int ee = 0;ee<3;ee++) {
    kExtEdgeFace[ee].resize(dOfs_x_edge_indices[ee].size(),dataH1.dataOnEntities[MBTRI][0].getIndices().size());
    Kext_edge_face[ee] = &*kExtEdgeFace[ee].data().begin();
  }
  kExtFaceFace.resize(dOfs_x_face_indices.size(),dOfs_x_face_indices.size());
  ierr = KExt_hh_hierarchical_face(
    r*eps,order_face,order_edge,
    N,N_face,N_edge,diffN,diffN_face,diffN_edge,
    t_loc,NULL,NULL,
    dofs_x,dofs_x_edge,dofs_x_face,
    &*kExtNodeFace.data().begin(),Kext_edge_face,&*kExtFaceFace.data().begin(),
    gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
  //cerr << "kExtNodeFace " << kExtNodeFace << endl;
  //cerr << "kExtFaceFace " << kExtFaceFace << endl;
  ierr = MatSetValues(snes_B,
    9,dofs_x_indices,
    kExtNodeFace.size2(),dofs_x_face_indices,
    &*kExtNodeFace.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  ierr = MatSetValues(snes_B,
    kExtFaceFace.size1(),dofs_x_face_indices,
    kExtFaceFace.size2(),dofs_x_face_indices,
    &*kExtFaceFace.data().begin(),ADD_VALUES); CHKERRQ(ierr);
  for(int ee = 0;ee<3;ee++) {
    //cerr << "kExtEdgeFace " << kExtEdgeFace[ee] << endl;
    ierr = MatSetValues(snes_B,
      kExtEdgeFace[ee].size1(),dofs_x_edge_indices[ee],
      kExtFaceFace.size2(),dofs_x_face_indices,
      Kext_edge_face[ee],ADD_VALUES); CHKERRQ(ierr);
  }

  kExtFaceEdge.resize(3);
  kExtNodeEdge.resize(3);
  kExtEdgeEdge.resize(3,3);
  for(int ee = 0;ee<3;ee++) {
    if(dOfs_x_edge_indices[ee].size()!=(unsigned int)(3*NBEDGE_H1(order_edge[ee]))) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    }
    kExtFaceEdge[ee].resize(dOfs_x_face_indices.size(),dOfs_x_edge_indices[ee].size());
    kExtNodeEdge[ee].resize(9,dOfs_x_edge_indices[ee].size());
    Kext_node_edge[ee] = &*kExtNodeEdge[ee].data().begin();
    Kext_face_edge[ee] = &*kExtFaceEdge[ee].data().begin();
    for(int EE = 0;EE<3;EE++) {
      kExtEdgeEdge(EE,ee).resize(dOfs_x_edge_indices[EE].size(),dOfs_x_edge_indices[ee].size());
      Kext_edge_edge[EE][ee] = &*kExtEdgeEdge(EE,ee).data().begin();
    }
  }
  ierr = KExt_hh_hierarchical_edge(
    r*eps,order_face,order_edge,
    N,N_face,N_edge,diffN,diffN_face,diffN_edge,
    t_loc,NULL,NULL,
    dofs_x,dofs_x_edge,dofs_x_face,
    Kext_node_edge,Kext_edge_edge,Kext_face_edge,
    gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
  for(int ee = 0;ee<3;ee++) {
    //cerr << "kExtFaceEdge: " << kExtFaceEdge[ee] << endl;
    //cerr << "kExtFaceEdge: " << kExtNodeEdge[ee] << endl;
    ierr = MatSetValues(snes_B,
      kExtFaceEdge[ee].size1(),dofs_x_face_indices,
      kExtFaceEdge[ee].size2(),dofs_x_edge_indices[ee],
      &*kExtFaceEdge[ee].data().begin(),ADD_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(snes_B,
      9,dofs_x_indices,
      kExtNodeEdge[ee].size2(),dofs_x_edge_indices[ee],
      &*kExtNodeEdge[ee].data().begin(),ADD_VALUES); CHKERRQ(ierr);
    for(int EE = 0;EE<3;EE++) {
      //cerr << kExtEdgeEdge(EE,ee) << endl;
      ierr = MatSetValues(snes_B,
	kExtEdgeEdge(EE,ee).size1(),dofs_x_edge_indices[EE],
	kExtEdgeEdge(EE,ee).size2(),dofs_x_edge_indices[ee],
	Kext_edge_edge[EE][ee],ADD_VALUES); CHKERRQ(ierr);
    }
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::reBaseToFaceLoocalCoordSystem(ublas::matrix<double> &t_glob_nodal) {
  PetscFunctionBegin;
  double s1[3],s2[3],normal[3],q[9];
  ierr = ShapeFaceBaseMBTRI(diffN,&*coords.data().begin(),normal,s1,s2); CHKERRQ(ierr);
  double nrm2_normal = cblas_dnrm2(3,normal,1);
  cblas_dscal(3,1./nrm2_normal,normal,1);
  cblas_dcopy(3,s1,1,&q[0],1);
  cblas_dcopy(3,s2,1,&q[3],1);
  cblas_dcopy(3,normal,1,&q[6],1);
  __CLPK_integer info;
  __CLPK_integer ipiv[3];
  info = lapack_dgesv(3,3,q,3,ipiv,&*t_glob_nodal.data().begin(),3);
  if(info != 0) {
    SETERRQ1(PETSC_COMM_SELF,1,"error solve dgesv info = %d",info);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::calcTraction() {
  PetscFunctionBegin;

  try {

  EntityHandle ent = fePtr->get_ent();
  map<int,bCPreassure>::iterator mip = mapPreassure.begin();
  tLoc.resize(3);
  tLoc[0] = tLoc[1] = tLoc[2] = 0;
  for(;mip!=mapPreassure.end();mip++) {
    if(mip->second.tRis.find(ent)!=mip->second.tRis.end()) {
      tLoc[2] -= mip->second.data.data.value1;
    }
  }
  tLocNodal.resize(3,3);
  for(int nn = 0;nn<3;nn++) {
    for(int dd = 0;dd<3;dd++) {
      tLocNodal(nn,dd) = tLoc[dd];
    }
  }

  map<int,bCForce>::iterator mif = mapForce.begin();
  for(;mif!=mapForce.end();mif++) {
    if(mif->second.tRis.find(ent)!=mif->second.tRis.end()) {
      tGlob.resize(3);
      tGlob[0] = mif->second.data.data.value3;
      tGlob[1] = mif->second.data.data.value4; 
      tGlob[2] = mif->second.data.data.value5;
      tGlob *= mif->second.data.data.value1;
      tGlobNodal.resize(3,3);
      for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  tGlobNodal(nn,dd) = tGlob[dd];
	}
      }
      ierr = reBaseToFaceLoocalCoordSystem(tGlobNodal); CHKERRQ(ierr);
      tLocNodal += tGlobNodal;
    }
  }

  //cerr << tLocNodal << endl;
  t_loc = &*tLocNodal.data().begin();
  //cerr << "tLocNodal: " << tLocNodal << endl;

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::operator()() {
  PetscFunctionBegin;

  //cerr << "MyTriangleSpatialFE::operator()()\n";

  try {

  try {

  dofs_X = &*coords.data().begin();
  for(int ee = 0;ee<3;ee++) {
    dofs_X_edge[ee] = NULL;
    idofs_X_edge[ee] = NULL;
  }
  dofs_X_face = NULL;
  idofs_X_face = NULL;

  dofs_x =  &*coords.data().begin();
  idofs_x = NULL;
  for(int ee = 0;ee<3;ee++) {
    order_edge[ee] = 0;
    N_edge[ee] = NULL;
    diffN_edge[ee] = NULL;
    dofs_x_edge[ee] = NULL;
    idofs_x_edge[ee] = NULL;
  }
  order_face = 0;
  N_face = NULL;
  diffN_face = NULL;
  dofs_x_face = NULL;
  idofs_x_face = NULL;

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  ierr = TriElementForcesAndSurcesCore::operator()(); CHKERRQ(ierr);
  ierr = calcTraction(); CHKERRQ(ierr);

  switch(snes_ctx) {
    case CTX_SNESNONE:
    case CTX_SNESSETFUNCTION: {
      tLocNodal *= *sCaleRhs;
      //cerr << "sCaleRhs " << *sCaleRhs << endl;
      //cerr << tLocNodal << endl;
      ierr = rHs(); CHKERRQ(ierr);
    }
    break;
    case CTX_SNESSETJACOBIAN: {
      tLocNodal *= *sCaleLhs;
      ierr = lHs(); CHKERRQ(ierr);
    }
    break;
  }

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::addForce(int ms_id) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_Cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_cubit_bc_data_structure(mapForce[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::addPreassure(int ms_id) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  const CubitMeshSets *cubit_meshset_ptr;
  ierr = mField.get_Cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
  ierr = cubit_meshset_ptr->get_cubit_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
  rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleMaterialFE::rHs() {
  PetscFunctionBegin;

  //cerr << "MyTriangleMaterialFE::rHs()\n";

  try {
    
  if(dOfs_X_indices.size()!=9) {
    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
  }
  /*cerr << "dOfs_X_indices: " << dOfs_X_indices << endl;
  cerr << "dOfs_X_indices: " << dOfs_X << endl;
  cerr << "order_face_material " << order_face_material << endl;
  cerr << "order_edge_material " 
    << order_edge_material[0] << " " 
    << order_edge_material[1] << " "
    << order_edge_material[2] << endl;*/

  fExtNode.resize(9);	
  ierr = Fext_H_hierarchical(
    order_face_material,order_edge_material,//2
    N,N_face,N_edge,diffN,diffN_face,diffN_edge,//8
    t_loc,NULL,NULL,//11
    dofs_X,dofs_X_edge,dofs_X_face,//14
    NULL,
    &*fExtNode.begin(),NULL,
    gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);

  Vec f = snes_f;
  if(uSeF) f = F; 
  ierr = VecSetValues(f,
    9,dofs_X_indices,
    &*fExtNode.data().begin(),ADD_VALUES); CHKERRQ(ierr);

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode NeummanForcesSurfaceComplexForLazy::MyTriangleMaterialFE::lHs() {
  PetscFunctionBegin;

  try {

  double center[3];
  tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
  cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
  double r = cblas_dnrm2(3,center,1);

  kExtNodeNode.resize(9,9);

  ierr = KExt_HH_hierarchical(
    r*eps,order_face_material,order_edge_material,
    N,N_face,N_edge,diffN,diffN_face,diffN_edge,
    t_loc,NULL,NULL,
    dofs_X,dofs_X_edge,dofs_X_face,
    &*kExtNodeNode.data().begin(),
    gaussPts.size2(),&gaussPts(2,0)); CHKERRQ(ierr);
  ierr = MatSetValues(snes_B,
    9,dofs_X_indices,
    9,dofs_X_indices,
    &*kExtNodeNode.data().begin(),ADD_VALUES); CHKERRQ(ierr);

  } catch (exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }

  PetscFunctionReturn(0);
}


}



