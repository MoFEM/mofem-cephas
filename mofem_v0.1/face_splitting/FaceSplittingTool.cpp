/* Copyright (C) 2014, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
 *
 * Implemented with MOLOKO in headphones https://www.youtube.com/watch?v=46IIAylUph0.
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

#include<FaceSplittingTool.hpp>
#include<FEM.h>
#include<complex_for_lazy.h>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace MoFEM {

PetscErrorCode FaceSplittingTools::buildKDTreeForCrackSurface(
    Range &entities,const BitRefLevel bit_mesh) {
  PetscFunctionBegin;

  if(bit_mesh.any()) {
    Range mesh_level_tris;
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr); 
    entities = intersect(entities,mesh_level_tris);
  }

  rval = moab_distance_from_crack_surface.delete_mesh(); CHKERR_PETSC(rval);

  Tag th_freez;
  rval = mField.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);

  Tag th_normal;
  double def_VAL[4] = { 0,0,0,0 };
  rval = moab_distance_from_crack_surface.tag_get_handle(
    "NORMAL",4,MB_TYPE_DOUBLE,th_normal,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

  //material prositions
  const DofMoFEMEntity_multiIndex *dofs_moabfield_ptr;
  ierr = mField.get_dofs(&dofs_moabfield_ptr); CHKERRQ(ierr);
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dof_iterator;

  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  //crack_front_edges = intersect(crack_front_edges,mesh_level_edges);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);

  Range nodes;
  rval = mField.get_moab().get_connectivity(entities,nodes,true); CHKERR_PETSC(rval);
  for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
    dof_iterator dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));
    dof_iterator hi_dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));	
    if(distance(dit,hi_dit)!=3) {
      SETERRQ(PETSC_COMM_SELF,1,"should three coordinates");
    }
    //get coords from material field
    double coords[3];
    for(;dit!=hi_dit;dit++) {
      coords[dit->get_dof_rank()] = dit->get_FieldData();
    }
    //create coords in moab for kdTree 
    rval = moab_distance_from_crack_surface.create_vertex(coords,map_nodes[*nit]); CHKERR_PETSC(rval);

    int freez;
    mField.get_moab().tag_get_data(th_freez,&*nit,1,&freez);

    if(freez == 0 && crack_front_edges_nodes.find(*nit)!=crack_front_edges_nodes.end() ) {

      dit = 
	dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("GRIFFITH_FORCE_TANGENT",*nit));
      hi_dit = 
	dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("GRIFFITH_FORCE_TANGENT",*nit));	
      if(distance(dit,hi_dit)!=3) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }

      double tangent[3] = { 0,0,0 };
      for(;dit!=hi_dit;dit++) {
	tangent[dit->get_dof_rank()] = dit->get_FieldData();
      }
      //double tangent_nrm2 =  cblas_dnrm2(3,tangent,1);
      //cblas_dscal(3,1./tangent_nrm2,tangent,1);
      //cerr << "tangent_nrm2 " << tangent_nrm2 << endl;

      dit = 
	dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MATERIAL_FORCE",*nit));
      hi_dit = 
	dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MATERIAL_FORCE",*nit));	
      if(distance(dit,hi_dit)!=3) {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }

      double force[3] = { 0,0,0 };
      for(;dit!=hi_dit;dit++) {
	force[dit->get_dof_rank()] = dit->get_FieldData();
      }
      //double force_nrm2 = cblas_dnrm2(3,force,1);
      //cblas_dscal(3,1./force_nrm2,force,1);
      //cerr << "force_nrm2 " << force_nrm2 << endl;

      double Spin_tangent[9] = { 0,0,0, 0,0,0, 0,0,0 };
      ierr = Spin(Spin_tangent,tangent); CHKERRQ(ierr);
      double normal[4] = { 0,0,0,1 };
      cblas_dgemv(CblasRowMajor,CblasNoTrans,3,3,-1,Spin_tangent,3,force,1,0,normal,1); 
      double nrm2 = cblas_dnrm2(3,normal,1);
      cblas_dscal(3,1./nrm2,normal,1);
      //cerr << "norm: " << normal[0] << " " << normal[1] << " " << normal[2] << endl;
    
      rval = moab_distance_from_crack_surface.tag_set_data(th_normal,&map_nodes[*nit],1,normal); CHKERR_PETSC(rval);

    }

  }

  double diffN[3*2];
  ierr = ShapeDiffMBTRI(diffN); CHKERRQ(ierr);

  Range new_entities;
  for(Range::iterator eit = entities.begin();eit!=entities.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    for(int nn = 0;nn<num_nodes;nn++) {
      new_conn[nn] = map_nodes[conn[nn]];
    }
    //create crack surface triangles for kdTree
    EntityHandle new_elem;
    rval = moab_distance_from_crack_surface.create_element(
      mField.get_moab().type_from_handle(*eit),new_conn,num_nodes,new_elem); CHKERR_PETSC(rval);
    new_entities.insert(new_elem);
    rval = moab_distance_from_crack_surface.get_connectivity(new_elem,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes]; 
    rval = moab_distance_from_crack_surface.get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    double normal[4] = { 0,0,0,1 };
    ierr = ShapeFaceNormalMBTRI(diffN,coords,normal); CHKERRQ(ierr);
    double nrm2 = cblas_dnrm2(3,normal,1);
    cblas_dscal(3,1./nrm2,normal,1);
    //calulate cracl surface normal
    rval = moab_distance_from_crack_surface.tag_set_data(th_normal,&new_elem,1,normal); CHKERR_PETSC(rval);
  }

  rval = moab_distance_from_crack_surface.create_meshset(MESHSET_SET,kdTree_rootMeshset_DistanceFromCrackSurface); CHKERR_PETSC(rval);
  rval = kdTree_DistanceFromCrackSurface.build_tree(new_entities,kdTree_rootMeshset_DistanceFromCrackSurface); CHKERR_PETSC(rval);

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::buildKDTreeForCrackSurface(const BitRefLevel bit_mesh) {
  PetscFunctionBegin;
  Range crack_surface;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface,true); CHKERRQ(ierr);
  if(bit_mesh.any()) {
    Range mesh_level_tris;
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr); 
    crack_surface = intersect(crack_surface,mesh_level_tris);
  }
  Tag th_interface_side;
  rval = mField.get_moab().tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
  Range entities;
  for(Range::iterator cit = crack_surface.begin();cit!=crack_surface.end();cit++) {
    int side;
    rval = mField.get_moab().tag_get_data(th_interface_side,&*cit,1,&side); CHKERR_PETSC(rval);
    if(side == 1) {
      entities.insert(*cit);
    }
  }
  ierr = buildKDTreeForCrackSurface(entities,0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::initBitLevelData(BitRefLevel bit_mesh) {
  PetscFunctionBegin;

  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr); 
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculateDistanceFromCrackSurface(Range &nodes,double alpha) {
  PetscFunctionBegin;	

  nodes = intersect(nodes,mesh_level_nodes);

  double def_VAL[3] = { 0,0,0 };
  rval = mField.get_moab().tag_get_handle(
    "DISTANCE_FROM_CRACK_SURFACE",1,MB_TYPE_DOUBLE,th_distance,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);
  rval = mField.get_moab().tag_get_handle(
    "PROJECTION_CRACK_SURFACE",3,MB_TYPE_DOUBLE,th_projection,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

  Tag th_normal;
  rval = moab_distance_from_crack_surface.tag_get_handle("NORMAL",th_normal); CHKERR_THROW(rval);

  //material prositions
  const DofMoFEMEntity_multiIndex *dofs_moabfield_ptr;
  ierr = mField.get_dofs(&dofs_moabfield_ptr); CHKERRQ(ierr);
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dof_iterator;


  for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
    double coords[3]; 
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    dof_iterator dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));
    dof_iterator hi_dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));	
    if(distance(dit,hi_dit)==3) {
      //get coords from material field
      for(;dit!=hi_dit;dit++) {
	coords[dit->get_dof_rank()] = dit->get_FieldData();
      }
    }

    double closest_point_out[3];
    EntityHandle triangle_out;
    rval = kdTree_DistanceFromCrackSurface.closest_triangle(
      kdTree_rootMeshset_DistanceFromCrackSurface,coords,closest_point_out,triangle_out); CHKERR_PETSC(rval);
    int num_nodes; 
    const EntityHandle* conn;
    rval = moab_distance_from_crack_surface.get_connectivity(triangle_out,conn,num_nodes,true); CHKERR_PETSC(rval);
    double face_coords[num_nodes*3];
    rval = moab_distance_from_crack_surface.get_coords(conn,num_nodes,face_coords); CHKERR_PETSC(rval);
    double normals[num_nodes*4];
    rval = moab_distance_from_crack_surface.tag_get_data(th_normal,conn,num_nodes,normals); CHKERR_PETSC(rval);

    double min_dist = -1;
    double distance[3];
    cblas_dcopy(3,coords,1,distance,1);
    for(int nn = 0;nn<3;nn++) {
      if(normals[nn*4+3]==0) continue;
      cblas_daxpy(3,-1,&face_coords[3*nn],1,distance,1);
      double dist0 = cblas_dnrm2(3,distance,1);
      if(min_dist<0) {
	min_dist = dist0;
      } else {
	min_dist = fmin(dist0,min_dist);
      }
      if(dist0 == min_dist) {
	double dot = -cblas_ddot(3,&normals[nn*4],1,distance,1);
	double projection[3];
	cblas_dcopy(3,distance,1,projection,1);
	cblas_daxpy(3,alpha*dot,&normals[nn*4],1,projection,1);
	double dist1 = cblas_dnrm2(3,projection,1);
	if(dist0>0) {
	  cblas_dscal(3,dist1/dist0,projection,1);
	}
	cblas_daxpy(3,1.,&face_coords[3*nn],1,projection,1);
	//cerr << "a dist1/dist0 " << dist0 << " " << dist1  << " " << dist1/dist0 << " " << dot << " ";
	//cerr << projection[0] << " " << projection[1] << " " << projection[2] << endl;
	rval = mField.get_moab().tag_set_data(th_distance,&*nit,1,&dot); CHKERR_PETSC(rval);
	rval = mField.get_moab().tag_set_data(th_projection,&*nit,1,projection); CHKERR_PETSC(rval);
      }
    }
    if(min_dist<0) {	
      cblas_daxpy(3,-1,closest_point_out,1,distance,1);
      double dist0 = cblas_dnrm2(3,distance,1);
      double normal[4];
      rval = moab_distance_from_crack_surface.tag_get_data(th_normal,&triangle_out,1,normal); CHKERR_PETSC(rval);
      double dot = -cblas_ddot(3,normal,1,distance,1);
      double projection[3];
      cblas_dcopy(3,distance,1,projection,1);
      cblas_daxpy(3,alpha*dot,normal,1,projection,1);
      double dist1 = cblas_dnrm2(3,projection,1);
      if(dist0>0) {
	cblas_dscal(3,dist1/dist0,projection,1);
      }
      cblas_daxpy(3,1.,closest_point_out,1,projection,1);
      //cerr << "b dist1/dist0 " << dist0 << " " << dist1  << " " << dist1/dist0 << " " << dot << " ";
      //cerr << projection[0] << " " << projection[1] << " " << projection[2] << endl;
      rval = mField.get_moab().tag_set_data(th_distance,&*nit,1,&dot); CHKERR_PETSC(rval);
      rval = mField.get_moab().tag_set_data(th_projection,&*nit,1,projection); CHKERR_PETSC(rval);
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculateDistanceCrackFrontNodesFromCrackSurface(double alpha) {
  PetscFunctionBegin;
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_front_edges = intersect(crack_front_edges,mesh_level_edges);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  ierr = calculateDistanceFromCrackSurface(crack_front_edges_nodes,alpha); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculateDistanceFromCrackSurface() {
  PetscFunctionBegin;
  Range nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets,nodes,true); CHKERR_PETSC(rval);
  ierr = calculateDistanceFromCrackSurface(nodes,1.); CHKERRQ(ierr);
  //if(verb>0) {
  /*  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,mesh_level_tets); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("distances.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }*/
  //}
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getOpositeForntEdges(bool createMeshset) {
  PetscFunctionBegin;

  //create meshset
  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,opositeFrontEdges); CHKERR_PETSC(rval);
  }

  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_front_edges = intersect(crack_front_edges,mesh_level_edges);

  //get tets adjacent to crack front
  Range crack_front_adj_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_adj_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_adj_tets;
  rval = mField.get_moab().get_adjacencies(crack_front_adj_nodes,3,false,crack_front_adj_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_adj_tets = intersect(crack_front_adj_tets,mesh_level_tets);

  //get edges
  Range crack_front_adj_tets_edges; 
  rval = mField.get_moab().get_adjacencies(crack_front_adj_tets,1,false,crack_front_adj_tets_edges,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_adj_tets_edges = subtract(crack_front_adj_tets_edges,crack_front_edges);

  /*//crack_front_test->skin faces->skin_edges;
  Skinner skin(&mField.get_moab());
  Range skin_faces; 
  rval = skin.find_skin(crack_front_adj_tets,false,skin_faces); CHKERR(rval);
  Range skin_edges;
  rval = mField.get_moab().get_adjacencies(skin_faces,1,false,skin_edges,Interface::UNION); CHKERR_PETSC(rval);*/

  //crack surface faces->nodes->edges;
  Range crack_surface;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface,true); CHKERRQ(ierr);
  Range crack_surface_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface,crack_surface_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_edges;
  rval = mField.get_moab().get_adjacencies(crack_surface_nodes,1,false,crack_surface_edges,Interface::UNION); CHKERR_PETSC(rval);
  //skin_edges = subtract(skin_edges,crack_surface_edges);
  crack_front_adj_tets_edges = subtract(crack_front_adj_tets_edges,crack_surface_edges);
 
  //for(Range::iterator sit = skin_edges.begin();sit!=skin_edges.end();sit++) {
  Range::iterator sit = crack_front_adj_tets_edges.begin();
  for(;sit!=crack_front_adj_tets_edges.end();sit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*sit,conn,num_nodes,true); CHKERR_PETSC(rval);
    //double coords[3*num_nodes];
    //rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    //cblas_daxpy(3,-1,coords,1,&coords[3],1);
    //double nrm2 = cblas_dnrm2(3,&coords[3],1);	
    double d[2];
    rval = mField.get_moab().tag_get_data(th_distance,conn,num_nodes,d); CHKERR_PETSC(rval);
    if(d[0]*d[1]<=0) {
      rval = mField.get_moab().add_entities(opositeFrontEdges,&*sit,1); CHKERR_PETSC(rval);
    }

  }
 
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getCrackFrontEntities(bool createMeshset,bool get_tets,bool get_faces,int verb) {
  PetscFunctionBegin;

  //get oposite crack foront faces, that is edges->nodes->tets
  Range opposite_crack_front_edges;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,opposite_crack_front_edges,false); CHKERR_PETSC(rval);
  opposite_crack_front_edges = intersect(opposite_crack_front_edges,mesh_level_edges); //only edges on problem mesh level
  Range opposite_crack_front_edges_nodes; 
  rval = mField.get_moab().get_connectivity(opposite_crack_front_edges,opposite_crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  //tets adjacent to opposite edges
  Range opposite_crack_front_edges_tets; 
  rval = mField.get_moab().get_adjacencies(opposite_crack_front_edges,3,false,opposite_crack_front_edges_tets,Interface::UNION); CHKERR_PETSC(rval);
  opposite_crack_front_edges_tets = intersect(opposite_crack_front_edges_tets,mesh_level_tets);

  //get crack front faces, that is edges->nodes->tets
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_front_edges = intersect(crack_front_edges,mesh_level_edges); //only edges on problem level
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_edges_nodes_faces;
  rval = mField.get_moab().get_adjacencies(crack_front_edges_nodes,2,false,crack_front_edges_nodes_faces,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_front_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(crack_front_edges_nodes,3,false,crack_front_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_edges_nodes_tets = intersect(crack_front_edges_nodes_tets,mesh_level_tets);
  Range crack_front_edges_faces;
  rval = mField.get_moab().get_adjacencies(crack_front_edges,2,false,crack_front_edges_faces,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_edges_faces = intersect(crack_front_edges_faces,mesh_level_tris);

  //common tets
  Range common_tets = intersect(crack_front_edges_nodes_tets,opposite_crack_front_edges_tets);

  Range crack_surface;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface,true); CHKERRQ(ierr);
  crack_surface = intersect(crack_surface,mesh_level_tris);
  Range crack_surface_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface,crack_surface_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_nodes_without_front = subtract(crack_surface_nodes,crack_front_edges_nodes);
  Range crack_surface_nodes_without_front_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_nodes_without_front,3,false,crack_surface_nodes_without_front_tets,Interface::UNION); CHKERR_PETSC(rval);
  common_tets = subtract(common_tets,crack_surface_nodes_without_front_tets);

  Range crack_front_edges_nodes_edges_faces;
  if(get_faces) {
    //get edges which are connecting crack front nodes but are not part of crack front itself
    //nodes of such edges, and node of crack front create face, which is part of crack.
    Range crack_front_edges_nodes_edges;
    rval = mField.get_moab().get_adjacencies(crack_front_edges_nodes,1,false,crack_front_edges_nodes_edges,Interface::UNION); CHKERR_PETSC(rval);
    crack_front_edges_nodes_edges = intersect(crack_front_edges_nodes_edges,mesh_level_edges);
    crack_front_edges_nodes_edges = subtract(crack_front_edges_nodes_edges,crack_front_edges);
    Range::iterator eit = crack_front_edges_nodes_edges.begin();
    for(;eit!=crack_front_edges_nodes_edges.end();) {
      Range edge_nodes;
      rval = mField.get_moab().get_connectivity(&*eit,1,edge_nodes,true); CHKERR_PETSC(rval);
      edge_nodes = subtract(edge_nodes,crack_front_edges_nodes);
      if(!edge_nodes.empty()) {
        eit = crack_front_edges_nodes_edges.erase(eit);
      } else {
        eit++;
      }
    }
    //faces of crack front adjacent edges connecting crack front nodes
    rval = mField.get_moab().get_adjacencies(
      crack_front_edges_nodes_edges,2,false,crack_front_edges_nodes_edges_faces,Interface::UNION); CHKERR_PETSC(rval);
    crack_front_edges_nodes_edges_faces = intersect(crack_front_edges_nodes_edges_faces,crack_front_edges_faces);
    crack_front_edges_nodes_edges_faces = subtract(crack_front_edges_nodes_edges_faces,crack_surface);
    //reject faces at corner which angles are smaller than PI
    //edges at corner should be not-conbex 
    Range common_tets_edges;
    rval = mField.get_moab().get_adjacencies(common_tets,1,false,common_tets_edges,Interface::UNION); CHKERR_PETSC(rval);
    Range::iterator fit = crack_front_edges_nodes_edges_faces.begin();
    for(;fit!=crack_front_edges_nodes_edges_faces.end();) {
      Range fit_edges;
      rval = mField.get_moab().get_adjacencies(&*fit,1,1,false,fit_edges); CHKERR_PETSC(rval);
      Range free_fit_edges;
      free_fit_edges = subtract(fit_edges,crack_front_edges);
      fit_edges = intersect(fit_edges,crack_front_edges);
      if(fit_edges.size()<2) {
        fit = crack_front_edges_nodes_edges_faces.erase(fit);
      } else {
        if(fit_edges.size()==2) {
  	if(intersect(free_fit_edges,common_tets_edges).empty()) {
  	  Range fit_tets;
  	  rval = mField.get_moab().get_adjacencies(&*fit,1,3,false,fit_tets); CHKERR_PETSC(rval);
  	  if(intersect(fit_tets,crack_surface_nodes_without_front_tets).size()>0) {
  	    fit = crack_front_edges_nodes_edges_faces.erase(fit); 
  	  } else {
  	    fit++;
  	  }
  	} else {
  	  fit++;
  	}
        }
      }
    }
  }

  if(get_tets) {
  
    //get faces
    Skinner skin(&mField.get_moab());
    //
    Range common_tets_faces;
    rval = mField.get_moab().get_adjacencies(common_tets,2,false,common_tets_faces,Interface::UNION); CHKERR_PETSC(rval);
    //remove faces on body surface
    Range outer_surface_skin;
    rval = skin.find_skin(mesh_level_tets,false,outer_surface_skin); CHKERR(rval);
    common_tets_faces = subtract(common_tets_faces,outer_surface_skin);
    //remove faces not adjacent to crack front
    common_tets_faces = intersect(common_tets_faces,crack_front_edges_nodes_faces);

    if(createMeshset) {
      rval = mField.get_moab().create_meshset(MESHSET_SET,crackFrontTests); CHKERR_PETSC(rval);
    }
    ierr = mField.get_moab().add_entities(crackFrontTests,common_tets); CHKERRQ(ierr);
    ierr = mField.get_moab().add_entities(crackFrontTests,common_tets_faces); CHKERRQ(ierr);

    //add nodes at new crack front
    Range common_tets_faces_nodes;
    rval = mField.get_moab().get_connectivity(common_tets_faces,common_tets_faces_nodes); CHKERR_PETSC(rval);
    common_tets_faces_nodes = subtract(common_tets_faces_nodes,crack_front_edges_nodes);
    ierr = mField.get_moab().add_entities(crackFrontTests,common_tets_faces_nodes); CHKERRQ(ierr);

    if(verb>0) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {
	EntityHandle out_meshset;
	rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,common_tets); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,common_tets_faces); CHKERR_PETSC(rval);
	//rval = mField.get_moab().add_entities(out_meshset,common_tets_edges); CHKERR_PETSC(rval);
	rval = mField.get_moab().write_file("crack_front_tets.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }
    }

  }

  if(get_faces) {
    ierr = mField.get_moab().add_entities(selectedCrackFaces,crack_front_edges_nodes_edges_faces); CHKERRQ(ierr);
    if(verb>0) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {
	EntityHandle out_meshset;
	rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = mField.get_moab().add_entities(out_meshset,crack_front_edges_nodes_edges_faces); CHKERR_PETSC(rval);
	rval = mField.get_moab().write_file("crack_front_faces.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getCrackFrontFaces(bool createMeshset,int verb) {
  PetscFunctionBegin;
  ierr = getCrackFrontEntities(createMeshset,false,true,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getCrackFrontTets(bool createMeshset,int verb) {
  PetscFunctionBegin;
  ierr = getCrackFrontEntities(createMeshset,true,false,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(bool createMeshset,int verb) {
  PetscFunctionBegin;

  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,chopTetsFaces); CHKERR_PETSC(rval);
  }

  //crack fornt tets
  Range crack_front_tets;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBTET,crack_front_tets,true); CHKERR_PETSC(rval);
  Range crack_front_tets_nodes;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBVERTEX,crack_front_tets_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_tets_faces;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBTRI,crack_front_tets_faces,true); CHKERR_PETSC(rval);

  Range crack_front_tets0 = crack_front_tets;
  Range crack_front_tets_nodes0 = crack_front_tets_nodes;
  Range crack_front_tets_faces0 = crack_front_tets_faces;

  //crack fornt edges
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  //front nodes
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  //nodes faces
  Range crack_front_edges_nodes_faces;
  rval = mField.get_moab().get_adjacencies(
    crack_front_edges_nodes,2,false,crack_front_edges_nodes_faces,Interface::UNION); CHKERR_PETSC(rval);
  //crack front edges nodes tets
  Range crack_front_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_front_edges_nodes,3,false,crack_front_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_edges_nodes_tets = intersect(crack_front_edges_nodes_tets,mesh_level_tets);
  //crack edges adjacent to crack front edges nodes
  Range crack_front_edges_nodes_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_front_edges_nodes,1,false,crack_front_edges_nodes_edges,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_edges_nodes_edges = intersect(crack_front_edges_nodes_edges,mesh_level_edges);

  //crack surface faces edges tets
  Range crack_surface;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface,true); CHKERRQ(ierr);
  Range crack_surface_edges;
  rval = mField.get_moab().get_adjacencies(crack_surface,1,false,crack_surface_edges,Interface::UNION); CHKERR_PETSC(rval);
  crack_surface_edges = subtract(crack_surface_edges,crack_front_edges);
  Range crack_surface_edges_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_edges,3,false,crack_surface_edges_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_surface_edges_tets = intersect(crack_surface_edges_tets,mesh_level_tets);

  if(verb>100) {
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_tets.size() %d\n",crack_front_tets.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_tets_nodes.size() %d\n",crack_front_tets_nodes.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_tets_faces.size() %d\n",crack_front_tets_faces.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_edges.size() %d\n",crack_front_edges.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_edges_nodes.size() %d\n",crack_front_edges_nodes.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_edges_nodes_tets.size() %d\n",crack_front_edges_nodes_tets.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_edges_nodes_faces.size() %d\n",crack_front_edges_nodes_faces.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_front_edges_nodes_edges.size() %d\n",crack_front_edges_nodes_edges.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_surface.size() %d\n",crack_surface.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_surface_edges.size() %d\n",crack_surface_edges.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"crack_surface_edges_tets.size() %d\n",crack_surface_edges_tets.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"mesh_level_tets.size() %d\n",mesh_level_tets.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"mesh_level_nodes.size() %d\n",mesh_level_nodes.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"mesh_level_edges.size() %d\n",mesh_level_edges.size());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"mesh_level_tris.size() %d\n",mesh_level_tris.size());
    PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  }

  Range _crack_front_tets_faces_;
  Range _crack_front_tets_edges_;
  Range _crack_front_free_nodes_;
  Range _crack_front_free_faces_;
  Range _crack_front_free_faces_edges_;
  Range _crack_front_free_faces_nodes_;
  Range _crack_front_body_skin_faces_;
  Range _crack_front_body_skin_faces_edges_;
  Range _crack_front_body_skin_edges_;
  Range _crack_front_body_skin_edges_nodes_;
  Range _crack_front_tets_skin_faces_;
  Range _other_crack_front_tets_faces_;
  Range _other_crack_front_tets_faces_edges_;
  Range _nodes_on_skin_surface_;

  int debug_ii = 0;

  if(verb>=0) {
    PetscPrintf(PETSC_COMM_WORLD,"Tets to chop: ");
  }

  Skinner skin(&mField.get_moab());
  Range mesh_level_tets_skin_faces;
  rval = skin.find_skin(mesh_level_tets,false,mesh_level_tets_skin_faces); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_faces_edges;
  rval = mField.get_moab().get_adjacencies(
    mesh_level_tets_skin_faces,1,false,mesh_level_tets_skin_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_faces_nodes;
  rval = mField.get_moab().get_connectivity(
    mesh_level_tets_skin_faces,mesh_level_tets_skin_faces_nodes,true); CHKERR_PETSC(rval);

  int ii = 0;

  do {

    //get faces adjacent to crack front tets
    _crack_front_tets_faces_.clear();
    rval = mField.get_moab().get_adjacencies(
      crack_front_tets,2,false,_crack_front_tets_faces_,Interface::UNION); CHKERR_PETSC(rval);

    //get free faces and their nodes
    //those nodes are excluded from chopping procedure
    //free faces are face of crack front faces minus faces of crack front tets
    _crack_front_free_faces_ = subtract(crack_front_tets_faces,_crack_front_tets_faces_);
    _crack_front_free_faces_nodes_.clear();
    rval = mField.get_moab().get_connectivity(_crack_front_free_faces_,_crack_front_free_faces_nodes_,true); CHKERR_PETSC(rval);

    //get faces on body surface adjacent to crack front tets
    //crack front faces on body surface is inetesection of skin of crack front tetst and body tets
    _crack_front_body_skin_faces_ = intersect(_crack_front_tets_faces_,mesh_level_tets_skin_faces); //faces on body skin
    //edges on body skin and adjacent to faces of crack front tets which are on body surface
    _crack_front_body_skin_faces_edges_.clear();
    rval = mField.get_moab().get_adjacencies(
      _crack_front_body_skin_faces_,1,false,_crack_front_body_skin_faces_edges_,Interface::UNION); CHKERR_PETSC(rval); 

    //Take skin of front tets. 
    _crack_front_tets_skin_faces_.clear();
    rval = skin.find_skin(crack_front_tets,false,_crack_front_tets_skin_faces_); CHKERR_PETSC(rval);
    _crack_front_tets_skin_faces_ = intersect(_crack_front_tets_skin_faces_,crack_front_edges_nodes_faces);
    _other_crack_front_tets_faces_ = subtract(_crack_front_tets_faces_,_crack_front_tets_skin_faces_);
    rval = mField.get_moab().get_adjacencies(
      _other_crack_front_tets_faces_,1,false,_other_crack_front_tets_faces_edges_,Interface::UNION); CHKERR_PETSC(rval); 

    //body surface edges
    _crack_front_tets_edges_.clear();
    rval = mField.get_moab().get_adjacencies(
      crack_front_tets,1,false,_crack_front_tets_edges_,Interface::UNION); CHKERR_PETSC(rval);
    //edges on body skin and adjacent to crack front tets
    //note that such edges are bigger set than _crack_front_body_skin_faces_edges_
    _crack_front_body_skin_edges_ = intersect(_crack_front_tets_edges_,mesh_level_tets_skin_faces_edges);
    //get the edges on body skin which are not adges on any face on body skin
    //removing node adjacent to such edge will crate gap in extenended crack surface
    _crack_front_body_skin_edges_ = subtract(_crack_front_body_skin_edges_,_other_crack_front_tets_faces_edges_);
    _crack_front_body_skin_edges_ = subtract(_crack_front_body_skin_edges_,_crack_front_body_skin_faces_edges_);
    //get nodes on body skin which can not be removed
    _crack_front_body_skin_edges_nodes_.clear();
    for(Range::iterator eit = _crack_front_body_skin_edges_.begin();
      eit!=_crack_front_body_skin_edges_.end();eit++) {
      Range adj_tets;
      rval = mField.get_moab().get_adjacencies(&*eit,1,3,false,adj_tets); CHKERR_PETSC(rval);
      adj_tets = intersect(adj_tets,crack_front_tets);
      if(adj_tets.size() == 1) {
	Range eit_conn;
	rval = mField.get_moab().get_connectivity(&*eit,1,eit_conn,true); CHKERR_PETSC(rval);
	_crack_front_body_skin_edges_nodes_.merge(eit_conn);
      }
    }

    //Get adjacent nodes to skin faces. Subtract from those nodes
    //nodes on crack front and nodes of free edges.
    _nodes_on_skin_surface_.clear();
    rval = mField.get_moab().get_connectivity(_crack_front_tets_skin_faces_,_nodes_on_skin_surface_,true); CHKERR_PETSC(rval);
    _nodes_on_skin_surface_ = intersect(_nodes_on_skin_surface_,crack_front_tets_nodes); 
    _nodes_on_skin_surface_ = subtract(_nodes_on_skin_surface_,crack_front_edges_nodes);
    _nodes_on_skin_surface_ = subtract(_nodes_on_skin_surface_,_crack_front_free_faces_nodes_);
    _nodes_on_skin_surface_ = subtract(_nodes_on_skin_surface_,_crack_front_body_skin_edges_nodes_);

    Range _nodes_on_skin_surface_edges_;
    rval = mField.get_moab().get_adjacencies(
      _nodes_on_skin_surface_,1,false,_nodes_on_skin_surface_edges_,Interface::UNION); CHKERR_PETSC(rval);
    _nodes_on_skin_surface_edges_ = intersect(_nodes_on_skin_surface_edges_,crack_front_edges_nodes_edges);
    _nodes_on_skin_surface_edges_ = subtract(_nodes_on_skin_surface_edges_,mesh_level_tets_skin_faces_edges);
    _nodes_on_skin_surface_edges_ = subtract(_nodes_on_skin_surface_edges_,crack_front_edges);

    for(Range::iterator eit = _nodes_on_skin_surface_edges_.begin();
	eit != _nodes_on_skin_surface_edges_.end(); eit++) {
    
      Range eit_tets;
      rval = mField.get_moab().get_adjacencies(&*eit,1,3,false,eit_tets); CHKERR_PETSC(rval);
      eit_tets = intersect(eit_tets,crack_front_tets);
      Range eit_tets_skin;
      rval = skin.find_skin(eit_tets,false,eit_tets_skin); CHKERR_PETSC(rval);
      Range eit_tets_skin_faces_edges;
      rval = mField.get_moab().get_adjacencies(
	eit_tets_skin,1,false,eit_tets_skin_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
      if(eit_tets_skin_faces_edges.find(*eit)!=eit_tets_skin_faces_edges.end()) {
	Range eit_faces;
	rval = mField.get_moab().get_adjacencies(&*eit,1,2,false,eit_faces); CHKERR_PETSC(rval);
	if(intersect(eit_faces,eit_tets_skin).size()>2) {
	  Range eit_nodes;
	  rval = mField.get_moab().get_connectivity(&*eit,1,eit_nodes,true); CHKERR_PETSC(rval);
	  _nodes_on_skin_surface_  = subtract(_nodes_on_skin_surface_,eit_nodes);
	}
      }

    }

    if(intersect(_nodes_on_skin_surface_,crack_front_edges_nodes).size()>0) {
      SETERRQ(PETSC_COMM_SELF,1,"should not happen");
    }

    if(verb>100) {
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n",_nodes_on_skin_surface_.size());
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"_nodes_on_skin_surface_.size() %d\n",_nodes_on_skin_surface_.size());
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"_crack_front_free_faces_nodes_.size() %d\n",_crack_front_free_faces_nodes_.size());
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"_crack_front_body_skin_edges_nodes_.size() %d\n",_crack_front_body_skin_edges_nodes_.size());
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
    }

    //if(_nodes_on_skin_surface_.empty()) {
      //PetscAttachDebugger();
    //}

    //check for tets which have all nodes on crack front old and/or new
    //such tets projected on crack surface can create zero volume if crack surface is plane,
    //in some case negative quality
    //*****************************
    if(_nodes_on_skin_surface_.empty()) {

        if(crack_front_tets.size()>0) {
      
	if(verb>=0) {
	  PetscPrintf(PETSC_COMM_WORLD," Zero Nodes");
	}

	_crack_front_free_faces_edges_.clear();
	rval = mField.get_moab().get_adjacencies(
	  _crack_front_free_faces_,1,false,_crack_front_free_faces_edges_,Interface::UNION); CHKERR_PETSC(rval);

	Range chop_tets;
	Range::iterator tit = crack_front_tets.begin();
	for(;tit!=crack_front_tets.end();tit++) {
	  Range tit_edges;
	  rval = mField.get_moab().get_adjacencies(&*tit,1,1,false,tit_edges); CHKERR_PETSC(rval);
	  tit_edges = subtract(tit_edges,crack_front_edges);
	  tit_edges = subtract(tit_edges,_crack_front_free_faces_edges_);
	  tit_edges = subtract(tit_edges,_crack_front_body_skin_edges_);
	  Range edge_to_remove;
	  Range::iterator eit = tit_edges.begin();
	  //select edge ajacent to one tet
	  for(;eit!=tit_edges.end();eit++) {
	    Range eit_tets;
	    rval = mField.get_moab().get_adjacencies(
	      &*eit,1,3,false,eit_tets); CHKERR_PETSC(rval);
	    eit_tets = intersect(eit_tets,crack_front_tets);
	    if(eit_tets.size()==1) {
	      edge_to_remove.insert(*eit);
	    }
	  }
	  if(edge_to_remove.empty()) {
	    SETERRQ(PETSC_COMM_SELF,1,"how this happen?");
	  }
	  eit = edge_to_remove.begin();
	  for(;eit!=edge_to_remove.end();eit++) {	
	    Range eit_nodes;
	    rval = mField.get_moab().get_connectivity(&*eit,1,eit_nodes,true); CHKERR_PETSC(rval);
	    //skip edges connecting nodes of crack front
	    if(intersect(eit_nodes,crack_front_edges_nodes).size()==0) {
	      continue;
	    }
	    Range eit_faces;
	    rval = mField.get_moab().get_adjacencies(
	      &*eit,1,2,false,eit_faces); CHKERR_PETSC(rval);
	    crack_front_tets_faces = subtract(crack_front_tets_faces,eit_faces);
	    Range eit_tets;
  	    rval = mField.get_moab().get_adjacencies(
	      &*eit,1,3,false,eit_tets); CHKERR_PETSC(rval);
	    eit_tets = intersect(eit_tets,crack_front_tets);
	    chop_tets.merge(eit_tets);
	    break;
	  }
	  if(eit == edge_to_remove.end()) {
	    Range tit_nodes;
	    rval = mField.get_moab().get_connectivity(&*tit,1,tit_nodes,true); CHKERR_PETSC(rval);
	    tit_nodes = intersect(tit_nodes,crack_front_edges_nodes);
	    Tag th_freez;
	    rval = mField.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
	    vector<int> freezed_nodes(tit_nodes.size());
	    rval = mField.get_moab().tag_get_data(th_freez,tit_nodes,&freezed_nodes[0]); CHKERR_PETSC(rval);
	    vector<int>::iterator vit;
	    for(vit = freezed_nodes.begin();vit != freezed_nodes.end(); vit++) {
	      if(*vit == 0) {
		SETERRQ(PETSC_COMM_SELF,1,"how this happen?");
	      }
	    }
	    chop_tets.insert(*tit);
	  }
	}
	crack_front_tets = subtract(crack_front_tets,chop_tets);
    
	debug_ii++;

	if(verb>1) {
	  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
	  if(pcomm->rank()==0) {
	    EntityHandle out_meshset;
	    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	    rval = mField.get_moab().add_entities(out_meshset,crack_front_tets); CHKERR_PETSC(rval);
	    rval = mField.get_moab().add_entities(out_meshset,crack_front_tets_faces); CHKERR_PETSC(rval);
	    ostringstream ss;
	    ss << "chop_debug_" << debug_ii << ".vtk";
	    rval = mField.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	  }
	}

	if(crack_front_tets.empty()) {
	  break;
	} else {
	  continue;
	}
      } else {
	SETERRQ(PETSC_COMM_SELF,1,"imposible case");
      }

    }
    //*****************************

    if(_nodes_on_skin_surface_.empty()) {
      SETERRQ(PETSC_COMM_SELF,1,"it should at least be one node");
    }
  
    double min_nit_q;
    EntityHandle min_node = 0;
    Range::iterator nit;
    nit = _nodes_on_skin_surface_.begin();
    for(;nit!=_nodes_on_skin_surface_.end();nit++) {
      Range check_nit;
      check_nit.insert(*nit);
      check_nit.merge(_crack_front_free_faces_nodes_);
      Range _chop_tets_;
      rval = mField.get_moab().get_adjacencies(&*nit,1,3,false,_chop_tets_); CHKERR_PETSC(rval);
      _chop_tets_ = intersect(_chop_tets_,crack_front_tets);
      double q_nit;
      ierr = calculate_qualityAfterProjectingNodes(
	check_nit,_chop_tets_,
	crack_front_edges_nodes,q_nit); CHKERRQ(ierr);
      if(nit == _nodes_on_skin_surface_.begin()) {
	min_nit_q = q_nit;
	min_node = *nit;
      } else {
	if(min_nit_q > q_nit) {
	  min_nit_q = q_nit;
	  min_node = *nit;
	}
      }
    }
    if(min_node == 0) {
      SETERRQ(PETSC_COMM_SELF,1,"no node to erase");
    }
    crack_front_tets_nodes.erase(min_node);

    //check if node connecting two or more parts of tets without common face only connected by edge
    //removing such node will create crack front node without all new front faces connected to it
    Range chop_tets;
    rval = mField.get_moab().get_adjacencies(&min_node,1,3,false,chop_tets); CHKERR_PETSC(rval);
    chop_tets = intersect(chop_tets,crack_front_tets);
    if(chop_tets.empty()) {
      SETERRQ(PETSC_COMM_SELF,1,"it is empty, algorithm is stack");
    }

    //remove choped tets
    //unsigned int nb_crack_front_tets = crack_front_tets.size();
    crack_front_tets = subtract(crack_front_tets,chop_tets);

    //get faces adjacent to removed node
    Range chop_faces;
    rval = mField.get_moab().get_adjacencies(&min_node,1,2,false,chop_faces); CHKERR_PETSC(rval);
    //remove all the faces adjacent to removed node, those faces not create crack surface
    chop_faces = intersect(crack_front_tets_faces,chop_faces);
    if(chop_faces.empty()) {
      SETERRQ(PETSC_COMM_SELF,1,"it is empty, algorithm is stack");
    }
    unsigned int nb_crack_front_tets_faces = crack_front_tets_faces.size();
    crack_front_tets_faces = subtract(crack_front_tets_faces,chop_faces);

    //remove faces with dangling nodes
    {
 
      Range::iterator fit = crack_front_tets_faces.begin();
      for(;fit!=crack_front_tets_faces.end();) {
	Range fit_nodes;
	rval = mField.get_moab().get_connectivity(&*fit,1,fit_nodes,true); CHKERR_PETSC(rval);
	fit_nodes = subtract(fit_nodes,crack_front_edges_nodes);
	fit_nodes = subtract(fit_nodes,mesh_level_tets_skin_faces_nodes);
	Range crack_front_tets_faces_but_one = crack_front_tets_faces;
	crack_front_tets_faces_but_one.erase(*fit);
	Range crack_front_tets_faces_but_one_nodes;
	rval = mField.get_moab().get_connectivity(
	  crack_front_tets_faces_but_one,crack_front_tets_faces_but_one_nodes,true); CHKERR_PETSC(rval);
	Range::iterator nit = fit_nodes.begin();
	for(;nit!=fit_nodes.end();nit++) {
	  if(crack_front_tets_faces_but_one_nodes.find(*nit) 
	      == crack_front_tets_faces_but_one_nodes.end()) {
	    crack_front_tets_faces.erase(*fit);
	    break;
	  }
	}
	if(nit!=fit_nodes.end()) {
	  fit = crack_front_tets_faces.begin();
	} else {
	  fit++;
	}
      }

    }

    //check if face at edge create T-connection
    {

      Range crack_front_tets_edges;
      rval = mField.get_moab().get_adjacencies(
	crack_front_tets,1,false,crack_front_tets_edges,Interface::UNION); CHKERR_PETSC(rval);
      Range crack_front_tets_faces_edges;
      rval = mField.get_moab().get_adjacencies(
	crack_front_tets_faces,1,false,crack_front_tets_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
      crack_front_tets_faces_edges = intersect(crack_front_tets_faces_edges,crack_front_edges_nodes_edges);
      crack_front_tets_faces_edges = subtract(crack_front_tets_faces_edges,crack_front_edges);
      crack_front_tets_faces_edges = subtract(crack_front_tets_faces_edges,crack_front_tets_edges);

      for(Range::iterator eit = crack_front_tets_faces_edges.begin();
	eit!=crack_front_tets_faces_edges.end();eit++) {
	Range eit_faces;
	rval = mField.get_moab().get_adjacencies(&*eit,1,2,false,eit_faces); CHKERR_PETSC(rval);
    	eit_faces = intersect(eit_faces,crack_front_tets_faces);
	unsigned int T_test = 2;
	if(mesh_level_tets_skin_faces_edges.find(*eit)!=mesh_level_tets_skin_faces_edges.end()) {
	  T_test = 1;
	} 
	if(eit_faces.size()>T_test) {
	  Range::iterator fit = eit_faces.begin();
	  for(;fit!=eit_faces.end();fit++) {
	    Range crack_front_tets_faces_but_one;
	    crack_front_tets_faces_but_one = crack_front_tets_faces;
	    crack_front_tets_faces_but_one.erase(*fit);
	    Range crack_front_tets_faces_but_one_edges;
	    rval = mField.get_moab().get_adjacencies(
	      crack_front_tets_faces_but_one,1,false,crack_front_tets_faces_but_one_edges,Interface::UNION); CHKERR_PETSC(rval);
	    crack_front_tets_faces_but_one_edges = intersect(crack_front_tets_faces_but_one_edges,crack_front_tets_faces_edges);
	    Range fit_edges;
	    rval = mField.get_moab().get_adjacencies(&*fit,1,1,false,fit_edges); CHKERR_PETSC(rval);
	    fit_edges = intersect(fit_edges,crack_front_tets_faces_edges);
	    fit_edges = subtract(fit_edges,crack_front_tets_faces_but_one_edges);
	    fit_edges = subtract(fit_edges,mesh_level_tets_skin_faces_edges);
	    if(fit_edges.size()>0) {
	      crack_front_tets_faces.erase(*fit);
	      break;
	    }
	  }
	}
      }

    }

    rval = mField.get_moab().get_connectivity(crack_front_tets,crack_front_tets_nodes,true); CHKERR_PETSC(rval);
    crack_front_tets_nodes = subtract(crack_front_tets_nodes,crack_front_edges_nodes);

    if(verb>=0) {
      PetscPrintf(PETSC_COMM_WORLD," (%d) %u (%u) <%u> [%4.4f]",
	ii++,crack_front_tets_nodes.size(),_nodes_on_skin_surface_.size(),crack_front_tets.size(),min_nit_q);
    }
    
    debug_ii++;

    if(verb>0) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {

      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,crack_front_tets); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,crack_front_tets_faces); CHKERR_PETSC(rval);
      ostringstream ss;
      ss << "chop_debug_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,chop_faces); CHKERR_PETSC(rval);
      ostringstream sss;
      sss << "chop_tets_debug_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(sss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,_crack_front_tets_skin_faces_); CHKERR_PETSC(rval);
      ostringstream ssss;
      ssss << "chop_crack_front_tets_skin_faces_debug_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(ssss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,_nodes_on_skin_surface_); CHKERR_PETSC(rval);
      ostringstream ssssss;
      ssssss << "chop_crack_nodes_on_skin_surface_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(ssssss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

      }

    }	

    //PetscBarrier(PETSC_NULL);

    if(crack_front_tets_faces.size() == nb_crack_front_tets_faces) { 
      SETERRQ(PETSC_COMM_SELF,1,"it is empty, algorithm is stack");
    }

  } while (!crack_front_tets.empty()); 

  if(verb>=0) {
    PetscPrintf(PETSC_COMM_WORLD,"\n",crack_front_tets_nodes.size());
  }

  //save meshset
  rval = mField.get_moab().add_entities(chopTetsFaces,crack_front_tets_faces); CHKERR_PETSC(rval);
  Range crack_front_tets_faces_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_tets_faces,crack_front_tets_faces_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_tets_faces_nodes_front;
  crack_front_tets_faces_nodes_front = subtract(crack_front_tets_faces_nodes,crack_front_edges_nodes);
  rval = mField.get_moab().add_entities(chopTetsFaces,crack_front_tets_faces_nodes_front); CHKERR_PETSC(rval);

  //calulate quality of choped faces
  double q;
  ierr = calculate_qualityAfterProjectingNodes(
    crack_front_tets_faces_nodes_front,mesh_level_tets,crack_front_edges_nodes,q); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"meshset quality = %4.3f\n",q);  

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::selectCrackFaces(bool createMeshset,int verb) {
  PetscFunctionBegin;
  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,selectedCrackFaces); CHKERR_PETSC(rval);
  }
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  //
  Range selecred_crack_surface_faces;
  //add faces form unfreezed noded
  Tag th_freez;
  rval = mField.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
  vector<int> freezed_nodes(crack_front_edges_nodes.size());
  rval = mField.get_moab().tag_get_data(th_freez,crack_front_edges_nodes,&freezed_nodes[0]); CHKERR_PETSC(rval);
  Range::iterator nit = crack_front_edges_nodes.begin();
  vector<int>::iterator vit = freezed_nodes.begin();
  for(;nit!=crack_front_edges_nodes.end();nit++,vit++) {
    if(*vit!=1) {
      Range faces;
      rval = mField.get_moab().get_adjacencies(
	&*nit,1,2,false,faces,Interface::UNION); CHKERR_PETSC(rval);
      selecred_crack_surface_faces.merge(faces);
    }
  }
  //
  Range chop_tets_faces;
  rval = mField.get_moab().get_entities_by_type(chopTetsFaces,MBTRI,chop_tets_faces,true); CHKERR_PETSC(rval);
  //
  selecred_crack_surface_faces = intersect(selecred_crack_surface_faces,chop_tets_faces);
  //
  rval = mField.get_moab().add_entities(selectedCrackFaces,selecred_crack_surface_faces); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::catMesh(const int verb) {
  PetscFunctionBegin;

  Range edges_to_refine;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,edges_to_refine,false); CHKERR_PETSC(rval);

  int current_ref_bit = meshRefineBitLevels.back();
  BitRefLevel current_ref = BitRefLevel().set(current_ref_bit);

  Range level_edges;
  ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
  edges_to_refine = intersect(edges_to_refine,mesh_level_edges);

  Range::iterator eit = edges_to_refine.begin();
  for(;eit!=edges_to_refine.end();) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double d,other_d;
    rval = mField.get_moab().tag_get_data(th_distance,&conn[0],1,&d); CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_data(th_distance,&conn[1],1,&other_d); CHKERR_PETSC(rval);
    double coords[3*num_nodes]; 
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,&coords[3],1,coords,1);
    double l = cblas_dnrm2(3,coords,1);
    if(fmin(fabs(d),fabs(other_d))/l<0.2) {
      eit = edges_to_refine.erase(*eit);
      continue;
    }
    //b = d => a*1 + d = other_d => a = other_d - d
    //a*s + b = 0 => s = -b/a => s = -d/(other_d - d)
    double s = -d/(other_d-d); 
    //cerr << s << endl;
    if( (s<0.2)||(s>0.8) ) {
      eit = edges_to_refine.erase(*eit);
      continue;
      //cerr << "EEEEEEEEEEEE\n";
    } else {
      eit++;
      continue;
    }
  }
  PetscPrintf(PETSC_COMM_WORLD,"nb. edges to refine = %u\n",edges_to_refine.size());	

  int last_ref_bit = meshRefineBitLevels.back();
  if(!meshIntefaceBitLevels.empty()) {
    if(last_ref_bit<meshIntefaceBitLevels.back()) {
      last_ref_bit = meshIntefaceBitLevels.back();
    }
  }
  last_ref_bit++;
  meshRefineBitLevels.push_back(last_ref_bit);     
  BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.add_verices_in_the_middel_of_edges(edges_to_refine,last_ref); CHKERRQ(ierr);

  Range level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);
  ierr = mField.refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
    EntityHandle cubit_meshset = cubit_it->meshset; 
    ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::meshRefine(const int verb) {
  PetscFunctionBegin;

  PetscBool flg = PETSC_TRUE;
  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    Tag th_set_ref_level;
    rval = mField.get_moab().tag_get_handle("_SET_REF_LEVEL",th_set_ref_level); CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_data(th_set_ref_level,&root_meshset,1,&nb_ref_levels); CHKERR_PETSC(rval);
    //cerr << "nb_ref_levels " << nb_ref_levels << endl;
  }

  BitRefLevel preserve_ref = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  const RefMoFEMEntity_multiIndex *refinedMoFemEntities_ptr;
  ierr = mField.get_ref_ents(&refinedMoFemEntities_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::index<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>::type::iterator refit,hi_refit;
  refit = refinedMoFemEntities_ptr->get<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>().lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  hi_refit = refinedMoFemEntities_ptr->get<Composite_EntType_mi_tag_and_ParentEntType_mi_tag>().upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  Range already_refined_edges;
  for(;refit!=hi_refit;refit++) {
    EntityHandle parent_ent = refit->get_parent_ent(); 
    already_refined_edges.insert(parent_ent);
    RefMoFEMEntity_multiIndex::iterator parent_rit;
    parent_rit = refinedMoFemEntities_ptr->find(parent_ent);
    if(parent_rit == refinedMoFemEntities_ptr->end()) {
      SETERRQ1(PETSC_COMM_SELF,1,
	  "data inconsistency, entity in database not found %lu",refit->get_parent_ent());
    }
    bool success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
	->modify(parent_rit,RefMoFEMEntity_change_add_bit(preserve_ref));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessfull");
    }
  }

  Range preserve_ref_tets;
  ierr = mField.get_entities_by_type_and_ref_level(
    preserve_ref,BitRefLevel().set(),MBTET,preserve_ref_tets); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(preserve_ref_tets); CHKERRQ(ierr);

  int current_ref_bit = meshRefineBitLevels.first();
  for(int ll = 1;ll<nb_ref_levels+1;ll++) {

      BitRefLevel current_ref = BitRefLevel().set(current_ref_bit);

      Range level_tets;
      ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
      ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);

      Range level_edges;
      ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

      Range crack_edges;
      ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_edges,true); CHKERRQ(ierr);
      Range crack_edges_nodes;
      rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);

      //PetscAttachDebugger();
      Range crack_edge_nodes_parents;
      for(Range::iterator nit = crack_edges_nodes.begin();nit!=crack_edges_nodes.end();nit++) {
	EntityHandle ent = *nit;
	do {
	  RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator refit;
	  refit = refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().find(ent);
	  if(refit == refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().end()) {
	    break;
	  }
	  if(refit->get_parent_ent() != 0) {
	    ent = refit->get_parent_ent();
	    crack_edge_nodes_parents.insert(ent);
	  }
	} while(refit->get_parent_ent() != 0);
      }
      Range crack_edge_nodes_parents_nodes;
      rval = mField.get_moab().get_connectivity(crack_edge_nodes_parents,crack_edge_nodes_parents_nodes,true); CHKERR_PETSC(rval);
      crack_edges_nodes.merge(crack_edge_nodes_parents_nodes);

      Range crack_edges_nodes_tets;
      rval = mField.get_moab().get_adjacencies(
	crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
      Range edges_to_refine;
      rval = mField.get_moab().get_adjacencies(crack_edges_nodes_tets,1,false,edges_to_refine,Interface::UNION); CHKERR_PETSC(rval);
      edges_to_refine.merge(already_refined_edges);
      edges_to_refine = intersect(edges_to_refine,level_edges);
      if(edges_to_refine.empty()) continue; 

      int last_ref_bit = meshRefineBitLevels.back();
      if(!meshIntefaceBitLevels.empty()) {
	if(last_ref_bit<meshIntefaceBitLevels.back()) {
	  last_ref_bit = meshIntefaceBitLevels.back();
	}
      }
      last_ref_bit++;
      meshRefineBitLevels.push_back(last_ref_bit);     
      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
      ierr = mField.add_verices_in_the_middel_of_edges(edges_to_refine,last_ref,2); CHKERRQ(ierr);
      ierr = mField.refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);
      ierr = mField.seed_ref_level_3D(level_tets,preserve_ref); CHKERRQ(ierr);
      
      for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
	EntityHandle cubit_meshset = cubit_it->meshset; 
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
      }

      current_ref_bit = meshRefineBitLevels.back();
  
  }

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("debug_mesh_refine.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::splitFaces(const int verb) {
  PetscFunctionBegin;

  const RefMoFEMEntity_multiIndex *refinedMoFemEntities_ptr;
  ierr = mField.get_ref_ents(&refinedMoFemEntities_ptr); CHKERRQ(ierr);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  BitRefLevel back_up_level;

  BitRefLevel current_ref = BitRefLevel().set(meshRefineBitLevels.back());
  BitRefLevel inheret_ents_from_level,inheret_ents_from_level_mask;
  if(!meshIntefaceBitLevels.empty()) {
    inheret_ents_from_level.set(meshIntefaceBitLevels.back());
    inheret_ents_from_level.set(BITREFLEVEL_SIZE-1);
  }
  inheret_ents_from_level_mask.set();

  Range interface_elements;
  ierr = mField.get_entities_by_type_and_ref_level(
    inheret_ents_from_level,inheret_ents_from_level_mask,MBTET,interface_elements); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(interface_elements); CHKERRQ(ierr);

  EntityHandle bit_meshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,bit_meshset); CHKERR_PETSC(rval);
  {

    ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,bit_meshset); CHKERRQ(ierr);
    ierr = mField.seed_finite_elements(bit_meshset); CHKERRQ(ierr);

    EntityHandle meshset_interface;
    ierr = mField.get_Cubit_msId_meshset(200,SideSet,meshset_interface); CHKERRQ(ierr);
    ierr = mField.get_msId_3dENTS_sides(meshset_interface,current_ref,true); CHKERRQ(ierr);

    int last_ref_bit = meshRefineBitLevels.back();
    if(!meshIntefaceBitLevels.empty()) {
      if(last_ref_bit<meshIntefaceBitLevels.back()) {
	last_ref_bit = meshIntefaceBitLevels.back();
      }
    }
    last_ref_bit++;
    meshIntefaceBitLevels.push_back(last_ref_bit);
    BitRefLevel last_ref = BitRefLevel().set(last_ref_bit);
  
    ierr = mField.get_msId_3dENTS_split_sides(
      bit_meshset,last_ref,
      inheret_ents_from_level,inheret_ents_from_level_mask,
      meshset_interface,false,true); CHKERRQ(ierr);

    //add refined ent to cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }

    //remove tets which have 4 nodes on crack surface, those are surce of problems. In case of planar crack, volume of such 
    //tets is negative. In some other cases quality could be negative.
    Range last_ref_tets,last_ref_tris;
    ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,last_ref_tets); CHKERRQ(ierr);
    ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTRI,last_ref_tris); CHKERRQ(ierr);
    Range interface_tris;
    rval = mField.get_moab().get_entities_by_type(meshset_interface,MBTRI,interface_tris,true); CHKERR_PETSC(rval);
    interface_tris = intersect(interface_tris,last_ref_tris);
    Range interface_tris_nodes;
    rval = mField.get_moab().get_connectivity(interface_tris,interface_tris_nodes,true); CHKERR_PETSC(rval);
    Range interface_tris_nodes_tets;
    rval = mField.get_moab().get_adjacencies(interface_tris_nodes,3,false,interface_tris_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
    interface_tris_nodes_tets = intersect(interface_tris_nodes_tets,last_ref_tets);
    Range interface_tris_nodes_tets_nodes;
    rval = mField.get_moab().get_connectivity(interface_tris_nodes_tets,interface_tris_nodes_tets_nodes,true); CHKERR_PETSC(rval);
    Range interface_tris_nodes_tets_nodes_minus_surface_nodes;
    interface_tris_nodes_tets_nodes_minus_surface_nodes = subtract(interface_tris_nodes_tets_nodes,interface_tris_nodes);
    Range interface_tris_nodes_tets_nodes_minus_surface_nodes_tets;
    rval = mField.get_moab().get_adjacencies(
      interface_tris_nodes_tets_nodes_minus_surface_nodes,
      3,false,interface_tris_nodes_tets_nodes_minus_surface_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
    interface_tris_nodes_tets_nodes_minus_surface_nodes_tets = intersect(
      interface_tris_nodes_tets_nodes_minus_surface_nodes_tets,last_ref_tets);

    Range tets_on_surface = subtract(interface_tris_nodes_tets,interface_tris_nodes_tets_nodes_minus_surface_nodes_tets);
    if(tets_on_surface.size()>0) {
      Range tets_on_surface_faces;
      rval = mField.get_moab().get_adjacencies(
	tets_on_surface,2,false,tets_on_surface_faces,Interface::UNION); CHKERR_PETSC(rval);
      tets_on_surface_faces = intersect(tets_on_surface_faces,interface_tris);
      
      if(pcomm->rank()>0) {
	if(verb>3) {    
	  EntityHandle out_meshset;
	  rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	  rval = mField.get_moab().add_entities(out_meshset,tets_on_surface); CHKERR_PETSC(rval);
	  //rval = mField.get_moab().add_entities(out_meshset,interface_tris); CHKERR_PETSC(rval);
	  //rval = mField.get_moab().add_entities(out_meshset,tets_on_surface_faces); CHKERR_PETSC(rval);
	  rval = mField.get_moab().write_file("tets_on_surface.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	  rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
	}
      }

      for(Range::iterator tit = tets_on_surface.begin();tit!=tets_on_surface.end();tit++) {

	RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator mit;
	mit = refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().find(*tit);
	if(mit == refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"no such tet in database");
	}
	bool success;
	success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
	  ->modify(mit,RefMoFEMEntity_change_set_nth_bit(last_ref_bit,false));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
	success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
	  ->modify(mit,RefMoFEMEntity_change_set_nth_bit(BITREFLEVEL_SIZE-1,true));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");

      }
      if(verb>=0) {    
	PetscPrintf(PETSC_COMM_WORLD,"number of block quasi-flat (on crack surface) tets: %u\n",tets_on_surface.size());
      }

      Range level_nodes;
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
      Range level_tets;
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);

      Range ents_to_set_level;

      Range level_tets_nodes;
      rval = mField.get_moab().get_connectivity(level_tets,level_tets_nodes,true); CHKERR_PETSC(rval);
      ents_to_set_level.merge(subtract(level_nodes,level_tets_nodes));
      if(verb>=0) {    
	PetscPrintf(PETSC_COMM_WORLD,"number of block entities of quasi-flat (on crack surface) tets: %u\n",ents_to_set_level.size());
      }

      for(Range::iterator eit = ents_to_set_level.begin();eit!=ents_to_set_level.end();eit++) {
	RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator mit;
	mit = refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().find(*eit);
	if(mit == refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().end()) {
	  SETERRQ(PETSC_COMM_SELF,1,"no such tet in database");
	}
	bool success;
	//success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
	//  ->modify(mit,RefMoFEMEntity_change_set_nth_bit(last_ref_bit,false));
	success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
	  ->modify(mit,RefMoFEMEntity_change_set_nth_bit(BITREFLEVEL_SIZE-1,true));
	if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
      }

    }

  }

  rval = mField.get_moab().delete_entities(&bit_meshset,1); CHKERR_PETSC(rval);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;
  bit_level0 = BitRefLevel().set(meshIntefaceBitLevels.back());

  if(verb>0) {    
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      BitRefLevel last_ref = BitRefLevel().set(meshIntefaceBitLevels.back());
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("debug_split_faces.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
  }

  BitRefLevel levels_to_preserve;
  {
    int *p = meshIntefaceBitLevels.begin();
    for(;p!=meshIntefaceBitLevels.end();p++) {
      levels_to_preserve.set(*p);
    }
  }
  Range ents_to_preserve;
  ierr = mField.get_entities_by_ref_level(
    levels_to_preserve,BitRefLevel().set(),ents_to_preserve); CHKERRQ(ierr);
  BitRefLevel preserve_ref = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  for(Range::iterator eit = ents_to_preserve.begin();
    eit!=ents_to_preserve.end();eit++) {
    RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator mit;
    mit = refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().find(*eit);
    if(mit == refinedMoFemEntities_ptr->get<MoABEnt_mi_tag>().end()) {
      SETERRQ1(
	PETSC_COMM_SELF,1,"no such ent in database, type %lu",
	mField.get_moab().type_from_handle(*eit));
    }
    bool success;
    success = const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)
      ->modify(mit,RefMoFEMEntity_change_add_bit(preserve_ref));
    if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::addNewSurfaceFaces_to_Cubit_msId200() {
  PetscFunctionBegin;
  EntityHandle meshset200;
  ierr = mField.get_Cubit_msId_meshset(200,SideSet,meshset200); CHKERRQ(ierr);
  Range new_crack_surface_faces;
  rval = mField.get_moab().get_entities_by_type(selectedCrackFaces,MBTRI,new_crack_surface_faces,false); CHKERR_PETSC(rval);
  rval = mField.get_moab().add_entities(meshset200,new_crack_surface_faces); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::addcrackFront_to_Cubit201() {
  PetscFunctionBegin;
  Range outer_surface_skin;
  Skinner skin(&mField.get_moab());
  rval = skin.find_skin(mesh_level_tets,false,outer_surface_skin); CHKERR(rval);
  Range outer_surface_skin_edges;
  rval = mField.get_moab().get_adjacencies(outer_surface_skin,1,false,outer_surface_skin_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_surface_tris;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_tris_skin_edges;
  rval = skin.find_skin(crack_surface_tris,false,crack_surface_tris_skin_edges); CHKERR(rval);
  crack_surface_tris_skin_edges = subtract(crack_surface_tris_skin_edges,outer_surface_skin_edges);
  if(mField.check_msId_meshset(201,SideSet)) {
    ierr = mField.delete_Cubit_msId(SideSet,201); CHKERRQ(ierr);
  }
  ierr = mField.add_Cubit_msId(SideSet,201); CHKERRQ(ierr);
  EntityHandle meshset201;
  ierr = mField.get_Cubit_msId_meshset(201,SideSet,meshset201); CHKERRQ(ierr);
  ierr = mField.get_moab().add_entities(meshset201,crack_surface_tris_skin_edges); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::projectCrackFrontNodes() {
  PetscFunctionBegin;

 //material prositions
  const DofMoFEMEntity_multiIndex *dofs_moabfield_ptr;
  ierr = mField.get_dofs(&dofs_moabfield_ptr); CHKERRQ(ierr);
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dof_iterator;
  //
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range::iterator nit = crack_front_edges_nodes.begin();
  for(;nit!=crack_front_edges_nodes.end();nit++) {
    ublas::vector<double,ublas::bounded_array<double,3> > new_coords;
    new_coords.resize(3);
    rval = mField.get_moab().tag_get_data(th_projection,&*nit,1,&*new_coords.data().begin()); CHKERR_PETSC(rval);
    dof_iterator dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));
    dof_iterator hi_dit = 
      dofs_moabfield_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("MESH_NODE_POSITIONS",*nit));
    if(distance(dit,hi_dit)!=3) {
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistencies");
    }
    for(;dit!=hi_dit;dit++) {
      dit->get_FieldData() = new_coords[dit->get_dof_rank()];
    }
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::squashIndices(const int verb) {
  PetscFunctionBegin;

  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"meshRefineBitLevels: ");
    int *p = meshRefineBitLevels.begin();
    for(;p!=meshRefineBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");

    PetscPrintf(PETSC_COMM_WORLD,"meshIntefaceBitLevels: ");
    p = meshIntefaceBitLevels.begin();
    for(;p!=meshIntefaceBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  }

  const RefMoFEMEntity_multiIndex *refinedMoFemEntities_ptr;
  ierr = mField.get_ref_ents(&refinedMoFemEntities_ptr); CHKERRQ(ierr);

  BitRefLevel maskPreserv;
  ierr = getMask(maskPreserv); CHKERRQ(ierr);
  BitRefLevel mask = maskPreserv;
  mask.flip();

  PetscBool flg = PETSC_TRUE;
  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    Tag th_set_ref_level;
    rval = mField.get_moab().tag_get_handle("_SET_REF_LEVEL",th_set_ref_level); CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_data(th_set_ref_level,&root_meshset,1,&nb_ref_levels); CHKERR_PETSC(rval);
  }

  vector<int> new_meshRefineBitLevels;
  vector<int> new_meshRefineBitLevels_map;
  new_meshRefineBitLevels.push_back(0);
  new_meshRefineBitLevels_map.push_back(0);
  new_meshRefineBitLevels.push_back(0);
  new_meshRefineBitLevels_map.push_back(meshRefineBitLevels.first());
  for(int nn = 1;nn<=nb_ref_levels;nn++) {
    new_meshRefineBitLevels.push_back(nn);
    int *p = (meshRefineBitLevels.end()-1)-(nb_ref_levels-nn);
    new_meshRefineBitLevels_map.push_back(*p);
  }
  new_meshRefineBitLevels[0] = new_meshRefineBitLevels.size()-1;

  vector<int> new_meshIntefaceBitLevels;
  vector<int> new_meshIntefaceBitLevels_map;
  new_meshIntefaceBitLevels.push_back(0);
  new_meshIntefaceBitLevels_map.push_back(0);
  new_meshIntefaceBitLevels.push_back(new_meshRefineBitLevels.size()-1);
  new_meshIntefaceBitLevels_map.push_back(meshIntefaceBitLevels.back());
  new_meshIntefaceBitLevels[0] = new_meshIntefaceBitLevels.size()-1;

  //set new memory
  bzero(meshRefineBitLevels.ptr,(BITREFLEVEL_SIZE-1)*sizeof(int));
  bcopy(&*new_meshRefineBitLevels.begin(),meshRefineBitLevels.ptr,new_meshRefineBitLevels.size()*sizeof(int));
  bzero(meshIntefaceBitLevels.ptr,(BITREFLEVEL_SIZE-1)*sizeof(int));
  bcopy(&*new_meshIntefaceBitLevels.begin(),meshIntefaceBitLevels.ptr,new_meshIntefaceBitLevels.size()*sizeof(int));

  if(verb>1) {
    PetscPrintf(PETSC_COMM_WORLD,"meshRefineBitLevels: ");
    int *p = meshRefineBitLevels.begin();
    for(;p!=meshRefineBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"meshIntefaceBitLevels: ");
    p = meshIntefaceBitLevels.begin();
    for(;p!=meshIntefaceBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  }

  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"meshRefineBitLevels (map): ");
    vector<int>::iterator vit,viit;
    vit = new_meshRefineBitLevels.begin();
    viit = new_meshRefineBitLevels_map.begin();
    for(;vit!=new_meshRefineBitLevels.end();vit++,viit++) {
      PetscPrintf(PETSC_COMM_WORLD," %d (%d)",*vit,*viit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"meshIntefaceBitLevels (map): ");
    vit = new_meshIntefaceBitLevels.begin();
    viit = new_meshIntefaceBitLevels_map.begin();
    for(;vit!=new_meshIntefaceBitLevels.end();vit++,viit++) {
      PetscPrintf(PETSC_COMM_WORLD," %d (%d)",*vit,*viit);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  }

  //squash bits
  RefMoFEMEntity_multiIndex::iterator mit = refinedMoFemEntities_ptr->begin();

  for(;mit!=refinedMoFemEntities_ptr->end();mit++) {
    
    if(mit->get_ent_type() == MBENTITYSET) continue;
    if(mit->get_BitRefLevel().none()) continue;

    if( (mask&mit->get_BitRefLevel()).none() ) {
      cerr << "mask:\n" << mask << endl;
      cerr << "mit->get_BitRefLevel():\n" << mit->get_BitRefLevel() << endl;
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistencies");
    }
    
    BitRefLevel new_bit;
    for(unsigned int nn = 1;nn<new_meshRefineBitLevels.size();nn++) {
      if((mit->get_BitRefLevel()&
	BitRefLevel().set(new_meshRefineBitLevels_map[nn])).any()) {
	new_bit.set(new_meshRefineBitLevels[nn]);
      }
    }

    for(unsigned int nn = 1;nn<new_meshIntefaceBitLevels.size();nn++) {
      if((mit->get_BitRefLevel()&
	BitRefLevel().set(new_meshIntefaceBitLevels_map[nn])).any()) {
	new_bit.set(new_meshIntefaceBitLevels[nn]);
      }
    }
    if( (mit->get_BitRefLevel()&
      BitRefLevel().set(BITREFLEVEL_SIZE-1)).any() ) {
      new_bit.set(BITREFLEVEL_SIZE-1);
    }
 
    if( (mit->get_BitRefLevel()&
      BitRefLevel().set(BITREFLEVEL_SIZE-2)).any() ) {
      new_bit.set(BITREFLEVEL_SIZE-2);
    } 

    /*cerr << *mit << endl;
    cerr << "mit->get_BitRefLevel():\n" << mit->get_BitRefLevel() << endl;
    cerr << "new_bit:\n" << new_bit << endl;*/

    if(new_bit.none()) {
      cerr << "mit->get_BitRefLevel():\n" << mit->get_BitRefLevel() << endl;
      SETERRQ(PETSC_COMM_SELF,1,"data inconsistencies");
    }

    //cerr << mit->get_BitRefLevel() << " : " << new_bit << endl;

    const_cast<RefMoFEMEntity_multiIndex*>(refinedMoFemEntities_ptr)->modify(mit,RefMoFEMEntity_change_set_bit(new_bit));
  }

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;
  bit_level0 = BitRefLevel().set(meshIntefaceBitLevels.back());

  if(verb>2) {
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
    rval = mField.get_moab().write_file("squash_mesh.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculate_qualityAfterProjectingNodes(
  Range &option_nodes,Range &intersect_tets,Range &crack_front_edges_nodes,double &current_q) {
  PetscFunctionBegin;
  Range opposite_crack_front_edges;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,opposite_crack_front_edges,false); CHKERR_PETSC(rval);
  double def_VAL[1] = { 0 };
  rval = mField.get_moab().tag_get_handle("MESHSET_PROJECTION_B",1,MB_TYPE_DOUBLE,th_b,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
  ierr = ShapeDiffMBTET(diffNTET); CHKERRQ(ierr);
  map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,3> > > projectedNodes;
  for(Range::iterator nit = option_nodes.begin();nit!=option_nodes.end();nit++) {
    ublas::vector<double,ublas::bounded_array<double,3> > 
      &closest_point_out = projectedNodes[*nit];
    closest_point_out.resize(3);
    rval = mField.get_moab().tag_get_data(
      th_projection,&*nit,1,&*closest_point_out.data().begin()); CHKERR_PETSC(rval);
  }
  map<EntityHandle,double> q_map;  
  Range adj_tets;
  rval = mField.get_moab().get_adjacencies(option_nodes,3,false,adj_tets,Interface::UNION); CHKERR_PETSC(rval);
  adj_tets = intersect(adj_tets,intersect_tets);
  for(Range::iterator tit = adj_tets.begin();tit!=adj_tets.end();tit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes]; 
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    double dofs_X[3*num_nodes];
    cblas_dcopy(3*num_nodes,coords,1,dofs_X,1);
    for(int nn = 0;nn<num_nodes;nn++) {
      if(crack_front_edges_nodes.find(conn[nn])!=crack_front_edges_nodes.end()) {
	continue;
      }
      map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,3> > >::iterator 
	mit = projectedNodes.find(conn[nn]);
      if(mit != projectedNodes.end() ) {
	cblas_dcopy(3,&*mit->second.data().begin(),1,&dofs_X[3*nn],1);
      }
    }
    double quality0 = 0,quality,b;
    {
      double coords_edges[2*3*6]; 
      ierr = get_edges_from_elem_coords(coords,coords_edges); CHKERRQ(ierr);
      double V =  Shape_intVolumeMBTET(diffNTET,&*coords); 
      /*cerr << "V " << V;
      double V1 =  Shape_intVolumeMBTET(diffNTET,dofs_X); 
      cerr << " V1 " << V1 << endl;*/
      double alpha[4] = {1,1,1,1};
      ierr = quality_volume_length_F(
	  V,alpha,0,
	  diffNTET,coords_edges,dofs_X,
	  NULL,NULL,NULL,
	  &quality0,&quality,&b,
	  NULL,NULL); 
      if(ierr != 0) {
	b = -2;
	quality = -2;
	ierr = 0;
	/*PetscSynchronizedPrintf(PETSC_COMM_WORLD,"dofs_X = [ %6.4e  %6.4e %6.4e %6.4e]\n",dofs_X[0],dofs_X[1],dofs_X[2],dofs_X[3]);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"dofs_X = [ %6.4e  %6.4e %6.4e %6.4e]\n",dofs_X[4],dofs_X[5],dofs_X[6],dofs_X[7]);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"dofs_X = [ %6.4e  %6.4e %6.4e %6.4e]\n",dofs_X[8],dofs_X[9],dofs_X[10],dofs_X[11]);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD,"dofs_X = [ %6.4e  %6.4e %6.4e %6.4e]\n",dofs_X[12],dofs_X[13],dofs_X[14],dofs_X[15]);
	PetscSynchronizedFlush(PETSC_COMM_WORLD);*/
	//CHKERRQ(ierr);
      }
    }
    Range new_front_nodes;
    new_front_nodes.insert(&conn[0],&conn[4]);
    new_front_nodes = intersect(new_front_nodes,option_nodes);
    for(Range::iterator nit = new_front_nodes.begin();nit!=new_front_nodes.end();nit++) {
      map<EntityHandle,double>::iterator mit = q_map.find(*nit);
      if(mit!=q_map.end()) {
	mit->second = fmin(mit->second,quality);
      } else {
	q_map[*nit] = quality;
      }
    }
    if(tit != adj_tets.begin()) {
      quality = fmin(quality,current_q);
    }
    current_q = quality;
  }
  for(map<EntityHandle,double>::iterator mit = q_map.begin();mit!=q_map.end();mit++) {
    //PetscPrintf(PETSC_COMM_WORLD,"node %u b = %3.2f\n",mit->first,mit->second);
    rval = mField.get_moab().tag_set_data(th_b,&mit->first,1,&mit->second); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getMask(BitRefLevel &maskPreserv,const int verb) {
  PetscFunctionBegin;

  PetscBool flg = PETSC_TRUE;
  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    const EntityHandle root_meshset = mField.get_moab().get_root_set();
    Tag th_set_ref_level;
    rval = mField.get_moab().tag_get_handle("_SET_REF_LEVEL",th_set_ref_level); CHKERR_PETSC(rval);
    rval = mField.get_moab().tag_get_data(th_set_ref_level,&root_meshset,1,&nb_ref_levels); CHKERR_PETSC(rval);
    //cerr << "nb_ref_levels " << nb_ref_levels << endl;
  }

  if(meshRefineBitLevels.size()<nb_ref_levels) {
    SETERRQ2(
      PETSC_COMM_SELF,1,"number of mesh levelsmeshRefineBitLevels.size()<nb_ref_levels %d < %d",
      meshRefineBitLevels.size(),nb_ref_levels);
  }

  if(verb>0) {
    PetscPrintf(PETSC_COMM_WORLD,"meshRefineBitLevels: ");
    int *p = meshRefineBitLevels.begin();
    for(;p!=meshRefineBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
    PetscPrintf(PETSC_COMM_WORLD,"meshIntefaceBitLevels: ");
    p = meshIntefaceBitLevels.begin();
    for(;p!=meshIntefaceBitLevels.end();p++) {
      PetscPrintf(PETSC_COMM_WORLD," %d",*p);
    }
    PetscPrintf(PETSC_COMM_WORLD,"\n");
  }

  maskPreserv.set();
  maskPreserv[meshRefineBitLevels.first()] = false;
  maskPreserv[meshIntefaceBitLevels.back()] = false;

  maskPreserv[BITREFLEVEL_SIZE-1] = false;
  maskPreserv[BITREFLEVEL_SIZE-2] = false;

  {
    int *p = meshRefineBitLevels.end();
    for(int nn = 0;nn<nb_ref_levels;nn++) {
      --p;
      maskPreserv[*p] = false;
    }
  }

  /*{
    int *p = meshIntefaceBitLevels.begin();
    for(;p!=meshIntefaceBitLevels.end();p++) {
      maskPreserv[*p] = false;
    }
  }*/

  if(verb>0) {
    ostringstream s;
    s << "maskPreserv: " << maskPreserv << endl;
    PetscPrintf(PETSC_COMM_WORLD,s.str().c_str());
  }

  PetscFunctionReturn(0);
}

PetscErrorCode main_refine_and_meshcat(FieldInterface& mField,FaceSplittingTools &face_splitting,bool cat_mesh,const int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;
 
  BitRefLevel bit_last_ref = BitRefLevel().set(face_splitting.meshRefineBitLevels.back());
  Range tets_on_last_ref_level;
  ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,tets_on_last_ref_level); CHKERRQ(ierr);
  ierr = mField.seed_ref_level_3D(tets_on_last_ref_level,bit_last_ref); CHKERRQ(ierr);
  ierr = face_splitting.buildKDTreeForCrackSurface(bit_level0); CHKERRQ(ierr);
  ierr = face_splitting.meshRefine(); CHKERRQ(ierr);
  bit_last_ref = BitRefLevel().set(face_splitting.meshRefineBitLevels.back());
  ierr = face_splitting.initBitLevelData(bit_last_ref);  CHKERRQ(ierr);
  ierr = face_splitting.calculateDistanceFromCrackSurface();  CHKERRQ(ierr);
  ierr = face_splitting.getOpositeForntEdges(true); CHKERRQ(ierr);
  if(cat_mesh) {
    ierr = face_splitting.catMesh(); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode main_select_faces_for_splitting(FieldInterface& mField,FaceSplittingTools &face_splitting,const int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  BitRefLevel bit_last_ref = BitRefLevel().set(face_splitting.meshRefineBitLevels.back());
  Range tets_on_last_ref_level;
  ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,tets_on_last_ref_level); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(tets_on_last_ref_level); CHKERRQ(ierr);

  ierr = face_splitting.buildKDTreeForCrackSurface(bit_level0); CHKERRQ(ierr);
  ierr = face_splitting.initBitLevelData(bit_last_ref);  CHKERRQ(ierr);
  ierr = face_splitting.calculateDistanceFromCrackSurface();  CHKERRQ(ierr);

  ierr = face_splitting.getOpositeForntEdges(true); CHKERRQ(ierr);

  if(verb>0) {

    if(pcomm->rank()==1) {

      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      //
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("out0.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      //
      rval = mField.get_moab().write_file("opositeFrontEdges.vtk","VTK","",&face_splitting.opositeFrontEdges,1); CHKERR_PETSC(rval);

    }

  }

  ierr = face_splitting.getCrackFrontTets(true,0); CHKERRQ(ierr);
  ierr = face_splitting.chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(true,0); CHKERRQ(ierr);
  ierr = face_splitting.selectCrackFaces(true); CHKERRQ(ierr);
  ierr = face_splitting.getCrackFrontFaces(true,10); CHKERRQ(ierr);

  if(verb>0) {

    if(pcomm->rank()==0) {

      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(bit_last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      //
      rval = mField.get_moab().write_file("opositeFrontEdges.vtk","VTK","",&face_splitting.opositeFrontEdges,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("crackFrontTests.vtk","VTK","",&face_splitting.crackFrontTests,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("chopTetsFaces.vtk","VTK","",&face_splitting.chopTetsFaces,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("selectedCrackFaces.vtk","VTK","",&face_splitting.selectedCrackFaces,1); CHKERR_PETSC(rval);

    }
  
  }

  PetscFunctionReturn(0);
}

PetscErrorCode main_split_faces_and_update_field_and_elements(FieldInterface& mField,FaceSplittingTools &face_splitting,const int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  BitRefLevel maskPreserv_flip;
  ierr = face_splitting.getMask(maskPreserv_flip); CHKERRQ(ierr);
  maskPreserv_flip.flip();
  ierr = mField.build_adjacencies(maskPreserv_flip); CHKERRQ(ierr);

  //remove all crack front elements
  ierr = mField.remove_ents_from_finite_element("C_TANGENT_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("C_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("CTC_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("dCT_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("C_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("CTC_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("CandCT_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_finite_element("CandCT_CRACK_SURFACE_ELEM_WITH_CRACK_FRONT",0,MBTRI); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACKFRONT_AREA",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACKFRONT_AREA",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACK_SURFACE",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("GRIFFITH_FORCE",0,MBVERTEX); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field("GRIFFITH_FORCE_TANGENT",0,MBVERTEX); CHKERRQ(ierr);

  BitRefLevel maskPreserv;
  ierr = face_splitting.getMask(maskPreserv,1); CHKERRQ(ierr);
  ierr = mField.delete_ents_by_bit_ref(maskPreserv,maskPreserv); CHKERRQ(ierr);
  ierr = face_splitting.squashIndices(0); CHKERRQ(ierr);
 
  BitRefLevel not_split_face_ref_level;
  not_split_face_ref_level.set(face_splitting.meshIntefaceBitLevels.back());
  not_split_face_ref_level.set(BITREFLEVEL_SIZE-1);
  not_split_face_ref_level.flip();
  ierr = mField.remove_ents_from_finite_element_by_bit_ref(not_split_face_ref_level,not_split_face_ref_level,1); CHKERRQ(ierr);
  ierr = mField.remove_ents_from_field_by_bit_ref(not_split_face_ref_level,not_split_face_ref_level,1); CHKERRQ(ierr);

  ierr = face_splitting.addNewSurfaceFaces_to_Cubit_msId200(); CHKERRQ(ierr);
  ierr = face_splitting.splitFaces(0); CHKERRQ(ierr);
  ierr = face_splitting.addcrackFront_to_Cubit201(); CHKERRQ(ierr);

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      Range tets;
      ierr = mField.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,tets); CHKERRQ(ierr); 
      rval = mField.get_moab().add_entities(out_meshset,tets); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("debug_meshs_after_delete.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
  }

  PetscFunctionReturn(0);
}

}

