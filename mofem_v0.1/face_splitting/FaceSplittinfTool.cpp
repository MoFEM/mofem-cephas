/* Copyright (C) 2014, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#include<FaceSplittinfTool.hpp>
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

  //material prositions
  const DofMoFEMEntity_multiIndex *dofs_moabfield_ptr;
  ierr = mField.get_dofs(&dofs_moabfield_ptr); CHKERRQ(ierr);
  typedef DofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dof_iterator;

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
    double coords[3];
    for(;dit!=hi_dit;dit++) {
      coords[dit->get_dof_rank()] = dit->get_FieldData();
    }
    //double coords[3]; 
    //rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    rval = moab_distance_from_crack_surface.create_vertex(coords,map_nodes[*nit]); CHKERR_PETSC(rval);
  }

  double diffN[3*2];
  ierr = ShapeDiffMBTRI(diffN); CHKERRQ(ierr);

  Tag th_normal;
  double def_VAL[3] = { 0,0,0 };
  rval = moab_distance_from_crack_surface.tag_get_handle("NORMAL",3,MB_TYPE_DOUBLE,th_normal,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

  Range new_entities;
  for(Range::iterator eit = entities.begin();eit!=entities.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    EntityHandle new_conn[num_nodes];
    for(int nn = 0;nn<num_nodes;nn++) {
      new_conn[nn] = map_nodes[conn[nn]];
    }
    EntityHandle new_elem;
    rval = moab_distance_from_crack_surface.create_element(
      mField.get_moab().type_from_handle(*eit),new_conn,num_nodes,new_elem); CHKERR_PETSC(rval);
    new_entities.insert(new_elem);
    rval = moab_distance_from_crack_surface.get_connectivity(new_elem,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes]; 
    rval = moab_distance_from_crack_surface.get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    double normal[3];
    ierr = ShapeFaceNormalMBTRI(diffN,coords,normal); CHKERRQ(ierr);
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
  const int def_side[] = {0};
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
PetscErrorCode FaceSplittingTools::calculateDistanceFromCrackSurface(Range &nodes) {
  PetscFunctionBegin;	

  nodes = intersect(nodes,mesh_level_nodes);

  double def_VAL[1] = { 0 };
  rval = mField.get_moab().tag_get_handle(
    "DISTANCE_FROM_CRACK_SURFACE",1,MB_TYPE_DOUBLE,th_distance,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
  Tag th_normal;
  rval = moab_distance_from_crack_surface.tag_get_handle("NORMAL",th_normal); CHKERR_THROW(rval);

  for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
    double coords[3]; 
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    double closest_point_out[3];
    EntityHandle triangle_out;
    rval = kdTree_DistanceFromCrackSurface.closest_triangle(
      kdTree_rootMeshset_DistanceFromCrackSurface,coords,closest_point_out,triangle_out); CHKERR_PETSC(rval);
    double distance[3];
    cblas_dcopy(3,coords,1,distance,1);
    cblas_daxpy(3,-1,closest_point_out,1,distance,1);
    double normal[3];
    rval = moab_distance_from_crack_surface.tag_get_data(th_normal,&triangle_out,1,normal); CHKERR_PETSC(rval);
    double dot  = cblas_ddot(3,normal,1,distance,1);
    double d = copysign(cblas_dnrm2(3,distance,1),dot);
    rval = mField.get_moab().tag_set_data(th_distance,&*nit,1,&d); CHKERR_PETSC(rval);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculateDistanceFromCrackSurface() {
  PetscFunctionBegin;
  Range nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets,nodes,true); CHKERR_PETSC(rval);
  ierr = calculateDistanceFromCrackSurface(nodes); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getOpositeForntEdges(bool createMeshset) {
  PetscFunctionBegin;

  //create meshset
  EntityHandle newFrontEdges_child;
  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,opositeFrontEdges); CHKERR_PETSC(rval);
    rval = mField.get_moab().create_meshset(MESHSET_SET,nodesOnCrackSurface); CHKERR_PETSC(rval);
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

  //crack_front_test->skin faces->skin_edges;
  Skinner skin(&mField.get_moab());
  Range skin_faces; 
  rval = skin.find_skin(crack_front_adj_tets,false,skin_faces); CHKERR(rval);
  Range skin_edges;
  rval = mField.get_moab().get_adjacencies(skin_faces,1,false,skin_edges,Interface::UNION); CHKERR_PETSC(rval);

  //crack surface faces->nodes->edges;
  Range crack_surface;
  ierr = mField.get_Cubit_msId_entities_by_dimension(200,SideSet,2,crack_surface,true); CHKERRQ(ierr);
  Range crack_surface_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface,crack_surface_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_edges;
  rval = mField.get_moab().get_adjacencies(crack_surface_nodes,1,false,crack_surface_edges,Interface::UNION); CHKERR_PETSC(rval);

  skin_edges = subtract(skin_edges,crack_surface_edges);
  for(Range::iterator sit = skin_edges.begin();sit!=skin_edges.end();sit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*sit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,coords,1,&coords[3],1);
    double nrm2 = cblas_dnrm2(3,&coords[3],1);	
    double d[2];
    rval = mField.get_moab().tag_get_data(th_distance,conn,num_nodes,d); CHKERR_PETSC(rval);
    double eps = 1e-2;
    if( ( fabs(d[0]/nrm2) < eps ) ) d[0] = 0;
    if( ( fabs(d[1]/nrm2) < eps ) ) d[1] = 0;
    if( (d[0]==0)||(d[1]==0) ) {
      if( d[0]==0 ) {
	rval = mField.get_moab().add_entities(nodesOnCrackSurface,&conn[0],1); CHKERR_PETSC(rval);
      }
      if( d[1]==0 ) {
	rval = mField.get_moab().add_entities(nodesOnCrackSurface,&conn[1],1); CHKERR_PETSC(rval);
      }
    } else {
      if(d[0]*d[1]<0) {
	rval = mField.get_moab().add_entities(opositeFrontEdges,&*sit,1); CHKERR_PETSC(rval);
      }
    }
  }
 
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getCrackSurfaceCorssingEdges(bool createMeshset) {
  PetscFunctionBegin;

  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,crackSurfaceCrossingEdges); CHKERR_PETSC(rval);
  }
  Range test_edges;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,test_edges,false); CHKERR_PETSC(rval);
  for(Range::iterator eit = test_edges.begin();eit!=test_edges.end();eit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,coords,1,&coords[3],1);
    double nrm2 = cblas_dnrm2(3,&coords[3],1);	
    cblas_dscal(3,1./nrm2,&coords[3],1);
    vector<EntityHandle> triangles_out;
    vector<double> distance_out;	
    const double tol = nrm2*1e-1;
    rval = kdTree_DistanceFromCrackSurface.ray_intersect_triangles(
      kdTree_rootMeshset_DistanceFromCrackSurface,tol,&coords[3],coords,
      triangles_out,distance_out); CHKERR_PETSC(rval);
    if(triangles_out.size()>1) {
      for(int nn = 0;nn<triangles_out.size();nn++) {
	cout << triangles_out[nn] << " " << distance_out[nn] << endl;
      }
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,&*eit,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,&triangles_out[0],triangles_out.size()); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("error.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      SETERRQ1(PETSC_COMM_SELF,1,"huston we have problem triangles_out.size() = %u",triangles_out.size());
    }
    if(!triangles_out.empty()) {
      ierr = mField.get_moab().add_entities(crackSurfaceCrossingEdges,&*eit,1); CHKERRQ(ierr);
    }
    /*cout << "Next:\n";
    ostream_iterator<double> out_it(cout, "\n");
    copy(distance_out.begin(),distance_out.end(), out_it);
    cout << "\n";*/
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::getCrackFrontTets(bool createMeshset) {
  PetscFunctionBegin;

  //get oposite crack foront faces, that is edges->nodes->tets
  Range opposite_crack_front_edges;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,opposite_crack_front_edges,false); CHKERR_PETSC(rval);
  opposite_crack_front_edges = intersect(opposite_crack_front_edges,mesh_level_edges);
  Range opposite_crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(opposite_crack_front_edges,opposite_crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range opposite_crack_front_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(opposite_crack_front_edges_nodes,3,false,opposite_crack_front_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  opposite_crack_front_edges_nodes_tets = intersect(opposite_crack_front_edges_nodes_tets,mesh_level_tets);

  //get crack front faces, that is edges->nodes->tets
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_front_edges = intersect(crack_front_edges,mesh_level_edges);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_edges_nodes_faces;
  rval = mField.get_moab().get_adjacencies(crack_front_edges_nodes,2,false,crack_front_edges_nodes_faces,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_front_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(crack_front_edges_nodes,3,false,crack_front_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_front_edges_nodes_tets = intersect(crack_front_edges_nodes_tets,mesh_level_tets);

  //common tets
  Range common_tets = intersect(crack_front_edges_nodes_tets,opposite_crack_front_edges_nodes_tets);

  //Nodes on crack surface
  Range nodes_on_crack_surface;
  rval = mField.get_moab().get_entities_by_type(nodesOnCrackSurface,MBVERTEX,nodes_on_crack_surface,false); CHKERR_PETSC(rval);
  Range nodes_on_crack_surface_tets;
  rval = mField.get_moab().get_adjacencies(nodes_on_crack_surface,3,false,nodes_on_crack_surface_tets,Interface::UNION); CHKERR_PETSC(rval);
  nodes_on_crack_surface_tets = intersect(nodes_on_crack_surface_tets,crack_front_edges_nodes_tets);
  common_tets.merge(nodes_on_crack_surface_tets);

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
  ierr = mField.get_moab().add_entities(crackFrontTests,nodes_on_crack_surface); CHKERRQ(ierr);
  ierr = mField.get_moab().add_entities(crackFrontTests,common_tets_faces_nodes); CHKERRQ(ierr);

  /*{
    EntityHandle out_meshset;
    rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    //rval = mField.get_moab().add_entities(out_meshset,common_tets); CHKERR_PETSC(rval);
    rval = mField.get_moab().add_entities(out_meshset,common_tets_faces); CHKERR_PETSC(rval);
    rval = mField.get_moab().write_file("debug.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
  }*/

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(bool createMeshset) {
  PetscFunctionBegin;

  Range crack_front_tets;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBTET,crack_front_tets,true); CHKERR_PETSC(rval);
  Range crack_front_tets_nodes;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBVERTEX,crack_front_tets_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_tets_faces;
  rval = mField.get_moab().get_entities_by_type(crackFrontTests,MBTRI,crack_front_tets_faces,true); CHKERR_PETSC(rval);

  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);

  if(createMeshset) {
    rval = mField.get_moab().create_meshset(MESHSET_SET,chopTetsFaces); CHKERR_PETSC(rval);
  }

  Range subtract_crack_front_tets_faces_chop_tets_faces;
  Range subtract_crack_front_tets_faces_chop_tets_faces_nodes;
  Range _crack_front_tets_faces_;
  Range _crack_front_tets_nodes_;

  int debug_ii = 0;

  do {

    cerr << "AAAAAAA " << crack_front_tets_nodes.size() << endl;

    double min_b;
    EntityHandle min_node = 0;
    for(Range::iterator nit = crack_front_tets_nodes.begin();nit!=crack_front_tets_nodes.end();nit++) {
      Range check_nodes;
      check_nodes.insert(*nit);
      check_nodes.merge(subtract_crack_front_tets_faces_chop_tets_faces_nodes);
      double b;
      ierr = calculate_qualityAfterProjectingNodes(check_nodes,b); CHKERRQ(ierr);
      if(nit == crack_front_tets_nodes.begin()) {
	min_b = b;
	min_node = *nit;
      } else {
	if(min_b > b) {
	  min_b = b;
	  min_node = *nit;
	}
      }
    }
    crack_front_tets_nodes.erase(min_node);

    Range chop_tets;
    rval = mField.get_moab().get_adjacencies(&min_node,1,3,false,chop_tets); CHKERR_PETSC(rval);
    chop_tets = intersect(chop_tets,crack_front_tets);
    //test if chop tets create single block
    if(chop_tets.size()>1) {
      Range seen_chop_tets;
      seen_chop_tets.insert(*chop_tets.begin());
      int nb_seen_chop_tets;
      do {
        nb_seen_chop_tets = seen_chop_tets.size();
	Range seen_chop_tets_faces;
	rval = mField.get_moab().get_adjacencies(seen_chop_tets,2,false,seen_chop_tets_faces,Interface::UNION); CHKERR_PETSC(rval);
	rval = mField.get_moab().get_adjacencies(seen_chop_tets_faces,3,false,seen_chop_tets,Interface::UNION); CHKERR_PETSC(rval);
	seen_chop_tets = intersect(seen_chop_tets,chop_tets);
      } while(nb_seen_chop_tets != seen_chop_tets.size());
      if(seen_chop_tets.size()!=chop_tets.size()) {
	//chop tets do not create single block
	continue;
      }
    }
    crack_front_tets = subtract(crack_front_tets,chop_tets);

    Range chop_faces;
    rval = mField.get_moab().get_adjacencies(&min_node,1,2,false,chop_faces,Interface::UNION); CHKERR_PETSC(rval);
    crack_front_tets_faces = subtract(crack_front_tets_faces,chop_faces);

    //get chop tets faces
    _crack_front_tets_faces_.clear();
    rval = mField.get_moab().get_adjacencies(crack_front_tets,2,false,_crack_front_tets_faces_,Interface::UNION); CHKERR_PETSC(rval);
    //from all faces subtract chop_tets_faces
    subtract_crack_front_tets_faces_chop_tets_faces = subtract(crack_front_tets_faces,_crack_front_tets_faces_);
    
    //get conncetivity of subtract_crack_front_tets_faces_chop_tets_faces
    rval = mField.get_moab().get_connectivity(
      subtract_crack_front_tets_faces_chop_tets_faces,subtract_crack_front_tets_faces_chop_tets_faces_nodes,true); CHKERR(rval);
    //get connectivity if crack_front_tets
    rval = mField.get_moab().get_connectivity(crack_front_tets,_crack_front_tets_nodes_,true); CHKERR_PETSC(rval);
    _crack_front_tets_nodes_ = intersect(_crack_front_tets_nodes_,crack_front_tets_nodes);
    crack_front_tets_nodes = subtract(_crack_front_tets_nodes_,subtract_crack_front_tets_faces_chop_tets_faces_nodes);

    subtract_crack_front_tets_faces_chop_tets_faces_nodes = subtract(subtract_crack_front_tets_faces_chop_tets_faces_nodes,crack_front_edges_nodes);
  
    debug_ii++;
    //if(debug_ii == 21) break; 

    {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,crack_front_tets); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,crack_front_tets_faces); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,subtract_crack_front_tets_faces_chop_tets_faces_nodes); CHKERR_PETSC(rval);
      ostringstream ss;
      ss << "chop_debug_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(out_meshset,chop_tets); CHKERR_PETSC(rval);
      ostringstream sss;
      sss << "chop_tets_debug_" << debug_ii << ".vtk";
      rval = mField.get_moab().write_file(sss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }	

  } while (!crack_front_tets.empty()); 


  rval = mField.get_moab().add_entities(chopTetsFaces,crack_front_tets); CHKERR_PETSC(rval);
  rval = mField.get_moab().add_entities(chopTetsFaces,crack_front_tets_faces); CHKERR_PETSC(rval);
  rval = mField.get_moab().add_entities(chopTetsFaces,subtract_crack_front_tets_faces_chop_tets_faces_nodes); CHKERR_PETSC(rval);

  //rval = mField.get_moab().add_entities(chopTetsFaces,subtract_crack_front_tets_faces_chop_tets_faces); CHKERR_PETSC(rval);
  //rval = mField.get_moab().add_entities(chopTetsFaces,_crack_front_tets_faces_); CHKERR_PETSC(rval);

  ierr = calculate_qualityAfterProjectingNodes(chopTetsFaces); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
/*PetscErrorCode FaceSplittingTools::getCrackFacesParts(EntityHandle meshset,EntityHandle meshset_edges) {
  PetscFunctionBegin;

  Range meshset_faces;
  rval = mField.get_moab().get_entities_by_type(meshset,MBTRI,meshset_faces,false); CHKERR_PETSC(rval);

  Range crack_front_edges,crack_front_edges_nodes;
  if(meshset_edges!=0) {
    rval = mField.get_moab().get_entities_by_type(meshset_edges,MBEDGE,crack_front_edges,true); CHKERR_PETSC(rval);
    rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes); CHKERR_PETSC(rval);
    //PetscPrintf(PETSC_COMM_WORLD,"crack_front_edges nb. edges = %u\n",crack_front_edges.size());
  }
  

  Range all_seen_faces;
  do {

    Range remaining_faces 
      = subtract(meshset_faces,all_seen_faces);
 
    unsigned int nb_option_faces;
    Range option_faces;
    option_faces.insert(*remaining_faces.begin());

    do {

      nb_option_faces = option_faces.size();

      Range option_faces_edges;
      rval = mField.get_moab().get_adjacencies(
	option_faces,1,false,option_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
      if(meshset_edges!=0) {
	option_faces_edges = subtract(option_faces_edges,crack_front_edges);
      }
      Range adj_faces;
      rval = mField.get_moab().get_adjacencies(
	option_faces_edges,2,false,adj_faces,Interface::UNION); CHKERR_PETSC(rval);
      adj_faces = intersect(adj_faces,meshset_faces);
      option_faces.merge(intersect(adj_faces,remaining_faces));

    } while (nb_option_faces != option_faces.size() );

    EntityHandle child;
    rval = mField.get_moab().create_meshset(MESHSET_SET,child); CHKERR_PETSC(rval);
    rval = mField.get_moab().add_child_meshset(meshset,child); CHKERR_PETSC(rval);
    rval = mField.get_moab().add_entities(child,option_faces); CHKERR_PETSC(rval);

    Range option_nodes;
    rval = mField.get_moab().get_connectivity(option_faces,option_nodes); CHKERR_PETSC(rval);
    option_nodes = subtract(option_nodes,crack_front_edges_nodes);
    rval = mField.get_moab().add_entities(child,option_nodes); CHKERR_PETSC(rval);

    all_seen_faces.merge(option_faces);

  } while (all_seen_faces.size()!=meshset_faces.size());

  vector<EntityHandle> children;
  rval = mField.get_moab().get_child_meshsets(meshset,children); CHKERR_PETSC(rval);
  PetscPrintf(PETSC_COMM_WORLD,"nb. of part = %u\n",children.size());

  PetscFunctionReturn(0);
}*/
/*PetscErrorCode FaceSplittingTools::getCrackFacesParts(const BitRefLevel bit_mesh) {
  PetscFunctionBegin;

  Range mesh_level_tets;
  Range mesh_level_tris;
  if(bit_mesh.any()) {
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr); 
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr); 
  }

  ierr = getCrackFacesParts(newFrontFaces,0); CHKERRQ(ierr);

  EntityHandle meshset_edges;
  ierr = mField.get_Cubit_msId_meshset(201,SideSet,meshset_edges); CHKERRQ(ierr);
  Range crack_front_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_front_edges,true); CHKERRQ(ierr);

  Skinner skin(&mField.get_moab());
  Range crack_front_nodes;
  rval = mField.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  Range crack_front_edges_faces;
  rval = mField.get_moab().get_adjacencies(crack_front_nodes,2,false,crack_front_edges_faces,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_front_tets;
  rval = mField.get_moab().get_adjacencies(crack_front_nodes,3,false,crack_front_tets,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_front_tets_skin_faces; 
  rval = skin.find_skin(crack_front_tets,false,crack_front_tets_skin_faces); CHKERR(rval);

  Range oposite_front_edge;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,oposite_front_edge,false); CHKERR_PETSC(rval);
  Range oposite_front_edge_nodes;
  rval = mField.get_moab().get_connectivity(oposite_front_edge,oposite_front_edge_nodes,true); CHKERR_PETSC(rval);

  vector<EntityHandle> children;
  rval = mField.get_moab().get_child_meshsets(newFrontFaces,children); CHKERR_PETSC(rval);
  for(unsigned int ii = 0;ii<children.size();ii++) {

    ierr = getCrackFacesParts(children[ii],meshset_edges); CHKERRQ(ierr);
    vector<EntityHandle> children_of_children;
    rval = mField.get_moab().get_child_meshsets(children[ii],children_of_children); CHKERR_PETSC(rval);
    if(children_of_children.size()!=2) {
      SETERRQ(PETSC_COMM_SELF,1,"should be two children meshsets");
    }

    //get tets iin between options
    Range option_faces_tets;
    vector<Range> option_faces;
    option_faces.resize(children_of_children.size());
    for(int jj = 0;jj<children_of_children.size();jj++) {

      //get option faces
      rval = mField.get_moab().get_entities_by_type(children_of_children[jj],MBTRI,option_faces[jj],false); CHKERR_PETSC(rval);
      //get option faces edges
      Range option_faces_edges;
      rval = mField.get_moab().get_adjacencies(option_faces[jj],1,false,option_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
      option_faces_edges = subtract(option_faces_edges,crack_front_edges);
      //get opetion faces edges test
      Range option_faces_edges_tets;
      rval = mField.get_moab().get_adjacencies(option_faces_edges,3,false,option_faces_edges_tets,Interface::UNION); CHKERR_PETSC(rval);
      option_faces_edges_tets = intersect(option_faces_edges_tets,mesh_level_tets);

      if(jj > 0) {
	option_faces_tets = intersect(option_faces_tets,option_faces_edges_tets);
      } else {
	option_faces_tets = option_faces_edges_tets;
      }
     
    }

    Range option_faces_tets_faces;
    rval = mField.get_moab().get_adjacencies(option_faces_tets,2,false,option_faces_tets_faces,Interface::UNION); CHKERR_PETSC(rval);
    vector<Range> chop_tets(children_of_children.size()),cat_out_faces(children_of_children.size());
    for(int jj = 0;jj<children_of_children.size();jj++) {
      //calulate quality 
      ierr = calculate_qualityAfterProjectingNodes(children_of_children[jj]); CHKERRQ(ierr);
      //get option nodes
      Range option_nodes;
      rval = mField.get_moab().get_entities_by_type(children_of_children[jj],MBVERTEX,option_nodes,false); CHKERR_PETSC(rval);
      //get tets adjacent to poor quality nodes
      for(Range::iterator nit = option_nodes.begin();nit!=option_nodes.end();nit++) {
	double b;
	rval = mField.get_moab().tag_get_data(th_b,&*nit,1,&b); CHKERR_PETSC(rval);
	if(b <= 0) {
	  Range tets;
	  rval = mField.get_moab().get_adjacencies(&*nit,1,3,false,tets); CHKERR_PETSC(rval);
	  if(bit_mesh.any()) {
	    tets = intersect(tets,mesh_level_tets);
	  }
	  chop_tets[jj].merge(tets);
	}
      }
      Range faces;
      rval = skin.find_skin(chop_tets[jj],false,faces); CHKERR(rval);

      cat_out_faces[jj] = intersect(faces,option_faces_tets_faces);
      cat_out_faces[jj] = subtract(cat_out_faces[jj],crack_front_tets_skin_faces);

      Range unattached_nodes;
      mField.get_moab().get_connectivity(cat_out_faces[jj],unattached_nodes,true);
      unattached_nodes = subtract(unattached_nodes,crack_front_nodes);
      unattached_nodes = subtract(unattached_nodes,oposite_front_edge_nodes);
      Range unattached_faces;
      rval = mField.get_moab().get_adjacencies(unattached_nodes,2,false,unattached_faces,Interface::UNION); CHKERR_PETSC(rval);
      cat_out_faces[jj] = subtract(cat_out_faces[jj],unattached_faces);

      //add choping faces tets
      option_faces_tets.merge(chop_tets[jj]);
    }

    Range skin_faces; 
    rval = skin.find_skin(option_faces_tets,false,skin_faces); CHKERR(rval);
    
    for(int jj = 0;jj<children_of_children.size();jj++) {
      
      Range _option_faces = intersect(option_faces[jj],skin_faces);
      _option_faces = subtract(option_faces[jj],_option_faces);
      rval = mField.get_moab().remove_entities(children_of_children[jj],_option_faces); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(children_of_children[jj],cat_out_faces[jj]); CHKERR_PETSC(rval);

      Range nodes;
      rval = mField.get_moab().get_entities_by_type(children_of_children[jj],MBVERTEX,nodes,false); CHKERR_PETSC(rval);
      rval = mField.get_moab().remove_entities(children_of_children[jj],nodes); CHKERR_PETSC(rval);

      Range faces;
      rval = mField.get_moab().get_entities_by_type(children_of_children[jj],MBTRI,faces,false); CHKERR_PETSC(rval);

      nodes.clear();
      rval = mField.get_moab().get_connectivity(faces,nodes,true); CHKERR_PETSC(rval);
      nodes = subtract(nodes,crack_front_nodes);
      rval = mField.get_moab().add_entities(children_of_children[jj],nodes); CHKERR_PETSC(rval);
     
    }


  }

  PetscFunctionReturn(0);
}*/
/*PetscErrorCode FaceSplittingTools::projectNodesAndEvaluateElementQuality() {
  PetscFunctionBegin;

  vector<EntityHandle> children;
  rval = mField.get_moab().get_child_meshsets(newFrontFaces,children); CHKERR_PETSC(rval);

  for(unsigned int ii = 0;ii<children.size();ii++) {

    Range part_faces;
    rval = mField.get_moab().get_entities_by_type(children[ii],MBTRI,part_faces,false); CHKERR_PETSC(rval);
    vector<EntityHandle> children_of_children;
    rval = mField.get_moab().get_child_meshsets(children[ii],children_of_children); CHKERR_PETSC(rval);

    for(unsigned int jj = 0;jj<children_of_children.size();jj++) {

      PetscPrintf(PETSC_COMM_WORLD,"AAA %u %u\n",ii,jj);
      ierr = calculate_qualityAfterProjectingNodes(children_of_children[jj]); CHKERRQ(ierr);

    }

  }
  
  PetscFunctionReturn(0);
}*/
PetscErrorCode FaceSplittingTools::catMesh() {
  PetscFunctionBegin;

  Range edges_to_refine;
  rval = mField.get_moab().get_entities_by_type(opositeFrontEdges,MBEDGE,edges_to_refine,false); CHKERR_PETSC(rval);

  int current_ref_bit = meshRefineBitLevels.back();
  BitRefLevel current_ref = BitRefLevel().set(current_ref_bit);
  Range level_edges;
  ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
  edges_to_refine = intersect(edges_to_refine,mesh_level_edges);
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
PetscErrorCode FaceSplittingTools::meshRefine() {
  PetscFunctionBegin;

  PetscBool flg = PETSC_TRUE;
  PetscInt nb_ref_levels;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_ref",&nb_ref_levels,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    nb_ref_levels = 0;
  }

  for(int ll = 1;ll<nb_ref_levels+1;ll++) {
 
      int current_ref_bit = meshRefineBitLevels.back();
      BitRefLevel current_ref = BitRefLevel().set(current_ref_bit);
  
      Range level_tets;
      ierr = mField.get_entities_by_ref_level(current_ref,BitRefLevel().set(),level_tets); CHKERRQ(ierr);

      Range crack_edges;
      ierr = mField.get_Cubit_msId_entities_by_dimension(201,SideSet,1,crack_edges,true); CHKERRQ(ierr);

      Range crack_edges_nodes;
      rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
      Range crack_edges_nodes_tets;
      rval = mField.get_moab().get_adjacencies(crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
      crack_edges_nodes_tets = intersect(level_tets.subset_by_type(MBTET),crack_edges_nodes_tets);

      Range edges_to_refine;
      rval = mField.get_moab().get_adjacencies(crack_edges_nodes_tets,1,false,edges_to_refine,Interface::UNION); CHKERR_PETSC(rval);
      int last_ref_bit = meshRefineBitLevels.back();
      if(!meshIntefaceBitLevels.empty()) {
	if(last_ref_bit<meshIntefaceBitLevels.back()) {
	  last_ref_bit = meshIntefaceBitLevels.back();
	}
      }
      last_ref_bit++;
      meshRefineBitLevels.push_back(last_ref_bit);     
      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
      cout << "AAAAAAA " << last_ref << endl;
      ierr = mField.add_verices_in_the_middel_of_edges(edges_to_refine,last_ref,2); CHKERRQ(ierr);
      ierr = mField.refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);
      
      for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
	EntityHandle cubit_meshset = cubit_it->meshset; 
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
	ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
      }
  
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::splitFaces() {
  PetscFunctionBegin;

  BitRefLevel current_ref = BitRefLevel().set(meshRefineBitLevels.back());
  cout << current_ref << endl;

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

    ierr = mField.get_msId_3dENTS_split_sides(bit_meshset,last_ref,meshset_interface,false,true); CHKERRQ(ierr);

    //add refined ent to cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset; 
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }
  }
  rval = mField.get_moab().delete_entities(&bit_meshset,1); CHKERR_PETSC(rval);

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::addNewSurfaceFaces_to_Cubit_msId200() {
  PetscFunctionBegin;
  EntityHandle meshset200;
  ierr = mField.get_Cubit_msId_meshset(200,SideSet,meshset200); CHKERRQ(ierr);
  Range new_crack_surface_faces;
  rval = mField.get_moab().get_entities_by_type(chopTetsFaces,MBTRI,new_crack_surface_faces,false); CHKERR_PETSC(rval);
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
PetscErrorCode FaceSplittingTools::calculate_qualityAfterProjectingNodes(EntityHandle meshset) {
  PetscFunctionBegin;
  Range option_nodes;
  rval = mField.get_moab().get_entities_by_type(meshset,MBVERTEX,option_nodes,false); CHKERR_PETSC(rval);
  double current_b;
  ierr = calculate_qualityAfterProjectingNodes(option_nodes,current_b); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"meshset b = %3.2f\n",current_b);
  rval = mField.get_moab().tag_set_data(th_b,&meshset,1,&current_b); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::calculate_qualityAfterProjectingNodes(Range &option_nodes,double &current_b) {
  PetscFunctionBegin;
  double def_VAL[1] = { 0 };
  rval = mField.get_moab().tag_get_handle("MESHSET_PROJECTION_B",1,MB_TYPE_DOUBLE,th_b,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
  ierr = ShapeDiffMBTET(diffNTET); CHKERRQ(ierr);
  map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,3> > > projectedNodes;
  for(Range::iterator nit = option_nodes.begin();nit!=option_nodes.end();nit++) {
    double coords[3]; 
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
      ublas::vector<double,ublas::bounded_array<double,3> > 
	&closest_point_out = projectedNodes[*nit];
      closest_point_out.resize(3);
      EntityHandle triangle_out;
      rval = kdTree_DistanceFromCrackSurface.closest_triangle(
      kdTree_rootMeshset_DistanceFromCrackSurface,
	coords,&*closest_point_out.data().begin(),triangle_out); CHKERR_PETSC(rval);
  }
  map<EntityHandle,double> b_map;  
  Range adj_tets;
  rval = mField.get_moab().get_adjacencies(option_nodes,3,false,adj_tets,Interface::UNION); CHKERR_PETSC(rval);
  for(Range::iterator tit = adj_tets.begin();tit!=adj_tets.end();tit++) {
    int num_nodes; 
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes]; 
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    double dofs_X[3*num_nodes];
    cblas_dcopy(3*num_nodes,coords,1,dofs_X,1);
    for(int nn = 0;nn<num_nodes;nn++) {
      map<EntityHandle,ublas::vector<double,ublas::bounded_array<double,3> > >::iterator 
	mit = projectedNodes.find(conn[nn]);
      if(mit != projectedNodes.end() ) {
	cblas_dcopy(3,&*mit->second.data().begin(),1,&dofs_X[3*nn],1);
      }
    }
    double coords_edges[2*3*6]; 
    ierr = get_edges_from_elem_coords(coords,coords_edges); CHKERRQ(ierr);
    double V =  Shape_intVolumeMBTET(diffNTET,&*coords); 
    double alpha[4] = {1,1,1,1};
    double quality0,quality,b;
    ierr = quality_volume_length_F(
	  V,alpha,0,
	  diffNTET,coords_edges,dofs_X,
	  NULL,NULL,NULL,
	  &quality0,&quality,&b,
	  NULL,NULL); CHKERRQ(ierr);
    Range new_front_nodes;
    new_front_nodes.insert(&conn[0],&conn[4]);
    new_front_nodes = intersect(new_front_nodes,option_nodes);
    for(Range::iterator nit = new_front_nodes.begin();nit!=new_front_nodes.end();nit++) {
      map<EntityHandle,double>::iterator mit = b_map.find(*nit);
      if(mit!=b_map.end()) {
	mit->second = fmin(mit->second,b);
      } else {
	b_map[*nit] = b;
      }
    }
    if(tit != adj_tets.begin()) {
      b = fmin(b,current_b);
    }
    current_b = b;
  }
  for(map<EntityHandle,double>::iterator mit = b_map.begin();mit!=b_map.end();mit++) {
    //PetscPrintf(PETSC_COMM_WORLD,"node %u b = %3.2f\n",mit->first,mit->second);
    rval = mField.get_moab().tag_set_data(th_b,&mit->first,1,&mit->second); CHKERR_PETSC(rval);
  }
  PetscFunctionReturn(0);
}

}

