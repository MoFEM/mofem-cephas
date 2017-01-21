/** \file CutMeshInterface.cpp
 * \brief Cut mesh by surface
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#include <version.h>
#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <NodeMerger.hpp>
#include <Core.hpp>

#include <CutMeshInterface.hpp>
#include <TetGenInterface.hpp>

namespace MoFEM {

  PetscErrorCode CutMeshInterface::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {
    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMCutMesh) {
      *iface = dynamic_cast<CutMeshInterface*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  CutMeshInterface::CutMeshInterface(const MoFEM::Core &core):
  cOre(const_cast<MoFEM::Core&>(core)) {
  }

  PetscErrorCode CutMeshInterface::setSurface(const Range &surface) {
    PetscFunctionBegin;
    sUrface = surface;
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::setVolume(const Range &volume) {
    PetscFunctionBegin;
    vOlume = volume;
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::mergeSurface(const Range &surface)  {
    PetscFunctionBegin;
    sUrface.merge(surface);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::mergeVolume(const Range &volume)  {
    PetscFunctionBegin;
    vOlume.merge(volume);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::buildTree() {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    treeSurfPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab,"ROOTSETSURF",true)
    );
    rval = treeSurfPtr->build(sUrface,rootSetSurf); CHKERRQ_MOAB(rval);
    Range faces;
    rval = moab.get_adjacencies(vOlume,2,false,faces,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::findEdgesToCut(int verb) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    // Range vol_edges;
    // rval = moab.get_adjacencies(
    //   vOlume,1,true,vol_edges,moab::Interface::UNION
    // ); CHKERRQ_MOAB(rval);
    edgesToCut.clear();
    cutEdges.clear();
    double ray_length;
    double ray_point[3],unit_ray_dir[3];
    VectorAdaptor vec_unit_ray_dir(
      3,ublas::shallow_array_adaptor<double>(3,unit_ray_dir)
    );
    VectorAdaptor vec_ray_point(
      3,ublas::shallow_array_adaptor<double>(3,ray_point)
    );

    Tag th_dist;
    double def_val[] = {0};
    rval = moab.tag_get_handle(
      "DIST",1,MB_TYPE_DOUBLE,th_dist,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val
    ); CHKERRQ_MOAB(rval);
    Tag th_dist_normal;
    rval = moab.tag_get_handle(
      "DIST_NORMAL",1,MB_TYPE_DOUBLE,th_dist_normal,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val
    ); CHKERRQ_MOAB(rval);

    Range vol_vertices;
    rval = moab.get_connectivity(vOlume,vol_vertices,true); CHKERRQ_MOAB(rval);
    for(Range::iterator vit = vol_vertices.begin();vit!=vol_vertices.end();vit++) {
      double coords[3];
      rval = moab.get_coords(&*vit,1,coords); CHKERR_MOAB(rval);
      VectorAdaptor point_in(3,ublas::shallow_array_adaptor<double>(3,coords));
      double p_out[3];
      EntityHandle facets_out;
      rval = treeSurfPtr->closest_to_location(&coords[0],rootSetSurf,p_out,facets_out); CHKERR_MOAB(rval);
      VectorAdaptor point_out(3,ublas::shallow_array_adaptor<double>(3,p_out));
      double normal[3];
      Util::normal(&moab,facets_out,normal[0],normal[1],normal[2]);
      VectorAdaptor n(3,ublas::shallow_array_adaptor<double>(3,normal));
      VectorDouble3 delta = point_out-point_in;
      double dist = norm_2(delta);
      double dist_normal = inner_prod(delta,n)/norm_2(n);
      rval = moab.tag_set_data(th_dist,&*vit,1,&dist); CHKERR_MOAB(rval);
      rval = moab.tag_set_data(th_dist_normal,&*vit,1,&dist_normal); CHKERRQ_MOAB(rval);
    }

    const double low_tol = 0;
    Range vol_edges;
    rval = moab.get_adjacencies(
      vOlume,1,true,vol_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    for(Range::iterator eit = vol_edges.begin();eit!=vol_edges.end(); eit++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double dist[num_nodes];
      rval = moab.tag_get_data(th_dist,conn,num_nodes,dist); CHKERRQ_MOAB(rval);
      double dist_normal[num_nodes];
      rval = moab.tag_get_data(th_dist_normal,conn,num_nodes,dist_normal); CHKERRQ_MOAB(rval);
      if(dist_normal[0]*dist_normal[1]<0) {
        ierr = getRayForEdge(
          *eit,vec_ray_point,vec_unit_ray_dir,ray_length
        ); CHKERRQ(ierr);
        // double s = fabs(dist_normal[0])/(fabs(dist_normal[0])+fabs(dist_normal[1]));
        // edgesToCut[*eit].dIst = s*ray_length;
        // edgesToCut[*eit].lEngth = ray_length;
        // edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
        // edgesToCut[*eit].rayPoint = vec_ray_point;
        // cutEdges.insert(*eit);
        std::vector< double > distances_out;
        std::vector< EntityHandle > facets_out;
        const double tol = ray_length*low_tol;
        rval = treeSurfPtr->ray_intersect_triangles(
          distances_out,facets_out,rootSetSurf,tol,ray_point,unit_ray_dir,&ray_length
        ); CHKERR_MOAB(rval);
        if(!distances_out.empty()) {
          if(distances_out[0]/ray_length < low_tol) {
            int num_nodes;
            const EntityHandle *conn;
            rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
            verticesOnSurface.insert(conn[0]);
          } else if(distances_out[0]/ray_length > 1-low_tol) {
            int num_nodes;
            const EntityHandle *conn;
            rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
            verticesOnSurface.insert(conn[1]);
          } else {
            edgesToCut[*eit].dIst = distances_out[0];
            edgesToCut[*eit].lEngth = ray_length;
            edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
            edgesToCut[*eit].rayPoint = vec_ray_point;
            cutEdges.insert(*eit);
          }
        }
      }
    }

    cutVolumes.clear();
    rval = moab.get_adjacencies(
      cutEdges,3,false,cutVolumes,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);

    Range edges;
    rval = moab.get_adjacencies(
      cutVolumes,1,false,edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    edges = subtract(edges,cutEdges);

    for(Range::iterator eit = edges.begin();eit!=edges.end();eit++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double dist[num_nodes];
      rval = moab.tag_get_data(th_dist,conn,num_nodes,dist); CHKERRQ_MOAB(rval);
      double dist_normal[num_nodes];
      rval = moab.tag_get_data(th_dist_normal,conn,num_nodes,dist_normal); CHKERRQ_MOAB(rval);
      if(dist_normal[0]*dist_normal[1]<0) {
        ierr = getRayForEdge(
          *eit,vec_ray_point,vec_unit_ray_dir,ray_length
        ); CHKERRQ(ierr);
        double s = fabs(dist_normal[0])/(fabs(dist_normal[0])+fabs(dist_normal[1]));
        edgesToCut[*eit].dIst = s*ray_length;
        edgesToCut[*eit].lEngth = ray_length;
        edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
        edgesToCut[*eit].rayPoint = vec_ray_point;
        cutEdges.insert(*eit);
      }
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::cutEdgesInMiddle(
    const BitRefLevel bit
  ) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    MeshRefinement *refiner;
    const RefEntity_multiIndex *ref_ents_ptr;
    PetscFunctionBegin;
    ierr = m_field.query_interface(refiner); CHKERRQ(ierr);
    ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    ierr = refiner->add_verices_in_the_middel_of_edges(cutEdges,bit); CHKERRQ(ierr);
    ierr = refiner->refine_TET(vOlume,bit,false); CHKERRQ(ierr);
    cutNewVolumes.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,BitRefLevel().set(),MBTET,cutNewVolumes
    ); CHKERRQ(ierr);
    cutNewSurfaces.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,bit,MBTRI,cutNewSurfaces
    ); CHKERRQ(ierr);
    // Find new vertices on catted edges
    verticecOnCutEdges.clear();
    cutNewVertices.clear();
    // Range opposite_edges_nodes;
    for(
      map<EntityHandle,TreeData>::iterator mit = edgesToCut.begin();
      mit!=edgesToCut.end();mit++
    ) {
      boost::shared_ptr<RefEntity> ref_ent = *(
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>
        ().find(boost::make_tuple(mit->first,MBVERTEX))
      );
      if((ref_ent->getBitRefLevel()&bit).any()) {
        EntityHandle vert = ref_ent->getRefEnt();
        cutNewVertices.insert(vert);
        verticecOnCutEdges[vert] = mit->second;
        // if(cutEdgesOutside.find(mit->first)!=cutEdgesOutside.end()) {
        //   opposite_edges_nodes.insert(vert);
        // }
      }
    }
    // At that point cutNewSurfaces has all newly created faces, not take all nodes
    // on those faces and subtract nodes on catted edges. Faces adjacent to nodes
    // which left are not part of surface.
    Range diff_verts;
    rval = moab.get_connectivity(cutNewSurfaces,diff_verts,true); CHKERRQ_MOAB(rval);
    diff_verts = subtract(diff_verts,unite(cutNewVertices,verticesOnSurface));
    // diff_verts = unite(diff_verts,opposite_edges_nodes);
    Range subtract_faces;
    rval = moab.get_adjacencies(
      diff_verts,2,false,subtract_faces,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    cutNewSurfaces = subtract(cutNewSurfaces,subtract_faces);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::trimEdgesInTheMiddle(const BitRefLevel bit) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    MeshRefinement *refiner;
    const RefEntity_multiIndex *ref_ents_ptr;
    PetscFunctionBegin;

    ierr = m_field.query_interface(refiner); CHKERRQ(ierr);
    ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    ierr = refiner->add_verices_in_the_middel_of_edges(trimEdges,bit); CHKERRQ(ierr);
    ierr = refiner->refine_TET(cutNewVolumes,bit,false); CHKERRQ(ierr);
    trimNewVolumes.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,bit,MBTET,trimNewVolumes
    ); CHKERRQ(ierr);
    // Get vertices which are on trim edges
    verticecOnTrimEdges.clear();
    trimNewVertices.clear();
    for(
      map<EntityHandle,TreeData>::iterator mit = edgesToTrim.begin();
      mit!=edgesToTrim.end();mit++
    ) {
      boost::shared_ptr<RefEntity> ref_ent = *(
        ref_ents_ptr->get<Composite_ParentEnt_And_EntType_mi_tag>
        ().find(boost::make_tuple(mit->first,MBVERTEX))
      );
      if((ref_ent->getBitRefLevel()&bit).any()) {
        EntityHandle vert = ref_ent->getRefEnt();
        trimNewVertices.insert(vert);
        verticecOnTrimEdges[vert] = mit->second;
      }
    }

    // Get faces which are trimmed
    trimNewSurfaces.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,bit,MBTRI,trimNewSurfaces
    ); CHKERRQ(ierr);
    Range trim_new_surfaces_nodes;
    rval = moab.get_connectivity(trimNewSurfaces,trim_new_surfaces_nodes,true); CHKERRQ_MOAB(rval);
    trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes,trimNewVertices);
    trim_new_surfaces_nodes = subtract(trim_new_surfaces_nodes,cutNewVertices);
    Range outside_edges_vertices;
    rval = moab.get_connectivity(outsideEdges,outside_edges_vertices,true); CHKERRQ_MOAB(rval);
    trim_new_surfaces_nodes = unite(trim_new_surfaces_nodes,outside_edges_vertices);
    Range faces_not_on_surface;
    rval = moab.get_adjacencies(trim_new_surfaces_nodes,2,false,faces_not_on_surface,moab::Interface::UNION);
    trimNewSurfaces = subtract(trimNewSurfaces,faces_not_on_surface);

    // Get surfaces which are not trimmed and add them to surface
    Range all_surfaces_on_bit_level;
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,BitRefLevel().set(),MBTRI,all_surfaces_on_bit_level
    ); CHKERRQ(ierr);
    all_surfaces_on_bit_level = intersect(all_surfaces_on_bit_level,cutNewSurfaces);
    trimNewSurfaces = unite(trimNewSurfaces,all_surfaces_on_bit_level);


    PetscFunctionReturn(0);

  }

  PetscErrorCode CutMeshInterface::moveMidNodesOnCutEdges() {
    PetscFunctionBegin;
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    // Range out_side_vertices;
    for(
      map<EntityHandle,TreeData>::iterator mit = verticecOnCutEdges.begin();
      mit!=verticecOnCutEdges.end(); mit++
    ) {
      double dist = mit->second.dIst;
      // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth << endl;
      VectorDouble3 new_coors = mit->second.rayPoint+dist*mit->second.unitRayDir;
      rval = moab.set_coords(&mit->first,1,&new_coors[0]); CHKERRQ_MOAB(rval);
    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode CutMeshInterface::moveMidNodesOnTrimedEdges() {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    // Range out_side_vertices;
    for(
      map<EntityHandle,TreeData>::iterator mit = verticecOnTrimEdges.begin();
      mit!=verticecOnTrimEdges.end(); mit++
    ) {
      double dist = mit->second.dIst;
      // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth << endl;
      VectorDouble3 new_coors = mit->second.rayPoint+dist*mit->second.unitRayDir;
      rval = moab.set_coords(&mit->first,1,&new_coors[0]); CHKERRQ_MOAB(rval);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::findEdgesToTrim(int verb) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    Range edges;
    rval = moab.get_adjacencies(
      cutNewSurfaces,1,false,edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);

    trimEdges.clear();
    outsideEdges.clear();
    edgesToTrim.clear();

    for(Range::iterator eit = edges.begin();eit!=edges.end();eit++) {
      // Get edge connectivity and coords
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double coords[3*num_nodes];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERR_MOAB(rval);
      // Put edges coords into boost vectors
      VectorAdaptor s0(3,ublas::shallow_array_adaptor<double>(3,&coords[0]));
      VectorAdaptor s1(3,ublas::shallow_array_adaptor<double>(3,&coords[3]));
      // Find point on surface closet to surface
      double point_out0[3];
      EntityHandle facets_out0;
      rval = treeSurfPtr->closest_to_location(
        &coords[0],rootSetSurf,point_out0,facets_out0
      ); CHKERR_MOAB(rval);
      double point_out1[3];
      EntityHandle facets_out1;
      rval = treeSurfPtr->closest_to_location(
        &coords[3],rootSetSurf,point_out1,facets_out1
      ); CHKERR_MOAB(rval);
      // Put closest point in boost vectors
      VectorAdaptor p0(3,ublas::shallow_array_adaptor<double>(3,point_out0));
      VectorAdaptor p1(3,ublas::shallow_array_adaptor<double>(3,point_out1));
      // Calculate deltas, i.e. vectors from edges to closet point on surface
      VectorDouble3 delta0,delta1;
      delta0 = p0-s0;
      delta1 = p1-s1;
      // moab.tag_set_data(th,&conn[0],1,&delta0[0]);
      // moab.tag_set_data(th,&conn[1],1,&delta1[0]);
      // Calculate distances
      double dist0 = norm_2(delta0);
      double dist1 = norm_2(delta1);
      double min_dist = fmin(dist0,dist1);
      double max_dist = fmax(dist0,dist1);
      // If one of nodes is on the surface and other is not, that edge is to trim
      const double tol = 1e-4;
      if(min_dist < tol && max_dist > tol) {
        trimEdges.insert(*eit);
        if(max_dist==dist0) {
          VectorDouble3 ray = p0-s0;
          double ray_length = norm_2(ray);
          edgesToTrim[*eit].dIst = dist0;
          edgesToTrim[*eit].lEngth = ray_length;
          edgesToTrim[*eit].unitRayDir = ray/ray_length;
          edgesToTrim[*eit].rayPoint = s0;
        } else {
          VectorDouble3 ray = p1-s1;
          double ray_length = norm_2(ray);
          edgesToTrim[*eit].dIst = dist1;
          edgesToTrim[*eit].lEngth = ray_length;
          edgesToTrim[*eit].unitRayDir = ray/ray_length;
          edgesToTrim[*eit].rayPoint = s1;
        }
      }
      if(min_dist > tol) {
        outsideEdges.insert(*eit);
      }
    }

    // Remove tris which are outside surface
    Range outside_edges_tris;
    rval = moab.get_adjacencies(
      outsideEdges,2,false,outside_edges_tris,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    for(Range::iterator tit = outside_edges_tris.begin();tit!=outside_edges_tris.end();tit++) {
      Range tit_edges;
      rval = moab.get_adjacencies(
        &*tit,1,1,false,tit_edges,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      if(intersect(tit_edges,outsideEdges).size()==3) {
        cutNewSurfaces.erase(*tit);
      }
    }
    cutNewVertices.clear();
    rval = moab.get_connectivity(cutNewSurfaces,cutNewVertices,true); CHKERRQ_MOAB(rval);


    // Remove form outside edges, those which are not part of cutNewSurfaces edges
    Range cut_new_surface_edges;
    rval = moab.get_adjacencies(
      cutNewSurfaces,1,false,cut_new_surface_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    outsideEdges = intersect(outsideEdges,cut_new_surface_edges);

    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::nodeMergeCutSurface(const BitRefLevel bit) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    // const RefEntity_multiIndex *ref_ents_ptr;
    Skinner skin(&m_field.get_moab());
    NodeMergerInterface *node_merger_ptr;
    PetscFunctionBegin;

    // Get access to ref entities multi-indices
    // ierr = m_field.get_ref_ents(&ref_ents_ptr); CHKERRQ(ierr);
    // Get interface to node merger
    ierr = m_field.query_interface(node_merger_ptr); CHKERRQ(ierr);

    // Get nodes on surface and surface edges
    Range surf_nodes;
    rval = moab.get_connectivity(cutNewSurfaces,surf_nodes,true);
    Range surf_skin;
    rval = skin.find_skin(0,cutNewSurfaces,false,surf_skin); CHKERR_MOAB(rval);
    Range surf_skin_nodes;
    rval = moab.get_connectivity(surf_skin,surf_skin_nodes,true);

    // Get tets adjacent to surface nodes
    Range adj_tets;
    rval = moab.get_adjacencies(
      surf_nodes,3,false,adj_tets,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    adj_tets = intersect(adj_tets,cutNewVolumes);
    Range adj_tets_edges;
    rval = moab.get_adjacencies(
      adj_tets,1,false,adj_tets_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    Range skin_adj_tets;
    rval = skin.find_skin(0,adj_tets,false,skin_adj_tets); CHKERR_MOAB(rval);

    Range skin_adj_tets_edges;
    rval = moab.get_adjacencies(skin_adj_tets,1,false,skin_adj_tets_edges,moab::Interface::UNION);
    adj_tets_edges = subtract(adj_tets_edges,skin_adj_tets_edges);

    map<double,EntityHandle> map_edges;
    for(
      Range::iterator
      eit = adj_tets_edges.begin();eit!=adj_tets_edges.end();eit++
     ) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double coords[3*num_nodes];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERR_MOAB(rval);
      VectorAdaptor s0(3,ublas::shallow_array_adaptor<double>(3,&coords[0]));
      VectorAdaptor s1(3,ublas::shallow_array_adaptor<double>(3,&coords[3]));
      map_edges[norm_2(s0-s1)] = *eit;
    }

    Range merged_edges;
    Range in_tets = adj_tets;
    for(map<double,EntityHandle>::iterator mit = map_edges.begin();mit!=map_edges.end();mit++) {

      if(merged_edges.find(mit->second)!=merged_edges.end()) {
        continue;
      }

      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(mit->second,conn,num_nodes,true); CHKERRQ_MOAB(rval);

      int fixed[] = {0,0};
      for(int nn = 0;nn!=2;nn++) {
        if(surf_skin_nodes.find(conn[nn])!=surf_skin_nodes.end()) {
          fixed[nn] = 2;
        } else if(surf_nodes.find(conn[nn])!=surf_nodes.end()) {
          fixed[nn] = 1;
        }
      }
      if(
        fixed[0] == 2 && fixed[1] == 2 &&
        surf_skin_nodes.find(mit->second)==surf_skin_nodes.end()
      ) {
        continue;
      }
      if(
        fixed[0] == 1 && fixed[1] == 1 &&
        surf_nodes.find(mit->second)!=surf_nodes.end()
      ) {
        continue;
      }


      EntityHandle father,mother;

      double move = -1;
      if(fixed[0] == fixed[1]) {
        father = conn[0];
        mother = conn[1];
        move = 0.5;
      } else if(fixed[0] > fixed[1]) {
        father = conn[0];
        mother = conn[1];
        move = 0;
      } else {
        father = conn[1];
        mother = conn[0];
        move = 0;
      }
      if(move == -1) {
        SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Don't know what to do");
      }


      Range out_tets;
      ierr = node_merger_ptr->mergeNodes(
        father,mother,out_tets,&in_tets,true,move
      ); CHKERRQ(ierr);
      if(!node_merger_ptr->getSucessMerge()) continue;

      Range adj_edge_faces;
      rval = moab.get_adjacencies(&mit->second,1,2,false,adj_edge_faces); CHKERRQ_MOAB(rval);
      Range adj_edge_faces_edges;
      rval = moab.get_adjacencies(
        adj_edge_faces,1,false,adj_edge_faces_edges,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      merged_edges.merge(adj_edge_faces_edges);

      // Update entities on surface
      Range adj_in_tets_faces;
      in_tets = subtract(out_tets,in_tets);
      rval = moab.get_adjacencies(
        in_tets,2,false,adj_in_tets_faces,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      adj_in_tets_faces = subtract(adj_in_tets_faces,cutNewSurfaces);
      for(
        Range::iterator fit = adj_in_tets_faces.begin();
        fit!=adj_in_tets_faces.end();fit++
      ) {
        EntityHandle parent_ent;
        rval = m_field.get_moab().tag_get_data(
          cOre.get_th_RefParentHandle(),&*fit,1,&parent_ent
        ); CHKERRQ_MOAB(rval);
        if(cutNewSurfaces.find(parent_ent)!=cutNewSurfaces.end()) {
          // Since parent of this face is part of the surface, add this child
          // entity.
          cutNewSurfaces.insert(*fit);
        }
      }
      cutNewSurfaces = subtract(cutNewSurfaces,adj_edge_faces);

      in_tets = out_tets;

    }

    ierr = m_field.seed_ref_level(in_tets,bit); CHKERRQ(ierr);
    Range bit_level_ents;
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,BitRefLevel().set(),MBTRI,bit_level_ents
    ); CHKERRQ(ierr);
    cutNewSurfaces = intersect(cutNewSurfaces,bit_level_ents);

    PetscFunctionReturn(0);
  }


  PetscErrorCode CutMeshInterface::getRayForEdge(
    const EntityHandle ent,
    VectorAdaptor ray_point,
    VectorAdaptor unit_ray_dir,
    double &ray_length
  ) const {
    MoABErrorCode rval;
    const MoFEM::Interface &m_field = cOre;
    const moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    int num_nodes;
    const EntityHandle *conn;
    rval = moab.get_connectivity(ent,conn,num_nodes,true); CHKERRQ_MOAB(rval);
    double coords[6];
    rval = moab.get_coords(conn,num_nodes,coords); CHKERR_MOAB(rval);
    VectorAdaptor s0(3,ublas::shallow_array_adaptor<double>(3,&coords[0]));
    VectorAdaptor s1(3,ublas::shallow_array_adaptor<double>(3,&coords[3]));
    noalias(ray_point) = s0;
    noalias(unit_ray_dir) = s1-s0;
    ray_length = norm_2(unit_ray_dir);
    unit_ray_dir /= ray_length;
    PetscFunctionReturn(0);
  }

  int CutMeshInterface::segmentPlane(
    VectorAdaptor s0,
    VectorAdaptor s1,
    VectorAdaptor x0,
    VectorAdaptor n,
    double &s
  ) const {
    VectorDouble3 u = s1 - s0;
    VectorDouble3 w = s0 - x0;
    double nu = inner_prod(n,u);
    double nw = -inner_prod(n,w);
    const double tol = 1e-4;
    if (fabs(nu) < tol) {           // segment is parallel to plane
        if (nw == 0)                      // segment lies in plane
            return 2;
        else
            return 0;                    // no intersection
    }
    // they are not parallel
    // compute intersect param
    s = nw / nu;
    if (s < 0 || s > 1)
        return 0;                        // no intersection
    return 1;
  }

}
