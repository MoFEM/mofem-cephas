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
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
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
    treeVolPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab,"ROOTSETVOl",true)
    );
    rval = treeVolPtr->build(faces,rootSetVol); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::findEdgesToCut(int verb) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    Range vol_edges;
    rval = moab.get_adjacencies(
      vOlume,1,true,vol_edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
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
    for(Range::iterator eit = vol_edges.begin();eit!=vol_edges.end();eit++) {
      ierr = getRayForEdge(
        *eit,vec_ray_point,vec_unit_ray_dir,ray_length
      ); CHKERRQ(ierr);
      std::vector< double > distances_out;
      std::vector< EntityHandle > facets_out;
      const double tol = ray_length*0.01;
      rval = treeSurfPtr->ray_intersect_triangles(
        distances_out,facets_out,rootSetSurf,tol,ray_point,unit_ray_dir,&ray_length
      ); CHKERR_MOAB(rval);
      if(!distances_out.empty()) {
        if(distances_out[0]/ray_length < 0.01) {
          int num_nodes;
          const EntityHandle *conn;
          rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
          verticesOnSurface.insert(conn[0]);
        } else if(distances_out[0]/ray_length > 0.99) {
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
      if(verb>0) {
        if(distances_out.size()>0) {
          std::cout << distances_out.size() << endl;
          for(int ii = 0;ii!=distances_out.size();ii++) {
            cout << "\t" << ii << " " << distances_out[ii] << endl;
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

    // Tag th;
    // double def_val[] = {0,0,0};
    // rval = moab.tag_get_handle("DIST",3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERRQ_MOAB(rval);

    cutEdgesOutside.clear();
    // Loop for all other edges, which are not cut by surface, but are in close
    // proximity to the surface
    for(Range::iterator eit = edges.begin();eit!=edges.end();eit++) {
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double coords[3*num_nodes];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERR_MOAB(rval);
      VectorAdaptor s0(3,ublas::shallow_array_adaptor<double>(3,&coords[0]));
      VectorAdaptor s1(3,ublas::shallow_array_adaptor<double>(3,&coords[3]));
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
      VectorAdaptor p0(3,ublas::shallow_array_adaptor<double>(3,point_out0));
      VectorAdaptor p1(3,ublas::shallow_array_adaptor<double>(3,point_out1));
      VectorDouble3 delta0,delta1;
      delta0 = p0-s0;
      delta1 = p1-s1;
      // moab.tag_set_data(th,&conn[0],1,&delta0[0]);
      // moab.tag_set_data(th,&conn[1],1,&delta1[0]);
      double dist0 = norm_2(delta0);
      double dist1 = norm_2(delta1);
      double normal[3];
      double *x0;
      if(dist0<dist1) {
        Util::normal(&moab,facets_out0,normal[0],normal[1],normal[2]);
        x0 = point_out0;
      } else {
        Util::normal(&moab,facets_out1,normal[0],normal[1],normal[2]);
        x0 = point_out1;
      }
      double s;
      int r = segmentPlane(
        s0,s1,
        VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,x0)),
        VectorAdaptor(3,ublas::shallow_array_adaptor<double>(3,&normal[0])),s
      );
      if(r == 1) {
        ierr = getRayForEdge(
          *eit,vec_ray_point,vec_unit_ray_dir,ray_length
        ); CHKERRQ(ierr);
        edgesToCut[*eit].dIst = s*ray_length;
        edgesToCut[*eit].lEngth = ray_length;
        edgesToCut[*eit].unitRayDir = vec_unit_ray_dir;
        edgesToCut[*eit].rayPoint = vec_ray_point;
        cutEdges.insert(*eit);
        cutEdgesOutside.insert(*eit);
      } else if(r == 2) {
        verticesOnSurface.insert(conn[0]);
        verticesOnSurface.insert(conn[1]);
      }
    }
    PetscFunctionReturn(0);
  }

  // PetscErrorCode CutMeshInterface::findTetOnTheFront(int verb) {
  //   MoABErrorCode rval;
  //   MoFEM::Interface &m_field = cOre;
  //   moab::Interface &moab = m_field.get_moab();
  //   PetscFunctionBegin;
  //   moab::Skinner skin(&moab);
  //   Range surf_skin;
  //   rval = skin.find_skin(0,sUrface,false,surf_skin); CHKERR_MOAB(rval);
  //   for(
  //     Range::iterator eit = surf_skin.begin();
  //     eit!=surf_skin.end();eit++
  //   ) {
  //     double ray_length;
  //     double ray_point[3],unit_ray_dir[3];
  //     VectorAdaptor vec_unit_ray_dir(
  //       3,ublas::shallow_array_adaptor<double>(3,unit_ray_dir)
  //     );
  //     VectorAdaptor vec_ray_point(
  //       3,ublas::shallow_array_adaptor<double>(3,ray_point)
  //     );
  //     ierr = getRayForEdge(*eit,vec_ray_point,vec_unit_ray_dir,ray_length); CHKERRQ(ierr);
  //     std::vector< double > distances_out;
  //     std::vector< EntityHandle > facets_out;
  //     const double tol = ray_length*0.01;
  //     rval = treeVolPtr->ray_intersect_triangles(
  //       distances_out,facets_out,rootSetVol,tol,ray_point,unit_ray_dir,&ray_length
  //     ); CHKERR_MOAB(rval);
  //     if(!distances_out.empty()) {
  //       rval = moab.get_adjacencies(
  //         &facets_out[0],facets_out.size(),3,false,frontTets,moab::Interface::UNION
  //       ); CHKERRQ_MOAB(rval);
  //       VectorAdaptor vec_ray_dir(
  //         3,ublas::shallow_array_adaptor<double>(3,unit_ray_dir)
  //       );
  //       VectorAdaptor vec_ray_point(
  //         3,ublas::shallow_array_adaptor<double>(3,ray_point)
  //       );
  //       for(int i = 0;i!=distances_out.size();i++) {
  //         // Only get the cloaset
  //         facesFrontToRemesh[facets_out[i]].dIst = distances_out[i];
  //         facesFrontToRemesh[facets_out[i]].lEngth = ray_length;
  //         facesFrontToRemesh[facets_out[i]].unitRayDir = vec_ray_dir;
  //         facesFrontToRemesh[facets_out[i]].rayPoint = vec_ray_point;
  //       }
  //       if(verb>0) {
  //         if(distances_out.size()>0) {
  //           std::cout << distances_out.size() << endl;
  //           for(int ii = 0;ii!=distances_out.size();ii++) {
  //             cout << "\t" << ii << " " << distances_out[ii] << endl;
  //           }
  //         }
  //       }
  //
  //     }
  //   }
  //   PetscFunctionReturn(0);
  // }

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
      bit,bit,MBTET,cutNewVolumes
    ); CHKERRQ(ierr);
    cutNewSurfaces.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,bit,MBTRI,cutNewSurfaces
    ); CHKERRQ(ierr);
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
    // cerr << cutNewVertices << endl;
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
    trimNewSurfaces.clear();
    ierr = m_field.get_entities_by_type_and_ref_level(
      bit,BitRefLevel().set(),MBTRI,trimNewSurfaces
    ); CHKERRQ(ierr);
    trimNewSurfaces = intersect(trimNewSurfaces,cutNewSurfaces);


    // cerr << trimNewVertices << endl;
    // Range diff_verts;
    // rval = moab.get_connectivity(cutNewSurfaces,diff_verts,true); CHKERRQ_MOAB(rval);
    // diff_verts = subtract(diff_verts,unite(cutNewVertices,verticesOnSurface));
    // Range subtract_faces;
    // rval = moab.get_adjacencies(
    //   diff_verts,2,false,subtract_faces,moab::Interface::UNION
    // ); CHKERRQ_MOAB(rval);
    // cutNewSurfaces = subtract(cutNewSurfaces,subtract_faces);
    PetscFunctionReturn(0);

  }

  PetscErrorCode CutMeshInterface::moveMidNodesOnCutEdges() {
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
      EntityHandle facets_out0;
      VectorDouble3 point_out0(3);
      rval = treeSurfPtr->closest_to_location(
        &new_coors[0],rootSetSurf,&point_out0[0],facets_out0
      ); CHKERR_MOAB(rval);
      dist = norm_2(point_out0-new_coors);
      const double tol = 1e-4;
      // if(dist>tol) {
      //   out_side_vertices.insert(&mit->first);
      // }
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
      int num_nodes;
      const EntityHandle *conn;
      rval = moab.get_connectivity(*eit,conn,num_nodes,true); CHKERRQ_MOAB(rval);
      double coords[3*num_nodes];
      rval = moab.get_coords(conn,num_nodes,coords); CHKERR_MOAB(rval);
      VectorAdaptor s0(3,ublas::shallow_array_adaptor<double>(3,&coords[0]));
      VectorAdaptor s1(3,ublas::shallow_array_adaptor<double>(3,&coords[3]));
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
      VectorAdaptor p0(3,ublas::shallow_array_adaptor<double>(3,point_out0));
      VectorAdaptor p1(3,ublas::shallow_array_adaptor<double>(3,point_out1));
      VectorDouble3 delta0,delta1;
      delta0 = p0-s0;
      delta1 = p1-s1;
      // moab.tag_set_data(th,&conn[0],1,&delta0[0]);
      // moab.tag_set_data(th,&conn[1],1,&delta1[0]);
      double dist0 = norm_2(delta0);
      double dist1 = norm_2(delta1);
      double min_dist = fmin(dist0,dist1);
      double max_dist = fmax(dist0,dist1);
      const double tol = 1e-4;
      if(min_dist < tol && max_dist > tol) {
        trimEdges.insert(*eit);
        edgesToTrim[*eit].dIst = dist0;
        edgesToTrim[*eit].rayPoint = s0;
        edgesToTrim[*eit].unitRayDir = (s1-s0);
        edgesToTrim[*eit].unitRayDir /= norm_2(edgesToTrim[*eit].unitRayDir);
      }
      if(min_dist > tol) {
        outsideEdges.insert(*eit);
      }
    }

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


  // #ifdef WITH_TETGEN
  //
  // PetscErrorCode CutMeshInterface::imprintFront(
  //   const BitRefLevel bit,
  //   const double tetgen_face_angle,
  //   int verb
  // ) {
  //   MoABErrorCode rval;
  //   MoFEM::Interface &m_field = cOre;
  //   moab::Interface &moab = m_field.get_moab();
  //   const RefEntity_multiIndex *refined_ents_ptr;
  //   PetscFunctionBegin;
  //   ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
  //   moab::Skinner skin(&moab);
  //   Range vol_skin;
  //   rval = skin.find_skin(0,vOlume,false,vol_skin); CHKERR_MOAB(rval);
  //   // {
  //   //   Range vertices;
  //   //   rval = moab.get_connectivity(frontTets,vertices,false); CHKERR_MOAB(rval);
  //   //   rval = moab.get_adjacencies(
  //   //     vertices,3,false,frontTets,moab::Interface::UNION
  //   //   ); CHKERR_MOAB(rval);
  //   //   frontTets = intersect(frontTets,vOlume);
  //   //   rval = moab.get_connectivity(frontTets,vertices,false); CHKERR_MOAB(rval);
  //   // }
  //   Range front_skin;
  //   rval = skin.find_skin(0,frontTets,false,front_skin); CHKERR_MOAB(rval);
  //   Range internal_front_skin = subtract(front_skin,vol_skin);
  //   Range front_tets_faces;
  //   rval = moab.get_adjacencies(
  //     frontTets,2,false,front_tets_faces,moab::Interface::UNION
  //   ); CHKERR_MOAB(rval);
  //   Range vertices;
  //   rval = moab.get_connectivity(frontTets,vertices,false); CHKERR_MOAB(rval);
  //   for(
  //     Range::iterator fit = front_tets_faces.begin();
  //     fit!=front_tets_faces.end();
  //     fit++
  //   ) {
  //     map<EntityHandle,TreeData>::iterator mit = facesFrontToRemesh.find(*fit);
  //     if(mit!=facesFrontToRemesh.end()) {
  //       double s = mit->second.dIst;
  //       VectorDouble3 new_coors = mit->second.rayPoint+s*mit->second.unitRayDir;
  //       EntityHandle vertex;
  //       rval = moab.create_vertex(&new_coors[0],vertex); CHKERRQ_MOAB(rval);
  //       vertices.insert(vertex);
  //       rval = moab.tag_set_data(cOre.get_th_RefParentHandle(),&vertex,1,&*fit); CHKERRQ_MOAB(rval);
  //       rval = moab.tag_set_data(cOre.get_th_RefBitLevel(),&vertex,1,&bit); CHKERRQ_MOAB(rval);
  //       std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
  //       const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
  //         boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),vertex))
  //       );
  //       if(!p_ent.second) {
  //         SETERRQ(
  //           m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"vertex not added"
  //         );
  //       }
  //     }
  //   }
  //
  //   map<int,Range> types_ents;
  //   //RIDGEVERTEX
  //   types_ents[TetGenInterface::RIDGEVERTEX].merge(vertices);
  //   //FREESEGVERTEX
  //   // types_ents[TetGenInterface::FREESEGVERTEX]
  //   //FREEFACETVERTEX
  //   // types_ents[TetGenInterface::FREEFACETVERTEX]
  //   //FREEVOLVERTEX
  //   // types_ents[TetGenInterface::FREEVOLVERTEX]
  //
  //   vector<pair<Range,int> > face_markers;
  //   for(
  //     Range::iterator fit = internal_front_skin.begin();
  //     fit!=internal_front_skin.end();fit++
  //   ) {
  //     Range facet;
  //     facet.insert(*fit);
  //     face_markers.push_back(pair<Range,int>(facet,1));
  //   }
  //
  //   Tag th_marker;
  //   int def_marker = 0;
  //   rval = moab.tag_get_handle(
  //     "TETGEN_MARKER",1,MB_TYPE_INTEGER,th_marker,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker
  //   ); CHKERRQ_MOAB(rval);
  //   vector<int> markers;
  //   markers.resize(vertices.size());
  //   fill(markers.begin(),markers.end(),0);
  //   rval = moab.tag_set_data(th_marker,vertices,&*markers.begin()); CHKERRQ_MOAB(rval);
  //   markers.resize(front_skin.size());
  //   fill(markers.begin(),markers.end(),0);
  //   rval = moab.tag_set_data(th_marker,front_skin,&*markers.begin()); CHKERRQ_MOAB(rval);
  //   Range front_skin_edges;
  //   rval = moab.get_adjacencies(
  //     front_skin,1,false,front_skin_edges,moab::Interface::UNION
  //   ); CHKERRQ_MOAB(rval);
  //   markers.resize(front_skin_edges.size());
  //   fill(markers.begin(),markers.end(),0);
  //   rval = moab.tag_set_data(th_marker,front_skin_edges,&*markers.begin()); CHKERRQ_MOAB(rval);
  //
  //   TetGenInterface *tetgen_iface;
  //   ierr = m_field.query_interface(tetgen_iface); CHKERRQ(ierr);
  //
  //   tetGenData.clear();
  //   tetGenData.push_back(new tetgenio);
  //   tetgenio &in = tetGenData.back();
  //
  //   Range ents_to_tetgen;
  //   ents_to_tetgen.merge(vertices);
  //   ents_to_tetgen.merge(internal_front_skin);
  //   rval = moab.get_adjacencies(
  //     internal_front_skin,1,false,ents_to_tetgen,moab::Interface::UNION
  //   ); CHKERRQ_MOAB(rval);
  //   ents_to_tetgen.merge(frontTets);
  //
  //   ierr = tetgen_iface->inData(ents_to_tetgen,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
  //   ierr = tetgen_iface->setGeomData(in,moabTetGenMap,tetGenMoabMap,types_ents); CHKERRQ(ierr);
  //   ierr = tetgen_iface->setFaceData(face_markers,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
  //
  //   if(verb>5) {
  //     if(m_field.getCommRank()==0) {
  //       char tetgen_in_file_name[] = "in";
  //       in.save_nodes(tetgen_in_file_name);
  //       in.save_elements(tetgen_in_file_name);
  //       in.save_faces(tetgen_in_file_name);
  //       in.save_edges(tetgen_in_file_name);
  //       in.save_poly(tetgen_in_file_name);
  //     }
  //   }
  //
  //   std::string tetgen_switches0;
  //   {
  //     ostringstream ss;
  //     ss << "piJYVV";
  //     // ss << "pYAVV" << endl;
  //     tetgen_switches0 = ss.str();
  //   }
  //   tetGenData.push_back(new tetgenio);
  //   {
  //     tetgenio &out = tetGenData.back();
  //     ierr = tetgen_iface->tetRahedralize(
  //       const_cast<char*>(tetgen_switches0.c_str()),in,out
  //     ); CHKERRQ(ierr);
  //     //save elems
  //     if(verb>5) {
  //       char tetgen_out_file_name[] = "out";
  //       out.save_nodes(tetgen_out_file_name);
  //       out.save_elements(tetgen_out_file_name);
  //       out.save_faces(tetgen_out_file_name);
  //       out.save_edges(tetgen_out_file_name);
  //       out.save_poly(tetgen_out_file_name);
  //     }
  //     //get mesh form TetGen and store it on second bit refined level
  //     ierr = tetgen_iface->outData(
  //       in,out,moabTetGenMap,tetGenMoabMap,bit,false,true
  //     ); CHKERRQ(ierr);
  //     ierr = m_field.seed_ref_level(subtract(vOlume,frontTets),bit); CHKERRQ(ierr);
  //   }
  //
  //
  //   PetscFunctionReturn(0);
  // }
  //
  // #endif //WITH_TETGEN

}
