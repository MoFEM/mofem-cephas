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
    obTreeSurfPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab,"ROOTSETSURF",true)
    );
    rval = obTreeSurfPtr->build(sUrface,rootSetSurf); CHKERRQ_MOAB(rval);
    Range faces;
    rval = moab.get_adjacencies(vOlume,2,false,faces,moab::Interface::UNION); CHKERRQ_MOAB(rval);
    obTreeVolPtr = boost::shared_ptr<OrientedBoxTreeTool>(
      new OrientedBoxTreeTool(&moab,"ROOTSETVOl",true)
    );
    rval = obTreeVolPtr->build(faces,rootSetVol); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::getRayForEdge(
    const EntityHandle ent,
    double *ray_point,
    double *unit_ray_dir,
    double *ray_length
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
    for(int d = 0;d!=3;d++) {
      ray_point[d] = coords[d];
      unit_ray_dir[d] = coords[3+d]-coords[d];
    }
    *ray_length = sqrt(
      unit_ray_dir[0]*unit_ray_dir[0]+
      unit_ray_dir[1]*unit_ray_dir[1]+
      unit_ray_dir[2]*unit_ray_dir[2]
    );
    for(int d = 0;d!=3;d++) {
      unit_ray_dir[d] /= *ray_length;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::findToCut(int verb) {
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
    cutVolumes.clear();
    for(Range::iterator eit = vol_edges.begin();eit!=vol_edges.end();eit++) {
      double ray_length;
      double ray_point[3],unit_ray_dir[3];
      ierr = getRayForEdge(*eit,ray_point,unit_ray_dir,&ray_length); CHKERRQ(ierr);
      std::vector< double > distances_out;
      std::vector< EntityHandle > facets_out;
      const double tol = ray_length*0.01;
      rval = obTreeSurfPtr->ray_intersect_triangles(
        distances_out,facets_out,rootSetSurf,tol,ray_point,unit_ray_dir,&ray_length
      ); CHKERR_MOAB(rval);
      if(!distances_out.empty()) {
        edgesToCut[*eit].dIst = distances_out[0];
        edgesToCut[*eit].lEngth = ray_length;
        VectorAdaptor adapt_vec_ray_dir(
          3,ublas::shallow_array_adaptor<double>(3,unit_ray_dir)
        );
        edgesToCut[*eit].unitRayDir = adapt_vec_ray_dir;
        VectorAdaptor adapt_vec_ray_point(
          3,ublas::shallow_array_adaptor<double>(3,ray_point)
        );
        edgesToCut[*eit].rayPoint = adapt_vec_ray_point;
        cutEdges.insert(*eit);
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
    rval = moab.get_adjacencies(
      cutEdges,3,false,cutVolumes,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::findTetOnTheFront(int verb) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    moab::Skinner skin(&moab);
    Range surf_skin;
    rval = skin.find_skin(0,sUrface,false,surf_skin); CHKERR_MOAB(rval);
    for(
      Range::iterator eit = surf_skin.begin();
      eit!=surf_skin.end();eit++
    ) {
      double ray_length;
      double ray_point[3],unit_ray_dir[3];
      ierr = getRayForEdge(*eit,ray_point,unit_ray_dir,&ray_length); CHKERRQ(ierr);
      std::vector< double > distances_out;
      std::vector< EntityHandle > facets_out;
      const double tol = ray_length*0.01;
      rval = obTreeVolPtr->ray_intersect_triangles(
        distances_out,facets_out,rootSetVol,tol,ray_point,unit_ray_dir,&ray_length
      ); CHKERR_MOAB(rval);
      if(!distances_out.empty()) {
        rval = moab.get_adjacencies(
          &facets_out[0],facets_out.size(),3,false,frontTets,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        VectorAdaptor adapt_vec_ray_dir(
          3,ublas::shallow_array_adaptor<double>(3,unit_ray_dir)
        );
        VectorAdaptor adapt_vec_ray_point(
          3,ublas::shallow_array_adaptor<double>(3,ray_point)
        );
        for(int i = 0;i!=distances_out.size();i++) {
          // Only get the cloaset
          facesFrontToRemesh[facets_out[i]].dIst = distances_out[i];
          facesFrontToRemesh[facets_out[i]].lEngth = ray_length;
          facesFrontToRemesh[facets_out[i]].unitRayDir = adapt_vec_ray_dir;
          facesFrontToRemesh[facets_out[i]].rayPoint = adapt_vec_ray_point;
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
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode CutMeshInterface::cutTets(
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
    for(
      map<EntityHandle,EdgeData>::iterator mit = edgesToCut.begin();
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
      }
    }
    cerr << cutNewVertices << endl;
    Range diff_verts;
    rval = moab.get_connectivity(cutNewSurfaces,diff_verts,true); CHKERRQ_MOAB(rval);
    diff_verts = subtract(diff_verts,cutNewVertices);
    Range subtract_faces;
    rval = moab.get_adjacencies(
      diff_verts,2,false,subtract_faces,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    cutNewSurfaces = subtract(cutNewSurfaces,subtract_faces);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CutMeshInterface::moveNodes() {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    for(
      map<EntityHandle,EdgeData>::iterator mit = verticecOnCutEdges.begin();
      mit!=verticecOnCutEdges.end(); mit++
    ) {
      double s = mit->second.dIst;
      // cout << s << " " << mit->second.dIst << " " << mit->second.lEngth << endl;
      VectorDouble3 new_coors = mit->second.rayPoint+s*mit->second.unitRayDir;
      rval = moab.set_coords(&mit->first,1,&new_coors[0]); CHKERRQ_MOAB(rval);
    }
    PetscFunctionReturn(0);
  }

  #ifdef WITH_TETGEN

  PetscErrorCode CutMeshInterface::imprintFront(
    const BitRefLevel bit,
    const double tetgen_face_angle,
    int verb
  ) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    const RefEntity_multiIndex *refined_ents_ptr;
    PetscFunctionBegin;
    ierr = m_field.get_ref_ents(&refined_ents_ptr); CHKERRQ(ierr);
    moab::Skinner skin(&moab);
    Range vol_skin;
    rval = skin.find_skin(0,vOlume,false,vol_skin); CHKERR_MOAB(rval);
    Range front_skin;
    rval = skin.find_skin(0,frontTets,false,front_skin); CHKERR_MOAB(rval);
    Range vertices;
    rval = moab.get_connectivity(frontTets,vertices,false); CHKERR_MOAB(rval);
    Range internal_front_skin = subtract(front_skin,vol_skin);
    for(
      Range::iterator fit = internal_front_skin.begin();
      fit!=internal_front_skin.end();
      fit++
    ) {
      map<EntityHandle,EdgeData>::iterator mit = facesFrontToRemesh.find(*fit);
      if(mit!=facesFrontToRemesh.end()) {
        double s = mit->second.dIst;
        VectorDouble3 new_coors = mit->second.rayPoint+s*mit->second.unitRayDir;
        EntityHandle vertex;
        rval = moab.create_vertex(&new_coors[0],vertex); CHKERRQ_MOAB(rval);
        vertices.insert(vertex);
        rval = moab.tag_set_data(
          cOre.get_th_RefParentHandle(),&vertex,1,&*fit
        ); CHKERRQ_MOAB(rval);
        rval = moab.tag_set_data(
          cOre.get_th_RefBitLevel(),&vertex,1,&bit
        ); CHKERRQ_MOAB(rval);
        std::pair<RefEntity_multiIndex::iterator,bool> p_ent =
        const_cast<RefEntity_multiIndex*>(refined_ents_ptr)->insert(
          boost::shared_ptr<RefEntity>(new RefEntity(m_field.get_basic_entity_data_ptr(),vertex))
        );
        if(!p_ent.second) {
          SETERRQ(
            m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"vertex not added"
          );
        }
      }

      map<int,Range> types_ents;
      //RIDGEVERTEX
      types_ents[TetGenInterface::RIDGEVERTEX].merge(vertices);
      //FREESEGVERTEX
      // types_ents[TetGenInterface::FREESEGVERTEX]
      //FREEFACETVERTEX
      // types_ents[TetGenInterface::FREEFACETVERTEX]
      //FREEVOLVERTEX
      // types_ents[TetGenInterface::FREEVOLVERTEX]

      Tag th_marker;
      int def_marker = 0;
      rval = moab.tag_get_handle(
        "TETGEN_MARKER",1,MB_TYPE_INTEGER,th_marker,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker
      ); CHKERRQ_MOAB(rval);
      vector<int> markers;
      fill(markers.begin(),markers.end(),0);
      rval = moab.tag_set_data(th_marker,vertices,&*markers.begin()); CHKERRQ_MOAB(rval);


      TetGenInterface *tetgen_iface;
      ierr = m_field.query_interface(tetgen_iface); CHKERRQ(ierr);

      tetGenData.clear();
      if(tetGenData.size()<1) {
        tetGenData.push_back(new tetgenio);
      }
      tetgenio &in = tetGenData.back();


      Range ents_to_tetgen;
      ents_to_tetgen.merge(vertices);
      ents_to_tetgen.merge(internal_front_skin);

      ierr = tetgen_iface->inData(ents_to_tetgen,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
      ierr = tetgen_iface->setGeomData(in,moabTetGenMap,tetGenMoabMap,types_ents); CHKERRQ(ierr);
      {
        vector<pair<Range,int> > markers;
        ierr = tetgen_iface->setFaceData(markers,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
      }

      if(verb>5) {
        if(m_field.getCommRank()==0) {
          char tetgen_in_file_name[] = "in";
          in.save_nodes(tetgen_in_file_name);
          in.save_elements(tetgen_in_file_name);
          in.save_faces(tetgen_in_file_name);
          in.save_edges(tetgen_in_file_name);
          in.save_poly(tetgen_in_file_name);
        }
      }

      std::string tetgen_switches;
      {
        ostringstream ss;
        ss << "rp" << tetgen_face_angle << "sqRS0JVV";
        tetgen_switches = ss.str();
      }

      tetGenData.push_back(new tetgenio);
      tetgenio &out = tetGenData.back();
      ierr = tetgen_iface->tetRahedralize(
        const_cast<char*>(tetgen_switches.c_str()),in,out
      ); CHKERRQ(ierr);

      //save elems
      if(verb>5) {
        char tetgen_out_file_name[] = "out";
        out.save_nodes(tetgen_out_file_name);
        out.save_elements(tetgen_out_file_name);
        out.save_faces(tetgen_out_file_name);
        out.save_edges(tetgen_out_file_name);
        out.save_poly(tetgen_out_file_name);
      }



    }

    PetscFunctionReturn(0);
  }

  #endif //WITH_TETGEN

}
