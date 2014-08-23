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

#ifdef WITH_TETGEM

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif

#include <MoFEM.hpp>
using namespace MoFEM;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>

#include <moab/Skinner.hpp>
//#include <moab/AdaptiveKDTree.hpp>

#include <moab/ParallelComm.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <FaceSplittingTool.hpp>
#include <TetGenInterface.hpp>

PetscErrorCode FaceSplittingTools::meshRefine(const int verb) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

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
  const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
  ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
  RefMoFEMEntity_multiIndex::index<Composite_EntType_and_ParentEntType_mi_tag>::type::iterator refit,hi_refit;
  refit = refinedEntitiesPtr_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().lower_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  hi_refit = refinedEntitiesPtr_ptr->get<Composite_EntType_and_ParentEntType_mi_tag>().upper_bound(boost::make_tuple(MBVERTEX,MBEDGE));
  Range already_refined_edges;
  for(;refit!=hi_refit;refit++) {
    EntityHandle parent_ent = refit->get_parent_ent(); 
    already_refined_edges.insert(parent_ent);
    RefMoFEMEntity_multiIndex::iterator parent_rit;
    parent_rit = refinedEntitiesPtr_ptr->find(parent_ent);
    if(parent_rit == refinedEntitiesPtr_ptr->end()) {
      SETERRQ1(PETSC_COMM_SELF,1,
	  "data inconsistency, entity in database not found %lu",refit->get_parent_ent());
    }
    bool success = const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)
	->modify(parent_rit,RefMoFEMEntity_change_add_bit(preserve_ref));
    if(!success) {
      SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessfull");
    }
  }

  int current_ref_bit = meshRefineBitLevels.first();
  for(int ll = 1;ll<nb_ref_levels+1;ll++) {

      BitRefLevel current_ref = BitRefLevel().set(current_ref_bit);

      Range level_tets;
      ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
      ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);

      Range level_edges;
      ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

      Range crack_edges;
      ierr = mField.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
      Range crack_edges_nodes;
      rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);

      //PetscAttachDebugger();
      Range crack_edge_nodes_parents;
      for(Range::iterator nit = crack_edges_nodes.begin();nit!=crack_edges_nodes.end();nit++) {
	EntityHandle ent = *nit;
	do {
	  RefMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator refit;
	  refit = refinedEntitiesPtr_ptr->get<MoABEnt_mi_tag>().find(ent);
	  if(refit == refinedEntitiesPtr_ptr->get<MoABEnt_mi_tag>().end()) {
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
      ierr = mField.query_interface(rEfiner); CHKERRQ(ierr);
      ierr = rEfiner->add_verices_in_the_middel_of_edges(edges_to_refine,last_ref,2); CHKERRQ(ierr);
      ierr = rEfiner->refine_TET(level_tets,last_ref,false); CHKERRQ(ierr);
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

  PetscErrorCode ierr;
  ErrorCode rval;

  const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
  ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
  BitRefLevel back_up_level;


  EntityHandle bit_meshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,bit_meshset); CHKERR_PETSC(rval);
  BitRefLevel current_ref = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,bit_meshset); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(bit_meshset); CHKERRQ(ierr);

  {

    EntityHandle meshset_interface;
    ierr = mField.get_Cubit_msId_meshset(200,SIDESET,meshset_interface); CHKERRQ(ierr);
    ierr = mField.query_interface(prismInterface); CHKERRQ(ierr);
    ierr = prismInterface->get_msId_3dENTS_sides(meshset_interface,current_ref,true,verb); CHKERRQ(ierr);

    int last_ref_bit = meshRefineBitLevels.back();
    if(!meshIntefaceBitLevels.empty()) {
      if(last_ref_bit<meshIntefaceBitLevels.back()) {
	last_ref_bit = meshIntefaceBitLevels.back();
      }
    }
    last_ref_bit++;
    meshIntefaceBitLevels.push_back(last_ref_bit);
    BitRefLevel last_ref = BitRefLevel().set(last_ref_bit);
  
    ierr = prismInterface->get_msId_3dENTS_split_sides(
      bit_meshset,last_ref,meshset_interface,false,true); CHKERRQ(ierr);

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

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = mField.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = mField.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = mField.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;
  bit_level0 = BitRefLevel().set(meshIntefaceBitLevels.back());

  if(verb>0) {    
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = mField.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      BitRefLevel last_ref = BitRefLevel().set(meshIntefaceBitLevels.back());
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,out_meshset); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("debug_split_faces.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
    }
  }

  PetscFunctionReturn(0);
}

#ifdef WITH_TETGEM

PetscErrorCode FaceSplittingTools::rebuildMeshWithTetGen(char switches[],const int verb) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;
  Skinner skin(&mField.get_moab());

  //get last refined bit level tets
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());//meshIntefaceBitLevels.back();
  Range mesh_level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);

  //crack front edges
  Range crack_edges;
  ierr = mField.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  //and nodes
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);

  //tets adj. to nodes
  Range crack_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  //get only tets on givem last bit ref level
  crack_edges_nodes_tets = intersect(crack_edges_nodes_tets,mesh_level_tets);

  //get second levelel tets
  Range crack_edges_nodes_tets_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges_nodes_tets,crack_edges_nodes_tets_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_tets_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes_tets_nodes,3,false,crack_edges_nodes_tets_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_nodes_tets_nodes_tets = intersect(crack_edges_nodes_tets_nodes_tets,mesh_level_tets);   
  Range crack_edges_nodes_tets_nodes_tets_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges_nodes_tets_nodes_tets,crack_edges_nodes_tets_nodes_tets_nodes,true); CHKERR_PETSC(rval);

  //get skins
  Range& tets1 = crack_edges_nodes_tets_nodes_tets;
  Range& nodes1 = crack_edges_nodes_tets_nodes_tets_nodes;
  Range tets1_faces;
  rval = mField.get_moab().get_adjacencies(
    tets1,2,true,tets1_faces,Interface::UNION); CHKERR_PETSC(rval);

  TetGenInterface *tetgen_iface;
  ierr = mField.query_interface(tetgen_iface); CHKERRQ(ierr);
  if(tetGenData.size()<1) {
    tetGenData.push_back(new tetgenio);
  }
  tetgenio &in = tetGenData.back();

  if(tetGenData.size()==1) {
    ierr = tetgen_iface->inData(nodes1,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
    //faces markers
    vector<pair<Range,int> > markers;
    Range merged_faces;
    //side set faces
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,sit)) {
      int id = sit->get_msId();
      Range faces;
      rval = mField.get_moab().get_entities_by_type(sit->meshset,MBTRI,faces,true); CHKERR_PETSC(rval);
      faces = intersect(faces,tets1_faces);
      Range::iterator it = faces.begin();
      for(;it!=faces.end();it++) {
	Range ent;
	ent.insert(*it);
	markers.push_back(pair<Range,int>(ent,id));
      }
      merged_faces.merge(faces);
    }
    //block skin faces
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,bit)) {
      Range block_tets;
      rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTET,block_tets,true); CHKERR_PETSC(rval);
      block_tets = intersect(block_tets,tets1);
      if(!block_tets.empty()) {
	Range block_tets_skin;
	rval = skin.find_skin(0,block_tets,false,block_tets_skin); CHKERR(rval);
	block_tets_skin = subtract(block_tets_skin,merged_faces);
	Range::iterator it = block_tets_skin.begin();
	for(;it!=block_tets_skin.end();it++) {
	  Range ent;
	  ent.insert(*it);
	  markers.push_back(pair<Range,int>(ent,-2));
	}
	merged_faces.merge(block_tets_skin);
      }
    }
    ierr = tetgen_iface->setFaceData(markers,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
    //make regions
    vector<pair<EntityHandle,int> > regions;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,bit)) {
      int id = bit->get_msId();
      Range block_tets;
      rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTET,block_tets,true); CHKERR_PETSC(rval);
      block_tets = intersect(block_tets,tets1);
      if(!block_tets.empty()) {
	Range::iterator it = block_tets.begin();
	for(;it!=block_tets.end();it++) {
	  regions.push_back(pair<EntityHandle,int>(*it,id));
	}
      }
    }
    ierr = tetgen_iface->setReginData(regions,in);  CHKERRQ(ierr);
  }

  //generate new mesh
  tetGenData.push_back(new tetgenio);
  tetgenio &out = tetGenData.back();
  char switches2[] = "pYAz";
  ierr = tetgen_iface->tetRahedralize(switches2,in,out); CHKERRQ(ierr);

  //get data from mesh
  int last_ref_bit = (meshRefineBitLevels.back()<meshIntefaceBitLevels.back()) ? meshIntefaceBitLevels.back() : meshRefineBitLevels.back();
  last_ref_bit++;
  meshRefineBitLevels.push_back(last_ref_bit);
  BitRefLevel last_ref = BitRefLevel().set(last_ref_bit);

  ierr = tetgen_iface->outData(in,out,moabTetGenMap,tetGenMoabMap,last_ref); CHKERRQ(ierr);
  
  map<int,Range> face_attributes_map;
  ierr = tetgen_iface->getTiangleAttributes(tetGenMoabMap,out,NULL,&face_attributes_map); CHKERRQ(ierr);
  map<int,Range> tet_region_map;
  ierr = tetgen_iface->getReginData(tetGenMoabMap,out,NULL,&tet_region_map); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,sit)) {
    int id = sit->get_msId();
    map<int,Range>::iterator mit = face_attributes_map.find(id);
    if(mit!=face_attributes_map.end()) {
      ierr = mField.get_moab().add_entities(sit->get_meshset(),mit->second); CHKERRQ(ierr);
    }
  }

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,bit)) {
    int id = bit->get_msId();  
    map<int,Range>::iterator mit = tet_region_map.find(id);
    if(mit!=tet_region_map.end()) {
      ierr = mField.get_moab().add_entities(bit->get_meshset(),mit->second); CHKERRQ(ierr);
    }
  }

  //add rest of elements to last bit level
  ierr = mField.seed_ref_level_3D(subtract(mesh_level_tets,tets1),last_ref); CHKERRQ(ierr);

  //split faces
  ierr = splitFaces(); CHKERRQ(ierr);

  if(verb>0) {
    EntityHandle meshset_level;
    rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_level); CHKERR_PETSC(rval);
    BitRefLevel last_int_bit_ref = last_ref;//BitRefLevel().set(meshIntefaceBitLevels.back());
    ierr = mField.get_entities_by_type_and_ref_level(last_int_bit_ref,BitRefLevel().set(),MBTET,meshset_level); CHKERRQ(ierr);
    //ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTRI,meshset_level); CHKERRQ(ierr);
    //rval = mField.get_moab().add_entities(meshset_level,tets1); CHKERR(rval);
    rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh.vtk","VTK","",&meshset_level,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&meshset_level,1); CHKERR_PETSC(rval);
    char tetgen_out_file_name[] = "out";
    out.save_elements(tetgen_out_file_name);
  }

  PetscFunctionReturn(0);
}

#endif


