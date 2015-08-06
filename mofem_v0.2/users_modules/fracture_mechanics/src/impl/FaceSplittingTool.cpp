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

#ifdef WITH_TETGEN

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif //WITH_TETGEN

#include <MoFEM.hpp>
using namespace MoFEM;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>

#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>

#include <moab/ParallelComm.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <FaceSplittingTool.hpp>
#include <TetGenInterface.hpp>
#include <NodeMerger.hpp>

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
      ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
      Range crack_edges_nodes;
      rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);

      Range crack_edges_nodes_tets;
      rval = mField.get_moab().get_adjacencies(
        crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION
      ); CHKERR_PETSC(rval);
      Range edges_to_refine;
      rval = mField.get_moab().get_adjacencies(crack_edges_nodes_tets,1,false,edges_to_refine,Interface::UNION); CHKERR_PETSC(rval);
      //edges_to_refine.merge(already_refined_edges);
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

      for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
        EntityHandle cubit_meshset = cubit_it->meshset;
        ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
      }

      current_ref_bit = meshRefineBitLevels.back();

  }

  {
    Range level_tets;
    BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
    ierr = mField.seed_ref_level(level_tets,BitRefLevel().set(BITREFLEVEL_SIZE-1)); CHKERRQ(ierr);
  }

  if(verb>5) {
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

  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);

  const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
  ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
  BitRefLevel back_up_level;

  Tag th_interface_side;
  rval = mField.get_moab().tag_get_handle("INTERFACE_SIDE",th_interface_side);
  if(rval == MB_SUCCESS) {
    rval = mField.get_moab().tag_delete(th_interface_side); CHKERR_PETSC(rval);
  }

  EntityHandle bit_meshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,bit_meshset); CHKERR_PETSC(rval);
  BitRefLevel current_ref = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.get_entities_by_type_and_ref_level(current_ref,BitRefLevel().set(),MBTET,bit_meshset); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(bit_meshset); CHKERRQ(ierr);

  {

    EntityHandle meshset_interface;
    ierr = mField.get_cubit_msId_meshset(200,SIDESET,meshset_interface); CHKERRQ(ierr);
    ierr = mField.query_interface(prismInterface); CHKERRQ(ierr);
    ierr = prismInterface->get_msId_3dENTS_sides(meshset_interface,current_ref,true,0); CHKERRQ(ierr);

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
      bit_meshset,last_ref,meshset_interface,true,true); CHKERRQ(ierr);

    //add refined ent to cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle cubit_meshset = cubit_it->meshset;
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,last_ref,cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }

    //remove tets which have 4 nodes on crack surface, those are source of problems. In case of planar crack, volume of such
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
    if(0) {//tets_on_surface.size()>0) { //swith off removing tets on surface
      Range tets_on_surface_faces;
      rval = mField.get_moab().get_adjacencies(
        tets_on_surface,2,false,tets_on_surface_faces,Interface::UNION
      ); CHKERR_PETSC(rval);
      tets_on_surface_faces = intersect(tets_on_surface_faces,interface_tris);

      if(pcomm->rank()>0) {
        if(verb>5) {
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

        RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator mit;
        mit = refinedEntitiesPtr_ptr->get<Ent_mi_tag>().find(*tit);
        if(mit == refinedEntitiesPtr_ptr->get<Ent_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,1,"no such tet in database");
        }
        bool success;
        success = const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)
        ->modify(mit,RefMoFEMEntity_change_set_nth_bit(last_ref_bit,false));
        if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsucceeded");
        success = const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)
        ->modify(mit,RefMoFEMEntity_change_set_nth_bit(BITREFLEVEL_SIZE-2,true));
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
        RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator mit;
        mit = refinedEntitiesPtr_ptr->get<Ent_mi_tag>().find(*eit);
        if(mit == refinedEntitiesPtr_ptr->get<Ent_mi_tag>().end()) {
          SETERRQ(PETSC_COMM_SELF,1,"no such tet in database");
        }
        bool success;
        //success = const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)
        //  ->modify(mit,RefMoFEMEntity_change_set_nth_bit(last_ref_bit,false));
        success = const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)
        ->modify(mit,RefMoFEMEntity_change_set_nth_bit(BITREFLEVEL_SIZE-2,true));
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

  if(verb>5) {
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

PetscErrorCode FaceSplittingTools::addCrackFront_to_Cubit201(int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Skinner skin(&mField.get_moab());
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());

  Range mesh_level_tets,mesh_level_tris;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);

  Range outer_surface_skin;
  rval = skin.find_skin(0,mesh_level_tets,false,outer_surface_skin); CHKERR_PETSC(rval);
  Range outer_surface_skin_edges;
  rval = mField.get_moab().get_adjacencies(outer_surface_skin,1,false,outer_surface_skin_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);

  Range crack_surface_tris_skin;
  rval = skin.find_skin(0,crack_surface_tris,false,crack_surface_tris_skin); CHKERR_PETSC(rval);

  crack_surface_tris_skin = subtract(crack_surface_tris_skin,outer_surface_skin_edges);
  EntityHandle meshset201;
  ierr = mField.get_cubit_msId_meshset(201,SIDESET,meshset201); CHKERRQ(ierr);
  rval = mField.get_moab().clear_meshset(&meshset201,1); CHKERR_PETSC(rval);
  ierr = mField.get_moab().add_entities(meshset201,crack_surface_tris_skin); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::roundCornersFillGaps_in_Cubit200(int nb,int verb) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  EntityHandle meshset200;
  ierr = mField.get_cubit_msId_meshset(200,SIDESET,meshset200); CHKERRQ(ierr);

  Skinner skin(&mField.get_moab());

  Range mesh_level_tris;
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);

  Range mesh_level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);

  Range crack_surface_tris_tets;
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);

  bool end;
  do {
    end = false;
    Range skin_edges;
    rval = skin.find_skin(0,crack_surface_tris,false,skin_edges); CHKERR(rval);
    Range skin_edges_tris;
    rval = mField.get_moab().get_adjacencies(
      skin_edges,2,false,skin_edges_tris,Interface::UNION); CHKERR_PETSC(rval);
    skin_edges_tris = intersect(skin_edges_tris,mesh_level_tris);
    skin_edges_tris = subtract(skin_edges_tris,crack_surface_tris);
    for(Range::iterator fit = skin_edges_tris.begin();fit!=skin_edges_tris.end();fit++) {
      Range adj_edges;
      rval = mField.get_moab().get_adjacencies(&*fit,1,1,false,adj_edges); CHKERR_PETSC(rval);
      adj_edges = intersect(adj_edges,skin_edges);
      if(adj_edges.size()>=(unsigned int)nb) {
        rval = mField.get_moab().get_adjacencies(
          crack_surface_tris,3,false,crack_surface_tris_tets,Interface::UNION
        ); CHKERR_PETSC(rval);
        crack_surface_tris_tets = intersect(crack_surface_tris_tets,mesh_level_tets);
        Range fit_tet;
        rval = mField.get_moab().get_adjacencies(&*fit,1,3,false,fit_tet); CHKERR_PETSC(rval);
        fit_tet = intersect(fit_tet,mesh_level_tets);
        //do not add faces which are adjacent to tets which have faces adjacent to existing crack surface
        if(intersect(fit_tet,crack_surface_tris_tets).size()==0) {
          if(verb>=0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"add face to meshset200\n");
          }
          ierr = mField.get_moab().add_entities(meshset200,&*fit,1); CHKERRQ(ierr);
          crack_surface_tris.insert(*fit);
          end = true;
        }
      }
    }
  } while(end);

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::crackFrontEdgeLengths(
  BitRefLevel bit_mesh,Range &to_split,Range &to_remove,int verb
) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Range mesh_level_tets,mesh_level_tris,mesh_level_edges,mesh_level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);

  //crack edges
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_edges;
  rval = mField.get_moab().get_adjacencies(crack_edges_nodes,1,false,crack_edges_nodes_edges,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_nodes_edges = intersect(crack_edges_nodes_edges,mesh_level_edges);
  Range crack_edges_nodes_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges_nodes_edges,crack_edges_nodes_edges_nodes,true); CHKERR_PETSC(rval);

  //crack surface
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris,crack_surface_tris_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_tris_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris,1,false,crack_surface_tris_edges,Interface::UNION
  ); CHKERR_PETSC(rval);
  Range crack_surface_tris_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris_nodes,3,false,crack_surface_tris_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_surface_tris_nodes_tets = intersect(crack_surface_tris_nodes_tets,mesh_level_tets);
  Range crack_surface_tris_nodes_tets_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris_nodes_tets,crack_surface_tris_nodes_tets_nodes,true); CHKERR_PETSC(rval);
  crack_surface_tris_nodes_tets_nodes = subtract(crack_surface_tris_nodes_tets_nodes,crack_surface_tris_nodes);

  //double def_VAL[1] = { 0 };
  //Tag th_length0;
  //rval = mField.get_moab().tag_get_handle(
    //"EDGE_LENGTH0",1,MB_TYPE_DOUBLE,th_length0,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

  map<EntityHandle,double> length_map1;
  Range::iterator eit;
  eit = crack_edges_nodes_edges.begin();
  for(;eit!=crack_edges_nodes_edges.end();eit++) {
    //double length0;
    //rval = mField.get_moab().tag_get_data(th_length0,&*eit,1,&length0); CHKERR(rval)
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[6];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,&coords[3],1,coords,1);
    double l = cblas_dnrm2(3,coords,1);
    length_map1[*eit] = l;
    //if(length0==0) {
      //rval = mField.get_moab().tag_set_data(th_length0,&*eit,1,&length0); CHKERR(rval)
    //}
  }

  EntityHandle kdTree_rootMeshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,kdTree_rootMeshset); CHKERR_PETSC(rval);
  AdaptiveKDTree kdTree(&mField.get_moab());
  rval = kdTree.build_tree(crack_surface_tris,&kdTree_rootMeshset); CHKERR_PETSC(rval);

  map<EntityHandle,double> length_map2;
  Range::iterator nit = crack_surface_tris_nodes_tets_nodes.begin();
  for(;nit!=crack_surface_tris_nodes_tets_nodes.end();nit++) {
    double coords[3];
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    EntityHandle triangle;
    double closest_point[3];
    rval = kdTree.closest_triangle(kdTree_rootMeshset,&coords[0],closest_point,triangle); CHKERR_PETSC(rval);
    cblas_daxpy(3,-1,closest_point,1,&coords[0],1);
    double l = cblas_dnrm2(3,&coords[0],1);
    length_map2[*nit] = l;
  }

  rval = mField.get_moab().delete_entities(&kdTree_rootMeshset,0); CHKERR_PETSC(rval);

  double ave_l1 = 0,max_l1 = 0,min_l1 = length_map1.begin()->second;
  for(
    map<EntityHandle,double>::iterator mit = length_map1.begin();
    mit!=length_map1.end();mit++
  ) {
    ave_l1 += mit->second;
    max_l1 = fmax(max_l1,mit->second);
    min_l1 = fmin(min_l1,mit->second);
  }
  ave_l1 /= length_map1.size();

  double sdev_l1 = 0;
  for(map<EntityHandle,double>::iterator mit = length_map1.begin();
    mit!=length_map1.end();mit++) {
    sdev_l1 += pow(mit->second-ave_l1,2);
  }
  sdev_l1 /= length_map1.size();
  sdev_l1 = sqrt(sdev_l1);

  double ave_l2 = 0,max_l2 = 0,min_l2 = length_map2.begin()->second;
  for(
    map<EntityHandle,double>::iterator mit = length_map2.begin();
    mit!=length_map2.end();mit++
  ) {
    ave_l2 += mit->second;
    max_l2 = fmax(max_l2,mit->second);
    min_l2 = fmin(min_l2,mit->second);
  }
  ave_l2 /= length_map2.size();

  double sdev_l2 = 0;
  for(map<EntityHandle,double>::iterator mit = length_map2.begin();
    mit!=length_map2.end();mit++) {
    sdev_l2 += pow(mit->second-ave_l2,2);
  }
  sdev_l2 /= length_map2.size();
  sdev_l2 = sqrt(sdev_l2);

  to_remove.clear();
  double diffNTET[12];
  ierr = ShapeDiffMBTET(diffNTET); CHKERRQ(ierr);
  eit = crack_edges_nodes.begin();
  for(;eit!=crack_edges_nodes.end();eit++) {
    Range::iterator eiit;
    Range adj_edges;
    rval = mField.get_moab().get_adjacencies(&*eit,1,1,false,adj_edges); CHKERR_PETSC(rval);
    adj_edges = intersect(adj_edges,crack_edges_nodes_edges);
    Range adj_tets;
    rval = mField.get_moab().get_adjacencies(&*eit,1,3,false,adj_tets); CHKERR_PETSC(rval);
    adj_tets = intersect(adj_tets,mesh_level_tets);
    Range bad_quality_edges_nodes;
    for(Range::iterator tit = adj_tets.begin();tit!=adj_tets.end();tit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = mField.get_moab().get_connectivity(*tit,conn,num_nodes,true); CHKERR_PETSC(rval);
      double coords[3*num_nodes];
      rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
      double coords_edges[2*3*6];
      ierr = get_edges_from_elem_coords(coords,coords_edges); CHKERRQ(ierr);
      double V = ShapeVolumeMBTET(diffNTET,&*coords);
      double alpha[4] = {1,1,1,1};
      double quality0,quality,b;
      ierr = quality_volume_length_F(
        V,alpha,0,
        diffNTET,coords_edges,coords,
        NULL,NULL,NULL,
        &quality0,&quality,&b,
        NULL,NULL
      );
      if(quality<0.1) {
        /*rval = mField.get_moab().get_adjacencies(
          &*tit,1,1,false,bad_quality_edges_nodes,Interface::UNION
        ); CHKERR_PETSC(rval);*/
        Range bad_quality_nodes;
        rval = mField.get_moab().get_adjacencies(
          &*tit,1,0,false,bad_quality_nodes,Interface::UNION
        ); CHKERR_PETSC(rval);
        bad_quality_edges_nodes.merge(bad_quality_nodes);
      }
    }

    Range adj_edges_on_crack_surface = intersect(adj_edges,crack_surface_tris_edges);
    adj_edges = subtract(adj_edges,crack_surface_tris_edges);
    Range adj_edges_nodes;
    rval = mField.get_moab().get_connectivity(adj_edges,adj_edges_nodes,true); CHKERR_PETSC(rval);
    adj_edges_nodes = subtract(adj_edges_nodes,crack_edges_nodes);
    if(intersect(crack_edges_nodes,adj_edges_nodes).size()>0) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    Range eit_to_remove;
    eiit = adj_edges_nodes.begin();
    for(;eiit != adj_edges_nodes.end();eiit++) {
      if(length_map2.find(*eiit)==length_map2.end()) {
        continue;
      }
      double l = length_map2[*eiit];
      if(l<fmax(ave_l2-sdev_l2,min_l2)) {
        eit_to_remove.insert(*eiit);
      }
    }

    to_remove.merge(intersect(eit_to_remove,bad_quality_edges_nodes));

    eiit = adj_edges_on_crack_surface.begin();
    for(;eiit != adj_edges_on_crack_surface.end();eiit++) {
      if(length_map1.find(*eiit)==length_map1.end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      double l = length_map1[*eiit];
      if(l>ave_l1+sdev_l1) {
        to_split.insert(*eiit);
      }
    }

  }

  if(intersect(crack_edges_nodes,to_remove).size()>0) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
  }

  /*Range quasi_flat_test;
  ierr = mField.get_entities_by_type_and_ref_level(BitRefLevel().set(BITREFLEVEL_SIZE-2),BitRefLevel().set(),MBTET,quasi_flat_test); CHKERRQ(ierr);
  Range quasi_flat_test_edges;
  rval = mField.get_moab().get_adjacencies(quasi_flat_test,1,false,quasi_flat_test_edges,Interface::UNION); CHKERR_PETSC(rval);
  to_split = subtract(to_split,quasi_flat_test_edges);
  to_split.clear();
  to_remove.clear();*/

  if(verb>=0) {
    ostringstream ss1;
    ss1 << "to split:\n " << to_split << endl << "to remove:\n " << to_remove << endl;
    ierr = PetscPrintf(PETSC_COMM_WORLD,ss1.str().c_str());
  }

  PetscFunctionReturn(0);
}
PetscErrorCode FaceSplittingTools::moveFrontNodesByVec(double v[]) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);


  Range mesh_level_edges;
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  if(crack_edges.empty()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"no crack front edges");
  }


  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);

  double coords[3];
  Range::iterator it = crack_edges_nodes.begin();
  for(;it!=crack_edges_nodes.end();it++) {
    rval = mField.get_moab().get_coords(&*it,1,coords); CHKERR_PETSC(rval);
    for(int dd = 0;dd<3;dd++) coords[dd] += v[dd];
    rval = mField.get_moab().set_coords(&*it,1,coords); CHKERR_PETSC(rval);
  }

  PetscFunctionReturn(0);
}

#ifdef WITH_TETGEN

PetscErrorCode FaceSplittingTools::rebuildMeshWithTetGen(vector<string> &switches,const int verb) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  {

    int back = meshRefineBitLevels.back();

    BitRefLevel maskPreserv;
    maskPreserv.set();
    maskPreserv[back] = false;

    //delete elements
    ierr = mField.remove_ents_from_finite_element("MATERIAL",0,MBTET); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("MATERIAL_COUPLED",0,MBTET); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("C_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("CTC_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("dCT_CRACKFRONT_AREA_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("GRIFFITH_FORCE_ELEMENT",0,MBTRI); CHKERRQ(ierr);
    //ierr = mField.remove_ents_from_finite_element("C_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    //ierr = mField.remove_ents_from_finite_element("CTC_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("CandCT_CRACK_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("BOTH_SIDE_OF_CRACK",0,MBPRISM); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("C_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("CTC_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_finite_element("CandCT_SURFACE_ELEM",0,MBTRI); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
      int msId = it->get_msId();
      if((msId < 10200)||(msId >= 10300)) continue;
      ostringstream ss0;
      ss0 << "C_SURFACE_ELEM_msId_" << msId;
      ostringstream ss1;
      ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
      ostringstream ss2;
      ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
      ierr = mField.remove_ents_from_finite_element(ss0.str(),0,MBTRI); CHKERRQ(ierr);
      ierr = mField.remove_ents_from_finite_element(ss1.str(),0,MBTRI); CHKERRQ(ierr);
      ierr = mField.remove_ents_from_finite_element(ss2.str(),0,MBTRI); CHKERRQ(ierr);
    }

    //delete field
    ierr = mField.remove_ents_from_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_CRACKFRONT_AREA",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_CRACKFRONT_AREA",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_CRACK_SURFACE",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_BOTH_SIDES",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("GRIFFITH_FORCE",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("GRIFFITH_FORCE_TANGENT",0,MBVERTEX); CHKERRQ(ierr);
    ierr = mField.remove_ents_from_field("LAMBDA_SURFACE",0,MBVERTEX); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
      int msId = it->get_msId();
      if((msId < 10200)||(msId >= 10300)) continue;
      ostringstream ss;
      ss << "LAMBDA_SURFACE_msId_" << msId;
      ierr = mField.remove_ents_from_field(ss.str(),0,MBVERTEX); CHKERRQ(ierr);
    }

    ierr = mField.delete_ents_by_bit_ref(maskPreserv,maskPreserv,true,0); CHKERRQ(ierr);

    //squash bits
    const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
    ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
    RefMoFEMEntity_multiIndex::iterator mit = refinedEntitiesPtr_ptr->begin();
    for(;mit!=refinedEntitiesPtr_ptr->end();mit++) {
      if(mit->get_ent_type() == MBENTITYSET) continue;
      if(mit->get_BitRefLevel().none()) continue;

      BitRefLevel new_bit;
      if((mit->get_BitRefLevel()&BitRefLevel().set(back)).any()) {
        new_bit.set(0);
        new_bit.set(BITREFLEVEL_SIZE-1);
      }

      const_cast<RefMoFEMEntity_multiIndex*>(refinedEntitiesPtr_ptr)->modify(mit,RefMoFEMEntity_change_set_bit(new_bit));
    }

    //set new memory
    bzero(meshRefineBitLevels.ptr,(BITREFLEVEL_SIZE-1)*sizeof(int));
    bzero(meshIntefaceBitLevels.ptr,(BITREFLEVEL_SIZE-1)*sizeof(int));
    meshRefineBitLevels.push_back(0);

    if(verb>5) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {
      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());

      EntityHandle meshset_out;
      //5
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTRI,meshset_out); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_5.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
      }
    }

  }

  Skinner skin(&mField.get_moab());

  Range to_split,to_remove;

  {
    BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
    ierr = crackFrontEdgeLengths(bit_mesh,to_split,to_remove,verb); CHKERRQ(ierr);

    // remove nodes on skin
    {

      Range mesh_level_edges;
      ierr = mField.get_entities_by_type_and_ref_level(
        bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges
      ); CHKERRQ(ierr);

      Range sideset_corners_edges;
      ierr = mField.get_cubit_msId_entities_by_dimension(
        100,SIDESET,1,sideset_corners_edges,true
      ); CHKERRQ(ierr);
      Range sideset_corners_nodes;
      rval = mField.get_moab().get_connectivity(
        sideset_corners_edges,sideset_corners_nodes,true
      ); CHKERR_PETSC(rval);
      sideset_corners_nodes = intersect(to_remove,sideset_corners_nodes);

      // get corners edges, between surface side-sets
      map<int,Range> surfaces_edges_map;
      Range surfaces;
      ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces,true); CHKERRQ(ierr);
      rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[102]); CHKERR_PETSC(rval);
      surfaces_edges_map[102] = intersect(surfaces_edges_map[102],mesh_level_edges);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
        int msId = it->get_msId();
        if((msId < 10200)||(msId >= 10300)) continue;
        EntityHandle meshset = it->get_meshset();
        surfaces.clear();
        rval = mField.get_moab().get_entities_by_type(meshset,MBTRI,surfaces,true); CHKERR_PETSC(rval);
        rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[msId]); CHKERR_PETSC(rval);
        surfaces_edges_map[msId] = intersect(surfaces_edges_map[msId],mesh_level_edges);
      }
      Range corners_edges;
      for(
        map<int,Range>::iterator mit = surfaces_edges_map.begin();
        mit!=surfaces_edges_map.end();mit++
      ) {
        corners_edges.merge(mit->second);
      }
      Range corners_nodes;
      rval = mField.get_moab().get_connectivity(
        corners_edges,corners_nodes,true
      ); CHKERR_PETSC(rval);
      corners_nodes = intersect(to_remove,corners_nodes);

      if(!corners_nodes.empty()) {

        // get nodes on crack front
        Range crack_edges;
        ierr = mField.get_cubit_msId_entities_by_dimension(
          201,SIDESET,1,crack_edges,true
        ); CHKERRQ(ierr);

        Range crack_edges_nodes;
        rval = skin.find_skin(0,crack_edges,false,crack_edges_nodes); CHKERR(rval);

        // get edges from side of skin of crack front
        Range crack_edges_nodes_edges;
        rval = mField.get_moab().get_adjacencies(
          crack_edges_nodes,1,true,crack_edges_nodes_edges,Interface::UNION
        ); CHKERR_PETSC(rval);
        crack_edges_nodes_edges = intersect(crack_edges_nodes_edges,mesh_level_edges);

        map<EntityHandle,EntityHandle> map_nodes_to_merge;
        for(
          Range::iterator eit = crack_edges_nodes_edges.begin();
          eit != crack_edges_nodes_edges.end();
          eit++
        ) {
          Range nodes;
          rval = mField.get_moab().get_connectivity(&*eit,1,nodes,true); CHKERR_PETSC(rval);
          EntityHandle front_node = intersect(crack_edges_nodes,nodes)[0];
          Range node_to_merge = intersect(corners_nodes,nodes);
          if(!node_to_merge.empty()) {
            map_nodes_to_merge[node_to_merge[0]] = front_node;
          }
        }
        NodeMergerInterface *node_merger_iface;
        ierr = mField.query_interface(node_merger_iface); CHKERRQ(ierr);
        for(
          map<EntityHandle,EntityHandle>::iterator mit = map_nodes_to_merge.begin();
          mit!=map_nodes_to_merge.end();mit++
        ) {
          int last_ref_bit = meshRefineBitLevels.back();
          last_ref_bit++;
          meshRefineBitLevels.push_back(last_ref_bit);
          BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
          to_remove.erase(mit->first);
          if(sideset_corners_nodes.find(mit->first)!=sideset_corners_nodes.end()) {
            ierr = node_merger_iface->mergeNodes(
              mit->first,mit->second,last_ref,bit_mesh
            ); CHKERRQ(ierr);
          } else {
            ierr = node_merger_iface->mergeNodes(
              mit->second,mit->first,last_ref,bit_mesh
            ); CHKERRQ(ierr);
          }
          for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
            EntityHandle meshset = cubit_it->meshset;
            ierr = mField.update_meshset_by_entities_children(meshset,last_ref,meshset,MBVERTEX,true); CHKERRQ(ierr);
            ierr = mField.update_meshset_by_entities_children(meshset,last_ref,meshset,MBEDGE,true); CHKERRQ(ierr);
            ierr = mField.update_meshset_by_entities_children(meshset,last_ref,meshset,MBTRI,true); CHKERRQ(ierr);
            ierr = mField.update_meshset_by_entities_children(meshset,last_ref,meshset,MBTET,true); CHKERRQ(ierr);
          }
          ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);
          if(verb>=0) {
            PetscPrintf(
              PETSC_COMM_WORLD,"Nodes merged %lu <- %lu\n",mit->second,mit->first
            );
          }
          bit_mesh = last_ref;
        }
        /*{
          EntityHandle meshset_out;
          rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
          BitRefLevel mask;
          mask.set();
          ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,mask,MBTET,meshset_out); CHKERRQ(ierr);
          rval = mField.get_moab().write_file("mesh_after_merge_node.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
          rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
        }*/
      }

    }

    Range level_tets,level_edges;
    ierr = mField.get_entities_by_type_and_ref_level(
      bit_mesh,BitRefLevel().set(),MBTET,level_tets
    ); CHKERRQ(ierr);
    rval = mField.get_moab().get_adjacencies(
      level_tets,1,true,level_edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    to_split = intersect(to_split,level_edges);
    if(to_split.size()>0) {
      int last_ref_bit = meshRefineBitLevels.back();
      if(!meshIntefaceBitLevels.empty()) {
        if(last_ref_bit<meshIntefaceBitLevels.back()) {
          last_ref_bit = meshIntefaceBitLevels.back();
        }
      }
      last_ref_bit++;
      meshRefineBitLevels.push_back(last_ref_bit);
      BitRefLevel new_ref = BitRefLevel().set(meshRefineBitLevels.back());
      ierr = mField.query_interface(rEfiner); CHKERRQ(ierr);
      ierr = rEfiner->add_verices_in_the_middel_of_edges(to_split,new_ref,2); CHKERRQ(ierr);
      ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);
      ierr = rEfiner->refine_TET(level_tets,new_ref,false); CHKERRQ(ierr);
      for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
        EntityHandle meshset = cubit_it->meshset;
        ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBVERTEX,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBEDGE,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTRI,true); CHKERRQ(ierr);
        ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTET,true); CHKERRQ(ierr);
      }
      ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);
    }
  }

  TetGenInterface *tetgen_iface;
  ierr = mField.query_interface(tetgen_iface); CHKERRQ(ierr);
  if(tetGenData.size()<1) {
    tetGenData.push_back(new tetgenio);
  }
  tetgenio &in = tetGenData.back();

  //get last refined bit level tets
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  Range mesh_level_tets,mesh_level_tris,mesh_level_edges,mesh_level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);
  Range mesh_level_tets_skin;
  rval = skin.find_skin(0,mesh_level_tets,false,mesh_level_tets_skin); CHKERR(rval);
  Range mesh_level_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(
    mesh_level_tets_skin,mesh_level_tets_skin_nodes,true
  ); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_edges;
  rval = mField.get_moab().get_adjacencies(
    mesh_level_tets_skin,1,false,mesh_level_tets_skin_edges,Interface::UNION
  ); CHKERR_PETSC(rval);

  //crack front edges
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  //tets edges
  Range crack_edges_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges,3,false,crack_edges_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_tets = intersect(crack_edges_tets,mesh_level_tets);
  //and nodes
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  crack_edges_nodes = intersect(crack_edges_nodes,mesh_level_nodes);
  Range crack_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_nodes_tets = intersect(crack_edges_nodes_tets,mesh_level_tets);
  Range crack_edges_nodes_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,1,false,crack_edges_nodes_edges,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_nodes_edges = intersect(crack_edges_nodes_edges,mesh_level_edges);

  //crack surfaces
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_tris_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris,1,false,crack_surface_tris_edges,Interface::UNION); CHKERR_PETSC(rval);
  //crack surface nodes, without crack front
  Range crack_surface_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris,crack_surface_tris_nodes,true); CHKERR_PETSC(rval);
  crack_surface_tris_nodes = subtract(crack_surface_tris_nodes,crack_edges_nodes);
  Range crack_surface_tris_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris_nodes,3,false,crack_surface_tris_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_surface_tris_nodes_tets = intersect(crack_surface_tris_nodes_tets,mesh_level_tets);
  Range crack_surface_tris_nodes_tets_nodes;
  rval = mField.get_moab().get_connectivity(
    crack_surface_tris_nodes_tets,crack_surface_tris_nodes_tets_nodes,true); CHKERR_PETSC(rval);
  //crack surface skin
  Range crack_surface_tris_skin;
  rval = skin.find_skin(0,crack_surface_tris,false,crack_surface_tris_skin); CHKERR(rval);
  Range crack_surface_tris_skin_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris_skin,crack_surface_tris_skin_nodes,true); CHKERR_PETSC(rval);

  Range corners_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  Range surfaces;
  ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces,true); CHKERRQ(ierr);
  map<int,Range> surfaces_edges_map;
  rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[102]); CHKERR_PETSC(rval);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    EntityHandle meshset = it->get_meshset();
    surfaces.clear();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTRI,surfaces,true); CHKERR_PETSC(rval);
    rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[msId]); CHKERR_PETSC(rval);
  }
  for(map<int,Range>::iterator mit = surfaces_edges_map.begin();mit!=surfaces_edges_map.end();mit++) {
    corners_edges.merge(mit->second);
  }
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);

  int crack_msId = -1;
  Range block_tets;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_ELASTICSET,it)) {
    block_tets.clear();
    EntityHandle meshset = it->get_meshset();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTET,block_tets,true); CHKERR_PETSC(rval);
    Range block_tets_edges;
    rval = mField.get_moab().get_adjacencies(
      block_tets,1,false,block_tets_edges,Interface::UNION); CHKERR_PETSC(rval);
    block_tets_edges = intersect(block_tets_edges,crack_edges);
    if(block_tets_edges.size()==crack_edges.size()) {
      crack_msId = it->get_msId();
      break;
    }
  }
  if(crack_msId == -1) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"can not find elastic block with crack front");
  }

  //get tets to remesh
  const int nb_tets_levels = 3;
  Range region_nodes = crack_edges_nodes;
  Range region_tets;
  for(int ii = 0;ii<nb_tets_levels;ii++) {
    rval = mField.get_moab().get_adjacencies(
      region_nodes,3,false,region_tets,Interface::UNION
    ); CHKERR_PETSC(rval);
    region_tets = intersect(region_tets,mesh_level_tets);
    region_tets = intersect(region_tets,block_tets);
    rval = mField.get_moab().get_connectivity(region_tets,region_nodes,true); CHKERR_PETSC(rval);
  }

  Range region_tets_faces;
  rval = mField.get_moab().get_adjacencies(
    region_tets,2,true,region_tets_faces,Interface::UNION
  ); CHKERR_PETSC(rval);
  Range region_tets_edges;
  rval = mField.get_moab().get_adjacencies(
    region_tets,1,true,region_tets_edges,Interface::UNION
  ); CHKERR_PETSC(rval);

  Range region_tets_skin;
  rval = skin.find_skin(0,region_tets,false,region_tets_skin); CHKERR(rval);
  region_tets_skin = intersect(region_tets_skin,region_tets_faces);
  Range region_tets_skin_without_boundary;
  region_tets_skin_without_boundary = subtract(region_tets_skin,mesh_level_tets_skin);
  region_tets_skin_without_boundary = intersect(region_tets_skin_without_boundary,region_tets_faces);
  Range region_tets_skin_boundary = subtract(region_tets_skin,region_tets_skin_without_boundary);

  Range region_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(
    region_tets_skin,region_tets_skin_nodes,true); CHKERR_PETSC(rval);

  Range region_tets_skin_without_boundary_nodes;
  rval = mField.get_moab().get_connectivity(
    region_tets_skin_without_boundary,region_tets_skin_without_boundary_nodes,true
  ); CHKERR_PETSC(rval);
  Range region_tets_skin_without_boundary_edges;
  rval = mField.get_moab().get_adjacencies(
    region_tets_skin_without_boundary,1,false,region_tets_skin_without_boundary_edges,Interface::UNION); CHKERR_PETSC(rval);
  region_tets_skin_without_boundary_edges = intersect(region_tets_skin_without_boundary_edges,region_tets_edges);

  Tag th_marker;
  int def_marker = 0;
  rval = mField.get_moab().tag_get_handle(
    "TETGEN_MARKER",1,MB_TYPE_INTEGER,th_marker,MB_TAG_CREAT|MB_TAG_SPARSE,&def_marker); CHKERR_PETSC(rval);
  vector<int> markers;

  crack_surface_tris = intersect(crack_surface_tris,region_tets_faces);
  crack_surface_tris_edges = intersect(crack_surface_tris_edges,region_tets_edges);

  map<int,Range> types_ents;
  //RIDGEVERTEX
  types_ents[TetGenInterface::RIDGEVERTEX].merge(region_tets_skin_without_boundary_nodes);
  //FREESEGVERTEX
  types_ents[TetGenInterface::FREESEGVERTEX].merge(crack_surface_tris_skin_nodes);
  //FREEFACETVERTEX
  types_ents[TetGenInterface::FREEFACETVERTEX].merge(region_tets_skin_nodes);
  types_ents[TetGenInterface::FREEFACETVERTEX] =  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEFACETVERTEX] =  subtract(types_ents[TetGenInterface::FREEFACETVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  //FREEVOLVERTEX
  types_ents[TetGenInterface::FREEVOLVERTEX].merge(region_nodes);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::RIDGEVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREESEGVERTEX]);
  types_ents[TetGenInterface::FREEVOLVERTEX] = subtract(types_ents[TetGenInterface::FREEVOLVERTEX],types_ents[TetGenInterface::FREEFACETVERTEX]);

  //boundary and nodes to remove

  markers.resize(region_tets_faces.size());
  fill(markers.begin(),markers.end(),0);
  rval = mField.get_moab().tag_set_data(th_marker,region_tets_faces,&*markers.begin()); CHKERR_PETSC(rval);
  markers.resize(region_tets_edges.size());
  fill(markers.begin(),markers.end(),0);
  rval = mField.get_moab().tag_set_data(th_marker,region_tets_edges,&*markers.begin()); CHKERR_PETSC(rval);
  markers.resize(region_nodes.size());
  fill(markers.begin(),markers.end(),0);
  rval = mField.get_moab().tag_set_data(th_marker,region_nodes,&*markers.begin()); CHKERR_PETSC(rval);

  markers.resize(region_tets_skin_without_boundary.size());
  fill(markers.begin(),markers.end(),-1);
  rval = mField.get_moab().tag_set_data(th_marker,region_tets_skin_without_boundary,&*markers.begin()); CHKERR_PETSC(rval);
  markers.resize(region_tets_skin_without_boundary_edges.size());
  fill(markers.begin(),markers.end(),-1);
  rval = mField.get_moab().tag_set_data(th_marker,region_tets_skin_without_boundary_edges,&*markers.begin()); CHKERR_PETSC(rval);

  int shift = 0;
  map<int,int> id_shift_map;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
    int msId = it->get_msId();
    id_shift_map[msId] = 1<<shift;
    shift++;
    Range surfaces_faces_msId;
    ierr = mField.get_cubit_msId_entities_by_dimension(msId,SIDESET,2,surfaces_faces_msId,true); CHKERRQ(ierr);
    surfaces_faces_msId = intersect(surfaces_faces_msId,region_tets_faces);
    markers.resize(surfaces_faces_msId.size());
    rval = mField.get_moab().tag_get_data(th_marker,surfaces_faces_msId,&*markers.begin()); CHKERR_PETSC(rval);
    for(unsigned int ii = 0;ii<markers.size();ii++) {
      markers[ii] |= id_shift_map[msId];
    }
    rval = mField.get_moab().tag_set_data(th_marker,surfaces_faces_msId,&*markers.begin()); CHKERR_PETSC(rval);
  }



  Range nodes_to_remove = subtract(to_remove.subset_by_type(MBVERTEX),

    unite(region_tets_skin_without_boundary_nodes,crack_surface_tris_nodes)
  );
  nodes_to_remove = subtract(nodes_to_remove,corners_edges_nodes);

  markers.resize(nodes_to_remove.size());
  fill(markers.begin(),markers.end(),-1);
  rval = mField.get_moab().tag_set_data(th_marker,nodes_to_remove,&*markers.begin()); CHKERR_PETSC(rval);

  //nodes
  if(tetGenData.size()==1) {

    Range ents_to_tetgen;
    ents_to_tetgen.merge(region_nodes);
    ents_to_tetgen.merge(region_tets);

    ents_to_tetgen.merge(crack_surface_tris);
    ents_to_tetgen.merge(crack_edges);
    ents_to_tetgen.merge(region_tets_skin_without_boundary);
    ents_to_tetgen.merge(region_tets_skin_without_boundary_edges);
    ents_to_tetgen.merge(region_tets_skin_boundary);

    if(verb > 10) {
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
      if(pcomm->rank()==0) {

      //BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());
      EntityHandle meshset_out;
      //0
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(meshset_out,ents_to_tetgen); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_0.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
      }
    }

    ierr = tetgen_iface->inData(ents_to_tetgen,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);
    ierr = tetgen_iface->setGeomData(in,moabTetGenMap,tetGenMoabMap,types_ents); CHKERRQ(ierr);

    vector<pair<Range,int> > markers;
    for(Range::iterator tit = crack_surface_tris.begin();tit!=crack_surface_tris.end();tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range,int>(facet,200));
    }
    for(Range::iterator tit = region_tets_skin_without_boundary.begin();tit!=region_tets_skin_without_boundary.end();tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range,int>(facet,1));
    }
    Range other_facets = subtract(region_tets_skin,region_tets_skin_without_boundary);
    for(Range::iterator tit = other_facets.begin();tit!=other_facets.end();tit++) {
      Range facet;
      facet.insert(*tit);
      markers.push_back(pair<Range,int>(facet,0));
    }
    ierr = tetgen_iface->setFaceData(markers,in,moabTetGenMap,tetGenMoabMap); CHKERRQ(ierr);

  }

  if(verb>5) {
    char tetgen_in_file_name[] = "in";
    in.save_nodes(tetgen_in_file_name);
    in.save_elements(tetgen_in_file_name);
    in.save_faces(tetgen_in_file_name);
    in.save_edges(tetgen_in_file_name);
    in.save_poly(tetgen_in_file_name);
  }

  //generate new mesh
  {
    vector<string>::iterator sw = switches.begin();
    for(int ii = 0;sw!=switches.end();sw++,ii++) {
      tetgenio &_in_ = tetGenData.back();
      tetGenData.push_back(new tetgenio);
      tetgenio &_out_ = tetGenData.back();
      char *s = const_cast<char*>(sw->c_str());
      ierr = tetgen_iface->tetRahedralize(s,_in_,_out_); CHKERRQ(ierr);
    }
  }
  tetgenio &out = tetGenData.back();
  //save elems
  if(verb>5) {
    char tetgen_out_file_name[] = "out";
    out.save_nodes(tetgen_out_file_name);
    out.save_elements(tetgen_out_file_name);
    out.save_faces(tetgen_out_file_name);
    out.save_edges(tetgen_out_file_name);
    out.save_poly(tetgen_out_file_name);
  }

  //get data from mesh
  int last_ref_bit = (meshRefineBitLevels.back()<meshIntefaceBitLevels.back()) ? meshIntefaceBitLevels.back() : meshRefineBitLevels.back();
  last_ref_bit++;
  meshRefineBitLevels.push_back(last_ref_bit);
  BitRefLevel last_ref = BitRefLevel().set(last_ref_bit);

  ierr = tetgen_iface->outData(in,out,moabTetGenMap,tetGenMoabMap,last_ref,false,false); CHKERRQ(ierr);

  Range last_ref_tets;
  ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,last_ref_tets); CHKERRQ(ierr);
  EntityHandle meshset_block;
  ierr = mField.get_cubit_msId_meshset(crack_msId,BLOCKSET,meshset_block); CHKERRQ(ierr);
  ierr = mField.get_moab().add_entities(meshset_block,last_ref_tets); CHKERRQ(ierr);

  Range tetgen_faces;
  map<int,Range> face_markers_map;
  ierr = tetgen_iface->getTriangleMarkers(tetGenMoabMap,out,&tetgen_faces,&face_markers_map); CHKERRQ(ierr);

  for(
    map<int,Range>::iterator mit = face_markers_map.begin();
    mit!=face_markers_map.end();mit++
  ) {
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
      int msId = it->get_msId();
      if(id_shift_map[msId]&mit->first) {
        EntityHandle meshset = it->get_meshset();
        ierr = mField.get_moab().add_entities(meshset,mit->second.subset_by_type(MBTRI)); CHKERRQ(ierr);
      }
    }
  }

  if(verb>5) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {

    BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());

    EntityHandle meshset_out;
    //1
    rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTRI,meshset_out); CHKERRQ(ierr);
    rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_1.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
    //2
    rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
    ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,meshset_out); CHKERRQ(ierr);
    rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_2.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
    rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);

    }
  }

  ierr = mField.seed_ref_level(
    subtract(mesh_level_tets,region_tets),last_ref
  ); CHKERRQ(ierr);
  //add rest of elements to last bit level
  ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);
  ierr = roundCornersFillGaps_in_Cubit200(2,2); CHKERRQ(ierr);
  ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);

  //split faces
  ierr = splitFaces(0); CHKERRQ(ierr);

  if(verb>5) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {

      BitRefLevel last_int_ref = BitRefLevel().set(meshIntefaceBitLevels.back());

      EntityHandle meshset_out;
      //3
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(last_int_ref,BitRefLevel().set(),MBTRI,meshset_out); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_3.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
      //4
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(last_int_ref,BitRefLevel().set(),MBTET,meshset_out); CHKERRQ(ierr);
      rval = mField.get_moab().write_file("rebuild_with_split_faces_tet_gen_mesh_4.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);

    }
  }

  PetscFunctionReturn(0);
}

#endif //WITH_TETGEN

PetscErrorCode FaceSplittingTools::cornerProblem(Range &new_nodes,int verb) {
  PetscFunctionBegin;

  if(new_nodes.empty()) PetscFunctionReturn(0);

  PetscErrorCode ierr;
  ErrorCode rval;

  // Get skin of level test, nodes and edges on the skin
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  Range mesh_level_tets,mesh_level_tris,mesh_level_edges,mesh_level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);
  Range mesh_level_tets_skin;
  Skinner skin(&mField.get_moab());
  rval = skin.find_skin(0,mesh_level_tets,false,mesh_level_tets_skin); CHKERR(rval);
  Range mesh_level_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets_skin,mesh_level_tets_skin_nodes,true); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_edges;
  rval = mField.get_moab().get_adjacencies(
    mesh_level_tets_skin,1,false,mesh_level_tets_skin_edges,Interface::UNION
  ); CHKERR_PETSC(rval);

  // Crack front edges, and nodes
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_on_skin;
  crack_edges_nodes_on_skin = intersect(crack_edges_nodes,mesh_level_tets_skin_nodes);

  Range crack_edges_nodes_tris;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,2,false,crack_edges_nodes_tris,Interface::UNION
  ); CHKERR_PETSC(rval);
  Range new_nodes_tris;
  rval = mField.get_moab().get_adjacencies(
    new_nodes,2,false,new_nodes_tris,Interface::UNION
  ); CHKERR_PETSC(rval);
  new_nodes_tris.merge(crack_edges_nodes_tris);


  Range new_nodes_tris_nodes;
  rval = mField.get_moab().get_connectivity(new_nodes_tris,new_nodes_tris_nodes,true); CHKERR_PETSC(rval);
  new_nodes_tris_nodes = subtract(new_nodes_tris_nodes,crack_edges_nodes);
  new_nodes_tris_nodes = subtract(new_nodes_tris_nodes,new_nodes);
  Range new_nodes_tris_nodes_tris;
  rval = mField.get_moab().get_adjacencies(
    new_nodes_tris_nodes,2,false,new_nodes_tris_nodes_tris,Interface::UNION
  ); CHKERR_PETSC(rval);
  new_nodes_tris = intersect(new_nodes_tris,mesh_level_tris);
  new_nodes_tris = subtract(new_nodes_tris,new_nodes_tris_nodes_tris);
  new_nodes_tris = subtract(new_nodes_tris,mesh_level_tets_skin);

  //cerr << "new_nodes_tris.size() " << new_nodes_tris.size() << endl;

  EntityHandle meshset200;
  ierr = mField.get_cubit_msId_meshset(200,SIDESET,meshset200); CHKERRQ(ierr);
  rval = mField.get_moab().add_entities(meshset200,new_nodes_tris); CHKERR_PETSC(rval);

  //add rest of elements to last bit level
  ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);
  ierr = roundCornersFillGaps_in_Cubit200(2,2); CHKERRQ(ierr);
  ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);

  new_nodes_tris_nodes.clear();
  rval = mField.get_moab().get_connectivity(new_nodes_tris,new_nodes_tris_nodes,true); CHKERR_PETSC(rval);

  for(
    map<EntityHandle,vector<double> >::iterator mit = nodesToMoveMap.begin();
    mit!=nodesToMoveMap.end();mit++
  ) {
    if(new_nodes_tris_nodes.find(mit->first)!=new_nodes_tris_nodes.end()) {
      rval = mField.get_moab().set_coords(&mit->first,1,&*mit->second.begin()); CHKERR_PETSC(rval);
    }
  }

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle meshset_out;
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(meshset_out,new_nodes_tris); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("new_nodes_tris.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::getCornerEdges(Range &edges_to_cut,int verb) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  // get entities at bit mesh level
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  Range mesh_level_tets,mesh_level_tris,mesh_level_edges,mesh_level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);
  // get entities on skin
  Range mesh_level_tets_skin;
  Skinner skin(&mField.get_moab());
  rval = skin.find_skin(0,mesh_level_tets,false,mesh_level_tets_skin); CHKERR(rval);
  Range mesh_level_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets_skin,mesh_level_tets_skin_nodes,true); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_edges;
  rval = mField.get_moab().get_adjacencies(mesh_level_tets_skin,1,false,mesh_level_tets_skin_edges,Interface::UNION); CHKERR_PETSC(rval);

  //corner edges and nodes
  Range corners_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  corners_edges = intersect(corners_edges,mesh_level_edges);
  Range corners_edges_nodes;
  rval = mField.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);

  // crack front edges, and nodes
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,1,false,crack_edges_nodes_edges,Interface::UNION
  ); CHKERR_PETSC(rval);
  // get crack front nodes on body skin
  Range crack_edges_nodes_on_skin;
  crack_edges_nodes_on_skin = intersect(crack_edges_nodes,mesh_level_tets_skin_nodes);
  Range crack_edges_nodes_on_skin_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes_on_skin,1,false,crack_edges_nodes_on_skin_edges,Interface::UNION
  ); CHKERR_PETSC(rval);

  // remove entities on skin
  Tag th_freez;
  const int def_order = 0;
  rval = mField.get_moab().tag_get_handle(
    "FROZEN_NODE",1,MB_TYPE_INTEGER,th_freez,MB_TAG_CREAT|MB_TAG_SPARSE,&def_order
  ); CHKERR_PETSC(rval);
  for(
    Range::iterator nit = crack_edges_nodes_on_skin.begin();
    nit!=crack_edges_nodes_on_skin.end();nit++
  ) {
    if(corners_edges_nodes.find(*nit)!=crack_edges_nodes.end()) {
      continue;
    }
    int freez = 0;
    rval = mField.get_moab().tag_get_data(th_freez,&*nit,1,&freez); CHKERR_PETSC(rval);
    if(freez != 0) {
      nit = crack_edges_nodes_on_skin.erase(nit);
      //cerr << *nit << endl;
    }
  }

  //cerr << "crack_edges_nodes_on_skin " << crack_edges_nodes_on_skin.size() << endl;

  // if no nodes on skin, do nothing
  if(crack_edges_nodes_on_skin.empty()) {
    PetscFunctionReturn(0);
  }

  // crack surfaces
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  // get edges on skin of crack surface
  Range crack_surface_skin;
  rval = skin.find_skin(0,crack_surface_tris,false,crack_surface_skin); CHKERR(rval);
  // get crack surface edges
  Range crack_surface_tris_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris,1,false,crack_surface_tris_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_surface_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris,crack_surface_tris_nodes,true); CHKERR_PETSC(rval);
  crack_surface_tris_nodes = subtract(crack_surface_tris_nodes,crack_edges_nodes);
  Range crack_surface_tris_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris_nodes,3,false,crack_surface_tris_nodes_tets,Interface::UNION
  ); CHKERR_PETSC(rval);

  //map surfaces
  Range surfaces;
  ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces,true); CHKERRQ(ierr);
  surfaces = intersect(surfaces,mesh_level_tris);
  map<int,Range> surfaces_edges_map;
  rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[102]); CHKERR_PETSC(rval);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    EntityHandle meshset = it->get_meshset();
    surfaces.clear();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTRI,surfaces,true); CHKERR_PETSC(rval);
    surfaces = intersect(surfaces,mesh_level_tris);
    rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[msId]); CHKERR_PETSC(rval);
  }
  // free edges are edges between two surface side sets, nodes can move along such edge
  Range free_edges;
  for(map<int,Range>::iterator mit = surfaces_edges_map.begin();mit!=surfaces_edges_map.end();mit++) {
    free_edges.merge(mit->second);
    //for(map<int,Range>::iterator miit = surfaces_edges_map.begin();miit!=surfaces_edges_map.end();miit++) {
      //free_edges.merge(intersect(mit->second,miit->second));
    //}
  }
  Range free_edges_nodes;
  rval = mField.get_moab().get_connectivity(free_edges,free_edges_nodes,true); CHKERR_PETSC(rval);
  free_edges_nodes = subtract(free_edges_nodes,corners_edges_nodes);

  Range front_nodes_on_skin = subtract(crack_edges_nodes_on_skin,free_edges_nodes);

  //cerr << "front_nodes_on_skin " << front_nodes_on_skin.size() << endl;

  Range nodes_on_corners;
  for(
    Range::iterator nit = front_nodes_on_skin.begin();nit!=front_nodes_on_skin.end();nit++
  ) {
    // get test adjacent to crack front node
    Range nit_tets;
    rval = mField.get_moab().get_adjacencies(
      &*nit,1,3,false,nit_tets,Interface::UNION
    ); CHKERR_PETSC(rval);
    nit_tets = intersect(nit_tets,mesh_level_tets);
    // not adjacent to crack surface
    nit_tets = subtract(nit_tets,crack_surface_tris_nodes_tets);
    // get edges
    Range nit_tets_edges;
    rval = mField.get_moab().get_adjacencies(
      nit_tets,1,false,nit_tets_edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    /*if(corners_edges_nodes.find(*nit)==corners_edges_nodes.end()) {
      // check if any node is on corner
      Range nit_tets_nodes;
      rval = mField.get_moab().get_connectivity(
        nit_tets,nit_tets_nodes,true
      ); CHKERR_PETSC(rval);
      Range intersecting_corners;;
      intersecting_corners = intersect(nit_tets_nodes,unite(corners_edges_nodes,free_edges_nodes));
      if(intersecting_corners.empty()) continue;
      nodes_on_corners.merge(intersecting_corners);
    }*/
    if(intersect(nit_tets_edges,unite(corners_edges,free_edges)).empty()) continue;
    // check if any edge is on the corner
    if(intersect(nit_tets_edges,unite(corners_edges,free_edges)).empty()) continue;
    nit_tets_edges = subtract(nit_tets_edges,crack_edges);
    nit_tets_edges = subtract(nit_tets_edges,crack_surface_tris_edges);
    // select edges to cut
    edges_to_cut.merge(nit_tets_edges);
  }

  if(nodes_on_corners.size()>0) {
    Range edges;
    rval = mField.get_moab().get_adjacencies(
      nodes_on_corners,1,false,edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    edges = subtract(edges,crack_edges);
    edges = subtract(edges,crack_surface_tris_edges);
    // select edges to cut
    edges_to_cut.merge(edges);
  }

  //cerr << "edges_to_cut " << edges_to_cut.size() << endl;

  if(edges_to_cut.size()>0) {
    // get edges to cut which intersect with corner edges
    Range corners_edges2;
    corners_edges2 = intersect(edges_to_cut,unite(corners_edges,free_edges));
    // get nodes which intersect with corner edges
    Range corners_edges_nodes2;
    rval = mField.get_moab().get_connectivity(
      corners_edges2,corners_edges_nodes2,true
    ); CHKERR_PETSC(rval);
    // get adjacent tets
    Range tets;
    rval = mField.get_moab().get_adjacencies(
      corners_edges_nodes2,3,false,tets,Interface::UNION
    ); CHKERR_PETSC(rval);
    tets = intersect(tets,mesh_level_tets);
    // edges of those tets
    Range tets_edges;
    rval = mField.get_moab().get_adjacencies(
      tets,1,false,tets_edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    // remove all edges of crack front
    tets_edges = subtract(tets_edges,crack_edges);
    tets_edges = subtract(tets_edges,crack_surface_tris_edges);
    edges_to_cut.merge(tets_edges);
  }

  // do not cut edges which coming from crack front itself
  edges_to_cut = subtract(edges_to_cut,crack_edges_nodes_edges);

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      EntityHandle meshset_out;
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      rval = mField.get_moab().add_entities(meshset_out,edges_to_cut); CHKERR_PETSC(rval);
      rval = mField.get_moab().write_file("edges_to_cut.vtk","VTK","",&meshset_out,1); CHKERR_PETSC(rval);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode FaceSplittingTools::propagateBySplit(Range &new_nodes,Range &edges_to_cut,int verb) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  // get bit level entities
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  Range mesh_level_tets,mesh_level_tris,mesh_level_edges;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  // get entities on skin
  Skinner skin(&mField.get_moab());
  Range mesh_level_tets_skin;
  rval = skin.find_skin(0,mesh_level_tets,false,mesh_level_tets_skin); CHKERR(rval);
  Range mesh_level_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets_skin,mesh_level_tets_skin_nodes,true); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_edges;
  rval = mField.get_moab().get_adjacencies(mesh_level_tets_skin,1,false,mesh_level_tets_skin_edges,Interface::UNION); CHKERR_PETSC(rval);

  // Crack front edges
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION); CHKERR_PETSC(rval);
  crack_edges_nodes_tets = intersect(crack_edges_nodes_tets,mesh_level_tets);
  Range crack_edges_tets_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes_tets,1,false,crack_edges_tets_edges,Interface::UNION
  ); CHKERR_PETSC(rval);

  // Crack surfaces
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  mesh_level_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_skin;
  rval = skin.find_skin(0,crack_surface_tris,false,crack_surface_skin); CHKERR(rval);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris,crack_surface_tris_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_tris_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris,1,false,crack_surface_tris_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_surface_skin_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_skin,crack_surface_skin_nodes,true); CHKERR_PETSC(rval);

  Range t1_side; // tets only on one side of crack surface
  {
    Range n3; //crack surface nodes without crack front nodes
    n3 = subtract(crack_surface_tris_nodes,crack_edges_nodes);
    n3 = subtract(n3,crack_surface_skin_nodes);
    Range t1; //crack tets adjacent to crack surface but not to crack front
    rval = mField.get_moab().get_adjacencies(
      n3,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    t1 = intersect(t1,mesh_level_tets);
    Range e0; //carck surface edges without crack front edges
    e0 = subtract(crack_surface_tris_edges,crack_edges);
    e0 = subtract(e0,crack_surface_skin);
    rval = mField.get_moab().get_adjacencies(
      e0,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    rval = mField.get_moab().get_adjacencies(
      crack_surface_tris,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    t1 = intersect(t1,mesh_level_tets);

    unsigned int nb_t1_side;
    t1_side.insert(*t1.begin());
    do {
      nb_t1_side = t1_side.size();
      Range f2;
      rval = mField.get_moab().get_adjacencies(
        t1_side,2,false,f2,Interface::UNION
      ); CHKERR_PETSC(rval);
      f2 = subtract(f2,crack_surface_tris);
      rval = mField.get_moab().get_adjacencies(
        f2,3,false,t1_side,Interface::UNION
      ); CHKERR_PETSC(rval);
      t1_side = intersect(t1_side,t1);
    } while(nb_t1_side != t1_side.size());
    Range t1_side_tris;
    rval = mField.get_moab().get_adjacencies(t1_side,2,false,t1_side_tris,Interface::UNION); CHKERR_PETSC(rval);
    crack_surface_tris = intersect(crack_surface_tris,t1_side_tris);
  }

  Range nodes_to_check;
  rval = mField.get_moab().get_connectivity(edges_to_cut,nodes_to_check,true); CHKERR_PETSC(rval);
  //cerr << nodes_to_check << endl;

  // calculate normals
  double diffN[3*2];
  ierr = ShapeDiffMBTRI(diffN); CHKERRQ(ierr);
  map<EntityHandle,vector<double> > normal_map;
  for(
    Range::iterator fit = crack_surface_tris.begin();
    fit!=crack_surface_tris.end();fit++
  ) {
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    normal_map[*fit].resize(3);
    double *normal = &*normal_map[*fit].begin();
    ierr = ShapeFaceNormalMBTRI(diffN,coords,normal); CHKERRQ(ierr);
    double nrm2 = cblas_dnrm2(3,normal,1);
    cblas_dscal(3,1./nrm2,normal,1);
    //set direction
    Range t3;
    rval = mField.get_moab().get_adjacencies(&*fit,1,3,false,t3,Interface::UNION); CHKERR_PETSC(rval);
    t3 = intersect(t3,t1_side);
    if(t3.size()!=1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    int side_number;
    int sense;
    int offset;
    rval = mField.get_moab().side_number(t3[0],*fit,side_number,sense,offset); CHKERR_PETSC(rval);
    if(sense == -1) {
      cblas_dscal(3,-1,normal,1);
    }
  }

  //build tree
  EntityHandle kdTree_rootMeshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,kdTree_rootMeshset); CHKERR_PETSC(rval);
  AdaptiveKDTree kdTree(&mField.get_moab());
  rval = kdTree.build_tree(crack_surface_tris,&kdTree_rootMeshset); CHKERR_PETSC(rval);

  Tag th_distance;
  if(verb>0) {
    double def_VAL[1] = { 0 };
    rval = mField.get_moab().tag_get_handle(
      "DISTANCE_FROM_CRACK_SURFACE",1,MB_TYPE_DOUBLE,th_distance,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL
    ); CHKERR_PETSC(rval);
  }

  // calculate distance from crack surface
  signedDistanceMap.clear();
  normalNodesMap.clear();
  for(Range::iterator nit = nodes_to_check.begin();nit!=nodes_to_check.end();nit++) {
    double coords[3];
    double closest_point[3];
    EntityHandle triangle;
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    rval = kdTree.closest_triangle(kdTree_rootMeshset,coords,closest_point,triangle); CHKERR_PETSC(rval);
    double delta[3];
    cblas_dcopy(3,closest_point,1,delta,1);
    cblas_daxpy(3,-1,coords,1,delta,1);
    map<EntityHandle,vector<double> >::iterator mit = normal_map.find(triangle);
    if(mit == normal_map.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    double *normal = &*mit->second.begin();
    double dot = cblas_ddot(3,normal,1,delta,1);
    signedDistanceMap[*nit] = dot;
    normalNodesMap[*nit] = mit->second;
    if(verb>0) {
      rval = mField.get_moab().tag_set_data(th_distance,&*nit,1,&dot); CHKERR_PETSC(rval);
    }
  }

  rval = mField.get_moab().delete_entities(&kdTree_rootMeshset,1); CHKERR_PETSC(rval);
  //cerr << "normalNodesMap.size() " << normalNodesMap.size() << endl;

  // mode nodes near to crack surface
  nodesToMoveMap.clear();
  map<EntityHandle,vector<double> >::iterator mit = normalNodesMap.begin();
  for(;mit!=normalNodesMap.end();mit++) {
    EntityHandle node = mit->first;
    double *normal = &*mit->second.begin();
    Range node_edges;
    rval = mField.get_moab().get_adjacencies(&node,1,1,false,node_edges,Interface::UNION); CHKERR_PETSC(rval);
    double max_l = 0;
    for(Range::iterator eit = node_edges.begin();eit!=node_edges.end();eit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
      double coords[3*num_nodes];
      rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
      cblas_daxpy(3,-1,&coords[3],1,coords,1);
      double l = fabs(cblas_ddot(3,normal,1,coords,1));
      max_l = fmax(l,max_l);
    }
    double dot = signedDistanceMap[node];
    if(fabs(dot/max_l) < 0.125) {
      //double coords[3];
      nodesToMoveMap[node].resize(3);
      double *coords = &*nodesToMoveMap[node].begin();
      rval = mField.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
      cblas_daxpy(3,dot,normal,1,coords,1);
      //rval = mField.get_moab().set_coords(&node,1,coords); CHKERR_PETSC(rval);
      signedDistanceMap[node] = 0;
      if(verb>0) {
        double dot0 = 0;
        rval = mField.get_moab().tag_set_data(th_distance,&node,1,&dot0); CHKERR_PETSC(rval);
      }
    }
  }

  for(
    map<EntityHandle,vector<double> >::iterator mit = nodesToMoveMap.begin();
    mit!=nodesToMoveMap.end();mit++
  ) {
    if(mesh_level_tets_skin_nodes.find(mit->first)!=mesh_level_tets_skin_nodes.end()) {
      //rval = mField.get_moab().set_coords(&mit->first,1,&*mit->second.begin()); CHKERR_PETSC(rval);
      new_nodes.insert(mit->first);
    }
  }

  //cerr << "nodesToMoveMap.size() " << nodesToMoveMap.size() << endl;
  //cerr << "new_nodes.size() " << new_nodes.size() << endl;

  Range edges_to_split;
  edgesToSplitMap.clear();
  mit = normalNodesMap.begin();
  for(;mit!=normalNodesMap.end();mit++) {
    EntityHandle node = mit->first;
    double dot = signedDistanceMap[node];
    if(dot>=0) continue;
    Range node_edges;
    rval = mField.get_moab().get_adjacencies(
      &node,1,1,false,node_edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    node_edges = intersect(node_edges,edges_to_cut);
    //cerr << node_edges.size() << endl;
    for(Range::iterator eit = node_edges.begin();eit!=node_edges.end();eit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
      EntityHandle node2 = (conn[0] == node) ? conn[1] : conn[0];
      if(signedDistanceMap.find(node2) == signedDistanceMap.end()) {
        continue;
      }
      double dot2 = signedDistanceMap[node2];
      //cerr << "\t" << dot*dot2 << endl;
      if(dot*dot2 >= 0) {
        continue;
      }
      double coords[3];
      rval = mField.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
      double coords2[3];
      rval = mField.get_moab().get_coords(&node2,1,coords2); CHKERR_PETSC(rval);
      cblas_daxpy(3,-1,coords,1,coords2,1);
      //double l = cblas_dnrm2(3,coords2,1);
      //dot2 = dot + a*l
      //a = (dot2-dot)/l
      //x = -dot*l/(dot2-dot)
      double x = -dot/(dot2-dot);
      cblas_daxpy(3,x,coords2,1,coords,1);
      edgesToSplitMap[*eit].resize(3);
      cblas_dcopy(3,coords,1,&*edgesToSplitMap[*eit].begin(),1);
      edges_to_split.insert(*eit);
    }
  }

  //cerr << edges_to_split << endl;

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {

      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());

      EntityHandle meshset_out;
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,meshset_out); CHKERRQ(ierr);
      rval = mField.get_moab().add_entities(meshset_out,edges_to_split); CHKERR_PETSC(rval);
      mField.get_moab().write_file("edges_to_split.vtk","VTK","",&meshset_out,1);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);

    }

  }

  if(edges_to_split.empty()) {
    PetscFunctionReturn(0);
  }

  for(
    map<EntityHandle,vector<double> >::iterator mit = nodesToMoveMap.begin();
    mit!=nodesToMoveMap.end();mit++
  ) {
    //rval = mField.get_moab().set_coords(&mit->first,1,&*mit->second.begin()); CHKERR_PETSC(rval);
    new_nodes.insert(mit->first);
  }

  //split edges
  Range level_tets;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  int last_ref_bit = meshRefineBitLevels.back();
  if(!meshIntefaceBitLevels.empty()) {
    if(last_ref_bit<meshIntefaceBitLevels.back()) {
      last_ref_bit = meshIntefaceBitLevels.back();
    }
  }
  last_ref_bit++;
  meshRefineBitLevels.push_back(last_ref_bit);
  BitRefLevel new_ref = BitRefLevel().set(meshRefineBitLevels.back());
  ierr = mField.query_interface(rEfiner); CHKERRQ(ierr);
  ierr = rEfiner->add_verices_in_the_middel_of_edges(edges_to_split,new_ref,2); CHKERRQ(ierr);
  ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);
  ierr = rEfiner->refine_TET(level_tets,new_ref,false); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
    EntityHandle meshset = cubit_it->meshset;
    ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBVERTEX,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBEDGE,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTRI,true); CHKERRQ(ierr);
    ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTET,true); CHKERRQ(ierr);
  }
  ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);

  Range nodes_on_split_edges;
  for(Range::iterator eit = edges_to_split.begin();eit!=edges_to_split.end();eit++) {
    const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
    ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
    RefMoFEMEntity_multiIndex::index<Composite_Ent_And_ParentEntType_mi_tag>::type::iterator it;
    it = refinedEntitiesPtr_ptr->get<Composite_Ent_And_ParentEntType_mi_tag>().find(boost::make_tuple(*eit,MBVERTEX));
    if(it == refinedEntitiesPtr_ptr->get<Composite_Ent_And_ParentEntType_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    EntityHandle node = it->get_ref_ent();
    nodes_on_split_edges.insert(node);
    rval = mField.get_moab().set_coords(&node,1,&*edgesToSplitMap[*eit].begin()); CHKERR_PETSC(rval);
    new_nodes.insert(node);
  }

  if(verb>0) {
    ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {

      BitRefLevel last_ref = BitRefLevel().set(meshRefineBitLevels.back());

      EntityHandle meshset_out;
      rval = mField.get_moab().create_meshset(MESHSET_SET,meshset_out); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(last_ref,BitRefLevel().set(),MBTET,meshset_out); CHKERRQ(ierr);
      mField.get_moab().write_file("split_sides.vtk","VTK","",&meshset_out,1);
      rval = mField.get_moab().delete_entities(&meshset_out,1); CHKERR_PETSC(rval);

    }

  }

  PetscFunctionReturn(0);
}

/*PetscErrorCode FaceSplittingTools::getNodesAndEdgesOnBodySurface(Range &edges_and_nodes,int verb = 0) {
  PetscFunctionBegin;

  // get entities at bit mesh level
  BitRefLevel bit_mesh = BitRefLevel().set(meshRefineBitLevels.back());
  Range mesh_level_tets,mesh_level_tris,mesh_level_edges,mesh_level_nodes;
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,mesh_level_tets); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTRI,mesh_level_tris); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBEDGE,mesh_level_edges); CHKERRQ(ierr);
  ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBVERTEX,mesh_level_nodes); CHKERRQ(ierr);
  // get entities on skin
  Range mesh_level_tets_skin;
  Skinner skin(&mField.get_moab());
  rval = skin.find_skin(0,mesh_level_tets,false,mesh_level_tets_skin); CHKERR(rval);
  Range mesh_level_tets_skin_nodes;
  rval = mField.get_moab().get_connectivity(mesh_level_tets_skin,mesh_level_tets_skin_nodes,true); CHKERR_PETSC(rval);
  Range mesh_level_tets_skin_edges;
  rval = mField.get_moab().get_adjacencies(mesh_level_tets_skin,1,false,mesh_level_tets_skin_edges,Interface::UNION); CHKERR_PETSC(rval);

  //map surfaces
  Range surfaces;
  ierr = mField.get_cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces,true); CHKERRQ(ierr);
  surfaces = intersect(surfaces,mesh_level_tris);
  map<int,Range> surfaces_edges_map;
  rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[102]); CHKERR_PETSC(rval);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    EntityHandle meshset = it->get_meshset();
    surfaces.clear();
    rval = mField.get_moab().get_entities_by_type(meshset,MBTRI,surfaces,true); CHKERR_PETSC(rval);
    surfaces = intersect(surfaces,mesh_level_tris);
    rval = skin.find_skin(0,surfaces,false,surfaces_edges_map[msId]); CHKERR_PETSC(rval);
  }
  // free edges are edges between two surface side sets, nodes can move along such edge
  Range free_edges;
  for(map<int,Range>::iterator mit = surfaces_edges_map.begin();mit!=surfaces_edges_map.end();mit++) {
    free_edges.merge(mit->second);
  }


  // Crack surfaces
  Range crack_surface_tris;
  ierr = mField.get_cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surface_tris,true); CHKERRQ(ierr);
  mesh_level_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_skin;
  rval = skin.find_skin(0,crack_surface_tris,false,crack_surface_skin); CHKERR(rval);
  crack_surface_tris = intersect(crack_surface_tris,mesh_level_tris);
  Range crack_surface_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_tris,crack_surface_tris_nodes,true); CHKERR_PETSC(rval);
  Range crack_surface_tris_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_surface_tris,1,false,crack_surface_tris_edges,Interface::UNION); CHKERR_PETSC(rval);
  Range crack_surface_skin_nodes;
  rval = mField.get_moab().get_connectivity(crack_surface_skin,crack_surface_skin_nodes,true); CHKERR_PETSC(rval);

  // get edges on body surface bridged by tets to body surface
  Range edges_on_surface;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes_tets,1,false,edges_on_surface,Interface::UNION
  ); CHKERR_PETSC(rval);
  edges_on_surface = intersect(edges_on_surface,mesh_level_tets_skin_edges);
  Range edges_on_surface_nodes;
  rval = mField.get_moab().get_connectivity(
    edges_on_surface,edges_on_surface_nodes,true
  ); CHKERR_PETSC(rval);
  edges_on_surface_nodes = subtract(edges_on_surface_nodes,crack_surface_tris_nodes);
  edges_on_surface_nodes = subtract(edges_on_surface_nodes,free_edges_nodes);

  // crack front edges, and nodes
  Range crack_edges;
  ierr = mField.get_cubit_msId_entities_by_dimension(201,SIDESET,1,crack_edges,true); CHKERRQ(ierr);
  crack_edges = intersect(crack_edges,mesh_level_edges);
  Range crack_edges_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges,crack_edges_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_tris;
  rval = mField.get_moab().get_adjacencies(
    crack_edges,2,false,crack_edges_tris,Interface::UNION
  ); CHKERR_PETSC(rval);
  crack_edges_tris = subtract(crack_edges_tris,crack_surface_tris);
  crack_edges_tris = intersect(crack_edges_tris,mesh_level_tris);
  Range crack_edges_tris_nodes;
  rval = mField.get_moab().get_connectivity(crack_edges_tris,crack_edges_tris_nodes,true); CHKERR_PETSC(rval);
  Range crack_edges_nodes_tets;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,3,false,crack_edges_nodes_tets,Interface::UNION
  ); CHKERR_PETSC(rval);
  crack_edges_nodes_tets = intersect(crack_edges_nodes_tets,mesh_level_tets);
  Range crack_edges_nodes_edges;
  rval = mField.get_moab().get_adjacencies(
    crack_edges_nodes,1,false,crack_edges_nodes_edges,Interface::UNION
  ); CHKERR_PETSC(rval);
  // this get edges on free edges on which crack is propagating
  free_edges = intersect(free_edges,crack_edges_nodes_edges);
  Range free_edges_nodes;
  rval = mField.get_moab().get_connectivity(free_edges,free_edges_nodes,true); CHKERR_PETSC(rval);

  Range t1_side; // tets only on one side of crack surface
  {
    Range n3; //crack surface nodes without crack front nodes
    n3 = subtract(crack_surface_tris_nodes,crack_edges_nodes);
    n3 = subtract(n3,crack_surface_skin_nodes);
    Range t1; //crack tets adjacent to crack surface but not to crack front
    rval = mField.get_moab().get_adjacencies(
      n3,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    t1 = intersect(t1,mesh_level_tets);
    Range e0; //carck surface edges without crack front edges
    e0 = subtract(crack_surface_tris_edges,crack_edges);
    e0 = subtract(e0,crack_surface_skin);
    rval = mField.get_moab().get_adjacencies(
      e0,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    rval = mField.get_moab().get_adjacencies(
      crack_surface_tris,3,false,t1,Interface::UNION
    ); CHKERR_PETSC(rval);
    t1 = intersect(t1,mesh_level_tets);

    unsigned int nb_t1_side;
    t1_side.insert(*t1.begin());
    do {
      nb_t1_side = t1_side.size();
      Range f2;
      rval = mField.get_moab().get_adjacencies(
        t1_side,2,false,f2,Interface::UNION
      ); CHKERR_PETSC(rval);
      f2 = subtract(f2,crack_surface_tris);
      rval = mField.get_moab().get_adjacencies(
        f2,3,false,t1_side,Interface::UNION
      ); CHKERR_PETSC(rval);
      t1_side = intersect(t1_side,t1);
    } while(nb_t1_side != t1_side.size());
    Range t1_side_tris;
    rval = mField.get_moab().get_adjacencies(t1_side,2,false,t1_side_tris,Interface::UNION); CHKERR_PETSC(rval);
    crack_surface_tris = intersect(crack_surface_tris,t1_side_tris);
  }

  // calculate normals
  double diffN[3*2];
  ierr = ShapeDiffMBTRI(diffN); CHKERRQ(ierr);
  map<EntityHandle,vector<double> > normal_map;
  for(
    Range::iterator fit = crack_surface_tris.begin();
    fit!=crack_surface_tris.end();fit++
  ) {
    int num_nodes;
    const EntityHandle* conn;
    rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    double coords[3*num_nodes];
    rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
    normal_map[*fit].resize(3);
    double *normal = &*normal_map[*fit].begin();
    ierr = ShapeFaceNormalMBTRI(diffN,coords,normal); CHKERRQ(ierr);
    double nrm2 = cblas_dnrm2(3,normal,1);
    cblas_dscal(3,1./nrm2,normal,1);
    //set direction
    Range t3;
    rval = mField.get_moab().get_adjacencies(&*fit,1,3,false,t3,Interface::UNION); CHKERR_PETSC(rval);
    t3 = intersect(t3,t1_side);
    if(t3.size()!=1) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    int side_number;
    int sense;
    int offset;
    rval = mField.get_moab().side_number(t3[0],*fit,side_number,sense,offset); CHKERR_PETSC(rval);
    if(sense == -1) {
      cblas_dscal(3,-1,normal,1);
    }
  }

  //build tree
  EntityHandle kdTree_rootMeshset;
  rval = mField.get_moab().create_meshset(MESHSET_SET,kdTree_rootMeshset); CHKERR_PETSC(rval);
  AdaptiveKDTree kdTree(&mField.get_moab());
  rval = kdTree.build_tree(crack_surface_tris,&kdTree_rootMeshset); CHKERR_PETSC(rval);

  Tag th_distance;
  if(verb>0) {
    double def_VAL[1] = { 0 };
    rval = mField.get_moab().tag_get_handle(
      "DISTANCE_FROM_CRACK_SURFACE",1,MB_TYPE_DOUBLE,th_distance,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL
    ); CHKERR_PETSC(rval);
  }

  // calculate distance from crack surface
  signedDistanceMap.clear();
  normalNodesMap.clear();
  for(
    Range::iterator nit = edges_on_surface_nodes.begin();nit!=edges_on_surface_nodes.end();nit++
  ) {
    double coords[3];
    double closest_point[3];
    EntityHandle triangle;
    rval = mField.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
    rval = kdTree.closest_triangle(kdTree_rootMeshset,coords,closest_point,triangle); CHKERR_PETSC(rval);
    double delta[3];
    cblas_dcopy(3,closest_point,1,delta,1);
    cblas_daxpy(3,-1,coords,1,delta,1);
    map<EntityHandle,vector<double> >::iterator mit = normal_map.find(triangle);
    if(mit == normal_map.end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    double *normal = &*mit->second.begin();
    double dot = cblas_ddot(3,normal,1,delta,1);
    signedDistanceMap[*nit] = dot;
    normalNodesMap[*nit] = mit->second;
    if(verb>0) {
      rval = mField.get_moab().tag_set_data(th_distance,&*nit,1,&dot); CHKERR_PETSC(rval);
    }
  }

  rval = mField.get_moab().delete_entities(&kdTree_rootMeshset,1); CHKERR_PETSC(rval);
  //cerr << "normalNodesMap.size() " << normalNodesMap.size() << endl;

  // get nodes near to crack surface extension
  nodesToMoveMap.clear();
  map<EntityHandle,vector<double> >::iterator mit = normalNodesMap.begin();
  for(;mit!=normalNodesMap.end();mit++) {
    EntityHandle node = mit->first;
    double *normal = &*mit->second.begin();
    Range node_edges;
    rval = mField.get_moab().get_adjacencies(&node,1,1,false,node_edges,Interface::UNION); CHKERR_PETSC(rval);
    double max_l = 0;
    for(Range::iterator eit = node_edges.begin();eit!=node_edges.end();eit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
      double coords[3*num_nodes];
      rval = mField.get_moab().get_coords(conn,num_nodes,coords); CHKERR_PETSC(rval);
      cblas_daxpy(3,-1,&coords[3],1,coords,1);
      double l = fabs(cblas_ddot(3,normal,1,coords,1));
      max_l = fmax(l,max_l);
    }
    double dot = signedDistanceMap[node];
    if(fabs(dot/max_l) < 0.125) {
      nodesToMoveMap[node].resize(3);
      double *coords = &*nodesToMoveMap[node].begin();
      rval = mField.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
      cblas_daxpy(3,dot,normal,1,coords,1);
      signedDistanceMap[node] = 0;
      if(verb>0) {
        double dot0 = 0;
        rval = mField.get_moab().tag_set_data(th_distance,&node,1,&dot0); CHKERR_PETSC(rval);
      }
    }
  }


  Range edges_to_split;
  edgesToSplitMap.clear();
  mit = normalNodesMap.begin();
  for(;mit!=normalNodesMap.end();mit++) {
    EntityHandle node = mit->first;
    double dot = signedDistanceMap[node];
    if(dot>=0) continue;
    Range node_edges;
    rval = mField.get_moab().get_adjacencies(
      &node,1,1,false,node_edges,Interface::UNION
    ); CHKERR_PETSC(rval);
    node_edges = intersect(node_edges,edges_to_cut);
    //cerr << node_edges.size() << endl;
    for(Range::iterator eit = node_edges.begin();eit!=node_edges.end();eit++) {
      int num_nodes;
      const EntityHandle* conn;
      rval = mField.get_moab().get_connectivity(*eit,conn,num_nodes,true); CHKERR_PETSC(rval);
      EntityHandle node2 = (conn[0] == node) ? conn[1] : conn[0];
      if(signedDistanceMap.find(node2) == signedDistanceMap.end()) {
        continue;
      }
      double dot2 = signedDistanceMap[node2];
      //cerr << "\t" << dot*dot2 << endl;
      if(dot*dot2 >= 0) {
        continue;
      }
      double coords[3];
      rval = mField.get_moab().get_coords(&node,1,coords); CHKERR_PETSC(rval);
      double coords2[3];
      rval = mField.get_moab().get_coords(&node2,1,coords2); CHKERR_PETSC(rval);
      cblas_daxpy(3,-1,coords,1,coords2,1);
      //double l = cblas_dnrm2(3,coords2,1);
      //dot2 = dot + a*l
      //a = (dot2-dot)/l
      //x = -dot*l/(dot2-dot)
      double x = -dot/(dot2-dot);
      cblas_daxpy(3,x,coords2,1,coords,1);
      edgesToSplitMap[*eit].resize(3);
      cblas_dcopy(3,coords,1,&*edgesToSplitMap[*eit].begin(),1);
      edges_to_split.insert(*eit);
    }
  }

  if(!edges_to_split.empty()) {

    //split edges
    Range level_tets;
    ierr = mField.get_entities_by_type_and_ref_level(bit_mesh,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
    int last_ref_bit = meshRefineBitLevels.back();
    if(!meshIntefaceBitLevels.empty()) {
      if(last_ref_bit<meshIntefaceBitLevels.back()) {
        last_ref_bit = meshIntefaceBitLevels.back();
      }
    }
    last_ref_bit++;
    meshRefineBitLevels.push_back(last_ref_bit);
    BitRefLevel new_ref = BitRefLevel().set(meshRefineBitLevels.back());
    ierr = mField.query_interface(rEfiner); CHKERRQ(ierr);
    ierr = rEfiner->add_verices_in_the_middel_of_edges(edges_to_split,new_ref,2); CHKERRQ(ierr);
    ierr = mField.seed_finite_elements(level_tets); CHKERRQ(ierr);
    ierr = rEfiner->refine_TET(level_tets,new_ref,false); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,cubit_it)) {
      EntityHandle meshset = cubit_it->meshset;
      ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(meshset,new_ref,meshset,MBTET,true); CHKERRQ(ierr);
    }
    ierr = addCrackFront_to_Cubit201(verb); CHKERRQ(ierr);

    for(
      Range::iterator eit = edges_to_split.begin();eit!=edges_to_split.end();eit++
    ) {
      const RefMoFEMEntity_multiIndex *refinedEntitiesPtr_ptr;
      ierr = mField.get_ref_ents(&refinedEntitiesPtr_ptr); CHKERRQ(ierr);
      RefMoFEMEntity_multiIndex::index<Composite_Ent_And_ParentEntType_mi_tag>::type::iterator it;
      it = refinedEntitiesPtr_ptr->get<Composite_Ent_And_ParentEntType_mi_tag>().find(boost::make_tuple(*eit,MBVERTEX));
      if(it == refinedEntitiesPtr_ptr->get<Composite_Ent_And_ParentEntType_mi_tag>().end()) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
      }
      EntityHandle node = it->get_ref_ent();
      nodesToMoveMap[node] = edgesToSplitMap[*eit];
      //rval = mField.get_moab().set_coords(&node,1,&*edgesToSplitMap[*eit].begin()); CHKERR_PETSC(rval);
      new_nodes.insert(node);
    }

  }

  {

    Ranges nodes_to_move;
    for(
      map<EntityHandle,vector<double> >::iterator mit = nodesToMoveMap.begin();
      mit!=nodesToMoveMap.end();mit++
    ) {
      nodes_to_move.insert(mit->first);
    }
    Range nodes_to_move_tris;
    rval = mField.get_moab().get_adjacencies(
      nodes_to_move,2,false,nodes_to_move_tris,Interface::UNION
    ); CHKERR_PETSC(rval);



  }

  PetscFunctionReturn(0);
}*/
