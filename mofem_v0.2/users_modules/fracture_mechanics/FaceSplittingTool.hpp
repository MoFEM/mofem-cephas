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

#ifndef __FACESPLITTINGTOOL_HPP__
#define __FACESPLITTINGTOOL_HPP__

namespace MoFEM {

struct FaceSplittingTools {

  FieldInterface& mField;
  MeshRefinment* rEfiner;
  PrismInterface* prismInterface; 
  
  Interface& moab_distance_from_crack_surface;
  moab::Core mb_instance_distance_from_crack_surface;
  AdaptiveKDTree kdTree_DistanceFromCrackSurface;

  FaceSplittingTools(FieldInterface& _mField): 
    mField(_mField),
    moab_distance_from_crack_surface(mb_instance_distance_from_crack_surface),
    kdTree_DistanceFromCrackSurface(&moab_distance_from_crack_surface),
    th_b(NULL),th_distance(NULL),th_projection(NULL) {

    kdTree_rootMeshset_DistanceFromCrackSurface = 0;
    opositeFrontEdges = 0;
    crackFrontTests = 0;
    chopTetsFaces = 0;
    selectedCrackFaces = 0;

    int def_bit_level_vec[BITREFLEVEL_SIZE];
    bzero(def_bit_level_vec,BITREFLEVEL_SIZE*sizeof(int));
    mField.get_moab().tag_get_handle(
      "_MESHREFINEBITLEVELS",BITREFLEVEL_SIZE*sizeof(int),MB_TYPE_OPAQUE,
      th_meshRefineBitLevels,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_vec); 
    mField.get_moab().tag_get_handle(
      "_MESHINTEFACEBITLEVELS",BITREFLEVEL_SIZE*sizeof(int),MB_TYPE_OPAQUE,
      th_meshIntefaceBitLevels,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_vec); 
    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    mField.get_moab().tag_get_by_ptr(th_meshRefineBitLevels,&root_meshset,1,(const void**)&ptr_meshRefineBitLevels);
    mField.get_moab().tag_get_by_ptr(th_meshIntefaceBitLevels,&root_meshset,1,(const void**)&ptr_meshIntefaceBitLevels);
    meshRefineBitLevels.ptr = ptr_meshRefineBitLevels;
    meshIntefaceBitLevels.ptr = ptr_meshIntefaceBitLevels;

    mField.query_interface(rEfiner);
    mField.query_interface(prismInterface);

  }

  ~FaceSplittingTools() {
    cleanMeshsets();
    /*if(kdTree_rootMeshset_DistanceFromCrackSurface!=0) {
      mField.get_moab().delete_entities(&kdTree_rootMeshset_DistanceFromCrackSurface,1);
      kdTree_rootMeshset_DistanceFromCrackSurface = 0;
    }*/
  }
 
  PetscErrorCode cleanMeshsets() {
    PetscFunctionBegin;

    if(opositeFrontEdges!=0) {
      rval = mField.get_moab().delete_entities(&opositeFrontEdges,1); CHKERR_PETSC(rval);
      opositeFrontEdges = 0;
    }
    if(crackFrontTests!=0) {
      rval = mField.get_moab().delete_entities(&crackFrontTests,1); CHKERR_PETSC(rval);
      crackFrontTests = 0;
    }
    if(chopTetsFaces!=0) {
      rval = mField.get_moab().delete_entities(&chopTetsFaces,1); CHKERR_PETSC(rval);
      chopTetsFaces=0;
    }
    if(selectedCrackFaces!=0) {
      rval = mField.get_moab().delete_entities(&selectedCrackFaces,1); CHKERR_PETSC(rval);
      selectedCrackFaces=0;
    }

    if(th_b != NULL) {
      rval = mField.get_moab().tag_delete(th_b); CHKERR_PETSC(rval);
      th_b = NULL;
    }

    if(th_distance != NULL) {
      rval = mField.get_moab().tag_delete(th_distance); CHKERR_PETSC(rval);
      th_distance = NULL;
    }

    if(th_projection != NULL) {
      rval = mField.get_moab().tag_delete(th_projection); CHKERR_PETSC(rval);
      th_projection = NULL;
    }
   
    PetscFunctionReturn(0);
  }

  //Cloasest face to crack surface, crack surface perspective 

  EntityHandle kdTree_rootMeshset_DistanceFromCrackSurface;
  map<EntityHandle,EntityHandle> map_nodes;
  PetscErrorCode buildKDTreeForCrackSurface(
    Range &entities,const BitRefLevel bit_mesh);
  PetscErrorCode buildKDTreeForCrackSurface(const BitRefLevel bit_mesh);

  //Init bit level data
  Range mesh_level_nodes;
  Range mesh_level_edges;
  Range mesh_level_tris;
  Range mesh_level_tets;
  PetscErrorCode initBitLevelData(const BitRefLevel bit_mesh);

  //Calulte distance on mesh
  PetscErrorCode calculateDistanceFromCrackSurface(Range &nodes,double alpha,int verb = 0);
  PetscErrorCode calculateDistanceCrackFrontNodesFromCrackSurface(double alpha,int verb = 0);
  PetscErrorCode calculateDistanceFromCrackSurface(int verb = 0);

  //Front edges

  EntityHandle opositeFrontEdges;

  PetscErrorCode getOpositeForntEdges(bool createMeshset,int verb = 0);

  //Front tets

  EntityHandle crackFrontTests;
  EntityHandle chopTetsFaces;
  EntityHandle selectedCrackFaces;

  PetscErrorCode getCrackFrontEntities(bool createMeshset,bool get_tets,bool get_faces,int verb = 0);
  PetscErrorCode getCrackFrontFaces(bool createMeshset,int verb = 0);
  PetscErrorCode getCrackFrontTets(bool createMeshset,int verb = 0);
  PetscErrorCode chopTetsUntilNonOneLeftOnlyCrackSurfaceFaces(bool createMeshset,int verb = 0);
  PetscErrorCode selectCrackFaces(bool createMeshset,int verb = 0);

  //Split new crack front faces
  Tag th_meshRefineBitLevels,th_meshIntefaceBitLevels;
  int *ptr_meshRefineBitLevels,*ptr_meshIntefaceBitLevels;

  struct  BitRefLevelVector {
    int* ptr;
    bool empty() { return !((bool)ptr[0]); }
    int size() { return ptr[0]; }
    void resize(int s) { ptr[0] = s; }
    int& first() { return ptr[1]; }
    int& back() { return ptr[ptr[0]]; }
    int* begin() { return &ptr[1]; };
    int* end() { return &ptr[ptr[0]+1]; }
    void push_back(int a) { 
      ptr[0]++;
      ptr[ptr[0]] = a; 
    }
  };

  BitRefLevelVector meshRefineBitLevels;
  BitRefLevelVector meshIntefaceBitLevels;

  PetscErrorCode catMesh(const int verb = -1);
  PetscErrorCode meshRefine(const int verb = -1);
  PetscErrorCode splitFaces(const int verb = -1);

  PetscErrorCode addNewSurfaceFaces_to_Cubit_msId200();
  PetscErrorCode addcrackFront_to_Cubit201();

  /** \brief project coords ibn crack surfec
   * do not use this function, it not respect body boundaries
   */
  PetscErrorCode projectCrackFrontNodes();

  PetscErrorCode getMask(BitRefLevel &maskPreserv,const int verb = -1);
  PetscErrorCode squashIndices(const int verb = -1);

  private:
  ErrorCode rval;
  PetscErrorCode ierr;

  double diffNTET[4*3];
  Tag th_b;
  Tag th_distance;
  Tag th_projection;

  PetscErrorCode calculate_qualityAfterProjectingNodes(
    Range &option_nodes,Range &intersect_tets,
    Range &crack_front_edges_nodes,double &current_q);

};

PetscErrorCode main_refine_and_meshcat(FieldInterface& mField,FaceSplittingTools &face_splitting,const int verb = 0);
PetscErrorCode main_select_faces_for_splitting(FieldInterface& mField,FaceSplittingTools &face_splitting,const int verb = 0);
PetscErrorCode main_split_faces_and_update_field_and_elements(FieldInterface& mField,FaceSplittingTools &face_splitting,const int verb = 0);

}

#endif // __FACESPLITTINGTOOL_HPP__

