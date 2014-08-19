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
  
  FaceSplittingTools(FieldInterface& _mField): 
    mField(_mField) {

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

  }

  ~FaceSplittingTools() {
    cleanMeshsets();
  }
 
  PetscErrorCode cleanMeshsets() {
    PetscFunctionBegin;

    PetscFunctionReturn(0);
  }

  //Init bit level data
  Range mesh_level_nodes;
  Range mesh_level_edges;
  Range mesh_level_tris;
  Range mesh_level_tets;
  PetscErrorCode initBitLevelData(const BitRefLevel bit_mesh);

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

  PetscErrorCode meshRefine(const int verb = -1);
  PetscErrorCode splitFaces(const int verb = -1);



};

}

#endif // __FACESPLITTINGTOOL_HPP__

