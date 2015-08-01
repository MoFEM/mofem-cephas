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

struct FaceSplittingTools {

  FieldInterface& mField;
  MeshRefinment* rEfiner;
  PrismInterface* prismInterface;

  FaceSplittingTools(FieldInterface& m_field):
    mField(m_field) {

    ErrorCode rval;

    int def_bit_level_vec[BITREFLEVEL_SIZE];
    bzero(def_bit_level_vec,BITREFLEVEL_SIZE*sizeof(int));
    rval = mField.get_moab().tag_get_handle(
      "_MESHREFINEBITLEVELS",BITREFLEVEL_SIZE*sizeof(int),MB_TYPE_OPAQUE,
      th_meshRefineBitLevels,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_vec);  CHKERR(rval);
    rval = mField.get_moab().tag_get_handle(
      "_MESHINTEFACEBITLEVELS",BITREFLEVEL_SIZE*sizeof(int),MB_TYPE_OPAQUE,
      th_meshIntefaceBitLevels,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_vec); CHKERR(rval);
    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    rval = mField.get_moab().tag_get_by_ptr(th_meshRefineBitLevels,&root_meshset,1,(const void**)&ptr_meshRefineBitLevels); CHKERR(rval);
    rval = mField.get_moab().tag_get_by_ptr(th_meshIntefaceBitLevels,&root_meshset,1,(const void**)&ptr_meshIntefaceBitLevels); CHKERR(rval);
    meshRefineBitLevels.ptr = ptr_meshRefineBitLevels;
    meshIntefaceBitLevels.ptr = ptr_meshIntefaceBitLevels;

    EntityHandle def_meshset = 0;
    rval = mField.get_moab().tag_get_handle(
      "_TETGEN_GEOMETRY",1,MB_TYPE_HANDLE,
      th_Polygons,MB_TAG_CREAT|MB_TAG_MESH,&def_meshset); CHKERR(rval);
    rval = mField.get_moab().tag_get_data(th_Polygons,&root_meshset,1,&polygonsMeshset); CHKERR(rval);
    if(polygonsMeshset == 0) {
      rval = mField.get_moab().create_meshset(MESHSET_SET,polygonsMeshset); CHKERR(rval);
      rval = mField.get_moab().tag_set_data(th_Polygons,&root_meshset,1,&polygonsMeshset); CHKERR(rval);
    }

  }

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

  PetscErrorCode meshRefine(const int verb = 0);
  PetscErrorCode splitFaces(const int verb = 0);

  //add edges to crack front
  PetscErrorCode addCrackFront_to_Cubit201(int verb = 0);
  PetscErrorCode roundCornersFillGaps_in_Cubit200(int nb,int verb = 0);

  //calculate length of edges adjacent to crack front
  PetscErrorCode crackFrontEdgeLengths(
    BitRefLevel bit_mesh,Range &to_split,Range &to_remove,int verb = 0
  );

  //move front nodes
  PetscErrorCode moveFrontNodesByVec(double v[]);

  EntityHandle polygonsMeshset;
  Tag th_Polygons;

  #ifdef WITH_TETGEN

  map<EntityHandle,unsigned long> moabTetGenMap;
  map<unsigned long,EntityHandle> tetGenMoabMap;

  boost::ptr_vector<tetgenio> tetGenData;
  PetscErrorCode rebuildMeshWithTetGen(
    vector<string> &switches,const int verb = -1);


  #endif

  PetscErrorCode getCornerEdges(Range &edges_to_cat,int verb = 0);
  PetscErrorCode propagateBySplit(Range &new_nodes,Range &edges_to_cat,int verb = 0);
  PetscErrorCode cornerProblem(Range &new_nodes,int verb = 0);

};


#endif // __FACESPLITTINGTOOL_HPP__
