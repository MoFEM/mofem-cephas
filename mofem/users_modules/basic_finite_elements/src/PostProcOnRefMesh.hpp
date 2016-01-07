/** \file PostProcOnRefMesh.hpp
 * \brief Post-process fields on refined mesh
 *
 * Create refined mesh, without enforcing continuity between element. Calculate
 * field values on nodes of that mesh.
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

#ifndef __POSTPROC_ON_REF_MESH_HPP
#define __POSTPROC_ON_REF_MESH_HPP

struct PostProcCommonOnRefMesh {

  struct CommonData {
    map<string,vector<ublas::vector<double> > > fieldMap;
    map<string,vector<ublas::matrix<double> > > gradMap;
  };

  struct OpGetFieldValues: public ForcesAndSurcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;
    CommonData &commonData;
    const string tagName;
    Vec V;

    OpGetFieldValues(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name,
      const string tag_name,
      CommonData &common_data,
      Vec v = PETSC_NULL
    ):
    ForcesAndSurcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts),
    commonData(common_data),
    tagName(tag_name),
    V(v) {}

    ublas::vector<double> vAlues;
    ublas::vector<double> *vAluesPtr;

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    );

  };

  struct OpGetFieldGradientValues: public ForcesAndSurcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;
    CommonData &commonData;
    const string tagName;
    Vec V;

    OpGetFieldGradientValues(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name,
      const string tag_name,
      CommonData &common_data,
      Vec v = PETSC_NULL
    ):
    ForcesAndSurcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts),
    commonData(common_data),
    tagName(tag_name),
    V(v)
    {}

    ublas::vector<double> vAlues;
    ublas::vector<double> *vAluesPtr;

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    );

  };

};

template<class ELEMENT>
struct PostProcTemplateOnRefineMesh: public ELEMENT {

  moab::Core coreMesh;
  Interface &postProcMesh;
  vector<EntityHandle> mapGaussPts;

  PostProcTemplateOnRefineMesh(FieldInterface &m_field):
  ELEMENT(m_field),
  postProcMesh(coreMesh) {
  }

  virtual PostProcCommonOnRefMesh::CommonData& getCommonData() {
    THROW_AT_LINE("not implemented");
  }

  /** \brief Add operator to post-process L2 or H1 field value

    \param field_name
    \param v If vector is given, values from vector are used to set tags on mesh

    Note:
    Name of the tag to store values on post-process mesh is the same as field name

  */
  PetscErrorCode addFieldValuesPostProc(const string field_name,Vec v = PETSC_NULL) {
    PetscFunctionBegin;
    ELEMENT::getOpPtrVector().push_back(
      new PostProcCommonOnRefMesh::OpGetFieldValues(
        postProcMesh,mapGaussPts,field_name,field_name,getCommonData(),v
      )
    );
    PetscFunctionReturn(0);
  }

  /** \brief Add operator to post-process L2 or H1 field value

    \param field_name
    \param tag_name to store results on post-process mesh
    \param v If vector is given, values from vector are used to set tags on mesh

  */
  PetscErrorCode addFieldValuesPostProc(const string field_name,const string tag_name,Vec v = PETSC_NULL) {
    PetscFunctionBegin;
    ELEMENT::getOpPtrVector().push_back(
      new PostProcCommonOnRefMesh::OpGetFieldValues(
        postProcMesh,mapGaussPts,field_name,tag_name,getCommonData(),v
      )
    );
    PetscFunctionReturn(0);
  }


  /** \brief Add operator to post-process L2 or H1 field gradient

    \param field_name
    \param v If vector is given, values from vector are used to set tags on mesh

    Note:
    Name of the tag to store values on post-process mesh is the same as field name

  */
  PetscErrorCode addFieldValuesGradientPostProc(const string field_name,Vec v = PETSC_NULL) {
    PetscFunctionBegin;
    ELEMENT::getOpPtrVector().push_back(
      new PostProcCommonOnRefMesh::OpGetFieldGradientValues(
        postProcMesh,mapGaussPts,field_name,field_name+"_GRAD",getCommonData(),v
      )
    );
    PetscFunctionReturn(0);
  }

  /** \brief Add operator to post-process L2 or H1 field gradient

    \param field_name
    \param tag_name to store results on post-process mesh
    \param v If vector is given, values from vector are used to set tags on mesh

  */
  PetscErrorCode addFieldValuesGradientPostProc(const string field_name,const string tag_name,Vec v = PETSC_NULL) {
    PetscFunctionBegin;
    ELEMENT::getOpPtrVector().push_back(
      new PostProcCommonOnRefMesh::OpGetFieldGradientValues(
        postProcMesh,mapGaussPts,field_name,tag_name,getCommonData(),v
      )
    );
    PetscFunctionReturn(0);
  }

};

/** \brief Post processing
  * \ingroup mofem_fs_post_proc
  */
struct PostProcVolumeOnRefinedMesh: public PostProcTemplateOnRefineMesh<VolumeElementForcesAndSourcesCore> {

  bool tenNodesPostProcTets;
  int nbOfRefLevels;

  PostProcVolumeOnRefinedMesh(
    FieldInterface &m_field,
    bool ten_nodes_post_proc_tets = true,
    int nb_ref_levels = -1
  ):
  PostProcTemplateOnRefineMesh<VolumeElementForcesAndSourcesCore>(m_field),
  tenNodesPostProcTets(ten_nodes_post_proc_tets),
  nbOfRefLevels(nb_ref_levels) {
  }

  virtual ~PostProcVolumeOnRefinedMesh() {
    ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
    if(pcomm_post_proc_mesh != NULL) {
      delete pcomm_post_proc_mesh;
    }
  }

  ublas::matrix<int> refTets;
  ublas::matrix<double> gaussPts_FirstOrder;

  // Gauss pts set on refined mesh
  int getRule(int order) { return -1; };

  struct CommonData: PostProcCommonOnRefMesh::CommonData {
    Range tEts;
  };
  CommonData commonData;

  virtual PostProcCommonOnRefMesh::CommonData& getCommonData() {
    return commonData;
  }

  /** \brief Generate reference mesh on single element

  Each element is subdivided on smaller elements, i.e. a reference mesh on
  single element is created. Nodes of such reference mesh are used as
  integration points at which field values are calculated and to
  each node a "moab" tag is attached to store those values.

  */
  PetscErrorCode generateReferenceElementMesh();

  /** \brief Set integration points

  If reference mesh is generated on single elements. This function maps
  reference coordinates into physical coordinates and create element
  on post-processing mesh.

  */
  PetscErrorCode setGaussPts(int order);


  /** \brief Clear operators list

  Clear operators list, user can use the same mesh instance to post-process
  different problem or the same problem with different set of post-processed
  fields.

  */
  PetscErrorCode clearOperators();

  PetscErrorCode preProcess();

  PetscErrorCode postProcess();

  /** \brief Add operator to post-process Hdiv field
  */
  PetscErrorCode addHdivFunctionsPostProc(const string field_name);

  struct OpHdivFunctions: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;

    OpHdivFunctions(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name
    ):
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    postProcMesh(post_proc_mesh),
    mapGaussPts(map_gauss_pts) {
    }

    PetscErrorCode doWork(
      int side,
      EntityType type,
      DataForcesAndSurcesCore::EntData &data
    );

  };

};

DEPRECATED typedef PostProcVolumeOnRefinedMesh PostPocOnRefinedMesh;

struct PostProcFatPrismOnRefinedMesh: public PostProcTemplateOnRefineMesh<FatPrismElementForcesAndSurcesCore> {

  bool tenNodesPostProcTets;

  PostProcFatPrismOnRefinedMesh(
    FieldInterface &m_field,
    bool ten_nodes_post_proc_tets = true
  ):
  PostProcTemplateOnRefineMesh<FatPrismElementForcesAndSurcesCore>(m_field),
  tenNodesPostProcTets(ten_nodes_post_proc_tets) {
  }

  virtual ~PostProcFatPrismOnRefinedMesh() {
    ParallelComm* pcomm_post_proc_mesh = ParallelComm::get_pcomm(&postProcMesh,MYPCOMM_INDEX);
    if(pcomm_post_proc_mesh != NULL) {
      delete pcomm_post_proc_mesh;
    }
  }

  int getRuleTrianglesOnly(int order) { return -1; };
  int getRuleThroughThickness(int order) { return -1; };

  PetscErrorCode setGaussPtsTrianglesOnly(int order_triangles_only);
  PetscErrorCode setGaussPtsThroughThickness(int order_thickness);
  PetscErrorCode generateReferenceElementMesh();

  map<EntityHandle,EntityHandle> elementsMap;

  PetscErrorCode preProcess();
  PetscErrorCode postProcess();

  struct CommonData: PostProcCommonOnRefMesh::CommonData {
  };
  CommonData commonData;

  virtual PostProcCommonOnRefMesh::CommonData& getCommonData() {
    return commonData;
  }

};

#endif //__POSTPROC_ON_REF_MESH_HPP

/***************************************************************************//**
 * \defgroup mofem_fs_post_proc Post Process
 * \ingroup user_modules
 ******************************************************************************/
