/** \file PostProcOnRefMesh.hpp
 * \brief Postprocess fields on refined mesh made for 10 Node tets
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

/** \brief Post processing
  * \ingroup mofem_fs_post_proc
  */
struct PostProcVolumeOnRefinedMesh: public VolumeElementForcesAndSourcesCore {

  moab::Core coreMesh;
  Interface &postProcMesh;

  bool tenNodesPostProcTets;
  int nbOfRefLevels;

  PostProcVolumeOnRefinedMesh(FieldInterface &m_field,
    bool ten_nodes_post_proc_tets = true,
    int nb_ref_levels = -1
  ):
  VolumeElementForcesAndSourcesCore(m_field),
  postProcMesh(coreMesh),
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
  vector<EntityHandle> mapGaussPts;

  // Gauss pts set on refined mesh
  int getRule(int order) { return -1; };

  struct CommonData {
    Range tEts;
    map<string,vector<ublas::vector<double> > > fieldMap;
    map<string,vector<ublas::matrix<double> > > gradMap;
  };
  CommonData commonData;

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
  on postprocessing mesh.

  */
  PetscErrorCode setGaussPts(int order);

  struct OpHdivFunctions: public VolumeElementForcesAndSourcesCore::UserDataOperator {

    Interface &postProcMesh;
    vector<EntityHandle> &mapGaussPts;

    OpHdivFunctions(
      Interface &post_proc_mesh,
      vector<EntityHandle> &map_gauss_pts,
      const string field_name):
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

  struct OpGetFieldValues: public VolumeElementForcesAndSourcesCore::UserDataOperator {

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
    VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
    postProcMesh(post_proc_mesh),mapGaussPts(map_gauss_pts),
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

  struct OpGetFieldGradientValues: public VolumeElementForcesAndSourcesCore::UserDataOperator {

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
      Vec v = PETSC_NULL):
      VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPCOL),
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

  /** \brief Add operator to postprocess Hdiv field
  */
  PetscErrorCode addHdivFunctionsPostProc(const string field_name);

  /** \brief Add operator to postprocess L2 or H1 field value

    \param field_name
    \param v If vector is given, values from vector are used to set tags on mesh

    Note:
    Name of the tag to store values on postprocess mesh is the same as field name

  */
  PetscErrorCode addFieldValuesPostProc(const string field_name,Vec v = PETSC_NULL);

  /** \brief Add operator to postprocess L2 or H1 field value

    \param field_name
    \param tag_name to store results on potprocess mesh
    \param v If vector is given, values from vector are used to set tags on mesh

  */
  PetscErrorCode addFieldValuesPostProc(const string field_name,const string tag_name,Vec v = PETSC_NULL);


  /** \brief Add operator to postprocess L2 or H1 field gradient

    \param field_name
    \param v If vector is given, values from vector are used to set tags on mesh

    Note:
    Name of the tag to store values on postprocess mesh is the same as field name

  */
  PetscErrorCode addFieldValuesGradientPostProc(const string field_name,Vec v = PETSC_NULL);

  /** \brief Add operator to postprocess L2 or H1 field gradient

    \param field_name
    \param tag_name to store results on potprocess mesh
    \param v If vector is given, values from vector are used to set tags on mesh

  */
  PetscErrorCode addFieldValuesGradientPostProc(const string field_name,const string tag_name,Vec v = PETSC_NULL);

  /** \brief Clear operators list

  Clear operators list, user can use the same mesh instance to postpocess
  different problem or the same problem with different set of postprocessed
  fields.

  */
  PetscErrorCode clearOperators();

  PetscErrorCode preProcess();

  PetscErrorCode postProcess();


};

DEPRECATED typedef PostProcVolumeOnRefinedMesh PostPocOnRefinedMesh;

#endif //__POSTPROC_ON_REF_MESH_HPP

/***************************************************************************//**
 * \defgroup mofem_fs_post_proc Post Process
 * \ingroup user_modules
 ******************************************************************************/
