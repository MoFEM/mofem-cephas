/** \file PostProc.hpp
 * \brief Post-process fields on refined mesh
 *
 * Methods for post-processing finite elements
 *
 * \ingroup mofem_fs_post_proc
 */

#ifndef __POST_PROC_HPP
#define __POST_PROC_HPP

namespace MoFEM {

/**
 * Each element is subdivided on smaller elements, i.e. a reference mesh on
 * single element is created. Nodes of such reference mesh are used as
 * integration points at which field values are calculated and to
 * each node a "moab" tag is attached to store those values.
 */
struct PostProcGenerateRefMeshBase;

using PostProcGenerateRefMeshPtr =
    boost::shared_ptr<PostProcGenerateRefMeshBase>;

/**
 * @brief Element for postprocessing. Uses MoAB to generate post-processing
 * mesh.
 *
 * @tparam T Finite Element Implementation
 */
template <EntityType T> struct PostProcGenerateRefMesh;

template <typename E> struct PostProcBrokenMeshInMoabBase : public E {

  PostProcBrokenMeshInMoabBase(MoFEM::Interface &m_field);
  virtual ~PostProcBrokenMeshInMoabBase();

  std::vector<EntityHandle> mapGaussPts;
  Range postProcElements;

  /**
   * @brief Get vector of vectors associated to integration points 
   * 
   * @return std::vector<EntityHandle>& 
   */
  inline auto &getMapGaussPts();

  /**
   * @brief Get postprocessing mesh
   *
   * @return moab::Interface&
   */
  inline auto &getPostProcMesh();

  /**
   * @brief Get postprocessing elements
   * 
   * @return auto& 
   */
  inline auto &getPostProcElements();

  /**
 * \brief wrote results in (MOAB) format, use "file_name.h5m"
 * @param  file_name file name (should always end with .h5m)
 * @return           error code

 * \ingroup mofem_fs_post_proc

 */
  MoFEMErrorCode writeFile(const std::string file_name);

  /**
   * @brief Generate vertices and elements
   *
   * @return MoFEMErrorCode
   */
  MoFEMErrorCode preProcess();

  MoFEMErrorCode postProcess();

  MoFEMErrorCode setGaussPts(int order);

protected:
  int getRule(int order);
  int getMaxLevel() const;

  moab::Core coreMesh;
  moab::Interface &postProcMesh;

  std::map<EntityType, PostProcGenerateRefMeshPtr> refElementsMap;
};

template <typename E> auto &PostProcBrokenMeshInMoabBase<E>::getMapGaussPts() {
  return mapGaussPts;
}

template <typename E> auto &PostProcBrokenMeshInMoabBase<E>::getPostProcMesh() {
  return postProcMesh;
}

template <typename E>
inline auto &PostProcBrokenMeshInMoabBase<E>::getPostProcElements() {
  return postProcElements;
}

template <typename E>
MoFEMErrorCode
PostProcBrokenMeshInMoabBase<E>::writeFile(const std::string file_name) {
  MoFEMFunctionBegin;
  ParallelComm *pcomm_post_proc_mesh =
      ParallelComm::get_pcomm(&postProcMesh, MYPCOMM_INDEX);
  if (pcomm_post_proc_mesh == NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "ParallelComm not allocated");
  CHKERR postProcMesh.write_file(file_name.c_str(), "MOAB",
                                 "PARALLEL=WRITE_PART");
  MoFEMFunctionReturn(0);
};

template <typename E> struct PostProcBrokenMeshInMoab;

template <>
struct PostProcBrokenMeshInMoab<VolumeElementForcesAndSourcesCore>
    : public PostProcBrokenMeshInMoabBase<VolumeElementForcesAndSourcesCore> {
  using PostProcBrokenMeshInMoabBase<
      VolumeElementForcesAndSourcesCore>::PostProcBrokenMeshInMoabBase;
};

template <>
struct PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>
    : public PostProcBrokenMeshInMoabBase<FaceElementForcesAndSourcesCore> {
  using PostProcBrokenMeshInMoabBase<
      FaceElementForcesAndSourcesCore>::PostProcBrokenMeshInMoabBase;
};

template <>
struct PostProcBrokenMeshInMoab<EdgeElementForcesAndSourcesCore>
    : public PostProcBrokenMeshInMoabBase<EdgeElementForcesAndSourcesCore> {
  using PostProcBrokenMeshInMoabBase<
      EdgeElementForcesAndSourcesCore>::PostProcBrokenMeshInMoabBase;
};

template <typename E>
boost::shared_ptr<E> make_post_proc_fe_in_moab(MoFEM::Interface &m_field);

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<VolumeElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field);

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<FaceElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field);

template <>
boost::shared_ptr<PostProcBrokenMeshInMoab<EdgeElementForcesAndSourcesCore>>
make_post_proc_fe_in_moab(MoFEM::Interface &m_field);

/**
 * @brief Post post-proc data at points from hash maps
 *
 * @tparam DIM1 dimension of vector in data_map_vec  and first column of
 * data_map_may
 * @tparam DIM2 dimension of second column in data_map_mat
 */
template <int DIM1, int DIM2>
struct OpPostProcMapInMoab : public ForcesAndSourcesCore::UserDataOperator {

  using DataMapVec = std::map<std::string, boost::shared_ptr<VectorDouble>>;
  using DataMapMat = std::map<std::string, boost::shared_ptr<MatrixDouble>>;

  /**
   * @brief Construct a new OpPostProcMapInMoab object
   *
   * @param post_proc_mesh postprocessing mesh
   * @param map_gauss_pts map of gauss points to nodes of postprocessing mesh
   * @param data_map_scalar hash map of scalar values (string is name of the
   * tag)
   * @param data_map_vec hash map of vector values
   * @param data_map_mat hash map of second order tensor values
   * @param data_symm_map_mat hash map of symmetric second order tensor values
   */
  OpPostProcMapInMoab(moab::Interface &post_proc_mesh,
                      std::vector<EntityHandle> &map_gauss_pts,
                      DataMapVec data_map_scalar, DataMapMat data_map_vec,
                      DataMapMat data_map_mat, DataMapMat data_symm_map_mat)
      : ForcesAndSourcesCore::UserDataOperator(
            NOSPACE, ForcesAndSourcesCore::UserDataOperator::OPSPACE),
        postProcMesh(post_proc_mesh), mapGaussPts(map_gauss_pts),
        dataMapScalar(data_map_scalar), dataMapVec(data_map_vec),
        dataMapMat(data_map_mat), dataMapSymmMat(data_symm_map_mat) {
    // Operator is only executed for vertices
    std::fill(&doEntities[MBEDGE], &doEntities[MBMAXTYPE], false);
  }
  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data);

private:
  moab::Interface &postProcMesh;
  std::vector<EntityHandle> &mapGaussPts;
  DataMapVec dataMapScalar;
  DataMapMat dataMapVec;
  DataMapMat dataMapMat;
  DataMapMat dataMapSymmMat;
};

template <int DIM1, int DIM2>
MoFEMErrorCode
OpPostProcMapInMoab<DIM1, DIM2>::doWork(int side, EntityType type,
                                        EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  std::array<double, 9> def;
  std::fill(def.begin(), def.end(), 0);

  auto get_tag = [&](const std::string name, size_t size) {
    Tag th;
    CHKERR postProcMesh.tag_get_handle(name.c_str(), size, MB_TYPE_DOUBLE, th,
                                       MB_TAG_CREAT | MB_TAG_SPARSE,
                                       def.data());
    return th;
  };

  MatrixDouble3by3 mat(3, 3);

  auto set_vector_3d = [&](auto &t) -> MatrixDouble3by3 & {
    mat.clear();
    for (size_t r = 0; r != DIM1; ++r)
      mat(0, r) = t(r);
    return mat;
  };

  auto set_matrix_3d = [&](auto &t) -> MatrixDouble3by3 & {
    mat.clear();
    for (size_t r = 0; r != DIM1; ++r)
      for (size_t c = 0; c != DIM2; ++c)
        mat(r, c) = t(r, c);
    return mat;
  };

  auto set_matrix_symm_3d = [&](auto &t) -> MatrixDouble3by3 & {
    mat.clear();
    for (size_t r = 0; r != DIM1; ++r)
      for (size_t c = 0; c != DIM1; ++c)
        mat(r, c) = t(r, c);
    return mat;
  };

  auto set_scalar = [&](auto t) -> MatrixDouble3by3 & {
    mat.clear();
    mat(0, 0) = t;
    return mat;
  };

  auto set_float_precision = [](const double x) {
    if (std::abs(x) < std::numeric_limits<float>::epsilon())
      return 0.;
    else
      return x;
  };

  auto set_tag = [&](auto th, auto gg, MatrixDouble3by3 &mat) {
    for (auto &v : mat.data())
      v = set_float_precision(v);
    return postProcMesh.tag_set_data(th, &mapGaussPts[gg], 1,
                                     &*mat.data().begin());
  };

  for (auto &m : dataMapScalar) {
    auto th = get_tag(m.first, 1);
    auto t_scl = getFTensor0FromVec(*m.second);
    auto nb_integration_pts = getGaussPts().size2();
    size_t gg = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      CHKERR set_tag(th, gg, set_scalar(t_scl));
      ++t_scl;
    }
  }

  for (auto &m : dataMapVec) {
    auto th = get_tag(m.first, 3);
    auto t_vec = getFTensor1FromMat<DIM1>(*m.second);
    auto nb_integration_pts = getGaussPts().size2();
    size_t gg = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      CHKERR set_tag(th, gg, set_vector_3d(t_vec));
      ++t_vec;
    }
  }

  for (auto &m : dataMapMat) {
    auto th = get_tag(m.first, 9);
    auto t_mat = getFTensor2FromMat<DIM1, DIM2>(*m.second);
    auto nb_integration_pts = getGaussPts().size2();
    size_t gg = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      CHKERR set_tag(th, gg, set_matrix_3d(t_mat));
      ++t_mat;
    }
  }

  for (auto &m : dataMapSymmMat) {
    auto th = get_tag(m.first, 9);
    auto t_mat = getFTensor2SymmetricFromMat<DIM1>(*m.second);
    auto nb_integration_pts = getGaussPts().size2();
    size_t gg = 0;
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      CHKERR set_tag(th, gg, set_matrix_symm_3d(t_mat));
      ++t_mat;
    }
  }

  MoFEMFunctionReturn(0);
}
} // namespace MoFEM

#endif //__POST_PROC_HPP