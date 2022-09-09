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
struct PostProcGenerateRefMeshBase {

  std::vector<MatrixDouble> levelShapeFunctions;
  std::vector<MatrixDouble> levelGaussPtsOnRefMesh;
  std::vector<ublas::matrix<int>> levelRef;

  EntityHandle startingVertEleHandle;
  std::vector<double *> verticesOnEleArrays;
  EntityHandle startingEleHandle;
  EntityHandle *eleConn;

  int countEle;
  int countVertEle;

  int nbVertices;
  int nbEles;

  PostProcGenerateRefMeshBase();
  MoFEMErrorCode getOptions();

  virtual MoFEMErrorCode generateReferenceElementMesh() = 0;

  PetscBool hoNodes;
  int defMaxLevel;
  std::string optPrefix;
};

using PostProcGenerateRefMeshPtr =
    boost::shared_ptr<PostProcGenerateRefMeshBase>;

/**
 * @brief Element for postprocessing. Uses MoAB to generate post-processing
 * mesh.
 *
 * @tparam T Finite Element Implementation
 */
template <EntityType T> struct PostProcGenerateRefMesh;

template <>
struct PostProcGenerateRefMesh<MBTET> : public PostProcGenerateRefMeshBase {
  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;
  MoFEMErrorCode generateReferenceElementMesh();
};

template <>
struct PostProcGenerateRefMesh<MBHEX> : public PostProcGenerateRefMeshBase {
  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;
  MoFEMErrorCode generateReferenceElementMesh();
};

template <>
struct PostProcGenerateRefMesh<MBTRI> : public PostProcGenerateRefMeshBase {
  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;
  MoFEMErrorCode generateReferenceElementMesh();
};

template <>
struct PostProcGenerateRefMesh<MBQUAD> : public PostProcGenerateRefMeshBase {
  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;
  MoFEMErrorCode generateReferenceElementMesh();
};

template <>
struct PostProcGenerateRefMesh<MBEDGE> : public PostProcGenerateRefMeshBase {
  using PostProcGenerateRefMeshBase::PostProcGenerateRefMeshBase;
  MoFEMErrorCode generateReferenceElementMesh();
};

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

template <typename E> int PostProcBrokenMeshInMoabBase<E>::getMaxLevel() const {
  auto get_element_max_dofs_order = [&]() {
    int max_order = 0;
    auto dofs_vec = E::getDataVectorDofsPtr();
    for (auto &dof : *dofs_vec) {
      const int dof_order = dof->getDofOrder();
      max_order = (max_order < dof_order) ? dof_order : max_order;
    };
    return max_order;
  };
  const auto dof_max_order = get_element_max_dofs_order();
  return (dof_max_order > 0) ? (dof_max_order - 1) / 2 : 0;
};

template <typename E>
PostProcBrokenMeshInMoabBase<E>::PostProcBrokenMeshInMoabBase(
    MoFEM::Interface &m_field)
    : E(m_field), postProcMesh(coreMesh) {}

template <typename E>
PostProcBrokenMeshInMoabBase<E>::~PostProcBrokenMeshInMoabBase() {
  ParallelComm *pcomm_post_proc_mesh =
      ParallelComm::get_pcomm(&postProcMesh, MYPCOMM_INDEX);
  if (pcomm_post_proc_mesh != NULL)
    delete pcomm_post_proc_mesh;
}

template <typename E> int PostProcBrokenMeshInMoabBase<E>::getRule(int order) {
  return -1;
};

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::setGaussPts(int order) {
  MoFEMFunctionBegin;

  auto type = type_from_handle(this->getFEEntityHandle());

  PostProcGenerateRefMeshPtr ref_ele;

  try {
    ref_ele = refElementsMap.at(type);
  } catch (const out_of_range &e) {
    SETERRQ1(
        PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
        "Generation of reference elements for type <%s> is not implemented",
        moab::CN::EntityTypeName(type));
  }

  auto set_gauss_pts = [&](auto &level_gauss_pts_on_ref_mesh, auto &level_ref,
                           auto &level_shape_functions,

                           auto start_vert_handle, auto start_ele_handle,
                           auto &verts_array, auto &conn, auto &ver_count,
                           auto &ele_count

                       ) {
    MoFEMFunctionBegin;

    size_t level = getMaxLevel();
    level = std::min(level, level_gauss_pts_on_ref_mesh.size() - 1);

    auto &level_ref_gauss_pts = level_gauss_pts_on_ref_mesh[level];
    auto &level_ref_ele = level_ref[level];
    auto &shape_functions = level_shape_functions[level];
    E::gaussPts.resize(level_ref_gauss_pts.size1(), level_ref_gauss_pts.size2(),
                       false);
    noalias(E::gaussPts) = level_ref_gauss_pts;

    const auto fe_ent = E::numeredEntFiniteElementPtr->getEnt();
    auto get_fe_coords = [&]() {
      const EntityHandle *conn;
      int num_nodes;
      CHK_MOAB_THROW(
          E::mField.get_moab().get_connectivity(fe_ent, conn, num_nodes, true),
          "error get connectivity");
      VectorDouble coords(num_nodes * 3);
      CHK_MOAB_THROW(
          E::mField.get_moab().get_coords(conn, num_nodes, &*coords.begin()),
          "error get coordinates");
      return coords;
    };

    auto coords = get_fe_coords();

    const int num_nodes = level_ref_gauss_pts.size2();
    mapGaussPts.resize(level_ref_gauss_pts.size2());

    FTensor::Index<'i', 3> i;
    FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_n(
        &*shape_functions.data().begin());
    FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_coords(
        &verts_array[0][ver_count], &verts_array[1][ver_count],
        &verts_array[2][ver_count]);
    for (int gg = 0; gg != num_nodes; ++gg, ++ver_count) {

      mapGaussPts[gg] = start_vert_handle + ver_count;

      auto set_float_precision = [](const double x) {
        if (std::abs(x) < std::numeric_limits<float>::epsilon())
          return 0.;
        else
          return x;
      };

      t_coords(i) = 0;
      auto t_ele_coords = getFTensor1FromArray<3, 3>(coords);
      for (int nn = 0; nn != CN::VerticesPerEntity(type); ++nn) {
        t_coords(i) += t_n * t_ele_coords(i);
        ++t_ele_coords;
        ++t_n;
      }

      for (auto ii : {0, 1, 2})
        t_coords(ii) = set_float_precision(t_coords(ii));

      ++t_coords;
    }

    Tag th;
    int def_in_the_loop = -1;
    CHKERR postProcMesh.tag_get_handle("NB_IN_THE_LOOP", 1, MB_TYPE_INTEGER, th,
                                       MB_TAG_CREAT | MB_TAG_SPARSE,
                                       &def_in_the_loop);

    postProcElements.clear();
    const int num_el = level_ref_ele.size1();
    const int num_nodes_on_ele = level_ref_ele.size2();
    auto start_e = start_ele_handle + ele_count;
    postProcElements = Range(start_e, start_e + num_el - 1);
    for (auto tt = 0; tt != level_ref_ele.size1(); ++tt, ++ele_count) {
      for (int nn = 0; nn != num_nodes_on_ele; ++nn) {
        conn[num_nodes_on_ele * ele_count + nn] =
            mapGaussPts[level_ref_ele(tt, nn)];
      }
    }

    const int n_in_the_loop = E::nInTheLoop;
    CHKERR postProcMesh.tag_clear_data(th, postProcElements, &n_in_the_loop);

    MoFEMFunctionReturn(0);
  };

  CHKERR set_gauss_pts(

      ref_ele->levelGaussPtsOnRefMesh, ref_ele->levelRef,
      ref_ele->levelShapeFunctions,

      ref_ele->startingVertEleHandle, ref_ele->startingEleHandle,
      ref_ele->verticesOnEleArrays, ref_ele->eleConn, ref_ele->countVertEle,
      ref_ele->countEle

  );

  MoFEMFunctionReturn(0);
};

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::preProcess() {
  moab::Interface &moab = coreMesh;
  MoFEMFunctionBegin;

  ParallelComm *pcomm_post_proc_mesh =
      ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm_post_proc_mesh != NULL)
    delete pcomm_post_proc_mesh;

  CHKERR postProcMesh.delete_mesh();

  auto get_ref_ele = [&](const EntityType type) {
    PostProcGenerateRefMeshPtr ref_ele_ptr;

    try {

      ref_ele_ptr = refElementsMap.at(type);

    } catch (const out_of_range &e) {

      switch (type) {
      case MBTET:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBTET>>();
        break;
      case MBHEX:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBHEX>>();
        break;
      case MBTRI:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBTRI>>();
        break;
      case MBQUAD:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBQUAD>>();
        break;
      case MBEDGE:
        ref_ele_ptr = boost::make_shared<PostProcGenerateRefMesh<MBEDGE>>();
        break;
      default:
        MOFEM_LOG("SELF", Sev::error)
            << "Generation of reference elements for type < "
            << moab::CN::EntityTypeName(type) << " > is not implemented";
        CHK_THROW_MESSAGE(MOFEM_NOT_IMPLEMENTED, "Element not implemented");
      }
    }

    CHK_THROW_MESSAGE(ref_ele_ptr->generateReferenceElementMesh(),
                      "Error when generating reference element");

    refElementsMap[type] = ref_ele_ptr;

    return ref_ele_ptr;
  };

  auto fe_ptr = this->problemPtr->numeredFiniteElementsPtr;

  auto miit =
      fe_ptr->template get<Composite_Name_And_Part_mi_tag>().lower_bound(
          boost::make_tuple(this->getFEName(), this->getLoFERank()));
  auto hi_miit =
      fe_ptr->template get<Composite_Name_And_Part_mi_tag>().upper_bound(
          boost::make_tuple(this->getFEName(), this->getHiFERank()));

  const int number_of_ents_in_the_loop = this->getLoopSize();
  if (std::distance(miit, hi_miit) != number_of_ents_in_the_loop) {
    SETERRQ(E::mField.get_comm(), MOFEM_DATA_INCONSISTENCY,
            "Wrong size of indicies. Inconsistent size number of iterated "
            "elements iterated by problem and from range.");
  }

  for (auto &m : refElementsMap) {
    m.second->nbVertices = 0;
    m.second->nbEles = 0;
    m.second->countEle = 0;
    m.second->countVertEle = 0;
  }

  for (; miit != hi_miit; ++miit) {
    auto type = (*miit)->getEntType();
    auto ref_ele = get_ref_ele(type);

    // Set pointer to element. So that getDataVectorDofsPtr in getMaxLevel
    // can work
    E::numeredEntFiniteElementPtr = *miit;
    bool add = true;
    if (E::exeTestHook) {
      add = E::exeTestHook(this);
    }

    if (add) {
      size_t level = getMaxLevel();
      level = std::min(level, ref_ele->levelGaussPtsOnRefMesh.size() - 1);
      ref_ele->nbVertices += ref_ele->levelGaussPtsOnRefMesh[level].size2();
      ref_ele->nbEles += ref_ele->levelRef[level].size1();
    }
  }

  auto alloc_vertices_and_elements_on_post_proc_mesh = [&]() {
    MoFEMFunctionBegin;

    ReadUtilIface *iface;
    CHKERR postProcMesh.query_interface(iface);

    for (auto &m : refElementsMap) {
      if (m.second->nbEles) {
        CHKERR iface->get_node_coords(3, m.second->nbVertices, 0,
                                      m.second->startingVertEleHandle,
                                      m.second->verticesOnEleArrays);
        CHKERR iface->get_element_connect(
            m.second->nbEles, m.second->levelRef[0].size2(), m.first, 0,
            m.second->startingEleHandle, m.second->eleConn);

        m.second->countEle = 0;
        m.second->countVertEle = 0;
      }
    }

    MoFEMFunctionReturn(0);
  };

  CHKERR alloc_vertices_and_elements_on_post_proc_mesh();

  MoFEMFunctionReturn(0);
}

template <typename E>
MoFEMErrorCode PostProcBrokenMeshInMoabBase<E>::postProcess() {
  MoFEMFunctionBeginHot;

  auto update_elements = [&]() {
    ReadUtilIface *iface;
    CHKERR this->postProcMesh.query_interface(iface);
    MoFEMFunctionBegin;

    for (auto &m : refElementsMap) {
      if (m.second->nbEles) {
        MOFEM_TAG_AND_LOG("SELF", Sev::noisy, "PostProc")
            << "Update < " << moab::CN::EntityTypeName(m.first)
            << m.second->countEle;
        CHKERR iface->update_adjacencies(
            m.second->startingEleHandle, m.second->countEle,
            m.second->levelRef[0].size2(), m.second->eleConn);
      }
    }

    MoFEMFunctionReturn(0);
  };

  auto resolve_shared_ents = [&]() {
    MoFEMFunctionBegin;
    auto get_lower_dimension = [&]() {
      int min_dim = 3;
      for (auto &m : refElementsMap) {
        const int dim = moab::CN::Dimension(m.first);
        if (m.second->nbEles) {
          min_dim = std::min(dim, min_dim);
        }
      }
      return min_dim;
    };

    auto remove_obsolete_entities = [&](auto &&min_dim) {
      Range edges;
      CHKERR postProcMesh.get_entities_by_type(0, MBEDGE, edges, false);

      Range faces;
      CHKERR postProcMesh.get_entities_by_dimension(0, 2, faces, false);

      Range vols;
      CHKERR postProcMesh.get_entities_by_dimension(0, 3, vols, false);

      Range ents;

      if (min_dim > 1)
        CHKERR postProcMesh.delete_entities(edges);
      else
        ents.merge(edges);

      if (min_dim > 2)
        CHKERR postProcMesh.delete_entities(faces);
      else
        ents.merge(faces);

      ents.merge(vols);

      return ents;
    };

    auto ents = remove_obsolete_entities(get_lower_dimension());

    ParallelComm *pcomm_post_proc_mesh =
        ParallelComm::get_pcomm(&(postProcMesh), MYPCOMM_INDEX);
    if (pcomm_post_proc_mesh == NULL) {
      // wrapRefMeshComm =
      // boost::make_shared<WrapMPIComm>(T::mField.get_comm(), false);
      pcomm_post_proc_mesh = new ParallelComm(
          &(postProcMesh),
          PETSC_COMM_WORLD /*(T::wrapRefMeshComm)->get_comm()*/);
    }

    int rank = E::mField.get_comm_rank();
    CHKERR postProcMesh.tag_clear_data(pcomm_post_proc_mesh->part_tag(), ents,
                                       &rank);

    CHKERR pcomm_post_proc_mesh->resolve_shared_ents(0);

    MoFEMFunctionReturn(0);
  };

  CHKERR update_elements();
  CHKERR resolve_shared_ents();

  MoFEMFunctionReturnHot(0);
}

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