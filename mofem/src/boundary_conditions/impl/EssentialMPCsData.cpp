/**
 * @file EssentialMPCsData.cpp
 * @brief MPC boundary conditions
 * @version 13.1
 * @date 2022-09-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <MoFEM.hpp>

namespace MoFEM {

MoFEMErrorCode setMPCParentAdjacency() {
  MoFEMFunctionBegin;

  MoFEMFunctionReturn(0);
}


EssentialPreProc<MPCsType>::EssentialPreProc(MoFEM::Interface &m_field,
                                             boost::shared_ptr<FEMethod> fe_ptr,
                                             bool is_spatial_positions)
    : mField(m_field), fePtr(fe_ptr), isSpatialPositions(is_spatial_positions) {
}

MoFEMErrorCode EssentialPreProc<MPCsType>::loadFileWithMPCs(
    Interface &m_field, const char *file_name, const char *options) {
  MoFEMFunctionBegin;

  auto simple = m_field.getInterface<Simple>();
  auto &moab = m_field.get_moab();

  // CHKERR moab.load_file(file_name, 0, options);
  CHKERR simple->loadFile("", file_name);

  constexpr bool is_debug = false;

  if (m_field.get_comm_size() == 1)
    MoFEMFunctionReturnHot(0);

  Range all_ents;
  CHKERR moab.get_entities_by_handle(0, all_ents, false);

  auto print_range_on_procs = [&](const Range &range, std::string name) {
    if(!is_debug) return;
    MOFEM_LOG("SYNC", Sev::inform)
        << name << " on proc [" << m_field.get_comm_rank() << "] : \n"
        << range;
    MOFEM_LOG_SYNCHRONISE(m_field.get_comm());
  };

  auto save_range_to_file = [&](const Range range, std::string name = "part") {
    if(!is_debug) return;
    int rr = m_field.get_comm_rank();
    ostringstream ss;
    ss << "out_" << name << "_" << rr << ".vtk";
    MOFEM_LOG("SELF", Sev::inform) << "Save debug part mesh " << ss.str();
    EntityHandle meshset;
    CHKERR moab.create_meshset(MESHSET_SET, meshset);
    CHKERR moab.add_entities(meshset, range);
    if (!range.empty())
      CHKERR moab.write_file(ss.str().c_str(), "VTK", "", &meshset, 1);
    CHKERR moab.delete_entities(&meshset, 1);
  };

  auto master_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
          std::regex((boost::format("%s(.*)") % "MPC_(.*)_LINKS").str()));
  Range mpc_ents, links_ents;
  for (auto m : master_meshset_ptr) {

    Range ents;
    CHKERR moab.get_entities_by_handle(m->getMeshset(), ents, true);

    links_ents.merge(ents);
    for (auto &link : links_ents) {
      Range verts;
      CHKERR moab.get_connectivity(&link, 1, verts, true);
      mpc_ents.merge(verts);
    }
  }

  // if(mpc_ents.empty()) {
  //   MoFEMFunctionReturnHot(0);
  // }
  
  save_range_to_file(all_ents, "all_ents");
  print_range_on_procs(mpc_ents, "mpc_ents");
  print_range_on_procs(links_ents, "links_ents");

  const int dim = simple->getDim();

  ParallelComm *pcomm =
      ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pcomm == NULL)
    pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

  Tag part_tag = pcomm->part_tag();
  Range tagged_sets, proc_ents, off_proc_ents, all_proc_ents;

  CHKERR moab.get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag, NULL, 1,
                                           tagged_sets, moab::Interface::UNION);
  print_range_on_procs(tagged_sets, "tagged_sets");
  for (auto &mit : tagged_sets) {
    int part;
    CHKERR moab.tag_get_data(part_tag, &mit, 1, &part);

    if (part == m_field.get_comm_rank()) {
      CHKERR moab.get_entities_by_dimension(mit, dim, proc_ents, true);
      CHKERR moab.get_entities_by_handle(mit, all_proc_ents, true);
    } else {
      CHKERR moab.get_entities_by_handle(mit, off_proc_ents, true);
    }
  }

  print_range_on_procs(proc_ents, "proc_ents");

  save_range_to_file(proc_ents);
  save_range_to_file(off_proc_ents, "all_proc_ents");
  all_proc_ents = subtract(all_proc_ents, links_ents);

  Range all_tets;
  CHKERR moab.get_entities_by_handle(0, all_tets, true);

  Range tets_skin;
  Skinner skin(&moab);
  CHKERR skin.find_skin(0, all_tets.subset_by_dimension(dim), false, tets_skin);


  // Find the skin of entities on this processor
  Range proc_ents_skin[4]; // stores only entities shared with other processors
  CHKERR skin.find_skin(0, proc_ents, false, proc_ents_skin[dim - 1]);


  Range skin_verts, verts_to_add, edges_to_add;
  CHKERR moab.get_connectivity(proc_ents_skin[dim - 1],
                               skin_verts); // FIXME: check for 3D
  for (auto &link : links_ents) {
    Range verts;
    CHKERR moab.get_connectivity(&link, 1, verts);
    for (auto &vert : verts)
      if (skin_verts.find(vert) != skin_verts.end()) {
        edges_to_add.insert(link);
        verts_to_add.merge(verts);
      }
  }

  // Remove shared entities from the skin
  proc_ents_skin[dim - 1] = subtract(proc_ents_skin[dim - 1], tets_skin);

  // Get adjacent entities and their connectivity

  // FIXME: these are NOT (!) only shared entities (includes outside)
  if (dim > 2) {
    CHKERR moab.get_adjacencies(proc_ents_skin[2], 1, false, proc_ents_skin[1],
                                moab::Interface::UNION);
  }
  Range adj_ents;
  CHKERR moab.get_connectivity(proc_ents_skin[1], adj_ents, true);
  proc_ents_skin[0].merge(adj_ents);

  proc_ents_skin[1].merge(edges_to_add);
  proc_ents_skin[0].merge(verts_to_add);

  print_range_on_procs(edges_to_add, "edges_to_add");
  print_range_on_procs(verts_to_add, "verts_to_add");

  // do it as early as possible
  // it slows down all get_connectivity and get_adjacencies

  // Remove entities not part of this processor
  save_range_to_file(proc_ents, "proc_ents");
  save_range_to_file(off_proc_ents, "off_proc_ents");

  all_tets = off_proc_ents;

  if (edges_to_add.empty())
    all_tets.merge(links_ents);

  save_range_to_file(proc_ents_skin[0], "proc_ents_skin0");
  save_range_to_file(proc_ents_skin[1], "proc_ents_skin1");

  for (int dd = dim; dd >= 0; --dd) {
    all_tets = subtract(all_tets, proc_ents_skin[dd]);
  }
  print_range_on_procs(proc_ents, "proc_ents");

  Range meshsets;
  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, true);
  for (auto m : meshsets) {
    CHKERR moab.remove_entities(m, all_tets);
  }

  print_range_on_procs(all_tets.subset_by_dimension(2), "all_tets to delete");
  save_range_to_file(all_tets.subset_by_dimension(2), "part_delete");

  // vertices do not need to be deleted
  for (int dd = dim; dd > 0; --dd) {
    CHKERR moab.delete_entities(all_tets.subset_by_dimension(dd));
  }

  // Update parallel status tag
  {
    Range all_ents;
    CHKERR moab.get_entities_by_handle(0, all_ents);
    std::vector<unsigned char> pstat_tag(all_ents.size(), 0);
    save_range_to_file(all_ents, "all_ents_part");
    CHKERR moab.tag_set_data(pcomm->pstatus_tag(), all_ents,
                             &*pstat_tag.begin());
  }

  Range all_ents_before;
  CHKERR moab.get_entities_by_handle(0, all_ents_before);
  save_range_to_file(all_ents_before, "all_ents_part_before");

  CHKERR pcomm->resolve_shared_ents(0, proc_ents, dim, -1, proc_ents_skin);

  Range all_ents_after;
  CHKERR moab.get_entities_by_handle(0, all_ents_after);
  save_range_to_file(all_ents_after, "all_ents_part_after");

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EssentialPreProc<MPCsType>::addLinksToPostProcMesh(
    Interface &m_field, moab::Interface &post_proc_mesh_interface,
    vector<std::string> field_names) {
  MoFEMFunctionBegin;

  // for each link vertex find entity in mfield database
  // and add it to post_proc_mesh_interface
  // for each field give in field_names save data on appropriate tag

  std::array<double, 9> def;
  std::fill(def.begin(), def.end(), 0);

  auto get_tag = [&](const std::string name, size_t size) {
    Tag th;
    size = size <= 1 ? 1 : (size <= 3 ? 3 : 9);
    CHKERR post_proc_mesh_interface.tag_get_handle(
        name.c_str(), size, MB_TYPE_DOUBLE, th, MB_TAG_CREAT | MB_TAG_SPARSE,
        def.data());
    return th;
  };

  // get links from meshsets
  auto master_meshset_ptr =
      m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
          std::regex((boost::format("%s(.*)") % "MPC_(.*)_LINKS").str()));
  Range links_ents;
  for (auto m : master_meshset_ptr) {
    Range ents;
    CHKERR m_field.get_moab().get_entities_by_type(m->getMeshset(), MBEDGE,
                                                           ents, true);
    links_ents.merge(ents);
  }

  if(links_ents.empty()) {
    MoFEMFunctionReturnHot(0);
  }

  // create edges in post_proc_mesh_interface
  std::vector<EntityHandle> p_links_edges(links_ents.size()); // post_proc edges
  std::vector<EntityHandle> o_nodes(links_ents.size() * 2); // original mesh nodes
  std::vector<EntityHandle> p_nodes(links_ents.size() * 2); // post_proc nodes

  int dd = 0;
  for (auto &link : links_ents) {
    Range verts;
    CHKERR m_field.get_moab().get_connectivity(&link, 1, verts, true);
    ublas::vector<EntityHandle> edge_conn(2);
    int gg = 0;
    for (auto &vert : verts) {
      VectorDouble coords(3);
      // fill o_nodes with vertices
      CHKERR m_field.get_moab().get_coords(&vert, 1, &coords[0]);
      CHKERR post_proc_mesh_interface.create_vertex(&coords[0], edge_conn[gg]);
      o_nodes[dd * 2 + gg] = vert;
      p_nodes[dd * 2 + gg] = edge_conn[gg];
      gg++;
    }

    CHKERR post_proc_mesh_interface.create_element(MBEDGE, &edge_conn[0], 2,
                                                   p_links_edges[dd++]);
  }

  const FieldEntity_multiIndex *field_ents;
  CHKERR m_field.get_field_ents(&field_ents);
  auto &field_ents_by_uid = field_ents->get<Unique_mi_tag>();

  for (auto field : field_names) {
    auto field_ptr = m_field.get_field_structure(field);
    const int nb_coefficients = field_ptr->getNbOfCoeffs();
    Tag tag = get_tag(field, nb_coefficients);

    int gg = 0;
    for (auto &node : o_nodes) {
      VectorDouble data(nb_coefficients);
      auto eit = field_ents_by_uid.find(FieldEntity::getLocalUniqueIdCalculate(
          m_field.get_field_bit_number(field), node));
      if (eit != field_ents_by_uid.end()) {
        noalias(data) = (*eit)->getEntFieldData();
      }
      auto ent = p_nodes[gg++];
      post_proc_mesh_interface.tag_set_data(tag, &ent, 1, &*data.begin());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EssentialPreProc<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (isSpatialPositions) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
  }

  MOFEM_LOG("BcMngWorld", Sev::warning)
      << "EssentialPreProc<MPCsType> has no effect here. ";

  MoFEMFunctionReturn(0);
}

EssentialPostProcLhs<MPCsType>::EssentialPostProcLhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Mat> lhs, SmartPetscObj<AO> ao)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vLhs(lhs), vAO(ao) {}

MoFEMErrorCode EssentialPostProcLhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    const auto problem_name = fe_ptr->problemPtr->getName();
    bool do_assembly = false;
    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          auto get_field_coeffs = [&](auto field_name) {
            auto field_ptr = mField.get_field_structure(field_name);
            return field_ptr->getNbOfCoeffs();
          };

          const auto nb_field_coeffs = get_field_coeffs(field_name);
          constexpr auto max_nb_dofs_per_node = 6;

          if (nb_field_coeffs > max_nb_dofs_per_node)
            MOFEM_LOG("WORLD", Sev::error)
                << "MultiPointConstraints PreProcLhs<MPCsType>: support only "
                   "up to 6 dofs per node.";
          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply MultiPointConstraints PreProc<MPCsType>: "
              << problem_name << "_" << field_name << "_" << block_name;

          auto mpc_type = mpc_bc->mpcType;
          switch (mpc_type) {
          case MPC::COUPLING:
          case MPC::TIE:
          case MPC::RIGID_BODY:
            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "MPC type not implemented");
          }

          // FIXME: there handle the vertices differently (block level)
          auto master_verts = mpc_bc->mpcMasterEnts.subset_by_type(MBVERTEX);
          auto slave_verts = mpc_bc->mpcSlaveEnts.subset_by_type(MBVERTEX);

          auto prb_name = fe_ptr->problemPtr->getName();
          auto get_flag = [&](int idx) { return (&mpc_bc->data.flag1)[idx]; };
          MOFEM_LOG("WORLD", Sev::verbose)
              << "Size of master nodes: " << mpc_bc->mpcMasterEnts.size()
              << " Size of slave nodes : " << mpc_bc->mpcSlaveEnts.size();

          if (mpc_bc->mpcMasterEnts.size() > mpc_bc->mpcSlaveEnts.size()) {
            MOFEM_LOG("WORLD", Sev::error)
                << "Size of master nodes < Size of slave nodes : "
                << mpc_bc->mpcMasterEnts.size() << " > " << mpc_bc->mpcSlaveEnts.size();
            // SETERRQ(PETSC_COMM_WORLD, MOFEM_OPERATION_UNSUCCESSFUL,
            //         "data inconsistency");
          }

          auto add_is = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          SmartPetscObj<IS> is_xyz_row_sum;
          SmartPetscObj<IS> is_m_row[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_m_col[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s_row[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s_col[max_nb_dofs_per_node];

          for (int dd = 0; dd != nb_field_coeffs; dd++) {
            if (get_flag(dd)) {
              // SmartPetscObj<IS> is_xyz_m;
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, ROW, field_name, dd, dd, is_m_row[dd],
                  &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, COL, field_name, dd, dd, is_m_col[dd],
                  &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, COL, field_name, dd, dd, is_s_col[dd],
                  &slave_verts);
              CHKERR is_mng->isCreateProblemFieldAndRank(
                  prb_name, ROW, field_name, dd, dd, is_s_row[dd],
                  &slave_verts);
              // ISView(is_s_row[dd], PETSC_VIEWER_STDOUT_WORLD);
              // ISView(is_s_col[dd], PETSC_VIEWER_STDOUT_WORLD);

              if (!mpc_bc->isReprocitical) {
                if (!is_xyz_row_sum) {
                  is_xyz_row_sum = is_s_row[dd];
                } else {
                  is_xyz_row_sum = add_is(is_xyz_row_sum, is_s_row[dd]);
                }
              }
            }
          }

          // if (is_xyz_row_sum) {
            SmartPetscObj<Mat> B =
                vLhs ? vLhs : SmartPetscObj<Mat>(fe_ptr->B, true);
            // The user is responsible for assembly if vLhs is provided
            MatType typem;
            CHKERR MatGetType(B, &typem);
            CHKERR MatSetOption(B, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
            if (typem == MATMPIBAIJ)
              CHKERR MatSetOption(B, MAT_USE_HASH_TABLE, PETSC_TRUE);

            if ((*fe_ptr->matAssembleSwitch) && !vLhs) {
              if (*fe_ptr->matAssembleSwitch) {
                CHKERR MatAssemblyBegin(B, MAT_FLUSH_ASSEMBLY);
                CHKERR MatAssemblyEnd(B, MAT_FLUSH_ASSEMBLY);
                *fe_ptr->matAssembleSwitch = false;
              }
            }

            if (vAO) {
              MOFEM_LOG("WORLD", Sev::error) << "No support for AO yet";
            }

            CHKERR MatZeroRowsIS(B, is_xyz_row_sum, 0, PETSC_NULL,
                                 PETSC_NULL);
          // }
            auto set_mat_values = [&](auto row_is, auto col_is, double val, double perturb = 0) {
              MoFEMFunctionBeginHot;
              const int *row_index_ptr;
              CHKERR ISGetIndices(row_is, &row_index_ptr);
              const int *col_index_ptr;
              CHKERR ISGetIndices(col_is, &col_index_ptr);
              int row_size, col_size;
              CHKERR ISGetLocalSize(row_is, &row_size);
              CHKERR ISGetLocalSize(col_is, &col_size);

              MatrixDouble locMat(row_size, col_size);
              fill(locMat.data().begin(), locMat.data().end(), 0.0);
              for (int i = 0; i != row_size; i++)
                for (int j = 0; j != col_size; j++)
                  if (i == j || col_size == 1)
                    locMat(i, j) = val;

              if (mpc_bc->isReprocitical) {
                locMat *= perturb;
                CHKERR ::MatSetValues(B, row_size, row_index_ptr, col_size,
                                      col_index_ptr, &*locMat.data().begin(),
                                      ADD_VALUES);
              } else {
                CHKERR ::MatSetValues(B, row_size, row_index_ptr, col_size,
                                      col_index_ptr, &*locMat.data().begin(),
                                      INSERT_VALUES);
              }

              CHKERR ISRestoreIndices(row_is, &row_index_ptr);
              CHKERR ISRestoreIndices(col_is, &col_index_ptr);

              MoFEMFunctionReturnHot(0);
            };

            for (int dd = 0; dd != nb_field_coeffs; dd++) {
              if (get_flag(dd)) {
                do_assembly = true;
                if (!mpc_bc->isReprocitical) {
                  CHKERR set_mat_values(is_s_row[dd], is_s_col[dd], vDiag);
                  CHKERR set_mat_values(is_s_row[dd], is_m_col[dd], -vDiag);
                } 
                else
                {
                  auto &pn = mpc_bc->mPenalty;
                  CHKERR set_mat_values(is_s_row[dd], is_s_col[dd], vDiag, pn);
                  CHKERR set_mat_values(is_s_row[dd], is_m_col[dd], -vDiag, pn);
                  CHKERR set_mat_values(is_m_row[dd], is_m_col[dd], vDiag, pn);
                  CHKERR set_mat_values(is_m_row[dd], is_s_col[dd], -vDiag, pn);
                }
              }
            }
        } // if regex
      } // mpc loop
    }   // bc loop
    if (do_assembly) {
      SmartPetscObj<Mat> B = vLhs ? vLhs : SmartPetscObj<Mat>(fe_ptr->B, true);
      CHKERR MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
      CHKERR MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    }

  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  } // if fe_ptr
  MoFEMFunctionReturn(0);
}

EssentialPostProcRhs<MPCsType>::EssentialPostProcRhs(
    MoFEM::Interface &m_field, boost::shared_ptr<FEMethod> fe_ptr, double diag,
    SmartPetscObj<Vec> rhs)
    : mField(m_field), fePtr(fe_ptr), vDiag(diag), vRhs(rhs) {}

MoFEMErrorCode EssentialPostProcRhs<MPCsType>::operator()() {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;

  if (auto fe_method_ptr = fePtr.lock()) {

    auto bc_mng = mField.getInterface<BcManager>();
    auto is_mng = mField.getInterface<ISManager>();
    auto vec_mng = mField.getInterface<VecManager>();
    const auto problem_name = fe_method_ptr->problemPtr->getName();
    // get connectivity from edge and assign master, slave ents
    for (auto bc : bc_mng->getBcMapByBlockName()) {
      if (auto mpc_bc = bc.second->mpcPtr) {

        auto &bc_id = bc.first;

        auto regex_str = (boost::format("%s_(.*)") % problem_name).str();
        if (std::regex_match(bc_id, std::regex(regex_str))) {

          auto [field_name, block_name] =
              BcManager::extractStringFromBlockId(bc_id, problem_name);

          auto get_field_coeffs = [&](auto field_name) {
            auto field_ptr = mField.get_field_structure(field_name);
            return field_ptr->getNbOfCoeffs();
          };

          const auto nb_field_coeffs = get_field_coeffs(field_name);
          constexpr auto max_nb_dofs_per_node = 6;

          if (nb_field_coeffs > max_nb_dofs_per_node)
            MOFEM_LOG("WORLD", Sev::error)
                << "MultiPointConstraints PreProcLhs<MPCsType>: support only "
                   "up to 6 dofs per node for now.";

          MOFEM_LOG("WORLD", Sev::noisy)
              << "Apply MultiPointConstraints PreProc<MPCsType>: "
              << problem_name << "_" << field_name << "_" << block_name;

          auto mpc_type = mpc_bc->mpcType;
          switch (mpc_type) {
          case MPC::COUPLING:
          case MPC::TIE:
          case MPC::RIGID_BODY:
            break;
          default:
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "MPC type not implemented");
          }

          auto master_verts = mpc_bc->mpcMasterEnts.subset_by_type(MBVERTEX);
          auto slave_verts = mpc_bc->mpcSlaveEnts.subset_by_type(MBVERTEX);

          auto prb_name = fe_method_ptr->problemPtr->getName();

          auto get_flag = [&](int idx) { return (&mpc_bc->data.flag1)[idx]; };

          auto add_is = [](auto is1, auto is2) {
            IS is;
            CHK_THROW_MESSAGE(ISExpand(is1, is2, &is), "is sum");
            return SmartPetscObj<IS>(is);
          };

          SmartPetscObj<IS> is_m[max_nb_dofs_per_node];
          SmartPetscObj<IS> is_s[max_nb_dofs_per_node];

          for (int dd = 0; dd != nb_field_coeffs; dd++) {
            if (get_flag(dd)) {
        
              CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                  prb_name, ROW, field_name, dd, dd, is_m[dd], &master_verts);
              CHKERR is_mng->isCreateProblemFieldAndRankLocal(
                  prb_name, ROW, field_name, dd, dd, is_s[dd], &slave_verts);
  
            }
          }

          {

            if (auto fe_ptr = fePtr.lock()) {
              auto snes_ctx = fe_ptr->snes_ctx;
              auto ts_ctx = fe_ptr->ts_ctx;
              const bool is_nonlinear = snes_ctx != FEMethod::CTX_SNESNONE ||
                                        ts_ctx != FEMethod::CTX_TSNONE;
              // is_nonlinear = is_nonlinear || mpc_bc->isReprocitical;
              const int contrb = mpc_bc->isReprocitical ? 1 : 0;
              SmartPetscObj<Vec> F =
                  vRhs ? vRhs : SmartPetscObj<Vec>(fe_ptr->f, true);

              if (fe_ptr->vecAssembleSwitch) {
                if ((*fe_ptr->vecAssembleSwitch) && !vRhs) {
                  CHKERR VecGhostUpdateBegin(F, ADD_VALUES, SCATTER_REVERSE);
                  CHKERR VecGhostUpdateEnd(F, ADD_VALUES, SCATTER_REVERSE);
                  CHKERR VecAssemblyBegin(F);
                  CHKERR VecAssemblyEnd(F);
                  *fe_ptr->vecAssembleSwitch = false;
                }
              }

              auto vec_set_values = [&](auto is_xyz_m, auto is_xyz_s, double val) {
                MoFEMFunctionBeginHot;
                const int *m_index_ptr;
                CHKERR ISGetIndices(is_xyz_m, &m_index_ptr);
                const int *s_index_ptr;
                CHKERR ISGetIndices(is_xyz_s, &s_index_ptr);
                int size_m;
                CHKERR ISGetLocalSize(is_xyz_m, &size_m);
                int size_s;
                CHKERR ISGetLocalSize(is_xyz_s, &size_s);

                double *f;
                CHKERR VecGetArray(F, &f);
                // if (size_m > size_s)
                //   MOFEM_LOG("WORLD", Sev::error)
                //       << "Size of master IS > Size of slave IS : " << size_m
                //       << " > " << size_s;
                // if (size_m == 0)
                //   MOFEM_LOG("WORLD", Sev::error)
                //       << "Size of master IS is " << size_m;
                
                  if (is_nonlinear) {
                    auto x = fe_ptr->x;
                    auto tmp_x = vectorDuplicate(F);

                    CHKERR vec_mng->setLocalGhostVector(problem_name, ROW,
                                                        tmp_x, INSERT_VALUES,
                                                        SCATTER_FORWARD);
                    const double *v;
                    const double *u;

                    CHKERR VecGetArrayRead(tmp_x, &u);
                    CHKERR VecGetArrayRead(x, &v);
                    if (size_m && size_s) {
                      for (auto i = 0; i != size_s; ++i) {
                          auto m_idx = size_m == 1 ? 0 : i;
                          f[s_index_ptr[i]] = val * (v[s_index_ptr[i]] - v[m_index_ptr[m_idx]]) + f[s_index_ptr[i]] * contrb ;
                      }
                    }
                    CHKERR VecRestoreArrayRead(x, &v);
                    CHKERR VecRestoreArrayRead(tmp_x, &u);
                  } else {
                    if (size_m && size_s) 
                    for (auto i = 0; i != size_s; ++i) {
                      f[s_index_ptr[i]] = 0;
                    }
                  }
                
                CHKERR VecRestoreArray(F, &f);
                CHKERR ISRestoreIndices(is_xyz_m, &m_index_ptr);
                CHKERR ISRestoreIndices(is_xyz_s, &s_index_ptr);

                MoFEMFunctionReturnHot(0);
              };

              for (int dd = 0; dd != nb_field_coeffs; dd++)
                if (get_flag(dd)) {

                  if (!mpc_bc->isReprocitical) {
                    CHKERR vec_set_values(is_m[dd], is_s[dd], vDiag);
                  } else {
                    auto &pn = mpc_bc->mPenalty;
                    CHKERR vec_set_values(is_m[dd], is_s[dd], pn);
                    CHKERR vec_set_values(is_s[dd], is_m[dd], pn);
                  }
                }
            };

            // User is responsible for assembly if vLhs is provided

            // ISView(is_xyz_m, PETSC_VIEWER_STDOUT_WORLD);
            // CHKERR MatZeroRowsColumnsIS(B, is_xyz_m, vDiag, PETSC_NULL,
            //                             PETSC_NULL);
          } 
        }
      }
    }
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Can not lock shared pointer");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM