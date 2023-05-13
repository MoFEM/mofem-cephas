/** \file BcManager.cpp
 * \brief Manages boundary conditions
 * \ingroup bc_manager
 */

namespace MoFEM {

MoFEMErrorCode
BcManager::query_interface(boost::typeindex::type_index type_index,
                           UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<BcManager *>(this);
  MoFEMFunctionReturnHot(0);
}

BcManager::BcManager(const Core &core) : cOre(const_cast<Core &>(core)) {

  if (!LogManager::checkIfChannelExist("BcMngWorld")) {
    auto core_log = logging::core::get();

    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmWorld(), "BcMngWorld"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSync(), "BcMngSync"));
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "BcMngSelf"));

    LogManager::setLog("BcMngWorld");
    LogManager::setLog("BcMngSync");
    LogManager::setLog("BcMngSelf");

    MOFEM_LOG_TAG("BcMngWorld", "BcMng");
    MOFEM_LOG_TAG("BcMngSync", "BcMng");
    MOFEM_LOG_TAG("BcMngSelf", "BcMng");
  }

  MOFEM_LOG("BcMngWorld", Sev::noisy) << "BC manager created";
}

MoFEMErrorCode BcManager::getOptions() {
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "BcManager options", "none");
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BcManager::removeBlockDOFsOnEntities(
    const std::string problem_name, const std::string block_name,
    const std::string field_name, int lo, int hi, bool get_low_dim_ents,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities(problem_name, block_name, field_name, lo, hi,
                                get_low_dim_ents);

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  for (auto m :
       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

           (boost::format("%s(.*)") % block_name).str()

               ))

  ) {
    const std::string bc_id =
        problem_name + "_" + field_name + "_" + m->getName();
    CHKERR remove_dofs_on_ents(bcMapByBlockName.at(bc_id)->bcEnts, lo, hi);
    bcMapByBlockName.at(bc_id)->bcMarkers = std::vector<unsigned char>();
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BcManager::pushMarkDOFsOnEntities(const std::string problem_name,
                                                 const std::string block_name,
                                                 const std::string field_name,
                                                 int lo, int hi,
                                                 bool get_low_dim_ents) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  auto get_dim = [&](const Range &ents) {
    for (auto d : {3, 2, 1})
      if (ents.num_of_dimension(d))
        return d;
    return 0;
  };

  auto get_adj_ents = [&](const Range &ents) {
    Range verts;
    CHKERR m_field.get_moab().get_connectivity(ents, verts, true);
    const auto dim = get_dim(ents);
    for (size_t d = 1; d < dim; ++d)
      CHKERR m_field.get_moab().get_adjacencies(ents, d, false, verts,
                                                moab::Interface::UNION);
    verts.merge(ents);
    CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(verts);
    return verts;
  };

  auto fix_disp = [&]() {
    MoFEMFunctionBegin;

    auto mark_fix_dofs = [&](std::vector<unsigned char> &marked_field_dofs,
                             const auto lo, const auto hi) {
      return prb_mng->modifyMarkDofs(problem_name, ROW, field_name, lo, hi,
                                     ProblemsManager::MarkOP::OR, 1,
                                     marked_field_dofs);
    };

    auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         bc->bcEnts, true);
        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
        CHKERR m->getAttributes(bc->bcAttributes);

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block " << m->getName() << " number of entities "
            << bc->bcEnts.size() << " number of attributes "
            << bc->bcAttributes.size() << " highest dim of entities "
            << get_dim(bc->bcEnts);
        CHKERR mark_fix_dofs(bc->bcMarkers, lo, hi);
        if (get_low_dim_ents) {
          auto low_dim_ents = get_adj_ents(bc->bcEnts);
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   low_dim_ents, bc->bcMarkers);
          bc->bcEnts.swap(low_dim_ents);
        } else
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEnts, bc->bcMarkers);

        const std::string bc_id =
            problem_name + "_" + field_name + "_" + m->getName();
        bcMapByBlockName[bc_id] = bc;
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

            (boost::format("%s(.*)") % block_name).str()

                ))

    );

    MoFEMFunctionReturn(0);
  };

  CHKERR fix_disp();

  MoFEMFunctionReturn(0);
}

boost::shared_ptr<BcManager::BCs>
BcManager::popMarkDOFsOnEntities(const std::string block_name) {
  auto bc_it = bcMapByBlockName.find(block_name);
  if (bc_it != bcMapByBlockName.end()) {
    auto bc = bc_it->second;
    bcMapByBlockName.erase(bc_it);
    return bc;
  }
  return boost::shared_ptr<BCs>();
}

BcManager::BcMarkerPtr
BcManager::getMergedBlocksMarker(std::vector<std::regex> bc_regex_vec) {
  BcManager::BcMarkerPtr boundary_marker_ptr;
  if (bcMapByBlockName.size()) {
    for (auto b : bcMapByBlockName) {
      for (auto &reg_name : bc_regex_vec) {
        if (std::regex_match(b.first, reg_name)) {
          if (!boundary_marker_ptr)
            boundary_marker_ptr =
                boost::make_shared<std::vector<char unsigned>>();
          boundary_marker_ptr->resize(b.second->bcMarkers.size(), 0);
          for (int i = 0; i != b.second->bcMarkers.size(); ++i) {
            (*boundary_marker_ptr)[i] |= b.second->bcMarkers[i];
          }
        }
      }
    }
  }
  return boundary_marker_ptr;
}

BcManager::BcMarkerPtr BcManager::getMergedBlocksMarker(
    const std::vector<BcManager::BcMarkerPtr> &boundary_markers_ptr_vec) {
  auto boundary_marker_ptr = boost::make_shared<std::vector<char unsigned>>();
  for (auto &bcm : boundary_markers_ptr_vec) {
    boundary_marker_ptr->resize(bcm->size(), 0);
    for (int i = 0; i != bcm->size(); ++i)
      (*boundary_marker_ptr)[i] |= (*bcm)[i];
  }
  return boundary_marker_ptr;
}

SmartPetscObj<IS> BcManager::getBlockIS(const std::string block_prefix,
                                        const std::string block_name,
                                        const std::string field_name,
                                        const std::string problem_name, int lo,
                                        int hi, SmartPetscObj<IS> is_expand) {
  Interface &m_field = cOre;

  const std::string bc_id =
      block_prefix + "_" + field_name + "_" + block_name + "(.*)";

  Range bc_ents;
  for (auto bc : getBcMapByBlockName()) {
    if (std::regex_match(bc.first, std::regex(bc_id))) {
      bc_ents.merge(*(bc.second->getBcEntsPtr()));
      MOFEM_LOG("BcMngWorld", Sev::verbose)
          << "Get entities from block and add to IS. Block name " << bc.first;
    }
  }

  SmartPetscObj<IS> is_bc;
  auto get_is = [&]() {
    MoFEMFunctionBegin;
    CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(bc_ents);
    CHKERR m_field.getInterface<ISManager>()->isCreateProblemFieldAndRank(
        problem_name, ROW, field_name, lo, hi, is_bc, &bc_ents);
    if (is_expand) {
      IS is_tmp;
      CHKERR ISExpand(is_bc, is_expand, &is_tmp);
      is_bc = SmartPetscObj<IS>(is_tmp);
    }
    CHKERR ISSort(is_bc);
    MoFEMFunctionReturn(0);
  };

  if (get_is())
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "IS is not created");

  return is_bc;
}

SmartPetscObj<IS> BcManager::getBlockIS(const std::string problem_name,
                                        const std::string block_name,
                                        const std::string field_name, int lo,
                                        int hi, SmartPetscObj<IS> is_expand) {
  return getBlockIS(problem_name, block_name, field_name, problem_name, lo, hi,
                    is_expand);
}

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  std::array<Range, 3> ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
           NODESET | DISPLACEMENTSET)) {

    const std::string bc_id =
        problem_name + "_" + field_name + "_DISPLACEMENTSET" +
        boost::lexical_cast<std::string>(m->getMeshsetId());

    auto bc = bcMapByBlockName.at(bc_id);

    if (bc->dispBcPtr) {
      if (bc->dispBcPtr->data.flag1) {
        ents_to_remove[0].merge(bc->bcEnts);
      }
      if (bc->dispBcPtr->data.flag2) {
        ents_to_remove[1].merge(bc->bcEnts);
      }
      if (bc->dispBcPtr->data.flag3) {
        ents_to_remove[2].merge(bc->bcEnts);
      }
      if (bc->dispBcPtr->data.flag4) {
        ents_to_remove[1].merge(bc->bcEnts);
        ents_to_remove[2].merge(bc->bcEnts);
      }
      if (bc->dispBcPtr->data.flag5) {
        ents_to_remove[0].merge(bc->bcEnts);
        ents_to_remove[2].merge(bc->bcEnts);
      }
      if (bc->dispBcPtr->data.flag6) {
        ents_to_remove[0].merge(bc->bcEnts);
        ents_to_remove[1].merge(bc->bcEnts);
      }
    }
    bc->bcMarkers = std::vector<unsigned char>();
  }

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  CHKERR remove_dofs_on_ents(ents_to_remove[0], 0, 0);
  CHKERR remove_dofs_on_ents(ents_to_remove[1], 1, 1);
  CHKERR remove_dofs_on_ents(ents_to_remove[2], 2, 2);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  Range ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
           NODESET | TEMPERATURESET)

  ) {
    const std::string bc_id =
        problem_name + "_" + field_name + "_TEMPERATURESET" +
        boost::lexical_cast<std::string>(m->getMeshsetId());
    auto bc = bcMapByBlockName.at(bc_id);
    ents_to_remove.merge(bc->bcEnts);
    bc->bcMarkers = std::vector<unsigned char>();
  }

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  CHKERR remove_dofs_on_ents(ents_to_remove, 0, MAX_DOFS_ON_ENTITY);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  Range ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(NODESET |
                                                                   HEATFLUXSET)

  ) {
    const std::string bc_id =
        problem_name + "_" + field_name + "_HEATFLUXSET" +
        boost::lexical_cast<std::string>(m->getMeshsetId());
    auto bc = bcMapByBlockName.at(bc_id);
    ents_to_remove.merge(bc->bcEnts);
    bc->bcMarkers = std::vector<unsigned char>();
  }

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  CHKERR remove_dofs_on_ents(ents_to_remove, 0, MAX_DOFS_ON_ENTITY);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcVectorMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcVectorMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  std::array<Range, 3> ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(BLOCKSET)) {

    const auto block_name = m->getName();

    std::string bc_id = problem_name + "_" + field_name + "_" + block_name;
    std::string regex_str;
    if (block_name_field_prefix) {
      regex_str = (boost::format("%s_%s_%s_((FIX_(ALL|X|Y|Z))|("
                                 "DISPLACEMENT|ROTATE))(.*)") %
                   problem_name % field_name % field_name)
                      .str();
    } else {
      regex_str = (boost::format("%s_%s_((FIX_(ALL|X|Y|Z))|("
                                 "DISPLACEMENT|ROTATE))(.*)") %
                   problem_name % field_name)
                      .str();
    }

    if (std::regex_match(bc_id, std::regex(regex_str))) {

      auto bc = bcMapByBlockName.at(bc_id);

      if (auto disp_bc = bc->dispBcPtr) {
        if (disp_bc->data.flag1) {
          ents_to_remove[0].merge(bc->bcEnts);
        }
        if (disp_bc->data.flag2) {
          ents_to_remove[1].merge(bc->bcEnts);
        }
        if (disp_bc->data.flag3) {
          ents_to_remove[2].merge(bc->bcEnts);
        }
        if (disp_bc->data.flag4) {
          ents_to_remove[1].merge(bc->bcEnts);
          ents_to_remove[2].merge(bc->bcEnts);
        }
        if (disp_bc->data.flag5) {
          ents_to_remove[0].merge(bc->bcEnts);
          ents_to_remove[2].merge(bc->bcEnts);
        }
        if (disp_bc->data.flag6) {
          ents_to_remove[0].merge(bc->bcEnts);
          ents_to_remove[1].merge(bc->bcEnts);
        }
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "BC type not implemented");
      }
    }
  }

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  CHKERR remove_dofs_on_ents(ents_to_remove[0], 0, 0);
  CHKERR remove_dofs_on_ents(ents_to_remove[1], 1, 1);
  CHKERR remove_dofs_on_ents(ents_to_remove[2], 2, 2);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
BcManager::removeBlockDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string block_name,
    const std::string field_name, bool get_low_dim_ents,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
      problem_name, block_name, field_name, get_low_dim_ents);

  Range ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(BLOCKSET)) {

    std::string bc_id = problem_name + "_" + field_name + "_" + m->getName();

    auto str = boost::format("%s_%s_%s(.*)")

               % problem_name % field_name % block_name;

    if (std::regex_match(bc_id, std::regex(str.str()))) {

      auto bc = bcMapByBlockName.at(bc_id);

      if (auto disp_bc = bc->tempBcPtr) {
        if (disp_bc->data.flag1) {
          ents_to_remove.merge(bc->bcEnts);
        }
      } else {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "BC type not implemented");
      }
    }
  }

  auto remove_dofs_on_ents = [&](const Range &ents, const int lo,
                                 const int hi) {
    if (is_distributed_mesh)
      return prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                           hi);
    else
      return prb_mng->removeDofsOnEntitiesNotDistributed(
          problem_name, field_name, ents, lo, hi);
  };

  CHKERR remove_dofs_on_ents(ents_to_remove, 0, MAX_DOFS_ON_ENTITY);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode
BcManager::pushMarkDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";

  auto get_dim = [&](const Range &ents) {
    for (auto d : {3, 2, 1})
      if (ents.num_of_dimension(d))
        return d;
    return 0;
  };

  auto get_adj_ents = [&](const Range &ents) {
    Range verts;
    CHKERR m_field.get_moab().get_connectivity(ents, verts, true);
    const auto dim = get_dim(ents);
    for (size_t d = 1; d < dim; ++d)
      CHKERR m_field.get_moab().get_adjacencies(ents, d, false, verts,
                                                moab::Interface::UNION);
    verts.merge(ents);
    CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(verts);
    return verts;
  };

  auto fix_disp = [&]() {
    MoFEMFunctionBegin;

    auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         bc->bcEnts, true);

        bc->dispBcPtr = boost::make_shared<DisplacementCubitBcData>();
        CHKERR m->getBcDataStructure(*(bc->dispBcPtr));

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block DISPLACEMENTSET number of entities "
            << bc->bcEnts.size() << " highest dim of entities "
            << get_dim(bc->bcEnts);
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->dispBcPtr;

        if (bc->dispBcPtr->data.flag1)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        if (bc->dispBcPtr->data.flag2)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        if (bc->dispBcPtr->data.flag3)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        if (bc->dispBcPtr->data.flag4) {

          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        }
        if (bc->dispBcPtr->data.flag5) {

          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        }
        if (bc->dispBcPtr->data.flag6) {

          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        }

        if (get_low_dim_ents) {
          auto low_dim_ents = get_adj_ents(bc->bcEnts);
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   low_dim_ents, bc->bcMarkers);
          bc->bcEnts.swap(low_dim_ents);
        } else
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEnts, bc->bcMarkers);

        const std::string bc_id =
            problem_name + "_" + field_name + "_DISPLACEMENTSET" +
            boost::lexical_cast<std::string>(m->getMeshsetId());
        bcMapByBlockName[bc_id] = bc;
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
            NODESET | DISPLACEMENTSET)

    );

    MoFEMFunctionReturn(0);
  };

  CHKERR fix_disp();

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";

  auto get_dim = [&](const Range &ents) {
    for (auto d : {3, 2, 1})
      if (ents.num_of_dimension(d))
        return d;
    return 0;
  };

  auto get_adj_ents = [&](const Range &ents) {
    Range verts;
    CHKERR m_field.get_moab().get_connectivity(ents, verts, true);
    const auto dim = get_dim(ents);
    for (size_t d = 1; d < dim; ++d)
      CHKERR m_field.get_moab().get_adjacencies(ents, d, false, verts,
                                                moab::Interface::UNION);
    verts.merge(ents);
    CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(verts);
    return verts;
  };

  auto fix_temp = [&]() {
    MoFEMFunctionBegin;

    auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         bc->bcEnts, true);
        bc->tempBcPtr = boost::make_shared<TemperatureCubitBcData>();
        CHKERR m->getBcDataStructure(*(bc->tempBcPtr));

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block TEMPERATURESET number of entities "
            << bc->bcEnts.size() << " highest dim of entities "
            << get_dim(bc->bcEnts);
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->tempBcPtr;

        CHKERR prb_mng->modifyMarkDofs(
            problem_name, ROW, field_name, 0, MAX_DOFS_ON_ENTITY,
            ProblemsManager::MarkOP::OR, 1, bc->bcMarkers);

        if (get_low_dim_ents) {
          auto low_dim_ents = get_adj_ents(bc->bcEnts);
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   low_dim_ents, bc->bcMarkers);
          bc->bcEnts.swap(low_dim_ents);
        } else
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEnts, bc->bcMarkers);

        const std::string bc_id =
            problem_name + "_" + field_name + "_TEMPERATURESET" +
            boost::lexical_cast<std::string>(m->getMeshsetId());
        bcMapByBlockName[bc_id] = bc;
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
            NODESET | TEMPERATURESET)

    );

    MoFEMFunctionReturn(0);
  };

  CHKERR fix_temp();

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";

  auto get_dim = [&](const Range &ents) {
    for (auto d : {3, 2, 1})
      if (ents.num_of_dimension(d))
        return d;
    return 0;
  };

  auto get_adj_ents = [&](const Range &ents) {
    Range verts;
    CHKERR m_field.get_moab().get_connectivity(ents, verts, true);
    const auto dim = get_dim(ents);
    for (size_t d = 1; d < dim; ++d)
      CHKERR m_field.get_moab().get_adjacencies(ents, d, false, verts,
                                                moab::Interface::UNION);
    verts.merge(ents);
    CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(verts);
    return verts;
  };

  auto fix_disp = [&]() {
    MoFEMFunctionBegin;

    auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         bc->bcEnts, true);
        bc->heatFluxBcPtr = boost::make_shared<HeatFluxCubitBcData>();
        CHKERR m->getBcDataStructure(*(bc->heatFluxBcPtr));

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block HEATFLUX number of entities " << bc->bcEnts.size()
            << " highest dim of entities " << get_dim(bc->bcEnts);
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->heatFluxBcPtr;

        CHKERR prb_mng->modifyMarkDofs(
            problem_name, ROW, field_name, 0, MAX_DOFS_ON_ENTITY,
            ProblemsManager::MarkOP::OR, 1, bc->bcMarkers);

        if (get_low_dim_ents) {
          auto low_dim_ents = get_adj_ents(bc->bcEnts);
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   low_dim_ents, bc->bcMarkers);
          bc->bcEnts.swap(low_dim_ents);
        } else
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEnts, bc->bcMarkers);

        const std::string bc_id =
            problem_name + "_" + field_name + "_HEATFLUXSET" +
            boost::lexical_cast<std::string>(m->getMeshsetId());
        bcMapByBlockName[bc_id] = bc;
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(SIDESET |
                                                                    HEATFLUXSET)

    );

    MoFEMFunctionReturn(0);
  };

  CHKERR fix_disp();

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcVectorMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  MoFEMFunctionBegin;

  auto mark_dofs = [&](const string block_name, const int &idx_0,
                            const int &idx_1) {
    MoFEMFunctionBeginHot;
    if (block_name_field_prefix) {
      const string field_block = field_name + "_" + block_name;
      CHKERR pushMarkDOFsOnEntities(problem_name, field_block, field_name,
                                    idx_0, idx_1, get_low_dim_ents);
    } else {

      CHKERR pushMarkDOFsOnEntities(problem_name, block_name, field_name, idx_0,
                                    idx_1, get_low_dim_ents);
    }
    MoFEMFunctionReturnHot(0);
  };

  // displacement
  CHKERR mark_dofs("FIX_X", 0, 0);
  CHKERR mark_dofs("FIX_Y", 1, 1);
  CHKERR mark_dofs("FIX_Z", 2, 2);
  CHKERR mark_dofs("FIX_ALL", 0, MAX_DOFS_ON_ENTITY);

  // rotation
  CHKERR mark_dofs("ROTATE_X", 1, 1);
  CHKERR mark_dofs("ROTATE_X", 2, 2);
  CHKERR mark_dofs("ROTATE_Y", 0, 0);
  CHKERR mark_dofs("ROTATE_Y", 2, 2);
  CHKERR mark_dofs("ROTATE_Z", 0, 0);
  CHKERR mark_dofs("ROTATE_Z", 1, 1);
  CHKERR mark_dofs("ROTATE_ALL", 0, MAX_DOFS_ON_ENTITY);

  std::string regex_str;
  if (block_name_field_prefix) {
    regex_str = (boost::format("%s_%s_%s_(.*)") % problem_name % field_name %
                 field_name)
                    .str();
  } else {
    regex_str = (boost::format("%s_%s_(.*)") % problem_name % field_name).str();
  }

  for (auto &m : bcMapByBlockName) {
    auto &bc_id = m.first;
    if (std::regex_match(bc_id, std::regex(regex_str))) {
      auto &bc = m.second;
      if (std::regex_match(bc_id, std::regex("(.*)_FIX_X(.*)"))) {
        bc->dispBcPtr = boost::make_shared<DisplacementCubitBcData>();
        bc->dispBcPtr->data.flag1 = 1;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value1 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block but have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() >= 1) {
          bc->dispBcPtr->data.value1 = bc->bcAttributes[0];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add X " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *bc->dispBcPtr;
      } else if (std::regex_match(bc_id, std::regex("(.*)_FIX_Y(.*)"))) {
        bc->dispBcPtr = boost::make_shared<DisplacementCubitBcData>();
        bc->dispBcPtr->data.flag2 = 1;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value2 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block but have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() == 1) {
          bc->dispBcPtr->data.value2 = bc->bcAttributes[0];
        } else if (bc->bcAttributes.size() >= 2) {
          bc->dispBcPtr->data.value2 = bc->bcAttributes[1];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add Y " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
      } else if (std::regex_match(bc_id, std::regex("(.*)_FIX_Z(.*)"))) {
        bc->dispBcPtr = boost::make_shared<DisplacementCubitBcData>();
        bc->dispBcPtr->data.flag3 = 1;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value3 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block but have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() == 1) {
          bc->dispBcPtr->data.value3 = bc->bcAttributes[0];
        } else if (bc->bcAttributes.size() == 3) {
          bc->dispBcPtr->data.value3 = bc->bcAttributes[2];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add Z " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
      } else if (std::regex_match(bc_id, std::regex("(.*)_FIX_ALL(.*)"))) {
        bc->dispBcPtr = boost::make_shared<DisplacementCubitBcData>();
        bc->dispBcPtr->data.flag1 = 1;
        bc->dispBcPtr->data.flag2 = 1;
        bc->dispBcPtr->data.flag3 = 1;
        if (bc->bcAttributes.size() >= 1) {
          bc->dispBcPtr->data.value1 = bc->bcAttributes[0];
        }
        if (bc->bcAttributes.size() >= 2) {
          bc->dispBcPtr->data.value2 = bc->bcAttributes[1];
        }
        if (bc->bcAttributes.size() >= 3) {
          bc->dispBcPtr->data.value3 = bc->bcAttributes[2];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add ALL " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
      } else if (std::regex_match(bc_id, std::regex("(.*)_ROTATE_X(.*)"))) {
        bc->dispBcPtr =
            boost::make_shared<DisplacementCubitBcDataWithRotation>();
        bc->dispBcPtr->data.flag4 = 1;
        bc->dispBcPtr->data.flag5 = 0;
        bc->dispBcPtr->data.flag6 = 0;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value4 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected at least one attribute on block (angle, center "
                 "coords) but have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() >= 1) {
          bc->dispBcPtr->data.value4 = bc->bcAttributes[0];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add X " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *bc->dispBcPtr;
        if (bc->bcAttributes.size() == 4 || bc->bcAttributes.size() == 6) {
          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation *>(
                      bc->dispBcPtr.get())) {
            auto &o = ext_disp_bc->rotOffset;
            for (int a = 0; a != 3; ++a)
              o[a] = bc->bcAttributes[bc->bcAttributes.size() - 3 + a];
            MOFEM_LOG("BcMngWorld", Sev::inform)
                << "Add Rotate X Center: " << o[0] << " " << o[1] << " "
                << o[2];
          }
        }
      } else if (std::regex_match(bc_id, std::regex("(.*)_ROTATE_Y(.*)"))) {
        bc->dispBcPtr =
            boost::make_shared<DisplacementCubitBcDataWithRotation>();
        bc->dispBcPtr->data.flag4 = 0;
        bc->dispBcPtr->data.flag5 = 1;
        bc->dispBcPtr->data.flag6 = 0;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value5 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block (angle, center coords) but "
                 "have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() == 1 ||
                   bc->bcAttributes.size() == 4) {
          bc->dispBcPtr->data.value5 = bc->bcAttributes[0];
        } else if (bc->bcAttributes.size() == 6) {
          bc->dispBcPtr->data.value5 = bc->bcAttributes[1];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add Y " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
        if (bc->bcAttributes.size() == 4 || bc->bcAttributes.size() == 6) {
          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation *>(
                      bc->dispBcPtr.get())) {
            auto &o = ext_disp_bc->rotOffset;
            for (int a = 0; a != 3; ++a)
              o[a] = bc->bcAttributes[bc->bcAttributes.size() - 3 + a];
            MOFEM_LOG("BcMngWorld", Sev::inform)
                << "Add Rotate Y Center: " << o[0] << " " << o[1] << " "
                << o[2];
          }
        }
      } else if (std::regex_match(bc_id, std::regex("(.*)_ROTATE_Z(.*)"))) {
        bc->dispBcPtr =
            boost::make_shared<DisplacementCubitBcDataWithRotation>();
        bc->dispBcPtr->data.flag4 = 0;
        bc->dispBcPtr->data.flag5 = 0;
        bc->dispBcPtr->data.flag6 = 1;
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value6 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block (angle, center coords) but "
                 "have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() == 1 ||
                   bc->bcAttributes.size() > 4) {
          bc->dispBcPtr->data.value6 = bc->bcAttributes[0];
        } else if (bc->bcAttributes.size() == 3) {
          bc->dispBcPtr->data.value6 = bc->bcAttributes[2];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add Z " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
        if (bc->bcAttributes.size() == 4 || bc->bcAttributes.size() == 6) {
          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation *>(
                      bc->dispBcPtr.get())) {
            auto &o = ext_disp_bc->rotOffset;
            for (int a = 0; a != 3; ++a)
              o[a] = bc->bcAttributes[bc->bcAttributes.size() - 3 + a];
            MOFEM_LOG("BcMngWorld", Sev::inform)
                << "Add Rotate Z Center: " << o[0] << " " << o[1] << " "
                << o[2];
          }
        }
      } else if (std::regex_match(bc_id, std::regex("(.*)_ROTATE_ALL(.*)"))) {
        bc->dispBcPtr =
            boost::make_shared<DisplacementCubitBcDataWithRotation>();
        bc->dispBcPtr->data.flag4 = 1;
        bc->dispBcPtr->data.flag5 = 1;
        bc->dispBcPtr->data.flag6 = 1;
        if (bc->bcAttributes.size() >= 1) {
          bc->dispBcPtr->data.value4 = bc->bcAttributes[0];
        }
        if (bc->bcAttributes.size() >= 2) {
          bc->dispBcPtr->data.value5 = bc->bcAttributes[1];
        }
        if (bc->bcAttributes.size() >= 3) {
          bc->dispBcPtr->data.value6 = bc->bcAttributes[2];
        }
        MOFEM_LOG("BcMngWorld", Sev::inform) << "Add ALL " << bc_id;
        MOFEM_LOG("BcMngWorld", Sev::inform) << *(bc->dispBcPtr);
        if (bc->bcAttributes.size() > 3) {
          if (auto ext_disp_bc =
                  dynamic_cast<DisplacementCubitBcDataWithRotation *>(
                      bc->dispBcPtr.get())) {
            auto &o = ext_disp_bc->rotOffset;
            for (int a = 0; a != 3; ++a)
              o[a] = bc->bcAttributes[3 + a];
            MOFEM_LOG("BcMngWorld", Sev::inform)
                << "Add Rotate ALL Center: " << o[0] << " " << o[1] << " "
                << o[2];
          }
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string block_name,
    const std::string field_name, bool get_low_dim_ents) {
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities(problem_name, block_name, field_name, 0,
                                MAX_DOFS_ON_ENTITY, get_low_dim_ents);

  auto regex_str =
      (boost::format("%s_%s_%s(.*)") % problem_name % field_name % block_name)
          .str();

  for (auto &m : bcMapByBlockName) {

    auto &bc_id = m.first;

    if (std::regex_match(bc_id, std::regex(regex_str))) {

      auto &bc = m.second;
      bc->tempBcPtr = boost::make_shared<TemperatureCubitBcData>();
      bc->tempBcPtr->data.flag1 = 1;
      if (bc->bcAttributes.empty()) {
        bc->tempBcPtr->data.value1 = 0;
        MOFEM_LOG("BcMngWorld", Sev::warning)
            << "Expected one attribute on block but have "
            << bc->bcAttributes.size();
      } else if (bc->bcAttributes.size() >= 1) {
        bc->tempBcPtr->data.value1 = bc->bcAttributes[0];
      }
    }
  }

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<DisplacementCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  MoFEMFunctionBegin;
  // that marks DOFs and create data when are set by cubit nodesets. 
  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);
  // that marks DOFs and create data when are set by blocsket.
  CHKERR pushMarkDOFsOnEntities<BcVectorMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<DisplacementCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  MoFEMFunctionBegin;
  // that remove DOFs when are set by cubit nodesets. 
  CHKERR removeBlockDOFsOnEntities<BcMeshsetType<DISPLACEMENTSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);
  // that remove DOFs when are by blocksets  
  CHKERR removeBlockDOFsOnEntities<BcVectorMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);
  // add more ways to remove bcs when appropiate
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<TemperatureCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  MoFEMFunctionBegin;
  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  auto get_block_name = [&]() {
    if (block_name_field_prefix)
      return (boost::format("%s_FIX_SCALAR") % field_name).str();
    else
      return field_name;
  };

  CHKERR pushMarkDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
      problem_name, get_block_name(), field_name, get_low_dim_ents);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<TemperatureCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  MoFEMFunctionBegin;
  CHKERR removeBlockDOFsOnEntities<BcMeshsetType<TEMPERATURESET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);

  auto get_block_name = [&]() {
    if (block_name_field_prefix)
      return (boost::format("%s_FIX_SCALAR") % field_name).str();
    else
      return field_name;
  };

  CHKERR removeBlockDOFsOnEntities<BcScalarMeshsetType<BLOCKSET>>(
      problem_name, get_block_name(), field_name, get_low_dim_ents,
      is_distributed_mesh);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<HeatFluxCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  MoFEMFunctionBegin;
  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);
  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<HeatFluxCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  MoFEMFunctionBegin;
  CHKERR removeBlockDOFsOnEntities<BcMeshsetType<HEATFLUXSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);
  MoFEMFunctionReturn(0);
}

std::pair<std::string, std::string>
BcManager::extractStringFromBlockId(const std::string block_id,
                                    const std::string prb_name) {

  // Assumes that field name is consist with letters and numbers.
  // No special characters.
  auto field_rgx_str =
      (boost::format("%s_([a-zA-Z0-9]*)_(.*)") % prb_name).str();
  std::regex field_rgx(field_rgx_str);
  std::smatch match_field_name;
  std::string field_name;
  std::string block_name;

  if (std::regex_search(block_id, match_field_name, field_rgx)) {
    field_name = match_field_name[1];
    block_name = match_field_name[2];
  } else {
    CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY,
                      "Field name and block name can not be resolved");
  }

  return std::make_pair(field_name, block_name);
}

} // namespace MoFEM
