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
  PetscBool flg = PETSC_TRUE;
  MoFEMFunctionBeginHot;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "BcManager options", "none");
  ierr = PetscOptionsEnd();
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BcManager::removeBlockDOFsOnEntities(
    const std::string problem_name, const std::string block_name,
    const std::string field_name, int lo, int hi, bool get_lod_dim_ents) {
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

  auto remove_dofs_on_ents = [&](const Range &&ents, const int lo,
                                 const int hi) {
    MoFEMFunctionBegin;
    CHKERR prb_mng->removeDofsOnEntities(problem_name, field_name, ents, lo,
                                         hi);
    MoFEMFunctionReturn(0);
  };

  auto get_block_ents = [&](const std::string block_name) {
    Range remove_ents;
    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
      if (it->getName().compare(0, block_name.length(), block_name) == 0) {
        CHKERR m_field.get_moab().get_entities_by_handle(it->meshset,
                                                         remove_ents, true);
        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block to remove " << block_name << " number of entities "
            << remove_ents.size() << " highest dim of entities "
            << get_dim(remove_ents);
      }
    }
    return remove_ents;
  };

  if (get_lod_dim_ents)
    CHKERR remove_dofs_on_ents(get_adj_ents(get_block_ents(block_name)), lo,
                               hi);
  else
    CHKERR remove_dofs_on_ents(get_block_ents(block_name), lo, hi);

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

    for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, BLOCKSET, it)) {
      if (it->getName().compare(0, block_name.length(), block_name) == 0) {

        const std::string bc_id =
            problem_name + "_" + field_name + "_" + it->getName();

        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(it->meshset,
                                                         bc->bcEnts, true);
        CHKERR it->getAttributes(bc->bcAttributes);

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block " << block_name << " number of entities "
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

        bcMapByBlockName[bc_id] = bc;
      }
    }
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
    boundary_marker_ptr = boost::make_shared<std::vector<char unsigned>>();
    for (auto b : bcMapByBlockName) {
      for (auto &reg_name : bc_regex_vec) {
        if (std::regex_match(b.first, reg_name)) {
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
      MOFEM_LOG("BcMngSelf", Sev::noisy)
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

} // namespace MoFEM
