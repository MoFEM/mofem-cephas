/** \file BcManager.cpp
 * \brief Manages boundary conditions
 * \ingroup bc_manager
 */

namespace MoFEM {

namespace BcManagerImplTools {

auto get_dim(const Range &ents) {
  for (auto d : {3, 2, 1})
    if (ents.num_of_dimension(d))
      return d;
  return 0;
};

auto get_adj_ents(moab::Interface &moab, const Range &ents) {
  Range verts;
  CHK_MOAB_THROW(moab.get_connectivity(ents, verts, true), "get verts");
  const auto dim = get_dim(ents);
  for (size_t d = 1; d < dim; ++d) {
    for (auto dd = d + 1; dd <= dim; ++dd) {
      CHK_MOAB_THROW(moab.get_adjacencies(ents.subset_by_dimension(dd), d,
                                          false, verts, moab::Interface::UNION),
                     "get adj");
    }
  }
  verts.merge(ents);
  return verts;
}
} // namespace BcManagerImplTools

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

  // if(problem_name.size())
  //   MOFEM_LOG("BcMngWorld", Sev::warning)
  //       << "Argument problem_name has no effect";

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
        CHKERR m->getAttributes(bc->bcAttributes);

        if (problem_name.size())
          CHKERR mark_fix_dofs(bc->bcMarkers, lo, hi);
        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block " << m->getName() << " number of attributes "
            << bc->bcAttributes.size();

        if (get_low_dim_ents) {
          auto low_dim_ents =
              BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
          bc->bcEnts.swap(low_dim_ents);
        }

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
        if (problem_name.size())
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

MoFEMErrorCode BcManager::addBlockDOFsToMPCs(const std::string problem_name,
                                             const std::string field_name,
                                             bool get_low_dim_ents,
                                             bool block_name_field_prefix,
                                             bool is_distributed_mesh) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";
  if (is_distributed_mesh)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument is_distributed_mesh=true has no effect";
  if (get_low_dim_ents)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument get_low_dim_ents=true has no effect";

  auto get_dim = [&](const Range &ents) {
    for (auto d : {3, 2, 1})
      if (ents.num_of_dimension(d))
        return d;
    return 0;
  };

  auto iterate_mpc_meshsets = [&]() {
    MoFEMFunctionBegin;

    auto mpc_meshset_ptr =
        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
            std::regex((boost::format("%s(.*)") % "MPC_(.*)").str()));

    for (auto m : mpc_meshset_ptr) {

      if (std::regex_match(m->getName(),
                           std::regex("(.*)COUPLING_LINKS(.*)"))) {

        auto bc = boost::make_shared<BCs>();
        bc->mpcPtr = boost::make_shared<MPCsType>();
        bc->mpcPtr->mpcType = MPC::COUPLING;

        std::string const corresponding_master_ms =
            std::regex_replace(m->getName(), std::regex("LINKS"), "MASTER");

        Range links_ents;
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         links_ents, true);

        Range master_nodes;
        if (m_field.getInterface<MeshsetsManager>()->checkMeshset(
                corresponding_master_ms)) {
          const CubitMeshSets *l;
          CHKERR m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
              corresponding_master_ms, &l);
          bc->mpcPtr->isReprocitical = false;

          CHKERR m_field.get_moab().get_entities_by_handle(l->getMeshset(),
                                                           master_nodes, true);
          // if (master_nodes.subset_by_dimension(0).size() <
          // links_ents.subset_by_dimension(1))
          {
            auto low_dim_ents = BcManagerImplTools::get_adj_ents(
                m_field.get_moab(), master_nodes);
            low_dim_ents = low_dim_ents.subset_by_dimension(0);
            master_nodes.swap(low_dim_ents);
          }

          MOFEM_LOG("BcMngWorld", Sev::verbose)
              << "Found block MASTER LINKS block: " << l->getName()
              << " Entities size: " << master_nodes.size();

        } else {
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "MASTER LINKS block not found: " << corresponding_master_ms
              << " setting reprocitical constraint. ("
              << bc->mpcPtr->isReprocitical << ").";
        }

        // if (std::regex_match(bc_id, std::regex("(.*)TIE(.*)")))
        //   bc->mpcPtr->mpcType = MPC::TIE;
        // if (std::regex_match(bc_id, std::regex("(.*)RIGID_BODY(.*)")))
        //   bc->mpcPtr->mpcType = MPC::RIGID_BODY;

        // std::cout << "links_ents range is: " << links_ents << "\n";
        for (auto &link : links_ents.subset_by_dimension(1)) {
          Range verts;
          CHKERR m_field.get_moab().get_connectivity(&link, 1, verts, true);
          // std::cout << "verts range is: " << verts << "\n";
          if (bc->mpcPtr->isReprocitical) {
            bc->mpcPtr->mpcMasterEnts.insert(verts[0]);
            bc->mpcPtr->mpcSlaveEnts.insert(verts[1]);
          } else {
            for (auto &m_node : verts)
              if (master_nodes.find(m_node) != master_nodes.end()) {
                // std::cout << "found ent: " << m_node << "\n";
                bc->mpcPtr->mpcMasterEnts.insert(m_node);
                bc->mpcPtr->mpcSlaveEnts.merge(
                    subtract(verts, Range(m_node, m_node)));
                break;
              }
          }
        }

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block MPC LINKS block: " << m->getName()
            << " Entities size (edges): " << links_ents.size()
            << " Entities size (nodes): "
            << bc->mpcPtr->mpcMasterEnts.size() +
                   bc->mpcPtr->mpcSlaveEnts.size()
            << " (" << bc->mpcPtr->mpcMasterEnts.size() << " "
            << bc->mpcPtr->mpcSlaveEnts.size() << ")";

        MOFEM_LOG("BcMngSync", Sev::noisy) << *bc->mpcPtr;
        MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

        // CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
        //     bc->mpcPtr->mpcMasterEnts);
        // CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
        //     bc->mpcPtr->mpcSlaveEnts);
        vector<double> mAttributes;
        CHKERR m->getAttributes(mAttributes);

        auto setFlags = [&](const auto &flags) {
          auto &d = bc->mpcPtr->data;
          if (flags.empty()) {
            d.flag1 = d.flag2 = d.flag3 = d.flag4 = d.flag5 = d.flag6 = true;
            return;
          }
          for (size_t i = 0; i < std::min(flags.size(), size_t(6)); ++i)
            (&d.flag1)[i] = flags[i] > 0.0;
        };

        setFlags(mAttributes);

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block " << m->getName() << " number of entities "
            << bc->mpcPtr->mpcMasterEnts.size() << " number of attributes "
            << mAttributes.size() << " highest dim of entities "
            << get_dim(bc->mpcPtr->mpcMasterEnts);

        // NOTE: we are not using markers at the moment
        // for (int i = 0; i != mAttributes.size(); ++i) {
        //   if (mAttributes[i] > 0.0)
        //     CHKERR mark_fix_dofs(bc->mpcPtr->bcMasterMarkers, i, i);
        // }

        // if (mAttributes.empty())
        //   CHKERR mark_fix_dofs(bc->mpcPtr->bcMasterMarkers, 0,
        //                        MAX_DOFS_ON_ENTITY);

        const std::string bc_id =
            problem_name + "_" + field_name + "_" + m->getName();
        bcMapByBlockName[bc_id] = bc;
      }
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR iterate_mpc_meshsets();

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

Range BcManager::getMergedBlocksRange(std::vector<std::regex> bc_regex_vec) {
  Range ents;
  if (bcMapByBlockName.size()) {
    for (auto b : bcMapByBlockName) {
      for (auto &reg_name : bc_regex_vec) {
        if (std::regex_match(b.first, reg_name)) {
          ents.merge(b.second->bcEnts);
        }
      }
    }
  }
  return ents;
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
BcManager::removeBlockDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  std::array<Range, 3> ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
           BLOCKSET | UNKNOWNNAME)) {

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

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
           BLOCKSET | UNKNOWNNAME)) {

    std::string bc_id = problem_name + "_" + field_name + "_" + m->getName();

    auto str = boost::format("%s_%s_%s(.*)")

               % problem_name % field_name % block_name;

    if (std::regex_match(bc_id, std::regex(str.str()))) {

      auto bc = bcMapByBlockName.at(bc_id);

      if (auto temp_bc = bc->tempBcPtr) {
        if (temp_bc->data.flag1) {
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

  // if(problem_name.size()==0)
  //   MOFEM_LOG("BcMngWorld", Sev::warning)
  //       << "Argument problem_name has no effect";

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";

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
            << "Found block DISPLACEMENTSET id = " << m->getMeshsetId();
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->dispBcPtr;

        MOFEM_LOG("BcMngSync", Sev::noisy)
            << "Found block DISPLACEMENTSET id = " << m->getMeshsetId()
            << " nb. of entities " << bc->bcEnts.size()
            << " highest dim of entities "
            << BcManagerImplTools::get_dim(bc->bcEnts);
        MOFEM_LOG("BcMngSync", Sev::noisy) << *bc->dispBcPtr;
        MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

        if (problem_name.size()) {

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
        }

        if (get_low_dim_ents) {
          auto low_dim_ents =
              BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
          bc->bcEnts.swap(low_dim_ents);
        }

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
        if (problem_name.size())
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
            << "Found block TEMPERATURESET id = " << m->getMeshsetId();
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->tempBcPtr;

        CHKERR prb_mng->modifyMarkDofs(
            problem_name, ROW, field_name, 0, MAX_DOFS_ON_ENTITY,
            ProblemsManager::MarkOP::OR, 1, bc->bcMarkers);

        if (get_low_dim_ents) {
          auto low_dim_ents =
              BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
          bc->bcEnts.swap(low_dim_ents);
        }

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
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

        CHKERR prb_mng->modifyMarkDofs(
            problem_name, ROW, field_name, 0, MAX_DOFS_ON_ENTITY,
            ProblemsManager::MarkOP::OR, 1, bc->bcMarkers);

        MOFEM_LOG("BcMngWorld", Sev::inform)
            << "Found block HEATFLUX id = " << m->getMeshsetId();
        MOFEM_LOG("BcMngWorld", Sev::inform) << *bc->heatFluxBcPtr;

        if (get_low_dim_ents) {
          auto low_dim_ents =
              BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
          bc->bcEnts.swap(low_dim_ents);
        }

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
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
MoFEMErrorCode
BcManager::pushMarkDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
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
        // for the ROTATE_X block the angles can be specified with either one or
        // three attributes, e.g. 1, coords or 1,0,0,coords
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value4 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block on block (angle (1 or 3), "
                 "center coords(3) but have "
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
        // for the ROTATE_Y block the angles can be specified with either one or
        // three attributes, e.g. 1, coords or 0,1,0,coords
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value5 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block on block (angle (1 or 3), "
                 "center coords(3) but have "
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
        // for the ROTATE_Z block the angles can be specified with either one or
        // three attributes, e.g. 1, coords or 0,0,1,coords
        if (bc->bcAttributes.empty()) {
          bc->dispBcPtr->data.value6 = 0;
          MOFEM_LOG("BcMngWorld", Sev::warning)
              << "Expected one attribute on block (angle (1 or 3), center "
                 "coords(3) but have "
              << bc->bcAttributes.size();
        } else if (bc->bcAttributes.size() == 1 ||
                   bc->bcAttributes.size() == 4) {
          bc->dispBcPtr->data.value6 = bc->bcAttributes[0];
        } else if (bc->bcAttributes.size() == 3 ||
                   bc->bcAttributes.size() == 6) {
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
  CHKERR pushMarkDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
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
  CHKERR removeBlockDOFsOnEntities<BcDisplacementMeshsetType<BLOCKSET>>(
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

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcMeshsetType<FORCESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  if (problem_name.empty())
    MOFEM_LOG("BcMngWorld", Sev::warning) << "Argument problem_name is empty";

  if (block_name_field_prefix)
    MOFEM_LOG("BcMngWorld", Sev::warning)
        << "Argument block_name_field_prefix=true has no effect";

  auto fix_force = [&]() {
    MoFEMFunctionBegin;

    auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
      MoFEMFunctionBegin;
      for (auto m : meshset_vec_ptr) {
        auto bc = boost::make_shared<BCs>();
        CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                         bc->bcEnts, true);
        bc->forceBcPtr = boost::make_shared<ForceCubitBcData>();
        CHKERR m->getBcDataStructure(*(bc->forceBcPtr));

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block FORCESET id = " << m->getMeshsetId();
        MOFEM_LOG("BcMngWorld", Sev::verbose) << *bc->forceBcPtr;

        MOFEM_LOG("BcMngSync", Sev::noisy)
            << "Found block FORCESET id = " << m->getMeshsetId()
            << " nb. of entities " << bc->bcEnts.size()
            << " highest dim of entities "
            << BcManagerImplTools::get_dim(bc->bcEnts);
        MOFEM_LOG("BcMngSync", Sev::noisy) << *bc->forceBcPtr;
        MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

        if (problem_name.size()) {

          if (bc->forceBcPtr->data.value2 > 0)
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
          if (bc->forceBcPtr->data.value3 > 0)
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
          if (bc->forceBcPtr->data.value4 > 0)
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);

          if (bc->forceBcPtr->data.value5 > 0) {

            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
          }
          if (bc->forceBcPtr->data.value5) {

            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
          }
          if (bc->forceBcPtr->data.value6) {

            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
            CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                           ProblemsManager::MarkOP::OR, 1,
                                           bc->bcMarkers);
          }
        }

        if (get_low_dim_ents) {
          auto low_dim_ents =
              BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
          bc->bcEnts.swap(low_dim_ents);
        }

        CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
            bc->bcEnts);
        if (problem_name.size())
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEnts, bc->bcMarkers);

        const std::string bc_id =
            problem_name + "_" + field_name + "_FORCESET" +
            boost::lexical_cast<std::string>(m->getMeshsetId());
        bcMapByBlockName[bc_id] = bc;
      }
      MoFEMFunctionReturn(0);
    };

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(NODESET |
                                                                    FORCESET)

    );

    MoFEMFunctionReturn(0);
  };

  CHKERR fix_force();

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  if (problem_name.size() == 0)
    MOFEM_LOG("BcMngWorld", Sev::warning) << "Argument problem_name is empty";

  auto get_force_block = [&](auto block_name) {
    MoFEMFunctionBegin;

    for (auto m :
         m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

             (boost::format("%s(.*)") % block_name).str()

                 ))

    ) {

      const auto block_name = m->getName();

      MOFEM_LOG("BcMngWorld", Sev::inform)
          << "Found force block " << block_name;

      auto bc = boost::make_shared<BCs>();
      CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                       bc->bcEnts, true);

      CHKERR m->getAttributes(bc->bcAttributes);
      if (bc->bcAttributes.size() < 3) {
        SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                "At least three block attributes for force block are expected");
      }

      bc->forceBcPtr = boost::make_shared<ForceCubitBcData>();
      // For details look at ForceCubitBcData in
      // mofem/src/multi_indices/BCData.hpp
      bc->forceBcPtr->data.value1 = 1;
      bc->forceBcPtr->data.value3 = bc->bcAttributes[0];
      bc->forceBcPtr->data.value4 = bc->bcAttributes[1];
      bc->forceBcPtr->data.value5 = bc->bcAttributes[2];

      MOFEM_LOG("BcMngWorld", Sev::inform) << *bc->forceBcPtr;
      MOFEM_LOG("BcMngSync", Sev::noisy)
          << "Found block FORCESET id = " << m->getMeshsetId()
          << " nb. of entities " << bc->bcEnts.size()
          << " highest dim of entities "
          << BcManagerImplTools::get_dim(bc->bcEnts);
      MOFEM_LOG("BcMngSync", Sev::noisy) << *bc->forceBcPtr;
      MOFEM_LOG_SEVERITY_SYNC(m_field.get_comm(), Sev::noisy);

      if (problem_name.size()) {

        if (bc->forceBcPtr->data.value2 > 0)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 0, 0,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        if (bc->forceBcPtr->data.value3 > 0)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 1, 1,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
        if (bc->forceBcPtr->data.value4 > 0)
          CHKERR prb_mng->modifyMarkDofs(problem_name, ROW, field_name, 2, 2,
                                         ProblemsManager::MarkOP::OR, 1,
                                         bc->bcMarkers);
      }

      if (get_low_dim_ents) {
        auto low_dim_ents =
            BcManagerImplTools::get_adj_ents(m_field.get_moab(), bc->bcEnts);
        bc->bcEnts.swap(low_dim_ents);
      }

      CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
          bc->bcEnts);
      if (problem_name.size())
        CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                 bc->bcEnts, bc->bcMarkers);

      const std::string bc_id =
          problem_name + "_" + field_name + "_" + block_name;
      bcMapByBlockName[bc_id] = bc;
    }
    MoFEMFunctionReturn(0);
  };

  CHKERR get_force_block("FORCE");

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::pushMarkDOFsOnEntities<ForceCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix) {
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<FORCESET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  CHKERR pushMarkDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix);

  MoFEMFunctionReturn(0);
}

template <>
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<BcMeshsetType<FORCESET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcMeshsetType<FORCESET>>(
      problem_name, field_name, get_low_dim_ents);

  std::array<Range, 3> ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(NODESET |
                                                                   FORCESET)) {

    const auto block_name = m->getName();
    std::string bc_id = problem_name + "_" + field_name + "_" + block_name;

    auto str = boost::format("%s_%s_%s(.*)")

               % problem_name % field_name % block_name;

    if (std::regex_match(bc_id, std::regex(str.str()))) {

      auto bc = bcMapByBlockName.at(bc_id);

      if (auto force_bc = bc->forceBcPtr) {
        if (force_bc->data.value3 > 0) {
          ents_to_remove[0].merge(bc->bcEnts);
        }
        if (force_bc->data.value4 > 0) {
          ents_to_remove[1].merge(bc->bcEnts);
        }
        if (force_bc->data.value5 > 0) {
          ents_to_remove[2].merge(bc->bcEnts);
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
BcManager::removeBlockDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  Interface &m_field = cOre;
  auto prb_mng = m_field.getInterface<ProblemsManager>();
  MoFEMFunctionBegin;

  CHKERR pushMarkDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents);

  std::array<Range, 3> ents_to_remove;

  for (auto m :

       m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(
           BLOCKSET | UNKNOWNNAME)) {

    const auto block_name = m->getName();
    std::string bc_id = problem_name + "_" + field_name + "_" + block_name;

    auto str = boost::format("%s_%s_%s(.*)")

               % problem_name % field_name % block_name;

    if (std::regex_match(bc_id, std::regex(str.str()))) {

      auto bc = bcMapByBlockName.at(bc_id);

      if (auto force_bc = bc->forceBcPtr) {
        if (force_bc->data.value3 > 0) {
          ents_to_remove[0].merge(bc->bcEnts);
        }
        if (force_bc->data.value4 > 0) {
          ents_to_remove[1].merge(bc->bcEnts);
        }
        if (force_bc->data.value5 > 0) {
          ents_to_remove[2].merge(bc->bcEnts);
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
MoFEMErrorCode BcManager::removeBlockDOFsOnEntities<ForceCubitBcData>(
    const std::string problem_name, const std::string field_name,
    bool get_low_dim_ents, bool block_name_field_prefix,
    bool is_distributed_mesh) {
  MoFEMFunctionBegin;

  CHKERR removeBlockDOFsOnEntities<BcMeshsetType<FORCESET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);

  CHKERR removeBlockDOFsOnEntities<BcForceMeshsetType<BLOCKSET>>(
      problem_name, field_name, get_low_dim_ents, block_name_field_prefix,
      is_distributed_mesh);

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode BcManager::pushMarkSideDofs(

    const std::string problem_name, const std::string block_name,
    const std::string field_name, int bridge_dim, int lo, int hi

) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  if (problem_name.empty())
    MoFEMFunctionReturnHot(0);

  auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
    auto prb_mng = m_field.getInterface<ProblemsManager>();
    MoFEMFunctionBegin;
    for (auto m : meshset_vec_ptr) {
      auto bc = boost::make_shared<BCs>();
      CHKERR m_field.get_moab().get_entities_by_handle(m->getMeshset(),
                                                       bc->bcEnts, true);
      CHKERR m->getAttributes(bc->bcAttributes);

      bc->dofsViewPtr = boost::make_shared<BCs::DofsView>();

      CHKERR m_field.getInterface<CommInterface>()->synchroniseEntities(
          bc->bcEnts);
      CHKERR prb_mng->getSideDofsOnBrokenSpaceEntities(
          *(bc->dofsViewPtr), problem_name, ROW, field_name, bc->bcEnts,
          bridge_dim, lo, hi);
      CHKERR prb_mng->markDofs(problem_name, ROW, *(bc->dofsViewPtr),
                               ProblemsManager::OR, bc->bcMarkers);

      MOFEM_LOG("BcMngWorld", Sev::inform)
          << "Found block " << m->getName() << " number of attributes "
          << bc->bcAttributes.size() << " number of entities "
          << bc->bcEnts.size();

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
}

MoFEMErrorCode BcManager::removeSideDOFs(const std::string problem_name,
                                         const std::string block_name,
                                         const std::string field_name,
                                         int bridge_dim, int lo, int hi,
                                         bool is_distributed_mesh) {
  Interface &m_field = cOre;
  MoFEMFunctionBegin;

  CHKERR pushMarkSideDofs(problem_name, block_name, field_name, bridge_dim, lo,
                          hi);

  auto iterate_meshsets = [&](auto &&meshset_vec_ptr) {
    MoFEMFunctionBegin;
    auto prb_mng = m_field.getInterface<ProblemsManager>();
    for (auto m : meshset_vec_ptr) {
      const std::string bc_id =
          problem_name + "_" + field_name + "_" + m->getName();
      auto &bc = bcMapByBlockName.at(bc_id);
      CHKERR prb_mng->removeDofs(problem_name, ROW, *(bc->dofsViewPtr), lo, hi);
      CHKERR prb_mng->removeDofs(problem_name, COL, *(bc->dofsViewPtr), lo, hi);
    }
    MoFEMFunctionReturn(0);
  };

  if (is_distributed_mesh) {

    CHKERR iterate_meshsets(

        m_field.getInterface<MeshsetsManager>()->getCubitMeshsetPtr(std::regex(

            (boost::format("%s(.*)") % block_name).str()

                ))

    );
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Not implemented");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
