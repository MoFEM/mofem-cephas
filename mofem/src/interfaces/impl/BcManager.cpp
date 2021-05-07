/** \file BcManager.cpp
 * \brief Manages boundary conditions
 * \ingroup bc_manager
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

namespace MoFEM {

MoFEMErrorCode BcManager::query_interface(const MOFEMuuid &uuid,
                                          UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMBcManager) {
    *iface = const_cast<BcManager *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
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
  auto prb_mng = m_field.getInterface<ProblemsManager>() MoFEMFunctionBegin;

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
                                                         bc->bcEdges, true);
        CHKERR it->getAttributes(bc->bcAttributes);

        MOFEM_LOG("BcMngWorld", Sev::verbose)
            << "Found block " << block_name << " number of entities "
            << bc->bcEdges.size() << " number of attributes "
            << bc->bcAttributes.size() << " highest dim of entities "
            << get_dim(bc->bcEdges);

        CHKERR mark_fix_dofs(bc->bcMarkers, lo, hi);
        if (get_low_dim_ents)
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   get_adj_ents(bc->bcEdges), bc->bcMarkers);
        else
          CHKERR prb_mng->markDofs(problem_name, ROW, ProblemsManager::AND,
                                   bc->bcEdges, bc->bcMarkers);

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

} // namespace MoFEM
