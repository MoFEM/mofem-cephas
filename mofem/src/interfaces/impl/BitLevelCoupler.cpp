/** \file BitLevelCoupler.cpp
 * \brief BitLevelCoupler interface implementation

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

MoFEMErrorCode
BitLevelCoupler::query_interface(const MOFEMuuid &uuid,
                                 UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMBitLevelCoupler) {
    *iface = const_cast<BitLevelCoupler *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::buildTree(const BitRefLevel &parent_level,
                                          int verb) {
  MoFEMFunctionBeginHot;
  Interface &m_field = cOre;
  treePtr.reset(new AdaptiveKDTree(&m_field.get_moab()));
  Range tets;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
      parent_level, BitRefLevel().set(), MBTET, tets);
  CHKERRG(ierr);
  rval = treePtr->build_tree(tets);
  CHKERRQ_MOAB(rval);
  if (verb > 2) {
    rval = treePtr->print();
    CHKERRQ_MOAB(rval);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::resetTree(const BitRefLevel &parent_level,
                                          int verb) {
  MoFEMFunctionBeginHot;
  treePtr->reset_tree();
  treePtr.reset();
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::getParent(const double *coords,
                                          EntityHandle &parent, bool tet_only,
                                          const double iter_tol,
                                          const double inside_tol, int verb) {
  MoFEMFunctionBeginHot;
  Interface &m_field = cOre;
  EntityHandle leaf_out;
  rval = treePtr->point_search(coords, leaf_out, iter_tol, inside_tol);
  CHKERRQ_MOAB(rval);
  bool is_in;
  Range tets;
  ierr = m_field.get_moab().get_entities_by_type(leaf_out, MBTET, tets);
  CHKERRG(ierr);
  Range::iterator tit = tets.begin();
  for (; tit != tets.end(); tit++) {
    ierr = getLocCoordsOnTet(*tit, coords, verb);
    CHKERRG(ierr);
    is_in = true;
    for (int nn = 0; nn < 4; nn++) {
      if (N[nn] < -inside_tol || N[nn] > 1 + inside_tol) {
        is_in = false;
        break;
      }
      if (!is_in)
        break;
    }
    parent = 0;
    if (is_in) {
      if (!tet_only) {
        // vertices
        if (fabs(N[0] - 1) < inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "node 0 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 0, 0, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1] - 1) < inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "node 1 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 0, 1, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2] - 2) < inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "node 2 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 0, 2, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3] - 1) < inside_tol) {
          if (verb > 1)
            std::cout << "node 3 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 0, 3, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        // edges
        if (fabs(N[0]) > inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "edge 0 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 0, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "edge 1 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 1, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) > inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "edge 2 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 2, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) > inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "edge 3 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 3, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "edge 4 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 4, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "edge 5 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 1, 5, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        // faces
        if (fabs(N[0]) > inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) < inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "face 0 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 2, 0, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) < inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "face 1 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 2, 1, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) > inside_tol && fabs(N[1]) < inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) > inside_tol) {
          if (verb > 1)
            std::cout << "face 2 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 2, 2, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        if (fabs(N[0]) > inside_tol && fabs(N[1]) > inside_tol &&
            fabs(N[2]) > inside_tol && fabs(N[3]) < inside_tol) {
          if (verb > 1)
            std::cout << "face 3 found " << std::endl;
          rval = m_field.get_moab().side_element(*tit, 2, 2, parent);
          CHKERRQ_MOAB(rval);
          MoFEMFunctionReturnHot(0);
        }
        // set parent
        if (parent != 0) {
          break;
        }
      }
      if (verb > 1)
        std::cout << "tet found " << std::endl;
      parent = *tit;
      break;
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::buildAdjacenciesVerticesOnTets(
    const BitRefLevel &parent_level, Range &children, bool vertex_elements,
    const double iter_tol, const double inside_tol, bool throw_error,
    int verb) {
  Interface &m_field = cOre;
  auto refined_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;
  // build Tree
  bool init_tree = false;

  // find parents of all nodes, if node has no parent then tetrahedral
  // containing that node is searched  node on tetrahedra my by part of face or
  // edge on that tetrahedral, this need to be verified
  auto it = refined_ptr->get<Ent_mi_tag>().lower_bound(
      get_id_for_min_type<MBVERTEX>());
  auto hi_it = refined_ptr->get<Ent_mi_tag>().upper_bound(
      get_id_for_max_type<MBVERTEX>());

  for (; it != hi_it; it++) {

    // entity on parent level, can be parent to yourself
    if (((*it)->getBitRefLevel() & parent_level).any())
      continue;

    if (verb > 1) {
      std::cout << *it << " " << (*it)->getBitRefLevel() << std::endl;
    }

    // that vertex is on parent bit level, no need to process
    // check if vertex has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = (*it)->getParentEnt();
    const auto ref_parent_ent = Core::makeSharedRefEntity(m_field, parent_ent);
    if ((ref_parent_ent->getBitRefLevel() & parent_level).any())
      continue;
    

    // check if vertex is on child entities set
    EntityHandle node = (*it)->getEnt();
    if (children.find(node) == children.end()) {
      continue;
    }

    // build a boundary volume Tree
    // if(!treePtr) {
    ierr = buildTree(parent_level, verb);
    CHKERRG(ierr);
    init_tree = true;
    //}

    double coords[3];
    rval = m_field.get_moab().get_coords(&node, 1, coords);
    CHKERRQ_MOAB(rval);
    EntityHandle parent = 0;
    ierr = getParent(coords, parent, false, iter_tol, inside_tol, verb);
    CHKERRG(ierr);
    ierr = chanegParent(refined_ptr->project<0>(it), parent);
    CHKERRG(ierr);
    if (throw_error && parent == 0) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "tets or any other entity for node not found");
    }
  }

  if (init_tree) {
    treePtr->reset_tree();
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::buildAdjacenciesEdgesFacesVolumes(
    const BitRefLevel &parent_level, Range &children, bool elements, int verb) {
  Interface &m_field = cOre;
  auto refined_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;

  if (verb > 2)
    std::cout << children << std::endl;

  std::vector<EntityHandle> conn_parents;

  Range::iterator eit, hi_eit;
  eit = children.begin();
  hi_eit = children.end();
  for (; eit != hi_eit; eit++) {

    // check entity type
    EntityType type;
    type = m_field.get_moab().type_from_handle(*eit);
    switch (type) {
    case MBEDGE:
    case MBTRI:
    case MBTET:
      break;
    default:
      continue;
    }

    // get ref entity iterator
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator it;
    it = refined_ptr->get<Ent_mi_tag>().find(*eit);
    // that entity is on parent bit level, no need to process
    if (((*it)->getBitRefLevel() & parent_level).any()) {
      continue;
    }

    if (verb > 1) {
      std::cout << "before: " << **it << std::endl;
    }

    // check if entity has a parent and parent is on parent bit level
    EntityHandle parent_ent;
    parent_ent = (*it)->getParentEnt();
    const auto ref_parent_ent = Core::makeSharedRefEntity(m_field, parent_ent);
    if ((ref_parent_ent->getBitRefLevel() & parent_level).any()) {
      if (!vErify)
        continue;
    }

    // connectivity
    int max_dim = m_field.get_moab().dimension_from_handle(*eit);
    const EntityHandle *conn;
    int num_nodes;
    ierr = m_field.get_moab().get_connectivity(*eit, conn, num_nodes);
    CHKERRG(ierr);
    conn_parents.resize(num_nodes);
    for (int nn = 0; nn < num_nodes; nn++) {
      conn_parents[nn] =
          Core::makeSharedRefEntity(m_field, conn[nn])->getParentEnt();
      RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator cit;
      cit = refined_ptr->get<Ent_mi_tag>().find(conn_parents[nn]);
      if (cit == refined_ptr->end()) {
        conn_parents[nn] = conn[nn];
        cit = refined_ptr->get<Ent_mi_tag>().find(conn_parents[nn]);
      }
      if (((*cit)->getBitRefLevel() & parent_level).none()) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "parent of vertex is not on parent bit level");
      }
      int ent_dim = m_field.get_moab().dimension_from_handle(conn_parents[nn]);
      max_dim = ent_dim > max_dim ? ent_dim : max_dim;
    }

    if (verb > 1)
      std::cout << "max_dim " << max_dim << std::endl;

    if (max_dim > 0) {

      for (; max_dim <= 3; max_dim++) {
        Range parent_ents;
        rval = m_field.get_moab().get_adjacencies(
            &*conn_parents.begin(), num_nodes, max_dim, false, parent_ents);
        CHKERRQ_MOAB(rval);
        parent_ents.erase((*it)->getEnt());
        if (!parent_ents.empty()) {
          ierr = chanegParent(refined_ptr->project<0>(it), *parent_ents.begin());
          CHKERRG(ierr);
          if (verb > 1) {
            std::cout << "after: " << **it << std::endl << std::endl;
          }
          break;
        }
      }

      if (!vErify && max_dim > 3) {
        ierr = chanegParent(refined_ptr->project<0>(it), 0);
        CHKERRG(ierr);
        if (verb > 1) {
          std::cout << "parent not found\n";
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::chanegParent(RefEntity_multiIndex::iterator it,
                                             EntityHandle parent) {
  Interface &m_field = cOre;
  auto refined_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;


  if (vErify) {
    ierr = verifyParent(it, parent);
    CHKERRG(ierr);
  }

  RefEntity_change_parent modifier(parent);
  bool success =
      const_cast<RefEntity_multiIndex *>(refined_ptr)->modify(it, modifier);
  if (!success) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "unsuccessful operation");
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::resetParents(Range &children, bool elements,
                                             int verb) {
  Interface &m_field = cOre;
  auto refined_ptr = m_field.get_ref_ents();
  MoFEMFunctionBeginHot;

  // access to ref dofs multi-index
  Range::iterator eit, hi_eit;
  eit = children.begin();
  hi_eit = children.end();
  for (; eit != hi_eit; eit++) {

    // get ref entity iterator
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator it;
    it = refined_ptr->get<Ent_mi_tag>().find(*eit);

    // resent entity parent
    ierr = chanegParent(refined_ptr->project<0>(it), 0);
    CHKERRG(ierr);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::verifyParent(RefEntity_multiIndex::iterator it,
                                             EntityHandle parent) {
  MoFEMFunctionBeginHot;

  if (parent != (*it)->getParentEnt()) {
    SETERRQ3(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "data inconsistency %lu != %lu for ent %lu", parent,
             (*it)->getParentEnt(), (*it)->getEnt());
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::getLocCoordsOnTet(EntityHandle tet,
                                                  const double *glob_coords,
                                                  int verb) {
  Interface &m_field = cOre;
  MoFEMFunctionBeginHot;

  int num_nodes;
  rval = m_field.get_moab().get_connectivity(tet, cOnn, num_nodes, true);
  CHKERRQ_MOAB(rval);
  rval = m_field.get_moab().get_coords(cOnn, num_nodes, cOords);
  CHKERRQ_MOAB(rval);
  double shifted_glob_coors[3];
  cblas_dcopy(3, glob_coords, 1, shifted_glob_coors, 1);
  for (int nn = 1; nn < 4; nn++) {
    cblas_daxpy(3, -1, cOords, 1, &cOords[nn * 3], 1);
  }
  cblas_daxpy(3, -1, cOords, 1, shifted_glob_coors, 1);
  for (int dd = 0; dd < 3; dd++) {
    cOords[dd] = 0;
  }
  ierr = ShapeDiffMBTET(diffN);
  CHKERRG(ierr);
  locCoords[0] = locCoords[1] = locCoords[2] = 0;
  ierr = ShapeMBTET(N, &locCoords[0], &locCoords[1], &locCoords[2], 1);
  CHKERRG(ierr);
  ierr = ShapeMBTET_inverse(N, diffN, cOords, shifted_glob_coors, locCoords);
  CHKERRG(ierr);
  ierr = ShapeMBTET(N, &locCoords[0], &locCoords[1], &locCoords[2], 1);
  CHKERRG(ierr);

  if (verb > 1) {
    std::cout << "N " << N[0] << " " << N[1] << " " << N[2] << " " << N[3]
              << std::endl;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode BitLevelCoupler::copyFieldDataFromParentToChildren(
    const std::vector<EntityHandle> &parents,
    const std::vector<EntityHandle> &children, const bool verify) {
  Interface &m_field = cOre;
  moab::Interface &moab = m_field.get_moab();
  auto fields_ptr = m_field.get_fields();
  MoFEMFunctionBegin;

  if (parents.size() != children.size()) {
    SETERRQ2(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
             "parents adn children vectors has to have the same size %d != %d",
             parents.size(), children.size());
  }

  const int nb_elems = parents.size();
  std::vector<double *> data(nb_elems);
  std::vector<int> data_size(nb_elems);

  for (auto fit = fields_ptr->begin(); fit != fields_ptr->end(); fit++) {

    // Verify consistency with database
    if (verify) {
      // Get pointer to multi-index with field entities
      auto *field_ents = m_field.get_field_ents();

      auto get_ents_max_order = [&](const std::vector<EntityHandle> &ents) {
        boost::shared_ptr<std::vector<const void *>> ents_max_order(
            new std::vector<const void *>());
        ents_max_order->resize(ents.size());
        CHKERR moab.tag_get_by_ptr((*fit)->th_AppOrder, &*ents.begin(),
                                   ents.size(), &*ents_max_order->begin());
        return ents_max_order;
      };

      auto max_order_parents = get_ents_max_order(parents);
      auto max_order_children = get_ents_max_order(children);

      auto pit = parents.begin();
      auto cit = children.begin();
      auto vit_parent_max_order = max_order_parents->begin();
      auto vit_child_max_order = max_order_children->begin();

      for (; pit != parents.end();
           ++pit, ++cit, ++vit_parent_max_order, ++vit_child_max_order) {
        // verify entity type
        if (moab.type_from_handle(*pit) != moab.type_from_handle(*cit)) {
          SETERRQ(m_field.get_comm(), MOFEM_DATA_INCONSISTENCY,
                  "inconsistent type");
        }
        // ref ent pointers
        auto ref_parent_ptr = Core::makeSharedRefEntity(m_field, *pit);
        auto ref_child_ptr = Core::makeSharedRefEntity(m_field, *cit);
        // create mofem entity objects
        boost::shared_ptr<FieldEntity> mofem_ent_parent(new FieldEntity(
            *fit, ref_parent_ptr,
            FieldEntity::makeSharedFieldDataAdaptorPtr(*fit, ref_parent_ptr),
            boost::shared_ptr<const int>(
                max_order_parents,
                static_cast<const int *>(*vit_parent_max_order))));
        boost::shared_ptr<FieldEntity> mofem_ent_child(new FieldEntity(
            *fit, ref_child_ptr,
            FieldEntity::makeSharedFieldDataAdaptorPtr(*fit, ref_child_ptr),
            boost::shared_ptr<const int>(
                max_order_children,
                static_cast<const int *>(*vit_child_max_order))));
        // check approximation order
        if (mofem_ent_parent->getMaxOrder() == mofem_ent_child->getMaxOrder()) {
          // approximation order is equal, simply copy data
          for (unsigned int dd = 0;
               dd != mofem_ent_child->getEntFieldData().size(); dd++) {
            mofem_ent_child->getEntFieldData()[dd] =
                mofem_ent_parent->getEntFieldData()[dd];
          }
        } else {
          // approximation orders is different
          FieldEntity_multiIndex::iterator fcit =
              field_ents->find(mofem_ent_child->getGlobalUniqueId());
          if (fcit == field_ents->end()) {
            // entity not in database, set order and copy data
            (FieldEntity_change_order(mofem_ent_parent->getMaxOrder()))(
                mofem_ent_child);
            for (unsigned int dd = 0;
                 dd != mofem_ent_child->getEntFieldData().size(); dd++) {
              mofem_ent_child->getEntFieldData()[dd] =
                  mofem_ent_parent->getEntFieldData()[dd];
            }
          } else {
            int parent_nb_dofs = mofem_ent_parent->getEntFieldData().size();
            int child_nb_dofs = mofem_ent_child->getEntFieldData().size();
            for (int dd = 0; dd != child_nb_dofs; dd++) {
              if (dd < parent_nb_dofs) {
                mofem_ent_child->getEntFieldData()[dd] =
                    mofem_ent_parent->getEntFieldData()[dd];
              } else {
                mofem_ent_child->getEntFieldData()[dd] = 0;
              }
            }
          }
        }
      }
    } else {

      auto copy_tag_data = [this, &moab, &data, &data_size](auto th, auto &p,
                                                            auto &c) {
        MoFEMFunctionBegin;
        // Get pointer and size of field values tag
        CHKERR moab.tag_get_by_ptr(th, p, (const void **)&data.front(),
                                   &data_size.front());
        // Set data
        CHKERR moab.tag_set_by_ptr(th, c, (void *const *)&data.front(),
                                   &data_size.front());
        MoFEMFunctionReturn(0);
      };

      Range p, c;
      p.insert_list(parents.begin(), parents.end());
      c.insert_list(children.begin(), children.end());

      Range parents_verts = p.subset_by_type(MBTET);
      Range children_verts = c.subset_by_type(MBTET);
      Range parents_without_verts = subtract(p, parents_verts);
      Range children_without_verts = subtract(c, children_verts);

      CHKERR copy_tag_data(fit->get()->th_FieldDataVerts, parents_verts,
                           children_verts);
      CHKERR copy_tag_data(fit->get()->th_FieldData, parents_verts,
                           children_verts);

    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode BitLevelCoupler::copyFieldDataFromParentToChildren(
    const BitRefLevel bit, const BitRefLevel mask, const bool verify) {
  Interface &m_field = cOre;
  // moab::Interface& moab = m_field.get_moab();
  MoFEMFunctionBeginHot;
  Range ents;
  ierr = m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(bit, mask,
                                                                      ents);
  CHKERRG(ierr);
  std::vector<EntityHandle> parents;
  std::vector<EntityHandle> children;
  for (Range::iterator eit = ents.begin(); eit != ents.end(); eit++) {
    const auto ref_ent = Core::makeSharedRefEntity(m_field, *eit);
    if (ref_ent->getParentEntType() == ref_ent->getEntType()) {
      children.push_back(*eit);
      parents.push_back(ref_ent->getParentEnt());
    }
  }
  ierr = copyFieldDataFromParentToChildren(parents, children, verify);
  CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
