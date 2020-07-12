/** \file FieldBlas.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_section_manager
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

MoFEMErrorCode FieldBlas::query_interface(const MOFEMuuid &uuid,
                                          UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMFieldBlas) {
    *iface = const_cast<FieldBlas *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

FieldBlas::FieldBlas(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}
FieldBlas::~FieldBlas() {}

MoFEMErrorCode FieldBlas::fieldLambda(FieldBlas::TwoFieldFunction lambda,
                                      const std::string field_name_x,
                                      const std::string field_name_y,
                                      bool error_if_missing,
                                      bool creat_if_missing) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto field_ents = m_field.get_field_ents();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;

  auto x_fit = fields_ptr->get<FieldName_mi_tag>().find(field_name_x);
  if (x_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "x field < %s > not found, (top tip: check spelling)",
             field_name_x.c_str());
  }
  auto y_fit = fields_ptr->get<FieldName_mi_tag>().find(field_name_y);
  if (y_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "y field < %s > not found, (top tip: check spelling)",
             field_name_y.c_str());
  }
  if ((*x_fit)->getSpace() != (*y_fit)->getSpace()) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "space for field < %s > and field <%s> are not compatible",
             field_name_x.c_str(), field_name_y.c_str());
  }
  if ((*x_fit)->getNbOfCoeffs() != (*y_fit)->getNbOfCoeffs()) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "rank for field < %s > and field <%s> are not compatible",
             field_name_x.c_str(), field_name_y.c_str());
  }

  typedef multi_index_container<
      boost::shared_ptr<DofEntity>,
      indexed_by<

          hashed_non_unique<
              tag<Composite_Ent_Order_And_CoeffIdx_mi_tag>,
              composite_key<

                  DofEntity,
                  const_mem_fun<DofEntity, EntityHandle, &DofEntity::getEnt>,
                  const_mem_fun<DofEntity, ApproximationOrder,
                                &DofEntity::getDofOrder>,
                  const_mem_fun<DofEntity, FieldCoefficientsNumber,
                                &DofEntity::getDofCoeffIdx>

                  >>

          >>
      DofEntity_multiIndex_composite_view;

  auto dof_lo_for_view = dofs_ptr->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*y_fit)->getBitNumber()));
  auto dof_hi_for_view = dofs_ptr->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*y_fit)->getBitNumber()));
      
  DofEntity_multiIndex_composite_view dof_composite_view;
  dof_composite_view.insert(dof_lo_for_view, dof_hi_for_view);

  auto x_eit = field_ents->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*x_fit)->getBitNumber()));
  auto x_eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*x_fit)->getBitNumber()));
  for (; x_eit != x_eit_hi; x_eit++) {
    int nb_dofs_on_x_entity = (*x_eit)->getNbDofsOnEnt();
    VectorAdaptor field_data = (*x_eit)->getEntFieldData();
    for (int dd = 0; dd != nb_dofs_on_x_entity; ++dd) {
      ApproximationOrder dof_order = (*x_eit)->getDofOrderMap()[dd];
      FieldCoefficientsNumber dof_rank = dd % (*x_eit)->getNbOfCoeffs();
      FieldData data = field_data[dd];
      auto dit = dof_composite_view.find(
          boost::make_tuple((*x_eit)->getEnt(), dof_order, dof_rank));
      if (dit == dof_composite_view.end()) {
        if (creat_if_missing) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                  "not yet implemented");
        } else {
          if (error_if_missing) {
            std::ostringstream ss;
            ss << "dof on ent " << (*x_eit)->getEnt() << " order " << dof_order
               << " rank " << dof_rank << " does not exist";
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    ss.str().c_str());
          } else {
            continue;
          }
        }
      }
      CHKERR lambda((*dit)->getFieldData(),data);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldAxpy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing,
                                    bool creat_if_missing) {
  MoFEMFunctionBegin;
  struct Axpy {
    const double alpha;
    Axpy(const double alpha) : alpha(alpha) {}
    inline MoFEMErrorCode operator()(double &fy, const double fx) {
      MoFEMFunctionBeginHot;
      fy += alpha * fx;
      MoFEMFunctionReturnHot(0);
    }
  };
  CHKERR fieldLambda(Axpy(alpha), field_name_x, field_name_y, error_if_missing,
                     creat_if_missing);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldCopy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing,
                                    bool creat_if_missing) {
  MoFEMFunctionBegin;
  struct Copy {
    const double alpha;
    Copy(const double alpha) : alpha(alpha) {}
    inline MoFEMErrorCode operator()(double &fy, const double fx) {
      MoFEMFunctionBeginHot;
      fy = alpha * fx;
      MoFEMFunctionReturnHot(0);
    }
  };
  CHKERR fieldLambda(Copy(alpha), field_name_x, field_name_y, error_if_missing,
                     creat_if_missing);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::setVertexDofs(FieldBlas::VertexCoordsFunction lambda,
                                        const std::string field_name,
                                        Range *sub_verts) {
  const MoFEM::Interface &m_field = cOre;
  MoFEMFunctionBegin;

  EntityHandle meshset = m_field.get_field_meshset(field_name);
  Range verts;
  CHKERR m_field.get_moab().get_entities_by_type(meshset, MBVERTEX, verts,
                                                 true);
  if (sub_verts)
    verts = intersect(*sub_verts, verts);

  struct LambdaMethod : EntityMethod {
    LambdaMethod(MoFEM::Interface &m_field, Range &verts,
                 FieldBlas::VertexCoordsFunction lambda)
        : EntityMethod(), mField(m_field), verts(verts), lambda(lambda),
          count(0), total(0) {}
    MoFEMErrorCode preProcess() {
      vit = verts.begin();
      return 0;
    }
    MoFEMErrorCode operator()() {
      MoFEMFunctionBegin;
      if (*vit != entPtr->getEnt())
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Inconsistent entity %ld != %ld", *vit, entPtr->getEnt());
      if(!count)
        CHKERR mField.get_moab().coords_iterate(vit, verts.end(), xCoord,
                                                yCoord, zCoord, count);
      CHKERR lambda(entPtr->getEntFieldData(), xCoord, yCoord, zCoord);
      ++xCoord;
      ++yCoord;
      ++zCoord;
      ++vit;
      ++total;
      --count;
      MoFEMFunctionReturn(0);
    }
    MoFEMErrorCode postProcess() {
      MoFEMFunctionBegin;
      if(total != verts.size())
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Inconsistent number of vertices in field meshset and in the "
                "field %d != %d",
                total, verts.size());
      MoFEMFunctionReturn(0);
    }

  private:
    MoFEM::Interface &mField;
    Range &verts;
    FieldBlas::VertexCoordsFunction lambda;
    int count;
    int total;
    Range::iterator vit;
    double *xCoord;
    double *yCoord;
    double *zCoord;
  };

  LambdaMethod lambda_method(const_cast<MoFEM::Interface &>(m_field), verts,
                             lambda);
  CHKERR const_cast<MoFEM::Interface &>(m_field).loop_entities(
      field_name, lambda_method, &verts, QUIET);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::setField(const double val, const EntityType type,
                                   const std::string field_name) {
  const MoFEM::Interface &m_field = cOre;
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBeginHot;

  auto dit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(
      boost::make_tuple(field_name, get_id_for_min_type(type)));
  auto hi_dit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(
      boost::make_tuple(field_name, get_id_for_max_type(type)));
  for (; dit != hi_dit; dit++) 
    (*dit)->getFieldData() = val;
  
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FieldBlas::setField(const double val, const EntityType type,
                                   const Range &ents,
                                   const std::string field_name) {
  const MoFEM::Interface &m_field = cOre;
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBeginHot;

  auto dit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(
      boost::make_tuple(field_name, get_id_for_min_type(type)));
  auto hi_dit = dofs_ptr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(
      boost::make_tuple(field_name, get_id_for_max_type(type)));
  EntityHandle ent, last = 0;
  bool cont = true;
  for (; dit != hi_dit; dit++) {
    ent = (*dit)->getEnt();
    if (ent != last) {
      if (ents.find(ent) == ents.end()) {
        cont = true;
      } else {
        cont = false;
      }
      last = ent;
    }
    if (cont)
      continue;
    (*dit)->getFieldData() = val;
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode FieldBlas::setField(const double val,
                                   const std::string field_name) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;
  auto fit = fields_ptr->get<FieldName_mi_tag>().find(field_name);
  if (fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             " field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }

  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*fit)->getBitNumber()));
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*fit)->getBitNumber()));
  for (; dit != hi_dit; dit++) {
    (*dit)->getFieldData() = val;
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldScale(const double alpha,
                                     const std::string field_name) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBeginHot;

  auto fit = fields_ptr->get<FieldName_mi_tag>().find(field_name);
  if (fit == fields_ptr->get<FieldName_mi_tag>().end()) 
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             " field < %s > not found, (top tip: check spelling)",
             field_name.c_str());

  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*fit)->getBitNumber()));
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*fit)->getBitNumber()));
  for (; dit != hi_dit; dit++) 
    (*dit)->getFieldData() *= alpha;
  
  MoFEMFunctionReturnHot(0);
}

} // namespace MoFEM
