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

MoFEMErrorCode
FieldBlas::query_interface(boost::typeindex::type_index type_index,
                           UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = const_cast<FieldBlas *>(this);
  MoFEMFunctionReturnHot(0);
}

FieldBlas::FieldBlas(const MoFEM::Core &core)
    : cOre(const_cast<MoFEM::Core &>(core)), dEbug(false) {}

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

  auto x_eit = field_ents->get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*x_fit)->getBitNumber()));
  auto x_eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*x_fit)->getBitNumber()));
  auto y_eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*y_fit)->getBitNumber()));

  for (; x_eit != x_eit_hi;) {

    const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
        (*y_fit)->getBitNumber(), (*x_eit)->getEnt());
    auto y_eit = field_ents->get<Unique_mi_tag>().find(lo_uid);

    if (y_eit == field_ents->end()) {

      if (creat_if_missing) {
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
                "Not implemented creation of DOFs on the fly");
      } else {
        if (error_if_missing) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Missing entity");
        }
      }

      ++x_eit;

    } else {

      auto check = [&]() {
        if (x_eit == x_eit_hi || x_eit == field_ents->end())
          return false;
        if (y_eit == y_eit_hi || y_eit == field_ents->end())
          return false;
        if ((*y_eit)->getEnt() != (*x_eit)->getEnt())
          return false;
        return true;
      };

      do {

        VectorAdaptor x_field_data = (*x_eit)->getEntFieldData();
        VectorAdaptor y_field_data = (*y_eit)->getEntFieldData();
        const auto size_x = x_field_data.size();
        const auto size_y = y_field_data.size();

        size_t dd = 0;
        for (; dd != std::min(size_x, size_y); ++dd)
          CHKERR lambda(y_field_data[dd], x_field_data[dd]);
        for (; dd < size_y; ++dd)
          y_field_data[dd] = 0;

        ++x_eit;
        ++y_eit;

      } while (check());
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
      if (!count)
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
      if (total != verts.size())
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

  const auto bit_number = m_field.get_field_bit_number(field_name);
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type));
  const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_max_type(type));

  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(lo_uid);
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(hi_uid);

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

  const auto bit_number = m_field.get_field_bit_number(field_name);
  const auto lo_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_min_type(type));
  const auto hi_uid = FieldEntity::getLocalUniqueIdCalculate(
      bit_number, get_id_for_max_type(type));

  auto dit = dofs_ptr->get<Unique_mi_tag>().lower_bound(lo_uid);
  auto hi_dit = dofs_ptr->get<Unique_mi_tag>().upper_bound(hi_uid);

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
    if (!cont)
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
