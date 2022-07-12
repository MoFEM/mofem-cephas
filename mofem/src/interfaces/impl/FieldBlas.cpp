/** \file FieldBlas.cpp
 * \brief Managing complexities for problem
 * \ingroup mofem_section_manager
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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

MoFEMErrorCode
FieldBlas::fieldLambdaOnValues(FieldBlas::OneFieldFunctionOnValues lambda,
                               const std::string field_name, Range *ents_ptr) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  MoFEMFunctionBegin;

  auto x_fit = fields_ptr->get<FieldName_mi_tag>().find(field_name);
  if (x_fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }

  auto wrap_lambda_on_enties =
      [lambda](boost::shared_ptr<FieldEntity> field_entity_ptr) {
        MoFEMFunctionBeginHot;

        auto field_data = field_entity_ptr->getEntFieldData();
        for (auto &v : field_data)
          v = lambda(v);

        MoFEMFunctionReturnHot(0);
      };

  CHKERR fieldLambdaOnEntities(wrap_lambda_on_enties, field_name, ents_ptr);
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode
FieldBlas::fieldLambdaOnEntities(FieldBlas::OneFieldFunctionOnEntities lambda,
                                 const std::string field_name,
                                 Range *ents_ptr) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
  auto field_ents = m_field.get_field_ents();
  auto dofs_ptr = m_field.get_dofs();
  MoFEMFunctionBegin;

  auto fit = fields_ptr->get<FieldName_mi_tag>().find(field_name);
  if (fit == fields_ptr->get<FieldName_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "field < %s > not found, (top tip: check spelling)",
             field_name.c_str());
  }

  // End of iterator for y_field
  const auto bit_number = (*fit)->getBitNumber();

  auto loop_ents = [&](const EntityHandle f, const EntityHandle s) {
    MoFEMFunctionBeginHot;

    const auto lo_uid = FieldEntity::getLoLocalEntityBitNumber(bit_number, f);
    const auto hi_uid = FieldEntity::getHiLocalEntityBitNumber(bit_number, s);

    // Start of iterator for x_field
    auto eit = field_ents->get<Unique_mi_tag>().lower_bound(lo_uid);
    // End of iterator for x_field
    auto eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(hi_uid);

    for (; eit != eit_hi; ++eit) {
      CHKERR lambda(*eit);
    }

    MoFEMFunctionReturnHot(0);
  };

  if (ents_ptr) {
    for (auto p = ents_ptr->const_pair_begin(); p != ents_ptr->const_pair_end();
         ++p) {
      CHKERR loop_ents(p->first, p->second);
    }
  } else {
    // we are looping from the very first possible entity handle (MBVERTEX) to
    // the very last (MBENTITYSET)
    CHKERR loop_ents(get_id_for_min_type<MBVERTEX>(),
                     get_id_for_max_type<MBENTITYSET>());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldLambdaOnValues(
    FieldBlas::TwoFieldFunctionOnValues lambda, const std::string field_name_x,
    const std::string field_name_y, bool error_if_missing, Range *ents_ptr) {
  const MoFEM::Interface &m_field = cOre;
  auto fields_ptr = m_field.get_fields();
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
             "space for field < %s > and field <%s> are not compatible for "
             "fieldblas",
             field_name_x.c_str(), field_name_y.c_str());
  }
  if ((*x_fit)->getNbOfCoeffs() != (*y_fit)->getNbOfCoeffs()) {
    SETERRQ2(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "rank for field < %s > and field <%s> are not compatible for fieldblas",
        field_name_x.c_str(), field_name_y.c_str());
  }

  auto wrap_lambda_on_enties =
      [lambda](boost::shared_ptr<FieldEntity> y_field_entity_ptr,
               const boost::shared_ptr<FieldEntity> x_field_entity_ptr) {
        MoFEMFunctionBeginHot;

        auto x_field_data = x_field_entity_ptr->getEntFieldData();
        auto y_field_data = y_field_entity_ptr->getEntFieldData();
        const auto size_x = x_field_data.size();
        const auto size_y = y_field_data.size();

        size_t dd = 0;
        for (; dd != std::min(size_x, size_y); ++dd)
          y_field_data[dd] = lambda(y_field_data[dd], x_field_data[dd]);
        for (; dd < size_y; ++dd)
          y_field_data[dd] = 0;

        MoFEMFunctionReturnHot(0);
      };

  CHKERR fieldLambdaOnEntities(wrap_lambda_on_enties, field_name_x,
                               field_name_y, error_if_missing, ents_ptr);
  MoFEMFunctionReturn(0);
};

MoFEMErrorCode
FieldBlas::fieldLambdaOnEntities(FieldBlas::TwoFieldFunctionOnEntities lambda,
                                 const std::string field_name_x,
                                 const std::string field_name_y,
                                 bool error_if_missing, Range *ents_ptr) {
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

  // End of iterator for y_field
  const auto x_bit_number = (*x_fit)->getBitNumber();
  const auto y_bit_number = (*y_fit)->getBitNumber();
  const auto y_eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId(y_bit_number));

  auto loop_ents = [&](const EntityHandle f, const EntityHandle s) {
    MoFEMFunctionBeginHot;

    const auto x_lo_uid =
        FieldEntity::getLoLocalEntityBitNumber(x_bit_number, f);
    const auto x_hi_uid =
        FieldEntity::getHiLocalEntityBitNumber(x_bit_number, s);

    // Start of iterator for x_field
    auto x_eit = field_ents->get<Unique_mi_tag>().lower_bound(x_lo_uid);
    // End of iterator for x_field
    auto x_eit_hi = field_ents->get<Unique_mi_tag>().upper_bound(x_hi_uid);

    for (; x_eit != x_eit_hi;) {

      const auto y_lo_uid = FieldEntity::getLocalUniqueIdCalculate(
          (*y_fit)->getBitNumber(), (*x_eit)->getEnt());
      auto y_eit = field_ents->get<Unique_mi_tag>().find(y_lo_uid);

      if (y_eit == field_ents->end()) {

        if (error_if_missing) {
          SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                  "Missing entity in y_field.");
        }

        ++x_eit;

      } else {

        auto check = [&]() {
          if (x_eit == x_eit_hi)
            return false;
          if (y_eit == y_eit_hi)
            return false;
          if ((*y_eit)->getEnt() != (*x_eit)->getEnt())
            return false;
          return true;
        };

        do {

          CHKERR lambda(*y_eit, *x_eit);

          ++x_eit;
          ++y_eit;

        } while (check());
      }
    }

    MoFEMFunctionReturnHot(0);
  };

  if (ents_ptr) {
    for (auto p = ents_ptr->const_pair_begin(); p != ents_ptr->const_pair_end();
         ++p) {
      CHKERR loop_ents(p->first, p->second);
    }
  } else {
    // we are looping from the very first possible entity handle (MBVERTEX) to
    // the very last (MBENTITYSET)
    CHKERR loop_ents(get_id_for_min_type<MBVERTEX>(),
                     get_id_for_max_type<MBENTITYSET>());
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldLambda(TwoFieldFunctionOnValues lambda,
                                      const std::string field_name_x,
                                      const std::string field_name_y,
                                      bool error_if_missing,
                                      bool create_if_missing) {
  if (create_if_missing)
    MOFEM_LOG("SELF", Sev::noisy)
        << "Option create_if_missing is set to true but have no effect";
  return fieldLambdaOnValues(lambda, field_name_x, field_name_y,
                             error_if_missing);
}

MoFEMErrorCode FieldBlas::fieldAxpy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing, Range *ents_ptr) {
  MoFEMFunctionBegin;
  auto axpy = [alpha](const double fy, const double fx) {
    return fy + alpha * fx;
  };
  CHKERR fieldLambdaOnValues(axpy, field_name_x, field_name_y, error_if_missing,
                             ents_ptr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldAxpy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing,
                                    bool create_if_missing) {
  return fieldAxpy(alpha, field_name_x, field_name_y, error_if_missing);
}

MoFEMErrorCode FieldBlas::fieldCopy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing, Range *ents_ptr) {
  MoFEMFunctionBegin;
  auto copy = [alpha](const double fy, const double fx) { return alpha * fx; };
  CHKERR fieldLambdaOnValues(copy, field_name_x, field_name_y, error_if_missing,
                             ents_ptr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode FieldBlas::fieldCopy(const double alpha,
                                    const std::string field_name_x,
                                    const std::string field_name_y,
                                    bool error_if_missing,
                                    bool create_if_missing) {
  return fieldCopy(alpha, field_name_x, field_name_y, error_if_missing);
}

MoFEMErrorCode FieldBlas::fieldScale(const double alpha,
                                     const std::string field_name,
                                     Range *ents_ptr) {
  MoFEMFunctionBeginHot;

  auto scale = [alpha](const double v) { return v * alpha; };
  CHKERR fieldLambdaOnValues(scale, field_name, ents_ptr);

  MoFEMFunctionReturnHot(0);
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

} // namespace MoFEM
