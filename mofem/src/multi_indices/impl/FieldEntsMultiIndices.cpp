/** \file FieldEntsMultiIndices.cpp
 * \brief Multi-index containers for entities
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

#define IS_BUILDING_MB
#include <moab/Error.hpp>

namespace MoFEM {

boost::shared_ptr<const FieldTmp<0, 0>> FieldEntityTmp<0, 0>::sFieldPtr;

FieldEntityTmp<0, 0>::FieldEntityTmp(
    const boost::shared_ptr<FieldTmp<0, 0>> field_ptr,
    const boost::shared_ptr<RefEntity> ref_ents_ptr,
    boost::shared_ptr<double *const> field_data_adaptor_ptr,
    boost::shared_ptr<const int> t_max_order_ptr)
    : interface_FieldImpl<FieldTmp<0, 0>, RefEntity>(field_ptr, ref_ents_ptr),
      tagMaxOrderPtr(t_max_order_ptr),
      fieldDataAdaptorPtr(field_data_adaptor_ptr) {

  localUId = getLocalUniqueIdCalculate(field_ptr->getBitNumber(),
                                       ref_ents_ptr->getEnt());

  if (PetscUnlikely(!fieldDataAdaptorPtr))
    THROW_MESSAGE("Pointer to field data adaptor not set");

  if (PetscUnlikely(!tagMaxOrderPtr))
    THROW_MESSAGE("Pointer to max order not set");
}

FieldEntityTmp<-1, -1>::FieldEntityTmp(
    const boost::shared_ptr<FieldTmp<-1, -1>> field_ptr,
    const boost::shared_ptr<RefEntity> ref_ents_ptr,
    boost::shared_ptr<double *const> field_data_adaptor_ptr,
    boost::shared_ptr<const int> t_max_order_ptr)
    : FieldEntityTmp<0, 0>(field_ptr, ref_ents_ptr, field_data_adaptor_ptr,
                           t_max_order_ptr),
      sFieldPtr(field_ptr) {}

boost::shared_ptr<FieldData *const>
FieldEntityTmp<0, 0>::makeSharedFieldDataAdaptorPtr(
    const boost::shared_ptr<Field> &field_ptr,
    const boost::shared_ptr<RefEntity> &ref_ents_ptr) {
  int size;
  FieldData *ptr;
  switch (ref_ents_ptr->getEntType()) {
  case MBVERTEX:
    size = field_ptr->getNbOfCoeffs();
    ptr = static_cast<FieldData *>(
        MoFEM::get_tag_ptr(field_ptr->moab, field_ptr->th_FieldDataVerts,
                           ref_ents_ptr->ent, &size));
    break;
  default:
    ptr = static_cast<FieldData *>(MoFEM::get_tag_ptr(
        field_ptr->moab, field_ptr->th_FieldData, ref_ents_ptr->ent, &size));
  }
  return boost::make_shared<FieldData *const>(ptr);
}

std::ostream &operator<<(std::ostream &os, const FieldEntity &e) {
  os << "ent_global_uid "
     << (UId)e.getGlobalUniqueId()
     // << " ent_local_uid " << (UId)e.get_local_unique_id()
     << " entity " << e.getEnt() << " type " << e.getEntType() << " pstatus "
     << std::bitset<8>(e.getPStatus()) << " owner handle " << e.getOwnerEnt()
     << " owner proc " << e.getOwnerProc() << " order " << e.getMaxOrder()
     << " " << *e.getFieldPtr();
  return os;
}

void FieldEntity_change_order::operator()(FieldEntity *e) {

  moab::Interface &moab = e->sPtr->getBasicDataPtr()->moab;
  const EntityHandle ent = e->getEnt();
  *const_cast<ApproximationOrder *>(e->getMaxOrderPtr()) = order;
  std::size_t nb_dofs = e->getOrderNbDofs(order) * e->getNbOfCoeffs();

  double *tag_field_data = nullptr;
  int tag_field_data_size = 0;

  auto set_verts = [&]() {
    if (e->getFieldPtr()->th_FieldDataVertsType == MB_TAG_SPARSE) {
      // Get pointer and size of field values tag
      rval = moab.tag_get_by_ptr(e->getFieldPtr()->th_FieldDataVerts, &ent, 1,
                                 (const void **)&tag_field_data,
                                 &tag_field_data_size);
      if (nb_dofs) {
        if (nb_dofs != tag_field_data_size) {
          rval = moab.tag_set_data(e->getFieldPtr()->th_FieldDataVerts, &ent, 1,
                                   &*data.begin());
          MOAB_THROW(rval);
        }
      } else if (rval == MB_SUCCESS) {
        rval = moab.tag_delete_data(e->getFieldPtr()->th_FieldDataVerts, &ent, 1);
        MOAB_THROW(rval);
      }
    } else {
      rval = moab.tag_get_by_ptr(e->getFieldPtr()->th_FieldDataVerts, &ent, 1,
                                 (const void **)&tag_field_data);
      MOAB_THROW(rval);
      rval = moab.tag_set_data(e->getFieldPtr()->th_FieldDataVerts, &ent, 1,
                               tag_field_data);
      MOAB_THROW(rval);
    }
  };

  auto set_default = [&]() {
    if (reduceTagSize || nb_dofs) {

      // Get pointer and size of field values tag
      rval = moab.tag_get_by_ptr(e->getFieldPtr()->th_FieldData, &ent, 1,
                                 (const void **)&tag_field_data,
                                 &tag_field_data_size);

      if ((reduceTagSize && nb_dofs != tag_field_data_size) ||
          nb_dofs > tag_field_data_size) {

        // Tag exist and are some data on it
        if (rval == MB_SUCCESS) {
          // Size of tag is different than new size, so copy data to new
          // container
          data.resize(tag_field_data_size);
          FieldData *ptr_begin = static_cast<FieldData *>(tag_field_data);
          FieldData *ptr_end =
              static_cast<FieldData *>(tag_field_data) + tag_field_data_size;
          std::copy(ptr_begin, ptr_end, data.begin());
        }

        if (rval != MB_SUCCESS || nb_dofs) {

          // Set field dof data
          data.resize(nb_dofs, 0);
          int tag_size[1];
          tag_size[0] = data.size();
          void const *tag_data[] = {&data[0]};
          rval = moab.tag_set_by_ptr(e->getFieldPtr()->th_FieldData, &ent, 1,
                                     tag_data, tag_size);
          MOAB_THROW(rval);
          rval = moab.tag_get_by_ptr(e->getFieldPtr()->th_FieldData, &ent, 1,
                                     (const void **)&tag_field_data,
                                     &tag_field_data_size);
          MOAB_THROW(rval);

        } else {

          rval = moab.tag_delete_data(e->getFieldPtr()->th_FieldData, &ent, 1);
          MOAB_THROW(rval);
        }
      }
    }
  };

  switch (e->getEntType()) {
  case MBVERTEX:
    set_verts();
    break;
  default:
    set_default();
  }

  if (nb_dofs)
    const_cast<double *&>(*(e->getEntFieldDataPtr())) = tag_field_data;
  else
    const_cast<double *&>(*(e->getEntFieldDataPtr())) = nullptr;
}

} // namespace MoFEM
