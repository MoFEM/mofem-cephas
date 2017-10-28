/** \file FieldCore.cpp
 * \brief Core interface methods for managing fields.
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

BitFieldId Core::getBitFieldId(const std::string& name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) {
    THROW_MESSAGE("field < "+name+" > not in database (top tip: check spelling)");
  }
  return (*miit)->getId();
}
std::string Core::getBitFieldIdName(const BitFieldId id) const {
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set.find(id);
  return (*miit)->getName();
}
EntityHandle Core::get_field_meshset(const BitFieldId id) const {
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set.find(id);
  if(miit==set.end()) THROW_MESSAGE("field not in database (top tip: check spelling)");
  return (*miit)->meshSet;
}
EntityHandle Core::get_field_meshset(const std::string& name) const {
  return get_field_meshset(getBitFieldId(name));
}

bool Core::check_field(const std::string &name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) return false;
  return true;
}

const Field* Core::get_field_structure(const std::string& name) {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) {
    throw MoFEMException(
      MOFEM_NOT_FOUND,
      std::string("field < "+name+" > not in databse (top tip: check spelling)").c_str()
    );
  }
  return miit->get();
}

PetscErrorCode Core::get_field_entities_by_dimension(const std::string name,int dim,Range &ents) const {

  MoFEMFunctionBeginHot;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_dimension(meshset,dim,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_field_entities_by_type(const std::string name,EntityType type,Range &ents) const {

  MoFEMFunctionBeginHot;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_type(meshset,type,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_field_entities_by_handle(const std::string name,Range &ents) const {

  MoFEMFunctionBeginHot;
  try {
    EntityHandle meshset = get_field_meshset(name);
    rval = moab.get_entities_by_handle(meshset,ents,true); CHKERRQ_MOAB(rval);
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::add_field(
  const std::string& name,
  const FieldSpace space,
  const FieldApproximationBase base,
  const FieldCoefficientsNumber nb_of_coefficients,
  const TagType tag_type,
  const enum MoFEMTypes bh,
  int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator fit;
  fit = fIelds.get<FieldName_mi_tag>().find(name);
  if(fit != fIelds.get<FieldName_mi_tag>().end() ) {
    if(bh == MF_EXCL) {
      SETERRQ1(cOmm,MOFEM_OPERATION_UNSUCCESSFUL,"field is <%s> in database",name.c_str());
    }
  } else {
    EntityHandle meshset;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
    //id
    BitFieldId id = getFieldShift();
    rval = moab.tag_set_data(th_FieldId,&meshset,1,&id); CHKERRQ_MOAB(rval);
    //space
    rval = moab.tag_set_data(th_FieldSpace,&meshset,1,&space); CHKERRQ_MOAB(rval);
    //base
    rval = moab.tag_set_data(th_FieldBase,&meshset,1,&base); CHKERRQ_MOAB(rval);
    //name
    void const* tag_data[] = { name.c_str() };
    int tag_sizes[1]; tag_sizes[0] = name.size();
    rval = moab.tag_set_by_ptr(th_FieldName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
    //name data prefix
    std::string name_data_prefix("_App_Data");
    void const* tag_prefix_data[] = { name_data_prefix.c_str() };
    int tag_prefix_sizes[1]; tag_prefix_sizes[0] = name_data_prefix.size();
    rval = moab.tag_set_by_ptr(th_FieldName_DataNamePrefix,&meshset,1,tag_prefix_data,tag_prefix_sizes); CHKERRQ_MOAB(rval);
    Tag th_AppOrder,th_FieldData,th_Rank;
    //data
    std::string Tag_data_name = name_data_prefix+name;
    const int def_len = 0;
    rval = moab.tag_get_handle(
      Tag_data_name.c_str(),
      def_len,
      MB_TYPE_OPAQUE,
      th_FieldData,
      MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
      NULL
    ); CHKERRQ_MOAB(rval);
    //order
    ApproximationOrder def_ApproximationOrder = -1;
    std::string Tag_ApproximationOrder_name = "_App_Order_"+name;
    rval = moab.tag_get_handle(
      Tag_ApproximationOrder_name.c_str(),
      sizeof(ApproximationOrder),
      MB_TYPE_OPAQUE,
      th_AppOrder,MB_TAG_CREAT|MB_TAG_BYTES|tag_type,
      &def_ApproximationOrder
    ); CHKERRQ_MOAB(rval);
    //rank
    int def_rank = 1;
    std::string Tag_rank_name = "_Field_Rank_"+name;
    rval = moab.tag_get_handle(
      Tag_rank_name.c_str(),
      sizeof(FieldCoefficientsNumber),
      MB_TYPE_OPAQUE,
      th_Rank,
      MB_TAG_CREAT|MB_TAG_BYTES|tag_type,
      &def_rank
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_set_data(th_Rank,&meshset,1,&nb_of_coefficients); CHKERRQ_MOAB(rval);
    //add meshset
    std::pair<Field_multiIndex::iterator,bool> p;
    try {
      CoordSystemsManager *cs_manger_ptr;
      ierr = getInterface(cs_manger_ptr); CHKERRQ(ierr);
      boost::shared_ptr<CoordSys > undefined_cs_ptr;
      ierr = cs_manger_ptr->getCoordSysPtr("UNDEFINED",undefined_cs_ptr); CHKERRQ(ierr);
      int sys_name_size[1];
      sys_name_size[0] = undefined_cs_ptr->getName().size();
      void const* sys_name[] = { &*undefined_cs_ptr->getNameRef().begin() };
      rval = moab.tag_set_by_ptr(
        cs_manger_ptr->get_th_CoordSysName(),&meshset,1,sys_name,sys_name_size
      ); CHKERRQ_MOAB(rval);
      EntityHandle coord_sys_id = undefined_cs_ptr->getMeshset();
      rval = moab.add_entities(coord_sys_id,&meshset,1); CHKERRQ_MOAB(rval);
      p = fIelds.insert(boost::make_shared<Field>(moab,meshset,undefined_cs_ptr));
      if(bh == MF_EXCL) {
        if(!p.second) SETERRQ1(
          cOmm,MOFEM_NOT_FOUND,
          "field not inserted %s (top tip, it could be already there)",
          Field(moab,meshset,undefined_cs_ptr).getName().c_str()
        );
      }
    } catch (MoFEMException const &e) {
      SETERRQ(cOmm,e.errorCode,e.errorMessage);
    }
    if(verbose > 0) {
      std::ostringstream ss;
      ss << "add: " << **p.first << std::endl;
      PetscPrintf(cOmm,ss.str().c_str());
    }
  }
  //unt
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::addEntsToFieldByDim(
  const Range &ents,const int dim,const std::string& name,int verb
) {

  *buildMoFEM = 0;
  EntityHandle idm = no_handle;
  if(verb==-1) verb = verbose;
  MoFEMFunctionBeginHot;
  try {
    idm = get_field_meshset(name);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  FieldSpace space;
  rval = moab.tag_get_data(th_FieldSpace,&idm,1,&space); CHKERRQ_MOAB(rval);
  Range nodes,faces,edges;
  switch(space) {
    case L2:
    rval = moab.add_entities(idm,ents); CHKERRQ_MOAB(rval);
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << name;
      ss << " nb. add ents " << ents.size();
      ss << std::endl;
      PetscSynchronizedPrintf(cOmm,ss.str().c_str());
    }
    break;
    case H1:
    rval = moab.add_entities(idm,ents); CHKERRQ_MOAB(rval);
    if(dim>0) {
      rval = moab.get_adjacencies(
        ents,0,false,nodes,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      {
        Range topo_nodes;
        rval = moab.get_connectivity(
          ents,topo_nodes,true
        ); CHKERRQ_MOAB(rval);
        Range mid_nodes;
        rval = moab.get_connectivity(
          ents,mid_nodes,false
        ); CHKERRQ_MOAB(rval);
        mid_nodes = subtract(mid_nodes,topo_nodes);
        nodes = subtract(nodes,mid_nodes);
      }
      rval = moab.add_entities(idm,nodes); CHKERRQ_MOAB(rval);
    }
    if(dim>2) {
      rval = moab.get_adjacencies(
        ents,2,false,faces,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(idm,faces); CHKERRQ_MOAB(rval);
    }
    if(dim>1) {
      rval = moab.get_adjacencies(
        ents,1,false,edges,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
    }
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << name;
      ss << " nb. add ents " << ents.size();
      ss << " nb. add faces " << faces.size();
      ss << " nb. add edges " << edges.size();
      ss << " nb. add nodes " << nodes.size();
      ss << std::endl;
      PetscSynchronizedPrintf(cOmm,ss.str().c_str());
    }
    break;
    case HCURL:
    rval = moab.add_entities(idm,ents); CHKERRQ_MOAB(rval);
    if(dim>2) {
      rval = moab.get_adjacencies(ents,2,false,faces,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(idm,faces); CHKERRQ_MOAB(rval);
    }
    if(dim>1) {
      rval = moab.get_adjacencies(
        ents,1,false,edges,moab::Interface::UNION
        ); CHKERRQ_MOAB(rval);
        rval = moab.add_entities(idm,edges); CHKERRQ_MOAB(rval);
      }
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << name;
      ss << " nb. add ents " << ents.size();
      ss << " nb. add faces " << faces.size();
      ss << " nb. add edges " << edges.size();
      ss << std::endl;
      PetscSynchronizedPrintf(cOmm,ss.str().c_str());
    }
    break;
    case HDIV:
    rval = moab.add_entities(idm,ents); CHKERRQ_MOAB(rval);
    if(dim>2) {
      rval = moab.get_adjacencies(
        ents,2,false,faces,moab::Interface::UNION
      ); CHKERRQ_MOAB(rval);
      rval = moab.add_entities(idm,faces); CHKERRQ_MOAB(rval);
    }
    if(verb>1) {
      std::ostringstream ss;
      ss << "add entities to field " << name;
      ss << " nb. add ents " << ents.size();
      ss << " nb. add faces " << faces.size();
      ss << std::endl;
      PetscSynchronizedPrintf(cOmm,ss.str().c_str());
    }
    break;
    default:
    SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"sorry, unknown space added to entity");
  }
  if(verb>1) {
    PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::add_ents_to_field_by_dim(
  const Range &ents,const int dim,const std::string& name,int verb
) {
  Range ents_dim = ents.subset_by_dimension(dim);
  return addEntsToFieldByDim(ents_dim,dim,name,verb);
}

PetscErrorCode Core::add_ents_to_field_by_type(
  const Range &ents,const EntityType type,const std::string& name,int verb
) {

  MoFEMFunctionBeginHot;
  Range ents_type = ents.subset_by_type(type);
  if(!ents_type.empty()) {
    const int dim = moab.dimension_from_handle(ents_type[0]);
    ierr = addEntsToFieldByDim(ents_type,dim,name,verb); CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::add_ents_to_field_by_dim(
  const EntityHandle meshset,const int dim,const std::string& name,const bool recursive,int verb
) {
  MoFEMFunctionBeginHot;
  Range ents;
  rval = moab.get_entities_by_dimension(meshset,dim,ents,recursive); CHKERRQ_MOAB(rval);
  ierr = addEntsToFieldByDim(ents,dim,name,verb); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::add_ents_to_field_by_type(
  const EntityHandle meshset,const EntityType type,const std::string& name,const bool recursive,int verb
) {
  MoFEMFunctionBeginHot;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents,recursive); CHKERRQ_MOAB(rval);
  if(!ents.empty()) {
    const int dim = moab.dimension_from_handle(ents[0]);
    ierr = addEntsToFieldByDim(ents,dim,name,verb); CHKERRQ(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::add_ents_to_field_by_EDGEs(const Range &edges,const std::string& name,int verb) {
  return add_ents_to_field_by_type(edges,MBEDGE,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_EDGEs(const EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBEDGE,name,true,verb);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBTRI,name,true,verb);
}
PetscErrorCode Core::add_ents_to_field_by_TRIs(const Range &tris,const std::string& name,int verb) {
  return add_ents_to_field_by_type(tris,MBTRI,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const Range &nodes,const std::string& name,int verb) {
  return add_ents_to_field_by_type(nodes,MBVERTEX,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_VERTICEs(const EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBVERTEX,name,true,verb);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const Range &tets,const std::string& name,int verb) {
  return add_ents_to_field_by_type(tets,MBTET,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_TETs(const EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBTET,name,true,verb);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(const Range &quads,const std::string& name,int verb) {
  return add_ents_to_field_by_type(quads,MBQUAD,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_QUADs(EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBQUAD,name,true,verb);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(const Range &prisms,const std::string& name,int verb) {
  return add_ents_to_field_by_type(prisms,MBPRISM,name,verb);
}
PetscErrorCode Core::add_ents_to_field_by_PRISMs(EntityHandle meshset,const std::string& name,int verb) {
  return add_ents_to_field_by_type(meshset,MBPRISM,name,true,verb);
}

PetscErrorCode Core::set_field_order(
  const Range &ents,const BitFieldId id,const ApproximationOrder order,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;

  //check field & meshset
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.find(id);
  if(miit==set_id.end()) SETERRQ(cOmm,MOFEM_NOT_FOUND,"no filed found");
  EntityHandle idm;
  try {
   idm = get_field_meshset(id);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }

  //intersection with field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle(idm,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  Range field_ents = intersect(ents,ents_of_id_meshset);
  if(verb>1) {
    PetscSynchronizedPrintf(
      cOmm,"nb. of ents for order change in the field <%s> %d\n",
      miit->get()->getName().c_str(),field_ents.size()
    );
  }

  //ent view by field id (in set all MoabEnts has the same FieldId)
  typedef FieldEntity_multiIndex::index<FieldName_mi_tag>::type EntsByName;
  EntsByName& set = entsFields.get<FieldName_mi_tag>();
  EntsByName::iterator eiit = set.lower_bound(miit->get()->getNameRef());
  FieldEntity_multiIndex_ent_view ents_id_view;
  if(eiit != set.end()) {
    EntsByName::iterator hi_eiit = set.upper_bound(miit->get()->getNameRef());
    std::copy(eiit,hi_eiit,std::back_inserter(ents_id_view));
  }
  if(verb>1) {
    PetscSynchronizedPrintf(
      cOmm,"nb. of ents in the multi index field <%s> %d\n",
      miit->get()->getName().c_str(),ents_id_view.size()
    );
  }


  //loop over ents
  int nb_ents_set_order_up = 0;
  int nb_ents_set_order_down = 0;
  int nb_ents_set_order_new = 0;

  Range new_ents;
  // for(Range::iterator eit = field_ents.begin();eit!=field_ents.end();eit++) {
  for(
    Range::const_pair_iterator pit = field_ents.const_pair_begin();
    pit!=field_ents.const_pair_end();pit++
  ) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;

    // Sanity check
    switch((*miit)->getSpace()) {
      case H1:
      if(moab.type_from_handle(first)==MBVERTEX) {
        if(order!=1) {
          SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,
            "approximation order for H1 space and vertex different than 1 makes not sense"
          );
        }
      }
      break;
      case HCURL:
      if(moab.type_from_handle(first)==MBVERTEX) {
        SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"HDIV space on vertices makes no sense");
      }
      break;
      case HDIV:
      if(moab.type_from_handle(first)==MBVERTEX) {
        SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"HDIV space on vertices makes no sense");
      }
      if(moab.type_from_handle(first)==MBEDGE) {
        SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"HDIV space on edges makes no sense");
      }
      break;
      default:
      break;
    }

    // Entity is in database, change order only if needed
    FieldEntity_multiIndex_ent_view::nth_index<1>::type::iterator vit,hi_vit;
    // vit = ents_id_view.get<1>().lower_bound(*eit);
    vit = ents_id_view.get<1>().lower_bound(first);
    hi_vit = ents_id_view.get<1>().upper_bound(second);
    for(;vit!=hi_vit;vit++,first++) {

      // It is a gap in field entities, those need to be added to the field
      while(first!=vit->get()->getEnt()&&first<=second) {
        new_ents.insert(first);
        first++;
      }

      // if(vit!=ents_id_view.get<1>().end()) {
      //entity is in database and order is changed or reset
      const ApproximationOrder old_approximation_order = (*vit)->getMaxOrder();
      if(old_approximation_order==order) continue;
      FieldEntity_multiIndex::iterator miit = entsFields.get<Unique_mi_tag>().
      find((*vit)->getGlobalUniqueId());

      if((*miit)->getMaxOrder()<order) nb_ents_set_order_up++;
      if((*miit)->getMaxOrder()>order) nb_ents_set_order_down++;

      {

      	// set dofs inactive if order is reduced, and set new order to entity if
      	// order is increased (note that dofs are not build if order is
      	// increased)

        DofEntityByNameAndEnt& dofs_by_name = dofsField.get<Composite_Name_And_Ent_mi_tag>();
        DofEntityByNameAndEnt::iterator dit = dofs_by_name.lower_bound(
          boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt())
        );
        if(dit!=dofs_by_name.end()) {
          DofEntityByNameAndEnt::iterator hi_dit = dofs_by_name.upper_bound(
            boost::make_tuple((*miit)->getNameRef(),(*miit)->getEnt())
          );
          for(;dit!=hi_dit;dit++) {
            if((*dit)->getDofOrder()<=order) continue;
            bool success = dofsField.modify(dofsField.project<0>(dit),DofEntity_active_change(false));
            if(!success) SETERRQ(cOmm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          }
        }

        bool success = entsFields.modify(entsFields.project<0>(miit),FieldEntity_change_order(order));
        if(!success) SETERRQ(cOmm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

      }

    }

    if(first<=second) {
      // This entity is not in databse, added to the vector of entities to which
      // tag with new order have to be set.
      new_ents.insert(first,second);
    }

    // else {
    //   // This entity is not in databse, added to the vector of entities to which
    //   // tag with new order have to be set.
    //   new_ents.insert(*eit);
    // }

  }


  // reserve memory for field  dofs
  boost::shared_ptr<std::vector<FieldEntity> > ents_array =
  boost::make_shared<std::vector<FieldEntity> >(std::vector<FieldEntity>());

  // Add sequence to field data structure. Note that entities are allocated
  // once into vector. This vector is passed into sequence as a weak_ptr.
  // Vector is destroyed at the point last entity inside that vector is
  // destroyed.
  miit->get()->getEntSequenceContainer()->push_back(ents_array);
  ents_array->reserve(new_ents.size());

  FieldEntity_change_order modify_order(order);

  //Range::iterator eit = new_ents.begin();
  //for(unsigned int ee = 0;ee!=new_ents.size();ee++,eit++) {
  for(
    Range::const_pair_iterator pit = new_ents.const_pair_begin();
    pit!=new_ents.const_pair_end();pit++
  ) {
    EntityHandle first = pit->first;
    EntityHandle second = pit->second;

    // get tags on entities
    Range new_pair(first,second);
    std::vector<ApproximationOrder*> tag_data_order(new_pair.size());
    rval = moab.tag_get_by_ptr(
      (*miit)->th_AppOrder,new_pair,(const void **)&tag_data_order[0]
    ); CHKERRQ_MOAB(rval);

    // Entity is not in database and order is changed or reset
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent,hi_miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().lower_bound(first);
    hi_miit_ref_ent = refinedEntities.get<Ent_mi_tag>().upper_bound(second);
    for(int ee = 0;miit_ref_ent!=hi_miit_ref_ent;miit_ref_ent++,first++,ee++) {
      // Set tag value
      *tag_data_order[ee] = order;
      // NOTE: This will work with newer compiler only, use push_back for back compatibility.
      // ents_array->emplace_back(*miit,*miit_ref_ent);
      ents_array->push_back(FieldEntity(*miit,*miit_ref_ent));
      modify_order(&(ents_array->back()));
      nb_ents_set_order_new++;
    }
    for(;first<=second;first++) {
      RefEntity ref_ent(basicEntityDataPtr,first);
      // FIXME: need some consistent policy in that case
      if(ref_ent.getBitRefLevel().none()) continue; // not on any mesh and not in database
      std::cerr << ref_ent << std::endl;
      std::cerr << "bit level " << ref_ent.getBitRefLevel() << std::endl;
      SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"Try to add entities which are not seeded or added to database");
    }

  }

  // Add entities to database
  std::vector<boost::shared_ptr<FieldEntity> > ents_shared_array;
  ents_shared_array.reserve(ents_array->size());
  for(
    std::vector<FieldEntity>::iterator
    vit = ents_array->begin();
    vit!=ents_array->end();vit++
  ) {
    // ents_shared_array.emplace_back(ents_array,&*vit);
    ents_shared_array.push_back(boost::shared_ptr<FieldEntity>(ents_array,&*vit));
  }

  // Add new ents to database
  entsFields.insert(ents_shared_array.begin(),ents_shared_array.end());


  if(verb>1) {
    PetscSynchronizedPrintf(
      cOmm,"nb. of entities in field <%s> for which order was increased %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_up,order
    );
    PetscSynchronizedPrintf(
      cOmm,"nb. of entities in field <%s> for which order was reduced %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_down,order
    );
    PetscSynchronizedPrintf(
      cOmm,"nb. of entities in field <%s> for which order set %d (order %d)\n",
      miit->get()->getName().c_str(),nb_ents_set_order_new,order
    );
    PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  }

  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field_order(
  const EntityHandle meshset,const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  rval = moab.get_entities_by_type(meshset,type,ents); CHKERRQ_MOAB(rval);
  if(verb>1) {
    PetscSynchronizedPrintf(cOmm,"nb. of ents for order change %d\n",ents.size());
  }
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  if(verb>1) {
    PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field_order(
  const EntityHandle meshset,const EntityType type,const std::string& name,const ApproximationOrder order,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try{
    ierr = set_field_order(meshset,type,getBitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field_order(const Range &ents,const std::string& name,const ApproximationOrder order,int verb) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  try{
    ierr = set_field_order(ents,getBitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,
  const EntityType type,const BitFieldId id,const ApproximationOrder order,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = BitRefManager(*this).getEntitiesByTypeAndRefLevel(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,id,order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::set_field_order_by_entity_type_and_bit_ref(
  const BitRefLevel &bit,const BitRefLevel &mask,
  const EntityType type,const std::string& name,const ApproximationOrder order,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  *buildMoFEM = 0;
  Range ents;
  ierr = BitRefManager(*this).getEntitiesByTypeAndRefLevel(bit,mask,type,ents,verb); CHKERRQ(ierr);
  try{
    ierr = set_field_order(ents,getBitFieldId(name),order,verb); CHKERRQ(ierr);
  } catch (MoFEMException const &e) {
    SETERRQ(cOmm,e.errorCode,e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::buildFieldForNoField(
  const BitFieldId id,std::map<EntityType,int> &dof_counter,int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  //field it
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  //find fiels
  FieldSetById::iterator miit = set_id.find(id);
  if(miit == set_id.end()) {
    SETERRQ(cOmm,MOFEM_NOT_FOUND,"field not found");
  }

  //ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle((*miit)->meshSet,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      cOmm,"ents in field %s meshset %d\n",(*miit)->getName().c_str(),ents_of_id_meshset.size()
    );
  }
  for(
    Range::iterator eit = ents_of_id_meshset.begin();
    eit != ents_of_id_meshset.end();eit++
  ) {
    //serch if field meshset is in database
    RefEntity_multiIndex::index<Ent_mi_tag>::type::iterator miit_ref_ent;
    miit_ref_ent = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(miit_ref_ent==refinedEntities.get<Ent_mi_tag>().end()) {
      SETERRQ(
        cOmm,
        MOFEM_DATA_INCONSISTENCY,
        "Entity is not in MoFEM databse, entities in field meshset need to be seeded (i.e. bit ref level add to them)"
      );
    }
    std::pair<FieldEntity_multiIndex::iterator,bool> e_miit;
    try {
      //create database entity
      e_miit = entsFields.insert(boost::make_shared<FieldEntity>(*miit,*miit_ref_ent));
    } catch (MoFEMException const &e) {
      SETERRQ(cOmm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(cOmm,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
    //this is nor real field in space (set order to zero)
    bool success = entsFields.modify(e_miit.first,FieldEntity_change_order(0));
    if(!success) SETERRQ(cOmm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    FieldCoefficientsNumber rank = 0;
    //create dofs on this entity (nb. of dofs is equal to rank)
    for(;rank<(*e_miit.first)->getNbOfCoeffs();rank++) {
      std::pair<DofEntity_multiIndex::iterator,bool> d_miit;
      //check if dof is in darabase
      d_miit.first = dofsField.project<0>(
        dofsField.get<Unique_mi_tag>().find(DofEntity::getGlobalUniqueIdCalculate(rank,*(e_miit.first)))
      );
      //if dof is not in databse
      if(d_miit.first==dofsField.end()) {
        //insert dof
        d_miit = dofsField.insert(boost::make_shared<DofEntity>(*(e_miit.first),0,rank,rank));
        if(d_miit.second) {
          dof_counter[MBENTITYSET]++; // Count entities in the meshset
        }
        bool success = dofsField.modify(d_miit.first,DofEntity_active_change(true));
        if(!success) SETERRQ(cOmm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
      }
      //check consistency
      assert((*d_miit.first)->getEntType()==(*e_miit.first)->getEntType());
      assert((*d_miit.first)->getId()==(*e_miit.first)->getId());
      assert((*d_miit.first)->getMaxOrder()==0);
    }
  }
  if(verb>2) {
    typedef DofEntity_multiIndex::index<FieldName_mi_tag>::type DofsByName;
    DofsByName &set = dofsField.get<FieldName_mi_tag>();
    DofsByName::iterator miit2 = set.lower_bound(miit->get()->getNameRef());
    DofsByName::iterator hi_miit2 = set.upper_bound(miit->get()->getNameRef());
    assert(miit2!=hi_miit2);
    for(;miit2!=hi_miit2;miit2++) {
      std::ostringstream ss;
      ss << *miit2 << std::endl;;
      PetscSynchronizedPrintf(cOmm,ss.str().c_str());
    }
    PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::buildFieldForL2H1HcurlHdiv(
  const BitFieldId id,
  std::map<EntityType,int> &dof_counter,
  std::map<EntityType,int> &inactive_dof_counter,
  int verb
) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;

  // Field by ID
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;

  // Find field
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator field_it = set_id.find(id);
  if(field_it == set_id.end()) {
    SETERRQ(cOmm,MOFEM_NOT_FOUND,"field not found");
  }
  const int rank = field_it->get()->getNbOfCoeffs();

  // Just check if there are any DOFs of given field are in database
  const bool dofs_on_field =
  dofsField.get<FieldName_mi_tag>().find(field_it->get()->getNameRef())!=
  dofsField.get<FieldName_mi_tag>().end();

  // Ents in the field meshset
  Range ents_of_id_meshset;
  rval = moab.get_entities_by_handle((*field_it)->meshSet,ents_of_id_meshset,false); CHKERRQ_MOAB(rval);
  if(verb>5) {
    PetscSynchronizedPrintf(
      cOmm,"ents in field %s meshset %d\n",(*field_it)->getName().c_str(),ents_of_id_meshset.size()
    );
  }

  // View of vertex entities to which dofs in given field are added
  FieldEntity_multiIndex_ent_view ents_view;
  std::vector<boost::shared_ptr<DofEntity> > dofs_shared_array;

  // Loop over all entities and insert dofs by Sequences on edges, faces and volumes
  Range::iterator eit = ents_of_id_meshset.begin();
  for(;eit!=ents_of_id_meshset.end();eit++) {

    // Find mofem entity
    FieldEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator e_miit;
    e_miit = entsFields.get<Composite_Name_And_Ent_mi_tag>().find(
      boost::make_tuple((*field_it)->getNameRef(),*eit)
    );
    if(e_miit == entsFields.get<Composite_Name_And_Ent_mi_tag>().end()) continue;

    // Get shared ptr to entity
    boost::shared_ptr<FieldEntity> field_ent = *e_miit;

    // Current dofs on entity
    const int current_nb_dofs_on_ent =
    !dofs_on_field ? 0 : dofsField.get<Composite_Name_And_Ent_mi_tag>().count(
      boost::make_tuple(field_it->get()->getNameRef(),*eit)
    );

    // Insert DOFs into databse
    const int nb_dofs_on_ent = field_ent->getNbDofsOnEnt();
    const int nb_active_dosf_on_ent = rank*field_ent->getOrderNbDofs(field_ent->getMaxOrder());

    if(
      field_ent->getEntType() == MBVERTEX &&
      current_nb_dofs_on_ent == 0
    ) {

      ents_view.push_back(field_ent);

      // // Only one DOF on entity, simply add it and job done
      // // boost::movelib::unique_ptr<DofEntity> mdof
      // // = boost::movelib::make_unique<DofEntity>(field_ent,0,0,0,true);
      // // dofsField.insert(boost::move(mdof));
      // dofsField.insert(
      //   boost::make_shared<DofEntity>(field_ent,0,0,0,true)
      // );

    } else if(current_nb_dofs_on_ent>nb_active_dosf_on_ent) {

      // This is when order is reduced, or no dofs on entity are deleted,
      // then some DOFs are set to inactive.

      // Get all DOFs which index is bigger than number of active DOFs
      DofEntity_multiIndex::index<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>::type::iterator
      dit,hi_dit;
      dit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().lower_bound(
        boost::make_tuple(
          field_it->get()->getNameRef(),*eit,nb_dofs_on_ent
        )
      );
      hi_dit = dofsField.get<Composite_Name_And_Ent_And_EntDofIdx_mi_tag>().upper_bound(
        boost::make_tuple(
          field_it->get()->getNameRef(),*eit,current_nb_dofs_on_ent
        )
      );

      // Modify DOFs as inactive
      for(;dit!=hi_dit;dit++) {
        bool success = dofsField.modify(
          dofsField.project<0>(dit),DofEntity_active_change(false)
        );
        if(!success) {
          SETERRQ(
            cOmm,
            MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful"
          );
        }
        ++inactive_dof_counter[dit->get()->getEntType()];
      }

    } else {

      // Allocate space for all dofs on this entity
      boost::shared_ptr<std::vector<DofEntity> > dofs_array
      = boost::make_shared<std::vector<DofEntity> >(std::vector<DofEntity>());
      dofs_array->reserve(nb_dofs_on_ent);

      // Set weak pointer on entity to vector/Sequence with allocated dofs on
      // this entity
      field_ent->getDofsSequence() = dofs_array;

      // Temerary vector to store shared pointers
      dofs_shared_array.clear();
      dofs_shared_array.reserve(nb_dofs_on_ent);

      // Create dofs instances and shared pointers
      int DD = 0;
      // Loop orders (loop until max entity order is set)
      for(int oo = 0;oo<=field_ent->getMaxOrder();oo++) {
        // Loop nb. dofs at order oo
        for(int dd = 0;dd<field_ent->getOrderNbDofsDiff(oo);dd++) {
          // Loop rank
          for(int rr = 0;rr<rank;rr++,DD++) {
            // push back dofs instanca
            // dofs_array->emplace_back(
            //   DofEntity(field_ent,oo,rr,DD,true)
            // );
            dofs_array->push_back(
              DofEntity(field_ent,oo,rr,DD,true)
            );
            // Push back shared_ptr for DoFS. Note shared pointer is aliased
            // to vector keeping all DOFs on the entity
            // dofs_shared_array.emplace_back(
            //   boost::shared_ptr<DofEntity>(dofs_array,&dofs_array->back())
            // );
            dofs_shared_array.push_back(
              boost::shared_ptr<DofEntity>(dofs_array,&dofs_array->back())
            );
          }
        }
      }
      if(DD != field_ent->getNbDofsOnEnt()) {
        std::ostringstream ss;
        ss << "rank " << rAnk << " ";
        ss << *field_ent << std::endl;
        SETERRQ3(
          cOmm,MOFEM_DATA_INCONSISTENCY,
          "Expected number of DOFs on entity not equal to number added to database (DD = %d != %d = field_ent->getNbDofsOnEnt())\n"
          "%s",
          DD,field_ent->getNbDofsOnEnt(),ss.str().c_str()
        );
      }

      // Finally add dofs to multi-index
      if(!dofs_shared_array.empty()) {

        // This is in case if some DOFs are already there, then need to
        // replace existing shared_ptr with new one. Then old sequence will be
        // released. If this is not done old and new sequence will exist in
        // memory, that will be waste of space.

        std::vector<boost::shared_ptr<DofEntity> >::iterator vit;
        vit = dofs_shared_array.begin();
        for(int ii = 0;ii!=current_nb_dofs_on_ent;ii++,vit++) {
          DofEntity_multiIndex::iterator d_miit;
          d_miit = dofsField.find(vit->get()->getGlobalUniqueId());
          if(d_miit == dofsField.end()) {
            SETERRQ(cOmm,MOFEM_DATA_INCONSISTENCY,"DOFs is not there, but should be");
          }
          bool success = dofsField.modify(
            d_miit,Dof_shared_ptr_change<DofEntity>(*vit)
          );
          if(!success) {
            SETERRQ(
              cOmm,
              MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful"
            );
          }
        }
        if(vit!=dofs_shared_array.end()) {
          // Those DOFs are added
          dof_counter[vit->get()->getEntType()] += std::distance(vit,dofs_shared_array.end());
          // Finally insert DOFs to database
          dofsField.insert(vit,dofs_shared_array.end());
        }

      }

    }

  }

  // Add vertices DOFs by bulk
  boost::shared_ptr<std::vector<DofEntity> > dofs_array
  = boost::make_shared<std::vector<DofEntity> >(std::vector<DofEntity>());
  // Add Sequence of DOFs to sequence container as weak_ptr
  dofs_array->reserve(rank*ents_view.size());
  // Add Sequence of DOFs to sequence container as weak_ptr
  field_it->get()->getDofSequenceContainer()->push_back(dofs_array);
  dofs_shared_array.clear();
  dofs_shared_array.reserve(dofs_array->size());
  for(
    FieldEntity_multiIndex_ent_view::iterator
    eit = ents_view.begin();eit!=ents_view.end();eit++
  ) {
    for(int r = 0;r!=rank;r++) {
      // Construct DOF
      dofs_array->push_back(DofEntity(*eit,1,r,r,true));
      // Construct aliased shared pointer
      dofs_shared_array.push_back(
        boost::shared_ptr<DofEntity>(dofs_array,&dofs_array->back())
      );
      ++dof_counter[MBVERTEX];
    }
  }
  // Insert into Multi-Index container
  dofsField.insert(dofs_shared_array.begin(),dofs_shared_array.end());

  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::build_fields(int verb) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::map<EntityType,int> dof_counter;
    std::map<EntityType,int> inactive_dof_counter;
    if (verb > 0) {
      PetscSynchronizedPrintf(
        cOmm,"Build Field %s (rank %d)\n",(*miit)->getName().c_str(),rAnk
      );
    }
    switch ((*miit)->getSpace()) {
      case NOFIELD:
      ierr = buildFieldForNoField(
        (*miit)->getId(),dof_counter,verb
      ); CHKERRQ(ierr);
      break;
      case L2:
      case H1:
      case HCURL:
      case HDIV:
      ierr = buildFieldForL2H1HcurlHdiv(
        (*miit)->getId(),
        dof_counter,
        inactive_dof_counter,
        verb
      ); CHKERRQ(ierr);
      break;
      default:
      SETERRQ(cOmm,MOFEM_NOT_IMPLEMENTED,"not implemented");
    }
    if (verb > 0) {
      int nb_added_dofs = 0;
      int nb_inactive_added_dofs = 0;
      for(std::map<EntityType,int>::iterator it = dof_counter.begin();it!=dof_counter.end();it++) {
        switch (it->first) {
          case MBVERTEX:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (vertices) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBEDGE:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (edges) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBTRI:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (triangles) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBQUAD:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (quads) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBTET:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (tets) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBPRISM:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (prisms) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          case MBENTITYSET:
          PetscSynchronizedPrintf(
            cOmm,
            "nb added dofs (meshsets) %d (inactive %d)\n",
            it->second,inactive_dof_counter[it->first]
          );
          break;
          default:
          SETERRQ(cOmm,MOFEM_NOT_IMPLEMENTED,"not implemented");
        }
        nb_added_dofs += it->second;
        nb_inactive_added_dofs += inactive_dof_counter[it->first];
      }
      if(verbose>0) {
        PetscSynchronizedPrintf(
          cOmm,
          "nb added dofs %d (number of inactive dofs %d)\n",
          nb_added_dofs,nb_inactive_added_dofs
        );
      }
    }
  }
  *buildMoFEM = 1<<0;
  if(verb>0) {
    PetscSynchronizedPrintf(cOmm,"Nb. dofs %u\n",dofsField.size());
  }
  PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
  //return 0;
}
PetscErrorCode Core::list_dofs_by_field_name(const std::string &field_name) const {
  MoFEMFunctionBeginHot;
  DofEntityByFieldName::iterator dit,hi_dit;
  dit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_dit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  for(;dit!=hi_dit;dit++) {
    std::ostringstream ss;
    ss << "rank " << rAnk << " ";
    ss << *dit << std::endl;
    PetscSynchronizedPrintf(cOmm,ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::list_fields() const {
  MoFEMFunctionBeginHot;
  typedef Field_multiIndex::index<BitFieldId_mi_tag>::type FieldSetById;
  const FieldSetById &set_id = fIelds.get<BitFieldId_mi_tag>();
  FieldSetById::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(cOmm,ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}


PetscErrorCode Core::list_adjacencies() const {
  MoFEMFunctionBeginHot;
  FieldEntityEntFiniteElementAdjacencyMap_multiIndex::iterator miit = entFEAdjacencies.begin();
  for(;miit!=entFEAdjacencies.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscSynchronizedPrintf(cOmm,ss.str().c_str());
  }
  PetscSynchronizedFlush(cOmm,PETSC_STDOUT);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_problem_finite_elements_entities(const std::string &problem_name,const std::string &fe_name,const EntityHandle meshset) {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(cOmm,1,"no such problem like < %s >",problem_name.c_str());
  NumeredEntFiniteElement_multiIndex &numeredFiniteElements = const_cast<NumeredEntFiniteElement_multiIndex&>(p_miit->numeredFiniteElements);
  NumeredEntFiniteElementbyName::iterator miit = numeredFiniteElements.get<FiniteElement_name_mi_tag>().lower_bound(fe_name);
  for(;miit!=numeredFiniteElements.get<FiniteElement_name_mi_tag>().upper_bound(fe_name);miit++) {
    EntityHandle ent = (*miit)->getEnt();
    rval = moab.add_entities(meshset,&ent,1); CHKERRQ_MOAB(rval);
    int part = (*miit)->getPart();
    rval = moab.tag_set_data(th_Part,&ent,1,&part); CHKERRQ_MOAB(rval);
  }
  MoFEMFunctionReturnHot(0);
}

FieldEntityByFieldName::iterator Core::get_ent_field_by_name_begin(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
}
FieldEntityByFieldName::iterator Core::get_ent_field_by_name_end(const std::string &field_name) const {
  return entsFields.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByFieldName::iterator Core::get_dofs_by_name_begin(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
}
DofEntityByFieldName::iterator Core::get_dofs_by_name_end(const std::string &field_name) const {
  return dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_begin(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_name,ent));
}
DofEntityByNameAndEnt::iterator
Core::get_dofs_by_name_and_ent_end(const std::string &field_name,const EntityHandle ent) const {
  return dofsField.get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_name,ent));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_begin(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().lower_bound(boost::make_tuple(field_name,type));
}
DofEntityByNameAndType::iterator
Core::get_dofs_by_name_and_type_end(const std::string &field_name,const EntityType type) const {
  return dofsField.get<Composite_Name_And_Type_mi_tag>().upper_bound(boost::make_tuple(field_name,type));
}
PetscErrorCode Core::check_number_of_ents_in_ents_field(const std::string& name) const {
  MoFEMFunctionBeginHot;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().find(name);
  if(it == fIelds.get<FieldName_mi_tag>().end()) {
    SETERRQ1(cOmm,1,"field not found < %s >",name.c_str());
  }
  EntityHandle meshset = (*it)->getMeshset();

  int num_entities;

  rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
  if(
    entsFields.get<FieldName_mi_tag>().count((*it)->getName()) > (unsigned int)num_entities
  ) {
    SETERRQ1(cOmm,1,"not equal number of entities in meshset and field multiindex < %s >",name.c_str());
  }
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::check_number_of_ents_in_ents_field() const {
  MoFEMFunctionBeginHot;
  Field_multiIndex::index<FieldName_mi_tag>::type::iterator it = fIelds.get<FieldName_mi_tag>().begin();
  for(;it!=fIelds.get<FieldName_mi_tag>().end();it++) {
    if((*it)->getSpace() == NOFIELD) continue; //FIXME: should be treated properly, not test is just skipped for this NOFIELD space
    EntityHandle meshset = (*it)->getMeshset();

    int num_entities;
    rval = moab.get_number_entities_by_handle(meshset,num_entities); CHKERRQ_MOAB(rval);
    if(entsFields.get<FieldName_mi_tag>().count((*it)->getName()) > (unsigned int)num_entities) {
      SETERRQ1(cOmm,1,"not equal number of entities in meshset and field multiindex < %s >",(*it)->getName().c_str());
    }
  }
  MoFEMFunctionReturnHot(0);
}


}
