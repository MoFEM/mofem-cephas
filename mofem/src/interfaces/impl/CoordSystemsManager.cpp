/** \file CoordSystemsManager.cpp
 * \brief Interface managing coordinate systems set to fields
 *
 */

/*
 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/


/** \file MeshsetsManager.cpp
 * \brief Interface to manage meshsets which carrying information about boundary conditions and material blocks
 *
 */

/**
 * The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
 * It can be freely used for educational and research purposes
 * by other institutions. If you use this softwre pleas cite my work.
 *
 * MoFEM is free software: you can redistribute it and/or modify it under
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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>
#include <fem_tools.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <CoordSystemsManager.hpp>

namespace MoFEM {

  MoFEMErrorCode CoordSystemsManager::query_interface(const MOFEMuuid& uuid, UnknownInterface** iface) const {

    MoFEMFunctionBeginHot;
    *iface = NULL;
    if(uuid == IDD_MOFEMMeshsetsManager) {
      *iface = const_cast<CoordSystemsManager*>(this);
      MoFEMFunctionReturnHot(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    MoFEMFunctionReturnHot(0);
  }

  CoordSystemsManager::CoordSystemsManager(const Core& core):
  cOre(const_cast<Core&>(core)) {
  }

  CoordSystemsManager::~CoordSystemsManager() {}

  MoFEMErrorCode CoordSystemsManager::getTags(int verb) {

    Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    MoFEMFunctionBeginHot;
    const int def_val_len = 0;
    //Coordinate systems
    const int def_coord_sys_dim[] = { 0,0,0,0 };
    rval = moab.tag_get_handle(
      "_CoordSysDim",4,MB_TYPE_INTEGER,th_CoordSysDim,MB_TAG_CREAT|MB_TAG_SPARSE,&def_coord_sys_dim
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_get_handle(
      "_CoordSysName",
      def_val_len,
      MB_TYPE_OPAQUE,
      th_CoordSysName,
      MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
      NULL
    ); CHKERRQ_MOAB(rval);
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode CoordSystemsManager::clearMap() {
    MoFEMFunctionBeginHot;
    coordinateSystems.clear();
    MoFEMFunctionReturnHot(0);
  }


  MoFEMErrorCode CoordSystemsManager::initialiseDatabaseFromMesh(int verb) {

    Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    MoFEMFunctionBeginHot;

    Range meshsets;
    rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERRQ_MOAB(rval);
    Range::iterator mit;
    //loop all meshsehset to find coordinate system
    mit = meshsets.begin();
    for(;mit!=meshsets.end();mit++) {
      try {
        const char *cs_name;
        int cs_name_size;
        rval = moab.tag_get_by_ptr(
          th_CoordSysName,&*mit,1,(const void **)&cs_name,&cs_name_size
        ) ;
        if(rval == MB_SUCCESS && cs_name_size) {
          int dim[4];
          rval = moab.tag_get_data(th_CoordSysDim,&*mit,1,dim); CHKERRQ_MOAB(rval);
          if((dim[0]+dim[1]+dim[2]+dim[3])!=0) {
            std::pair<CoordSys_multiIndex::iterator,bool> p = coordinateSystems.insert(
              boost::make_shared<CoordSys>(moab,*mit)
            );
            if(!p.second) {
              SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"meshset to coord system not inserted");
            }
          }
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
    }
    { // Create cartesian coordinate system if not exist
      CoordSys_multiIndex::index<CoordSysName_mi_tag >::type::iterator csit;
      csit = coordinateSystems.get<CoordSysName_mi_tag >().find("CARTESIAN3D");
      if(csit==coordinateSystems.get<CoordSysName_mi_tag >().end()) {
        EntityHandle meshset;
        rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
        int dim[] = { 3,0,0,0 };
        rval = moab.tag_set_data(th_CoordSysDim,&meshset,1,dim); CHKERRQ_MOAB(rval);
        std::string sys_name_str = "CARTESIAN3D";
        void const* sys_name[] = { sys_name_str.c_str() };
        int sys_name_size[1];
        sys_name_size[0] = sys_name_str.size();
        rval = moab.tag_set_by_ptr(
          th_CoordSysName,&meshset,1,sys_name,sys_name_size
        ); CHKERRQ_MOAB(rval);
        std::pair<CoordSys_multiIndex ::iterator,bool> p = coordinateSystems.insert(
          boost::make_shared<CoordSys>(moab,meshset)
        );
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"MeshSet to coord system not inserted");
        }
      }
      csit = coordinateSystems.get<CoordSysName_mi_tag>().find("UNDEFINED");
      if(csit==coordinateSystems.get<CoordSysName_mi_tag>().end()) {
        EntityHandle meshset;
        rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
        int dim[] = { -1,0,0,0 };
        rval = moab.tag_set_data(th_CoordSysDim,&meshset,1,dim); CHKERRQ_MOAB(rval);
        std::string sys_name_str = "UNDEFINED";
        void const* sys_name[] = { sys_name_str.c_str() };
        int sys_name_size[1];
        sys_name_size[0] = sys_name_str.size();
        rval = moab.tag_set_by_ptr(
          th_CoordSysName,&meshset,1,sys_name,sys_name_size
        ); CHKERRQ_MOAB(rval);
        std::pair<CoordSys_multiIndex ::iterator,bool> p = coordinateSystems.insert(
          boost::make_shared<CoordSys>(moab,meshset)
        );
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"MeshSet to coord system not inserted");
        }
      }
    }
    //PetscSynchronizedFlush(comm,PETSC_STDOUT);
    CoordSys_multiIndex::index<CoordSysName_mi_tag>::type::iterator undefined_cs_it;
    undefined_cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find("UNDEFINED");
    if(undefined_cs_it==coordinateSystems.get<CoordSysName_mi_tag>().end()) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Undefined system not found");
    }

    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode CoordSystemsManager::addCoordinateSystem(const int cs_dim[],const std::string name) {

    Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    EntityHandle meshset;
    MoFEMFunctionBeginHot;
    rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
    rval = moab.tag_set_data(th_CoordSysDim,&meshset,1,cs_dim); CHKERRQ_MOAB(rval);
    void const* sys_name[] = { name.c_str() };
    int sys_name_size[1];
    sys_name_size[0] = name.size();
    rval = moab.tag_set_by_ptr(
      th_CoordSysName,&meshset,1,sys_name,sys_name_size
    ); CHKERRQ_MOAB(rval);
    std::pair<CoordSys_multiIndex ::iterator,bool> p = coordinateSystems.insert(
      boost::make_shared<CoordSys>(moab,meshset)
    );
    if(!p.second) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"MeshSet to coord system not inserted");
    }
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode CoordSystemsManager::setFieldCoordinateSystem(const std::string field_name,const std::string cs_name) {

    Interface &m_field = cOre;
    const Field_multiIndex *fields_ptr;
    MoFEMFunctionBeginHot;
    ierr = m_field.get_fields(&fields_ptr); CHKERRG(ierr);
    Field_multiIndex::index<FieldName_mi_tag>::type::iterator field_it;
    field_it = fields_ptr->get<FieldName_mi_tag>().find(field_name);
    if(field_it==fields_ptr->get<FieldName_mi_tag>().end()) {
      SETERRQ1(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Field < %s > not found",field_name.c_str()
      );
    }
    CoordSys_multiIndex::index<CoordSysName_mi_tag>::type::iterator cs_it;
    cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find(cs_name);
    if(cs_it==coordinateSystems.get<CoordSysName_mi_tag>().end()) {
      SETERRQ1(
        PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Coord system < %s > not found",cs_name.c_str()
      );
    }
    int dim = 1;
    for(int alpha = 0;alpha<4;alpha++) {
      if((*cs_it)->getDim(alpha)>0) {
        dim *= (*cs_it)->getDim(alpha);
      }
    }
    switch((*field_it)->getSpace()) {
      case NOSPACE:
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"No space given");
      case H1:
      if((*field_it)->getNbOfCoeffs()!=dim) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "dimension mismatch of field and coordinate system"
          "cs dim %d field rank %d",
          dim,(*field_it)->getNbOfCoeffs()
        );
      }
      break;
      case HDIV:
      case HCURL:
      if(3*(*field_it)->getNbOfCoeffs()!=dim) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "dimension mismatch of field and coordinate system"
          "cs dim %d field rank %d",
          dim,(*field_it)->getNbOfCoeffs()
        );
      }
      break;
      case L2:
      if((*field_it)->getNbOfCoeffs()!=dim) {
        SETERRQ2(
          PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
          "dimension mismatch of field and coordinate system"
          "cs dim %d field rank %d",
          dim,(*field_it)->getNbOfCoeffs()
        );
      }
      case NOFIELD:
      case LASTSPACE:
      {};
    }
    bool success = const_cast<Field_multiIndex*>(fields_ptr)->modify(
      fields_ptr->project<0>(field_it),FieldChangeCoordinateSystem(*cs_it)
    );
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode CoordSystemsManager::getCoordSysPtr(const EntityHandle id,boost::shared_ptr<CoordSys> &cs_ptr) {
    MoFEMFunctionBeginHot;
    CoordSys_multiIndex::index<Meshset_mi_tag>::type::iterator cs_it;
    cs_it = coordinateSystems.get<Meshset_mi_tag>().find(id);
    if(cs_it==coordinateSystems.get<Meshset_mi_tag>().end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "Unknown Coordinate System ms_id %lu",
        id
      );
    }
    cs_ptr = *cs_it;
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode CoordSystemsManager::getCoordSysPtr(const string name,boost::shared_ptr<CoordSys> &cs_ptr) {
    MoFEMFunctionBeginHot;
    CoordSys_multiIndex::index<CoordSysName_mi_tag>::type::iterator cs_it;
    cs_it = coordinateSystems.get<CoordSysName_mi_tag>().find(name);
    if(cs_it==coordinateSystems.get<CoordSysName_mi_tag>().end()) {
      SETERRQ1(
        PETSC_COMM_SELF,
        MOFEM_DATA_INCONSISTENCY,
        "Unknown Coordinate System <%s>",
        name.c_str()
      );
    }
    cs_ptr = *cs_it;
    MoFEMFunctionReturnHot(0);
  }

}
