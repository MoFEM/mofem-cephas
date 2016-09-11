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
#include <FEMMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <MeshsetsManager.hpp>

namespace MoFEM {

  PetscErrorCode MeshsetsManager::queryInterface(const MOFEMuuid& uuid, UnknownInterface** iface) {

    PetscFunctionBegin;
    *iface = NULL;
    if(uuid == IDD_MOFEMMeshsetsManager) {
      *iface = dynamic_cast<MeshsetsManager*>(this);
      PetscFunctionReturn(0);
    }
    if(uuid == IDD_MOFEMUnknown) {
      *iface = dynamic_cast<UnknownInterface*>(this);
      PetscFunctionReturn(0);
    }
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
    PetscFunctionReturn(0);
  }

  MeshsetsManager::MeshsetsManager(const MoFEM::Core& core):
  cOre(const_cast<MoFEM::Core&>(core)) {
  }

  PetscErrorCode MeshsetsManager::clearMap() {
    PetscFunctionBegin;
    cubitMeshsets.clear();
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::initialiseDatabseInformationFromMesh(int verb) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    Range meshsets;
    rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERRQ_MOAB(rval);
    for(Range::iterator mit = meshsets.begin();mit!=meshsets.end();mit++) {
      try {
        //check if meshset is cubit meshset
        CubitMeshSets base_meshset(moab,*mit);
        if((base_meshset.cubitBcType&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
          std::pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(base_meshset);
          if(!p.second) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"meshset not inserted");
          }
          if(verb > 0) {
            std::ostringstream ss;
            ss << "read cubit " << base_meshset << std::endl;
            //PetscSynchronizedPrintf(comm,ss.str().c_str());
            PetscPrintf(m_field.get_comm(),ss.str().c_str());
          }
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
    }
    PetscFunctionReturn(0);
  }


  PetscErrorCode MeshsetsManager::getTags(int verb) {
    MoABErrorCode rval;
    PetscFunctionBegin;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    int default_val = -1;
    rval = moab.tag_get_handle(
      DIRICHLET_SET_TAG_NAME,1, MB_TYPE_INTEGER,
      nsTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,1, MB_TYPE_INTEGER,
      ssTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERRQ_MOAB(rval);
    const int def_bc_data_len = 0;
    std::string tag_name = std::string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA";
    rval = moab.tag_get_handle(
      tag_name.c_str(),
      def_bc_data_len,
      MB_TYPE_OPAQUE,
      nsTag_data,
      MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,
      NULL
    ); CHKERRQ_MOAB(rval);
    tag_name = std::string(NEUMANN_SET_TAG_NAME)+"__BC_DATA";
    rval = moab.tag_get_handle(
      tag_name.c_str(),
      def_bc_data_len,
      MB_TYPE_OPAQUE,
      ssTag_data,
      MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER,
      bhTag,MB_TAG_SPARSE|MB_TAG_CREAT,&default_val); CHKERRQ_MOAB(rval);
    std::vector<unsigned int> def_uint_zero(3,0);
    rval= moab.tag_get_handle(BLOCK_HEADER,3*sizeof(unsigned int),MB_TYPE_INTEGER,
      bhTag_header,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_uint_zero[0]
    ); CHKERRQ_MOAB(rval);
    Tag block_attribs;
    int def_Block_Attributes_length = 0;
    rval = moab.tag_get_handle(BLOCK_ATTRIBUTES,def_Block_Attributes_length,MB_TYPE_DOUBLE,
      block_attribs,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL
    ); CHKERRQ_MOAB(rval);
    Tag entity_name_tag;
    rval = moab.tag_get_handle(
      NAME_TAG_NAME,NAME_TAG_SIZE,MB_TYPE_OPAQUE,entity_name_tag,MB_TAG_SPARSE|MB_TAG_CREAT
    ); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printDisplacementSet() const {
    PetscErrorCode ierr;
    DisplacementCubitBcData mydata;
    PetscFunctionBegin;
    ierr = printBcSet(mydata,NODESET|mydata.tYpe.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printPressureSet() const {
    PetscErrorCode ierr;
    PressureCubitBcData mydata;
    PetscFunctionBegin;
    ierr = printBcSet(mydata,SIDESET|mydata.tYpe.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printForceSet() const {
    PetscErrorCode ierr;
    ForceCubitBcData mydata;
    PetscFunctionBegin;
    ierr = printBcSet(mydata,NODESET|mydata.tYpe.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printTemperatureSet() const {
    PetscErrorCode ierr;
    TemperatureCubitBcData mydata;
    PetscFunctionBegin;
    ierr = printBcSet(mydata,NODESET|mydata.tYpe.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printHeatFluxSet() const {
    PetscErrorCode ierr;
    HeatFluxCubitBcData mydata;
    PetscFunctionBegin;
    ierr = printBcSet(mydata,SIDESET|mydata.tYpe.to_ulong()); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::printMaterialsSet() const {
    MoABErrorCode rval;
    PetscErrorCode ierr;
    PetscFunctionBegin;
    const MoFEM::Interface& m_field = cOre;
    const moab::Interface& moab = m_field.get_moab();
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_((*this),BLOCKSET|MAT_ELASTICSET,it)) {
      Mat_Elastic data;
      ierr = it->getAttributeDataStructure(data); CHKERRQ(ierr);
      std::ostringstream ss;
      ss << *it << std::endl;
      ss << data;
      Range tets;
      rval = moab.get_entities_by_type(it->meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
      ss << "MAT_ELATIC msId "<< it->getMeshSetId() << " nb. tets " << tets.size() << std::endl;
      ss << std::endl;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|MAT_THERMALSET,it)) {
        Mat_Thermal data;
        ierr = it->getAttributeDataStructure(data); CHKERRQ(ierr);
        std::ostringstream ss;
        ss << *it << std::endl;
        ss << data;
        PetscPrintf(m_field.get_comm(),ss.str().c_str());
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|MAT_MOISTURESET,it)) {
      Mat_Moisture data;
      ierr = it->getAttributeDataStructure(data); CHKERRQ(ierr);
      std::ostringstream ss;
      ss << *it << std::endl;
      ss << data;
      PetscPrintf(m_field.get_comm(),ss.str().c_str());
    }

    PetscFunctionReturn(0);
  }


  bool MeshsetsManager::checkMeshset(const int ms_id,const CubitBCType cubit_bc_type) {
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      return true;
    }
    return false;
  }

  PetscErrorCode MeshsetsManager::addMeshset(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    if(checkMeshset(ms_id,cubit_bc_type)) {
      SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
    }
    try {
      CubitMeshSets cmeshset(moab,cubit_bc_type,ms_id);
      if((cmeshset.cubitBcType&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
        std::pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(cmeshset);
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"meshset not inserted");
        }
        if(name.size()>0) {
          bool success  = cubitMeshsets.modify(p.first,CubitMeshSets_change_name(moab,name));
          if(!success) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"name to cubit meshset can not be set");
          }
        }
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::setAttribites(
    const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name
  ) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    cit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(cit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
    }
    if(name.size()>0) {
      bool success  = cubitMeshsets.modify(cubitMeshsets.project<0>(cit),CubitMeshSets_change_name(moab,name));
      if(!success) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"name to cubit meshset can not be set");
      }
    }
    bool success = cubitMeshsets.modify(
      cubitMeshsets.project<0>(cit),CubitMeshSets_change_attributes(moab,attributes)
    );
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::setAttribitesByDataStructure(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name
  ) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    cit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(cit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
    }
    if(name.size()>0) {
      bool success  = cubitMeshsets.modify(cubitMeshsets.project<0>(cit),CubitMeshSets_change_name(moab,name));
      if(!success) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"name to cubit meshset can not be set");
      }
    }
    bool success = cubitMeshsets.modify(
      cubitMeshsets.project<0>(cit),CubitMeshSets_change_attributes_data_structure(moab,data)
    );
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::setBcData(
    const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
  ) {
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    cit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(cit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
    }
    bool success = cubitMeshsets.modify(
      cubitMeshsets.project<0>(cit),CubitMeshSets_change_bc_data_structure(moab,data)
    );
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::deleteMeshset(const CubitBCType cubit_bc_type,const int ms_id) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(miit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(PETSC_COMM_SELF,1,"such cubit meshset is already there",ms_id);
    }
    EntityHandle meshset = miit->getMeshSet();
    cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().erase(miit);
    rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getCubitMeshsetPtr(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      *cubit_meshset_ptr = &*miit;
    } else {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"msId = %d is not there",ms_id);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getEntitiesByDimension(
    const int msId,const unsigned int cubit_bc_type, const int dimension,Range &entities,const bool recursive
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(msId,cubit_bc_type));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      ierr = miit->getMeshSetIdEntitiesByDimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
    } else {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"msId = %d is not there",msId);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getEntitiesByDimension(
    const int ms_id,const unsigned int cubit_bc_type, Range &entities,const bool recursive
  ) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      ierr = miit->getMeshSetIdEntitiesByDimension(moab,entities,recursive); CHKERRQ(ierr);
    } else {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"ms_id = %d is not there",ms_id);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getMeshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      meshset = miit->meshset;
    } else {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"ms_id = %d is not there",ms_id);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getMeshsetsByType(const unsigned int cubit_bc_type,Range &meshsets) {
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
      miit = cubitMeshsets.get<CubitMeshSets_mi_tag>().lower_bound(cubit_bc_type);
    CubitMeshSet_multiIndex::index<CubitMeshSets_mi_tag>::type::iterator
      hi_miit = cubitMeshsets.get<CubitMeshSets_mi_tag>().upper_bound(cubit_bc_type);
    for(;miit!=hi_miit;miit++) {
      meshsets.insert(miit->meshset);
    }
    PetscFunctionReturn(0);
  }




}
