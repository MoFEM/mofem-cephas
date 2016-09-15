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

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

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
      ss << "MAT_ELATIC msId "<< it->getMeshsetId() << " nb. tets " << tets.size() << std::endl;
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

  PetscErrorCode MeshsetsManager::addEntitiesToMeshset(const CubitBCType cubit_bc_type,const int ms_id,Range &ents) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    cit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
      boost::make_tuple(ms_id,cubit_bc_type.to_ulong())
    );
    if(cit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
    }
    EntityHandle meshset = cit->getMeshset();
    rval = moab.add_entities(meshset,ents); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::addEntitiesToMeshset(const CubitBCType cubit_bc_type,const int ms_id,const EntityHandle *ents,const int nb_ents) {
    MoABErrorCode rval;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
    cit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(
      boost::make_tuple(ms_id,cubit_bc_type.to_ulong())
    );
    if(cit==cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
    }
    EntityHandle meshset = cit->getMeshset();
    rval = moab.add_entities(meshset,ents,nb_ents); CHKERRQ_MOAB(rval);
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
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
    }
    if(name.size()>0) {
      bool success  = cubitMeshsets.modify(cubitMeshsets.project<0>(cit),CubitMeshSets_change_name(moab,name));
      if(!success) {
        SETERRQ(m_field.get_comm(),MOFEM_OPERATION_UNSUCCESSFUL,"name to cubit meshset can not be set");
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
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
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
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
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
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"such cubit meshset is already there",ms_id);
    }
    EntityHandle meshset = miit->getMeshset();
    cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().erase(miit);
    rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getCubitMeshsetPtr(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type.to_ulong()));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      *cubit_meshset_ptr = &*miit;
    } else {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"msId = %d is not there",ms_id);
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
      ierr = miit->getMeshsetIdEntitiesByDimension(moab,dimension,entities,recursive); CHKERRQ(ierr);
    } else {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"msId = %d is not there",msId);
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
      ierr = miit->getMeshsetIdEntitiesByDimension(moab,entities,recursive); CHKERRQ(ierr);
    } else {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"ms_id = %d is not there",ms_id);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::getMeshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
    MoFEM::Interface &m_field = cOre;
    PetscFunctionBegin;
    CubitMeshSet_multiIndex::index<Composite_Cubit_msId_And_MeshSetType_mi_tag>::type::iterator
      miit = cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().find(boost::make_tuple(ms_id,cubit_bc_type));
    if(miit!=cubitMeshsets.get<Composite_Cubit_msId_And_MeshSetType_mi_tag>().end()) {
      meshset = miit->meshset;
    } else {
      SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"ms_id = %d is not there",ms_id);
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

  struct BlockData {

    const CubitMeshSets *cubitMeshsetPtr;

    int iD;
    string addType;
    string nAme;
    CubitBC bcType;

    Mat_Elastic matElastic;
    DisplacementCubitBcData dispBc;
    ForceCubitBcData forceBc;
    PressureCubitBcData pressureBc;
    TemperatureCubitBcData temperatureBc;
    HeatFluxCubitBcData heatFluxBc;


    std::vector<double> aTtr;
    BlockData():
    aTtr(10) {
      strncpy(dispBc.data.name,"Displacement",12);
      strncpy(forceBc.data.name,"Force",5);
      strncpy(pressureBc.data.name,"Pressure",8);
      strncpy(temperatureBc.data.name,"Temperature",11);
      strncpy(heatFluxBc.data.name,"HeatFlux",8);
    }

  };

  PetscErrorCode MeshsetsManager::setMeshsetFromFile(const string file_name) {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscFunctionBegin;
    try {
      std::ifstream ini_file(file_name.c_str(),std::ifstream::in);
      po::variables_map vm;
      po::options_description config_file_options;
      map<int,BlockData> block_lists;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        block_lists[it->getMeshsetId()].cubitMeshsetPtr = &*it;
        std::ostringstream str_add;
        str_add << "block_" << it->getMeshsetId() << ".add";
        std::ostringstream str_id;
        str_id << "block_" << it->getMeshsetId() << ".id";
        std::ostringstream str_name;
        str_name << "block_" << it->getMeshsetId() << ".name";
        config_file_options.add_options()
        (str_add.str().c_str(),po::value<string>(&block_lists[it->getMeshsetId()].addType)->default_value("UNKNOWNSET"),"Add block set")
        (str_id.str().c_str(),po::value<int>(&block_lists[it->getMeshsetId()].iD)->default_value(-1),"Id of meshset")
        (str_name.str().c_str(),po::value<string>(&block_lists[it->getMeshsetId()].nAme)->default_value(""),"Name of the meshset");
        // Block attributes
        for(int ii = 1;ii<=10;ii++) {
          std::ostringstream str_user;
          str_user << "block_" << it->getMeshsetId() << ".user" << ii;
          config_file_options.add_options()
          (str_user.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].aTtr[ii-1])->default_value(0.0),"Add block attribute");
        }
        // Mat elastic
        {
          // double Young; 			///< Young's modulus
          // double Poisson; 		///< Poisson's ratio
          // double ThermalExpansion;	///< Thermal expansion
          std::ostringstream str_young;
          str_young << "block_" << it->getMeshsetId() << ".young";
          std::ostringstream str_poisson;
          str_poisson << "block_" << it->getMeshsetId() << ".poisson";
          std::ostringstream str_thermal_expansion;
          str_thermal_expansion << "block_" << it->getMeshsetId() << ".thermalexpansion";
          config_file_options.add_options()
          (str_young.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].matElastic.data.Young)->default_value(-1),"Young modulus")
          (str_poisson.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].matElastic.data.Poisson)->default_value(-2),"Poisson ratio")
          (str_thermal_expansion.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].matElastic.data.ThermalExpansion)->default_value(-1),"Thermal expansion");
        }
        // Displacement bc
        {
          // char flag1; //< Flag for X-Translation (0: N/A, 1: specified)
          // char flag2; //< Flag for Y-Translation (0: N/A, 1: specified)
          // char flag3; //< Flag for Z-Translation (0: N/A, 1: specified)
          // char flag4; //< Flag for X-Rotation (0: N/A, 1: specified)
          // char flag5; //< Flag for Y-Rotation (0: N/A, 1: specified)
          // char flag6; //< Flag for Z-Rotation (0: N/A, 1: specified)
          // double value1; //< Value for X-Translation
          // double value2; //< Value for Y-Translation
          // double value3; //< Value for Z-Translation
          // double value4; //< Value for X-Rotation
          // double value5; //< Value for Y-Rotation
          // double value6; //< Value for Z-Rotation
          std::ostringstream str_flag1;
          str_flag1 << "block_" << it->getMeshsetId() << ".disp_flag1";
          std::ostringstream str_flag2;
          str_flag2 << "block_" << it->getMeshsetId() << ".disp_flag2";
          std::ostringstream str_flag3;
          str_flag3 << "block_" << it->getMeshsetId() << ".disp_flag3";
          std::ostringstream str_flag4;
          str_flag4 << "block_" << it->getMeshsetId() << ".disp_flag4";
          std::ostringstream str_flag5;
          str_flag5 << "block_" << it->getMeshsetId() << ".disp_flag5";
          std::ostringstream str_flag6;
          str_flag6 << "block_" << it->getMeshsetId() << ".disp_flag6";
          std::ostringstream str_value1;
          str_value1 << "block_" << it->getMeshsetId() << ".disp_ux";
          std::ostringstream str_value2;
          str_value2 << "block_" << it->getMeshsetId() << ".disp_uy";
          std::ostringstream str_value3;
          str_value3 << "block_" << it->getMeshsetId() << ".disp_uz";
          std::ostringstream str_value4;
          str_value4 << "block_" << it->getMeshsetId() << ".disp_rx";
          std::ostringstream str_value5;
          str_value5 << "block_" << it->getMeshsetId() << ".disp_ry";
          std::ostringstream str_value6;
          str_value6 << "block_" << it->getMeshsetId() << ".disp_rz";
          config_file_options.add_options()
          (str_flag1.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag1)->default_value(0),"flag1")
          (str_flag2.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag2)->default_value(0),"flag2")
          (str_flag3.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag3)->default_value(0),"flag3")
          (str_flag4.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag4)->default_value(0),"flag4")
          (str_flag5.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag5)->default_value(0),"flag5")
          (str_flag6.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].dispBc.data.flag6)->default_value(0),"flag6")
          (str_value1.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value1)->default_value(0),"value1")
          (str_value2.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value2)->default_value(0),"value2")
          (str_value3.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value3)->default_value(0),"value3")
          (str_value4.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value4)->default_value(0),"value4")
          (str_value5.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value5)->default_value(0),"value5")
          (str_value6.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].dispBc.data.value6)->default_value(0),"value6");
        }
        // Fore BC data
        {
          // char zero[3]; //< 3 zeros
          // double value1; //< Force magnitude
          // double value2; //< Moment magnitude
          // double value3; //< X-component of force direction vector
          // double value4; //< Y-component of force direction vector
          // double value5; //< Z-component of force direction vector
          // double value6; //< X-component of moment direction vector
          // double value7; //< Y-component of moment direction vector
          // double value8; //< Z-component of moment direction vector
          // char zero2; // 0
          std::ostringstream str_value1;
          str_value1 << "block_" << it->getMeshsetId() << ".force_magnitude";
          std::ostringstream str_value2;
          str_value2 << "block_" << it->getMeshsetId() << ".moment_magnitude";
          std::ostringstream str_value3;
          str_value3 << "block_" << it->getMeshsetId() << ".moment_fx";
          std::ostringstream str_value4;
          str_value4 << "block_" << it->getMeshsetId() << ".moment_fy";
          std::ostringstream str_value5;
          str_value5 << "block_" << it->getMeshsetId() << ".moment_fz";
          std::ostringstream str_value6;
          str_value6 << "block_" << it->getMeshsetId() << ".moment_mx";
          std::ostringstream str_value7;
          str_value7 << "block_" << it->getMeshsetId() << ".moment_my";
          std::ostringstream str_value8;
          str_value8 << "block_" << it->getMeshsetId() << ".moment_mz";
          config_file_options.add_options()
          (str_value1.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value1)->default_value(0),"value1")
          (str_value2.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value2)->default_value(0),"value2")
          (str_value3.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value3)->default_value(0),"value3")
          (str_value4.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value4)->default_value(0),"value4")
          (str_value5.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value5)->default_value(0),"value5")
          (str_value6.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value6)->default_value(0),"value6")
          (str_value7.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value7)->default_value(0),"value7")
          (str_value8.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].forceBc.data.value8)->default_value(0),"value8");
        }
        {
          // char name[11]; //< 11 characters for "Temperature"
          // char pre1; //< This is always zero
          // char pre2; //< 0: temperature is not applied on thin shells (default); 1: temperature is applied on thin shells
          // char flag1; //< 0: N/A, 1: temperature value applied (not on thin shells)
          // char flag2; //< 0: N/A, 1: temperature applied on thin shell middle
          // char flag3; //< 0: N/A, 1: thin shell temperature gradient specified
          // char flag4; //< 0: N/A, 1: top thin shell temperature
          // char flag5; //< 0: N/A, 1: bottom thin shell temperature
          // char flag6; //< This is always zero
          // double value1; //< Temperature (default case - no thin shells)
          // double value2; //< Temperature for middle of thin shells
          // double value3; //< Temperature gradient for thin shells
          // double value4; //< Temperature for top of thin shells
          // double value5; //< Temperature for bottom of thin shells
          // double value6; //< This is always zero, i.e. ignore
          std::ostringstream str_flag1;
          str_flag1 << "block_" << it->getMeshsetId() << ".temperature_flag1";
          std::ostringstream str_value1;
          str_value1 << "block_" << it->getMeshsetId() << ".temperature_t";
          config_file_options.add_options()
          (str_flag1.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].temperatureBc.data.flag1)->default_value(0),"flag1")
          (str_value1.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].temperatureBc.data.value1)->default_value(0),"value1");
          // TODO: Add more cases, see above
        }
        // Sideset
        {
          // char name[8];   //< 8 characters for "Pressure"
          // char zero;      //< This is always zero
          // char flag2;     //< 0: Pressure is interpreted as pure pressure 1: pressure is interpreted as total force
          // double value1;  //< Pressure value
          std::ostringstream str_flag2;
          str_flag2 << "block_" << it->getMeshsetId() << ".pressure_flag2";
          std::ostringstream str_value1;
          str_value1 << "block_" << it->getMeshsetId() << ".pressure_magnitude";
          config_file_options.add_options()
          (str_flag2.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].pressureBc.data.flag2)->default_value(0),"flag2")
          (str_value1.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].pressureBc.data.value1)->default_value(0),"value1");
        }
        {
          // char name[8]; //< 8 characters for "HeatFlux" (no space)
          // char pre1; //< This is always zero
          // char pre2; //< 0: heat flux is not applied on thin shells (default); 1: heat flux is applied on thin shells
          // char flag1; //< 0: N/A, 1: normal heat flux case (i.e. single value, case without thin shells)
          // char flag2; //< 0: N/A, 1: Thin shell top heat flux specified
          // char flag3; //< 0: N/A, 1: Thin shell bottom heat flux specidied
          // double value1; //< Heat flux value for default case (no thin shells)
          // double value2; //< Heat flux (thin shell top)
          // double value3; //< Heat flux (thin shell bottom)
          std::ostringstream str_flag1;
          str_flag1 << "block_" << it->getMeshsetId() << ".heatflux_flag1";
          std::ostringstream str_value1;
          str_value1 << "block_" << it->getMeshsetId() << ".heatflux_magnitude";
          config_file_options.add_options()
          (str_flag1.str().c_str(),po::value<char>(&block_lists[it->getMeshsetId()].heatFluxBc.data.flag1)->default_value(0),"flag1")
          (str_value1.str().c_str(),po::value<double>(&block_lists[it->getMeshsetId()].heatFluxBc.data.value1)->default_value(0),"value1");
        }
      }
      po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
      store(parsed,vm);
      po::notify(vm);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {

        CubitBCType bc_type;
        unsigned jj = 0;
        while(1<<jj != LASTSET_BC) {
          if(string(CubitBCNames[jj+1])==block_lists[it->getMeshsetId()].addType) {
            cerr << CubitBCNames[jj+1] << " ";
            bc_type = 1<<jj;
          }
          ++jj;
        }
        if(bc_type.none()) {
          SETERRQ1(
            m_field.get_comm(),
            MOFEM_DATA_INCONSISTENCY,
            "Unrecognized type %s\n",block_lists[it->getMeshsetId()].addType.c_str()
          );
        }
        std::cerr << bc_type.to_ulong() << " ";
        std::cerr << it->getMeshsetId() << " ";
        std::cerr << block_lists[it->getMeshsetId()].addType << endl;

        if(bc_type.to_ulong()==BLOCKSET) block_lists[it->getMeshsetId()].bcType = BLOCKSET;
        else if(bc_type.to_ulong()==NODESET) block_lists[it->getMeshsetId()].bcType = NODESET;
        else if(bc_type.to_ulong()==SIDESET) block_lists[it->getMeshsetId()].bcType = SIDESET;
        else {
          SETERRQ1(
            m_field.get_comm(),
            MOFEM_DATA_INCONSISTENCY,
            "Not yet implemented type %s\n",
            block_lists[it->getMeshsetId()].addType.c_str()
          );
        }
        if(block_lists[it->getMeshsetId()].iD==-1) {
          SETERRQ1(
            m_field.get_comm(),
            MOFEM_DATA_INCONSISTENCY,
            "Unset iD number %d\n",
            block_lists[it->getMeshsetId()].iD
          );
        }

      }
      std::vector<std::string> additional_parameters;
      additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
      for(std::vector<std::string>::iterator vit = additional_parameters.begin();
      vit!=additional_parameters.end();vit++) {
        ierr = PetscPrintf(m_field.get_comm(),"** WARNING Unrecognized option %s\n",vit->c_str()); CHKERRQ(ierr);
      }
      for(map<int,BlockData>::iterator mit = block_lists.begin();mit!=block_lists.end();mit++) {
        switch(mit->second.bcType) {
          case UNKNOWNSET:
          break;
          case BLOCKSET: {
            if(
              (CubitBCType(mit->second.bcType)&mit->second.cubitMeshsetPtr->getBcType()).any()&&
              mit->second.iD==mit->second.cubitMeshsetPtr->getMeshsetId()
            ) {
              // Meshset is the same, only modification
            } else {
              ierr = addMeshset(mit->second.bcType,mit->second.iD,mit->second.nAme); CHKERRQ(ierr);
              EntityHandle meshset = mit->second.cubitMeshsetPtr->getMeshset();
              ierr = addEntitiesToMeshset(mit->second.bcType,mit->second.iD,&meshset,1); CHKERRQ(ierr);
            }
            //Add attributes
            ierr = setAttribites(mit->second.bcType,mit->second.iD,mit->second.aTtr); CHKERRQ(ierr);
            //Add material elastic data if value are physical (i.e. Young > 0, Poisson in (-1.0.5) and ThermalExpansion>0)
            if(mit->second.matElastic.data.Young!=-1) {
              ierr = setAttribitesByDataStructure(mit->second.bcType,mit->second.iD,mit->second.matElastic); CHKERRQ(ierr);
            }
          }
          break;
          case NODESET: {
            if(
              (CubitBCType(mit->second.bcType)&mit->second.cubitMeshsetPtr->getBcType()).any()&&
              mit->second.iD==mit->second.cubitMeshsetPtr->getMeshsetId()
            ) {
              // Meshset is the same, only modification
            } else {
              ierr = addMeshset(mit->second.bcType,mit->second.iD); CHKERRQ(ierr);
              EntityHandle meshset = mit->second.cubitMeshsetPtr->getMeshset();
              ierr = addEntitiesToMeshset(mit->second.bcType,mit->second.iD,&meshset,1); CHKERRQ(ierr);
            }
            //Add displacement bc
            if(
              mit->second.dispBc.data.flag1||
              mit->second.dispBc.data.flag2||
              mit->second.dispBc.data.flag3||
              mit->second.dispBc.data.flag4||
              mit->second.dispBc.data.flag5||
              mit->second.dispBc.data.flag6
            ) {
              if(mit->second.dispBc.data.flag1=='0') mit->second.dispBc.data.flag1=0;
              if(mit->second.dispBc.data.flag1=='N') mit->second.dispBc.data.flag1=0;
              if(mit->second.dispBc.data.flag1) mit->second.dispBc.data.flag1=1;
              if(mit->second.dispBc.data.flag2=='0') mit->second.dispBc.data.flag2=0;
              if(mit->second.dispBc.data.flag2=='N') mit->second.dispBc.data.flag2=0;
              if(mit->second.dispBc.data.flag2) mit->second.dispBc.data.flag2=1;
              if(mit->second.dispBc.data.flag3=='0') mit->second.dispBc.data.flag3=0;
              if(mit->second.dispBc.data.flag3=='N') mit->second.dispBc.data.flag3=0;
              if(mit->second.dispBc.data.flag3) mit->second.dispBc.data.flag3=1;
              if(mit->second.dispBc.data.flag4=='0') mit->second.dispBc.data.flag4=0;
              if(mit->second.dispBc.data.flag4=='N') mit->second.dispBc.data.flag4=0;
              if(mit->second.dispBc.data.flag4) mit->second.dispBc.data.flag4=1;
              if(mit->second.dispBc.data.flag5=='0') mit->second.dispBc.data.flag5=0;
              if(mit->second.dispBc.data.flag5=='N') mit->second.dispBc.data.flag5=0;
              if(mit->second.dispBc.data.flag5) mit->second.dispBc.data.flag5=1;
              if(mit->second.dispBc.data.flag6=='0') mit->second.dispBc.data.flag6=0;
              if(mit->second.dispBc.data.flag6=='N') mit->second.dispBc.data.flag6=0;
              if(mit->second.dispBc.data.flag6) mit->second.dispBc.data.flag6=1;
              ierr = setBcData(mit->second.bcType,mit->second.iD,mit->second.dispBc); CHKERRQ(ierr);
            }
            if(mit->second.forceBc.data.value1!=0||mit->second.forceBc.data.value2!=0) {
              ierr = setBcData(mit->second.bcType,mit->second.iD,mit->second.forceBc); CHKERRQ(ierr);
            }
            // Add temperature boundary condition
            if(mit->second.temperatureBc.data.flag1) {
              if(mit->second.temperatureBc.data.flag1=='0') mit->second.temperatureBc.data.flag1=0;
              if(mit->second.temperatureBc.data.flag1=='N') mit->second.temperatureBc.data.flag1=0;
              if(mit->second.temperatureBc.data.flag1) mit->second.temperatureBc.data.flag1=1;
              ierr = setBcData(mit->second.bcType,mit->second.iD,mit->second.temperatureBc); CHKERRQ(ierr);
            }
          }
          break;
          case SIDESET: {
            if(
              (CubitBCType(mit->second.bcType)&mit->second.cubitMeshsetPtr->getBcType()).any()&&
              mit->second.iD==mit->second.cubitMeshsetPtr->getMeshsetId()
            ) {
              // Meshset is the same, only modification
            } else {
              ierr = addMeshset(mit->second.bcType,mit->second.iD); CHKERRQ(ierr);
              EntityHandle meshset = mit->second.cubitMeshsetPtr->getMeshset();
              ierr = addEntitiesToMeshset(mit->second.bcType,mit->second.iD,&meshset,1); CHKERRQ(ierr);
            }
            // Add pressure
            if(mit->second.pressureBc.data.value1!=0) {
              if(mit->second.pressureBc.data.flag2=='0') mit->second.pressureBc.data.flag2=0;
              if(mit->second.pressureBc.data.flag2=='N') mit->second.pressureBc.data.flag2=0;
              if(mit->second.pressureBc.data.flag2) mit->second.pressureBc.data.flag2=1;
              ierr = setBcData(mit->second.bcType,mit->second.iD,mit->second.pressureBc); CHKERRQ(ierr);
            }
            // Add heat flux
            if(mit->second.heatFluxBc.data.value1!=0) {
              if(mit->second.heatFluxBc.data.flag1=='0') mit->second.heatFluxBc.data.flag1=0;
              if(mit->second.heatFluxBc.data.flag1=='N') mit->second.heatFluxBc.data.flag1=0;
              if(mit->second.heatFluxBc.data.flag1) mit->second.heatFluxBc.data.flag1=1;
              ierr = setBcData(mit->second.bcType,mit->second.iD,mit->second.heatFluxBc); CHKERRQ(ierr);
            }
          }
          break;
          default:
          SETERRQ(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"Not yet implemented type\n");
        }
      }
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode MeshsetsManager::setMeshsetFromFile() {
    PetscErrorCode ierr;
    MoFEM::Interface &m_field = cOre;
    moab::Interface &moab = m_field.get_moab();
    PetscBool flg_file;
    char meshset_file_name[255];
    PetscFunctionBegin;
    ierr = PetscOptionsBegin(m_field.get_comm(),"","Set meshsets form file","none"); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-meshsets_config",
      "meshsets config  file name","",
      "add_cubit_meshsets.in",
      meshset_file_name,
      255,
      &flg_file
    ); CHKERRQ(ierr);
    if(flg_file==PETSC_TRUE) {
      ifstream f(meshset_file_name);
      if(!f.good()) {
        SETERRQ1(m_field.get_comm(),MOFEM_DATA_INCONSISTENCY,"File configuring meshsets ( %s ) can not be open\n",meshset_file_name);
      }
      ierr = setMeshsetFromFile(string(meshset_file_name)); CHKERRQ(ierr);
    }
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


}
