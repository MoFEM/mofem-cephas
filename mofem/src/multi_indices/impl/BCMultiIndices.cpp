/** \file BCMultiIndices.cpp
 * \brief Structures for managing boundary conditions
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

#include <Includes.hpp>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <BCMultiIndices.hpp>

namespace MoFEM {

//moab base meshsets
MoFEMErrorCode CubitMeshSets::getTagsHanlders(Interface &moab) {
  MoFEMFunctionBeginHot;
  ErrorCode rval;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,nsTag); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,ssTag); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle((
    std::string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA").c_str(),nsTag_data
  ); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle(
    (std::string(NEUMANN_SET_TAG_NAME)+"__BC_DATA").c_str(),ssTag_data
  ); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME,bhTag); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle(BLOCK_HEADER,bhTag_header); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  rval = moab.tag_get_handle(BLOCK_ATTRIBUTES,thBlockAttribs); CHKERRQ_MOAB(rval); MOAB_THROW(rval);
  rval = moab.tag_get_handle(NAME_TAG_NAME,entityNameTag); CHKERRQ_MOAB(rval);MOAB_THROW(rval);
  MoFEMFunctionReturnHot(0);
}
CubitMeshSets::CubitMeshSets(Interface &moab,const EntityHandle _meshset):
  meshset(_meshset),
  cubitBcType(UNKNOWNSET),
  msId(NULL),
  tag_bc_data(NULL),
  tag_bc_size(0),
  tag_block_header_data(NULL),
  tag_block_attributes(NULL),
  tag_block_attributes_size(0),
  tagName(NULL),
  meshsets_mask(NODESET|SIDESET|BLOCKSET) {
  
  ierr = getTagsHanlders(moab); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ErrorCode rval;
  rval = moab.tag_get_tags_on_entity(meshset,tag_handles); MOAB_THROW(rval);
  std::vector<Tag>::iterator tit = tag_handles.begin();
  for(;tit!=tag_handles.end();tit++) {
    if(
      *tit == nsTag ||
      *tit == ssTag ||
      *tit == bhTag) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&msId); MOAB_THROW(rval);
    }
    if(*tit == nsTag) {
      if(*msId != -1) {
        cubitBcType = NODESET;
      }
    }
    if(*tit == ssTag) {
      if(*msId != -1) {
        cubitBcType = SIDESET;
      }
    }
    if(*tit == bhTag) {
      if(*msId != -1) {
        cubitBcType = BLOCKSET;
      }
    }
    if((*tit == nsTag_data)||(*tit == ssTag_data)) {
      rval = moab.tag_get_by_ptr(
        *tit,&meshset,1,(const void **)&tag_bc_data,&tag_bc_size
      ); MOAB_THROW(rval);
      
      ierr = getTypeFromBcData(cubitBcType); if(ierr>0) THROW_MESSAGE("unrecognized bc_data type");
    }
    if(*tit == bhTag_header) {
      int tag_length;
      rval = moab.tag_get_length(*tit,tag_length);
      rval = moab.tag_get_by_ptr(
        *tit,&meshset,1,(const void **)&tag_block_header_data
      ); MOAB_THROW(rval);
      if(tag_block_header_data[1]>0) cubitBcType |= MATERIALSET;
    }
    if(*tit == thBlockAttribs) {
      rval = moab.tag_get_by_ptr(
        *tit,&meshset,1,(const void **)&tag_block_attributes,&tag_block_attributes_size
      ); MOAB_THROW(rval);
    }
    if(*tit == entityNameTag) {
      rval = moab.tag_get_by_ptr(entityNameTag,&meshset,1,(const void **)&tagName); MOAB_THROW(rval);
      
      ierr = getTypeFromName(cubitBcType); if(ierr>0) THROW_MESSAGE("unrecognized Cubit name type");
    }
  }

  //If BC set has name, unset UNKNOWNNAME
  if(
    cubitBcType.to_ulong() & (
      DISPLACEMENTSET|
      FORCESET|
      PRESSURESET|
      VELOCITYSET|
      ACCELERATIONSET|
      TEMPERATURESET|
      HEATFLUXSET|
      INTERFACESET
    )
  ) {
    if( (cubitBcType & CubitBCType(UNKNOWNNAME)).any() ) {
      cubitBcType = cubitBcType & (~CubitBCType(UNKNOWNNAME));
    }
  }

}
CubitMeshSets::CubitMeshSets(Interface &moab,const CubitBCType _cubit_bc_type,const int _ms_id):
cubitBcType(_cubit_bc_type),msId(NULL),
tag_bc_data(NULL),
tag_bc_size(0),
tag_block_header_data(NULL),
tag_block_attributes(NULL),
tag_block_attributes_size(0),
tagName(NULL),
meshsets_mask(NODESET|SIDESET|BLOCKSET)  {
  
  ierr = getTagsHanlders(moab); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ErrorCode rval;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); MOAB_THROW(rval);
  switch(cubitBcType.to_ulong()) {
    case NODESET:
    rval = moab.tag_set_data(nsTag,&meshset,1,&_ms_id); MOAB_THROW(rval);
    rval = moab.tag_get_by_ptr(nsTag,&meshset,1,(const void **)&msId); MOAB_THROW(rval);
    break;
    case SIDESET:
    rval = moab.tag_set_data(ssTag,&meshset,1,&_ms_id); MOAB_THROW(rval);
    rval = moab.tag_get_by_ptr(ssTag,&meshset,1,(const void **)&msId); MOAB_THROW(rval);
    break;
    case BLOCKSET:
    rval = moab.tag_set_data(bhTag,&meshset,1,&_ms_id); MOAB_THROW(rval);
    rval = moab.tag_get_by_ptr(bhTag,&meshset,1,(const void **)&msId); MOAB_THROW(rval);
    break;
    default:
    {
      PetscTraceBackErrorHandler(
        PETSC_COMM_WORLD,
        __LINE__,PETSC_FUNCTION_NAME,__FILE__,
        MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,"data inconstancy",PETSC_NULL
      );
      PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,
        __LINE__,PETSC_FUNCTION_NAME,__FILE__,
        MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,"data inconstancy",PETSC_NULL
      );
    }
  }
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByDimension(Interface &moab,const int dimension,Range &entities,const bool recursive)  const {
  MoFEMFunctionBeginHot;
  ErrorCode rval;
  //rval = moab.list_entity(meshset); CHKERRQ_MOAB(rval);
  rval = moab.get_entities_by_dimension(meshset,dimension,entities,recursive);
  if(rval !=  MB_SUCCESS) {
    std::ostringstream ss;
    ss << "bc set " << *this << std::endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByDimension(Interface &moab,Range &entities,const bool recursive)  const {
  MoFEMFunctionBeginHot;
  if((cubitBcType&CubitBCType(BLOCKSET)).any()) {
    if(tag_block_header_data!=NULL) {
      return getMeshsetIdEntitiesByDimension(moab,tag_block_header_data[2],entities,recursive);
    } else {
      SETERRQ(PETSC_COMM_SELF,1,"dimension unknown");
    }
  }
  if((cubitBcType&CubitBCType(NODESET)).any()) {
    return getMeshsetIdEntitiesByDimension(moab,0,entities,recursive);
  }
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode CubitMeshSets::getMeshsetIdEntitiesByType(Interface &moab,const EntityType type,Range &entities,const bool recursive)  const {
  MoFEMFunctionBeginHot;
  ErrorCode rval;
  //rval = moab.list_entity(meshset); CHKERRQ_MOAB(rval);
  rval = moab.get_entities_by_type(meshset,type,entities,recursive);
  if(rval !=  MB_SUCCESS) {
    std::ostringstream ss;
    ss << "bc set " << *this << std::endl;
    PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());
  }
  CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getBcData(std::vector<char>& bc_data) const {
  MoFEMFunctionBeginHot;
  bc_data.resize(tag_bc_size);
  copy(&tag_bc_data[0],&tag_bc_data[tag_bc_size],bc_data.begin());
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getBlockHeaderData(std::vector<unsigned int>& material_data) const {
    MoFEMFunctionBeginHot;
    material_data.resize(3);
    copy(&tag_block_header_data[0],&tag_block_header_data[3],material_data.begin());
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::printBlockHeaderData(std::ostream& os) const {
    MoFEMFunctionBeginHot;
    if(tag_block_header_data!=NULL) {
      std::vector<unsigned int> material_data;
      getBlockHeaderData(material_data);
      os << "block_header_data = ";
      std::vector<unsigned int>::iterator vit = material_data.begin();
      for(;vit!=material_data.end();vit++) {
        os << std::hex << (int)((unsigned int)*vit) << " ";
      }
      os << ": ";
      vit = material_data.begin();
      for(;vit!=material_data.end();vit++) {
        os << *vit;
      }
      os << std::endl;
    } else {
      os << "no block header data" << std::endl;
    }
    MoFEMFunctionReturnHot(0);
}

std::string CubitMeshSets::getName() const {
  if(tagName!=NULL) {
    return std::string(tagName);
  } else {
    return "NoNameSet";
  }
}

MoFEMErrorCode CubitMeshSets::printName(std::ostream& os) const {
    MoFEMFunctionBeginHot;
    std::string name = getName();
    os << std::endl;
    os << "Block name:  " << name << std::endl;
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getTypeFromBcData(const std::vector<char> &bc_data,CubitBCType &type) const {
    MoFEMFunctionBeginHot;

    //See CubitBCType in common.hpp
    if(bc_data.size()==0) {
      MoFEMFunctionReturnHot(0);
    }

    if (strcmp (&bc_data[0],"Displacement") == 0) type |= DISPLACEMENTSET;
    else if (strcmp (&bc_data[0],"Force") == 0) type |= FORCESET;
    else if (strcmp (&bc_data[0],"Velocity") == 0) type |= VELOCITYSET;
    else if (strcmp (&bc_data[0],"Acceleration") == 0) type |= ACCELERATIONSET;
    else if (strcmp (&bc_data[0],"Temperature") == 0) type |= TEMPERATURESET;
    else if (strcmp (&bc_data[0],"Pressure") == 0) type |= PRESSURESET;
    else if (strcmp (&bc_data[0],"HeatFlux") == 0) type |= HEATFLUXSET;
    else if (strcmp (&bc_data[0],"cfd_bc") == 0) type |= INTERFACESET;
    else type |= UNKNOWNNAME;


    MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode CubitMeshSets::getTypeFromBcData(CubitBCType &type) const {
  MoFEMFunctionBeginHot;
  
  std::vector<char> bc_data;
  ierr = getBcData(bc_data); CHKERRG(ierr);
  ierr = getTypeFromBcData(bc_data,type); CHKERRG(ierr);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::printBcData(std::ostream& os) const {
  MoFEMFunctionBeginHot;
  std::vector<char> bc_data;
  getBcData(bc_data);
  os << "bc_data = ";
  std::vector<char>::iterator vit = bc_data.begin();
  for(;vit!=bc_data.end();vit++) {
    os << std::hex << (int)((unsigned char)*vit) << " ";
  }
  os << ": ";
  vit = bc_data.begin();
  for(;vit!=bc_data.end();vit++) {
    os << *vit;
  }
  os << std::endl;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getAttributes(std::vector<double>& attributes) const {
  MoFEMFunctionBeginHot;
  attributes.resize(tag_block_attributes_size);
  if(tag_block_attributes_size>0) {
    copy(&tag_block_attributes[0],&tag_block_attributes[tag_block_attributes_size],attributes.begin());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::setAttributes(moab::Interface &moab,const std::vector<double>& attributes) {
  
  MoFEMFunctionBeginHot;
  int tag_size[] = { (int)attributes.size() };
  void const* tag_data[] = { &*attributes.begin() };
  rval = moab.tag_set_by_ptr(thBlockAttribs,&meshset,1,tag_data,tag_size); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_by_ptr(
    thBlockAttribs,&meshset,1,(const void **)&tag_block_attributes,&tag_block_attributes_size
  ); CHKERRQ_MOAB(rval);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::printAttributes(std::ostream& os) const {
    MoFEMFunctionBeginHot;
    std::vector<double> attributes;
    getAttributes(attributes);
    os << std::endl;
    os << "Block attributes" << std::endl;
    os << "----------------" << std::endl;
    for(unsigned int ii = 0;ii<attributes.size();ii++) {
      os << "attr. no: " << ii+1 << "   value: " << attributes[ii] << std::endl;
    }
    os << std::endl;
    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getTypeFromName(const std::string &name,CubitBCType &type) const {
    MoFEMFunctionBeginHot;

    //See CubitBCType in common.hpp
    if (name.compare(0,11,"MAT_ELASTIC") == 0) {
      type |= MAT_ELASTICSET;
    } else if (name.compare(0,11,"MAT_THERMAL") == 0) {
      type |= MAT_THERMALSET;
    } else if (name.compare(0,12,"MAT_MOISTURE") == 0) {
      type |= MAT_MOISTURESET;
    } else if (name.compare(0,10,"MAT_INTERF") == 0) {
      type |= MAT_INTERFSET;
    } else if (name.compare(0,11,"BODY_FORCES") == 0) {
      type |= BODYFORCESSET;
    }
    //To be extended as appropriate
    else { type |= UNKNOWNNAME; }

    MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode CubitMeshSets::getTypeFromName(CubitBCType &type) const {
    MoFEMFunctionBeginHot;
    
    std::string name = getName();
    ierr = getTypeFromName(name,type); CHKERRG(ierr);
    MoFEMFunctionReturnHot(0);
}

std::ostream& operator<<(std::ostream& os,const CubitMeshSets& e) {

  // get name of cubit meshset
  std::ostringstream ss;
  unsigned jj = 0;
  while(1<<jj != LASTSET_BC) {
    const CubitBCType jj_bc_type = 1<<jj;
    if((e.cubitBcType&jj_bc_type).any()) {
      string bc_type_name;
      ss << " " << string(CubitBCNames[jj+1]);
    }
    ++jj;
  }

  // push data to stream
  os << "meshset " << e.meshset << " type" << ss.str();
  if(e.msId != NULL) os << " msId " << *(e.msId);
  if(e.tagName!=NULL) {
    os << " name " << e.getName();
  }
  if(e.tag_block_header_data != NULL) {
    os << " block header: ";
    os << " blockCol = " << e.tag_block_header_data[0];
    os << " blockMat = " << e.tag_block_header_data[1];
    os << " blockDimension = " << e.tag_block_header_data[2];
  }

  return os;
}

std::ostream& operator<<(std::ostream& os,const DisplacementCubitBcData& e) {
    os << "\n";
    os << "D i s p l a c e m e n t \n \n";
    os << "Flag for X-Translation (0/1): " << (int)e.data.flag1 << "\n";
    os << "Flag for Y-Translation (0/1): " << (int)e.data.flag2 << "\n";
    os << "Flag for Z-Translation (0/1): " << (int)e.data.flag3 << "\n";
    os << "Flag for X-Rotation (0/1): " << (int)e.data.flag4 << "\n";
    os << "Flag for Y-Rotation (0/1): " << (int)e.data.flag5 << "\n";
    os << "Flag for Z-Rotation (0/1): " << (int)e.data.flag6 << "\n \n";

    if (e.data.flag1 == 1)
        os << "Displacement magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Displacement magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Displacement magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Displacement magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Displacement magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Displacement magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Displacement magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Displacement magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Displacement magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Displacement magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Displacement magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
    else os << "Displacement magnitude (Z-Rotation): N/A" << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const ForceCubitBcData& e) {
    os << "\n";
    os << "F o r c e \n \n";
    os << "Force magnitude: " << e.data.value1 << "\n";
    os << "Moment magnitude: " << e.data.value2 << "\n";
    os << "Force direction vector (X-component): " << e.data.value3 << "\n";
    os << "Force direction vector (Y-component): " << e.data.value4 << "\n";
    os << "Force direction vector (Z-component): " << e.data.value5 << "\n";
    os << "Moment direction vector (X-component): " << e.data.value6 << "\n";
    os << "Moment direction vector (Y-component): " << e.data.value7 << "\n";
    os << "Moment direction vector (Z-component): " << e.data.value8 << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const VelocityCubitBcData& e) {
    os << "\n";
    os << "V e l o c i t y \n \n";
    if (e.data.flag1 == 1)
        os << "Velocity magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Velocity magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Velocity magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Velocity magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Velocity magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Velocity magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Velocity magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Velocity magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Velocity magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Velocity magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Velocity magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
    else os << "Velocity magnitude (Z-Rotation): N/A" << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const AccelerationCubitBcData& e) {
    os << "\n";
    os << "A c c e l e r a t i o n \n \n";
    if (e.data.flag1 == 1)
        os << "Acceleration magnitude (X-Translation): " << e.data.value1 << "\n";
    else os << "Acceleration magnitude (X-Translation): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Acceleration magnitude (Y-Translation): " << e.data.value2 << "\n";
    else os << "Acceleration magnitude (Y-Translation): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Acceleration magnitude (Z-Translation): " << e.data.value3 << "\n";
    else os << "Acceleration magnitude (Z-Translation): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Acceleration magnitude (X-Rotation): " << e.data.value4 << "\n";
    else os << "Acceleration magnitude (X-Rotation): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Acceleration magnitude (Y-Rotation): " << e.data.value5 << "\n";
    else os << "Acceleration magnitude (Y-Rotation): N/A" << "\n";
    if (e.data.flag6 == 1)
        os << "Acceleration magnitude (Z-Rotation): " << e.data.value6 << "\n \n";
    else os << "Acceleration magnitude (Z-Rotation): N/A" << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const TemperatureCubitBcData& e) {
    os << "\n";
    os << "T e m p e r a t u r e \n \n";
    if (e.data.flag1 == 1)
        os << "Temperature: " << e.data.value1 << "\n";
    else os << "Temperature (default case): N/A" << "\n";
    if (e.data.flag2 == 1)
        os << "Temperature (thin shell middle): " << e.data.value2 << "\n";
    else os << "Temperature (thin shell middle): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Temperature (thin shell gradient): " << e.data.value3 << "\n";
    else os << "Temperature (thin shell gradient): N/A" << "\n";
    if (e.data.flag4 == 1)
        os << "Temperature (thin shell top): " << e.data.value4 << "\n";
    else os << "Temperature (thin shell top): N/A" << "\n";
    if (e.data.flag5 == 1)
        os << "Temperature (thin shell bottom): " << e.data.value5 << "\n \n";
    else os << "Temperature (thin shell bottom): N/A" << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const PressureCubitBcData& e) {
    os << "\n";
    os << "P r e s s u r e \n \n";
    os << "Pressure flag2: " << (int)e.data.flag2 << "\n";
    os << "Pressure value: " << e.data.value1 << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const HeatFluxCubitBcData& e) {
    os << "\n";
    os << "H e a t  F l u x \n \n";
    if (e.data.flag1 == 1)
        os << "Heat flux value: " << e.data.value1 << "\n";
    else os << "Heat flux is applied on thin shells" << "\n";
    if (e.data.flag2 == 1)
        os << "Heat flux value (thin shell top): " << e.data.value2 << "\n";
    else os << "Heat flux value (thin shell top): N/A" << "\n";
    if (e.data.flag3 == 1)
        os << "Heat flux value (thin shell bottom): " << e.data.value3 << "\n \n";
    else os << "Heat flux value (thin shell bottom): N/A" << "\n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const CfgCubitBcData& e) {
    os << "\n";
    os << "CFD BC \n \n";
    return os;
}

std::ostream& operator<<(std::ostream& os,const BlockSetAttributes& e)
  {
    os << std::endl << "Blcok attributes" << std::endl;
    os << "-------------------" << std::endl;
    os << "User attribute 1 = " << e.data.User1 << std::endl;
    os << "User attribute 2 = " << e.data.User2 << std::endl;
    os << "User attribute 3 = " << e.data.User3 << std::endl;
    os << "User attribute 4 = " << e.data.User4 << std::endl;
    os << "User attribute 5 = " << e.data.User5 << std::endl;
    os << "User attribute 6 = " << e.data.User6 << std::endl;
    os << "User attribute 7 = " << e.data.User7 << std::endl;
    os << "User attribute 8 = " << e.data.User7 << std::endl;
    os << "User attribute 9 = " << e.data.User7 << std::endl;
    os << "User attribute 10 = " << e.data.User10 << std::endl << std::endl;
    return os;
  }

std::ostream& operator<<(std::ostream& os,const Mat_Elastic& e)
    {
        os << std::endl << "Material Properties" << std::endl;
        os << "-------------------" << std::endl;
        os << "Young's modulus  = " << e.data.Young << std::endl;
        os << "Poisson's ratio  = " << e.data.Poisson << std::endl;
        os << "Thermal expansion = " << e.data.ThermalExpansion << std::endl;
        os << "User attribute 1 = " << e.data.User1 << std::endl;
        os << "User attribute 2 = " << e.data.User2 << std::endl;
        os << "User attribute 3 = " << e.data.User3 << std::endl;
        os << "User attribute 4 = " << e.data.User4 << std::endl;
        os << "User attribute 5 = " << e.data.User5 << std::endl;
        os << "User attribute 6 = " << e.data.User6 << std::endl;
        os << "User attribute 7 = " << e.data.User7 << std::endl << std::endl;
        return os;
    }

std::ostream& operator<<(std::ostream& os,const Mat_Elastic_EberleinHolzapfel1& e)
    {
        os << std::endl << "Material Properties" << std::endl;
        os << "-------------------" << std::endl;
        os << "Young's modulus  = " << e.data.Young << std::endl;
        os << "Poisson's ratio  = " << e.data.Poisson << std::endl;
        os << "k1 = " << e.data.k1 << std::endl;
        os << "k2 = " << e.data.k2 << std::endl;
        os << "a0_x = " << e.data.a0x << std::endl;
        os << "a0_y = " << e.data.a0y << std::endl;
        os << "a0_z = " << e.data.a0z << std::endl;
        os << "a1_x = " << e.data.a1x << std::endl;
        os << "a1_y = " << e.data.a1y << std::endl;
        os << "a1_Z = " << e.data.a1z << std::endl << std::endl;
        return os;
    }


std::ostream& operator<<(std::ostream& os,const Mat_Thermal& e)
    {
        os << std::endl << "Material Properties" << std::endl;
        os << "-------------------" << std::endl;
        os << "Conductivity  = " << e.data.Conductivity << std::endl;
        os << "User attribute 1 = " << e.data.HeatCapacity << std::endl;
        os << "User attribute 2 = " << e.data.User2 << std::endl;
        os << "User attribute 3 = " << e.data.User3 << std::endl;
        os << "User attribute 4 = " << e.data.User4 << std::endl;
        os << "User attribute 5 = " << e.data.User5 << std::endl;
        os << "User attribute 6 = " << e.data.User6 << std::endl;
        os << "User attribute 7 = " << e.data.User7 << std::endl;
        os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
        return os;
    }



std::ostream& operator<<(std::ostream& os,const Mat_Moisture& e) {
  os << std::endl << "Material Properties" << std::endl;
  os << "-------------------" << std::endl;
  os << "Diffusivity  = " << e.data.Diffusivity << std::endl;
  os << "Viscosity = " << e.data.Viscosity << std::endl;
  os << "Permeability = " << e.data.Permeability << std::endl;
  os << "User attribute 3 = " << e.data.User3 << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
  return os;
}



std::ostream& operator<<(std::ostream& os,const Block_BodyForces& e) {
  os << std::endl << "Block Body Forces" << std::endl;
  os << "-------------------" << std::endl;
  os << "density  = " << e.data.density << std::endl;
  os << "acceleration_x = " << e.data.acceleration_x << std::endl;
  os << "acceleration_y = " << e.data.acceleration_y << std::endl;
  os << "acceleration_z = " << e.data.acceleration_z << std::endl;
  os << "User attribute 4 = " << e.data.User4 << std::endl;
  os << "User attribute 5 = " << e.data.User5 << std::endl;
  os << "User attribute 6 = " << e.data.User6 << std::endl;
  os << "User attribute 7 = " << e.data.User7 << std::endl;
  os << "User attribute 8 = " << e.data.User8 << std::endl << std::endl;
  return os;
}


std::ostream& operator<<(std::ostream& os,const Mat_Elastic_TransIso& e) {
    os << std::endl << "Material Properties" << std::endl;
    os << "-------------------" << std::endl;
    os << "Young's modulus in xy plane (Ep)     = " << e.data.Youngp << std::endl;
    os << "Young's modulus in z-direction (Ez)  = " << e.data.Youngz << std::endl;
    os << "Poisson's ratio in xy plane (vp)     = " << e.data.Poissonp << std::endl;
    os << "Poisson's ratio in z-direction (vpz) = " << e.data.Poissonpz << std::endl;
    os << "Shear modulus in z-direction (Gzp)   = " << e.data.Shearzp << std::endl << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os,const Mat_Interf& e) {
    os << std::endl << "Material Properties" << std::endl;
    os << "-------------------" << std::endl;
    os << "Elastic module	= " << e.data.alpha << std::endl << std::endl;
    os << "Damage coupling	= " << e.data.beta << std::endl << std::endl;
    os << "Strengh		= " << e.data.ft << std::endl << std::endl;
    os << "Fracture energy	= " << e.data.Gf << std::endl << std::endl;

    return os;
}

void CubitMeshSets_change_name::operator()(CubitMeshSets &e) {
  
  switch (e.cubitBcType.to_ulong()) {
    case BLOCKSET:
    {
      nAme.resize(NAME_TAG_SIZE);
      rval = mOab.tag_set_data(e.entityNameTag,&e.meshset,1,nAme.c_str()); MOAB_THROW(rval);
      rval = mOab.tag_get_by_ptr(e.entityNameTag,&e.meshset,1,(const void **)&e.tagName); MOAB_THROW(rval);
      
      CubitBCType type;
      ierr = e.getTypeFromName(type); if(ierr>0) THROW_MESSAGE("unrecognized Cubit name type");
      e.cubitBcType |= type;
    }; break;
    case NODESET:
    case SIDESET:
    {
      nAme.resize(NAME_TAG_SIZE);
      rval = mOab.tag_set_data(e.entityNameTag,&e.meshset,1,nAme.c_str()); MOAB_THROW(rval);
      rval = mOab.tag_get_by_ptr(e.entityNameTag,&e.meshset,1,(const void **)&e.tagName); MOAB_THROW(rval);
    }; break;
    default:
    THROW_MESSAGE("not implemented for this CubitBC type");
  }
}

void CubitMeshSets_change_add_bit_to_cubit_bc_type::operator()(CubitMeshSets &e) {
  e.cubitBcType |= bIt;
}

void CubitMeshSets_change_attributes::operator()(CubitMeshSets &e) {
  
  ierr = e.setAttributes(mOab,aTtr);
  if(ierr>0) THROW_MESSAGE("Attributes not changed");
}

void CubitMeshSets_change_attributes_data_structure::operator()(CubitMeshSets &e) {
  
  
  // Need to run this to set tag size in number of doubles, don;t know nothing about structure
  int tag_size[] = { (int)(aTtr.getSizeOfData()/sizeof(double)) };
  void const* tag_data[] = { aTtr.getDataPtr() };
  rval = mOab.tag_set_by_ptr(e.thBlockAttribs,&e.meshset,1,tag_data,tag_size); MOAB_THROW(rval);
  rval = mOab.tag_get_by_ptr(
    e.thBlockAttribs,&e.meshset,1,(const void **)&e.tag_block_attributes,&e.tag_block_attributes_size
  ); MOAB_THROW(rval);
  // Here I know about structure
  ierr = e.setAttributeDataStructure(aTtr);
  if(ierr>0) THROW_MESSAGE("Attributes not changed");
}

void CubitMeshSets_change_bc_data_structure::operator()(CubitMeshSets &e) {
  
  
  // Need to run this to set tag size, don;t know nothing about structure
  int tag_size[] = { (int)bcData.getSizeOfData() };
  void const* tag_data[] = { bcData.getDataPtr() };
  if((e.cubitBcType&CubitBCType(NODESET)).any()) {
    rval = mOab.tag_set_by_ptr( e.nsTag_data,&e.meshset,1,tag_data,tag_size
    ); MOAB_THROW(rval);
    rval = mOab.tag_get_by_ptr(
      e.nsTag_data,&e.meshset,1,(const void **)&e.tag_bc_data,&e.tag_bc_size
    ); MOAB_THROW(rval);
  } else if((e.cubitBcType&CubitBCType(SIDESET)).any()) {
    rval = mOab.tag_set_by_ptr(
      e.ssTag_data,&e.meshset,1,tag_data,tag_size
    ); MOAB_THROW(rval);
    rval = mOab.tag_get_by_ptr(
      e.ssTag_data,&e.meshset,1,(const void **)&e.tag_bc_data,&e.tag_bc_size
    ); MOAB_THROW(rval);
  } else {
    THROW_MESSAGE("You have to have NODESET or SIDESET to apply BC data on it");
  }
  // Here I know about structure
  ierr = e.setBcDataStructure(bcData);
  if(ierr>0) THROW_MESSAGE("Attributes not changed");
  // Get Type form BC data
  ierr = e.getTypeFromBcData(e.cubitBcType); if(ierr>0) THROW_MESSAGE("unrecognized bc_data type");
}



}
