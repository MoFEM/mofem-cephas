/** \file CoreDataStructures.cpp
 * \brief Myltindex containes, data structures and other low-level functions 
 * 
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
 *
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


#include <CoreDataStructures.hpp>

/*
   Defines the function where the compiled source is located; used 
   in printing error messages. This is defined here in case the user
   does not declare it.
*/
#ifndef __SDIR__
#define __SDIR__ "unknown file source"
#endif

namespace MoFEM {

//moab base meshsets
PetscErrorCode CubitMeshSets::get_tags_hanlders(Interface &moab) {
  PetscFunctionBegin;
  ErrorCode rval;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,nsTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,ssTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle((string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA").c_str(),nsTag_data); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle((string(NEUMANN_SET_TAG_NAME)+"__BC_DATA").c_str(),ssTag_data); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME,bhTag); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle("BLOCK_HEADER",bhTag_header); CHKERR(rval);CHKERR_THROW(rval);
  rval = moab.tag_get_handle("Block_Attributes",block_attribs); CHKERR(rval); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(NAME_TAG_NAME,entityNameTag); CHKERR(rval);CHKERR_THROW(rval);
  PetscFunctionReturn(0);
}
CubitMeshSets::CubitMeshSets(Interface &moab,const EntityHandle _meshset): 
  meshset(_meshset),CubitBCType(UNKNOWNSET),msId(NULL),tag_bc_data(NULL),tag_bc_size(0),
  tag_block_header_data(NULL),tag_block_attributes(NULL),tag_block_attributes_size(0),tag_name_data(NULL),
  meshsets_mask(NODESET|SIDESET|BLOCKSET) {
  PetscErrorCode ierr;
  ierr = get_tags_hanlders(moab); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ErrorCode rval;
  rval = moab.tag_get_tags_on_entity(meshset,tag_handles); CHKERR(rval);CHKERR_THROW(rval);
  vector<Tag>::iterator tit = tag_handles.begin();
  for(;tit!=tag_handles.end();tit++) {
    if(
      *tit == nsTag ||
      *tit == ssTag ||
      *tit == bhTag) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&msId); CHKERR(rval);CHKERR_THROW(rval);
    }
    if(
      (*tit == nsTag_data)||
      (*tit == ssTag_data)) {
    }
    if(*tit == nsTag) {
      if(*msId != -1) {
	CubitBCType = NODESET;
      }
    }
    if(*tit == ssTag) {
      if(*msId != -1) {
	CubitBCType = SIDESET;
      }
    }
    if(*tit == bhTag) {
      if(*msId != -1) {
	CubitBCType = BLOCKSET;
      }
    }
    if(
      (*tit == nsTag_data) ||
      (*tit == ssTag_data)) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&tag_bc_data,&tag_bc_size); CHKERR(rval);CHKERR_THROW(rval);
      PetscErrorCode ierr;
      ierr = get_type_from_bc_data(CubitBCType); if(ierr>0) throw("unrecognised bc_data type");
    }
    if(*tit == bhTag_header) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&tag_block_header_data); CHKERR(rval);CHKERR_THROW(rval);
      if(tag_block_header_data[9]>0) CubitBCType |= MATERIALSET;
    }
    if(*tit == block_attribs) {
      rval = moab.tag_get_by_ptr(*tit,&meshset,1,(const void **)&tag_block_attributes,&tag_block_attributes_size); CHKERR(rval); CHKERR_THROW(rval);
      //for(int ii = 0;ii<tag_block_attributes_size;ii++) {
	//cerr << "RRRRR " << tag_block_attributes[ii] << endl;
      //}
    }
    if(*tit == entityNameTag) {
      rval = moab.tag_get_by_ptr(entityNameTag,&meshset,1,(const void **)&tag_name_data); CHKERR(rval); CHKERR_THROW(rval);
      PetscErrorCode ierr;
      ierr = get_type_from_Cubit_name(CubitBCType); if(ierr>0) throw("unrecognised Cubit name type");
    }
  }

  //If BC set has name, unset UNKNOWNCUBITNAME
  if(CubitBCType.to_ulong() & (	  
      DISPLACEMENTSET|
      FORCESET|
      PRESSURESET|
      VELOCITYSET|
      ACCELERATIONSET|
      TEMPERATURESET|
      HEATFLUXSET|
      INTERFACESET) ) {
     
      if( (CubitBCType & CubitBC_BitSet(UNKNOWNCUBITNAME)).any() ) {
	CubitBCType = CubitBCType & (~CubitBC_BitSet(UNKNOWNCUBITNAME));  
      }
  }

}
CubitMeshSets::CubitMeshSets(Interface &moab,const CubitBC_BitSet _CubitBCType,const int _msId): 
  CubitBCType(_CubitBCType),msId(NULL),
  tag_bc_data(NULL),tag_bc_size(0),
  tag_block_header_data(NULL),
  tag_block_attributes(NULL),
  tag_block_attributes_size(0),
  tag_name_data(NULL),
  meshsets_mask(NODESET|SIDESET|BLOCKSET)  {
  PetscErrorCode ierr;
  ierr = get_tags_hanlders(moab); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ErrorCode rval;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_THROW(rval);

  switch(_CubitBCType.to_ulong()) {
    case NODESET:
      rval = moab.tag_set_data(nsTag,&meshset,1,&_msId); CHKERR(rval);CHKERR_THROW(rval);
      rval = moab.tag_get_by_ptr(nsTag,&meshset,1,(const void **)&msId); CHKERR(rval);CHKERR_THROW(rval);
    break;
    case SIDESET:
      rval = moab.tag_set_data(ssTag,&meshset,1,&_msId); CHKERR(rval);CHKERR_THROW(rval);
      rval = moab.tag_get_by_ptr(ssTag,&meshset,1,(const void **)&msId); CHKERR(rval);CHKERR_THROW(rval);
    break;
    case BLOCKSET:
      rval = moab.tag_set_data(bhTag,&meshset,1,&_msId); CHKERR(rval);CHKERR_THROW(rval);
      rval = moab.tag_get_by_ptr(bhTag,&meshset,1,(const void **)&msId); CHKERR(rval);CHKERR_THROW(rval);
    break;
    default: {
      PetscTraceBackErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	"not implemented yet",PETSC_NULL);
      PetscMPIAbortErrorHandler(PETSC_COMM_WORLD,__LINE__,PETSC_FUNCTION_NAME,__FILE__,__SDIR__,1,PETSC_ERROR_INITIAL,
	"not implemented yet",PETSC_NULL);

    }
  }


}
PetscErrorCode CubitMeshSets::get_Cubit_msId_entities_by_dimension(Interface &moab,const int dimension,Range &entities,const bool recursive)  const {
  PetscFunctionBegin;
  ErrorCode rval;
  rval = moab.get_entities_by_dimension(meshset,dimension,entities,recursive); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::get_Cubit_msId_entities_by_dimension(Interface &moab,Range &entities,const bool recursive)  const {
  PetscFunctionBegin;
  if((CubitBCType&CubitBC_BitSet(BLOCKSET)).any()) {
    if(tag_block_header_data!=NULL) {
      return get_Cubit_msId_entities_by_dimension(moab,tag_block_header_data[11],entities,recursive);
    } else {
      SETERRQ(PETSC_COMM_SELF,1,"dimension unknown");
    }
  }
  if((CubitBCType&CubitBC_BitSet(NODESET)).any()) {
    return get_Cubit_msId_entities_by_dimension(moab,0,entities,recursive);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode CubitMeshSets::get_Cubit_bc_data(vector<char>& bc_data) const {
  PetscFunctionBegin;
  bc_data.resize(tag_bc_size);
  copy(&tag_bc_data[0],&tag_bc_data[tag_bc_size],bc_data.begin());
  PetscFunctionReturn(0);
}

PetscErrorCode CubitMeshSets::get_Cubit_block_header_data(vector<unsigned int>& material_data) const {
    PetscFunctionBegin;
    copy(&tag_block_header_data[0],&tag_block_header_data[12],material_data.begin());
    PetscFunctionReturn(0);
}

PetscErrorCode CubitMeshSets::print_Cubit_block_header_data(ostream& os) const {
    PetscFunctionBegin;
    vector<unsigned int> material_data;
    get_Cubit_block_header_data(material_data);
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
    PetscFunctionReturn(0);
}

string CubitMeshSets::get_Cubit_name() const {
  if(tag_name_data!=NULL) {
    return string(tag_name_data);
  } else {
    return "NoNameSet";
  }
}
    
PetscErrorCode CubitMeshSets::print_Cubit_name(ostream& os) const {
    PetscFunctionBegin;
    string name = get_Cubit_name();
    os << endl;
    os << "Block name:  " << name << endl;
    PetscFunctionReturn(0);
}
           
PetscErrorCode CubitMeshSets::get_type_from_bc_data(const vector<char> &bc_data,CubitBC_BitSet &type) const {
    PetscFunctionBegin;
    
    //See CubitBC_BitSet in common.hpp
    if(bc_data.size()==0) {
      PetscFunctionReturn(0);
    }
    
    if (strcmp (&bc_data[0],"Displacement") == 0)
        type |= DISPLACEMENTSET;
    else if (strcmp (&bc_data[0],"Force") == 0)
        type |= FORCESET;
    else if (strcmp (&bc_data[0],"Velocity") == 0)
        type |= VELOCITYSET;
    else if (strcmp (&bc_data[0],"Acceleration") == 0)
        type |= ACCELERATIONSET;
    else if (strcmp (&bc_data[0],"Temperature") == 0)
        type |= TEMPERATURESET;
    else if (strcmp (&bc_data[0],"Pressure") == 0)
        type |= PRESSURESET;
    else if (strcmp (&bc_data[0],"HeatFlux") == 0)
        type |= HEATFLUXSET;
    else if (strcmp (&bc_data[0],"cfd_bc") == 0)
        type |= INTERFACESET;
    else SETERRQ(PETSC_COMM_SELF,1,"this bc_data is unknown");
    
    PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::get_type_from_bc_data(CubitBC_BitSet &type) const {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  vector<char> bc_data;
  ierr = get_Cubit_bc_data(bc_data); CHKERRQ(ierr);
  ierr = get_type_from_bc_data(bc_data,type); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::print_Cubit_bc_data(ostream& os) const {
  PetscFunctionBegin;
  vector<char> bc_data;
  get_Cubit_bc_data(bc_data);
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
  PetscFunctionReturn(0);
}
PetscErrorCode CubitMeshSets::get_Cubit_attributes(vector<double>& attributes) const {
  PetscFunctionBegin;
  attributes.resize(tag_block_attributes_size);
  if(tag_block_attributes_size>0) {
    copy(&tag_block_attributes[0],&tag_block_attributes[tag_block_attributes_size],attributes.begin());
  }
  PetscFunctionReturn(0);
}
    
PetscErrorCode CubitMeshSets::print_Cubit_attributes(ostream& os) const {
    PetscFunctionBegin;
    vector<double> attributes;
    get_Cubit_attributes(attributes);
    os << endl;
    os << "Block attributes" << endl;
    os << "----------------" << endl;
    for(unsigned int ii = 0;ii<attributes.size();ii++)
        {
            os << "attr. no: " << ii+1 << "   value: " << attributes[ii] << endl;
        }
    os << endl;
    PetscFunctionReturn(0);
}

PetscErrorCode CubitMeshSets::get_type_from_Cubit_name(const string &name,CubitBC_BitSet &type) const {
    PetscFunctionBegin;

    //See CubitBC_BitSet in common.hpp
    if (name.compare(0,11,"MAT_ELASTIC") == 0) {
        type |= MAT_ELASTICSET; }
    else if (name.compare(0,11,"MAT_THERMAL") == 0) {
        type |= MAT_THERMALSET; }
    else if (name.compare(0,12,"MAT_MOISTURE") == 0) {
      type |= MAT_MOISTURESET; }
    else if (name.compare(0,10,"MAT_INTERF") == 0) {
        type |= MAT_INTERFSET; }
    else if (name.compare(0,11,"BODY_FORCES") == 0) {
	type |= BLOCK_BODYFORCESSET; }
    
        //To be extended as appropriate
    
    else { type |= UNKNOWNCUBITNAME; }

    PetscFunctionReturn(0);
}

PetscErrorCode CubitMeshSets::get_type_from_Cubit_name(CubitBC_BitSet &type) const {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    string name = get_Cubit_name();
    ierr = get_type_from_Cubit_name(name,type); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
    
ostream& operator<<(ostream& os,const CubitMeshSets& e) {
  os << "meshset " << e.meshset << " type " << e.CubitBCType;
  if(e.msId != NULL) os << " msId " << *(e.msId);
  if(e.tag_name_data!=NULL) {
    os << " name " << e.get_Cubit_name();
  }
  if(e.tag_block_header_data != NULL) {
    os << " block header: ";
    os << " blockID = " << e.tag_block_header_data[0];
    os << " blockElemType = " << e.tag_block_header_data[1];
    os << " blockMat = " << e.tag_block_header_data[9];
    os << " blockAttributeOrder = " << e.tag_block_header_data[5];
    os << " blockDimension = " << e.tag_block_header_data[11];
  }
  return os;
}

ostream& operator<<(ostream& os,const DisplacementCubitBcData& e) {
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

ostream& operator<<(ostream& os,const ForceCubitBcData& e) {
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

ostream& operator<<(ostream& os,const VelocityCubitBcData& e) {
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
 
ostream& operator<<(ostream& os,const AccelerationCubitBcData& e) {
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

ostream& operator<<(ostream& os,const TemperatureCubitBcData& e) {
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

ostream& operator<<(ostream& os,const PressureCubitBcData& e) {
    os << "\n";
    os << "P r e s s u r e \n \n";
    os << "Pressure value: " << e.data.value1 << "\n \n";
    return os;
}

ostream& operator<<(ostream& os,const HeatfluxCubitBcData& e) {
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

ostream& operator<<(ostream& os,const CfgCubitBcData& e) {
    os << "\n";
    os << "CFD BC \n \n";
    return os;   
}
 
ostream& operator<<(ostream& os,const BlockSetAttributes& e)
  {
    os << endl << "Blcok attributes" << endl;
    os << "-------------------" << endl;
    os << "User attribute 1 = " << e.data.User1 << endl;
    os << "User attribute 2 = " << e.data.User2 << endl;
    os << "User attribute 3 = " << e.data.User3 << endl;
    os << "User attribute 4 = " << e.data.User4 << endl;
    os << "User attribute 5 = " << e.data.User5 << endl;
    os << "User attribute 6 = " << e.data.User6 << endl;
    os << "User attribute 7 = " << e.data.User7 << endl;
    os << "User attribute 8 = " << e.data.User7 << endl;
    os << "User attribute 9 = " << e.data.User7 << endl;
    os << "User attribute 10 = " << e.data.User10 << endl << endl;
    return os;
  }
       
ostream& operator<<(ostream& os,const Mat_Elastic& e)
    {
        os << endl << "Material Properties" << endl;
        os << "-------------------" << endl;
        os << "Young's modulus  = " << e.data.Young << endl;
        os << "Poisson's ratio  = " << e.data.Poisson << endl;
        os << "Thermal expansion = " << e.data.ThermalExpansion << endl;
        os << "User attribute 1 = " << e.data.User1 << endl;
        os << "User attribute 2 = " << e.data.User2 << endl;
        os << "User attribute 3 = " << e.data.User3 << endl;
        os << "User attribute 4 = " << e.data.User4 << endl;
        os << "User attribute 5 = " << e.data.User5 << endl;
        os << "User attribute 6 = " << e.data.User6 << endl;
        os << "User attribute 7 = " << e.data.User7 << endl << endl;
        return os;
    }

ostream& operator<<(ostream& os,const Mat_Elastic_EberleinHolzapfel1& e)
    {
        os << endl << "Material Properties" << endl;
        os << "-------------------" << endl;
        os << "Young's modulus  = " << e.data.Young << endl;
        os << "Poisson's ratio  = " << e.data.Poisson << endl;
        os << "k1 = " << e.data.k1 << endl;
        os << "k2 = " << e.data.k2 << endl;
        os << "a0_x = " << e.data.a0x << endl;
        os << "a0_y = " << e.data.a0y << endl;
        os << "a0_z = " << e.data.a0z << endl;
        os << "a1_x = " << e.data.a1x << endl;
        os << "a1_y = " << e.data.a1y << endl;
        os << "a1_Z = " << e.data.a1z << endl << endl;
        return os;
    }

    
ostream& operator<<(ostream& os,const Mat_Thermal& e)
    {
        os << endl << "Material Properties" << endl;
        os << "-------------------" << endl;
        os << "Conductivity  = " << e.data.Conductivity << endl;
        os << "User attribute 1 = " << e.data.HeatCapacity << endl;
        os << "User attribute 2 = " << e.data.User2 << endl;
        os << "User attribute 3 = " << e.data.User3 << endl;
        os << "User attribute 4 = " << e.data.User4 << endl;
        os << "User attribute 5 = " << e.data.User5 << endl;
        os << "User attribute 6 = " << e.data.User6 << endl;
        os << "User attribute 7 = " << e.data.User7 << endl;
        os << "User attribute 8 = " << e.data.User8 << endl << endl;
        return os;
    }
  
  
  
  ostream& operator<<(ostream& os,const Mat_Moisture& e)
  {
    os << endl << "Material Properties" << endl;
    os << "-------------------" << endl;
    os << "Diffusivity  = " << e.data.Diffusivity << endl;
    os << "User attribute 1 = " << e.data.User1 << endl;
    os << "User attribute 2 = " << e.data.User2 << endl;
    os << "User attribute 3 = " << e.data.User3 << endl;
    os << "User attribute 4 = " << e.data.User4 << endl;
    os << "User attribute 5 = " << e.data.User5 << endl;
    os << "User attribute 6 = " << e.data.User6 << endl;
    os << "User attribute 7 = " << e.data.User7 << endl;
    os << "User attribute 8 = " << e.data.User8 << endl << endl;
    return os;
  }


  
ostream& operator<<(ostream& os,const Block_BodyForces& e)
    {
        os << endl << "Block Body Forces" << endl;
        os << "-------------------" << endl;
        os << "density  = " << e.data.density << endl;
        os << "acceleration_x = " << e.data.acceleration_x << endl;
        os << "acceleration_y = " << e.data.acceleration_y << endl;
        os << "acceleration_z = " << e.data.acceleration_z << endl;
        os << "User attribute 4 = " << e.data.User4 << endl;
        os << "User attribute 5 = " << e.data.User5 << endl;
        os << "User attribute 6 = " << e.data.User6 << endl;
        os << "User attribute 7 = " << e.data.User7 << endl;
        os << "User attribute 8 = " << e.data.User8 << endl << endl;
        return os;
    }

    
ostream& operator<<(ostream& os,const Mat_Elastic_TransIso& e)
{
    os << endl << "Material Properties" << endl;
    os << "-------------------" << endl;
    os << "Young's modulus in xy plane (Ep)     = " << e.data.Youngp << endl;
    os << "Young's modulus in z-direction (Ez)  = " << e.data.Youngz << endl;
    os << "Poisson's ratio in xy plane (vp)     = " << e.data.Poissonp << endl;
    os << "Poisson's ratio in z-direction (vpz) = " << e.data.Poissonpz << endl;
    os << "Shear modulus in z-direction (Gzp)   = " << e.data.Shearzp << endl << endl;
    return os;
}
    
ostream& operator<<(ostream& os,const Mat_Interf& e)
{
    os << endl << "Material Properties" << endl;
    os << "-------------------" << endl;
    os << "Elastic modulus multiplier        = " << e.data.alpha << endl << endl;
	  os << "Damage coupling multiplier        = " << e.data.beta << endl << endl;
	  os << "Maximum resisting stress of crack = " << e.data.ft		<< endl << endl;
	  os << "Fracture energy                   = " << e.data.Gf << endl << endl;

    return os;
}
    
}
