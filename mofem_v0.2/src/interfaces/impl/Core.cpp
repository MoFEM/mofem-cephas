/** \file Core.cpp
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


#include <moab/ParallelComm.hpp>

#include <boost/ptr_container/ptr_map.hpp>

#include <petscsys.h>
#include <petscvec.h> 
#include <petscmat.h> 
#include <petscsnes.h> 
#include <petscts.h>
#include <petscconfiginfo.h> 

#include <version.h>
#include <definitions.h>
#include <h1_hdiv_hcurl_l2.h>

#include <Common.hpp>

#include <LoopMethods.hpp>
#include <Core.hpp>

#include <CoreDataStructures.hpp>

#include <TetGenInterface.hpp>

#ifdef WITH_NETGEN
  namespace nglib {
  #include <nglib.h>
  }
  using namespace nglib;
  #include <NetGenInterface.hpp>
#endif

#include <NodeMerger.hpp>
#include <BitLevelCoupler.hpp>

namespace MoFEM {

const static int debug = 1;

PetscErrorCode print_MoFem_verison(MPI_Comm comm) {
  PetscFunctionBegin;
  PetscPrintf(comm,"version %d.%d.%d\n",MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
  PetscPrintf(comm,"git commit id %s\n",GIT_SHA1_NAME);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::queryInterface(const MOFEMuuid& uuid,FieldUnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMPrismInterface) {
    *iface = dynamic_cast<PrismInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMMeshRefine) {
    *iface = dynamic_cast<MeshRefinment*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMSeriesRecorder) {
    *iface = dynamic_cast<SeriesRecorder*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMFieldInterface) {
    *iface = dynamic_cast<FieldInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<FieldInterface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown inteface");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::query_interface_type(const std::type_info& type,void*& ptr) {
  PetscFunctionBegin;
  
  #ifdef WITH_TETGEN
  if(type == typeid(TetGenInterface)) {
    if(iFaces.find(IDD_MOFEMTetGegInterface.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMTetGegInterface.uUId.to_ulong()] = new TetGenInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMTetGegInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  #ifdef WITH_NETGEN
  if(type == typeid(NetGenInterface)) {
    if(iFaces.find(IDD_MOFEMNetGegInterface.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMNetGegInterface.uUId.to_ulong()] = new NetGenInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMNetGegInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  //Node merger
  if(type == typeid(NodeMergerInterface)) {
    if(iFaces.find(IDD_MOFENNodeMerger.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFENNodeMerger.uUId.to_ulong()] = new NodeMergerInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFENNodeMerger.uUId.to_ulong());
    PetscFunctionReturn(0);
  } 

  if(type == typeid(MeshRefinment)) {
    ptr = static_cast<MeshRefinment*>(this);
  } else if(type == typeid(SeriesRecorder)) {
    ptr = static_cast<SeriesRecorder*>(this);
  } else if(type == typeid(PrismInterface)) {
    ptr = static_cast<PrismInterface*>(this);
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"unknown inteface");
  }
  PetscFunctionReturn(0);
}

bool Core::isGloballyInitialised = false;

static void error_printf_hilight(void) {
#if defined(PETSC_HAVE_UNISTD_H) && defined(PETSC_USE_ISATTY)
  if (PetscErrorPrintf == PetscErrorPrintfDefault) {
    if (isatty(fileno(PETSC_STDERR))) fprintf(PETSC_STDERR,"\033[1;32m");
  }
#endif
}

static void error_printf_normal(void) {
#if defined(PETSC_HAVE_UNISTD_H) && defined(PETSC_USE_ISATTY)
  if (PetscErrorPrintf == PetscErrorPrintfDefault) {
    if (isatty(fileno(PETSC_STDERR))) fprintf(PETSC_STDERR,"\033[0;39m\033[0;49m");
  }
#endif
}

PetscErrorCode mofem_error_handler(MPI_Comm comm,int line,const char *fun,const char *file,PetscErrorCode n,PetscErrorType p,const char *mess,void *ctx) {
  PetscFunctionBegin; 

  int rank = 0;
  if (comm != PETSC_COMM_SELF) MPI_Comm_rank(comm,&rank);

  if(!rank) {

    PetscBool ismain,isunknown;
  
    PetscStrncmp(fun,"main",4,&ismain); 
    PetscStrncmp(fun,"unknown",7,&isunknown); 

    if(ismain || isunknown) { 

      stringstream strs_version;
      strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;


      error_printf_hilight();
      (*PetscErrorPrintf)("--------------------- MoFEM Error Message---------------------------------------------------------------------------\n"); 
      (*PetscErrorPrintf)("MoFEM version %d.%d.%d\n",MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD); 
      (*PetscErrorPrintf)("MoFEM git commit id %s\n",GIT_SHA1_NAME); 
      (*PetscErrorPrintf)("See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/guidelines_bug_reporting.html for bug reporting.\n");
      (*PetscErrorPrintf)("See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/faq_and_bugs.html for trouble shooting.\n");
      error_printf_normal(); 

      PetscTraceBackErrorHandler(PETSC_COMM_SELF,line,fun,file,n,p,mess,ctx);

      error_printf_hilight();
      (*PetscErrorPrintf)("----------MoFEM End of Error Message -------send entire error message to CMatGU <cmatgu@googlegroups.com> ----------\n"); 
      error_printf_normal(); 

    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before aborting so will not accidently kill process 0*/
    PetscSleep(10.0);
    abort();

  }

  PetscFunctionReturn(0);
}

Core::Core(Interface& _moab,MPI_Comm _comm,int _verbose): 
  moab(_moab),comm(_comm),verbose(_verbose) {

  if(!isGloballyInitialised) {
    PetscPushErrorHandler(mofem_error_handler,PETSC_NULL);
    isGloballyInitialised = true;
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,comm);
  MPI_Comm_size(comm,&sIze);
  MPI_Comm_rank(comm,&rAnk);

  const EntityHandle root_meshset = moab.get_root_set();
  if(verbose>0) {
    print_MoFem_verison(comm);
  }

  // Version
  Tag th_version;
  stringstream strs_version;
  strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;
  string version = strs_version.str();
  rval = moab.tag_get_handle("_MoFEM_VERSION",version.size()*sizeof(char),MB_TYPE_OPAQUE,
    th_version,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,NULL); CHKERR_THROW(rval);
  const char *ptr_version = version.c_str();
  rval = moab.tag_set_data(th_version,&root_meshset,1,ptr_version); CHKERR_THROW(rval);

  //tags saved in vtk-files
  const int def_part = -1;
  rval = moab.tag_get_handle("PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part); CHKERR_THROW(rval);

  //Tags Ref
  EntityHandle def_handle = 0;
  rval = moab.tag_get_handle("_RefParentHandle",1,MB_TYPE_HANDLE,th_RefParentHandle,MB_TAG_CREAT|MB_TAG_SPARSE,&def_handle); CHKERR_THROW(rval);
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle("_RefType",2,MB_TYPE_INTEGER,th_RefType,MB_TAG_CREAT|MB_TAG_SPARSE,def_type); CHKERR_THROW(rval);
  BitRefLevel def_bit_level = 0;
  rval = moab.tag_get_handle("_RefBitLevel",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); CHKERR_THROW(rval);
  BitRefLevel def_bit_level_mask = BitRefLevel().set();
  rval = moab.tag_get_handle("_RefBitLevelMask",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel_Mask,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_mask); CHKERR_THROW(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle("_RefBitEdge",sizeof(BitRefEdges),MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_egde); CHKERR_THROW(rval);
    
  //Tags Field
  const unsigned long int def_id = 0;
  rval = moab.tag_get_handle("_FieldId",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FieldId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  FieldSpace def_space = LASTSPACE;
  rval = moab.tag_get_handle("_FieldSpace",sizeof(FieldSpace),MB_TYPE_OPAQUE,
    th_FieldSpace,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_space); CHKERR_THROW(rval);
  const int def_val_len = 0;
  rval = moab.tag_get_handle("_FieldName",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FieldName_DataNamePrefix",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName_DataNamePrefix,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);

  //Tags FE
  rval = moab.tag_get_handle("_FEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_FEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEName",def_val_len,MB_TYPE_OPAQUE,
    th_FEName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdCol",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdRow",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_FEIdData",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);

  //Tags Problem
  rval = moab.tag_get_handle("_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
    th_ProblemId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_ProblemFEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemName",def_val_len,MB_TYPE_OPAQUE,
    th_ProblemName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  DofIdx def_nbdofs = 0;
  rval = moab.tag_get_handle("_ProblemNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemLocalNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemGhostNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs); CHKERR_THROW(rval);

  //Global Variables
  //Fields
  int def_shift = 1;
  rval = moab.tag_get_handle("_FieldShift",1,MB_TYPE_INTEGER,th_FieldShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  const void* tag_data[1];
  rval = moab.tag_get_by_ptr(th_FieldShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  f_shift = (int*)tag_data[0];
  //FE
  rval = moab.tag_get_handle("_FEShift",1,MB_TYPE_INTEGER,th_FEShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  MoFEMFiniteElement_shift = (int*)tag_data[0];
  //Problem
  rval = moab.tag_get_handle("_ProblemShift",1,MB_TYPE_INTEGER,th_ProblemShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  p_shift = (int*)tag_data[0];
  //SaftyNets
  int def_bool = 0;
  rval = moab.tag_get_handle("_MoFEMBuild",1,MB_TYPE_INTEGER,th_MoFEMBuild,MB_TAG_CREAT|MB_TAG_MESH,&def_bool); 
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  rval = moab.tag_get_by_ptr(th_MoFEMBuild,&root_meshset,1,(const void **)&build_MoFEM); CHKERR_THROW(rval);
  //Meshsets
  int default_val = -1;
  rval = moab.tag_get_handle(DIRICHLET_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    nsTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(NEUMANN_SET_TAG_NAME,1, MB_TYPE_INTEGER, 
    ssTag, MB_TAG_SPARSE|MB_TAG_CREAT, &default_val); CHKERR_THROW(rval);
  const int def_bc_data_len = 0;
  std::string tag_name = std::string(DIRICHLET_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    nsTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  tag_name = std::string(NEUMANN_SET_TAG_NAME)+"__BC_DATA";
  rval = moab.tag_get_handle(tag_name.c_str(),def_bc_data_len,MB_TYPE_OPAQUE,
    ssTag_data,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER,
    bhTag,MB_TAG_SPARSE|MB_TAG_CREAT,&default_val); CHKERR_THROW(rval);
  std::vector<unsigned int> def_uint_zero(3,0);
  rval= moab.tag_get_handle(BLOCK_HEADER,3*sizeof(unsigned int),MB_TYPE_INTEGER,
    bhTag_header,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_uint_zero[0]); CHKERR_THROW(rval); 
  Tag block_attribs;
  int def_Block_Attributes_length = 0;
  rval = moab.tag_get_handle(BLOCK_ATTRIBUTES,def_Block_Attributes_length,MB_TYPE_DOUBLE,
    block_attribs,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval); 
  Tag entity_name_tag;
  rval = moab.tag_get_handle(
    NAME_TAG_NAME,NAME_TAG_SIZE,MB_TYPE_OPAQUE,entity_name_tag,MB_TAG_SPARSE|MB_TAG_CREAT); CHKERR_THROW(rval);
  //Series
  rval = moab.tag_get_handle("_SeriesName",def_val_len,MB_TYPE_OPAQUE,
    th_SeriesName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL); CHKERR_THROW(rval);
  //For VTK files
  int def_elem_type = MBMAXTYPE;
  rval = moab.tag_get_handle("ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type); CHKERR_THROW(rval); 
  //
  clear_map();
  initialiseDatabseInformationFromMesh(verbose); 
  // Petsc Logs
  PetscLogEventRegister("FE_preProcess",0,&USER_EVENT_preProcess);
  PetscLogEventRegister("FE_operator",0,&USER_EVENT_operator);
  PetscLogEventRegister("FE_postProcess",0,&USER_EVENT_postProcess);
  PetscLogEventRegister("FielCore_createMat",0,&USER_EVENT_createMat);

}
Core::~Core() {
}
Interface& Core::get_moab() {
  return moab;
}
MPI_Comm Core::get_comm() {
  return comm;
}
BitFieldId Core::get_BitFieldId(const string& name) const {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) THROW_AT_LINE("field < "+name+" > not in databse (top tip: check spelling)");
  return miit->get_id();
}
string Core::get_BitFieldId_name(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  return miit->get_name();
}
EntityHandle Core::get_field_meshset(const BitFieldId id) const {
  typedef MoFEMField_multiIndex::index<BitFieldId_mi_tag>::type field_set_by_id;
  const field_set_by_id &set = moabFields.get<BitFieldId_mi_tag>();
  field_set_by_id::iterator miit = set.find(id);
  if(miit==set.end()) THROW_AT_LINE("field not in databse (top tip: check spelling)");
  return miit->meshset;
}
EntityHandle Core::get_field_meshset(const string& name) const {
  return get_field_meshset(get_BitFieldId(name));
}

bool Core::check_field(const string &name) const {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) return false;
  return true;
}
const MoFEMField* Core::get_field_structure(const string& name) {
  typedef MoFEMField_multiIndex::index<FieldName_mi_tag>::type field_set_by_name;
  const field_set_by_name &set = moabFields.get<FieldName_mi_tag>();
  field_set_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    throw MofemException(MOFEM_NOT_FOUND,
      string("field < "+name+" > not in databse (top tip: check spelling)").c_str());
  }
  return &*miit;
}
BitFieldId Core::get_field_shift() {
  if(*f_shift >= BITFIELDID_SIZE) {
    char msg[] = "number of fields exceeded";
    PetscTraceBackErrorHandler(
      comm,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCT,PETSC_ERROR_INITIAL,msg,PETSC_NULL);
    PetscMPIAbortErrorHandler(comm,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCT,PETSC_ERROR_INITIAL,msg,PETSC_NULL);
  }
  return BitFieldId().set(((*f_shift)++)-1);
}
BitFEId Core::get_BitFEId() {
  assert((unsigned int)*MoFEMFiniteElement_shift<BitFEId().set().to_ulong());
  return BitFEId(1<<(((*MoFEMFiniteElement_shift)++)-1)); 
}
BitProblemId Core::get_problem_shift() {
  assert((unsigned int)*p_shift<BitProblemId().set().to_ulong());
  return BitProblemId(1<<(((*p_shift)++)-1)); 
}
PetscErrorCode Core::clear_map() {
  PetscFunctionBegin;
  refinedEntities.clear();
  refinedFiniteElements.clear();
  moabFields.clear();
  entsMoabField.clear();
  dofsMoabField.clear();
  finiteElements.clear();
  finiteElementsMoFEMEnts.clear();
  entFEAdjacencies.clear();
  moFEMProblems.clear();
  cubit_meshsets.clear();
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_prism_to_mofem_database(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(RefMoFEMEntity(moab,prism));
    if(p_ent.second) {
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      p_MoFEMFiniteElement = refinedFiniteElements.insert(
	ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first)));
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(prism,conn,num_nodes,true); CHKERR_THROW(rval);
      Range face_side3,face_side4;
      rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERR_PETSC(rval);
      rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERR_PETSC(rval);
      if(face_side3.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prism don't have side face 3");
      if(face_side4.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prims don't have side face 4");
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side3.begin());
      p_MoFEMFiniteElement.first->get_side_number_ptr(moab,*face_side4.begin());
    } 
  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }
  PetscFunctionReturn(0);
}

}
