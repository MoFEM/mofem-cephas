/** \file Core.cpp
 * \brief Mylti-index containers, data structures and other low-level functions
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
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <MaterialBlocks.hpp>
#include <CubitBCData.hpp>
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
#include <FieldInterface.hpp>
#include <MeshRefinment.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

// Interfaces
#include <TetGenInterface.hpp>
#ifdef WITH_NETGEN
  namespace nglib {
  #include <nglib.h>
  }
  using namespace nglib;
  #include <NetGenInterface.hpp>
#endif

#include <NodeMerger.hpp>
#include <PrismsFromSurfaceInterface.hpp>

#include <boost/scoped_ptr.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <BitLevelCoupler.hpp>

namespace MoFEM {

//const static int debug = 1;

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
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown inteface");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::query_interface_type(const std::type_info& type,void*& ptr) {
  PetscFunctionBegin;

  // TetGen
  #ifdef WITH_TETGEN
  if(type == typeid(TetGenInterface)) {
    if(iFaces.find(IDD_MOFEMTetGegInterface.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMTetGegInterface.uUId.to_ulong()] = new TetGenInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMTetGegInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  // NetGen
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
    if(iFaces.find(IDD_MOFEMNodeMerger.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMNodeMerger.uUId.to_ulong()] = new NodeMergerInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMNodeMerger.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //BitLevelCoupler
  if(type == typeid(BitLevelCouplerInterface)) {
    if(iFaces.find(IDD_MOFEMBitLevelCoupler.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMBitLevelCoupler.uUId.to_ulong()] = new BitLevelCouplerInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMBitLevelCoupler.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Create prism elements from surface Elements
  if(type == typeid(PrismsFromSurfaceInterface)) {
    if(iFaces.find(IDD_MOFEMPrismsFromSurface.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMPrismsFromSurface.uUId.to_ulong()] = new PrismsFromSurfaceInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMPrismsFromSurface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(MeshRefinment)) {
    ptr = static_cast<MeshRefinment*>(this);
  } else if(type == typeid(SeriesRecorder)) {
    ptr = static_cast<SeriesRecorder*>(this);
  } else if(type == typeid(PrismInterface)) {
    ptr = static_cast<PrismInterface*>(this);
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown inteface");
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

    if(p == PETSC_ERROR_INITIAL) {
      error_printf_hilight();
      (*PetscErrorPrintf)("--------------------- MoFEM Error Message---------------------------------------------------------------------------\n");
      (*PetscErrorPrintf)("MoFEM version %d.%d.%d\n",MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
      (*PetscErrorPrintf)("MoFEM git commit id %s\n",GIT_SHA1_NAME);
      (*PetscErrorPrintf)("See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/guidelines_bug_reporting.html for bug reporting.\n");
      (*PetscErrorPrintf)("See http://userweb.eng.gla.ac.uk/lukasz.kaczmarczyk/MoFem/html/faq_and_bugs.html for trouble shooting.\n");
      error_printf_normal();

    }

    PetscTraceBackErrorHandler(PETSC_COMM_SELF,line,fun,file,n,p,mess,ctx);

    PetscBool ismain,isunknown;

    PetscStrncmp(fun,"main",4,&ismain);
    PetscStrncmp(fun,"unknown",7,&isunknown);

    if(ismain || isunknown) {

      stringstream strs_version;
      strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;

      error_printf_hilight();
      (*PetscErrorPrintf)("----------MoFEM End of Error Message -------send entire error message to CMatGU <cmatgu@googlegroups.com> ----------\n");
      error_printf_normal();

    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before aborting so will not accidently kill process 0*/
    PetscSleep(10.0);
    abort();

  }

  PetscFunctionReturn(n);
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
  rval = moab.tag_get_handle(
    "_MoFEM_VERSION",version.size()*sizeof(char),MB_TYPE_OPAQUE,
    th_version,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,NULL
  ); CHKERR_THROW(rval);
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
  rval = moab.tag_get_handle(
    "_RefBitLevel",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level
  ); CHKERR_THROW(rval);
  BitRefLevel def_bit_level_mask = BitRefLevel().set();
  rval = moab.tag_get_handle(
    "_RefBitLevelMask",sizeof(BitRefLevel),MB_TYPE_OPAQUE,
    th_RefBitLevel_Mask,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level_mask
  ); CHKERR_THROW(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle(
    "_RefBitEdge",sizeof(BitRefEdges),MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_egde
  ); CHKERR_THROW(rval);

  //Tags Field
  const unsigned long int def_id = 0;
  rval = moab.tag_get_handle(
    "_FieldId",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FieldId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  FieldSpace def_space = LASTSPACE;
  rval = moab.tag_get_handle(
    "_FieldSpace",sizeof(FieldSpace),MB_TYPE_OPAQUE,
    th_FieldSpace,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_space
  ); CHKERR_THROW(rval);
  const int def_val_len = 0;
  rval = moab.tag_get_handle(
    "_FieldName",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FieldName_DataNamePrefix",def_val_len,MB_TYPE_OPAQUE,
    th_FieldName_DataNamePrefix,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);

  //Tags FE
  rval = moab.tag_get_handle(
    "_FEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_FEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FEName",def_val_len,MB_TYPE_OPAQUE,
    th_FEName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FEIdCol",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FEIdRow",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FEIdData",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdData,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);

  //Tags Problem
  rval = moab.tag_get_handle("_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
    th_ProblemId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_ProblemFEId,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_id
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle("_ProblemName",def_val_len,MB_TYPE_OPAQUE,
    th_ProblemName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);
  DofIdx def_nbdofs = 0;
  rval = moab.tag_get_handle(
    "_ProblemNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_ProblemNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_ProblemLocalNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_ProblemGhostNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofRow,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_ProblemLocalNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_ProblemGhostNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofCol,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_nbdofs
  ); CHKERR_THROW(rval);

  //Coordinate systems
  const int def_coord_system = NO_CORD_SYSTEM_ON_THE_MESHSET;
  rval = moab.tag_get_handle(
    "_CoordSysId",1,MB_TYPE_INTEGER,th_CoordSystem,MB_TAG_CREAT|MB_TAG_SPARSE,&def_coord_system
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_FieldCoordSysId",1,MB_TYPE_INTEGER,th_FieldCoordSystem,MB_TAG_CREAT|MB_TAG_SPARSE,&def_coord_system
  ); CHKERR_THROW(rval);
  rval = moab.tag_get_handle(
    "_CoordSysName",def_val_len,MB_TYPE_OPAQUE,
    th_CoordSysName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);

  //Global Variables
  //Fields
  int def_shift = 1;
  rval = moab.tag_get_handle("_FieldShift",1,MB_TYPE_INTEGER,th_FieldShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  const void* tag_data[1];
  rval = moab.tag_get_by_ptr(th_FieldShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  fShift = (int*)tag_data[0];
  //FE
  rval = moab.tag_get_handle("_FEShift",1,MB_TYPE_INTEGER,th_FEShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_FEShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  feShift = (int*)tag_data[0];
  //Problem
  rval = moab.tag_get_handle("_ProblemShift",1,MB_TYPE_INTEGER,th_ProblemShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERR_THROW(rval);
  rval = moab.tag_get_by_ptr(th_ProblemShift,&root_meshset,1,tag_data); CHKERR_THROW(rval);
  pShift = (int*)tag_data[0];
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
    bhTag_header,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_uint_zero[0]
  ); CHKERR_THROW(rval);
  Tag block_attribs;
  int def_Block_Attributes_length = 0;
  rval = moab.tag_get_handle(BLOCK_ATTRIBUTES,def_Block_Attributes_length,MB_TYPE_DOUBLE,
    block_attribs,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);
  Tag entity_name_tag;
  rval = moab.tag_get_handle(
    NAME_TAG_NAME,NAME_TAG_SIZE,MB_TYPE_OPAQUE,entity_name_tag,MB_TAG_SPARSE|MB_TAG_CREAT
  ); CHKERR_THROW(rval);
  //Series
  rval = moab.tag_get_handle("_SeriesName",def_val_len,MB_TYPE_OPAQUE,
    th_SeriesName,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES|MB_TAG_VARLEN,NULL
  ); CHKERR_THROW(rval);
  //For VTK files
  int def_elem_type = MBMAXTYPE;
  rval = moab.tag_get_handle(
    "ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type
  ); CHKERR_THROW(rval);
  //
  clearMap();
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
  if(miit==set.end()) {
    THROW_AT_LINE("field < "+name+" > not in database (top tip: check spelling)");
  }
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
  if(miit==set.end()) THROW_AT_LINE("field not in database (top tip: check spelling)");
  return miit->meshSet;
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
    throw MoFEMException(
      MOFEM_NOT_FOUND,
      string("field < "+name+" > not in databse (top tip: check spelling)").c_str()
    );
  }
  return &*miit;
}
BitFieldId Core::getFieldShift() {
  if(*fShift >= BITFIELDID_SIZE) {
    char msg[] = "number of fields exceeded";
    PetscTraceBackErrorHandler(
      comm,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,msg,PETSC_NULL);
    PetscMPIAbortErrorHandler(comm,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,msg,PETSC_NULL);
  }
  return BitFieldId().set(((*fShift)++)-1);
}
BitFEId Core::getFEShift() {
  assert((unsigned int)*feShift<BitFEId().set().to_ulong());
  return BitFEId(1<<(((*feShift)++)-1));
}
BitProblemId Core::getProblemShift() {
  assert((unsigned int)*pShift<BitProblemId().set().to_ulong());
  return BitProblemId(1<<(((*pShift)++)-1));
}
PetscErrorCode Core::clearMap() {
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
  cubitMeshsets.clear();
  PetscFunctionReturn(0);
}

PetscErrorCode Core::addPrismToDatabase(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(RefMoFEMEntity(moab,prism));
    if(p_ent.second) {
      pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      p_MoFEMFiniteElement = refinedFiniteElements.insert(
	      ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ent.first))
      );
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
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::synchronise_entities(Range &ents,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  //ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

  //make a buffer
  vector<vector<EntityHandle> > sbuffer(sIze);

  Range::iterator eit = ents.begin();
  for(;eit!=ents.end();eit++) {

    RefMoFEMEntity_multiIndex::iterator meit;
    meit = refinedEntities.get<Ent_mi_tag>().find(*eit);
    if(meit == refinedEntities.get<Ent_mi_tag>().end()) {
      continue;
      SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
        "rank %d entity %lu not exist on database, local entity can not be found for this owner",
        rAnk,*eit);
    }

    unsigned char pstatus = meit->get_pstatus();

    if(pstatus == 0) continue;

    if(verb>1) {
      ostringstream zz;
      zz << "pstatus " <<  bitset<8>(pstatus) << " ";
      PetscSynchronizedPrintf(comm,"%s",zz.str().c_str());
    }

    for(int proc = 0; proc<MAX_SHARING_PROCS && -1 != meit->get_sharing_procs_ptr()[proc]; proc++) {
      if(meit->get_sharing_procs_ptr()[proc] == -1) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_IMPOSIBLE_CASE,"sharing processor not set");
      }
      if(meit->get_sharing_procs_ptr()[proc] == rAnk) {
        continue;
      }
      EntityHandle handle_on_sharing_proc = meit->get_sharing_handlers_ptr()[proc];
      sbuffer[meit->get_sharing_procs_ptr()[proc]].push_back(handle_on_sharing_proc);
      if(verb>1) {
	PetscSynchronizedPrintf(comm,"send %lu (%lu) to %d at %d\n",
	  meit->get_ref_ent(),handle_on_sharing_proc,meit->get_sharing_procs_ptr()[proc],rAnk);
      }
      if(!(pstatus&PSTATUS_MULTISHARED)) {
	break;
      }
    }


  }

  int nsends = 0; 			// number of messages to send
  vector<int> sbuffer_lengths(sIze); 	// length of the message to proc
  const size_t block_size = sizeof(EntityHandle)/sizeof(int);
  for(int proc  = 0;proc<sIze;proc++) {

    if(!sbuffer[proc].empty()) {

      sbuffer_lengths[proc] = sbuffer[proc].size()*block_size;
      nsends++;

    } else {

      sbuffer_lengths[proc] = 0;

    }

  }

  // Make sure it is a PETSc comm
  ierr = PetscCommDuplicate(comm,&comm,NULL); CHKERRQ(ierr);

  vector<MPI_Status> status(sIze);

  // Computes the number of messages a node expects to receive
  int nrecvs;	// number of messages received
  ierr = PetscGatherNumberOfMessages(comm,NULL,&sbuffer_lengths[0],&nrecvs); CHKERRQ(ierr);

  // Computes info about messages that a MPI-node will receive, including (from-id,length) pairs for each message.
  int *onodes;		// list of node-ids from which messages are expected
  int *olengths;	// corresponding message lengths
  ierr = PetscGatherMessageLengths(comm,nsends,nrecvs,&sbuffer_lengths[0],&onodes,&olengths);  CHKERRQ(ierr);

  // Gets a unique new tag from a PETSc communicator. All processors that share
  // the communicator MUST call this routine EXACTLY the same number of times.
  // This tag should only be used with the current objects communicator; do NOT
  // use it with any other MPI communicator.
  int tag;
  ierr = PetscCommGetNewTag(comm,&tag); CHKERRQ(ierr);

  // Allocate a buffer sufficient to hold messages of size specified in
  // olengths. And post Irecvs on these buffers using node info from onodes
  int **rbuf;		// must bee freed by user
  MPI_Request *r_waits; // must bee freed by user
  // rbuf has a pointers to messages. It has size of of nrecvs (number of
  // messages) +1. In the first index a block is allocated,
  // such that rbuf[i] = rbuf[i-1]+olengths[i-1].
  ierr = PetscPostIrecvInt(comm,tag,nrecvs,onodes,olengths,&rbuf,&r_waits); CHKERRQ(ierr);

  MPI_Request *s_waits; // status of sens messages
  ierr = PetscMalloc1(nsends,&s_waits); CHKERRQ(ierr);

  // Send messeges
  for(int proc=0,kk=0; proc<sIze; proc++) {
    if(!sbuffer_lengths[proc]) continue; // no message to send to this proc
    ierr = MPI_Isend(
      &(sbuffer[proc])[0], 	// buffer to send
      sbuffer_lengths[proc], 	// message length
      MPIU_INT,proc,       	// to proc
      tag,comm,s_waits+kk); CHKERRQ(ierr);
    kk++;
  }

  // Wait for received
  if(nrecvs) {
    ierr = MPI_Waitall(nrecvs,r_waits,&status[0]);CHKERRQ(ierr);
  }
  // Wait for send messages
  if(nsends) {
    ierr = MPI_Waitall(nsends,s_waits,&status[0]);CHKERRQ(ierr);
  }

  if(verb>0) {
    PetscSynchronizedPrintf(comm,"nb. before ents %u\n",ents.size());
  }

  // synchronise range
  for(int kk = 0;kk<nrecvs;kk++) {

    int len = olengths[kk];
    int *data_from_proc = rbuf[kk];

    for(int ee = 0;ee<len;ee+=block_size) {

      EntityHandle ent;
      bcopy(&data_from_proc[ee],&ent,sizeof(EntityHandle));
      RefMoFEMEntity_multiIndex::index<Ent_mi_tag>::type::iterator meit;
      meit = refinedEntities.get<Ent_mi_tag>().find(ent);
      if(meit == refinedEntities.get<Ent_mi_tag>().end()) {
	SETERRQ2(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,
	  "rank %d entity %lu not exist on database, local entity can not be found for this owner",rAnk,ent);
      }
      if(verb>2) {
	PetscSynchronizedPrintf(comm,"received %ul (%ul) from %d at %d\n",meit->get_ref_ent(),ent,onodes[kk],rAnk);
      }
      ents.insert(meit->get_ref_ent());

    }

  }

  if(verb>0) {
    PetscSynchronizedPrintf(comm,"nb. after ents %u\n",ents.size());
  }


  // Cleaning
  ierr = PetscFree(s_waits); CHKERRQ(ierr);
  ierr = PetscFree(rbuf[0]); CHKERRQ(ierr);
  ierr = PetscFree(rbuf); CHKERRQ(ierr);
  ierr = PetscFree(r_waits); CHKERRQ(ierr);
  ierr = PetscFree(onodes); CHKERRQ(ierr);
  ierr = PetscFree(olengths); CHKERRQ(ierr);

  if(verb>0) {
    PetscSynchronizedFlush(comm,PETSC_STDOUT);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode Core::rebuild_database(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clearMap(); CHKERRQ(ierr);
  ierr = initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::initialiseDatabseInformationFromMesh(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range meshsets;
  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERR_PETSC(rval);
  //loop all meshsets in moab to find meshets with boundary conditions
  Range::iterator mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    try {
      //check if meshset is cubit meshset
      CubitMeshSets base_meshset(moab,*mit);
      if((base_meshset.cubit_bc_type&CubitBCType(NODESET|SIDESET|BLOCKSET)).any()) {
        pair<CubitMeshSet_multiIndex::iterator,bool> p = cubitMeshsets.insert(base_meshset);
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"meshset not inserted");
        }
        if(verb > 0) {
          ostringstream ss;
          ss << "read cubit" << base_meshset << endl;
          //PetscSynchronizedPrintf(comm,ss.str().c_str());
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
  }
  //loop all meshsehset to find coordinate system
  mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    try {
      int cs_id = CARTESIAN_COORD_SYSTEM;
      rval = moab.tag_get_data(th_CoordSystem,&*mit,1,&cs_id); CHKERR_PETSC(rval);
      if(cs_id!=NO_CORD_SYSTEM_ON_THE_MESHSET) {
        CoordSys coord_sys(moab,*mit);
        pair<CoordSys_multiIndex ::iterator,bool> p = coordinateSystems.insert(coord_sys);
        if(!p.second) {
          SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"meshset to coord system not inserted");
        }
      }
    } catch (MoFEMException const &e) {
      SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
    }
  }
  { // Create cartesian coordinate system if not exist
    CoordSys_multiIndex::index<CoordSysID_mi_tag>::type::iterator csit;
    csit = coordinateSystems.get<CoordSysID_mi_tag>().find(CARTESIAN_COORD_SYSTEM);
    if(csit==coordinateSystems.get<CoordSysID_mi_tag>().end()) {
      EntityHandle meshset;
      rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERR_PETSC(rval);
      int id = CARTESIAN_COORD_SYSTEM;
      rval = moab.tag_set_data(th_CoordSystem,&meshset,1,&id); CHKERR_PETSC(rval);
      string sys_name_str = "CARTESIAN_GENERIC";
      void const* sys_name[] = { sys_name_str.c_str() };
      int sys_name_size[1];
      sys_name_size[0] = sys_name_str.size();
      rval = moab.tag_set_by_ptr(
        th_CoordSysName,&meshset,1,sys_name,sys_name_size
      ); CHKERR_PETSC(rval);
      CoordSys coord_sys(moab,meshset);
      pair<CoordSys_multiIndex ::iterator,bool> p = coordinateSystems.insert(coord_sys);
      if(!p.second) {
        SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"MeshSet to coord system not inserted");
      }
    }
  }
  //PetscSynchronizedFlush(comm,PETSC_STDOUT);
  CoordSys_multiIndex::index<CoordSysID_mi_tag>::type::iterator cartesian_it;
  cartesian_it = coordinateSystems.get<CoordSysID_mi_tag>().find(CARTESIAN_COORD_SYSTEM);
  if(cartesian_it==coordinateSystems.get<CoordSysID_mi_tag>().end()) {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Generic Cartesian system not found");
  }
  mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    BitFieldId field_id;
    //get bit id form field tag
    rval = moab.tag_get_data(th_FieldId,&*mit,1,&field_id); CHKERR_PETSC(rval);
    //check if meshset if field meshset
    if(field_id!=0) {
      pair<MoFEMField_multiIndex::iterator,bool> p;
      try {
        int coord_sys_id;
        rval = moab.tag_get_data(th_FieldCoordSystem,&*mit,1,&coord_sys_id); CHKERR_PETSC(rval);
        CoordSys_multiIndex::index<CoordSysID_mi_tag>::type::iterator cs_it = cartesian_it;
        if(coord_sys_id) {
          cs_it = coordinateSystems.get<CoordSysID_mi_tag>().find(coord_sys_id);
          if(cs_it==coordinateSystems.get<CoordSysID_mi_tag>().end()) {
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Unknown Coordinate System");
          }
        }
        p = moabFields.insert(MoFEMField(moab,*mit,&*cs_it));
        if(verb > 0) {
          ostringstream ss;
          ss << "read field " << *p.first << endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
      if(p.first->get_space()==NOFIELD) {
        assert(p.first->meshSet == *mit);
        //add field to ref ents
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*mit));
        NOT_USED(p_ref_ent);
      } else {
        Range ents;
        rval = moab.get_entities_by_handle(*mit,ents,false); CHKERR_PETSC(rval);
        if(verb > 1) {
          ostringstream ss;
          ss << "read field ents " << ents.size() << endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*eit));
          try {
            MoFEMEntity moabent(moab,&*p.first,&*p_ref_ent.first);
            pair<MoFEMEntity_multiIndex::iterator,bool> p_ent = entsMoabField.insert(moabent);
            NOT_USED(p_ent);
          } catch (const std::exception& ex) {
            ostringstream ss;
            ss << ex.what() << endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
          }
        }
      }
    }
    BitFieldId fe_id;
    //get bit id from fe tag
    rval = moab.tag_get_data(th_FEId,&*mit,1,&fe_id); CHKERR_PETSC(rval);
    //check if meshset is finite element meshset
    if(fe_id!=0) {
      pair<MoFEMFiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(MoFEMFiniteElement(moab,*mit));
      if(verb > 0) {
        ostringstream ss;
        ss << "read finite element " << *p.first << endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
      NOT_USED(p);
      assert(p.first->meshset == *mit);
      Range ents;
      rval = moab.get_entities_by_type(*mit,MBENTITYSET,ents,false); CHKERR_PETSC(rval);
      rval = moab.get_entities_by_handle(*mit,ents,true); CHKERR_PETSC(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        pair<RefMoFEMEntity_multiIndex::iterator,bool> p_ref_ent = refinedEntities.insert(RefMoFEMEntity(moab,*eit));
        pair<RefMoFEMElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
        try {
          switch (moab.type_from_handle(*eit)) {
            case MBVERTEX:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_VERTEX(moab,&*p_ref_ent.first)));
            break;
            case MBEDGE:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_EDGE(moab,&*p_ref_ent.first)));
            break;
            case MBTRI:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TRI(moab,&*p_ref_ent.first)));
            break;
            case MBTET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_TET(moab,&*p_ref_ent.first)));
            break;
            case MBPRISM:
            ierr = addPrismToDatabase(*eit,verb); CHKERRQ(ierr);
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_PRISM(moab,&*p_ref_ent.first)));
            break;
            case MBENTITYSET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefMoFEMElement(new RefMoFEMElement_MESHSET(moab,&*p_ref_ent.first)));
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
          }
          if(p_MoFEMFiniteElement.second) {
            //PetscPrintf(comm,"Warrning: this entity should be already in refined finite elements database");
            //SETERRQ(PETSC_COMM_SELF,1,"data inconsistency, this entity should be already in refined finite elements database");
          }
        } catch (MoFEMException const &e) {
          SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
        }
      }
    }
    BitProblemId problem_id;
    //get bit id form problem tag
    rval = moab.tag_get_data(th_ProblemId,&*mit,1,&problem_id); CHKERR_PETSC(rval);
    //check if meshset if problem meshset
    if(problem_id!=0) {
      pair<MoFEMProblem_multiIndex::iterator,bool> p = moFEMProblems.insert(MoFEMProblem(moab,*mit));
      if(verb > 0) {
        ostringstream ss;
        ss << "read problem " << *p.first << endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
    }
    //check if meshset is Series meshset
    {
      const void* tag_name_data;
      int tag_name_size;
      rval = moab.tag_get_by_ptr(th_SeriesName,&*mit,1,(const void **)&tag_name_data,&tag_name_size);
      if(rval == MB_SUCCESS) {
        pair<Series_multiIndex::iterator,bool> p = sEries.insert(MoFEMSeries(moab,*mit));
        if(verb > 0) {
          ostringstream ss;
          ss << "read series " << *p.first << endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    }
  }
  //build ref entities meshset
  for(int dd = 0;dd<=3;dd++) {
    Range ents;
    rval = moab.get_entities_by_dimension(0,dd,ents,false); CHKERR_PETSC(rval);
    Range::iterator eit = ents.begin();
    for(;eit!=ents.end();eit++) {
      switch (moab.type_from_handle(*eit)) {
        case MBVERTEX:
        case MBEDGE:
        case MBTRI:
        case MBTET:
        case MBPRISM:
        break;
        default:
        continue;
      }
      RefMoFEMEntity mofem_ent(moab,*eit);
      BitRefLevel bit = mofem_ent.get_BitRefLevel();
      if(bit.none()) {
        continue;
      }
      pair<RefMoFEMEntity_multiIndex::iterator,bool> p;
      p = refinedEntities.insert(mofem_ent);
    }
  }
  //build series steps
  for(Series_multiIndex::iterator sit = sEries.begin();sit!=sEries.end();sit++) {
    int nb_steps;
    ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
    int ss = 0;
    for(;ss<nb_steps;ss++) {
      pair<SeriesStep_multiIndex::iterator,bool> p = seriesSteps.insert(MoFEMSeriesStep(moab,&*sit,ss));
      if(verb > 0) {
        ostringstream ss;
        ss << "add series step " << *p.first << endl;
        PetscPrintf(comm,ss.str().c_str());
      }
    }
  }
  if(verb > 2) {
    list_fields();
    list_finite_elements();
    list_problem();
  }
  PetscFunctionReturn(0);
}

}
