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

#include <version.h>
#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

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

#include <UnknownInterface.hpp>
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

// Interfaces
#include <MeshRefinement.hpp>
#include <MeshsetsManager.hpp>
#include <CoordSystemsManager.hpp>
#include <TetGenInterface.hpp>
#include <MedInterface.hpp>
#include <NodeMerger.hpp>
#include <PrismsFromSurfaceInterface.hpp>

#include <boost/scoped_ptr.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <BitLevelCoupler.hpp>

extern "C" {
  void macro_is_depracted_using_deprecated_function() {}
}

namespace MoFEM {

//const static int debug = 1;

PetscErrorCode print_MoFem_verison(MPI_Comm comm) {
  PetscFunctionBegin;
  PetscPrintf(comm,"version %d.%d.%d\n",MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD);
  PetscPrintf(comm,"git commit id %s\n",GIT_SHA1_NAME);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::queryInterface(const MOFEMuuid& uuid,UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
  if(uuid == IDD_MOFEMPrismInterface) {
    *iface = dynamic_cast<PrismInterface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMMeshRefine) {
    *iface = dynamic_cast<MeshRefinement*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMSeriesRecorder) {
    *iface = dynamic_cast<SeriesRecorder*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMInterface) {
    *iface = dynamic_cast<Interface*>(this);
    PetscFunctionReturn(0);
  }
  if(uuid == IDD_MOFEMUnknown) {
    *iface = dynamic_cast<Interface*>(this);
    PetscFunctionReturn(0);
  }
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::query_interface_type(const std::type_info& type,void*& ptr) const {
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

  // MedInterface
  #ifdef WITH_MED
  if(type == typeid(MedInterface)) {
    if(iFaces.find(IDD_MOFEMMedInterface.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMMedInterface.uUId.to_ulong()] = new MedInterface(*this);
    }
    ptr = iFaces.at(IDD_MOFEMMedInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  //Meshsets manager
  if(type == typeid(MeshsetsManager)) {
    if(iFaces.find(IDD_MOFEMMeshsetsManager.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMMeshsetsManager.uUId.to_ulong()] = new MeshsetsManager(*this);
    }
    ptr = iFaces.at(IDD_MOFEMMeshsetsManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Cooordinate systems manager
  if(type == typeid(CoordSystemsManager)) {
    if(iFaces.find(IDD_MOFEMCoordsSystemsManager.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMCoordsSystemsManager.uUId.to_ulong()] = new CoordSystemsManager(*this);
    }
    ptr = iFaces.at(IDD_MOFEMCoordsSystemsManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

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

  if(type == typeid(MeshRefinement)) {
    if(iFaces.find(IDD_MOFEMMeshRefine.uUId.to_ulong()) == iFaces.end()) {
      iFaces[IDD_MOFEMMeshRefine.uUId.to_ulong()] = new MeshRefinement(*this);
    }
    ptr = iFaces.at(IDD_MOFEMMeshRefine.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(SeriesRecorder)) {
    ptr = static_cast<SeriesRecorder*>(const_cast<Core*>(this));
  } else if(type == typeid(PrismInterface)) {
    ptr = static_cast<PrismInterface*>(const_cast<Core*>(this));
  } else {
    SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
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
      (*PetscErrorPrintf)("See http://mofem.eng.gla.ac.uk/mofem/html/guidelines_bug_reporting.html for bug reporting.\n");
      (*PetscErrorPrintf)("See http://mofem.eng.gla.ac.uk/mofem/html/faq_and_bugs.html for trouble shooting.\n");
      error_printf_normal();

    }

    PetscTraceBackErrorHandler(PETSC_COMM_SELF,line,fun,file,n,p,mess,ctx);

    PetscBool ismain,isunknown;

    PetscStrncmp(fun,"main",4,&ismain);
    PetscStrncmp(fun,"unknown",7,&isunknown);

    if(ismain || isunknown) {

      std::stringstream strs_version;
      strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;

      error_printf_hilight();
      (*PetscErrorPrintf)("----------MoFEM End of Error Message -------send entire error message to CMatGU <cmatgu@googlegroups.com> ----------\n");
      error_printf_normal();

    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before aborting so will not accidentally kill process 0*/
    PetscSleep(10.0);
    abort();

  }

  PetscFunctionReturn(n);
}

Core::Core(moab::Interface& _moab,MPI_Comm _comm,int _verbose):
moab(_moab),
comm(_comm),
verbose(_verbose) {

  if(!isGloballyInitialised) {
    PetscPushErrorHandler(mofem_error_handler,PETSC_NULL);
    isGloballyInitialised = true;
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,comm);

  MPI_Comm_size(comm,&sIze);
  MPI_Comm_rank(comm,&rAnk);

  if(verbose>0) {
    print_MoFem_verison(comm);
  }

  ierr = query_interface(meshsetsManagerPtr); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  ierr = getTags(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = clearMap(); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  basicEntityDataPtr = boost::shared_ptr<BasicEntityData>(new BasicEntityData(moab));
  ierr = initialiseDatabseInformationFromMesh(verbose); CHKERRABORT(PETSC_COMM_WORLD,ierr);

  // Petsc Logs
  PetscLogEventRegister("FE_preProcess",0,&USER_EVENT_preProcess);
  PetscLogEventRegister("FE_operator",0,&USER_EVENT_operator);
  PetscLogEventRegister("FE_postProcess",0,&USER_EVENT_postProcess);
  PetscLogEventRegister("MoFEMCreateMat",0,&USER_EVENT_createMat);
  PetscLogEventRegister("MoFEMBuildProblem",0,&USER_EVENT_buildProblem);

}
Core::~Core() {
}
moab::Interface& Core::get_moab() {
  return moab;
}
const moab::Interface& Core::get_moab() const {
  return moab;
}
MPI_Comm Core::get_comm() const {
  return comm;
}
BitFieldId Core::get_BitFieldId(const std::string& name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) {
    THROW_MESSAGE("field < "+name+" > not in database (top tip: check spelling)");
  }
  return (*miit)->getId();
}
std::string Core::get_BitFieldId_name(const BitFieldId id) const {
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
  return get_field_meshset(get_BitFieldId(name));
}

bool Core::check_field(const std::string &name) const {
  typedef Field_multiIndex::index<FieldName_mi_tag>::type FieldSetByName;
  const FieldSetByName &set = fIelds.get<FieldName_mi_tag>();
  FieldSetByName::iterator miit = set.find(name);
  if(miit==set.end()) return false;
  return true;
}

bool Core::check_finite_element(const std::string &name) const {
  typedef FiniteElement_multiIndex::index<FiniteElement_name_mi_tag>::type FeSetByName;
  const FeSetByName &set = finiteElements.get<FiniteElement_name_mi_tag>();
  FeSetByName::iterator miit = set.find(name);
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
  PetscErrorCode ierr;
  PetscFunctionBegin;
  refinedEntities.clear();
  refinedFiniteElements.clear();
  fIelds.clear();
  entsFields.clear();
  dofsField.clear();
  finiteElements.clear();
  entsFiniteElements.clear();
  entFEAdjacencies.clear();
  pRoblems.clear();
  sEries.clear();
  seriesSteps.clear();

  MeshsetsManager *m_manger_ptr;
  ierr = query_interface(m_manger_ptr); CHKERRQ(ierr);
  ierr = m_manger_ptr->clearMap(); CHKERRQ(ierr);

  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->clearMap(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode Core::addPrismToDatabase(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    std::pair<RefEntity_multiIndex::iterator,bool> p_ent;
    p_ent = refinedEntities.insert(boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,prism)));
    if(p_ent.second) {
      std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      p_MoFEMFiniteElement = refinedFiniteElements.insert(
	      ptrWrapperRefElement(boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*p_ent.first)))
      );
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(prism,conn,num_nodes,true); MOAB_THROW(rval);
      Range face_side3,face_side4;
      rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERRQ_MOAB(rval);
      rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERRQ_MOAB(rval);
      if(face_side3.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prism don't have side face 3");
      if(face_side4.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prims don't have side face 4");
      p_MoFEMFiniteElement.first->getSideNumberPtr(moab,*face_side3.begin());
      p_MoFEMFiniteElement.first->getSideNumberPtr(moab,*face_side4.begin());
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::getTags(int verb) {
  // PetscErrorCode ierr;
  MoABErrorCode rval;

  PetscFunctionBegin;

  const EntityHandle root_meshset = moab.get_root_set();
  if(root_meshset) {
    THROW_MESSAGE("Root meshset should be 0");
  }
  // Version
  Tag th_version;
  std::stringstream strs_version;
  strs_version << "MoFEM_version_" << MoFEM_VERSION_MAJOR << "." << MoFEM_VERSION_MINOR << "." << MoFEM_VERSION_BUILD;
  std::string version = strs_version.str();
  rval = moab.tag_get_handle(
    "_MoFEM_VERSION",
    version.size()*sizeof(char),
    MB_TYPE_OPAQUE,
    th_version,
    MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,
    NULL
  );
  if(rval==MB_ALREADY_ALLOCATED) {
    rval = MB_SUCCESS;
  } else {
    CHKERRQ_MOAB(rval);
    const char *ptr_version = version.c_str();
    rval = moab.tag_set_data(th_version,&root_meshset,1,ptr_version); CHKERRQ_MOAB(rval);
  }

  //tags saved in vtk-files
  const int def_part = -1;
  rval = moab.tag_get_handle(
    "PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part
  ); CHKERRQ_MOAB(rval);

  //Tags Ref
  EntityHandle def_handle = 0;
  rval = moab.tag_get_handle(
    "_RefParentHandle",
    1,
    MB_TYPE_HANDLE,
    th_RefParentHandle,
    MB_TAG_CREAT|MB_TAG_SPARSE,
    &def_handle
  ); CHKERRQ_MOAB(rval);
  BitRefLevel def_bit_level = 0;
  rval = moab.tag_get_handle(
    "_RefBitLevel",
    sizeof(BitRefLevel),
    MB_TYPE_OPAQUE,
    th_RefBitLevel,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,
    &def_bit_level
  ); CHKERRQ_MOAB(rval);
  BitRefLevel def_bit_level_mask = BitRefLevel().set();
  rval = moab.tag_get_handle(
    "_RefBitLevelMask",
    sizeof(BitRefLevel),
    MB_TYPE_OPAQUE,
    th_RefBitLevel_Mask,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,
    &def_bit_level_mask
  ); CHKERRQ_MOAB(rval);
  BitRefEdges def_bit_egde = 0;
  rval = moab.tag_get_handle(
    "_RefBitEdge",
    sizeof(BitRefEdges),
    MB_TYPE_OPAQUE,
    th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,
    &def_bit_egde
  ); CHKERRQ_MOAB(rval);
  const int def_type[] = {0,0};
  rval = moab.tag_get_handle(
    "_RefType",
    2,
    MB_TYPE_INTEGER,
    th_RefType,
    MB_TAG_CREAT|MB_TAG_SPARSE,
    def_type
  ); CHKERRQ_MOAB(rval);

  //Tags Field
  const unsigned long int def_id = 0;
  rval = moab.tag_get_handle(
    "_FieldId",
    sizeof(BitFieldId),
    MB_TYPE_OPAQUE,
    th_FieldId,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,
    &def_id
  ); CHKERRQ_MOAB(rval);
  FieldSpace def_space = LASTSPACE;
  rval = moab.tag_get_handle(
    "_FieldSpace",
    sizeof(FieldSpace),
    MB_TYPE_OPAQUE,
    th_FieldSpace,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,
    &def_space
  ); CHKERRQ_MOAB(rval);
  FieldApproximationBase def_base = LASTBASE;
  rval = moab.tag_get_handle(
    "_FieldBase",
    sizeof(FieldApproximationBase),
    MB_TYPE_OPAQUE,
    th_FieldBase,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,
    &def_base
  ); CHKERRQ_MOAB(rval);

  const int def_val_len = 0;
  rval = moab.tag_get_handle(
    "_FieldName",
    def_val_len,
    MB_TYPE_OPAQUE,
    th_FieldName,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
    NULL
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_FieldName_DataNamePrefix",
    def_val_len,
    MB_TYPE_OPAQUE,
    th_FieldName_DataNamePrefix,
    MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,
    NULL
  ); CHKERRQ_MOAB(rval);

  //Tags FE
  rval = moab.tag_get_handle(
    "_FEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_FEId,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_FEName",def_val_len,MB_TYPE_OPAQUE,
    th_FEName,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,NULL
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_FEIdCol",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdCol,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_FEIdRow",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdRow,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_FEIdData",sizeof(BitFieldId),MB_TYPE_OPAQUE,
    th_FEIdData,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);

  //Tags Problem
  rval = moab.tag_get_handle("_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
    th_ProblemId,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle("_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
    th_ProblemFEId,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle("_ProblemName",def_val_len,MB_TYPE_OPAQUE,
    th_ProblemName,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,NULL
  ); CHKERRQ_MOAB(rval);
  DofIdx def_nbdofs = 0;
  rval = moab.tag_get_handle(
    "_ProblemNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsRow,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_ProblemNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemNbDofsCol,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_ProblemLocalNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofRow,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_ProblemGhostNbDofsRow",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofRow,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_ProblemLocalNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemLocalNbDofCol,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);
  rval = moab.tag_get_handle(
    "_ProblemGhostNbDofsCol",sizeof(DofIdx),MB_TYPE_OPAQUE,
    th_ProblemGhostNbDofCol,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_nbdofs
  ); CHKERRQ_MOAB(rval);

  //Global Variables
  //Fields
  int def_shift = 1;
  rval = moab.tag_get_handle("_FieldShift",1,MB_TYPE_INTEGER,th_FieldShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERRQ_MOAB(rval);
  const void* tag_data[1];
  rval = moab.tag_get_by_ptr(th_FieldShift,&root_meshset,1,tag_data); CHKERRQ_MOAB(rval);
  fShift = (int*)tag_data[0];
  //FE
  rval = moab.tag_get_handle("_FEShift",1,MB_TYPE_INTEGER,th_FEShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERRQ_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_FEShift,&root_meshset,1,tag_data); CHKERRQ_MOAB(rval);
  feShift = (int*)tag_data[0];
  //Problem
  rval = moab.tag_get_handle("_ProblemShift",1,MB_TYPE_INTEGER,th_ProblemShift,MB_TAG_CREAT|MB_TAG_MESH,&def_shift);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  CHKERRQ_MOAB(rval);
  rval = moab.tag_get_by_ptr(th_ProblemShift,&root_meshset,1,tag_data); CHKERRQ_MOAB(rval);
  pShift = (int*)tag_data[0];
  //SaftyNets
  int def_bool = 0;
  rval = moab.tag_get_handle("_MoFEMBuild",1,MB_TYPE_INTEGER,th_MoFEMBuild,MB_TAG_CREAT|MB_TAG_MESH,&def_bool);
  if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
  rval = moab.tag_get_by_ptr(th_MoFEMBuild,&root_meshset,1,(const void **)&buildMoFEM); CHKERRQ_MOAB(rval);
  //Series
  rval = moab.tag_get_handle("_SeriesName",def_val_len,MB_TYPE_OPAQUE,
    th_SeriesName,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_VARLEN|MB_TAG_SPARSE,NULL
  ); CHKERRQ_MOAB(rval);

  //Meshsets with boundary conditions and material sets
  MeshsetsManager *meshsets_manager_ptr;
  ierr = query_interface(meshsets_manager_ptr); CHKERRQ(ierr);
  ierr = meshsets_manager_ptr->getTags(verb); CHKERRQ(ierr);

  //Coordinate systems
  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->getTags(verb); CHKERRQ(ierr);

  //For VTK files
  int def_elem_type = MBMAXTYPE;
  rval = moab.tag_get_handle(
    "ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type
  ); CHKERRQ_MOAB(rval);

  PetscFunctionReturn(0);
}

PetscErrorCode Core::clear_database(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = clearMap(); CHKERRQ(ierr);
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
  MeshsetsManager *m_manger_ptr;
  ierr = query_interface(m_manger_ptr); CHKERRQ(ierr);
  ierr = m_manger_ptr->initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);
  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);

  Range meshsets;
  rval = moab.get_entities_by_type(0,MBENTITYSET,meshsets,true);  CHKERRQ_MOAB(rval);
  Range::iterator mit;
  mit = meshsets.begin();
  for(;mit!=meshsets.end();mit++) {
    BitFieldId field_id;
    // Get bit id form field tag
    rval = moab.tag_get_data(th_FieldId,&*mit,1,&field_id); CHKERRQ_MOAB(rval);
    // Check if meshset if field meshset
    if(field_id!=0) {
      std::pair<Field_multiIndex::iterator,bool> p;
      try {
        EntityHandle coord_sys_id;
        rval = moab.tag_get_data(
          cs_manger_ptr->get_th_CoordSysMeshset(),&*mit,1,&coord_sys_id
        ); CHKERRQ_MOAB(rval);
        boost::shared_ptr<CoordSys> cs_ptr;
        if(coord_sys_id!=0) {
          ierr = cs_manger_ptr->getCoordSysPtr(coord_sys_id,cs_ptr); CHKERRQ(ierr);
        } else {
          ierr = cs_manger_ptr->getCoordSysPtr("UNDEFINED",cs_ptr); CHKERRQ(ierr);
        }
        p = fIelds.insert(boost::shared_ptr<Field>(new Field(moab,*mit,cs_ptr)));
        if(verb > 0) {
          std::ostringstream ss;
          ss << "read field " << **p.first << std::endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
      if((*p.first)->getSpace()==NOFIELD) {
        assert((*p.first)->meshSet == *mit);
        //add field to ref ents
        std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent;
        p_ref_ent = refinedEntities.insert(boost::shared_ptr<RefEntity>(
          new RefEntity(basicEntityDataPtr,*mit))
        );
        NOT_USED(p_ref_ent);
      } else {
        Range ents;
        rval = moab.get_entities_by_handle(*mit,ents,false); CHKERRQ_MOAB(rval);
        if(verb > 1) {
          std::ostringstream ss;
          ss << "read field ents " << ents.size() << std::endl;;
          PetscPrintf(comm,ss.str().c_str());
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent;
          p_ref_ent = refinedEntities.insert(boost::shared_ptr<
            RefEntity>(new RefEntity(basicEntityDataPtr,*eit))
          );
          try {
            boost::shared_ptr<MoFEMEntity> moabent(new MoFEMEntity(*p.first,*p_ref_ent.first));
            std::pair<MoFEMEntity_multiIndex::iterator,bool> p_ent = entsFields.insert(moabent);
            NOT_USED(p_ent);
          } catch (const std::exception& ex) {
            std::ostringstream ss;
            ss << ex.what() << std::endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
          }
        }
      }
    }
    // Check for finite elements
    BitFieldId fe_id;
    // Get bit id from fe tag
    rval = moab.tag_get_data(th_FEId,&*mit,1,&fe_id); CHKERRQ_MOAB(rval);
    //check if meshset is finite element meshset
    if(fe_id!=0) {
      std::pair<FiniteElement_multiIndex::iterator,bool> p = finiteElements.insert(
        boost::shared_ptr<FiniteElement>(new FiniteElement(moab,*mit))
      );
      if(verb > 0) {
        std::ostringstream ss;
        ss << "read finite element " << **p.first << std::endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
      NOT_USED(p);
      assert((*p.first)->meshset == *mit);
      Range ents;
      rval = moab.get_entities_by_type(*mit,MBENTITYSET,ents,false); CHKERRQ_MOAB(rval);
      rval = moab.get_entities_by_handle(*mit,ents,true); CHKERRQ_MOAB(rval);
      Range::iterator eit = ents.begin();
      for(;eit!=ents.end();eit++) {
        std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent;
        p_ref_ent = refinedEntities.insert(boost::shared_ptr<RefEntity>(
          new RefEntity(basicEntityDataPtr,*eit))
        );
        std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
        try {
          switch (moab.type_from_handle(*eit)) {
            case MBVERTEX:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_VERTEX(moab,*p_ref_ent.first)))
            );
            break;
            case MBEDGE:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_EDGE(moab,*p_ref_ent.first)))
            );
            break;
            case MBTRI:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TRI(moab,*p_ref_ent.first)))
            );
            break;
            case MBTET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TET(moab,*p_ref_ent.first)))
            );
            break;
            case MBPRISM:
            ierr = addPrismToDatabase(*eit,verb); CHKERRQ(ierr);
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*p_ref_ent.first)))
            );
            break;
            case MBENTITYSET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_MESHSET(moab,*p_ref_ent.first)))
            );
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
    rval = moab.tag_get_data(th_ProblemId,&*mit,1,&problem_id); CHKERRQ_MOAB(rval);
    //check if meshset if problem meshset
    if(problem_id!=0) {
      std::pair<MoFEMProblem_multiIndex::iterator,bool> p = pRoblems.insert(MoFEMProblem(moab,*mit));
      if(verb > 0) {
        std::ostringstream ss;
        ss << "read problem " << *p.first << std::endl;;
        PetscPrintf(comm,ss.str().c_str());
      }
    }
    //check if meshset is Series meshset
    {
      const void* tag_name_data;
      int tag_name_size;
      rval = moab.tag_get_by_ptr(th_SeriesName,&*mit,1,(const void **)&tag_name_data,&tag_name_size);
      if(rval == MB_SUCCESS) {
        std::pair<Series_multiIndex::iterator,bool> p = sEries.insert(MoFEMSeries(moab,*mit));
        if(verb > 0) {
          std::ostringstream ss;
          ss << "read series " << *p.first << std::endl;
          PetscPrintf(comm,ss.str().c_str());
        }
      }
    }
  }
  //build ref entities meshset
  for(int dd = 0;dd<=3;dd++) {
    Range ents;
    rval = moab.get_entities_by_dimension(0,dd,ents,false); CHKERRQ_MOAB(rval);
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
      boost::shared_ptr<RefEntity> mofem_ent(
        new RefEntity(basicEntityDataPtr,*eit)
      );
      BitRefLevel bit = mofem_ent->getBitRefLevel();
      if(bit.none()) {
        continue;
      }
      std::pair<RefEntity_multiIndex::iterator,bool> p;
      p = refinedEntities.insert(mofem_ent);
    }
  }
  //build series steps
  for(Series_multiIndex::iterator sit = sEries.begin();sit!=sEries.end();sit++) {
    int nb_steps;
    ierr = sit->get_nb_steps(moab,nb_steps); CHKERRQ(ierr);
    int ss = 0;
    for(;ss<nb_steps;ss++) {
      std::pair<SeriesStep_multiIndex::iterator,bool> p = seriesSteps.insert(MoFEMSeriesStep(moab,&*sit,ss));
      if(verb > 0) {
        std::ostringstream ss;
        ss << "add series step " << *p.first << std::endl;
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

// PetscErrorCode Core::add_coordinate_system(const int cs_dim[],const std::string name) {
//   CoordSystemsManager *cs_manger_ptr;
//   PetscFunctionBegin;
//   ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
//   ierr = cs_manger_ptr->addCoordinateSystem(cs_dim,name); CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }
//
// PetscErrorCode Core::set_field_coordinate_system(const std::string field_name,const std::string cs_name) {
//   CoordSystemsManager *cs_manger_ptr;
//   PetscFunctionBegin;
//   ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
//   ierr = cs_manger_ptr->setFieldCoordinateSystem(field_name,cs_name); CHKERRQ(ierr);
//   PetscFunctionReturn(0);
// }

// cubit meshsets

PetscErrorCode Core::print_cubit_displacement_set() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printDisplacementSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_pressure_set() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printPressureSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_force_set() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printForceSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_temperature() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printTemperatureSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_heat_flux_set() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printHeatFluxSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_materials_set() const {
  PetscErrorCode ierr;
  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printMaterialsSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

bool Core::check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type) {
  return meshsetsManagerPtr->checkMeshset(ms_id,cubit_bc_type);
}

PetscErrorCode Core::add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
  return meshsetsManagerPtr->addMeshset(cubit_bc_type,ms_id,name);
}

PetscErrorCode Core::set_cubit_msId_attribites(
  const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name
) {
  return meshsetsManagerPtr->setAttribites(cubit_bc_type,ms_id,attributes,name);
}
PetscErrorCode Core::set_cubit_msId_attribites_data_structure(
  const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name
) {
  return meshsetsManagerPtr->setAttribitesByDataStructure(cubit_bc_type,ms_id,data,name);
}
PetscErrorCode Core::set_cubit_msId_bc_data_structure(
  const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
) {
  return meshsetsManagerPtr->setBcData(cubit_bc_type,ms_id,data);
}
PetscErrorCode Core::delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id) {
  return meshsetsManagerPtr->deleteMeshset(cubit_bc_type,ms_id);
}
PetscErrorCode Core::get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
  return meshsetsManagerPtr->getCubitMeshsetPtr(ms_id,cubit_bc_type,cubit_meshset_ptr);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(
  const int msId,const CubitBCType cubit_bc_type,const int dimension,Range &entities,const bool recursive
) {
  return meshsetsManagerPtr->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),dimension,entities,recursive);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,Range &entities,const bool recursive) {
  return meshsetsManagerPtr->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),entities,recursive);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(
  const int ms_id,const unsigned int cubit_bc_type,const int dimension,Range &entities,const bool recursive
) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),dimension,entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int ms_id,const unsigned int cubit_bc_type,
  Range &entities,const bool recursive) {
  PetscFunctionBegin;
  ierr = get_cubit_msId_entities_by_dimension(ms_id,CubitBCType(cubit_bc_type),entities,recursive); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_cubit_msId_meshset(const int ms_id,const unsigned int cubit_bc_type,EntityHandle &meshset) {
  return meshsetsManagerPtr->getMeshset(ms_id,cubit_bc_type,meshset);
}

PetscErrorCode Core::get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) {
  return meshsetsManagerPtr->getMeshsetsByType(cubit_bc_type,meshsets);
}

PetscErrorCode Core::get_fields(const Field_multiIndex **fields_ptr) const {
  PetscFunctionBegin;
  *fields_ptr = &fIelds;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const {
  PetscFunctionBegin;
  *refined_entities_ptr = &refinedEntities;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_ref_finite_elements(const RefElement_multiIndex **refined_finite_elements_ptr) const {
  PetscFunctionBegin;
  *refined_finite_elements_ptr = &refinedFiniteElements;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_problem(const std::string &problem_name,const MoFEMProblem **problem_ptr) const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems.find(problem_name);
  if(p_miit == problems.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "problem < %s > not found, (top tip: check spelling)",problem_name.c_str()
    );
  }
  *problem_ptr = &*p_miit;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_field_ents(const MoFEMEntity_multiIndex **field_ents) const {
  PetscFunctionBegin;
  *field_ents = &entsFields;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::get_dofs(const DofEntity_multiIndex **dofs_ptr) const {
  PetscFunctionBegin;
  *dofs_ptr = &dofsField;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::seed_ref_level_2D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range ents2d;
  rval = moab.get_entities_by_dimension(meshset,2,ents2d,false); CHKERRQ_MOAB(rval);
  ierr = seed_ref_level(ents2d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::seed_ref_level(const Range &ents,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  Range seeded_ents;
  try {
    if(verb > 1) {
      PetscSynchronizedPrintf(comm,"nb. entities for seed %d\n",ents.size());
    }
    Range::iterator tit = ents.begin();
    for(;tit!=ents.end();tit++) {
      boost::shared_ptr<RefEntity> ref_ent(new RefEntity(basicEntityDataPtr,*tit));
      std::bitset<8> ent_pstat(ref_ent->getPStatus());
      ent_pstat.flip(0);
      std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(ref_ent);
      if(debug > 0) {
        ierr = test_moab(moab,*tit); CHKERRQ(ierr);
      }
      bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
      if(!success) {
        SETERRQ(
          comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful"
        );
      }
      if(verb>2) {
        std::ostringstream ss;
        ss << **p_ent.first;
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
      std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      switch((*p_ent.first)->getEntType()) {
        case MBVERTEX:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_VERTEX(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBEDGE:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_EDGE(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBTRI:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TRI(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBTET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_TET(moab,*p_ent.first)))
        );
        seeded_ents.insert(*tit);
        break;
        case MBPRISM:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_PRISM(moab,*p_ent.first)))
        );
        break;
        case MBENTITYSET:
        p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
          boost::shared_ptr<RefElement>(new RefElement_MESHSET(moab,*p_ent.first)))
        );
        break;
        default:
        SETERRQ(comm,MOFEM_NOT_IMPLEMENTED,"not implemented");
      }
      if(verb>3) {
        std::ostringstream ss;
        ss << *(p_MoFEMFiniteElement.first->getRefElement());
        PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
      }
    }
    if(!seeded_ents.empty()) {
      int dim = moab.dimension_from_handle(seeded_ents[0]);
      for(int dd = 0;dd<dim;dd++) {
        Range ents;
        rval = moab.get_adjacencies(seeded_ents,dd,true,ents,moab::Interface::UNION); CHKERRQ_MOAB(rval);
        if(dd == 2) {
          // currently only works with triangles
          ents = ents.subset_by_type(MBTRI);
        }
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
            boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,*eit))
          );
          bool success = refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
          if(!success) SETERRQ(comm,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
          if(verb>2) {
            std::ostringstream ss;
            ss << *(p_ent.first);
            PetscSynchronizedPrintf(comm,"%s\n",ss.str().c_str());
          }
        }
      }
    }
    if(verb>2) {
      PetscSynchronizedPrintf(comm,"\n");
      PetscSynchronizedFlush(comm,PETSC_STDOUT);
    }
  } catch (MoFEMException const &e) {
    SETERRQ(comm,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  Range ents3d;
  rval = moab.get_entities_by_dimension(meshset,3,ents3d,false); CHKERRQ_MOAB(rval);
  ierr = seed_ref_level(ents3d,bit,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  std::pair<RefEntity_multiIndex::iterator,bool> p_ent = refinedEntities.insert(
    boost::shared_ptr<RefEntity>(new RefEntity(basicEntityDataPtr,meshset))
  );
  refinedEntities.modify(p_ent.first,RefEntity_change_add_bit(bit));
  ptrWrapperRefElement pack_fe(
    boost::shared_ptr<RefElement>(new RefElement_MESHSET(moab,*p_ent.first))
  );
  std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement = refinedFiniteElements.insert(pack_fe);
  if(verbose > 0) {
    std::ostringstream ss;
    ss << "add meshset as ref_ent " << *(p_MoFEMFiniteElement.first->getRefElement()) << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}


}
