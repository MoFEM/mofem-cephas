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
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <UnknownInterface.hpp>
#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

// Interfaces
#include <ProblemsManager.hpp>
#include <Simple.hpp>
#include <ISManager.hpp>
#include <BitRefManager.hpp>
#include <VecManager.hpp>
#include <FieldBlas.hpp>
#include <MeshRefinement.hpp>
#include <SeriesRecorder.hpp>
#include <PrismInterface.hpp>
#include <CutMeshInterface.hpp>
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

PetscErrorCode Core::queryInterface(const MOFEMuuid& uuid,UnknownInterface** iface) {
  PetscFunctionBegin;
  *iface = NULL;
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

  ptr = NULL;

  // TetGen
  #ifdef WITH_TETGEN
  if(type == typeid(TetGenInterface)) {
    if(iFaces.find(IDD_MOFEMTetGegInterface.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMTetGegInterface.uUId.to_ulong();
      iFaces.insert(uid,new TetGenInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMTetGegInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  // MedInterface
  #ifdef WITH_MED
  if(type == typeid(MedInterface)) {
    if(iFaces.find(IDD_MOFEMMedInterface.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMMedInterface.uUId.to_ulong();
      iFaces.insert(uid,new MedInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMMedInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }
  #endif

  // Problems manager
  if(type == typeid(ProblemsManager)) {
    if(iFaces.find(IDD_MOFEMProblemsManager.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMProblemsManager.uUId.to_ulong();
      iFaces.insert(uid,new ProblemsManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMProblemsManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  // Simple interface
  if(type == typeid(Simple)) {
    if(iFaces.find(IDD_MOFEMSimple.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMSimple.uUId.to_ulong();
      iFaces.insert(uid,new Simple(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMSimple.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(ISManager)) {
    if(iFaces.find(IDD_MOFEMISManager.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMISManager.uUId.to_ulong();
      iFaces.insert(uid,new ISManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMISManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(VecManager)) {
    if(iFaces.find(IDD_MOFEMVEC.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMVEC.uUId.to_ulong();
      iFaces.insert(uid,new VecManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMVEC.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(FieldBlas)) {
    if(iFaces.find(IDD_MOFEMFieldBlas.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMFieldBlas.uUId.to_ulong();
      iFaces.insert(uid,new FieldBlas(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMFieldBlas.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(BitRefManager)) {
    if(iFaces.find(IDD_MOFEMBitRefManager.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMBitRefManager.uUId.to_ulong();
      iFaces.insert(uid,new BitRefManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMBitRefManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Meshsets manager
  if(type == typeid(MeshsetsManager)) {
    if(iFaces.find(IDD_MOFEMMeshsetsManager.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMMeshsetsManager.uUId.to_ulong();
      iFaces.insert(uid,new MeshsetsManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMMeshsetsManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Coordinate systems manager
  if(type == typeid(CoordSystemsManager)) {
    if(iFaces.find(IDD_MOFEMCoordsSystemsManager.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMCoordsSystemsManager.uUId.to_ulong();
      iFaces.insert(uid,new CoordSystemsManager(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMCoordsSystemsManager.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Node merger
  if(type == typeid(NodeMergerInterface)) {
    if(iFaces.find(IDD_MOFEMNodeMerger.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMNodeMerger.uUId.to_ulong();
      iFaces.insert(uid,new NodeMergerInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMNodeMerger.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //BitLevelCoupler
  if(type == typeid(BitLevelCoupler)) {
    if(iFaces.find(IDD_MOFEMBitLevelCoupler.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMBitLevelCoupler.uUId.to_ulong();
      iFaces.insert(uid,new BitLevelCoupler(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMBitLevelCoupler.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  //Create prism elements from surface Elements
  if(type == typeid(PrismsFromSurfaceInterface)) {
    if(iFaces.find(IDD_MOFEMPrismsFromSurface.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMPrismsFromSurface.uUId.to_ulong();
      iFaces.insert(uid,new PrismsFromSurfaceInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMPrismsFromSurface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(MeshRefinement)) {
    if(iFaces.find(IDD_MOFEMMeshRefine.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMMeshRefine.uUId.to_ulong();
      iFaces.insert(uid,new MeshRefinement(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMMeshRefine.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(PrismInterface)) {
    if(iFaces.find(IDD_MOFEMPrismInterface.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMPrismInterface.uUId.to_ulong();
      iFaces.insert(uid,new PrismInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMPrismInterface.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(CutMeshInterface)) {
    if(iFaces.find(IDD_MOFEMCutMesh.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMCutMesh.uUId.to_ulong();
      iFaces.insert(uid,new CutMeshInterface(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMCutMesh.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  if(type == typeid(SeriesRecorder)) {
    if(iFaces.find(IDD_MOFEMSeriesRecorder.uUId.to_ulong()) == iFaces.end()) {
      unsigned long int uid = IDD_MOFEMSeriesRecorder.uUId.to_ulong();
      iFaces.insert(uid,new SeriesRecorder(*this));
    }
    ptr = &iFaces.at(IDD_MOFEMSeriesRecorder.uUId.to_ulong());
    PetscFunctionReturn(0);
  }

  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");

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

static PetscErrorCode mofem_error_handler(
  MPI_Comm comm,int line,const char *fun,const char *file,PetscErrorCode n,PetscErrorType p,const char *mess,void *ctx
) {
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
      (*PetscErrorPrintf)("----------MoFEM End of Error Message -------send entire error message to mofem-group@googlegroups.com ----------\n");
      error_printf_normal();

    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before aborting so will not accidentally kill process 0*/
    PetscSleep(10.0);
    abort();

  }

  PetscFunctionReturn(n);
}

PetscErrorCode print_verison() {
  //
  PetscFunctionBegin;
  PetscFunctionReturn(0);
}

Core::Core(moab::Interface& moab,MPI_Comm comm,int verbose):
moab(moab),
verbose(verbose),
cOmm(0) {
  if(!isGloballyInitialised) {
    PetscPushErrorHandler(mofem_error_handler,PETSC_NULL);
    isGloballyInitialised = true;
  }
  // Duplicate petsc communicator
  ierr = PetscCommDuplicate(comm,&cOmm,NULL); CHKERRABORT(comm,ierr);
  MPI_Comm_size(cOmm,&sIze);
  MPI_Comm_rank(cOmm,&rAnk);
  // CHeck if moab has set communicator if not set communicator interbally
  ParallelComm* pComm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pComm == NULL) {
    pComm =  new ParallelComm(&moab,comm);
  }
  ierr = getTags(); CHKERRABORT(cOmm,ierr);
  ierr = clearMap(); CHKERRABORT(cOmm,ierr);
  basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  ierr = initialiseDatabseInformationFromMesh(verbose); CHKERRABORT(cOmm,ierr);
  // Petsc Logs
   PetscLogEventRegister("FE_preProcess",0,&USER_EVENT_preProcess);
   PetscLogEventRegister("FE_operator",0,&USER_EVENT_operator);
   PetscLogEventRegister("FE_postProcess",0,&USER_EVENT_postProcess);
   PetscLogEventRegister("MoFEMCreateMat",0,&USER_EVENT_createMat);
   PetscLogEventRegister("MoFEMBuildProblem",0,&USER_EVENT_buildProblem);

   // Print version
  if(verbose>0) {
    ierr = PetscPrintf(
      cOmm,
      "lib version %d.%d.%d\n",
      MoFEM_VERSION_MAJOR,MoFEM_VERSION_MINOR,MoFEM_VERSION_BUILD
    ); CHKERRABORT(cOmm,ierr);
    ierr = PetscPrintf(cOmm,"git commit id %s\n",GIT_SHA1_NAME); CHKERRABORT(cOmm,ierr);
  }

}
Core::~Core() {
  int flg;
  MPI_Finalized(&flg);
  // Destroy interfaces
  iFaces.clear();
  // Pop error hanlder
  if(isGloballyInitialised) {
    if(!flg) {
      PetscPopErrorHandler();
    }
    isGloballyInitialised = false;
  }
  // Destroy communictaor
  if(!flg) {
    ierr = PetscCommDestroy(&cOmm); CHKERRABORT(cOmm,ierr);
  }
}

BitFieldId Core::getFieldShift() {
  if(*fShift >= BITFIELDID_SIZE) {
    char msg[] = "number of fields exceeded";
    PetscTraceBackErrorHandler(
      cOmm,
      __LINE__,PETSC_FUNCTION_NAME,__FILE__,
      MOFEM_DATA_INCONSISTENCY,PETSC_ERROR_INITIAL,msg,PETSC_NULL);
    PetscMPIAbortErrorHandler(cOmm,
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

  // Cleaning databases in iterfaces
  SeriesRecorder *series_recorder_ptr;
  ierr = query_interface(series_recorder_ptr); CHKERRQ(ierr);
  ierr = series_recorder_ptr->clearMap(); CHKERRQ(ierr);
  MeshsetsManager *m_manger_ptr;
  ierr = query_interface(m_manger_ptr); CHKERRQ(ierr);
  ierr = m_manger_ptr->clearMap(); CHKERRQ(ierr);
  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->clearMap(); CHKERRQ(ierr);

  // Cleaning databases
  refinedEntities.clear();
  refinedFiniteElements.clear();
  fIelds.clear();
  entsFields.clear();
  dofsField.clear();
  finiteElements.clear();
  entsFiniteElements.clear();
  entFEAdjacencies.clear();
  pRoblems.clear();


  PetscFunctionReturn(0);
}

PetscErrorCode Core::addPrismToDatabase(const EntityHandle prism,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  try {
    std::pair<RefEntity_multiIndex::iterator,bool> p_ent;
    p_ent = refinedEntities.insert(boost::make_shared<RefEntity>(basicEntityDataPtr,prism));
    if(p_ent.second) {
      std::pair<RefElement_multiIndex::iterator,bool> p_MoFEMFiniteElement;
      p_MoFEMFiniteElement = refinedFiniteElements.insert(
	      ptrWrapperRefElement(boost::shared_ptr<RefElement>(new RefElement_PRISM(*p_ent.first)))
      );
      int num_nodes;
      const EntityHandle* conn;
      rval = moab.get_connectivity(prism,conn,num_nodes,true); MOAB_THROW(rval);
      Range face_side3,face_side4;
      rval = moab.get_adjacencies(conn,3,2,false,face_side3); CHKERRQ_MOAB(rval);
      rval = moab.get_adjacencies(&conn[3],3,2,false,face_side4); CHKERRQ_MOAB(rval);
      if(face_side3.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prism don't have side face 3");
      if(face_side4.size()!=1) SETERRQ(PETSC_COMM_SELF,1,"prims don't have side face 4");
      p_MoFEMFiniteElement.first->getSideNumberPtr(*face_side3.begin());
      p_MoFEMFiniteElement.first->getSideNumberPtr(*face_side4.begin());
    }
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::getTags(int verb) {
  //


  PetscFunctionBegin;

  const EntityHandle root_meshset = moab.get_root_set();
  if(root_meshset) {
    THROW_MESSAGE("Root meshset should be 0");
  }

  // Set version
  {
    Version version;
    ierr = getFileVersion(moab,version); CHKERRQ(ierr);
    // PetscPrintf(
    //   comm,"file version %d.%d.%d\n",
    //   version.majorVersion,version.minorVersion,version.buildVersion
    // );
  }

  // Global Variables
  {
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
  }

  // Tags saved in vtk-files
  {
    const int def_part = -1;
    rval = moab.tag_get_handle(
      "PARTITION",1,MB_TYPE_INTEGER,th_Part,MB_TAG_CREAT|MB_TAG_SPARSE,&def_part
    ); CHKERRQ_MOAB(rval);
    int def_elem_type = MBMAXTYPE;
    rval = moab.tag_get_handle(
      "ElemType",1,MB_TYPE_INTEGER,th_ElemType,MB_TAG_CREAT|MB_TAG_SPARSE,&def_elem_type
    ); CHKERRQ_MOAB(rval);
  }

  // Tags Ref
  {
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
  }

  // Tags Field
  {
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
  }

  //Tags FE
  {
    const unsigned long int def_id = 0;
    const int def_val_len = 0;
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
  }

  //Tags Problem
  {
    const unsigned long int def_id = 0;
    const int def_val_len = 0;
    rval = moab.tag_get_handle(
      "_ProblemId",sizeof(BitProblemId),MB_TYPE_OPAQUE,
      th_ProblemId,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_get_handle(
      "_ProblemFEId",sizeof(BitFEId),MB_TYPE_OPAQUE,
      th_ProblemFEId,MB_TAG_CREAT|MB_TAG_BYTES|MB_TAG_SPARSE,&def_id
    ); CHKERRQ_MOAB(rval);
    rval = moab.tag_get_handle(
      "_ProblemName",def_val_len,MB_TYPE_OPAQUE,
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
  }

  //Meshsets with boundary conditions and material sets
  MeshsetsManager *meshsets_manager_ptr;
  ierr = query_interface(meshsets_manager_ptr); CHKERRQ(ierr);
  ierr = meshsets_manager_ptr->getTags(verb); CHKERRQ(ierr);

  // Series recorder
  SeriesRecorder *series_recorder_ptr;
  ierr = query_interface(series_recorder_ptr); CHKERRQ(ierr);
  ierr = series_recorder_ptr->getTags(verb); CHKERRQ(ierr);

  //Coordinate systems
  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->getTags(verb); CHKERRQ(ierr);


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

  CoordSystemsManager *cs_manger_ptr;
  ierr = query_interface(cs_manger_ptr); CHKERRQ(ierr);

  // Initialize coordinate systems
  ierr = cs_manger_ptr->initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);

  // Initialize database
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
        const char *cs_name;
        int cs_name_size;
        rval = moab.tag_get_by_ptr(
          cs_manger_ptr->get_th_CoordSysName(),&*mit,1,(const void **)&cs_name,&cs_name_size
        );
        boost::shared_ptr<CoordSys> cs_ptr;
        if(rval == MB_SUCCESS && cs_name_size) {
          ierr = cs_manger_ptr->getCoordSysPtr(std::string(cs_name,cs_name_size),cs_ptr); CHKERRQ(ierr);
        } else {
          ierr = cs_manger_ptr->getCoordSysPtr("UNDEFINED",cs_ptr); CHKERRQ(ierr);
        }
        p = fIelds.insert(boost::shared_ptr<Field>(new Field(moab,*mit,cs_ptr)));
        if(verb > 0) {
          std::ostringstream ss;
          ss << "read field " << **p.first << std::endl;;
          PetscPrintf(cOmm,ss.str().c_str());
        }
      } catch (MoFEMException const &e) {
        SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
      }
      if((*p.first)->getSpace()==NOFIELD) {
        assert((*p.first)->meshSet == *mit);
        //add field to ref ents
        std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent;
        p_ref_ent = refinedEntities.insert(
          boost::make_shared<RefEntity>(basicEntityDataPtr,*mit)
        );
        NOT_USED(p_ref_ent);
      } else {
        Range ents;
        rval = moab.get_entities_by_handle(*mit,ents,false); CHKERRQ_MOAB(rval);
        if(verb > 1) {
          std::ostringstream ss;
          ss << "read field ents " << ents.size() << std::endl;;
          PetscPrintf(cOmm,ss.str().c_str());
        }
        boost::shared_ptr<std::vector<FieldEntity> > ents_array =
        boost::make_shared<std::vector<FieldEntity> >(std::vector<FieldEntity>());
        // Add sequence to field data structure. Note that entities are allocated
        // once into vector. This vector is passed into sequence as a weak_ptr.
        // Vector is destroyed at the point last entity inside that vector is
        // destroyed.
        p.first->get()->getEntSeqenceContainer()->push_back(ents_array);
        ents_array->reserve(ents.size());
        std::vector<boost::shared_ptr<FieldEntity> > ents_shared_array;
        ents_shared_array.reserve(ents.size());
        Range::iterator eit = ents.begin();
        for(;eit!=ents.end();eit++) {
          std::pair<RefEntity_multiIndex::iterator,bool> p_ref_ent;
          p_ref_ent = refinedEntities.insert(
            boost::make_shared<RefEntity>(basicEntityDataPtr,*eit)
          );
          try {
            // NOTE: This will work with newer compiler only, use push_back for back compatibility.
            // ents_array->emplace_back(*p.first,*p_ref_ent.first);
            // ents_shared_array.emplace_back(ents_array,&ents_array->back());
            ents_array->push_back(FieldEntity(*p.first,*p_ref_ent.first));
            ents_shared_array.push_back(
              boost::shared_ptr<FieldEntity>(ents_array,&ents_array->back())
            );
          } catch (const std::exception& ex) {
            std::ostringstream ss;
            ss << ex.what() << std::endl;
            SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
          }
        }
        entsFields.insert(ents_shared_array.begin(),ents_shared_array.end());
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
        PetscPrintf(cOmm,ss.str().c_str());
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
              boost::shared_ptr<RefElement>(new RefElement_VERTEX(*p_ref_ent.first)))
            );
            break;
            case MBEDGE:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_EDGE(*p_ref_ent.first)))
            );
            break;
            case MBTRI:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TRI(*p_ref_ent.first)))
            );
            break;
            case MBTET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_TET(*p_ref_ent.first)))
            );
            break;
            case MBPRISM:
            ierr = addPrismToDatabase(*eit,verb); CHKERRQ(ierr);
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_PRISM(*p_ref_ent.first)))
            );
            break;
            case MBENTITYSET:
            p_MoFEMFiniteElement = refinedFiniteElements.insert(ptrWrapperRefElement(
              boost::shared_ptr<RefElement>(new RefElement_MESHSET(*p_ref_ent.first)))
            );
            break;
            default:
            SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"Only finite elements of type MBTET, MBPRISM and MBENTITYSET are implemented");
          }
          if(p_MoFEMFiniteElement.second) {
            //PetscPrintf(cOmm,"Warrning: this entity should be already in refined finite elements database");
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
      std::pair<Problem_multiIndex::iterator,bool> p = pRoblems.insert(Problem(moab,*mit));
      if(verb > 0) {
        std::ostringstream ss;
        ss << "read problem " << *p.first << std::endl;;
        PetscPrintf(cOmm,ss.str().c_str());
      }
    }
    //check if meshset is Series meshset
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

  if(verb > 2) {
    list_fields();
    list_finite_elements();
    list_problem();
  }

  // Initialize interfaces
  MeshsetsManager *m_manger_ptr;
  ierr = query_interface(m_manger_ptr); CHKERRQ(ierr);
  ierr = m_manger_ptr->initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);
  SeriesRecorder *series_recorder_ptr;
  ierr = query_interface(series_recorder_ptr); CHKERRQ(ierr);
  ierr = series_recorder_ptr->initialiseDatabseInformationFromMesh(verb); CHKERRQ(ierr);

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

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printDisplacementSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_pressure_set() const {

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printPressureSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_force_set() const {

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printForceSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_temperature() const {

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printTemperatureSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_heat_flux_set() const {

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printHeatFluxSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::print_cubit_materials_set() const {

  MeshsetsManager *meshsets_manager;
  PetscFunctionBegin;
  ierr = query_interface(meshsets_manager); CHKERRQ(ierr);
  ierr = meshsets_manager->printMaterialsSet(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

bool Core::check_msId_meshset(const int ms_id,const CubitBCType cubit_bc_type) {
  return get_meshsets_manager_ptr()->checkMeshset(ms_id,cubit_bc_type);
}

PetscErrorCode Core::add_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id,const std::string name) {
  return get_meshsets_manager_ptr()->addMeshset(cubit_bc_type,ms_id,name);
}

PetscErrorCode Core::set_cubit_msId_attribites(
  const CubitBCType cubit_bc_type,const int ms_id,const std::vector<double> &attributes,const std::string name
) {
  return get_meshsets_manager_ptr()->setAttribites(cubit_bc_type,ms_id,attributes,name);
}
PetscErrorCode Core::set_cubit_msId_attribites_data_structure(
  const CubitBCType cubit_bc_type,const int ms_id,const GenericAttributeData &data,const std::string name
) {
  return get_meshsets_manager_ptr()->setAttribitesByDataStructure(cubit_bc_type,ms_id,data,name);
}
PetscErrorCode Core::set_cubit_msId_bc_data_structure(
  const CubitBCType cubit_bc_type,const int ms_id,const GenericCubitBcData &data
) {
  return get_meshsets_manager_ptr()->setBcData(cubit_bc_type,ms_id,data);
}
PetscErrorCode Core::delete_cubit_msId(const CubitBCType cubit_bc_type,const int ms_id) {
  return get_meshsets_manager_ptr()->deleteMeshset(cubit_bc_type,ms_id);
}
PetscErrorCode Core::get_cubit_msId(const int ms_id,const CubitBCType cubit_bc_type,const CubitMeshSets **cubit_meshset_ptr) {
  return get_meshsets_manager_ptr()->getCubitMeshsetPtr(ms_id,cubit_bc_type,cubit_meshset_ptr);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(
  const int msId,const CubitBCType cubit_bc_type,const int dimension,Range &entities,const bool recursive
) {
  return get_meshsets_manager_ptr()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),dimension,entities,recursive);
}
PetscErrorCode Core::get_cubit_msId_entities_by_dimension(const int msId,const CubitBCType cubit_bc_type,Range &entities,const bool recursive) {
  return get_meshsets_manager_ptr()->getEntitiesByDimension(msId,cubit_bc_type.to_ulong(),entities,recursive);
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
  return get_meshsets_manager_ptr()->getMeshset(ms_id,cubit_bc_type,meshset);
}

PetscErrorCode Core::get_cubit_meshsets(const unsigned int cubit_bc_type,Range &meshsets) {
  return get_meshsets_manager_ptr()->getMeshsetsByType(cubit_bc_type,meshsets);
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

PetscErrorCode Core::get_problem(const std::string &problem_name,const Problem **problem_ptr) const {
  PetscFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
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

PetscErrorCode Core::get_problems(const Problem_multiIndex **problems_ptr) const {
  PetscFunctionBegin;
  *problems_ptr = &pRoblems;
  PetscFunctionReturn(0);
}

PetscErrorCode Core::get_field_ents(const FieldEntity_multiIndex **field_ents) const {
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
  return BitRefManager(*this).setBitRefLevelByDim(meshset,2,bit,verb);
}

PetscErrorCode Core::seed_ref_level(const Range &ents,const BitRefLevel &bit,const bool only_tets,int verb) {
  return BitRefManager(*this).setBitRefLevel(ents,bit,only_tets,verb);
}

PetscErrorCode Core::seed_ref_level_3D(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  return BitRefManager(*this).setBitRefLevelByDim(meshset,3,bit,verb);
}

PetscErrorCode Core::seed_ref_level_MESHSET(const EntityHandle meshset,const BitRefLevel &bit,int verb) {
  return BitRefManager(*this).setBitLevelToMeshset(meshset,bit,verb);
}

MeshsetsManager* Core::get_meshsets_manager_ptr() {
  MeshsetsManager* meshsets_manager_ptr;
  query_interface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

const MeshsetsManager* Core::get_meshsets_manager_ptr() const {
  MeshsetsManager* meshsets_manager_ptr;
  query_interface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

}
