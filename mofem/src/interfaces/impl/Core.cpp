/** \file Core.cpp
 * \brief Multi-index containers, data structures and other low-level functions
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

extern "C" {
  void macro_is_depracted_using_deprecated_function() {}
}

namespace MoFEM {

PetscErrorCode Core::query_interface(const MOFEMuuid& uuid,UnknownInterface** iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if(uuid == IDD_MOFEMCoreInterface) {
    *iface = static_cast<CoreInterface*>(const_cast<Core*>(this));
    MoFEMFunctionReturnHot(0);
  } else if(uuid == IDD_MOFEMDeprecatedCoreInterface) {
    *iface = static_cast<DeprecatedCoreInterface*>(const_cast<Core*>(this));
    MoFEMFunctionReturnHot(0);
  }
  // Get sub-interface
  unsigned long int id = uuid.uUId.to_ulong();
  boost::ptr_map<unsigned long,UnknownInterface>::iterator it;
  it = iFaces.find(id);
  if(it != iFaces.end()) {
    *iface = it->second;
    MoFEMFunctionReturnHot(0);
  }
  *iface = NULL; 
  SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"unknown interface");
  MoFEMFunctionReturnHot(0);
}

bool Core::isGloballyInitialised = false;

static void error_printf_highlight(void) {
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
  MoFEMFunctionBeginHot;

  int rank = 0;
  if (comm != PETSC_COMM_SELF) MPI_Comm_rank(comm,&rank);

  if(!rank) {

    if(p == PETSC_ERROR_INITIAL) {
      error_printf_highlight();
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

      error_printf_highlight();
      (*PetscErrorPrintf)("----------MoFEM End of Error Message -------send entire error message to mofem-group@googlegroups.com ----------\n");
      error_printf_normal();

    }

  } else {

    /* do not print error messages since process 0 will print them, sleep before
     * aborting so will not accidentally kill process 0*/
    PetscSleep(10.0);
    abort();

  }

  MoFEMFunctionReturnHot(n);
}

template<class IFACE>
PetscErrorCode Core::regSubInterface(const MOFEMuuid& uid) {
  MoFEMFunctionBeginHot;
  ierr = registerInterface<IFACE>(uid); CHKERRQ(ierr);
  unsigned long int id = uid.uUId.to_ulong();
  iFaces.insert(id,new IFACE(*this));
  MoFEMFunctionReturnHot(0);
}

Core::Core(moab::Interface& moab,MPI_Comm comm,int verbose):
moab(moab),
cOmm(0),
verbose(verbose) {

  if(!isGloballyInitialised) {
    PetscPushErrorHandler(mofem_error_handler,PETSC_NULL);
    isGloballyInitialised = true;
  }

  // Register interfaces for this implementation
  ierr = registerInterface<UnknownInterface>(IDD_MOFEMUnknown); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = registerInterface<CoreInterface>(IDD_MOFEMCoreInterface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = registerInterface<DeprecatedCoreInterface>(IDD_MOFEMDeprecatedCoreInterface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  // Register sub interfaces
  ierr = regSubInterface<Simple>(IDD_MOFEMSimple); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<ProblemsManager>(IDD_MOFEMProblemsManager); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<ISManager>(IDD_MOFEMISManager); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<VecManager>(IDD_MOFEMVEC); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<FieldBlas>(IDD_MOFEMFieldBlas); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<BitRefManager>(IDD_MOFEMBitRefManager); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<Tools>(IDD_MOFEMTools); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<MeshsetsManager>(IDD_MOFEMMeshsetsManager); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<CoordSystemsManager>(IDD_MOFEMCoordsSystemsManager); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<NodeMergerInterface>(IDD_MOFEMNodeMerger); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<BitLevelCoupler>(IDD_MOFEMBitLevelCoupler); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<PrismsFromSurfaceInterface>(IDD_MOFEMPrismsFromSurface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<MeshRefinement>(IDD_MOFEMMeshRefine); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<PrismInterface>(IDD_MOFEMPrismInterface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<CutMeshInterface>(IDD_MOFEMCutMesh); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = regSubInterface<SeriesRecorder>(IDD_MOFEMSeriesRecorder); CHKERRABORT(PETSC_COMM_SELF,ierr);
  #ifdef WITH_TETGEN
  ierr = regSubInterface<TetGenInterface>(IDD_MOFEMTetGegInterface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  #endif
  #ifdef WITH_MED
  ierr = regSubInterface<MedInterface>(IDD_MOFEMMedInterface); CHKERRABORT(PETSC_COMM_SELF,ierr);
  #endif

  // Register MOFEM events in PETSc
  PetscLogEventRegister("FE_preProcess",0,&MOFEM_EVENT_preProcess);
  PetscLogEventRegister("FE_operator",0,&MOFEM_EVENT_operator);
  PetscLogEventRegister("FE_postProcess",0,&MOFEM_EVENT_postProcess);
  PetscLogEventRegister("MoFEMCreateMat",0,&MOFEM_EVENT_createMat);

  // Duplicate PETSc communicator
  ierr = PetscCommDuplicate(comm,&cOmm,NULL); CHKERRABORT(comm,ierr);
  MPI_Comm_size(cOmm,&sIze);
  MPI_Comm_rank(cOmm,&rAnk);
  // CHeck if moab has set communicator if not set communicator interbally
  ParallelComm* pComm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pComm == NULL) {
    pComm =  new ParallelComm(&moab,comm);
  }

  // Initialize database
  ierr = getTags(); CHKERRABORT(cOmm,ierr);
  ierr = clearMap(); CHKERRABORT(cOmm,ierr);
  basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  ierr = initialiseDatabaseFromMesh(verbose); CHKERRABORT(cOmm,ierr);

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
  // Pop error handler
  if(isGloballyInitialised) {
    if(!flg) {
      PetscPopErrorHandler();
    }
    isGloballyInitialised = false;
  }
  // Destroy communicator
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
  MoFEMFunctionBeginHot;
  // Cleaning databases in interfaces
  ierr = getInterface<SeriesRecorder>()->clearMap(); CHKERRQ(ierr);
  ierr = getInterface<MeshsetsManager>()->clearMap(); CHKERRQ(ierr);
  ierr = getInterface<CoordSystemsManager>()->clearMap(); CHKERRQ(ierr);
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
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::addPrismToDatabase(const EntityHandle prism,int verb) {
  MoFEMFunctionBeginHot;
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
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::getTags(int verb) {
  MoFEMFunctionBeginHot;

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
    //Safety nets
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
    BitRefEdges def_bit_edge = 0;
    rval = moab.tag_get_handle(
      "_RefBitEdge",
      sizeof(BitRefEdges),
      MB_TYPE_OPAQUE,
      th_RefBitEdge,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,
      &def_bit_edge
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
  ierr = getInterface(meshsets_manager_ptr); CHKERRQ(ierr);
  ierr = meshsets_manager_ptr->getTags(verb); CHKERRQ(ierr);

  // Series recorder
  SeriesRecorder *series_recorder_ptr;
  ierr = getInterface(series_recorder_ptr); CHKERRQ(ierr);
  ierr = series_recorder_ptr->getTags(verb); CHKERRQ(ierr);

  //Coordinate systems
  CoordSystemsManager *cs_manger_ptr;
  ierr = getInterface(cs_manger_ptr); CHKERRQ(ierr);
  ierr = cs_manger_ptr->getTags(verb); CHKERRQ(ierr);


  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::clear_database(int verb) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  ierr = clearMap(); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::rebuild_database(int verb) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;
  ierr = clearMap(); CHKERRQ(ierr);
  ierr = initialiseDatabaseFromMesh(verb); CHKERRQ(ierr);
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::initialiseDatabaseFromMesh(int verb) {
  MoFEMFunctionBeginHot;
  if(verb==-1) verb = verbose;

  CoordSystemsManager *cs_manger_ptr;
  ierr = getInterface(cs_manger_ptr); CHKERRQ(ierr);

  // Initialize coordinate systems
  ierr = cs_manger_ptr->initialiseDatabaseFromMesh(verb); CHKERRQ(ierr);

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
        p.first->get()->getEntSequenceContainer()->push_back(ents_array);
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
            //PetscPrintf(cOmm,"Warring: this entity should be already in refined finite elements database");
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
  ierr = getInterface(m_manger_ptr); CHKERRQ(ierr);
  ierr = m_manger_ptr->initialiseDatabaseFromMesh(verb); CHKERRQ(ierr);
  SeriesRecorder *series_recorder_ptr;
  ierr = getInterface(series_recorder_ptr); CHKERRQ(ierr);
  ierr = series_recorder_ptr->initialiseDatabaseFromMesh(verb); CHKERRQ(ierr);

  MoFEMFunctionReturnHot(0);
}

// cubit meshsets

PetscErrorCode Core::get_fields(const Field_multiIndex **fields_ptr) const {
  MoFEMFunctionBeginHot;
  *fields_ptr = &fIelds;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const {
  MoFEMFunctionBeginHot;
  *refined_entities_ptr = &refinedEntities;
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::get_ref_finite_elements(const RefElement_multiIndex **refined_finite_elements_ptr) const {
  MoFEMFunctionBeginHot;
  *refined_finite_elements_ptr = &refinedFiniteElements;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_problem(const std::string &problem_name,const Problem **problem_ptr) const {
  MoFEMFunctionBeginHot;
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
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_problems(const Problem_multiIndex **problems_ptr) const {
  MoFEMFunctionBeginHot;
  *problems_ptr = &pRoblems;
  MoFEMFunctionReturnHot(0);
}

PetscErrorCode Core::get_field_ents(const FieldEntity_multiIndex **field_ents) const {
  MoFEMFunctionBeginHot;
  *field_ents = &entsFields;
  MoFEMFunctionReturnHot(0);
}
PetscErrorCode Core::get_dofs(const DofEntity_multiIndex **dofs_ptr) const {
  MoFEMFunctionBeginHot;
  *dofs_ptr = &dofsField;
  MoFEMFunctionReturnHot(0);
}

MeshsetsManager* Core::get_meshsets_manager_ptr() {
  MeshsetsManager* meshsets_manager_ptr;
  getInterface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

const MeshsetsManager* Core::get_meshsets_manager_ptr() const {
  MeshsetsManager* meshsets_manager_ptr;
  getInterface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

}
