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

MoFEMErrorCode Core::query_interface(const MOFEMuuid &uuid,
                                     UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMCoreInterface) {
    *iface = static_cast<CoreInterface *>(const_cast<Core *>(this));
    MoFEMFunctionReturnHot(0);
  } else if (uuid == IDD_MOFEMDeprecatedCoreInterface) {
    *iface = static_cast<DeprecatedCoreInterface *>(const_cast<Core *>(this));
    MoFEMFunctionReturnHot(0);
  }
  // Get sub-interface
  unsigned long int id = uuid.uUId.to_ulong();
  boost::ptr_map<unsigned long, UnknownInterface>::iterator it;
  it = iFaces.find(id);
  if (it != iFaces.end()) {
    *iface = it->second;
    MoFEMFunctionReturnHot(0);
  }
  *iface = NULL;
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

bool Core::isGloballyInitialised = false;

MoFEMErrorCode Core::Initialize(int *argc, char ***args, const char file[],
                                const char help[]) {
  ierr = PetscInitialize(argc, args, file, help);
  CHKERRG(ierr);
  ierr = PetscPushErrorHandler(mofem_error_handler, PETSC_NULL);
  CHKERRG(ierr);
  isGloballyInitialised = true;
  return MOFEM_SUCESS;
}

MoFEMErrorCode Core::Finalize() {
   CHKERRQ(ierr);
  ierr = PetscPopErrorHandler();
  CHKERRG(ierr);
  isGloballyInitialised = false;
  return PetscFinalize();
}

template<class IFACE>
MoFEMErrorCode Core::regSubInterface(const MOFEMuuid& uid) {
  MoFEMFunctionBeginHot;
  ierr = registerInterface<IFACE>(uid); CHKERRG(ierr);
  unsigned long int id = uid.uUId.to_ulong();
  iFaces.insert(id,new IFACE(*this));
  MoFEMFunctionReturnHot(0);
}

Core::Core(moab::Interface &moab, MPI_Comm comm, const int verbose,
           const bool distributed_mesh)
    : moab(moab), cOmm(0), verbose(verbose),
      initaliseAndBuildField(PETSC_FALSE),
      initaliseAndBuildFiniteElements(PETSC_FALSE) {

  // This is deprecated ONE should use MoFEM::Core::Initialize
  if (!isGloballyInitialised) {
    PetscPushErrorHandler(mofem_error_handler, PETSC_NULL);
    isGloballyInitialised = true;
  }

  // Register interfaces for this implementation
  ierr = registerInterface<UnknownInterface>(IDD_MOFEMUnknown);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = registerInterface<CoreInterface>(IDD_MOFEMCoreInterface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = registerInterface<DeprecatedCoreInterface>(
      IDD_MOFEMDeprecatedCoreInterface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // Register sub interfaces
  ierr = regSubInterface<Simple>(IDD_MOFEMSimple);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<ProblemsManager>(IDD_MOFEMProblemsManager);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<ISManager>(IDD_MOFEMISManager);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<VecManager>(IDD_MOFEMVEC);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<FieldBlas>(IDD_MOFEMFieldBlas);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<BitRefManager>(IDD_MOFEMBitRefManager);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<Tools>(IDD_MOFEMTools);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<MeshsetsManager>(IDD_MOFEMMeshsetsManager);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<CoordSystemsManager>(IDD_MOFEMCoordsSystemsManager);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<NodeMergerInterface>(IDD_MOFEMNodeMerger);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<BitLevelCoupler>(IDD_MOFEMBitLevelCoupler);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr =
      regSubInterface<PrismsFromSurfaceInterface>(IDD_MOFEMPrismsFromSurface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<MeshRefinement>(IDD_MOFEMMeshRefine);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<PrismInterface>(IDD_MOFEMPrismInterface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<CutMeshInterface>(IDD_MOFEMCutMesh);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  ierr = regSubInterface<SeriesRecorder>(IDD_MOFEMSeriesRecorder);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
#ifdef WITH_TETGEN
  ierr = regSubInterface<TetGenInterface>(IDD_MOFEMTetGegInterface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
#endif
#ifdef WITH_MED
  ierr = regSubInterface<MedInterface>(IDD_MOFEMMedInterface);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
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
  if (distributed_mesh)
    basicEntityDataPtr->setDistributedMesh();
  else
    basicEntityDataPtr->unSetDistributedMesh();

  ierr = getOptions(verbose); CHKERRABORT(cOmm,ierr);
  ierr = initialiseDatabaseFromMesh(verbose); CHKERRABORT(cOmm,ierr);

  // Print version
  if (verbose > QUIET) {
    char petsc_version[255];
    ierr = PetscGetVersion(petsc_version, 255);
    CHKERRABORT(cOmm, ierr);
    ierr = PetscPrintf(cOmm, "MoFEM version %d.%d.%d (%s %s) \n",
                       MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR,
                       MoFEM_VERSION_BUILD, MOAB_VERSION_STRING, petsc_version);
    CHKERRABORT(cOmm, ierr);
    ierr = PetscPrintf(cOmm, "git commit id %s\n", GIT_SHA1_NAME);
    CHKERRABORT(cOmm, ierr);
  }
}

Core::~Core() {
  PetscBool is_finalized;
  PetscFinalized(&is_finalized);
  // Destroy interfaces
  iFaces.clear();
  // This is deprecated ONE should use MoFEM::Core::Initialize
  if (isGloballyInitialised && is_finalized) {
    isGloballyInitialised = false;
  }
  // Destroy communicator
  if (!is_finalized) {
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

MoFEMErrorCode Core::clearMap() {
  MoFEMFunctionBeginHot;
  // Cleaning databases in interfaces
  ierr = getInterface<SeriesRecorder>()->clearMap(); CHKERRG(ierr);
  ierr = getInterface<MeshsetsManager>()->clearMap(); CHKERRG(ierr);
  ierr = getInterface<CoordSystemsManager>()->clearMap(); CHKERRG(ierr);
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

MoFEMErrorCode Core::addPrismToDatabase(const EntityHandle prism, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  std::pair<RefEntity_multiIndex::iterator, bool> p_ent;
  p_ent = refinedEntities.insert(
      boost::make_shared<RefEntity>(basicEntityDataPtr, prism));
  if (p_ent.second) {
    std::pair<RefElement_multiIndex::iterator, bool> p;
    p = refinedFiniteElements.insert(
        boost::shared_ptr<RefElement>(new RefElement_PRISM(*p_ent.first)));
    int num_nodes;
    const EntityHandle *conn;
    CHKERR moab.get_connectivity(prism, conn, num_nodes, true);
    Range face_side3, face_side4;
    CHKERR moab.get_adjacencies(conn, 3, 2, false, face_side3);
    CHKERR moab.get_adjacencies(&conn[3], 3, 2, false, face_side4);
    if (face_side3.size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "prism don't have side face 3");
    if (face_side4.size() != 1)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "prims don't have side face 4");
    p.first->get()->getSideNumberPtr(*face_side3.begin());
    p.first->get()->getSideNumberPtr(*face_side4.begin());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::getTags(int verb) {
  MoFEMFunctionBeginHot;

  const EntityHandle root_meshset = moab.get_root_set();
  if (root_meshset) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Root meshset should be 0");
  }

  // Set version
  {
    Version version;
    ierr = getFileVersion(moab, version);
    CHKERRG(ierr);
  }

  // Global Variables
  {
    // Fields
    int def_shift = 1;
    rval = moab.tag_get_handle("_FieldShift", 1, MB_TYPE_INTEGER, th_FieldShift,
                               MB_TAG_CREAT | MB_TAG_MESH, &def_shift);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    CHKERRG(rval);
    const void *tag_data[1];
    rval = moab.tag_get_by_ptr(th_FieldShift, &root_meshset, 1, tag_data);
    CHKERRG(rval);
    fShift = (int *)tag_data[0];
    // FE
    rval = moab.tag_get_handle("_FEShift", 1, MB_TYPE_INTEGER, th_FEShift,
                               MB_TAG_CREAT | MB_TAG_MESH, &def_shift);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    CHKERRG(rval);
    rval = moab.tag_get_by_ptr(th_FEShift, &root_meshset, 1, tag_data);
    CHKERRG(rval);
    feShift = (int *)tag_data[0];
    // Problem
    rval = moab.tag_get_handle("_ProblemShift", 1, MB_TYPE_INTEGER,
                               th_ProblemShift, MB_TAG_CREAT | MB_TAG_MESH,
                               &def_shift);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    CHKERRG(rval);
    rval = moab.tag_get_by_ptr(th_ProblemShift, &root_meshset, 1, tag_data);
    CHKERRG(rval);
    pShift = (int *)tag_data[0];
    // Safety nets
    int def_bool = 0;
    rval = moab.tag_get_handle("_MoFEMBuild", 1, MB_TYPE_INTEGER, th_MoFEMBuild,
                               MB_TAG_CREAT | MB_TAG_MESH, &def_bool);
    if (rval == MB_ALREADY_ALLOCATED)
      rval = MB_SUCCESS;
    rval = moab.tag_get_by_ptr(th_MoFEMBuild, &root_meshset, 1,
                               (const void **)&buildMoFEM);
    CHKERRG(rval);
  }

  // Tags saved in vtk-files
  {
    const int def_part = -1;
    rval = moab.tag_get_handle("PARTITION", 1, MB_TYPE_INTEGER, th_Part,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_part);
    CHKERRG(rval);
    int def_elem_type = MBMAXTYPE;
    rval = moab.tag_get_handle("ElemType", 1, MB_TYPE_INTEGER, th_ElemType,
                               MB_TAG_CREAT | MB_TAG_SPARSE, &def_elem_type);
    CHKERRG(rval);
  }

  // Tags Ref
  {
    EntityHandle def_handle = 0;

    rval = moab.tag_get_handle("_RefParentHandle", 1, MB_TYPE_HANDLE,
                               th_RefParentHandle, MB_TAG_CREAT | MB_TAG_SPARSE,
                               &def_handle);
    CHKERRG(rval);
    BitRefLevel def_bit_level = 0;
    rval = moab.tag_get_handle(
        "_RefBitLevel", sizeof(BitRefLevel), MB_TYPE_OPAQUE, th_RefBitLevel,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_bit_level);
    CHKERRG(rval);
    BitRefLevel def_bit_level_mask = BitRefLevel().set();
    rval = moab.tag_get_handle("_RefBitLevelMask", sizeof(BitRefLevel),
                               MB_TYPE_OPAQUE, th_RefBitLevel_Mask,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_bit_level_mask);
    CHKERRG(rval);
    BitRefEdges def_bit_edge = 0;
    rval = moab.tag_get_handle(
        "_RefBitEdge", sizeof(BitRefEdges), MB_TYPE_OPAQUE, th_RefBitEdge,
        MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES, &def_bit_edge);
    CHKERRG(rval);
    const int def_type[] = {0, 0};
    rval = moab.tag_get_handle("_RefType", 2, MB_TYPE_INTEGER, th_RefType,
                               MB_TAG_CREAT | MB_TAG_SPARSE, def_type);
    CHKERRG(rval);
  }

  // Tags Field
  {
    const unsigned long int def_id = 0;

    rval = moab.tag_get_handle(
        "_FieldId", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FieldId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    FieldSpace def_space = LASTSPACE;
    rval                 = moab.tag_get_handle(
        "_FieldSpace", sizeof(FieldSpace), MB_TYPE_OPAQUE, th_FieldSpace,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_space);
    CHKERRG(rval);
    FieldApproximationBase def_base = LASTBASE;
    rval                            = moab.tag_get_handle(
        "_FieldBase", sizeof(FieldApproximationBase), MB_TYPE_OPAQUE,
        th_FieldBase, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_base);
    CHKERRG(rval);
    const int def_val_len = 0;
    rval = moab.tag_get_handle(
        "_FieldName", def_val_len, MB_TYPE_OPAQUE, th_FieldName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_FieldName_DataNamePrefix", def_val_len, MB_TYPE_OPAQUE,
        th_FieldName_DataNamePrefix,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERRG(rval);
  }

  // Tags FE
  {
    const unsigned long int def_id = 0;
    const int def_val_len          = 0;

    rval = moab.tag_get_handle(
        "_FEId", sizeof(BitFEId), MB_TYPE_OPAQUE, th_FEId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_FEName", def_val_len, MB_TYPE_OPAQUE, th_FEName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_FEIdCol", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdCol,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_FEIdRow", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdRow,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_FEIdData", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdData,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
  }

  // Tags Problem
  {
    const unsigned long int def_id = 0;
    const int def_val_len          = 0;

    rval = moab.tag_get_handle(
        "_ProblemId", sizeof(BitProblemId), MB_TYPE_OPAQUE, th_ProblemId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_ProblemFEId", sizeof(BitFEId), MB_TYPE_OPAQUE, th_ProblemFEId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERRG(rval);
    rval = moab.tag_get_handle(
        "_ProblemName", def_val_len, MB_TYPE_OPAQUE, th_ProblemName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERRG(rval);
    DofIdx def_nbdofs = 0;
    rval = moab.tag_get_handle("_ProblemNbDofsRow", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemNbDofsRow,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
    rval = moab.tag_get_handle("_ProblemNbDofsCol", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemNbDofsCol,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
    rval = moab.tag_get_handle("_ProblemLocalNbDofsRow", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemLocalNbDofRow,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
    rval = moab.tag_get_handle("_ProblemGhostNbDofsRow", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemGhostNbDofRow,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
    rval = moab.tag_get_handle("_ProblemLocalNbDofsCol", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemLocalNbDofCol,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
    rval = moab.tag_get_handle("_ProblemGhostNbDofsCol", sizeof(DofIdx),
                               MB_TYPE_OPAQUE, th_ProblemGhostNbDofCol,
                               MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
                               &def_nbdofs);
    CHKERRG(rval);
  }

  // Meshsets with boundary conditions and material sets
  MeshsetsManager *meshsets_manager_ptr;
  ierr = getInterface(meshsets_manager_ptr);
  CHKERRG(ierr);
  ierr = meshsets_manager_ptr->getTags(verb);
  CHKERRG(ierr);

  // Series recorder
  SeriesRecorder *series_recorder_ptr;
  ierr = getInterface(series_recorder_ptr);
  CHKERRG(ierr);
  ierr = series_recorder_ptr->getTags(verb);
  CHKERRG(ierr);

  // Coordinate systems
  CoordSystemsManager *cs_manger_ptr;
  ierr = getInterface(cs_manger_ptr);
  CHKERRG(ierr);
  ierr = cs_manger_ptr->getTags(verb);
  CHKERRG(ierr);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_database(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clearMap();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::rebuild_database(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  CHKERR clearMap();
  CHKERR initialiseDatabaseFromMesh(verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::getOptions(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  CHKERR PetscOptionsBegin(cOmm, optionsPrefix.c_str(), "Mesh cut options",
                           "See MoFEM documentation");

  CHKERR PetscOptionsBool(
      "-mofem_init_fields", "Initialise fields on construction", "",
      initaliseAndBuildField, &initaliseAndBuildField, NULL);

  CHKERR PetscOptionsBool(
      "-mofem_init_fields", "Initialise fields on construction", "",
      initaliseAndBuildFiniteElements, &initaliseAndBuildFiniteElements, NULL);

  // TODO: Add read verbosity level
  // TODO: Add option to initalise problems ??? -  DO WE REALLY NEED THAT

  ierr = PetscOptionsEnd();
  CHKERRG(ierr);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::initialiseDatabaseFromMesh(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  CoordSystemsManager *cs_manger_ptr;
  CHKERR getInterface(cs_manger_ptr);

  // Initialize coordinate systems
  CHKERR cs_manger_ptr->initialiseDatabaseFromMesh(verb);

  Range ref_elems_to_add;

  // Initialize database
  Range meshsets;
  CHKERR moab.get_entities_by_type(0, MBENTITYSET, meshsets, false);
  Range special_meshsets;
  for (Range::iterator mit = meshsets.begin(); mit != meshsets.end(); mit++) {
    BitFieldId field_id;
    // Get bit id form field tag
    CHKERR moab.tag_get_data(th_FieldId, &*mit, 1, &field_id);
    // Check if meshset if field meshset
    if (field_id != 0) {
      std::pair<Field_multiIndex::iterator, bool> p;
      const char *cs_name;
      int cs_name_size;
      boost::shared_ptr<CoordSys> cs_ptr;
      rval = moab.tag_get_by_ptr(cs_manger_ptr->get_th_CoordSysName(), &*mit, 1,
                                 (const void **)&cs_name, &cs_name_size);
      if (rval == MB_SUCCESS && cs_name_size) {
        CHKERR cs_manger_ptr->getCoordSysPtr(std::string(cs_name, cs_name_size),
                                             cs_ptr);
      } else {
        CHKERR cs_manger_ptr->getCoordSysPtr("UNDEFINED", cs_ptr);
      }
      p = fIelds.insert(
          boost::shared_ptr<Field>(new Field(moab, *mit, cs_ptr)));
      if (verb >= VERBOSE) {
        std::ostringstream ss;
        ss << "read field " << **p.first << std::endl;
        PetscPrintf(cOmm, ss.str().c_str());
      }
      special_meshsets.insert(*mit);
    }
    // Check for finite elements
    BitFieldId fe_id;
    // Get bit id from fe tag
    CHKERR moab.tag_get_data(th_FEId, &*mit, 1, &fe_id);
    // check if meshset is finite element meshset
    if (fe_id != 0) {
      std::pair<FiniteElement_multiIndex::iterator, bool> p =
          finiteElements.insert(
              boost::shared_ptr<FiniteElement>(new FiniteElement(moab, *mit)));
      if (verb >= VERBOSE) {
        std::ostringstream ss;
        ss << "read finite element " << **p.first << std::endl;
        PetscPrintf(cOmm, ss.str().c_str());
      }
      NOT_USED(p);
      Range ents;
      CHKERR moab.get_entities_by_type(*mit, MBENTITYSET, ents, false);
      CHKERR moab.get_entities_by_handle(*mit, ents, true);
      ref_elems_to_add.merge(ents);
      special_meshsets.insert(*mit);
    }
    BitProblemId problem_id;
    // get bit id form problem tag
    CHKERR moab.tag_get_data(th_ProblemId, &*mit, 1, &problem_id);
    // check if meshset if problem meshset
    if (problem_id != 0) {
      std::pair<Problem_multiIndex::iterator, bool> p =
          pRoblems.insert(Problem(moab, *mit));
      if (verb >= VERBOSE) {
        std::ostringstream ss;
        ss << "read problem " << *p.first << " bit ref "
           << p.first->getBitRefLevel() << " bit mask "
           << p.first->getMaskBitRefLevel() << std::endl;
        PetscPrintf(cOmm, ss.str().c_str());
      }
      special_meshsets.insert(*mit);
    }
  }

  // Add entities to database
  Range bit_ref_ents;
  CHKERR moab.get_entities_by_handle(0, bit_ref_ents, false);
  bit_ref_ents = subtract(bit_ref_ents,special_meshsets);
  CHKERR getInterface<BitRefManager>()->filterEntitiesByRefLevel(
      BitRefLevel().set(), BitRefLevel().set(), bit_ref_ents);
  CHKERR getInterface<BitRefManager>()->setEntitiesBitRefLevel(bit_ref_ents);
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(
      ref_elems_to_add);

  // Build field entities
  for (Field_multiIndex::iterator fit = fIelds.begin(); fit != fIelds.end();
       ++fit) {
      Range ents_of_id_meshset;
      CHKERR moab.get_entities_by_handle(fit->get()->getMeshset(),
                                         ents_of_id_meshset, false);
      CHKERR set_field_order(ents_of_id_meshset, fit->get()->getId(), -1);
  }

  if (initaliseAndBuildField || initaliseAndBuildFiniteElements) {
    CHKERR build_fields(verb);
    if (initaliseAndBuildFiniteElements) {
      CHKERR build_finite_elements(verb);
    }
  }

  if (verb > VERY_NOISY) {
    list_fields();
    list_finite_elements();
    list_problem();
  }

  // Initialize interfaces
  MeshsetsManager *m_manger_ptr;
  CHKERR getInterface(m_manger_ptr);
  CHKERR m_manger_ptr->initialiseDatabaseFromMesh(verb);
  SeriesRecorder *series_recorder_ptr;
  CHKERR getInterface(series_recorder_ptr);
  CHKERR series_recorder_ptr->initialiseDatabaseFromMesh(verb);

  MoFEMFunctionReturn(0);
}

// cubit meshsets

MoFEMErrorCode Core::get_fields(const Field_multiIndex **fields_ptr) const {
  MoFEMFunctionBeginHot;
  *fields_ptr = &fIelds;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::get_ref_ents(const RefEntity_multiIndex **refined_entities_ptr) const {
  MoFEMFunctionBeginHot;
  *refined_entities_ptr = &refinedEntities;
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode Core::get_ref_finite_elements(
    const RefElement_multiIndex **refined_finite_elements_ptr) const {
  MoFEMFunctionBeginHot;
  *refined_finite_elements_ptr = &refinedFiniteElements;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::get_problem(const std::string &problem_name,
                                 const Problem **problem_ptr) const {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  const ProblemsByName &problems = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = problems.find(problem_name);
  if (p_miit == problems.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "problem < %s > not found, (top tip: check spelling)",
             problem_name.c_str());
  }
  *problem_ptr = &*p_miit;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::get_problems(const Problem_multiIndex **problems_ptr) const {
  MoFEMFunctionBeginHot;
  *problems_ptr = &pRoblems;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::get_field_ents(const FieldEntity_multiIndex **field_ents) const {
  MoFEMFunctionBeginHot;
  *field_ents = &entsFields;
  MoFEMFunctionReturnHot(0);
}
MoFEMErrorCode Core::get_dofs(const DofEntity_multiIndex **dofs_ptr) const {
  MoFEMFunctionBeginHot;
  *dofs_ptr = &dofsField;
  MoFEMFunctionReturnHot(0);
}

MeshsetsManager *Core::get_meshsets_manager_ptr() {
  MeshsetsManager *meshsets_manager_ptr;
  getInterface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

const MeshsetsManager *Core::get_meshsets_manager_ptr() const {
  MeshsetsManager *meshsets_manager_ptr;
  getInterface(meshsets_manager_ptr);
  return meshsets_manager_ptr;
}

} // namespace MoFEM
