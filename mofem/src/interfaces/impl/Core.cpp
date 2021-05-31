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

#include <MoFEM.hpp>

#include "impl/ErrorHandler.cpp"

extern "C" {
void macro_is_deprecated_using_deprecated_function() {}
}

namespace MoFEM {

Core::WrapMPIComm::WrapMPIComm(MPI_Comm &comm, MPI_Comm &duplicated_comm,
                               bool petsc)
    : comm(comm), duplicatedComm(duplicated_comm), isPetscComm(petsc) {
  if (isPetscComm) {
    ierr = PetscCommDuplicate(comm, &duplicated_comm, NULL);
    CHKERRABORT(comm, ierr);
  } else {
    int ierr = MPI_Comm_dup(comm, &duplicated_comm);
    if(ierr) {
      THROW_MESSAGE("MPI_Comm_dup not working");
    }
  }
}
Core::WrapMPIComm::~WrapMPIComm() {
  if (isPetscComm) {
    ierr = PetscCommDestroy(&duplicatedComm);
    CHKERRABORT(comm, ierr);
  } else {
    int ierr = MPI_Comm_free(&duplicatedComm);
    if (ierr) {
      THROW_MESSAGE("MPI_Comm_free not working");
    }
  }
}

constexpr const int CoreTmp<0>::value;
constexpr const int CoreTmp<-1>::value;

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
int Core::mpiInitialised;
PetscBool Core::isInitialized;

MoFEMErrorCode Core::Initialize(int *argc, char ***args, const char file[],
                                const char help[]) {

  MPI_Initialized(&mpiInitialised);
  if (!mpiInitialised)
    MPI_Init(argc, args);

  PetscInitialized(&isInitialized);
  if (isInitialized == PETSC_FALSE) {
    PetscInitialize(argc, args, file, help);
    PetscPushErrorHandler(mofem_error_handler, PETSC_NULL);
  }

  LogManager::createDefaultSinks(MPI_COMM_WORLD);
  PetscVFPrintf = LogManager::logPetscFPrintf;
  CHKERR LogManager::getOptions();
  isGloballyInitialised = true;

  MOFEM_LOG_CHANNEL("WORLD");
  char petsc_version[255];
  CHKERR PetscGetVersion(petsc_version, 255);
  MOFEM_LOG_C("WORLD", Sev::inform, "MoFEM version %d.%d.%d (%s %s)",
              MoFEM_VERSION_MAJOR, MoFEM_VERSION_MINOR, MoFEM_VERSION_BUILD,
              MOAB_VERSION_STRING, petsc_version);
  MOFEM_LOG_C("WORLD", Sev::inform, "git commit id %s", GIT_SHA1_NAME);

  auto log_time = [&](const auto perefix, auto time) {
    MOFEM_LOG("WORLD", Sev::inform)
        << perefix << time.date().year() << "-" << time.date().month() << "-"
        << time.date().day() << " " << time.time_of_day().hours() << ":"
        << time.time_of_day().minutes() << ":" << time.time_of_day().seconds();
  };

  // Get current system time
  log_time("Local time: ", boost::posix_time::second_clock::local_time());
  log_time("UTC time: ", boost::posix_time::second_clock::universal_time());

  return MOFEM_SUCCESS;
}

MoFEMErrorCode Core::Finalize() {
  if (isGloballyInitialised) {
    PetscPopErrorHandler();
    isGloballyInitialised = false;

    if (isInitialized == PETSC_FALSE) {
      PetscBool is_finalized;
      PetscFinalized(&is_finalized);
      if (!is_finalized)
        PetscFinalize();
    }

    if (!mpiInitialised) {
      int mpi_finalized;
      MPI_Finalized(&mpi_finalized);
      if (!mpi_finalized)
        MPI_Finalize();
    }
  }

  return 0;
}

// Use SFINAE to decide which template should be run,
// if exist getSubInterfaceOptions run this one.
template <class T>
static auto get_sub_iface_options_imp(T *const ptr, int)
    -> decltype(ptr->getSubInterfaceOptions()) {
  return ptr->getSubInterfaceOptions();
};

// Use SFINAE to decide which template should be run,
// if getSubInterfaceOptions not exist run this one.
template <class T>
static auto get_sub_iface_options_imp(T *const ptr, long) -> MoFEMErrorCode {
  return 0;
};

template <class IFACE>
MoFEMErrorCode Core::regSubInterface(const MOFEMuuid &uid) {
  MoFEMFunctionBegin;
  CHKERR registerInterface<IFACE>(uid, false);
  unsigned long int id = uid.uUId.to_ulong();
  IFACE *ptr = new IFACE(*this);

  // If sub interface has function getSubInterfaceOptions run
  // it after construction. getSubInterfaceOptions is used to
  // get parameters from command line.
  // See SFINAE:
  // https://stackoverflow.com/questions/257288/is-it-possible-to-write-a-template-to-check-for-a-functions-existence
  // https://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error
  auto get_sub_iface_options = [](auto *const ptr) {
    return get_sub_iface_options_imp(ptr, 0);
  };
  CHKERR get_sub_iface_options(ptr);

  iFaces.insert(id, ptr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::coreGenericConstructor(moab::Interface &moab,
                                            MPI_Comm comm, const int verbose) {
  MoFEMFunctionBegin;

  // This is deprecated ONE should use MoFEM::Core::Initialize
  if (!isGloballyInitialised)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "MoFEM globally is not initialised, call MoFEM::Core::Initialize");

  // Create duplicate communicator
  wrapMPIMOABComm = boost::make_shared<WrapMPIComm>(comm, moabComm, false);
  MPI_Comm_size(mofemComm, &sIze);
  MPI_Comm_rank(mofemComm, &rAnk);

  // CHeck if moab has set communicator if not set communicator interbally
  ParallelComm *pComm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
  if (pComm == NULL)
    pComm = new ParallelComm(&moab, moabComm);

  // Register interfaces for this implementation
  CHKERR registerInterface<UnknownInterface>(IDD_MOFEMUnknown);
  CHKERR registerInterface<CoreInterface>(IDD_MOFEMCoreInterface);
  CHKERR registerInterface<DeprecatedCoreInterface>(
      IDD_MOFEMDeprecatedCoreInterface);

  // Register MOFEM events in PETSc
  PetscLogEventRegister("FE_preProcess", 0, &MOFEM_EVENT_preProcess);
  PetscLogEventRegister("FE_operator", 0, &MOFEM_EVENT_operator);
  PetscLogEventRegister("FE_postProcess", 0, &MOFEM_EVENT_postProcess);
  PetscLogEventRegister("MoFEMCreateMat", 0, &MOFEM_EVENT_createMat);

  MoFEMFunctionReturn(0);
}

Core::CoreTmp(moab::Interface &moab, ///< MoAB interface
              MPI_Comm comm,         ///< MPI communicator
              const int verbose      ///< Verbosity level

              )
    : CoreTmp(moab, comm, verbose, CoreValue<0>()) {

  // Register sub-interfaces
  ierr = this->registerSubInterfaces();
  CHKERRABORT(comm, ierr);
  ierr = this->clearMap();
  CHKERRABORT(comm, ierr);
  ierr = this->getTags();
  CHKERRABORT(comm, ierr);
  ierr = this->getOptions(verbose);
  CHKERRABORT(comm, ierr);

  this->basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  setRefEntBasicDataPtr(*this, this->basicEntityDataPtr);

  ierr = this->initialiseDatabaseFromMesh(verbose);
  CHKERRABORT(comm, ierr);
}

CoreTmp<-1>::CoreTmp(
    moab::Interface &moab,      ///< MoAB interface
    MPI_Comm comm,              ///< MPI communicator
    const int verbose          ///< Verbosity level

    )
    : CoreTmp<0>(moab, comm, verbose, CoreValue<-1>()) {

  // Register sub-interfaces
  ierr = this->registerSubInterfaces();
  CHKERRABORT(comm, ierr);
  ierr = this->clearMap();
  CHKERRABORT(comm, ierr);
  ierr = this->getTags();
  CHKERRABORT(comm, ierr);
  ierr = this->getOptions(verbose);
  CHKERRABORT(comm, ierr);

  this->basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  setRefEntBasicDataPtr(*this, this->basicEntityDataPtr);

  ierr = this->initialiseDatabaseFromMesh(verbose);
  CHKERRABORT(comm, ierr);
}

Core::~CoreTmp() {
  PetscBool is_finalized;
  PetscFinalized(&is_finalized);
  // Destroy interfaces
  iFaces.clear();
  // This is deprecated ONE should use MoFEM::Core::Initialize
  if (isGloballyInitialised && is_finalized) {
    isGloballyInitialised = false;
  }
}

MoFEMErrorCode Core::initialiseDatabaseFromMesh(int verb) {
  MOFEM_LOG_CHANNEL("WORLD");
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  Range ref_elems_to_add;

  auto m_moab = &get_moab();

  // Initialize database
  Range meshsets;
  CHKERR get_moab().get_entities_by_type(0, MBENTITYSET, meshsets, false);
  Range special_meshsets;
  for (auto mit : meshsets) {
    BitFieldId field_id;
    // Get bit id form field tag
    CHKERR get_moab().tag_get_data(th_FieldId, &mit, 1, &field_id);
    // Check if meshset if field meshset
    if (field_id != 0) {

      auto p = fIelds.insert(boost::make_shared<Field>(moab, mit));

      if (verb > QUIET)
        MOFEM_LOG("WORLD", Sev::verbose) << "Read field " << **p.first;

      if (!p.second) {
        // Field meshset exists, remove duplicate meshsets from other
        // processors.
        Range ents;
        CHKERR get_moab().get_entities_by_handle(mit, ents, true);
        CHKERR get_moab().add_entities((*p.first)->getMeshset(), ents);
        CHKERR get_moab().delete_entities(&mit, 1);
      } else {
        special_meshsets.insert(mit);
      }
    }
    // Check for finite elements
    BitFieldId fe_id;
    // Get bit id from fe tag
    CHKERR get_moab().tag_get_data(th_FEId, &mit, 1, &fe_id);
    // check if meshset is finite element meshset
    if (fe_id != 0) {
      std::pair<FiniteElement_multiIndex::iterator, bool> p =
          finiteElements.insert(
              boost::shared_ptr<FiniteElement>(new FiniteElement(moab, mit)));
      if (verb > QUIET)
        MOFEM_LOG("WORLD", Sev::verbose) << "Read finite element " << **p.first;

      Range ents;
      CHKERR get_moab().get_entities_by_type(mit, MBENTITYSET, ents, false);
      CHKERR get_moab().get_entities_by_handle(mit, ents, true);
      ref_elems_to_add.merge(ents);
      if (!p.second) {
        // Finite element mesh set exist, could be created on other processor.
        // Remove duplicate.
        CHKERR get_moab().add_entities((*p.first)->getMeshset(), ents);
        CHKERR get_moab().delete_entities(&mit, 1);
      } else {
        special_meshsets.insert(mit);
      }
    }
    BitProblemId problem_id;
    // get bit id form problem tag
    CHKERR get_moab().tag_get_data(th_ProblemId, &mit, 1, &problem_id);
    // check if meshset if problem meshset
    if (problem_id != 0) {
      std::pair<Problem_multiIndex::iterator, bool> p =
          pRoblems.insert(Problem(moab, mit));
      if (verb > QUIET) {
        MOFEM_LOG("WORLD", Sev::verbose) << "Read problem " << *p.first;
        MOFEM_LOG("WORLD", Sev::noisy)
            << "\tBitRef " << p.first->getBitRefLevel() << " BitMask "
            << p.first->getMaskBitRefLevel();
      }

      if (!p.second) {
        // Problem meshset exists, could be created on other processor.
        // Remove duplicate.
        Range ents;
        CHKERR get_moab().get_entities_by_handle(mit, ents, true);
        CHKERR get_moab().get_entities_by_type(mit, MBENTITYSET, ents, true);
        CHKERR get_moab().add_entities(p.first->meshset, ents);
        CHKERR get_moab().delete_entities(&mit, 1);
      } else {
        special_meshsets.insert(mit);
      }
    }
  }

  // Add entities to database
  Range bit_ref_ents;
  CHKERR get_moab().get_entities_by_handle(0, bit_ref_ents, false);
  bit_ref_ents = subtract(bit_ref_ents, special_meshsets);
  CHKERR getInterface<BitRefManager>()->filterEntitiesByRefLevel(
      BitRefLevel().set(), BitRefLevel().set(), bit_ref_ents);
  CHKERR getInterface<BitRefManager>()->setEntitiesBitRefLevel(bit_ref_ents);
  CHKERR getInterface<BitRefManager>()->setElementsBitRefLevel(
      ref_elems_to_add);

  // Build field entities
  for (auto field : fIelds) {
    if (field->getSpace() != NOSPACE) {
      Range ents_of_id_meshset;
      CHKERR get_moab().get_entities_by_handle(field->getMeshset(),
                                               ents_of_id_meshset, false);
      CHKERR set_field_order(ents_of_id_meshset, field->getId(), -1, verb);
    }
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

MoFEMErrorCode Core::setMoabInterface(moab::Interface &new_moab, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  // clear moab database
  CHKERR clearMap();

  // set new reference
  moab = std::ref(new_moab);

  // check if moab has set communicator if not set communicator internally
  ParallelComm *pComm = ParallelComm::get_pcomm(&new_moab, MYPCOMM_INDEX);
  if (pComm == NULL) {
    pComm = new ParallelComm(&new_moab, moabComm);
  }

  // create MoFEM tags
  CHKERR getTags();

  // Create basic entity data struture
  basicEntityDataPtr = boost::make_shared<BasicEntityData>(moab);
  setRefEntBasicDataPtr(*this, this->basicEntityDataPtr);

  // Initalise database
  CHKERR this->initialiseDatabaseFromMesh(verb);

  MoFEMFunctionReturn(0);
};

MoFEMErrorCode Core::registerSubInterfaces() {
  MoFEMFunctionBegin;

  iFaces.clear();

  // Register sub interfaces
  CHKERR regSubInterface<LogManager>(IDD_MOFEMLogManager);
  CHKERR regSubInterface<Simple>(IDD_MOFEMSimple);
  CHKERR regSubInterface<PipelineManager>(IDD_MOFEMBasic);
  CHKERR regSubInterface<ProblemsManager>(IDD_MOFEMProblemsManager);
  CHKERR regSubInterface<MatrixManager>(IDD_MOFEMMatrixManager);
  CHKERR regSubInterface<ISManager>(IDD_MOFEMISManager);
  CHKERR regSubInterface<VecManager>(IDD_MOFEMVEC);
  CHKERR regSubInterface<FieldBlas>(IDD_MOFEMFieldBlas);
  CHKERR regSubInterface<BitRefManager>(IDD_MOFEMBitRefManager);
  CHKERR regSubInterface<Tools>(IDD_MOFEMTools);
  CHKERR regSubInterface<CommInterface>(IDD_MOFEMComm);
  CHKERR regSubInterface<MeshsetsManager>(IDD_MOFEMMeshsetsManager);
  CHKERR regSubInterface<NodeMergerInterface>(IDD_MOFEMNodeMerger);
  CHKERR regSubInterface<BitLevelCoupler>(IDD_MOFEMBitLevelCoupler);
  CHKERR regSubInterface<PrismsFromSurfaceInterface>(
      IDD_MOFEMPrismsFromSurface);
  CHKERR regSubInterface<MeshRefinement>(IDD_MOFEMMeshRefine);
  CHKERR regSubInterface<PrismInterface>(IDD_MOFEMPrismInterface);
  CHKERR regSubInterface<CutMeshInterface>(IDD_MOFEMCutMesh);
  CHKERR regSubInterface<SeriesRecorder>(IDD_MOFEMSeriesRecorder);
#ifdef WITH_TETGEN
  CHKERR regSubInterface<TetGenInterface>(IDD_MOFEMTetGegInterface);
#endif
#ifdef WITH_MED
  CHKERR regSubInterface<MedInterface>(IDD_MOFEMMedInterface);
#endif
  CHKERR regSubInterface<FieldEvaluatorInterface>(IDD_MOFEMFieldEvaluator);
  CHKERR regSubInterface<BcManager>(IDD_MOFEMBcManager);

  MoFEMFunctionReturn(0);
};

BitFieldId Core::getFieldShift() {
  if (*fShift >= BITFIELDID_SIZE)
    THROW_MESSAGE("Number of field elements exceeded");
  return BitFieldId().set(((*fShift)++) - 1);
}
BitFEId Core::getFEShift() {
  if (*feShift >= BitFEId().set().to_ulong())
    THROW_MESSAGE("Number of finite elements exceeded");
  return BitFEId(1 << (((*feShift)++) - 1));
}

BitProblemId Core::getProblemShift() {
  if (*pShift >= BitProblemId().set().to_ulong())
    THROW_MESSAGE("Number of problems exceeded");
  return BitProblemId(1 << (((*pShift)++) - 1));
}

MoFEMErrorCode Core::clearMap() {
  MoFEMFunctionBegin;
  // Cleaning databases in interfaces
  CHKERR getInterface<SeriesRecorder>()->clearMap();
  CHKERR getInterface<MeshsetsManager>()->clearMap();
  CHKERR getInterface<CutMeshInterface>()->clearMap();
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
  MoFEMFunctionReturn(0);
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
    CHKERR get_moab().get_connectivity(prism, conn, num_nodes, true);
    Range face_side3, face_side4;
    CHKERR get_moab().get_adjacencies(conn, 3, 2, false, face_side3);
    CHKERR get_moab().get_adjacencies(&conn[3], 3, 2, false, face_side4);
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
  MoFEMFunctionBegin;

  const EntityHandle root_meshset = get_moab().get_root_set();
  if (root_meshset) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Root meshset should be 0");
  }

  // Set version
  {
    Version version;
    CHKERR getFileVersion(moab, version);
  }

  // Global Variables
  {

    auto check_tag_allocated = [](auto &rval) {
      MoFEMFunctionBeginHot;
      if (rval == MB_ALREADY_ALLOCATED)
        rval = MB_SUCCESS;
      else
        CHKERRG(rval);
      MoFEMFunctionReturnHot(0);
    };

    // Fields
    int def_shift = 1;
    rval = get_moab().tag_get_handle("_FieldShift", 1, MB_TYPE_INTEGER,
                                     th_FieldShift, MB_TAG_CREAT | MB_TAG_MESH,
                                     &def_shift);
    CHKERR check_tag_allocated(rval);

    const void *tag_data[1];
    CHKERR get_moab().tag_get_by_ptr(th_FieldShift, &root_meshset, 1, tag_data);
    fShift = (int *)tag_data[0];
    // FE
    rval = get_moab().tag_get_handle("_FEShift", 1, MB_TYPE_INTEGER, th_FEShift,
                                     MB_TAG_CREAT | MB_TAG_MESH, &def_shift);
    CHKERR check_tag_allocated(rval);

    CHKERR get_moab().tag_get_by_ptr(th_FEShift, &root_meshset, 1, tag_data);
    feShift = (int *)tag_data[0];
    // Problem
    rval = get_moab().tag_get_handle("_ProblemShift", 1, MB_TYPE_INTEGER,
                                     th_ProblemShift,
                                     MB_TAG_CREAT | MB_TAG_MESH, &def_shift);
    CHKERR check_tag_allocated(rval);

    CHKERR get_moab().tag_get_by_ptr(th_ProblemShift, &root_meshset, 1,
                                     tag_data);
    pShift = (int *)tag_data[0];
    // Safety nets
    int def_bool = 0;
    rval = get_moab().tag_get_handle("_MoFEMBuild", 1, MB_TYPE_INTEGER,
                                     th_MoFEMBuild, MB_TAG_CREAT | MB_TAG_MESH,
                                     &def_bool);
    CHKERR check_tag_allocated(rval);

    CHKERR get_moab().tag_get_by_ptr(th_MoFEMBuild, &root_meshset, 1,
                                     (const void **)&buildMoFEM);
  }

  // Tags saved in vtk-files
  {
    const int def_part = -1;
    CHKERR get_moab().tag_get_handle("PARTITION", 1, MB_TYPE_INTEGER, th_Part,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, &def_part);
  }

  // Tags Ref
  {
    const int def_part = -1;
    CHKERR get_moab().tag_get_handle("_MeshsetPartition", 1, MB_TYPE_INTEGER,
                                     th_Part, MB_TAG_CREAT | MB_TAG_SPARSE,
                                     &def_part);
    EntityHandle def_handle = 0;
    CHKERR get_moab().tag_get_handle("_RefParentHandle", 1, MB_TYPE_HANDLE,
                                     th_RefParentHandle,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, &def_handle);
    BitRefLevel def_bit_level = 0;
    CHKERR get_moab().tag_get_handle(
        "_RefBitLevel", sizeof(BitRefLevel), MB_TYPE_OPAQUE, th_RefBitLevel,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_bit_level);
    BitRefLevel def_bit_level_mask = BitRefLevel().set();
    CHKERR get_moab().tag_get_handle(
        "_RefBitLevelMask", sizeof(BitRefLevel), MB_TYPE_OPAQUE,
        th_RefBitLevel_Mask, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_bit_level_mask);
    BitRefEdges def_bit_edge = 0;
    CHKERR get_moab().tag_get_handle(
        "_RefBitEdge", sizeof(BitRefEdges), MB_TYPE_OPAQUE, th_RefBitEdge,
        MB_TAG_CREAT | MB_TAG_SPARSE | MB_TAG_BYTES, &def_bit_edge);
    const int def_type[] = {0, 0};
    CHKERR get_moab().tag_get_handle("_RefType", 2, MB_TYPE_INTEGER, th_RefType,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, def_type);
  }

  // Tags Field
  {
    const unsigned long int def_id = 0;
    CHKERR get_moab().tag_get_handle(
        "_FieldId", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FieldId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    FieldSpace def_space = LASTSPACE;
    CHKERR get_moab().tag_get_handle(
        "_FieldSpace", sizeof(FieldSpace), MB_TYPE_OPAQUE, th_FieldSpace,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_space);
    FieldApproximationBase def_base = LASTBASE;
    CHKERR get_moab().tag_get_handle(
        "_FieldBase", sizeof(FieldApproximationBase), MB_TYPE_OPAQUE,
        th_FieldBase, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_base);
    const int def_val_len = 0;
    CHKERR get_moab().tag_get_handle(
        "_FieldName", def_val_len, MB_TYPE_OPAQUE, th_FieldName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERR get_moab().tag_get_handle(
        "_FieldName_DataNamePrefix", def_val_len, MB_TYPE_OPAQUE,
        th_FieldName_DataNamePrefix,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
  }

  // Tags FE
  {
    const unsigned long int def_id = 0;
    const int def_val_len = 0;
    CHKERR get_moab().tag_get_handle(
        "_FEId", sizeof(BitFEId), MB_TYPE_OPAQUE, th_FEId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERR get_moab().tag_get_handle(
        "_FEName", def_val_len, MB_TYPE_OPAQUE, th_FEName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    CHKERR get_moab().tag_get_handle(
        "_FEIdCol", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdCol,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERR get_moab().tag_get_handle(
        "_FEIdRow", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdRow,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERR get_moab().tag_get_handle(
        "_FEIdData", sizeof(BitFieldId), MB_TYPE_OPAQUE, th_FEIdData,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
  }

  // Tags Problem
  {
    const unsigned long int def_id = 0;
    const int def_val_len = 0;
    CHKERR get_moab().tag_get_handle(
        "_ProblemId", sizeof(BitProblemId), MB_TYPE_OPAQUE, th_ProblemId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERR get_moab().tag_get_handle(
        "_ProblemFEId", sizeof(BitFEId), MB_TYPE_OPAQUE, th_ProblemFEId,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE, &def_id);
    CHKERR get_moab().tag_get_handle(
        "_ProblemName", def_val_len, MB_TYPE_OPAQUE, th_ProblemName,
        MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_VARLEN | MB_TAG_SPARSE, NULL);
    DofIdx def_nbdofs = 0;
    CHKERR get_moab().tag_get_handle(
        "_ProblemNbDofsRow", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemNbDofsRow, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
    CHKERR get_moab().tag_get_handle(
        "_ProblemNbDofsCol", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemNbDofsCol, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
    CHKERR get_moab().tag_get_handle(
        "_ProblemLocalNbDofsRow", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemLocalNbDofRow, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
    CHKERR get_moab().tag_get_handle(
        "_ProblemGhostNbDofsRow", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemGhostNbDofRow, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
    CHKERR get_moab().tag_get_handle(
        "_ProblemLocalNbDofsCol", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemLocalNbDofCol, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
    CHKERR get_moab().tag_get_handle(
        "_ProblemGhostNbDofsCol", sizeof(DofIdx), MB_TYPE_OPAQUE,
        th_ProblemGhostNbDofCol, MB_TAG_CREAT | MB_TAG_BYTES | MB_TAG_SPARSE,
        &def_nbdofs);
  }

  // Meshsets with boundary conditions and material sets
  MeshsetsManager *meshsets_manager_ptr;
  CHKERR getInterface(meshsets_manager_ptr);
  CHKERR meshsets_manager_ptr->getTags(verb);

  // Series recorder
  SeriesRecorder *series_recorder_ptr;
  CHKERR getInterface(series_recorder_ptr);
  CHKERR series_recorder_ptr->getTags(verb);

  MoFEMFunctionReturn(0);
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
  CHKERR this->clearMap();
  CHKERR this->getTags(verb);
  CHKERR this->initialiseDatabaseFromMesh(verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::set_moab_interface(moab::Interface &new_moab, int verb) {
  return this->setMoabInterface(new_moab, verb);
};

MoFEMErrorCode CoreTmp<-1>::set_moab_interface(moab::Interface &new_moab,
                                               int verb) {
  return this->setMoabInterface(new_moab, verb);
};

MoFEMErrorCode Core::getOptions(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  CHKERR PetscOptionsBegin(mofemComm, optionsPrefix.c_str(), "Mesh cut options",
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

MoFEMErrorCode
Core::get_finite_elements(const FiniteElement_multiIndex **fe_ptr) const {
  MoFEMFunctionBeginHot;
  *fe_ptr = &finiteElements;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::get_ents_finite_elements(
    const EntFiniteElement_multiIndex **fe_ent_ptr) const {
  MoFEMFunctionBeginHot;
  *fe_ent_ptr = &entsFiniteElements;
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

MoFEMErrorCode Core::get_ents_elements_adjacency(
    const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
        *dofs_elements_adjacency) const {
  MoFEMFunctionBeginHot;
  *dofs_elements_adjacency = &entFEAdjacencies;
  MoFEMFunctionReturnHot(0);
}

const FieldEntityEntFiniteElementAdjacencyMap_multiIndex *
Core::get_ents_elements_adjacency() const {
  return &entFEAdjacencies;
}

const Field_multiIndex *Core::get_fields() const { return &fIelds; }
const RefEntity_multiIndex *Core::get_ref_ents() const {
  return &refinedEntities;
}
const RefElement_multiIndex *Core::get_ref_finite_elements() const {
  return &refinedFiniteElements;
}
const FiniteElement_multiIndex *Core::get_finite_elements() const {
  return &finiteElements;
}
const EntFiniteElement_multiIndex *Core::get_ents_finite_elements() const {
  return &entsFiniteElements;
}
const FieldEntity_multiIndex *Core::get_field_ents() const {
  return &entsFields;
}
const DofEntity_multiIndex *Core::get_dofs() const { return &dofsField; }
const Problem *Core::get_problem(const std::string problem_name) const {
  const Problem *prb;
  CHKERR get_problem(problem_name, &prb);
  return prb;
}
const Problem_multiIndex *Core::get_problems() const { return &pRoblems; }

template <int V, typename std::enable_if<(V >= 0), int>::type * = nullptr>
void set_ref_ent_basic_data_ptr_impl(boost::shared_ptr<BasicEntityData> &ptr) {
  RefEntityTmp<V>::basicDataPtr = ptr;
};

template <int V, typename std::enable_if<(V < 0), int>::type * = nullptr>
void set_ref_ent_basic_data_ptr_impl(boost::shared_ptr<BasicEntityData> &ptr) {
  return;
};

void Core::setRefEntBasicDataPtr(MoFEM::Interface &m_field,
                                      boost::shared_ptr<BasicEntityData> &ptr) {

  switch (m_field.getValue()) {
  case -1:
    set_ref_ent_basic_data_ptr_impl<-1>(ptr);
    break;
  case 0:
    set_ref_ent_basic_data_ptr_impl<0>(ptr);
    break;
  case 1:
    set_ref_ent_basic_data_ptr_impl<1>(ptr);
    break;
  default:
    THROW_MESSAGE("Core index can vary from -1 to MAX_CORE_TMP");
  }

};

boost::shared_ptr<RefEntityTmp<0>>
Core::makeSharedRefEntity(MoFEM::Interface &m_field, const EntityHandle ent) {

  boost::shared_ptr<RefEntityTmp<0>> ref_ent_ptr;

  switch (m_field.getValue()) {
  case -1:
    ref_ent_ptr = boost::shared_ptr<RefEntityTmp<0>>(

        new RefEntityTmp<-1>(m_field.get_basic_entity_data_ptr(), ent)

    );
    break;
  case 0:
    ref_ent_ptr = boost::shared_ptr<RefEntityTmp<0>>(

        new RefEntityTmp<0>(m_field.get_basic_entity_data_ptr(), ent)

    );
    break;
  case 1:
    ref_ent_ptr = boost::shared_ptr<RefEntityTmp<0>>(

        new RefEntityTmp<1>(m_field.get_basic_entity_data_ptr(), ent)

    );
    break;
  default:
    THROW_MESSAGE("Core index can vary from -1 to MAX_CORE_TMP");
  }

  return ref_ent_ptr;
}

boost::shared_ptr<RefEntityTmp<0>>
Core::make_shared_ref_entity(const EntityHandle ent) {
  return this->makeSharedRefEntity(*this, ent);
}

boost::shared_ptr<RefEntityTmp<0>>
CoreTmp<-1>::make_shared_ref_entity(const EntityHandle ent) {
  return this->makeSharedRefEntity(*this, ent);
}

} // namespace MoFEM
