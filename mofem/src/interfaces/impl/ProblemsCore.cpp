/** \file ProblemsCore.cpp
 * \brief Managing complexities for problem
 */


#include <MoFEM.hpp>

#define ProblemCoreFunctionBegin                                               \
  MoFEMFunctionBegin;                                                          \
  MOFEM_LOG_CHANNEL("WORLD");                                                  \
  MOFEM_LOG_CHANNEL("SYNC");                                                   \
  MOFEM_LOG_FUNCTION();                                                        \
  MOFEM_LOG_TAG("SYNC", "ProblemCore");                                        \
  MOFEM_LOG_TAG("WORLD", "ProblemCore")

namespace MoFEM {

bool Core::check_problem(const string name) {
  Problem_multiIndex::index<Problem_mi_tag>::type::iterator pit;
  pit = pRoblems.get<Problem_mi_tag>().find(name);
  if (pit == pRoblems.get<Problem_mi_tag>().end()) {
    return false;
  }
  return true;
}

MoFEMErrorCode Core::addProblem(const BitProblemId id, const std::string &name,
                                int verb) {
  ProblemCoreFunctionBegin;

  if (verb == -1)
    verb = verbose;
  EntityHandle meshset;
  CHKERR get_moab().create_meshset(MESHSET_SET, meshset);
  CHKERR get_moab().tag_set_data(th_ProblemId, &meshset, 1, &id);

  // Add problem meshset to partion meshset. In case of no elements
  // on processor part, when mesh file is read, finite element meshset is
  // prevented from deletion by moab reader.
  auto add_meshset_to_partition = [&](auto meshset) {
    MoFEMFunctionBegin;
    const void *tag_vals[] = {&rAnk};
    ParallelComm *pcomm = ParallelComm::get_pcomm(
        &get_moab(), get_basic_entity_data_ptr()->pcommID);
    Tag part_tag = pcomm->part_tag();
    Range tagged_sets;
    CHKERR get_moab().get_entities_by_type_and_tag(0, MBENTITYSET, &part_tag,
                                                   tag_vals, 1, tagged_sets,
                                                   moab::Interface::UNION);
    for (auto s : tagged_sets)
      CHKERR get_moab().add_entities(s, &meshset, 1);
    MoFEMFunctionReturn(0);
  };
  CHKERR add_meshset_to_partition(meshset);

  void const *tag_data[] = {name.c_str()};
  int tag_sizes[1];
  tag_sizes[0] = name.size();
  CHKERR get_moab().tag_set_by_ptr(th_ProblemName, &meshset, 1, tag_data,
                                   tag_sizes);
  // create entry
  auto p = pRoblems.insert(Problem(moab, meshset));
  if (!p.second) {
    MOFEM_LOG_FUNCTION();
    MOFEM_LOG_ATTRIBUTES("SELF", LogManager::BitScope);
    MOFEM_LOG("SELF", Sev::error) << "Following problem can not be added:";
    MOFEM_LOG("SELF", Sev::error)
        << "Problem " << name << " id(" << id.to_ulong() << ") " << id
        << " added meshset " << meshset;
    MOFEM_LOG("SELF", Sev::error) << "List of problems already in databse:";
    for (auto &p : pRoblems) {
      MOFEM_LOG("SELF", Sev::error)
          << p.getName() << " id(" << p.getId().to_ulong() << ") " << id
          << " on meshset " << p.meshset;
    }
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Problem not added");
  }

  MOFEM_LOG("WORLD", Sev::inform) << "Add problem " << name;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_problem(const std::string &name, enum MoFEMTypes bh,
                                 int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  auto miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (miit == pRoblems.get<Problem_mi_tag>().end()) {

    int p_shift = 0;
    for (; pRoblems.get<BitProblemId_mi_tag>().find(BitProblemId().set(
               p_shift)) != pRoblems.get<BitProblemId_mi_tag>().end();
         ++p_shift) {
    }

    auto id = BitProblemId().set(p_shift);
    CHKERR addProblem(id, name, verb);
    
  } else if (bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "problem is in database %s",
             name.c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::delete_problem(const std::string name) {
  MoFEMFunctionBegin;
  auto p_miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (p_miit == pRoblems.get<Problem_mi_tag>().end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "no such problem like < %s >",
             name.c_str());
  }
  const EntityHandle meshset = p_miit->meshset;
  pRoblems.get<Problem_mi_tag>().erase(p_miit);
  CHKERR get_moab().delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

BitProblemId Core::getBitProblemId(const std::string &name) const {
  auto p_miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (p_miit == pRoblems.get<Problem_mi_tag>().end()) {
    THROW_MESSAGE("no such problem like " + name + " >");
  }
  return p_miit->getId();
}

MoFEMErrorCode Core::list_problem() const {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<BitProblemId_mi_tag>::type ProblemById;
  const ProblemById &set_id = pRoblems.get<BitProblemId_mi_tag>();
  ProblemById::iterator miit = set_id.begin();
  for (; miit != set_id.end(); miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscPrintf(mofemComm, ss.str().c_str());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_add_finite_element(const std::string &name_problem,
                                        const std::string &fe_name) {
  MoFEMFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  if (miit == set.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is not there",
             name_problem.c_str());
  }
  BitFEId f_id = getBitFEId(fe_name);
  bool success = set.modify(miit, ProblemFiniteElementChangeBitAdd(f_id));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_problem_unset_finite_element(const std::string &name_problem,
                                          const std::string &fe_name) {
  MoFEMFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  if (miit == set.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is not there",
             name_problem.c_str());
  }
  BitFEId f_id = getBitFEId(fe_name);
  bool success = set.modify(miit, ProblemFiniteElementChangeBitUnSet(f_id));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_problem_ref_level_add_bit(const std::string &name_problem,
                                       const BitRefLevel &bit) {
  MoFEMFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if (miit == set.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is there",
             ss.str().c_str());
  bool success = set.modify(miit, ProblemChangeRefLevelBitAdd(bit));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, 1, "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_problem_ref_level_set_bit(const std::string &name_problem,
                                       const BitRefLevel &bit) {
  MoFEMFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if (miit == set.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is there",
             ss.str().c_str());
  bool success = set.modify(miit, ProblemChangeRefLevelBitSet(bit));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_problem_mask_ref_level_add_bit(const std::string &name_problem,
                                            const BitRefLevel &bit) {
  MoFEMFunctionBegin;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  if (miit == set.end()) {
    std::ostringstream ss;
    ss << name_problem;
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is there",
             ss.str().c_str());
  }
  bool success = set.modify(miit, ProblemChangeRefLevelBitDofMaskAdd(bit));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::modify_problem_mask_ref_level_set_bit(const std::string &name_problem,
                                            const BitRefLevel &bit) {
  MoFEMFunctionBegin typedef Problem_multiIndex::index<Problem_mi_tag>::type
      ProblemsByName;
  ProblemsByName &set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  if (miit == set.end()) {
    std::ostringstream ss;
    ss << name_problem;
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is there",
             ss.str().c_str());
  }
  bool success = set.modify(miit, ProblemChangeRefLevelBitDofMaskSet(bit));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_problem_on_distributed_mesh(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  Problem_multiIndex::iterator p_miit = pRoblems.begin();
  for (; p_miit != pRoblems.end(); p_miit++) {
    CHKERR getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(
        const_cast<Problem *>(&*p_miit), verb);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_problem(const std::string problem_name, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &prob_by_name = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = prob_by_name.find(problem_name);
  if (p_miit == prob_by_name.end()) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "problem < %s > not found, (top tip: check spelling)",
             problem_name.c_str());
  }
  // zero rows
  bool success = prob_by_name.modify(p_miit, ProblemZeroNbRowsChange());
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  // zero cols
  success = prob_by_name.modify(p_miit, ProblemZeroNbColsChange());
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  // clear finite elements
  success =
      prob_by_name.modify(p_miit, ProblemClearNumeredFiniteElementsChange());
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  // clear data structures
  success = prob_by_name.modify(p_miit, ProblemClearSubProblemData());
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  success = prob_by_name.modify(p_miit, ProblemClearComposedProblemData());
  if (!success)
    SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
            "modification unsuccessful");
  if (p_miit->getRowDofsSequence())
    p_miit->getRowDofsSequence()->clear();
  if (p_miit->getColDofsSequence())
    p_miit->getColDofsSequence()->clear();
  if (p_miit->getSubData())
    p_miit->getSubData().reset();
  if (p_miit->getComposedProblemsData())
    p_miit->getComposedProblemsData().reset();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::build_problems(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  if (!((*buildMoFEM) & BUILD_FIELD))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "fields not build");
  if (!((*buildMoFEM) & BUILD_FE))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "FEs not build");
  if (!((*buildMoFEM) & BUILD_ADJ))
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "adjacencies not build");
  // iterate problems
  Problem_multiIndex::iterator p_miit = pRoblems.begin();
  for (; p_miit != pRoblems.end(); p_miit++) {
    Problem *problem_ptr = const_cast<Problem *>(&*p_miit);
    CHKERR getInterface<ProblemsManager>()->buildProblem(problem_ptr, false,
                                                         verb);
  }
  *buildMoFEM |= BUILD_PROBLEM;
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::clear_problems(int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  // iterate problems
  for (auto p_miit = pRoblems.begin(); p_miit != pRoblems.end(); p_miit++)
    CHKERR clear_problem(p_miit->getName(), verb);
  MoFEMFunctionReturn(0);
}

#define SET_BASIC_METHOD(METHOD, PROBLEM_PTR)                                  \
  {                                                                            \
    METHOD.rAnk = rAnk;                                                        \
    METHOD.sIze = sIze;                                                        \
    METHOD.problemPtr = PROBLEM_PTR;                                           \
    METHOD.fieldsPtr = &fIelds;                                                \
    METHOD.refinedEntitiesPtr = &refinedEntities;                              \
    METHOD.entitiesPtr = &entsFields;                                          \
    METHOD.dofsPtr = &dofsField;                                               \
    METHOD.refinedFiniteElementsPtr = &refinedFiniteElements;                  \
    METHOD.finiteElementsPtr = &finiteElements;                                \
    METHOD.finiteElementsEntitiesPtr = &entsFiniteElements;                    \
    METHOD.adjacenciesPtr = &entFEAdjacencies;                                 \
  }

MoFEMErrorCode Core::problem_basic_method_preProcess(const Problem *problem_ptr,
                                                     BasicMethod &method,
                                                     int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  // finite element
  SET_BASIC_METHOD(method, problem_ptr)
  PetscLogEventBegin(MOFEM_EVENT_preProcess, 0, 0, 0, 0);
  CHKERR method.preProcess();
  PetscLogEventEnd(MOFEM_EVENT_preProcess, 0, 0, 0, 0);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::problem_basic_method_preProcess(const std::string &problem_name,
                                      BasicMethod &method, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "problem is not in database %s",
             problem_name.c_str());
  CHKERR problem_basic_method_preProcess(&*p_miit, method, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::problem_basic_method_postProcess(const Problem *problem_ptr,
                                       BasicMethod &method, int verb) {
  MoFEMFunctionBegin;
  SET_BASIC_METHOD(method, problem_ptr)

  PetscLogEventBegin(MOFEM_EVENT_postProcess, 0, 0, 0, 0);
  CHKERR method.postProcess();
  PetscLogEventEnd(MOFEM_EVENT_postProcess, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Core::problem_basic_method_postProcess(const std::string &problem_name,
                                       BasicMethod &method, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;

  // find p_miit
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end())
    SETERRQ1(mofemComm, 1, "problem is not in database %s",
             problem_name.c_str());

  CHKERR problem_basic_method_postProcess(&*p_miit, method, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(
    const Problem *problem_ptr, const std::string &fe_name, FEMethod &method,
    int lower_rank, int upper_rank,
    boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr, MoFEMTypes bh,
    CacheTupleWeakPtr cache_ptr, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  CacheTupleSharedPtr tmp_cache_ptr;
  if (!cache_ptr.use_count()) {
    tmp_cache_ptr = boost::make_shared<CacheTuple>();
    CHKERR cache_problem_entities(problem_ptr->getName(), tmp_cache_ptr);
    method.cacheWeakPtr = tmp_cache_ptr;
  } else {
    method.cacheWeakPtr = cache_ptr;
  }

  if (!fe_ptr)
    fe_ptr = problem_ptr->numeredFiniteElementsPtr;

  auto miit = fe_ptr->get<Composite_Name_And_Part_mi_tag>().lower_bound(
      boost::make_tuple(fe_name, lower_rank));
  auto hi_miit = fe_ptr->get<Composite_Name_And_Part_mi_tag>().upper_bound(
      boost::make_tuple(fe_name, upper_rank));

  if (miit == hi_miit && (bh & MF_EXIST)) {
    if (!check_finite_element(fe_name)) {
      SETERRQ1(mofemComm, MOFEM_NOT_FOUND, "finite element < %s > not found",
               fe_name.c_str());
    }
  }

  method.feName = fe_name;
  method.loopSize =
      std::distance(miit, hi_miit); // Set numbers of elements in the loop
  method.loHiFERank = std::make_pair(lower_rank, upper_rank);

  SET_BASIC_METHOD(method, &*problem_ptr)
  
  PetscLogEventBegin(MOFEM_EVENT_preProcess, 0, 0, 0, 0);
  CHKERR method.preProcess();
  PetscLogEventEnd(MOFEM_EVENT_preProcess, 0, 0, 0, 0);

  PetscLogEventBegin(MOFEM_EVENT_operator, 0, 0, 0, 0);
  for (int nn = 0; miit != hi_miit; miit++, nn++) {

    method.nInTheLoop = nn; // Index of element in the loop
    method.numeredEntFiniteElementPtr = *miit;

    if (method.exeTestHook) {
      if (method.exeTestHook(&method)) {
        CHKERR method();
      }
    } else {
      CHKERR method();
    }

  }
  PetscLogEventEnd(MOFEM_EVENT_operator, 0, 0, 0, 0);

  PetscLogEventBegin(MOFEM_EVENT_postProcess, 0, 0, 0, 0);
  CHKERR method.postProcess();
  PetscLogEventEnd(MOFEM_EVENT_postProcess, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(
    const std::string &problem_name, const std::string &fe_name,
    FEMethod &method,
    boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr, MoFEMTypes bh,
    CacheTupleWeakPtr cache_ptr, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  CHKERR loop_finite_elements(problem_name, fe_name, method, rAnk, rAnk, fe_ptr,
                              bh, cache_ptr, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(
    const std::string &problem_name, const std::string &fe_name,
    FEMethod &method, int lower_rank, int upper_rank,
    boost::shared_ptr<NumeredEntFiniteElement_multiIndex> fe_ptr, MoFEMTypes bh,
    CacheTupleWeakPtr cache_ptr, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;

  auto &prb_by_name = pRoblems.get<Problem_mi_tag>();
  auto p_miit = prb_by_name.find(problem_name);
  if (p_miit == prb_by_name.end())
    SETERRQ1(mofemComm, MOFEM_INVALID_DATA, "Problem <%s> is not in database",
             problem_name.c_str());

  CHKERR loop_finite_elements(&*p_miit, fe_name, method, lower_rank, upper_rank,
                              fe_ptr, bh, cache_ptr, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(const Problem *problem_ptr,
                               const std::string &field_name, RowColData rc,
                               DofMethod &method, int lower_rank,
                               int upper_rank, int verb) {
  MoFEMFunctionBegin;
  SET_BASIC_METHOD(method, &*problem_ptr);
  typedef NumeredDofEntity_multiIndex::index<Unique_mi_tag>::type
      NumeredDofsByUId;
  NumeredDofsByUId *dofs;
  switch (rc) {
  case ROW:
    dofs = &problem_ptr->numeredRowDofsPtr->get<Unique_mi_tag>();
    break;
  case COL:
    dofs = &problem_ptr->numeredColDofsPtr->get<Unique_mi_tag>();
    break;
  default:
    SETERRQ(mofemComm, MOFEM_DATA_INCONSISTENCY, "Not implemented");
  }

  auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
    method.fieldPtr = *field_it;
  } else {
    SETERRQ(mofemComm, MOFEM_NOT_FOUND,
            ("Field not found " + field_name).c_str());
  }

  auto miit = dofs->lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_miit = dofs->upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));

  method.loopSize = std::distance(miit, hi_miit);
  method.loHiFERank = std::make_pair(lower_rank, upper_rank);      

  CHKERR method.preProcess();

  int nn = 0;
  for (; miit != hi_miit; miit++, nn++) {
    if ((*miit)->getPart() >= lower_rank && (*miit)->getPart() <= upper_rank) {
      method.nInTheLoop = nn; // Index of element in the loop
      method.dofPtr = miit->get()->getDofEntityPtr();
      method.dofNumeredPtr = *miit;
      CHKERR method();
    }
  }

  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(
    const std::string &problem_name, const std::string &field_name,
    RowColData rc,     // ROW or COL
    DofMethod &method, // Finite element instance processed on each DOF
    int lower_rank,    // Only DOFs on processor higher or equal to this are
                       // processed
    int upper_rank,    // Only DOFs lowest or higher to this are processed
    int verb           // verbosity level
) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "problem not in database %s",
             problem_name.c_str());
  CHKERR loop_dofs(&*p_miit, field_name, rc, method, lower_rank, upper_rank,
                   verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(const std::string &problem_name,
                               const std::string &field_name, RowColData rc,
                               DofMethod &method, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  CHKERR loop_dofs(problem_name, field_name, rc, method, 0, sIze, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(const std::string &field_name, DofMethod &method,
                               int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  SET_BASIC_METHOD(method, nullptr);

  auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
    method.fieldPtr = *field_it;
  } else {
    SETERRQ(mofemComm, MOFEM_NOT_FOUND,
            ("Field not found " + field_name).c_str());
  }

  auto miit = dofsField.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_miit = dofsField.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));

  method.loopSize = std::distance(miit, hi_miit);
  method.loHiFERank = std::make_pair(0, get_comm_size());

  CHKERR method.preProcess();
  for (int nn = 0; miit != hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr = *miit;
    CHKERR method();
  }
  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_entities(const Problem *problem_ptr,
                                   const std::string field_name, RowColData rc,
                                   EntityMethod &method, int lower_rank,
                                   int upper_rank, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  decltype(problem_ptr->numeredRowDofsPtr) dofs;
  switch (rc) {
  case ROW:
    dofs = problem_ptr->numeredRowDofsPtr;
    break;
  case COL:
    dofs = problem_ptr->numeredColDofsPtr;
    break;
  default:
    SETERRQ(mofemComm, MOFEM_DATA_INCONSISTENCY,
            "It works only with rows or columns");
  }

  auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
    method.fieldPtr = *field_it;
  } else {
    SETERRQ(mofemComm, MOFEM_NOT_FOUND,
            ("Field not found " + field_name).c_str());
  }

  auto miit = dofs->lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_miit = dofs->upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));

  using FieldEntity_view_multiIndex = multi_index_container<

      boost::shared_ptr<FieldEntity>,
      indexed_by<

          ordered_unique<

              tag<Ent_mi_tag>,
              const_mem_fun<FieldEntity::interface_type_RefEntity, EntityHandle,
                            &FieldEntity::getEnt>>

          >>;

  FieldEntity_view_multiIndex ents_view;
  auto hint = ents_view.begin();
  for (; miit != hi_miit; ++miit)
    if ((*miit)->getPart() >= lower_rank && (*miit)->getPart() <= upper_rank)
      ents_view.emplace_hint(hint, (*miit)->getFieldEntityPtr());

  SET_BASIC_METHOD(method, problem_ptr);

  method.loopSize = ents_view.size();
  method.loHiFERank = std::make_pair(lower_rank, upper_rank);

  CHKERR method.preProcess();
  method.nInTheLoop = 0;
  for (auto &field_ent : ents_view) {
    method.entPtr = field_ent;
    CHKERR method();
    ++method.nInTheLoop;
  }
  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_entities(const std::string problem_name,
                                   const std::string field_name, RowColData rc,
                                   EntityMethod &method, int lower_rank,
                                   int upper_rank, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  // find p_miit
  auto &prb = pRoblems.get<Problem_mi_tag>();
  auto p_miit = prb.find(problem_name);
  if (p_miit == prb.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "problem not in database %s",
             problem_name.c_str());
  CHKERR loop_entities(&*p_miit, field_name, rc, method, lower_rank, upper_rank,
                       verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_entities(const std::string problem_name,
                                   const std::string field_name, RowColData rc,
                                   EntityMethod &method, int verb) {
  return loop_entities(problem_name, field_name, rc, method, rAnk, rAnk, verb);
}

MoFEMErrorCode Core::loop_entities(const std::string field_name,
                                   EntityMethod &method,
                                   Range const *const ents, int verb) {
  MoFEMFunctionBegin;
  if (verb == DEFAULT_VERBOSITY)
    verb = verbose;
  SET_BASIC_METHOD(method, nullptr);

  auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
  if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
    method.fieldPtr = *field_it;
  } else {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "Field not found %s",
             field_name.c_str());
  }

  auto lo_eit = entsFields.get<Unique_mi_tag>().lower_bound(
      FieldEntity::getLoBitNumberUId((*field_it)->getBitNumber()));
  auto hi_eit = entsFields.get<Unique_mi_tag>().upper_bound(
      FieldEntity::getHiBitNumberUId((*field_it)->getBitNumber()));

  typedef multi_index_container<
      boost::shared_ptr<FieldEntity>,
      indexed_by<ordered_unique<
          tag<Ent_mi_tag>, const_mem_fun<FieldEntity::interface_RefEntity,
                                         EntityHandle, &FieldEntity::getEnt>>>>
      FieldEntity_view_multiIndex;

  FieldEntity_view_multiIndex ents_view;
  ents_view.insert(lo_eit, hi_eit);

  method.loopSize = ents_view.size();
  method.loHiFERank = std::make_pair(0, get_comm_size());

  CHKERR method.preProcess();
  method.nInTheLoop = 0;

  if (ents)
    for (auto p = ents->const_pair_begin(); p != ents->const_pair_end(); ++p)
      for (auto feit = ents_view.lower_bound(p->first);
           feit != ents_view.upper_bound(p->second); ++feit) {
        method.entPtr = *feit;
        CHKERR method();
        ++method.nInTheLoop;
      }
  else
    for (auto &field_ent : ents_view) {
      method.entPtr = field_ent;
      CHKERR method();
      ++method.nInTheLoop;
    }

  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::cache_problem_entities(const std::string prb_name,
                                            CacheTupleWeakPtr cache_weak_ptr) {
  ProblemCoreFunctionBegin;

  if (auto cache_ptr = cache_weak_ptr.lock()) {
    auto p_miit = pRoblems.get<Problem_mi_tag>().find(prb_name);
    if (p_miit == pRoblems.get<Problem_mi_tag>().end())
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "problem not in database %s",
               prb_name.c_str());

    const BitRefLevel &prb_bit = p_miit->getBitRefLevel();
    const BitRefLevel &prb_mask = p_miit->getBitRefLevelMask();
    const BitFEId &prb_fe_id = p_miit->getBitFEId();

    auto &row_dofs = p_miit->numeredRowDofsPtr;
    auto &col_dofs = p_miit->numeredColDofsPtr;

    auto &cache_data = std::get<0>(*cache_ptr);
    auto &cache_row = std::get<1>(*cache_ptr);
    auto &cache_col = std::get<2>(*cache_ptr);

    cache_row.resize(entsFields.size());
    if (row_dofs != col_dofs)
      cache_col.resize(entsFields.size());
    cache_data.resize(entsFields.size());

    size_t idx = 0;
    for (auto it = entsFields.begin(); it != entsFields.end(); ++it, ++idx) {

      const auto uid = (*it)->getLocalUniqueId();
      auto r = entFEAdjacencies.get<Unique_mi_tag>().equal_range(uid);
      for (auto lo = r.first; lo != r.second; ++lo) {

        if ((lo->getBitFEId() & prb_fe_id).any()) {

          const BitRefLevel &fe_bit = lo->entFePtr->getBitRefLevel();

          // if entity is not problem refinement level
          if (((fe_bit & prb_mask) != fe_bit) || ((fe_bit & prb_bit).none()))
            continue;

          auto cache_numered_dofs = [&](auto &numered_dofs, auto &cache_vec,
                                        auto &ent_cache) {
            auto dit = numered_dofs->lower_bound(uid);
            decltype(dit) hi_dit;
            if (dit != numered_dofs->end())
              hi_dit = numered_dofs->upper_bound(
                  uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
            else
              hi_dit = dit;

            ent_cache = boost::shared_ptr<EntityCacheNumeredDofs>(
                cache_ptr, &(cache_vec[idx]));
            cache_vec[idx].loHi = {dit, hi_dit};
          };

          auto cache_dofs = [&](auto &dofs, auto &cache_vec, auto &ent_cache) {
            auto dit = dofs.lower_bound(uid);
            decltype(dit) hi_dit;
            if (dit != dofs.end())
              hi_dit = dofs.upper_bound(
                  uid | static_cast<UId>(MAX_DOFS_ON_ENTITY - 1));
            else
              hi_dit = dit;

            ent_cache = boost::shared_ptr<EntityCacheDofs>(cache_ptr,
                                                           &(cache_vec[idx]));
            cache_vec[idx].loHi = {dit, hi_dit};
          };

          cache_numered_dofs(row_dofs, cache_row, (*it)->entityCacheRowDofs);
          if (row_dofs != col_dofs) {
            if (cache_col.size() != entsFields.size())
              cache_col.resize(entsFields.size());
            cache_numered_dofs(col_dofs, cache_col, (*it)->entityCacheColDofs);
          } else {
            (*it)->entityCacheColDofs = (*it)->entityCacheRowDofs;
          }

          cache_dofs(dofsField, cache_data, (*it)->entityCacheDataDofs);

          break;
        }
      }
    }
  } else {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Cache not allocated");
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
