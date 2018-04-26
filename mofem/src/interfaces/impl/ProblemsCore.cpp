/** \file ProblemsCore.cpp
 * \brief Managing complexities for problem
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
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  EntityHandle meshset;
  CHKERR moab.create_meshset(MESHSET_SET | MESHSET_TRACK_OWNER, meshset);
  CHKERR moab.tag_set_data(th_ProblemId, &meshset, 1, &id);
  void const *tag_data[] = {name.c_str()};
  int tag_sizes[1];
  tag_sizes[0] = name.size();
  CHKERR moab.tag_set_by_ptr(th_ProblemName, &meshset, 1, tag_data, tag_sizes);
  // create entry
  std::pair<Problem_multiIndex::iterator, bool> p =
      pRoblems.insert(Problem(moab, meshset));
  NOT_USED(p);
  assert(p.second);
  if (verb > 0) {
    std::ostringstream ss;
    ss << "add problem: " << name << std::endl;
    PetscPrintf(cOmm, ss.str().c_str());
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::add_problem(const std::string &name, enum MoFEMTypes bh,
                                 int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  auto miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (miit == pRoblems.get<Problem_mi_tag>().end()) {
    BitProblemId id = getProblemShift();
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
  CHKERR moab.delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

BitProblemId Core::getBitProblemId(const std::string &name) const {
  auto p_miit = pRoblems.get<Problem_mi_tag>().find(name);
  if (p_miit == pRoblems.get<Problem_mi_tag>().end()) {
    THROW_MESSAGE(
      "no such problem like " + name + " >"
    );
  }
  return p_miit->getId();
}

MoFEMErrorCode Core::list_problem() const {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<BitProblemId_mi_tag>::type ProblemById;
  const ProblemById &set_id  = pRoblems.get<BitProblemId_mi_tag>();
  ProblemById::iterator miit = set_id.begin();
  for (; miit != set_id.end(); miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscPrintf(cOmm, ss.str().c_str());
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_add_finite_element(const std::string &name_problem,
                                        const std::string &fe_name) {
  MoFEMFunctionBeginHot;
  try {
    typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
    ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
    ProblemsByName::iterator miit = set.find(name_problem);
    if (miit == set.end()) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
               "this problem <%s> is not there", name_problem.c_str());
    }
    BitFEId f_id = getBitFEId(fe_name);
    bool success = set.modify(miit, ProblemFiniteElementChangeBitAdd(f_id));
    if (!success)
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_unset_finite_element(const std::string &name_problem,
                                          const std::string &fe_name) {
  MoFEMFunctionBeginHot;
  try {
    typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
    ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
    ProblemsByName::iterator miit = set.find(name_problem);
    if (miit == set.end()) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND,
               "this problem <%s> is not there", name_problem.c_str());
    }
    BitFEId f_id = getBitFEId(fe_name);
    bool success = set.modify(miit, ProblemFiniteElementChangeBitUnSet(f_id));
    if (!success)
      SETERRQ(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
              "modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_ref_level_add_bit(const std::string &name_problem,
                                       const BitRefLevel &bit) {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if (miit == set.end())
    SETERRQ1(PETSC_COMM_SELF, MOFEM_NOT_FOUND, "this problem <%s> is there",
             ss.str().c_str());
  bool success = set.modify(miit, ProblemChangeRefLevelBitAdd(bit));
  if (!success)
    SETERRQ(PETSC_COMM_SELF, 1, "modification unsuccessful");
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_ref_level_set_bit(const std::string &name_problem,
                                       const BitRefLevel &bit) {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_mask_ref_level_add_bit(const std::string &name_problem,
                                            const BitRefLevel &bit) {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
Core::modify_problem_mask_ref_level_set_bit(const std::string &name_problem,
                                            const BitRefLevel &bit) {
  MoFEMFunctionBeginHot;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &set           = pRoblems.get<Problem_mi_tag>();
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_problem_on_distributed_mesh(int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  DofEntity_multiIndex_active_view dofs_rows, dofs_cols;
  Problem_multiIndex::iterator p_miit = pRoblems.begin();
  for (; p_miit != pRoblems.end(); p_miit++) {
    ierr = getInterface<ProblemsManager>()->buildProblemOnDistributedMesh(
        const_cast<Problem *>(&*p_miit), verb);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_problem(const std::string &problem_name, int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  ProblemsByName &prob_by_name    = pRoblems.get<Problem_mi_tag>();
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
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::build_problems(int verb) {
  MoFEMFunctionBeginHot;
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
    ierr                 = build_problem(problem_ptr, false, verb);
    CHKERRG(ierr);
  }
  *buildMoFEM |= BUILD_PROBLEM;
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Core::clear_problems(int verb) {
  MoFEMFunctionBeginHot;
  if (verb == -1)
    verb = verbose;
  // iterate problems
  for (Problem_multiIndex::iterator p_miit = pRoblems.begin();
       p_miit != pRoblems.end(); p_miit++) {
    ierr = clear_problem(p_miit->getName(), verb);
    CHKERRG(ierr);
  }
  MoFEMFunctionReturnHot(0);
}

#define SET_BASIC_METHOD(METHOD, PROBLEM_PTR)                                  \
  {                                                                            \
    METHOD.rAnk                      = rAnk;                                   \
    METHOD.sIze                      = sIze;                                   \
    METHOD.problemPtr                = PROBLEM_PTR;                            \
    METHOD.fieldsPtr                 = &fIelds;                                \
    METHOD.refinedEntitiesPtr        = &refinedEntities;                       \
    METHOD.entitiesPtr               = &entsFields;                            \
    METHOD.dofsPtr                   = &dofsField;                             \
    METHOD.refinedFiniteElementsPtr  = &refinedFiniteElements;                 \
    METHOD.finiteElementsPtr         = &finiteElements;                        \
    METHOD.finiteElementsEntitiesPtr = &entsFiniteElements;                    \
    METHOD.adjacenciesPtr            = &entFEAdjacencies;                      \
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
  ProblemsByName &pRoblems_set    = pRoblems.get<Problem_mi_tag>();
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
  ProblemsByName &pRoblems_set    = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end())
    SETERRQ1(cOmm, 1, "problem is not in database %s", problem_name.c_str());

  CHKERR problem_basic_method_postProcess(&*p_miit, method, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(const std::string &problem_name,
                                          const std::string &fe_name,
                                          FEMethod &method, MoFEMTypes bh,
                                          int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  CHKERR loop_finite_elements(problem_name, fe_name, method, rAnk, rAnk, bh,
                              verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(
    const Problem *problem_ptr, const std::string &fe_name,
    FEMethod &method, // reference to finite element implementation
    int lower_rank,   // only elements on part between low and up rank are
                      // processed
    int upper_rank,
    MoFEMTypes bh, // is set to MF_EXIST, throw error if element is not declared
                   // in database
    int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;

  method.feName = fe_name;
  SET_BASIC_METHOD(method, &*problem_ptr)
  PetscLogEventBegin(MOFEM_EVENT_preProcess, 0, 0, 0, 0);
  CHKERR method.preProcess();
  PetscLogEventEnd(MOFEM_EVENT_preProcess, 0, 0, 0, 0);

  NumeredEntFiniteElementbyNameAndPart &numered_fe =
      problem_ptr->numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>();
  NumeredEntFiniteElementbyNameAndPart::iterator miit =
      numered_fe.lower_bound(boost::make_tuple(fe_name, lower_rank));
  NumeredEntFiniteElementbyNameAndPart::iterator hi_miit =
      numered_fe.upper_bound(boost::make_tuple(fe_name, upper_rank));

  if (miit == hi_miit && (bh & MF_EXIST)) {
    if (!check_finite_element(fe_name)) {
      SETERRQ1(cOmm, MOFEM_NOT_FOUND, "finite element < %s > not found",
               fe_name.c_str());
    }
  }

  method.loopSize = distance(miit, hi_miit);
  for (int nn = 0; miit != hi_miit; miit++, nn++) {

    // back_miit--;

    method.nInTheLoop                 = nn;
    method.numeredEntFiniteElementPtr = *miit;
    method.dataPtr                    = (*miit)->sPtr->data_dofs;
    method.rowPtr                     = (*miit)->rows_dofs;
    method.colPtr                     = (*miit)->cols_dofs;

    PetscLogEventBegin(MOFEM_EVENT_operator, 0, 0, 0, 0);
    CHKERR method();
    PetscLogEventEnd(MOFEM_EVENT_operator, 0, 0, 0, 0);

  }

  PetscLogEventBegin(MOFEM_EVENT_postProcess, 0, 0, 0, 0);
  CHKERR method.postProcess();
  PetscLogEventEnd(MOFEM_EVENT_postProcess, 0, 0, 0, 0);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_finite_elements(const std::string &problem_name,
                                          const std::string &fe_name,
                                          FEMethod &method, int lower_rank,
                                          int upper_rank, MoFEMTypes bh,
                                          int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &pRoblems_set    = pRoblems.get<Problem_mi_tag>();
  ProblemsByName::iterator p_miit = pRoblems_set.find(problem_name);
  if (p_miit == pRoblems_set.end())
    SETERRQ1(cOmm, 1, "problem is not in database %s", problem_name.c_str());

  CHKERR loop_finite_elements(&*p_miit, fe_name, method, lower_rank, upper_rank,
                              bh, verb);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(const Problem *problem_ptr,
                               const std::string &field_name, RowColData rc,
                               DofMethod &method, int lower_rank,
                               int upper_rank, int verb) {
  MoFEMFunctionBegin;
  SET_BASIC_METHOD(method, &*problem_ptr);
  typedef NumeredDofEntity_multiIndex::index<
      Composite_Name_And_Part_mi_tag>::type NumeredDofsByNameAndPart;
  NumeredDofsByNameAndPart *dofs;
  switch (rc) {
  case ROW:
    dofs = &problem_ptr->numeredDofsRows->get<Composite_Name_And_Part_mi_tag>();
    break;
  case COL:
    dofs = &problem_ptr->numeredDofsCols->get<Composite_Name_And_Part_mi_tag>();
    break;
  default:
    SETERRQ(cOmm, MOFEM_DATA_INCONSISTENCY, "not implemented");
  }
  NumeredDofsByNameAndPart::iterator miit =
      dofs->lower_bound(boost::make_tuple(field_name, lower_rank));
  NumeredDofsByNameAndPart::iterator hi_miit =
      dofs->upper_bound(boost::make_tuple(field_name, upper_rank));
  if (miit != hi_miit) {
    method.fieldPtr = miit->get()->getFieldPtr();
  } else {
    Field_multiIndex::index<FieldName_mi_tag>::type::iterator field_it =
        fIelds.get<FieldName_mi_tag>().find(field_name);
    if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
      method.fieldPtr = *field_it;
    }
  }
  CHKERR method.preProcess();
  for (; miit != hi_miit; miit++) {
    method.dofPtr        = miit->get()->getDofEntityPtr();
    method.dofNumeredPtr = *miit;
    CHKERR method();
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
  if (verb == -1)
    verb = verbose;
  typedef Problem_multiIndex::index<Problem_mi_tag>::type ProblemsByName;
  // find p_miit
  ProblemsByName &pRoblems_set    = pRoblems.get<Problem_mi_tag>();
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
  if (verb == -1)
    verb = verbose;
  CHKERR loop_dofs(problem_name, field_name, rc, method, 0, sIze, verb);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_dofs(const std::string &field_name, DofMethod &method,
                               int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  SET_BASIC_METHOD(method, NULL);
  auto miit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  auto hi_miit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  if (miit != hi_miit) {
    method.fieldPtr = miit->get()->getFieldPtr();
  } else {
    auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
    if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
      method.fieldPtr = *field_it;
    }
  }
  method.loopSize = distance(miit, hi_miit);
  CHKERR method.preProcess();
  for (int nn = 0; miit != hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr     = *miit;
    CHKERR method();
  }
  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Core::loop_entities(const std::string &field_name,
                                   EntityMethod &method, int verb) {
  MoFEMFunctionBegin;
  if (verb == -1)
    verb = verbose;
  SET_BASIC_METHOD(method, NULL);
  auto miit = entsFields.get<FieldName_mi_tag>().lower_bound(field_name);
  auto hi_miit = entsFields.get<FieldName_mi_tag>().upper_bound(field_name);
  if (miit != hi_miit) {
    method.fieldPtr = miit->get()->getFieldPtr();
  } else {
    auto field_it = fIelds.get<FieldName_mi_tag>().find(field_name);
    if (field_it != fIelds.get<FieldName_mi_tag>().end()) {
      method.fieldPtr = *field_it;
    }
  }
  method.loopSize = distance(miit, hi_miit);
  CHKERR method.preProcess();
  for (int nn = 0; miit != hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr     = *miit;
    CHKERR method();
  }
  CHKERR method.postProcess();
  MoFEMFunctionReturn(0);
}



}
