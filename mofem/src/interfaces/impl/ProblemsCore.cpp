/** \file ProblemCore.cpp
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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

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

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <MeshRefinement.hpp>
#include <PrismInterface.hpp>
#include <SeriesRecorder.hpp>
#include <Core.hpp>

#include <boost/scoped_array.hpp>

#include <moab/MeshTopoUtil.hpp>

namespace MoFEM {

struct __attribute__ ((__packed__)) IdxDataType {
  int globalDof;
  char uId[sizeof(UId)];
  IdxDataType(const GlobalUId &uid,int global_dof):
    globalDof(global_dof) {
    bcopy(&uid,uId,sizeof(UId));
  }
};

bool Core::check_problem(const string name) {
  MoFEMProblem_multiIndex::index<Problem_mi_tag>::type::iterator pit;
  pit = pRoblems.get<Problem_mi_tag>().find(name);
  if(pit==pRoblems.get<Problem_mi_tag>().end()) {
    return false;
  }
  return true;
}

PetscErrorCode Core::add_problem(const BitProblemId id,const std::string& name) {
  PetscFunctionBegin;
  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset); CHKERRQ_MOAB(rval);
  rval = moab.tag_set_data(th_ProblemId,&meshset,1,&id); CHKERRQ_MOAB(rval);
  void const* tag_data[] = { name.c_str() };
  int tag_sizes[1]; tag_sizes[0] = name.size();
  rval = moab.tag_set_by_ptr(th_ProblemName,&meshset,1,tag_data,tag_sizes); CHKERRQ_MOAB(rval);
  //create entry
  std::pair<MoFEMProblem_multiIndex::iterator,bool> p = pRoblems.insert(MoFEMProblem(moab,meshset));
  NOT_USED(p);
  assert(p.second);
  if(verbose>0) {
    std::ostringstream ss;
    ss << "add problem: " << name << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::add_problem(const std::string& name,enum MoFEMTypes bh,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  if(miit==set.end()) {
    BitProblemId id = getProblemShift();
    ierr = add_problem(id,name); CHKERRQ(ierr);
  } else if(bh == MF_EXCL) {
    SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"problem is in database %s",name.c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::delete_problem(const std::string name) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name &mofem_problems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = mofem_problems_set.find(name);
  if(p_miit == mofem_problems_set.end()) {
    SETERRQ1(PETSC_COMM_SELF,1,"no such problem like < %s >",name.c_str());
  }
  EntityHandle meshset = p_miit->meshset;
  mofem_problems_set.erase(p_miit);
  rval = moab.delete_entities(&meshset,1); CHKERRQ_MOAB(rval);
  PetscFunctionReturn(0);
}

BitProblemId Core::get_BitProblemId(const std::string& name) const {
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  const mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name);
  return miit->getId();
}

PetscErrorCode Core::list_problem() const {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<BitProblemId_mi_tag>::type problem_set_by_id;
  const problem_set_by_id &set_id = pRoblems.get<BitProblemId_mi_tag>();
  problem_set_by_id::iterator miit = set_id.begin();
  for(;miit!=set_id.end();miit++) {
    std::ostringstream ss;
    ss << *miit << std::endl;
    PetscPrintf(comm,ss.str().c_str());
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::modify_problem_add_finite_element(const std::string &name_problem,const std::string &fe_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
    mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
    mofem_problems_by_name::iterator miit = set.find(name_problem);
    if(miit==set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",name_problem.c_str());
    }
    BitFEId f_id = getBitFEId(fe_name);
    bool success = set.modify(miit,ProblemFiniteElementChangeBitAdd(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::modify_problem_unset_finite_element(const std::string &name_problem,const std::string &fe_name) {
  PetscFunctionBegin;
  try {
    typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
    mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
    mofem_problems_by_name::iterator miit = set.find(name_problem);
    if(miit==set.end()) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is not there",name_problem.c_str());
    }
    BitFEId f_id = getBitFEId(fe_name);
    bool success = set.modify(miit,ProblemFiniteElementChangeBitUnSet(f_id));
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode Core::modify_problem_ref_level_add_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitAdd(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,1,"modification unsuccessful");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::modify_problem_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitSet(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::modify_problem_mask_ref_level_set_bit(const std::string &name_problem,const BitRefLevel &bit) {
  PetscFunctionBegin;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  mofem_problems_by_name& set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator miit = set.find(name_problem);
  std::ostringstream ss;
  ss << name_problem;
  if(miit==set.end()) SETERRQ1(PETSC_COMM_SELF,MOFEM_NOT_FOUND,"this problem <%s> is there",ss.str().c_str());
  bool success = set.modify(miit,ProblemChangeRefLevelBitDofMaskSet(bit));
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  PetscFunctionReturn(0);
}

PetscErrorCode Core::build_problem_on_distributed_mesh(
  const std::string &name,const bool square_matrix,int verb
) {
  ProblemsManager *problems_manager_ptr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = query_interface(problems_manager_ptr); CHKERRQ(ierr);
  ierr = problems_manager_ptr->buildProblemOnDistributedMesh(
    name,square_matrix,verb
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem_on_distributed_mesh(
  MoFEMProblem *problem_ptr,const bool square_matrix,int verb
) {
  ProblemsManager *problems_manager_ptr;
  PetscFunctionBegin;
  ierr = query_interface(problems_manager_ptr); CHKERRQ(ierr);
  ierr = problems_manager_ptr->buildProblemOnDistributedMesh(
    problem_ptr,square_matrix,verb
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem_on_distributed_mesh(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  DofEntity_multiIndex_active_view dofs_rows,dofs_cols;
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  for(;p_miit!=pRoblems.end();p_miit++) {
    ierr = build_problem_on_distributed_mesh(const_cast<MoFEMProblem*>(&*p_miit),verb); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_mesh(
  const Range &ents,const int dim,const int adj_dim,const int n_parts,int verb
) {
  ProblemsManager *prb_mng_ptr;
  PetscFunctionBegin;
  ierr = query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionMesh(ents,dim,adj_dim,n_parts,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problem(const std::string &problem_name,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type pRoblemsByName;
  pRoblemsByName &prob_by_name = pRoblems.get<Problem_mi_tag>();
  pRoblemsByName::iterator p_miit = prob_by_name.find(problem_name);
  if(p_miit == prob_by_name.end()) {
    SETERRQ1(
      PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,
      "problem < %s > not found, (top tip: check spelling)",problem_name.c_str()
    );
  }
  //zero rows
  bool success = prob_by_name.modify(p_miit,ProblemZeroNbRowsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //zero cols
  success = prob_by_name.modify(p_miit,ProblemZeroNbColsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  //clear finite elements
  success = prob_by_name.modify(p_miit,ProblemClearNumeredFiniteElementsChange());
  if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");

  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem(MoFEMProblem *problem_ptr,const bool square_matrix,int verb) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  // Note: Only allowe changes on problem_ptr structure which not influence multindex
  // indexing are allowd.
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->buildProblem(problem_ptr,square_matrix,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problem(const std::string &problem_name,const bool square_matrix,int verb) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  // Note: Only allowe changes on problem_ptr structure which not influence multindex
  // indexing are allowd.
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->buildProblem(problem_name,square_matrix,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::build_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  if(!((*buildMoFEM)&BUILD_FIELD)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"fields not build");
  if(!((*buildMoFEM)&BUILD_FE)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"FEs not build");
  if(!((*buildMoFEM)&BUILD_ADJ)) SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"adjacencies not build");
  //iterate problems
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  for(;p_miit!=pRoblems.end();p_miit++) {
    MoFEMProblem *problem_ptr =  const_cast<MoFEMProblem*>(&*p_miit);
    ierr = build_problem(problem_ptr,false,verb); CHKERRQ(ierr);
  }
  *buildMoFEM |= BUILD_PROBLEM;
  PetscFunctionReturn(0);
}
PetscErrorCode Core::clear_problems(int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  MoFEMProblem_multiIndex::iterator p_miit = pRoblems.begin();
  //iterate problems
  for(;p_miit!=pRoblems.end();p_miit++) {
    //zero rows
    bool success = pRoblems.modify(p_miit,ProblemZeroNbRowsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //zero cols
    success = pRoblems.modify(p_miit,ProblemZeroNbColsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
    //clear finite elements
    success = pRoblems.modify(p_miit,ProblemClearNumeredFiniteElementsChange());
    if(!success) SETERRQ(PETSC_COMM_SELF,MOFEM_OPERATION_UNSUCCESSFUL,"modification unsuccessful");
  }
  PetscFunctionReturn(0);
}
PetscErrorCode Core::partition_simple_problem(const std::string &name,int verb) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->partitionSimpleProblem(name,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_compose_problem(
  const std::string &name,
  const std::string &problem_for_rows,
  bool copy_rows,
  const std::string &problem_for_cols,
  bool copy_cols,
  int verb
) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->inheretPartition(
    name,problem_for_rows,copy_rows,problem_for_cols,copy_cols,verb
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::build_sub_problem(
  const std::string &out_name,
  const std::vector<std::string> &fields_row,
  const std::vector<std::string> &fields_col,
  const std::string &main_problem,
  const bool square_matrix,
  int verb
) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);

  if(verb==-1) verb = verbose;
  ierr = problem_manager_ptr->buildSubProblem(
    out_name,fields_row,fields_col,main_problem,square_matrix,verb
  ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode Core::printPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->printPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::debugPartitionedProblem(const MoFEMProblem *problem_ptr,int verb) {
  ProblemsManager *problem_manager_ptr;
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = query_interface(problem_manager_ptr); CHKERRQ(ierr);
  ierr = problem_manager_ptr->debugPartitionedProblem(problem_ptr,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_problem(const std::string &name,int verb) {
  ProblemsManager *prb_mng_ptr;
  PetscFunctionBegin;
  ierr = query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionProblem(name,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_finite_elements(
  const std::string &name,
  bool part_from_moab,
  int low_proc,
  int hi_proc,
  int verb
) {
  PetscFunctionBegin;
  ProblemsManager *prb_mng_ptr;
  PetscFunctionBegin;
  ierr = query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements(
    name,part_from_moab,low_proc,hi_proc,verb
  ); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::partition_ghost_dofs(const std::string &name,int verb) {
  PetscFunctionBegin;
  ProblemsManager *prb_mng_ptr;
  PetscFunctionBegin;
  ierr = query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionGhostDofs(name,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#define SET_BASIC_METHOD(METHOD,PROBLEM_PTR) \
  { \
    METHOD.rAnk = rAnk; \
    METHOD.sIze = sIze; \
    METHOD.problemPtr = PROBLEM_PTR; \
    METHOD.fieldsPtr = &fIelds; \
    METHOD.refinedEntitiesPtr = &refinedEntities; \
    METHOD.entitiesPtr = &entsFields; \
    METHOD.dofsPtr = &dofsField; \
    METHOD.refinedFiniteElementsPtr = &refinedFiniteElements; \
    METHOD.finiteElementsPtr = &finiteElements; \
    METHOD.finiteElementsEntitiesPtr = &entsFiniteElements; \
    METHOD.adjacenciesPtr = &entFEAdjacencies; \
  }

PetscErrorCode Core::problem_basic_method_preProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element
  SET_BASIC_METHOD(method,problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);
  PetscFunctionReturn(0);
}

PetscErrorCode Core::problem_basic_method_preProcess(const std::string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());
  ierr = problem_basic_method_preProcess(&*p_miit,method,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const MoFEMProblem *problem_ptr,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(method,problem_ptr)

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::problem_basic_method_postProcess(const std::string &problem_name,BasicMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;

  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());

  ierr = problem_basic_method_postProcess(&*p_miit,method,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const std::string &problem_name,
  const std::string &fe_name,
  FEMethod &method,MoFEMTypes bh,
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;

  ierr = loop_finite_elements(problem_name,fe_name,method,rAnk,rAnk,bh,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const MoFEMProblem *problem_ptr,
  const std::string &fe_name,
  FEMethod &method, // reference to finite element implementation
  int lower_rank, // only elements on part between low and up rank are processed
  int upper_rank,
  MoFEMTypes bh, // is set to MF_EXIST, throw error if element is not declared in databse
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  // finite element

  method.feName = fe_name;
  SET_BASIC_METHOD(method,&*problem_ptr)
  PetscLogEventBegin(USER_EVENT_preProcess,0,0,0,0);
  ierr = method.preProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_preProcess,0,0,0,0);

  NumeredEntFiniteElementbyNameAndPart &numered_fe =
  problem_ptr->numeredFiniteElements.get<Composite_Name_And_Part_mi_tag>();
  NumeredEntFiniteElementbyNameAndPart::iterator miit =
  numered_fe.lower_bound(boost::make_tuple(fe_name,lower_rank));
  NumeredEntFiniteElementbyNameAndPart::iterator hi_miit =
  numered_fe.upper_bound(boost::make_tuple(fe_name,upper_rank));

  if(miit==hi_miit && bh&MF_EXIST) {
    if(!check_finite_element(fe_name)) {
      SETERRQ1(comm,MOFEM_NOT_FOUND,"finite element < %s > not found",fe_name.c_str());
    }
  }

  // NumeredEntFiniteElementbyNameAndPart::iterator back_miit = hi_miit;
  method.loopSize = distance(miit,hi_miit);
  for(int nn = 0;miit!=hi_miit;miit++,nn++) {

    // back_miit--;

    method.nInTheLoop = nn;
    method.numeredEntFiniteElementPtr = &*(*miit);
    method.dataPtr = &((*miit)->sPtr->data_dofs);
    method.rowPtr = (*miit)->rows_dofs.get();
    method.colPtr = (*miit)->cols_dofs.get();

    try {
      PetscLogEventBegin(USER_EVENT_operator,0,0,0,0);
      ierr = method(); CHKERRQ(ierr);
      PetscLogEventEnd(USER_EVENT_operator,0,0,0,0);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "FE method " << typeid(method).name() //boost::core::demangle(typeid(method).name())
      << "   throw in method: " << ex.what()
      << " at line " << __LINE__
      << " in file " << __FILE__ << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }

  }

  PetscLogEventBegin(USER_EVENT_postProcess,0,0,0,0);
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscLogEventEnd(USER_EVENT_postProcess,0,0,0,0);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_finite_elements(
  const std::string &problem_name,
  const std::string &fe_name,
  FEMethod &method,
  int lower_rank,
  int upper_rank,
  MoFEMTypes bh,int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem is not in database %s",problem_name.c_str());

  ierr = loop_finite_elements(&*p_miit,fe_name,method,lower_rank,upper_rank,bh,verb); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const MoFEMProblem *problem_ptr,const std::string &field_name,RowColData rc,EntMethod &method,int lower_rank,int upper_rank,int verb
) {
  PetscFunctionBegin;
  SET_BASIC_METHOD(method,&*problem_ptr);
  typedef NumeredDofEntity_multiIndex::index<Composite_Name_And_Part_mi_tag>::type NumeredDofsByNameAndPart;
  NumeredDofsByNameAndPart *dofs;
  switch (rc) {
    case ROW:
      dofs = &problem_ptr->numered_dofs_rows->get<Composite_Name_And_Part_mi_tag>();
      break;
    case COL:
      dofs = &problem_ptr->numered_dofs_cols->get<Composite_Name_And_Part_mi_tag>();
      break;
    default:
     SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,"not implemented");
  }
  NumeredDofsByNameAndPart::iterator miit = dofs->lower_bound(boost::make_tuple(field_name,lower_rank));
  NumeredDofsByNameAndPart::iterator hi_miit = dofs->upper_bound(boost::make_tuple(field_name,upper_rank));
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(;miit!=hi_miit;miit++) {
    method.dofPtr = &(*(*miit)->getDofEntityPtr());
    method.dofNumeredPtr = &*(*miit);
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const std::string &problem_name,
  const std::string &field_name,
  RowColData rc,                     // ROW or COL
  EntMethod &method,                 // Finite element instance proceesd on each DOF
  int lower_rank,                    // Only DOFs on processor higher or equal to this are processed
  int upper_rank,                    // Only DOFs lowest or higher to this are processed
  int verb                           // verbosity level
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  typedef MoFEMProblem_multiIndex::index<Problem_mi_tag>::type mofem_problems_by_name;
  // find p_miit
  mofem_problems_by_name &pRoblems_set = pRoblems.get<Problem_mi_tag>();
  mofem_problems_by_name::iterator p_miit = pRoblems_set.find(problem_name);
  if(p_miit == pRoblems_set.end()) SETERRQ1(comm,1,"problem not in database %s",problem_name.c_str());
  ierr = loop_dofs(&*p_miit,field_name,rc,method,lower_rank,upper_rank,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(
  const std::string &problem_name,
  const std::string &field_name,
  RowColData rc,
  EntMethod &method,
  int verb
) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  ierr = loop_dofs(problem_name,field_name,rc,method,0,sIze,verb); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode Core::loop_dofs(const std::string &field_name,EntMethod &method,int verb) {
  PetscFunctionBegin;
  if(verb==-1) verb = verbose;
  SET_BASIC_METHOD(method,NULL);
  DofEntityByFieldName::iterator miit,hi_miit;
  miit = dofsField.get<FieldName_mi_tag>().lower_bound(field_name);
  hi_miit = dofsField.get<FieldName_mi_tag>().upper_bound(field_name);
  method.loopSize = distance(miit,hi_miit);
  ierr = method.preProcess(); CHKERRQ(ierr);
  for(int nn = 0; miit!=hi_miit; miit++, nn++) {
    method.nInTheLoop = nn;
    method.dofPtr = &*(*miit);
    method.dofNumeredPtr = NULL;
    try {
      ierr = method(); CHKERRQ(ierr);
    } catch (MoFEMException const &e) {
      SETERRQ(comm,e.errorCode,e.errorMessage);
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << "throw in method: " << ex.what() << std::endl;
      SETERRQ(comm,MOFEM_DATA_INCONSISTENCY,ss.str().c_str());
    }
  }
  ierr = method.postProcess(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

}
