/* This file is part of MoFEM.
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

#include <MoFEM.hpp>
using namespace MoFEM;

#include <MethodForForceScaling.hpp>
#include <DirichletBC.hpp>
#include <PostProcOnRefMesh.hpp>
#include <ThermalElement.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <AnalyticalDirichlet.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

static int debug = 1;
static char help[] = "...\n\n";

struct AnaliticalFunction {

  std::vector<ublas::vector<double> > val;

  std::vector<ublas::vector<double> >& operator()(double x,double y,double z) {
    val.resize(1);
    val[0].resize(1);
    (val[0])[0] = pow(x,1);
    return val;
  }
};

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm)
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  BARRIER_RANK_END(pcomm)

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("TEMP",H1,1); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.add_problem("BC_PROBLEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_add_bit("BC_PROBLEM",bit_level0); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 2;
  }
  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ThermalElement thermal_elements(m_field);
  ierr = thermal_elements.addThermalElements("TEMP"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","THERMAL_FE"); CHKERRQ(ierr);

  Range bc_tris;
  for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(m_field,"ANALYTICAL_BC",it)) {
    rval = moab.get_entities_by_type(it->getMeshSet(),MBTRI,bc_tris,true); CHKERRQ_MOAB(rval);
  }

  AnalyticalDirichletBC analytical_bc(m_field);
  ierr = analytical_bc.initializeProblem(m_field,"BC_FE","TEMP",bc_tris); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("BC_PROBLEM","BC_FE"); CHKERRQ(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  ierr = m_field.partition_simple_problem("BC_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("BC_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("BC_PROBLEM"); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.VecCreateGhost("TEST_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec T;
  ierr = VecDuplicate(F,&T); CHKERRQ(ierr);
  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays("TEST_PROBLEM",&A); CHKERRQ(ierr);

  ierr = thermal_elements.setThermalFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFiniteElementLhsOperators("TEMP",A); CHKERRQ(ierr);
  ierr = thermal_elements.setThermalFluxFiniteElementRhsOperators("TEMP",F); CHKERRQ(ierr);

  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  //analytical Dirichlet bc
  AnalyticalDirichletBC::DirichletBC analytical_ditihlet_bc(m_field,"TEMP",A,T,F);

  //solve for ditihlet bc dofs
  ierr = analytical_bc.setProblem(m_field,"BC_PROBLEM"); CHKERRQ(ierr);

  boost::shared_ptr<AnaliticalFunction> testing_function = boost::shared_ptr<AnaliticalFunction>(new AnaliticalFunction);

  ierr = analytical_bc.setApproxOps(m_field,"TEMP",bc_tris,testing_function,0); CHKERRQ(ierr);
  ierr = analytical_bc.solveProblem(m_field,"BC_PROBLEM","BC_FE",analytical_ditihlet_bc,bc_tris); CHKERRQ(ierr);

  ierr = analytical_bc.destroyProblem(); CHKERRQ(ierr);

  //preproc
  ierr = m_field.problem_basic_method_preProcess("TEST_PROBLEM",analytical_ditihlet_bc); CHKERRQ(ierr);
  //ierr = m_field.set_global_ghost_vector("TEST_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = m_field.loop_finite_elements("TEST_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeRhs()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","THERMAL_FE",thermal_elements.getLoopFeLhs()); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","THERMAL_FLUX_FE",thermal_elements.getLoopFeFlux()); CHKERRQ(ierr);

  //postproc
  ierr = m_field.problem_basic_method_postProcess("TEST_PROBLEM",analytical_ditihlet_bc); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = VecScale(F,-1); CHKERRQ(ierr);
  //std::string wait;
  //std::cout << "\n matrix is coming = \n" << std::endl;
  //ierr = MatView(A,PETSC_VIEWER_DRAW_WORLD);
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,T); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  //ierr = m_field.problem_basic_method_preProcess("TEST_PROBLEM",analytical_ditihlet_bc); CHKERRQ(ierr);
  ierr = m_field.set_global_ghost_vector("TEST_PROBLEM",ROW,T,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = m_field.set_local_ghost_vector("TEST_PROBLEM",ROW,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //ierr = VecView(T,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PetscReal pointwisenorm;
  ierr = VecMax(T,NULL,&pointwisenorm);
  std::cout << "\n The Global Pointwise Norm of error for this problem is : " << pointwisenorm << std::endl;

  // PetscViewer viewer;
  // PetscViewerASCIIOpen(PETSC_COMM_WORLD,"thermal_with_analytical_bc.txt",&viewer);
  // ierr = VecChop(T,1e-4); CHKERRQ(ierr);
  // ierr = VecView(T,viewer); CHKERRQ(ierr);
  // ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


  double sum = 0;
  ierr = VecSum(T,&sum); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"sum  = %9.8e\n",sum); CHKERRQ(ierr);
  double fnorm;
  ierr = VecNorm(T,NORM_2,&fnorm); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"fnorm  = %9.8e\n",fnorm); CHKERRQ(ierr);
  if(fabs(sum+6.46079983e-01)>1e-7) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
  }
  if(fabs(fnorm-4.26080052e+00)>1e-6) {
    SETERRQ(PETSC_COMM_WORLD,MOFEM_ATOM_TEST_INVALID,"Failed to pass test");
  }


  if(debug) {

    PostProcVolumeOnRefinedMesh post_proc(m_field);
    ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("TEMP"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = post_proc.addFieldValuesGradientPostProc("TEMP"); CHKERRQ(ierr);
    ierr = m_field.loop_finite_elements("TEST_PROBLEM","THERMAL_FE",post_proc); CHKERRQ(ierr);
    ierr = post_proc.writeFile("out.h5m"); CHKERRQ(ierr);

  }

  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
