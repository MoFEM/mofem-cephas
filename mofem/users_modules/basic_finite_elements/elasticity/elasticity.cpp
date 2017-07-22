/** \file elasticity.cpp
 * \ingroup nonlinear_elastic_elem
 * \example elasticity.cpp

 The example shows how to solve the linear elastic problem. An example can read
 file with temperature field, then thermal stresses are included.

 What example can do:
 - take into account temperature field, i.e. calculate thermal stresses and deformation
 - stationary and time depend field is considered
 - take into account gravitational body forces
 - take in account fluid pressure
 - can work with higher order geometry definition
 - works on distributed meshes
 - multi-grid solver where each grid level is approximation order level
 - each mesh block can have different material parameters and approximation order

See example how code can be used \cite jordi:2017,
 \image html SquelaDamExampleByJordi.png "Example what you can do with this code. Analysis of the arch dam of Susqueda, located in Catalonia (Spain)" width=800px

 This is an example of application code; it does not show how elements are implemented. Example presents how to:
 - read mesh
 - set-up problem
 - run finite elements on the problem
 - assemble matrices and vectors
 - solve the linear system of equations
 - save results


 If you like to see how to implement finite elements, material, are other parts of the code, look here;
 - Hooke material, see \ref Hooke
 - Thermal-stress assembly, see \ref  ThermalElement
 - Body forces element, see \ref BodyFroceConstantField
 - Fluid pressure element, see \ref FluidPressure
 - The general implementation of an element for arbitrary Lagrangian-Eulerian
 aformulation for a nonlinear elastic problem is here \ref
 NonlinearElasticElement. Here we limit ourselves to Hooke equation and fix
 mesh, so the problem becomes linear. Not that elastic element is implemented
 with automatic differentiation.

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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include <Hooke.hpp>

using namespace boost::numeric;




static char help[] =
  "-my_block_config set block data\n"
  "\n";

//const double young_modulus = 1;
//const double poisson_ratio = 0.0;

struct BlockOptionData {
  int oRder;
  double yOung;
  double pOisson;
  double initTemp;
  BlockOptionData():
    oRder(-1),
    yOung(-1),
    pOisson(-2),
    initTemp(0) {}
};

int main(int argc, char *argv[]) {

  // Initialize PETSCc
  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

    PetscBool flg_block_config,flg_file;
    char mesh_file_name[255];
    char block_config_file[255];
    PetscInt order = 2;
    PetscBool is_partitioned = PETSC_FALSE;

    // Read options from command line
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","Elastic Config","none"); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-my_file",
      "mesh file name","",
      "mesh.h5m",mesh_file_name,255,&flg_file
    ); CHKERRQ(ierr);
    ierr = PetscOptionsInt(
      "-my_order",
      "default approximation order","",
      2,&order,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsBool(
      "-my_is_partitioned",
      "set if mesh is partitioned (this result that each process keeps only part of the mes","",
      PETSC_FALSE,&is_partitioned,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsString("-my_block_config",
    "elastic configure file name","",
    "block_conf.in",block_config_file,255,&flg_block_config
  ); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Throw error if file with mesh is not give
  if(flg_file != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // Create mesh databse
  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;

  // Create moab communicator
  // Create separate MOAB communicator, it is duplicate of PETSc communicator.
  // NOTE That this should eliminate potential communication problems between
  // MOAB and PETSC functions.
  MPI_Comm moab_comm_world;
  MPI_Comm_dup(PETSC_COMM_WORLD,&moab_comm_world);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,moab_comm_world);

  // Read whole mesh or part of is if partitioned
  if(is_partitioned == PETSC_TRUE) {
    // This is a case of distributed mesh and algebra. In that case each processor
    // keep only part of the problem.
    const char *option;
    option = "PARALLEL=READ_PART;"
    "PARALLEL_RESOLVE_SHARED_ENTS;"
    "PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  } else {
    // If that case we have distributed algebra, i.e. assembly of vectors and
    // matrices is in parallel, but whole mesh is stored on all processors.
    // Solver and matrix scales well, however problem set-up of problem is
    // not fully parallel.
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  }

  // Create MoFEM databse and link it to MoAB
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  // Print boundary conditions and material parameters
  MeshsetsManager *mmanager_ptr;
  ierr = m_field.query_interface(mmanager_ptr); CHKERRQ(ierr);
  ierr = mmanager_ptr->printDisplacementSet(); CHKERRQ(ierr);
  ierr = mmanager_ptr->printForceSet(); CHKERRQ(ierr);
  ierr = mmanager_ptr->printMaterialsSet(); CHKERRQ(ierr);

  // Set bit refinement level to all entities (we do not refine mesh in this example
  // so all entities have the same bit refinement level)
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  // Declare approximation fields
  ierr = m_field.add_field("DISPLACEMENT",H1,AINSWORTH_LOBBATO_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  // We can use higher oder geometry to define body
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);

  // Declare problem

  // Add entities (by tets) to the field ( all entities in the mesh, root_set = 0 )
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  // Set apportion order.
  // See Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes.
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  // Set order of approximation of geometry.
  // Apply 2nd order only on skin (or in in whole body)
  {
    Skinner skin(&m_field.get_moab());
    Range faces,tets;
    rval = m_field.get_moab().get_entities_by_type(0,MBTET,tets); CHKERRQ_MOAB(rval);
    // rval = skin.find_skin(0,tets,false,faces); CHKERRQ_MOAB(rval);
    Range edges;
    rval = m_field.get_moab().get_adjacencies(
      tets,1,false,edges,moab::Interface::UNION
    ); CHKERRQ_MOAB(rval);
    // rval = m_field.get_moab().get_adjacencies(
    //   faces,1,false,edges,moab::Interface::UNION
    // ); CHKERRQ_MOAB(rval);
    // ierr = m_field.synchronise_entities(edges); CHKERRQ(ierr);
    ierr = m_field.set_field_order(edges,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  }
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  // Configure blocks by parsing config file. It allows setting approximation
  // order for each block independently.
  std::map<int,BlockOptionData> block_data;
  if(flg_block_config) {
    try {
      ifstream ini_file(block_config_file);
      //std::cerr << block_config_file << std::endl;
      po::variables_map vm;
      po::options_description config_file_options;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        std::ostringstream str_order;
        str_order << "block_" << it->getMeshsetId() << ".displacement_order";
        config_file_options.add_options()
        (str_order.str().c_str(),po::value<int>(&block_data[it->getMeshsetId()].oRder)->default_value(order));
        std::ostringstream str_cond;
        str_cond << "block_" << it->getMeshsetId() << ".young_modulus";
        config_file_options.add_options()
        (str_cond.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].yOung)->default_value(-1));
        std::ostringstream str_capa;
        str_capa << "block_" << it->getMeshsetId() << ".poisson_ratio";
        config_file_options.add_options()
        (str_capa.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].pOisson)->default_value(-2));
        std::ostringstream str_init_temp;
        str_init_temp << "block_" << it->getMeshsetId() << ".initial_temperature";
        config_file_options.add_options()
        (str_init_temp.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].initTemp)->default_value(0));
      }
      po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
      store(parsed,vm);
      po::notify(vm);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        if(block_data[it->getMeshsetId()].oRder == -1) continue;
        if(block_data[it->getMeshsetId()].oRder == order) continue;
        PetscPrintf(PETSC_COMM_WORLD,"Set block %d order to %d\n",it->getMeshsetId(),block_data[it->getMeshsetId()].oRder);
        Range block_ents;
        rval = moab.get_entities_by_handle(it->getMeshset(),block_ents,true); CHKERRQ_MOAB(rval);
        Range ents_to_set_order;
        rval = moab.get_adjacencies(block_ents,3,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ_MOAB(rval);
        ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
        rval = moab.get_adjacencies(block_ents,2,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ_MOAB(rval);
        rval = moab.get_adjacencies(block_ents,1,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ_MOAB(rval);
        ierr = m_field.synchronise_entities(ents_to_set_order); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents_to_set_order,"DISPLACEMENT",block_data[it->getMeshsetId()].oRder); CHKERRQ(ierr);
      }
      std::vector<std::string> additional_parameters;
      additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
      for(std::vector<std::string>::iterator vit = additional_parameters.begin();
      vit!=additional_parameters.end();vit++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARNING Unrecognized option %s\n",vit->c_str()); CHKERRQ(ierr);
      }
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
  }

  // Add elastic element
  boost::shared_ptr<Hooke<adouble> > hooke_adouble_ptr(new Hooke<adouble>());
  boost::shared_ptr<Hooke<double> > hooke_double_ptr(new Hooke<double>());
  NonlinearElasticElement elastic(m_field,2);
  ierr = elastic.setBlocks(hooke_double_ptr,hooke_adouble_ptr); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = elastic.setOperators("DISPLACEMENT","MESH_NODE_POSITIONS",false,true); CHKERRQ(ierr);

  // Update material parameters. Set material parameters block by block.
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    int id = it->getMeshsetId();
    if(block_data[id].yOung>0) {
      elastic.setOfBlocks[id].E = block_data[id].yOung;
      ierr = PetscPrintf(
        PETSC_COMM_WORLD,"Block %d set Young modulus %3.4g\n",id,elastic.setOfBlocks[id].E
      ); CHKERRQ(ierr);
    }
    if(block_data[id].pOisson>=-1) {
      elastic.setOfBlocks[id].PoissonRatio = block_data[id].pOisson;
      ierr = PetscPrintf(
        PETSC_COMM_WORLD,"Block %d set Poisson ratio %3.4g\n",id,elastic.setOfBlocks[id].PoissonRatio
      ); CHKERRQ(ierr);
    }
  }

  // Add body force element. This is only declaration of element. not its implementation.
  ierr = m_field.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    Range tets;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_type(tets,MBTET,"BODY_FORCE"); CHKERRQ(ierr);
  }

  // Add Neumann forces, i.e. pressure or traction forces applied on body surface. This
  // is only declaration not implementation.
  ierr = MetaNeummanForces::addNeumannBCElements(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addElement(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaEdgeForces::addElement(m_field,"DISPLACEMENT"); CHKERRQ(ierr);

  // Add fluid pressure finite elements. This is special pressure on the surface from
  // fluid, i.e. preseeure which linearly change with the depth.
  FluidPressure fluid_pressure_fe(m_field);
  // This function only declare element. Element is implemented by operators in
  // class FluidPressure.
  fluid_pressure_fe.addNeumannFluidPressureBCElements("DISPLACEMENT");

  // Add elements for thermo elasticity if temperature field is defined.
  ThermalStressElement thermal_stress_elem(m_field);
  // Check if TEMP field exist, and then add element.
  if(!m_field.check_field("TEMP")) {
    bool add_temp_field = false;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      if(block_data[it->getMeshsetId()].initTemp!=0) {
        add_temp_field = true;
        break;
      }
    }
    if(add_temp_field) {
      ierr = m_field.add_field("TEMP",H1,AINSWORTH_LEGENDRE_BASE,1,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_type(0,MBTET,"TEMP"); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"TEMP",1); CHKERRQ(ierr);
    }
  }
  if(m_field.check_field("TEMP")) {
    ierr = thermal_stress_elem.addThermalSterssElement("ELASTIC","DISPLACEMENT","TEMP"); CHKERRQ(ierr);
  }

  // All is declared, at that point build fields first,
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  // If 10-node test are on the mesh, use mid nodes to set HO-geometry. Class Projection10NodeCoordsOnField
  // do the trick.
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  if(m_field.check_field("TEMP")) {
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      if(block_data[it->getMeshsetId()].initTemp!=0) {
        PetscPrintf(PETSC_COMM_WORLD,"Set block %d temperature to %3.2g\n",
        it->getMeshsetId(),block_data[it->getMeshsetId()].initTemp);
        Range block_ents;
        rval = moab.get_entities_by_handle(it->meshset,block_ents,true); CHKERRQ_MOAB(rval);
        Range vertices;
        rval = moab.get_connectivity(block_ents,vertices,true); CHKERRQ_MOAB(rval);
        ierr = m_field.set_field(block_data[it->getMeshsetId()].initTemp,MBVERTEX,vertices,"TEMP"); CHKERRQ(ierr);
      }
    }
  }

  // Build database for elements. Actual implementation of element is not need here,
  // only elements has to be declared.
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  // Build adjacencies between elements and field entities
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  // Register MOFEM DM implementation in PETSc
  ierr = DMRegister_MGViaApproxOrders("MOFEM"); CHKERRQ(ierr);

  // Create DM manager
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,"MOFEM");CHKERRQ(ierr);
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,"ELASTIC_PROB",bit_level0); CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMMoFEMSetIsPartitioned(dm,is_partitioned); CHKERRQ(ierr);
  // Add elements to DM manager
  ierr = DMMoFEMAddElement(dm,"ELASTIC"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"BODY_FORCE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"FLUID_PRESSURE_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"FORCE_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"PRESSURE_FE"); CHKERRQ(ierr);
  ierr = DMSetUp(dm); CHKERRQ(ierr);

  // Create matrices & vectors. Note that native PETSc DM interface is used, but
  // ubder the PETSc interface MOFEM implementation is running.
  Vec F,D,D0;
  ierr = DMCreateGlobalVector(dm,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D0); CHKERRQ(ierr);
  Mat Aij;
  ierr = DMCreateMatrix(dm,&Aij); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

  // Zero vectors and matrices
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecZeroEntries(D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  // This controls how kinematic constrains are set, by blockset or nodeset. Cubit
  // sets kinetic booundary conditions by blockset.
  bool flag_cubit_disp = false;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|DISPLACEMENTSET,it)) {
    flag_cubit_disp = true;
  }

  // Below particular implementations of finite elements are used to assemble
  // problem matrixes and vectors.  Implementation of element does not change
  // how element is declared.

  // Assemble Aij and F. Define dirihlet_bc element, which stets constrains to MatrixDouble
  // and the right hand side vector.
  boost::shared_ptr<FEMethod> dirihlet_bc_ptr;

  // if normally defined boundary conditions are not found, try to use DISPLACEMENT blockset.
  // To implementations available here, depending how BC is defined on mesh file.
  if(!flag_cubit_disp){
    dirihlet_bc_ptr = boost::shared_ptr<FEMethod>(new
      DirichletBCFromBlockSetFEMethodPreAndPostProcWithFlags(
        m_field,"DISPLACEMENT","DISPLACEMENT",Aij,D0,F
      )
    );
  } else {
    dirihlet_bc_ptr = boost::shared_ptr<FEMethod>(
      new DisplacementBCFEMethodPreAndPostProc(m_field,"DISPLACEMENT",Aij,D0,F)
    );
  }
  // That sets dirihlet_bc objects that problem is linear, i.e. no newton (SNES) solver
  // is run for this problem.
  dirihlet_bc_ptr->snes_ctx = FEMethod::CTX_SNESNONE;
  dirihlet_bc_ptr->ts_ctx = FEMethod::CTX_TSNONE;

  // D0 vector will store initial displacements
  ierr = VecZeroEntries(D0); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = DMoFEMMeshToLocalVector(dm,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  // Run dirihlet_bc, from that on the mesh set values in vector D0. Run implementation
  // of particular dirihlet_bc.
  ierr = DMoFEMPreProcessFiniteElements(dm,dirihlet_bc_ptr.get()); CHKERRQ(ierr);
  // Set values from D0 on the field (on the mesh)
  ierr = DMoFEMMeshToLocalVector(dm,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  // Calculate residual forces as result of applied kinematic constrains. Run implementation
  // of particular finite element implementation. Look how NonlinearElasticElement is implemented,
  // in that case. We run NonlinearElasticElement with hook material.
  // Calculate right hand side vector
  elastic.getLoopFeRhs().snes_f = F;
  ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&elastic.getLoopFeRhs()); CHKERRQ(ierr);
  // Assemble matrix
  elastic.getLoopFeLhs().snes_B = Aij;
  ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&elastic.getLoopFeLhs()); CHKERRQ(ierr);

  // Assemble pressure and traction forces. Run particular implemented for do this, see
  // MetaNeummanForces how this is implemented.
  boost::ptr_map<std::string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setMomentumFluxOperators(m_field,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  {
    boost::ptr_map<std::string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
    for(;mit!=neumann_forces.end();mit++) {
      ierr = DMoFEMLoopFiniteElements(dm,mit->first.c_str(),&mit->second->getLoopFe()); CHKERRQ(ierr);
    }
  }
  // Assemble forces applied to nodes, see implementation in NodalForce
  boost::ptr_map<std::string,NodalForce> nodal_forces;
  ierr = MetaNodalForces::setOperators(m_field,nodal_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  {
    boost::ptr_map<std::string,NodalForce>::iterator fit = nodal_forces.begin();
    for(;fit!=nodal_forces.end();fit++) {
      ierr = DMoFEMLoopFiniteElements(dm,fit->first.c_str(),&fit->second->getLoopFe()); CHKERRQ(ierr);
    }
  }
  // Assemble edge forces
  boost::ptr_map<std::string,EdgeForce> edge_forces;
  ierr = MetaEdgeForces::setOperators(m_field,edge_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  {
    boost::ptr_map<std::string,EdgeForce>::iterator fit = edge_forces.begin();
    for(;fit!=edge_forces.end();fit++) {
      ierr = DMoFEMLoopFiniteElements(dm,fit->first.c_str(),&fit->second->getLoopFe()); CHKERRQ(ierr);
    }
  }
  // Assemble body forces, implemented in BodyFroceConstantField
  BodyFroceConstantField body_forces_methods(m_field);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->getMeshsetId()); CHKERRQ(ierr);
  }
  ierr = DMoFEMLoopFiniteElements(dm,"BODY_FORCE",&body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  // Assemble fluid pressure forces
  ierr = fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators("DISPLACEMENT",F,false,true); CHKERRQ(ierr);
  ierr = DMoFEMLoopFiniteElements(dm,"FLUID_PRESSURE_FE",&fluid_pressure_fe.getLoopFe()); CHKERRQ(ierr);

  // Apply kinematic constrain to right hand side vector and matrix
  ierr = DMoFEMPostProcessFiniteElements(dm,dirihlet_bc_ptr.get()); CHKERRQ(ierr);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  // Set matrix positive defined and symmetric for Cholesky and icc pre-conditioner
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = VecScale(F,-1); CHKERRQ(ierr);

  // Create solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetDM(solver,dm); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  // Setup multi-grid pre-conditioner if set from command line
  {
    //from PETSc example ex42.c
    PetscBool same = PETSC_FALSE;
    PC pc;
    ierr = KSPGetPC(solver,&pc); CHKERRQ(ierr);
    PetscObjectTypeCompare((PetscObject)pc,PCMG,&same);
    if (same) {
      PCMGSetUpViaApproxOrdersCtx pc_ctx(dm,Aij,true);
      ierr = PCMGSetUpViaApproxOrders(pc,&pc_ctx); CHKERRQ(ierr);
      ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
    } else {
      // Operators are already set, do not use DM for doing that
      ierr = KSPSetDMActive(solver,PETSC_FALSE); CHKERRQ(ierr);
    }
  }
  ierr = KSPSetInitialGuessKnoll(solver,PETSC_FALSE); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(solver,PETSC_TRUE); CHKERRQ(ierr);
  // Set up solver
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  // Set up post-procesor. It is some generic implementation of finite element.
  PostProcVolumeOnRefinedMesh post_proc(m_field);
  // Add operators to the elements, strating with some generic
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  // Add problem specific operator on element to post-process stresses
  post_proc.getOpPtrVector().push_back(
    new PostPorcHookStress(
      m_field,
      post_proc.postProcMesh,
      post_proc.mapGaussPts,
      "DISPLACEMENT",
      post_proc.commonData,
      &elastic.setOfBlocks
    )
  );

  // Temperature field is defined on the mesh
  if(m_field.check_field("TEMP")) {

    // Create thermal vector
    Vec F_thermal;
    ierr = VecDuplicate(F,&F_thermal); CHKERRQ(ierr);

    // Set up implementation for calculation of thermal stress vector. Look how
    // thermal stresses and vector is assembled in ThermalStressElement.
    ierr = thermal_stress_elem.setThermalStressRhsOperators("DISPLACEMENT","TEMP",F_thermal); CHKERRQ(ierr);

    SeriesRecorder *recorder_ptr;
    ierr = m_field.query_interface(recorder_ptr); CHKERRQ(ierr);
    // Read time series and do thermo-elastic analysis, this is when time dependent
    // temperature problem was run before on the mesh. It means that before non-stationary
    // problem was solved for temperature and filed "TEMP" is stored for subsequent time
    // steps in the recorder.
    if( recorder_ptr->check_series("THEMP_SERIES") ) {
      // This is time dependent case, so loop of data series stored by tape recorder.
      // Loop over time steps
      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {
        PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
        // Load field data for this time step
        ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
        ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        // Calculate the right-hand side vector as result of thermal stresses. It MetaNodalForces
        // that on "ELASTIC" element data structure the element implementation in thermal_stress_elem
        // is executed.
        ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
        // Assemble vector
        ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
        // Accumulate ghost dofs
        ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

        // Calculate norm of vector and print values
        PetscReal nrm_F;
        ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
        PetscReal nrm_F_thremal;
        ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);

        ierr = VecScale(F_thermal,-1); CHKERRQ(ierr); //check this !!!
        ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);

        // Set dirihlet boundary to thermal stresses vector
        dirihlet_bc_ptr->snes_x = D;
        dirihlet_bc_ptr->snes_f = F_thermal;
        ierr = DMoFEMPostProcessFiniteElements(dm,dirihlet_bc_ptr.get()); CHKERRQ(ierr);

        // Solve problem
        ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
        // Add boundary conditions vector
        ierr = VecAXPY(D,1.,D0); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

        // Sava data on the mesh
        ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

        // Save data on mesh
        ierr = DMoFEMPreProcessFiniteElements(dm,dirihlet_bc_ptr.get()); CHKERRQ(ierr);
        // Post-process results
        ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
        std::ostringstream o1;
        o1 << "out_" << sit->step_number << ".h5m";
        ierr = post_proc.writeFile(o1.str().c_str()); CHKERRQ(ierr);
      }
    } else {

      // This is a case when stationary problem for temperature was solved.
      ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      // Calculate the right-hand side vector with thermal stresses
      ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
      // Assemble vector
      ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
      // Accumulate ghost dofs
      ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

      // Calculate norm
      PetscReal nrm_F;
      ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
      PetscReal nrm_F_thremal;
      ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);

      // Add thermal stress vector and other forces vector
      ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);
      ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);

      // Apply kinetic boundary conditions
      dirihlet_bc_ptr->snes_x = D;
      dirihlet_bc_ptr->snes_f = F_thermal;
      ierr = DMoFEMPostProcessFiniteElements(dm,dirihlet_bc_ptr.get()); CHKERRQ(ierr);

      // Solve problem
      ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
      ierr = VecAXPY(D,1.,D0); CHKERRQ(ierr);
      // Update ghost values for solution vector
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      // Save data on mesh
      ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
      // Save results to file
      ierr = post_proc.writeFile("out.h5m"); CHKERRQ(ierr);

    }

    // Destroy vector, no needed any more
    ierr = VecDestroy(&F_thermal); CHKERRQ(ierr);

  } else {
    // Elastic analysis no temperature field
    // Solve for vector D
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    // Add kinetic boundary conditions
    ierr = VecAXPY(D,1.,D0); CHKERRQ(ierr);
    // Update ghost values
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    // Save data from vector on mesh
    ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    // Porticoes results by running post_proc implementation of element
    ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
    // Weite mesh in parallel (uisng h5m MOAB format, writing is in parallel)
    ierr = post_proc.writeFile("out.h5m"); CHKERRQ(ierr);
  }

  // Calculate elastic energy
  elastic.getLoopFeEnergy().snes_ctx = SnesMethod::CTX_SNESNONE;
  elastic.getLoopFeEnergy().eNergy = 0;
  ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&elastic.getLoopFeEnergy()); CHKERRQ(ierr);
  // Print elastic energy
  PetscPrintf(PETSC_COMM_WORLD,"Elastic energy %6.4e\n",elastic.getLoopFeEnergy().eNergy);

  // Destroy matrices, vecors, solver and DM
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&D0); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = DMDestroy(&dm); CHKERRQ(ierr);

  MPI_Comm_free(&moab_comm_world);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

  return 0;
}
