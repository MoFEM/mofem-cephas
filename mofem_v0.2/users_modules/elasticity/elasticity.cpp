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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include <DirichletBC.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>

#include <PotsProcOnRefMesh.hpp>
#include <PostProcHookStresses.hpp>

#include <adolc/adolc.h> 
#include <NonLienarElasticElement.hpp>
#include <Hooke.hpp>

#include <PCMGSetUpViaApproxOrders.hpp>

using namespace boost::numeric;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = 
  "-my_block_config set block data\n"
  "\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;

struct BlockOptionData {
  int oRder;
  double yOung;
  double pOisson;
  double initTemp;
  BlockOptionData():
    oRder(-1),
    yOung(-1),
    pOisson(-1),
    initTemp(0) {}
};

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;

  PetscBool flg_block_config,flg_file;
  char mesh_file_name[255];
  char block_config_file[255];
  PetscInt order = 2;
  PetscBool is_partitioned = PETSC_FALSE;

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","Elastic Config","none"); CHKERRQ(ierr);
  ierr = PetscOptionsString("-my_file",
    "mesh file name","",
    "mesh.h5m",mesh_file_name,255,&flg_file); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-my_order",
    "default approximation order","",
      2,&order,PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsBool("-my_is_partitioned", 
    "set if mesh is partitioned (this result that each process keeps only part of the mes","",
    PETSC_TRUE,&is_partitioned,PETSC_NULL); CHKERRQ(ierr);
  ierr = PetscOptionsString("-my_block_config",
    "elastic configure file name","",
    "block_conf.in",block_config_file,255,&flg_block_config); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //Reade parameters from line command
  if(flg_file != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;"
      "PARALLEL_RESOLVE_SHARED_ENTS;"
      "PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  Range meshset_level0;
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"meshset_level0 %d\n",meshset_level0.size());
  PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT); 

  //Define problem

  //Fields
  ierr = m_field.add_field("DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

  //Declare problem

  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT",2); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);

  //ierr = m_field.synchronise_field_entities("DISPLACEMENT"); CHKERRQ(ierr);
  //ierr = m_field.synchronise_field_entities("MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  // configure blocks by parsing config file
  // it allow to set approximation order for each block independently
  map<int,BlockOptionData> block_data;
  if(flg_block_config) {
    try {
      ifstream ini_file(block_config_file);  
      //cerr << block_config_file << endl;
      po::variables_map vm;
      po::options_description config_file_options;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        ostringstream str_order;
        str_order << "block_" << it->get_msId() << ".displacemet_order";
        config_file_options.add_options()
	 (str_order.str().c_str(),po::value<int>(&block_data[it->get_msId()].oRder)->default_value(order));
	ostringstream str_cond;
        str_cond << "block_" << it->get_msId() << ".young_modulus";
        config_file_options.add_options()
	 (str_cond.str().c_str(),po::value<double>(&block_data[it->get_msId()].yOung)->default_value(-1));
	ostringstream str_capa;
        str_capa << "block_" << it->get_msId() << ".poisson_ratio";
        config_file_options.add_options()
	 (str_capa.str().c_str(),po::value<double>(&block_data[it->get_msId()].pOisson)->default_value(-1));
	ostringstream str_init_temp;
        str_init_temp << "block_" << it->get_msId() << ".initial_temperature";
        config_file_options.add_options()
	 (str_init_temp.str().c_str(),po::value<double>(&block_data[it->get_msId()].initTemp)->default_value(0));
      }
      po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
      store(parsed,vm);
      po::notify(vm); 
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
	if(block_data[it->get_msId()].oRder == -1) continue;
        if(block_data[it->get_msId()].oRder == order) continue;
	PetscPrintf(PETSC_COMM_WORLD,"Set block %d order to %d\n",it->get_msId(),block_data[it->get_msId()].oRder);
	Range block_ents;
	rval = moab.get_entities_by_handle(it->get_meshset(),block_ents,true); CHKERR(rval);
	Range ents_to_set_order;
	ierr = moab.get_adjacencies(block_ents,3,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
	ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
	ierr = moab.get_adjacencies(block_ents,2,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
	ierr = moab.get_adjacencies(block_ents,1,false,ents_to_set_order,Interface::UNION); CHKERRQ(ierr);
	ierr = m_field.synchronise_entities(ents_to_set_order,2); CHKERRQ(ierr);
	ierr = m_field.synchronise_entities(ents_to_set_order,2); CHKERRQ(ierr);
	ierr = m_field.synchronise_entities(ents_to_set_order,2); CHKERRQ(ierr);

        ierr = m_field.set_field_order(ents_to_set_order,"DISPLACEMENT",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
      }
      vector<string> additional_parameters;
      additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
      for(vector<string>::iterator vit = additional_parameters.begin();
	vit!=additional_parameters.end();vit++) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARNING Unrecognised option %s\n",vit->c_str()); CHKERRQ(ierr);
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
  }

  //define eleatic element
  Hooke<adouble> hooke_adouble;
  Hooke<double> hooke_double;
  NonlinearElasticElement elastic(m_field,2);
  ierr = elastic.setBlocks(&hooke_double,&hooke_adouble); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = elastic.setOperators("DISPLACEMENT","MESH_NODE_POSITIONS",false,true); CHKERRQ(ierr);

  ierr = m_field.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("BODY_FORCE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    Range tets;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    ierr = m_field.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
  }

  //neuman forces
  ierr = MetaNeummanForces::addNeumannBCElements(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = MetaNodalForces::addNodalForceElement(m_field,"DISPLACEMENT"); CHKERRQ(ierr);

  //define fluid pressure finite elements
  FluidPressure fluid_pressure_fe(m_field);
  fluid_pressure_fe.addNeumannFluidPressureBCElements("DISPLACEMENT");
  //define elements for thermo elasticity if themperature field avelible
  ThermalStressElement thermal_stress_elem(m_field);

  if(!m_field.check_field("TEMP")) {
    bool add_temp_field = false;
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      if(block_data[it->get_msId()].initTemp!=0) {
	add_temp_field = true;
	break;
      }
    }
    if(add_temp_field) {
      ierr = m_field.add_field("TEMP",H1,1,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(0,"TEMP"); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"TEMP",1); CHKERRQ(ierr);
    }
  }

  if(m_field.check_field("TEMP")) {
    ierr = thermal_stress_elem.addThermalSterssElement("ELASTIC","DISPLACEMENT","TEMP"); CHKERRQ(ierr);
  } 

  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //get HO gemetry for 10 node tets
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  if(m_field.check_field("TEMP")) {
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
      if(block_data[it->get_msId()].initTemp!=0) {
	PetscPrintf(PETSC_COMM_WORLD,"Set block %d temperature to %3.2g\n",
	  it->get_msId(),block_data[it->get_msId()].initTemp);
	Range block_ents;
	rval = moab.get_entities_by_handle(it->meshset,block_ents,true); CHKERR(rval);
	Range vertices;
	rval = moab.get_connectivity(block_ents,vertices,true); CHKERR_PETSC(rval);
	ierr = m_field.set_field(block_data[it->get_msId()].initTemp,MBVERTEX,vertices,"TEMP"); CHKERRQ(ierr);
      }
    }
  }

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //print bcs
  ierr = m_field.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = m_field.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = m_field.print_cubit_materials_set(); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","BODY_FORCE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","FLUID_PRESSURE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_PROB","PRESSURE_FE"); CHKERRQ(ierr);

  DMType dm_name = "ELASTIC_PROB";
  ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);
  //craete dm instance
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,dm_name);CHKERRQ(ierr);
  //set dm datastruture whict created mofem datastructures
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  ierr = DMMoFEMSetIsPartitioned(dm,is_partitioned); CHKERRQ(ierr);
  //add elements to dm
  ierr = DMMoFEMAddElement(dm,"ELASTIC"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"BODY_FORCE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"FLUID_PRESSURE_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"FORCE_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"PRESSURE_FE"); CHKERRQ(ierr);

  ierr = DMSetUp(dm); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = DMCreateGlobalVector_MoFEM(dm,&F); CHKERRQ(ierr);
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = DMCreateMatrix_MoFEM(dm,&Aij); CHKERRQ(ierr);
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  //assemble Aij and F
 
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"DISPLACEMENT",Aij,D,F);
  //set kinematic boundary conditions
  ierr = DMoFEMPreProcessFiniteElements(dm,&my_dirichlet_bc); CHKERRQ(ierr);

  //elastic element matrix
  elastic.getLoopFeLhs().snes_B = Aij;
  ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&elastic.getLoopFeLhs()); CHKERRQ(ierr);

  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = DMoFEMLoopFiniteElements(dm,mit->first.c_str(),&mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //noadl forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  ierr = MetaNodalForces::setNodalForceElementOperators(m_field,nodal_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    ierr = DMoFEMLoopFiniteElements(dm,fit->first.c_str(),&fit->second->getLoopFe()); CHKERRQ(ierr);
  }
  //body forces
  BodyFroceConstantField body_forces_methods(m_field);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|BODYFORCESSET,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = DMoFEMLoopFiniteElements(dm,"BODY_FORCE",&body_forces_methods.getLoopFe()); CHKERRQ(ierr);
  //fluid pressure
  ierr = fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators("DISPLACEMENT",F,false,true); CHKERRQ(ierr);
  ierr = DMoFEMLoopFiniteElements(dm,"FLUID_PRESSURE_FE",&fluid_pressure_fe.getLoopFe()); CHKERRQ(ierr);

  //postproc
  ierr = DMoFEMPostProcessFiniteElements(dm,&my_dirichlet_bc); CHKERRQ(ierr);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(Aij,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = VecScale(F,-1); CHKERRQ(ierr);

  //PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  {
    //from PETSc example ex42.c
    PetscBool same = PETSC_FALSE;
    PC pc;
    ierr = KSPGetPC(solver,&pc); CHKERRQ(ierr);
    PetscObjectTypeCompare((PetscObject)pc,PCMG,&same);
    if (same) {
      ierr = PCMGSetUpViaApproxOrders(pc,&m_field,"ELASTIC_PROB"); CHKERRQ(ierr);
    }
  }
  ierr = KSPSetUp(solver); CHKERRQ(ierr);


  PostPocOnRefinedMesh post_proc(m_field);
  ierr = post_proc.generateRefereneElemenMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISPLACEMENT"); CHKERRQ(ierr);
  //add postpocessing for sresses
  post_proc.get_op_to_do_Rhs().push_back(
	  new PostPorcStress(
	    m_field,
	    post_proc.postProcMesh,
	    post_proc.mapGaussPts,
	    "DISPLACEMENT",
	    post_proc.commonData));

  if(m_field.check_field("TEMP")) {
    //read time series and do thermo elastci analysis
    Vec F_thermal;
    ierr = VecDuplicate(F,&F_thermal); CHKERRQ(ierr);
    ierr = thermal_stress_elem.setThermalStressRhsOperators("DISPLACEMENT","TEMP",F_thermal); CHKERRQ(ierr);

    SeriesRecorder *recorder_ptr;
    ierr = m_field.query_interface(recorder_ptr); CHKERRQ(ierr);
    if( recorder_ptr->check_series("THEMP_SERIES") ) {
      for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {
	PetscPrintf(PETSC_COMM_WORLD,"Process step %d\n",sit->get_step_number());
	ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
	ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);

	PetscReal nrm_F;
	ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
	PetscReal nrm_F_thremal;
	ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);
	ierr = VecScale(F_thermal,-1); CHKERRQ(ierr); //check this !!!
	ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
        my_dirichlet_bc.snes_x = D;
        my_dirichlet_bc.snes_f = F_thermal;
	ierr = DMoFEMPostProcessFiniteElements(dm,&my_dirichlet_bc); CHKERRQ(ierr);
	ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

	//Save data on mesh
	ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
	ostringstream o1;
	o1 << "out_" << sit->step_number << ".h5m";
	rval = post_proc.postProcMesh.write_file(o1.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
      }
    } else {
      ierr = VecZeroEntries(F_thermal); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&thermal_stress_elem.getLoopThermalStressRhs()); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F_thermal,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F_thermal); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F_thermal); CHKERRQ(ierr);
 
      PetscReal nrm_F;
      ierr = VecNorm(F,NORM_2,&nrm_F); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F = %6.4e\n",nrm_F);
      PetscReal nrm_F_thremal;
      ierr = VecNorm(F_thermal,NORM_2,&nrm_F_thremal); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"norm2 F_thernal = %6.4e\n",nrm_F_thremal);
      ierr = VecScale(F_thermal,-1); CHKERRQ(ierr);  // check this !!!
      ierr = VecAXPY(F_thermal,1,F); CHKERRQ(ierr);
      my_dirichlet_bc.snes_x = D;
      my_dirichlet_bc.snes_f = F_thermal;
      ierr = DMoFEMPostProcessFiniteElements(dm,&my_dirichlet_bc); CHKERRQ(ierr);

      ierr = KSPSolve(solver,F_thermal,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

      //Save data on mesh
      ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
      rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
    }
    ierr = VecDestroy(&F_thermal); CHKERRQ(ierr);
  } else {
    // elastic analys
    ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  
    //Save data on mesh
    ierr = DMoFEMMeshToLocalVector(dm,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&post_proc); CHKERRQ(ierr);
    rval = post_proc.postProcMesh.write_file("out.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
  }

  elastic.getLoopFeEnergy().snes_ctx = SnesMethod::CTX_SNESNONE;
  ierr = DMoFEMLoopFiniteElements(dm,"ELASTIC",&elastic.getLoopFeEnergy()); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Elastic energy %6.4e\n",elastic.getLoopFeEnergy().eNergy);
      
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PetscFinalize();

}

