/** \file thermal_unsteady.cpp
 \ingroup mofem_thermal_elem
 \brief Example of thermal unsteady analyze.
*/

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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#ifdef __GROUNDSURFACETEMERATURE_HPP

  #include <GenricClimateModel.hpp>
  #include <GroundSurfaceTemerature.hpp>

  #include <time.h>
  extern "C" {
    #include <spa.h>
  }
  #include <CrudeClimateModel.hpp>

#endif // __GROUNDSURFACETEMERATURE_HPP

static char help[] =
  "-my_file mesh file\n"
  "-order set approx. order to all blocks\n"
  "-my_block_config set block data\n"
  "-my_ground_analysis_data data for crude climate model\n"
  "\n";

struct BlockOptionData {
  int oRder;
  double cOnductivity;
  double cApacity;
  double initTemp;
  BlockOptionData():
    oRder(-1),
    cOnductivity(-1),
    cApacity(-1),
    initTemp(0) {}
};

struct MonitorPostProc: public FEMethod {

  MoFEM::Interface &mField;
  PostProcVolumeOnRefinedMesh postProc;

  bool iNit;
  int pRT;

  MonitorPostProc(MoFEM::Interface &m_field):
    FEMethod(),mField(m_field),postProc(m_field),iNit(false) {
    
    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 1;
    }
  }

  MoFEMErrorCode preProcess() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode operator()() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode postProcess() {
    MoFEMFunctionBeginHot;
    
    if(!iNit) {
      ierr = postProc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("TEMP"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("TEMP_RATE"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesGradientPostProc("TEMP"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      iNit = true;
    }
    int step;
    ierr = TSGetTimeStepNumber(ts,&step); CHKERRQ(ierr);
    
    if((step)%pRT==0) {
      ierr = mField.loop_finite_elements("DMTHERMAL","THERMAL_FE",postProc); CHKERRQ(ierr);
      std::ostringstream sss;
      sss << "out_thermal_" << step << ".h5m";
      ierr = postProc.writeFile(sss.str().c_str()); CHKERRQ(ierr);
    }
    MoFEMFunctionReturnHot(0);
  }

};


int main(int argc, char *argv[]) {

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
  const char *option;
  option = "";

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  char time_data_file_for_ground_surface[255];
  PetscBool ground_temperature_analys = PETSC_FALSE;
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_ground_analysis_data",
    time_data_file_for_ground_surface,255,&ground_temperature_analys); CHKERRQ(ierr);
  if(ground_temperature_analys) {
    #ifndef __GROUNDSURFACETEMERATURE_HPP
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR to do ground thermal analys MoFEM need to be complilet wiith ADOL-C");
    #endif // __GROUNDSURFACETEMERATURE_HPP
  }

  DMType dm_name = "DMTHERMAL";
  ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);
  //craete dm instance
  DM dm;
  ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
  ierr = DMSetType(dm,dm_name);CHKERRQ(ierr);

  //create MoAB
  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);
  //create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level (this allow to set refinement levels for h-adaptivity)
  //onlt one level is used in this example
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields H1 space rank 1
  ierr = m_field.add_field("TEMP",H1,AINSWORTH_LEGENDRE_BASE,1,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("TEMP_RATE",H1,AINSWORTH_LEGENDRE_BASE,1,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);

  //Add field H1 space rank 3 to approximate gemetry using heierachical basis
  //For 10 node tets, before use, gemetry is projected on that field (see below)
  ierr = m_field.add_field(
    "MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO
  ); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field (root_mesh, i.e. on all mesh etities fields are approx.)
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"TEMP"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"TEMP_RATE"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //for simplicity of example to all entities is applied the same order
  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(root_set,MBTET,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP_RATE",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP_RATE",1); CHKERRQ(ierr);

  //gemetry approximation is set to 2nd oreder
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  // configure blocks by parsing config file
  // it allow to set approximation order for each block independettly
  PetscBool block_config;
  char block_config_file[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_block_config",block_config_file,255,&block_config); CHKERRQ(ierr);
  std::map<int,BlockOptionData> block_data;
  bool solar_radiation = false;
  if(block_config) {
    try {
      ifstream ini_file(block_config_file);
      //std::cerr << block_config_file << std::endl;
      po::variables_map vm;
      po::options_description config_file_options;
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        std::ostringstream str_order;
        str_order << "block_" << it->getMeshsetId() << ".temperature_order";
        config_file_options.add_options()
        (str_order.str().c_str(),po::value<int>(&block_data[it->getMeshsetId()].oRder)->default_value(order));
        std::ostringstream str_cond;
        str_cond << "block_" << it->getMeshsetId() << ".heat_conductivity";
        config_file_options.add_options()
        (str_cond.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].cOnductivity)->default_value(-1));
        std::ostringstream str_capa;
        str_capa << "block_" << it->getMeshsetId() << ".heat_capacity";
        config_file_options.add_options()
        (str_capa.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].cApacity)->default_value(-1));
        std::ostringstream str_init_temp;
        str_init_temp << "block_" << it->getMeshsetId() << ".initail_temperature";
        config_file_options.add_options()
        (str_init_temp.str().c_str(),po::value<double>(&block_data[it->getMeshsetId()].initTemp)->default_value(0));
      }
      config_file_options.add_options()
      ("climate_model.solar_radiation",po::value<bool>(&solar_radiation)->default_value(false));
      po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
      store(parsed,vm);
      po::notify(vm);
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
        if(block_data[it->getMeshsetId()].oRder == -1) continue;
        if(block_data[it->getMeshsetId()].oRder == order) continue;
        PetscPrintf(PETSC_COMM_WORLD,"Set block %d oRder to %d\n",it->getMeshsetId(),block_data[it->getMeshsetId()].oRder);
        Range block_ents;
        rval = moab.get_entities_by_handle(it->meshset,block_ents,true); CHKERRG(rval);
        Range ents_to_set_order;
        ierr = moab.get_adjacencies(block_ents,3,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ(ierr);
        ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
        ierr = moab.get_adjacencies(block_ents,2,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ(ierr);
        ierr = moab.get_adjacencies(block_ents,1,false,ents_to_set_order,moab::Interface::UNION); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents_to_set_order,"TEMP",block_data[it->getMeshsetId()].oRder); CHKERRQ(ierr);
        ierr = m_field.set_field_order(ents_to_set_order,"TEMP_RATE",block_data[it->getMeshsetId()].oRder); CHKERRQ(ierr);
      }
      std::vector<std::string> additional_parameters;
      additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
      for(std::vector<std::string>::iterator vit = additional_parameters.begin();
      vit!=additional_parameters.end();vit++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARRNING Unrecognised option %s\n",vit->c_str()); CHKERRQ(ierr);
      }
    } catch (const std::exception& ex) {
      std::ostringstream ss;
      ss << ex.what() << std::endl;
      SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
    }
  }

  //this default class to calculate thermal elements
  ThermalElement thermal_elements(m_field);
  ierr = thermal_elements.addThermalElements("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalFluxElement("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalConvectionElement("TEMP"); CHKERRQ(ierr);
  ierr = thermal_elements.addThermalRadiationElement("TEMP"); CHKERRQ(ierr);
  //add rate of temerature to data field of finite element
  ierr = m_field.modify_finite_element_add_field_data("THERMAL_FE","TEMP_RATE"); CHKERRQ(ierr);
  //and temperature element default element operators at integration (gauss) points
  ierr = thermal_elements.setTimeSteppingProblem("TEMP","TEMP_RATE"); CHKERRQ(ierr);

  //set block material data from opetion file
  std::map<int,ThermalElement::BlockData>::iterator mit;
  mit = thermal_elements.setOfBlocks.begin();
  for(;mit!=thermal_elements.setOfBlocks.end();mit++) {
    //std::cerr << mit->first << std::endl;
    //std::cerr << block_data[mit->first].cOnductivity  << " " << block_data[mit->first].cApacity << std::endl;
    if(block_data[mit->first].cOnductivity != -1) {
      PetscPrintf(PETSC_COMM_WORLD,"Set block %d heat conductivity to %3.2e\n",
      mit->first,block_data[mit->first].cOnductivity);
      for(int dd = 0;dd<3;dd++) {
        mit->second.cOnductivity_mat(dd,dd) = block_data[mit->first].cOnductivity;
      }
    }
    if(block_data[mit->first].cApacity != -1) {
      PetscPrintf(PETSC_COMM_WORLD,"Set block %d heat capacity to %3.2e\n",
      mit->first,block_data[mit->first].cApacity);
      mit->second.cApacity = block_data[mit->first].cApacity;
    }
  }

  #ifdef __GROUNDSURFACETEMERATURE_HPP
  GroundSurfaceTemerature ground_surface(m_field);
  CrudeClimateModel time_data(time_data_file_for_ground_surface);
  GroundSurfaceTemerature::PreProcess exectuteGenericClimateModel(&time_data);
  if(ground_temperature_analys) {
    ierr = ground_surface.addSurfaces("TEMP");   CHKERRQ(ierr);
    ierr = ground_surface.setOperators(&time_data,"TEMP"); CHKERRQ(ierr);
  }
  #endif //__GROUNDSURFACETEMERATURE_HPP

  //build database, i.e. declare dofs, elements and ajacencies

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //priject 10 node tet approximation of gemetry on hierarhical basis
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    if(block_data[it->getMeshsetId()].initTemp!=0) {
      Range block_ents;
      rval = moab.get_entities_by_handle(it->meshset,block_ents,true); CHKERRG(rval);
      Range vertices;
      rval = moab.get_connectivity(block_ents,vertices,true); CHKERRG(rval);
      ierr = m_field.set_field(block_data[it->getMeshsetId()].initTemp,MBVERTEX,vertices,"TEMP"); CHKERRQ(ierr);
    }
  }

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  // delete old temerature recorded series
  SeriesRecorder *recorder_ptr;
  ierr = m_field.getInterface(recorder_ptr); CHKERRQ(ierr);
  if(recorder_ptr->check_series("THEMP_SERIES")) {
    /*for(_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder_ptr,"THEMP_SERIES",sit)) {
      ierr = recorder_ptr->load_series_data("THEMP_SERIES",sit->get_step_number()); CHKERRQ(ierr);
    }*/
    ierr = recorder_ptr->delete_recorder_series("THEMP_SERIES"); CHKERRQ(ierr);
  }

  //set dm data structure which created mofem data structures
  ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRQ(ierr);
  ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
  //add elements to dm
  ierr = DMMoFEMAddElement(dm,"THERMAL_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_FLUX_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_CONVECTION_FE"); CHKERRQ(ierr);
  ierr = DMMoFEMAddElement(dm,"THERMAL_RADIATION_FE"); CHKERRQ(ierr);
  #ifdef __GROUNDSURFACETEMERATURE_HPP
  if(ground_temperature_analys) {
    ierr = DMMoFEMAddElement(dm,"GROUND_SURFACE_FE"); CHKERRQ(ierr);
  }
  #endif //__GROUNDSURFACETEMERATURE_HPP

  ierr = DMSetUp(dm); CHKERRQ(ierr);

  //create matrices
  Vec T,F;
  ierr = DMCreateGlobalVector_MoFEM(dm,&T); CHKERRQ(ierr);
  ierr = VecDuplicate(T,&F); CHKERRQ(ierr);
  Mat A;
  ierr = DMCreateMatrix_MoFEM(dm,&A); CHKERRQ(ierr);

  DirichletTemperatureBc dirichlet_bc(m_field,"TEMP",A,T,F);
  ThermalElement::UpdateAndControl update_velocities(m_field,"TEMP","TEMP_RATE");
  ThermalElement::TimeSeriesMonitor monitor(m_field,"THEMP_SERIES","TEMP");
  MonitorPostProc post_proc(m_field);

  //Initialize data with values save of on the field
  ierr = VecZeroEntries(T); CHKERRQ(ierr);
  ierr = DMoFEMMeshToLocalVector(dm,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = DMoFEMPreProcessFiniteElements(dm,&dirichlet_bc); CHKERRQ(ierr);
  ierr = DMoFEMMeshToGlobalVector(dm,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  //preprocess
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&update_velocities,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&dirichlet_bc,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,&dirichlet_bc,NULL); CHKERRQ(ierr);
  #ifdef __GROUNDSURFACETEMERATURE_HPP
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&exectuteGenericClimateModel,NULL); CHKERRQ(ierr);
  { // add preporcessor, calculating angle on which sun ray on the surface
    if(solar_radiation) {
      boost::ptr_vector<GroundSurfaceTemerature::SolarRadiationPreProcessor>::iterator it,hi_it;
      it = ground_surface.preProcessShade.begin();
      hi_it = ground_surface.preProcessShade.end();
      for(;it!=hi_it;it++) {
        ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&*it,NULL); CHKERRQ(ierr);
      }
    }
  }
  #endif //__GROUNDSURFACETEMERATURE_HPP

  //loops rhs
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_FE",&thermal_elements.feRhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_FLUX_FE",&thermal_elements.feFlux,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_CONVECTION_FE",&thermal_elements.feConvectionRhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIFunction(dm,"THERMAL_RADIATION_FE",&thermal_elements.feRadiationRhs,NULL,NULL); CHKERRQ(ierr);
  #ifdef __GROUNDSURFACETEMERATURE_HPP
  if(ground_temperature_analys) {
    ierr = DMMoFEMTSSetIFunction(dm,"GROUND_SURFACE_FE",&ground_surface.getFeGroundSurfaceRhs(),NULL,NULL); CHKERRQ(ierr);
  }
  #endif //__GROUNDSURFACETEMERATURE_HPP


  //loops lhs
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_FE",&thermal_elements.feLhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_CONVECTION_FE",&thermal_elements.feConvectionLhs,NULL,NULL); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,"THERMAL_RADIATION_FE",&thermal_elements.feRadiationLhs,NULL,NULL); CHKERRQ(ierr);
  #ifdef __GROUNDSURFACETEMERATURE_HPP
  if(ground_temperature_analys) {
    ierr = DMMoFEMTSSetIJacobian(dm,"GROUND_SURFACE_FE",&ground_surface.getFeGroundSurfaceLhs(),NULL,NULL); CHKERRQ(ierr);
  }
  #endif //__GROUNDSURFACETEMERATURE_HPP

  //postprocess
  ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,NULL,&dirichlet_bc); CHKERRQ(ierr);
  ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,NULL,&dirichlet_bc); CHKERRQ(ierr);

  TsCtx *ts_ctx;
  DMMoFEMGetTsCtx(dm,&ts_ctx);
  //add monitor operator
  ts_ctx->get_postProcess_to_do_Monitor().push_back(&monitor);
  ts_ctx->get_postProcess_to_do_Monitor().push_back(&post_proc);

  //create time solver
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  ierr = TSSetIFunction(ts,F,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,A,A,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  //add monitor to TS solver
  ierr = TSMonitorSet(ts,f_TSMonitorSet,ts_ctx,PETSC_NULL); CHKERRQ(ierr); // !!!

  ierr = recorder_ptr->add_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);
  //start to record
  ierr = recorder_ptr->initialize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
  ierr = TSSetDM(ts,dm); CHKERRQ(ierr);

  ierr = TSSolve(ts,T); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  //end recoder
  ierr = recorder_ptr->finalize_series_recorder("THEMP_SERIES"); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);

  // save solution, if boundary conditions are defined you can use that file in mechanical problem
  // to calculate thermal stresses
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-dm_is_partitioned",&is_partitioned,PETSC_NULL); CHKERRQ(ierr);
  if(is_partitioned) {
    rval = moab.write_file("solution.h5m"); CHKERRG(rval);
  } else {
    if(pcomm->rank()==0) {
      rval = moab.write_file("solution.h5m"); CHKERRG(rval);
    }
  }


  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&T); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  ierr = TSDestroy(&ts); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
