/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <moab/AdaptiveKDTree.hpp>
extern "C" {
  #include <spa.h>
}

#include <time.h>
#include <adolc/adolc.h> 
//#include <ThermalElement.hpp>
#include<moab/Skinner.hpp>

#include <GenricClimateModel.hpp>
#include <GroundSurfaceTemerature.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct MyTimeData: public GenricClimateModel {

  MyTimeData(): GenricClimateModel() {

    spaData.delta_ut1 = 0;    	// Fractional second difference between UTC and UT which is used
				// to adjust UTC for earth's irregular rotation rate and is derived
				// from observation only and is reported in this bulletin:
				// http://maia.usno.navy.mil/ser7/ser7.dat,
				// where delta_ut1 = DUT1
				// valid range: -1 to 1 second (exclusive), error code 17

    spaData.delta_t = 0;     	// Difference between earth rotation time and terrestrial time
				// It is derived from observation only and is reported in this
				// bulletin: http://maia.usno.navy.mil/ser7/ser7.dat,
				// where delta_t = 32.184 + (TAI-UTC) - DUT1
				// valid range: -8000 to 8000 seconds, error code: 7

    spaData.timezone = 0;     	// Observer time zone (negative west of Greenwich)
				// valid range: -18   to   18 hours,   error code: 8

    spaData.elevation = 10;    	// Observer elevation [meters]
				// valid range: -6500000 or higher meters,    error code: 11

    spaData.pressure = 1013.25; // Annual average local pressure [millibars]
				// valid range:    0 to 5000 millibars,       error code: 12

    spaData.temperature = 20;  	// Annual average local temperature [degrees Celsius]
				// valid range: -273 to 6000 degrees Celsius, error code; 13


    spaData.slope = 0;        	// Surface slope (measured from the horizontal plane)
				// valid range: -360 to 360 degrees, error code: 14

    spaData.azm_rotation = 0; 	// Surface azimuth rotation (measured from south to projection of
				// surface normal on horizontal plane, negative east)
				// valid range: -360 to 360 degrees, error code: 15

    spaData.atmos_refract = 0.5667; // Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
				// valid range: -5   to   5 degrees, error code: 16

    //Longest day (Solstice)
    spaData.year = 2014;	// 4-digit year,      valid range: -2000 to 6000, error code: 1
    spaData.month = 6;          // 2-digit month,         valid range: 1 to  12,  error code: 2
    spaData.day = 21;           // 2-digit day,           valid range: 1 to  31,  error code: 3

    spaData.hour = 0;        	// Observer local hour,   valid range: 0 to  24,  error code: 4
    spaData.minute = 0;         // Observer local minute, valid range: 0 to  59,  error code: 5
    spaData.second = 0;         // Observer local second, valid range: 0 to <60,  error code: 6	

    //This is London
    spaData.longitude = 0.1275;   // Observer longitude (negative west of Greenwich)
				  // valid range: -180  to  180 degrees, error code: 9

    spaData.latitude = 51.5072;    // Observer latitude (negative south of equator)
				  // valid range: -90   to   90 degrees, error code: 10
  }

  spa_data spaData;

  PetscErrorCode set(double t = 0) {
    PetscFunctionBegin;

    spaData.function = SPA_ZA_RTS;

    int r;
    r = spa_calculate(&spaData);
    if(r) {
      SETERRQ1(PETSC_COMM_SELF,1,"wrong input data for solar position calulator error codde %d",r);
    }
    
    zenith = spaData.zenith;
    azimuth = spaData.azimuth;

    PetscPrintf(PETSC_COMM_WORLD,
      "Suntransit %3.2f Sunrise %3.2f Sunset %3.2f\n" ,
      spaData.suntransit,spaData.sunrise,spaData.sunset);

    PetscFunctionReturn(0);
  }

};

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  const char *option;
  option = "";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields H1 space rank 1
  ierr = m_field.add_field("TEMP",H1,1); CHKERRQ(ierr);

  //Add field H1 space rank 3 to approximate gemetry using heierachical basis
  //For 10 node tets, before use, gemetry is projected on that field (see below)
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set(); 
  //add entities to field (root_mesh, i.e. on all mesh etities fields are approx.)
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"TEMP"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //for simplicity of example to all entities is applied the same order
  ierr = m_field.set_field_order(root_set,MBTET,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"TEMP",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"TEMP",1); CHKERRQ(ierr);

  //gemetry approximation is set to 2nd oreder
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  GroundSurfaceTemerature ground_surface(m_field);
  ground_surface.addSurfaces("TEMP"); 

  //build database, i.e. declare dofs, elements and ajacencies

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //priject 10 node tet approximation of gemetry on hierarhical basis 
  Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
  ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  MyTimeData time_data;
  ierr = ground_surface.setOperators(1,&time_data,"TEMP"); CHKERRQ(ierr);
  GroundSurfaceTemerature::Shade *shade_ptr = &*ground_surface.preProcessShade.begin();

  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,tets);  CHKERRQ(ierr);
  shade_ptr->getSkin(tets);

  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET,meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(meshset,shade_ptr->sKin); CHKERR_PETSC(rval);

  struct tm start_time;
  start_time.tm_sec = time_data.spaData.second;
  start_time.tm_min = time_data.spaData.minute;
  start_time.tm_hour = time_data.spaData.hour;
  start_time.tm_mday = time_data.spaData.day;
  start_time.tm_mon = time_data.spaData.month-1;
  start_time.tm_year = time_data.spaData.year-1900;
  time_t t0 = mktime(&start_time);
  time_t t = t0;
  for(;t<t0+60*60*24;t+=4*60*60) {

    struct tm current_time;
    current_time = *gmtime(&t);//localtime(&t);
    time_data.spaData.second = current_time.tm_sec;
    time_data.spaData.minute = current_time.tm_min;
    time_data.spaData.hour = current_time.tm_hour;
    time_data.spaData.day = current_time.tm_mday;
    time_data.spaData.month = current_time.tm_mon+1;
    time_data.spaData.year = current_time.tm_year+1900;

    PetscPrintf(PETSC_COMM_WORLD,"%s",asctime(&current_time));
        
    ierr = time_data.set(); CHKERRQ(ierr);
    ierr = shade_ptr->preProcess(); CHKERRQ(ierr);

    ostringstream ss;
    ss << "out_" << t << ".vtk";

    rval = moab.write_file(ss.str().c_str(),"VTK","",&meshset,1); CHKERR_PETSC(rval); CHKERR_PETSC(rval);

  }

  //define problems
  ierr = m_field.add_problem("GROUND_SURFACE"); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("GROUND_SURFACE","GROUND_SURFACE_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("GROUND_SURFACE",bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);
  //partition
  ierr = m_field.partition_problem("GROUND_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("GROUND_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("GROUND_SURFACE"); CHKERRQ(ierr);

  for(_IT_GET_DOFS_FIELD_BY_NAME_AND_TYPE_FOR_LOOP_(m_field,"TEMP",MBVERTEX,dof)) {
    dof->get_FieldData() = 20;
  }

  //create matrices
  Vec F;
  ierr = m_field.VecCreateGhost("GROUND_SURFACE",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("GROUND_SURFACE",&Aij); CHKERRQ(ierr);

  ground_surface.feGroundSurfaceRhs.ts_F = F;
  ground_surface.feGroundSurfaceRhs.ts_a = 1;
  ground_surface.feGroundSurfaceRhs.ts_B = Aij;
  ground_surface.feGroundSurfaceRhs.ts_a = 1;

  ground_surface.feGroundSurfaceLhs.ts_F = F;
  ground_surface.feGroundSurfaceLhs.ts_a = 1;
  ground_surface.feGroundSurfaceLhs.ts_B = Aij;
  ground_surface.feGroundSurfaceLhs.ts_a = 1;

  ierr = m_field.loop_finite_elements("GROUND_SURFACE","GROUND_SURFACE_FE",ground_surface.getFeGroundSurfaceRhs()); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

  ierr = m_field.loop_finite_elements("GROUND_SURFACE","GROUND_SURFACE_FE",ground_surface.getFeGroundSurfaceLhs()); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"ground_surface_temperature.txt",&viewer); CHKERRQ(ierr);

  //ierr = VecChop(F,1e-4); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecView(F,viewer); CHKERRQ(ierr);
  
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);
  MatChop(Aij,1e-4);
  //MatView(Aij,PETSC_VIEWER_STDOUT_WORLD);
  MatView(Aij,viewer);
  //std::string wait;
  //std::cin >> wait;

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);


  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}



