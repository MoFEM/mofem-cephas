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

void tetcircumcenter_tp(double a[3],double b[3],double c[3], double d[3],
  double circumcenter[3],double *xi,double *eta,double *zeta);
void tricircumcenter3d_tp(double a[3],double b[3],double c[3],
  double circumcenter[3],double *xi,double *eta);

#include <spa.h>

}

#include <time.h>

#include <adolc/adolc.h> 
#include <ThermalElement.hpp>
#include <GroundSurfaceTemerature.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct MyTimeData: public GroundSurfaceTemerature::TimeDependendData {

  MyTimeData(): GroundSurfaceTemerature::TimeDependendData() {

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

    spaData.longitude = 4.25;   // Observer longitude (negative west of Greenwich)
				// valid range: -180  to  180 degrees, error code: 9

    spaData.latitude = 55.8;    // Observer latitude (negative south of equator)
				// valid range: -90   to   90 degrees, error code: 10

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

    spaData.atmos_refract;	// Atmospheric refraction at sunrise and sunset (0.5667 deg is typical)
				// valid range: -5   to   5 degrees, error code: 16


    spaData.year = 2015;	// 4-digit year,      valid range: -2000 to 6000, error code: 1
    spaData.month = 7;          // 2-digit month,         valid range: 1 to  12,  error code: 2
    spaData.day = 10;           // 2-digit day,           valid range: 1 to  31,  error code: 3

    spaData.hour = 0;        	// Observer local hour,   valid range: 0 to  24,  error code: 4
    spaData.minute = 0;         // Observer local minute, valid range: 0 to  59,  error code: 5
    spaData.second = 0;         // Observer local second, valid range: 0 to <60,  error code: 6	

  }

  spa_data spaData;

  PetscErrorCode set() {
    PetscFunctionBegin;

    int r;
    r = spa_calculate(&spaData);
    if(r) {
      SETERRQ1(PETSC_COMM_SELF,1,"wrong input data for solar position calulator error codde %d",r);
    }
    
    zenith = spaData.zenith;
    azimuth = spaData.azimuth;

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

  ThermalElement thermal_elements(m_field);
  GroundSurfaceTemerature ground_surface(thermal_elements);

  MyTimeData time_data;

  Range bare_soil_tris;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
    if(it->get_Cubit_name().compare(0,8,"BARESOIL") == 0) {
      Range tris;
      rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      bare_soil_tris.merge(tris);
    }
  }

  GroundSurfaceTemerature::BareSoil bare_soil(bare_soil_tris);
  GroundSurfaceTemerature::Shade shade(m_field,&time_data,bare_soil);


  Range tets;
  ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTET,tets);  CHKERRQ(ierr);
  shade.getSkin(tets);

  EntityHandle meshset;
  rval = moab.create_meshset(MESHSET_SET,meshset); CHKERR_PETSC(rval);
  rval = moab.add_entities(meshset,shade.sKin); CHKERR_PETSC(rval);

  struct tm start_time;
  start_time.tm_sec = time_data.spaData.second;
  start_time.tm_min = time_data.spaData.minute;
  start_time.tm_hour = time_data.spaData.hour;
  start_time.tm_mday = time_data.spaData.day;
  start_time.tm_mon = time_data.spaData.month-1;
  start_time.tm_year = time_data.spaData.year-1900;

  time_t t0 = mktime(&start_time);
  time_t t = t0;
  for(;t<t0+60*60*24;t+=2*60) {

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
    ierr = shade.preProcess(); CHKERRQ(ierr);

    ostringstream ss;
    ss << "out_" << t << ".vtk";

    rval = moab.write_file(ss.str().c_str(),"VTK","",&meshset,1); CHKERR_PETSC(rval); CHKERR_PETSC(rval);

  }
  
  PetscFinalize();

  return 0;
}



