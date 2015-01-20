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

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iterator>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

#include <CrudeClimateModel.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = 
  "-nb_days [days]\n"
  "-my_step_size [days]\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  PetscBool flg;
  PetscScalar nb_days;
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-my_nb_days",&nb_days,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    nb_days = 365;
  }
  PetscScalar time_step;
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-my_step_size",&time_step,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    time_step = 0.2; // 5th part of day 
  }

  CrudeClimateModel time_data("parameters.in");
  ierr = time_data.testCode(60*60*24*nb_days,60*60*24*time_step); CHKERRQ(ierr);

  PetscFinalize();
  return 0;
}



