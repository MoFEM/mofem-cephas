/** \file gel_constitutive_equation_test.cpp
  \brief Atom test verifying implementation of gel constitutive equations
  \ingroup gel
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
#include <PostProcOnRefMesh.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adolc/adolc.h>
#include <Gels.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

ErrorCode rval;
PetscErrorCode ierr;
static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  // Create gel instance
  Gel gel(m_field);

  Gel::BlockMaterialData material_data;

  // Set material parameters
  material_data.gAlpha = 1;
  material_data.vAlpha = 0.3;
  material_data.gBeta = 1;
  material_data.vBeta = 0.3;
  material_data.gBetaHat = 1;
  material_data.vBetaHat = 0.3;
  material_data.oMega = 1;
  material_data.vIscosity = 1;
  material_data.pErmeability = 2.;

  Gel::ConstitutiveEquation<double> ce(material_data);

  // Set state variables
  ublas::matrix<double> &F = ce.F;
  F.resize(3,3);
  F.clear();
  F(0,0) = 0.01;  F(0,1) = 0.02;  F(0,2) = 0.03;
  F(1,0) = 0.04;  F(1,1) = 0.05;  F(1,2) = 0.06;
  F(2,0) = 0.07;  F(2,1) = 0.08;  F(2,2) = 0.09;

  ublas::matrix<double> &F_dot = ce.FDot;
  F_dot.resize(3,3);
  F_dot.clear();
  F_dot(0,0) = 0.01;  F_dot(0,1) = 0.02;  F_dot(0,2) = 0.03;
  F_dot(1,0) = 0.04;  F_dot(1,1) = 0.05;  F_dot(1,2) = 0.06;
  F_dot(2,0) = 0.07;  F_dot(2,1) = 0.08;  F_dot(2,2) = 0.09;

  ublas::matrix<double> &strainHat = ce.strainHat;
  strainHat.resize(3,3);
  strainHat.clear();
  strainHat(0,0) = 0.01;  strainHat(0,1) = 0.02;  strainHat(0,2) = 0.03;
  strainHat(1,0) = 0.04;  strainHat(1,1) = 0.05;  strainHat(1,2) = 0.06;
  strainHat(2,0) = 0.07;  strainHat(2,1) = 0.08;  strainHat(2,2) = 0.09;

  ublas::matrix<double> &strainHatDot = ce.strainHatDot;
  strainHatDot.resize(3,3);
  strainHatDot.clear();
  strainHatDot(0,0) = 0.01;  strainHatDot(0,1) = 0.02;  strainHatDot(0,2) = 0.03;
  strainHatDot(1,0) = 0.04;  strainHatDot(1,1) = 0.05;  strainHatDot(1,2) = 0.06;
  strainHatDot(2,0) = 0.07;  strainHatDot(2,1) = 0.08;  strainHatDot(2,2) = 0.09;

  ce.mU = 1;
  ublas::vector<double> &gradient_mu = ce.gradientMu;
  gradient_mu.resize(3);
  gradient_mu.clear();
  gradient_mu[0] = 1;
  gradient_mu[1] = 0;
  gradient_mu[2] = 0;

  // Creeate tee like output, printing results on screen and simultaneously
  // writing results to file
  typedef tee_device<ostream, ofstream> TeeDevice;
  typedef stream<TeeDevice> TeeStream;
  ofstream ofs("gel_constitutive_model_test.txt");
  TeeDevice my_tee(cout,ofs);
  TeeStream my_split(my_tee);

  // Do calculations

  ierr = ce.calculateCauchyDefromationTensor(); CHKERRQ(ierr);
  my_split << "C\n" << ce.C << endl << endl;

  ierr = ce.calculateStrainTotal(); CHKERRQ(ierr);
  my_split << "strainTotal\n" << ce.strainTotal << endl << endl;

  ierr = ce.calculateTraceStrainTotalDot(); CHKERRQ(ierr);
  my_split << "traceStrianTotalDot\n" << ce.traceStrainTotalDot << endl << endl;

  ierr = ce.calculateVolumeDot(); CHKERRQ(ierr);
  my_split << "volumeDot\n" << ce.volumeDot << endl << endl;

  ierr = ce.calculateStressAlpha(); CHKERRQ(ierr);
  my_split << "stressAlpha\n" << ce.stressAlpha << endl << endl;

  ierr = ce.calculateStressBeta(); CHKERRQ(ierr);
  my_split << "stressBeta\n" << ce.stressBeta << endl << endl;

  ierr = ce.calculateStressBetaHat(); CHKERRQ(ierr);
  my_split << "stressBeta\n" << ce.stressBetaHat << endl << endl;

  ierr = ce.calculateStressTotal(); CHKERRQ(ierr);
  my_split << "stressTotal\n" << ce.stressTotal << endl << endl;

  ierr = ce.calculateSolventFlux(); CHKERRQ(ierr);
  my_split << "solventFlux\n" << ce.solventFlux << endl << endl;

  ierr = ce.calculateStrainHatFlux(); CHKERRQ(ierr);
  my_split << "strainHatDot\n" << ce.strainHatDot << endl << endl;

  ierr = ce.calculateResidualStrainHat(); CHKERRQ(ierr);
  my_split << "residualStrainHat\n" << ce.residualStrainHat << endl << endl;

  PetscFinalize();

  return 0;
}
