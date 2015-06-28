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

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adolc/adolc.h>
#include <Gels.hpp>

int main(int argc, char *argv[]) {

  /*cout<<"from Gel main "<<endl;
  Gel<double> Gel_double;

  double vAlpha, gAlpha;

  vAlpha=0.3; gAlpha=1;
   ublas::matrix<double> strainMain;
  strainMain.resize(3,3); strainMain.clear();
  strainMain(0,0)=0.01;  strainMain(0,1)=0.02;  strainMain(0,2)=0.03;
  strainMain(1,0)=0.04;  strainMain(1,1)=0.05;  strainMain(1,2)=0.06;
  strainMain(2,0)=0.07;  strainMain(2,1)=0.08;  strainMain(2,2)=0.09;

  cout<<"strainMain = "<<strainMain<<endl;
  Gel_double.calcualteStressAlpha(strainMain, vAlpha, gAlpha);

  ublas::matrix<double> strainHat;
  strainHat=0.5*strainMain;
  cout<<"strainHat = "<<strainHat<<endl;

  double vBeta, gBeta;
  vBeta=0.3; gBeta=1;
  Gel_double.calcualteStressBeta(strainMain, strainHat, vBeta, gBeta);


  double vBetaHat, gBetaHat;
  vBetaHat=0.3; gBetaHat=1;
  Gel_double.calcualteStrainHatDot(vBetaHat, gBetaHat);

  double mu, omega;
  mu=0.2; omega=3;
  Gel_double.calcualteStressBetaHat(mu, omega);


  double permeability, viscosity;
  permeability=1;
  viscosity=2;
  ublas::vector<double> gradientMu;
  gradientMu.resize(3);
  gradientMu(0)=0.4;  gradientMu(1)=0.3;  gradientMu(2)=0.5;
  Gel_double.calcualteFlux(permeability, viscosity, omega, gradientMu);*/

  return 0;
}
