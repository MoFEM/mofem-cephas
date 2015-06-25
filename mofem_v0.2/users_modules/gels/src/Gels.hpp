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

#ifndef __GEL_HPP__
#define __GEL_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

template<typename TYPE>
struct Gel {

  ublas::matrix<TYPE> stressAlpha, stressBeta, strainHatDot, stressBetaHat;
  ublas::vector<TYPE> fLux;
  ublas::matrix<TYPE> traceStrainMainMat;
  

  
  virtual PetscErrorCode calcualteStressAlpha(ublas::matrix<TYPE> strainMain, double vAlpha, double gAlpha) {
    PetscFunctionBegin;
    TYPE traceStrainMain;
    traceStrainMain=strainMain(0,0)+strainMain(1,1)+strainMain(2,2);
//    cout<<"traceStrainMain  "<<traceStrainMain<<endl;
    traceStrainMainMat.resize(3,3); traceStrainMainMat.clear();
    for(int ii=0; ii<3; ii++){
      traceStrainMainMat(ii,ii)=traceStrainMain;
    }
//    cout<<"traceStrainMainMat1  "<<traceStrainMainMat<<endl;
    stressAlpha=2*gAlpha*(strainMain+(vAlpha/(1-2*vAlpha))*traceStrainMainMat);
    cout<<"stressAlpha  "<<stressAlpha<<endl;
    PetscFunctionReturn(0);
  }


  virtual PetscErrorCode calcualteStressBeta(ublas::matrix<TYPE> strainMain, ublas::matrix<TYPE> strainHat, double vBeta, double gBeta) {
    PetscFunctionBegin;
    TYPE traceStrainHat;
    ublas::matrix<TYPE> traceStrianHatMat;
    traceStrainHat=strainHat(0,0)+strainHat(1,1)+strainHat(2,2);
    traceStrianHatMat.resize(3,3); traceStrianHatMat.clear();
    for(int ii=0; ii<3; ii++){
      traceStrianHatMat(ii,ii)=traceStrainHat;
    }
//    cout<<"traceStrainMainMat2  "<<traceStrainMainMat<<endl;
    stressBeta=2*gBeta*((strainMain-strainHat)+(vBeta/(1-2*vBeta))*(traceStrainMainMat-traceStrianHatMat));
    cout<<"stressBeta  "<<stressBeta<<endl;
    PetscFunctionReturn(0);
  }


  virtual PetscErrorCode calcualteStrainHatDot(double vBetaHat, double gBetaHat) {
    PetscFunctionBegin
    cout<<"\n\n\n\nstressBeta1  "<<stressBeta<<endl;
    TYPE traceStressBeta;
    ublas::matrix<TYPE> traceStressBetaMat;
    traceStressBeta=stressBeta(0,0)+stressBeta(1,1)+stressBeta(2,2);
    cout<<"traceStressBeta  "<<traceStressBeta<<endl;
    traceStressBetaMat.resize(3,3); traceStressBetaMat.clear();
    for(int ii=0; ii<3; ii++){
      traceStressBetaMat(ii,ii)=traceStressBeta;
    }
    strainHatDot=(1/(2*gBetaHat))*(stressBeta-(vBetaHat/(1+vBetaHat))*traceStressBetaMat);
    cout<<"strainHatDot  "<<strainHatDot<<endl;
    PetscFunctionReturn(0);
  }
  
  
  virtual PetscErrorCode calcualteStressBetaHat(double mu, double omega) {
    PetscFunctionBegin
    stressBetaHat.resize(3,3);  stressBetaHat.clear();
    for(int ii=0; ii<3; ii++){
      stressBetaHat(ii,ii)=1;
    }
    stressBetaHat=(mu/omega)*stressBetaHat;
    cout<<"\n\n\nstressBetaHat  "<<stressBetaHat<<endl;
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode calcualteFlux(double permeability, double viscosity, double omega, ublas::vector<TYPE> gradientMu) {
    PetscFunctionBegin
    fLux=-(permeability/(viscosity*omega*omega))*gradientMu;
    cout<<"\n\n\nfLux  "<<fLux<<endl;
    PetscFunctionReturn(0);
  }

  
};

#endif //__GEL_HPP__
