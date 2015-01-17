/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * This file is part of MoFEM.
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

#ifndef __GENRICCLIMATEMODEL_HPP
#define __GENRICCLIMATEMODEL_HPP

struct GenricClimateModel {

  double T0; // reference temperature (K)
  double e0; // reference saturation vapor pressure (es at a certain temp, usually 0 deg C) (Pa)
  double Rv; // gas constant for water vapor (J*K/Kg)
  double Lv; // latent heat of vaporization of water (J)

  double u10;		//< wind at high 10m (m/s)	
  double CR; 		//< cloudness factor (0–1, dimensionless)
  double Ta;		//< air temperature (C)
  double Td;		//< dew point temperature (C)
  double P;		//< pressure
  double Rs; 		//< observed solar radiation (W/m2)

  double zenith;       //topocentric zenith angle [degrees]
  double azimuth;      //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]

  /** \brief Clausius-Clapeyron equation
    */
  template <typename TYPE> 
  TYPE calulateVapourPressureClausiusClapeyron(TYPE T) { 
    return e0*exp((Lv/Rv)*((1./T0)-(1./(T+T0))));
  }

  template <typename TYPE> 
  TYPE calulateVapourPressureTetenFormula(TYPE T) {
    const double b = 17.2694;
    const double T1 = 273.15;
    const double T2 = 35.86;
    return e0*exp(b*(T+T0-T1)/(T+T0-T2));
  }

  template <typename TYPE> 
  TYPE calulateVapourPressure(TYPE T) {
    //This use Tetent's formula by default
    return calulateVapourPressureTetenFormula(T);
  }

  template <typename TYPE>
  TYPE calulateMixingRatio(TYPE T,double P) {
    const double eps = 0.622;
    TYPE e = calulateVapourPressure(T);
    return eps*T/(P-T);
  }

  template <typename TYPE>   
  double calculateAbsoluteVirtualTempertaure(TYPE T,double P) {
    const double c = 0.379;
    double e = calulateVapourPressure(T);
    return (T+T0)/(1-c*e/P);
  }

  GenricClimateModel() {

    T0 = 273.15; // reference temperature (K)
    e0 = 611; // reference saturation vapor pressure (es at a certain temp, usually 0 deg C) (Pa)
    Rv = 461.5; // gas constant for water vapor (J*K/Kg)
    Lv = 2.5e6; // latent heat of vaporization of water (J)

    u10 = 0;			//< wind at high 10m (m/s)	
    CR = 0; 			//< cloudness factor (0–1, dimensionless)
    Ta = 10;			//< air temperature (C)
    Td = 5;			//< dew point temperature (C)
    P = 101325;		//< pressure
    Rs = 1361; 		//< observed solar radiation (W/m2)

    zenith = 0;       //topocentric zenith angle [degrees]
    azimuth = 0;      //topocentric azimuth angle (eastward from north) [for navigators and solar radiation]

  }

  virtual PetscErrorCode set(double t = 0) = 0;
};

#endif //__GENRICCLIMATEMODEL_HPP
