/** \file GroundSurfaceTemerature.hpp 
 * \brief Operators and data structures for thermal analys
 *
 * Implementation of boundary conditions for ground temerature
 * 
 */

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

#ifndef __GROUNDSURFACETEMERATURE_HPP
#define __GROUNDSURFACETEMERATURE_HPP

struct GroundSurfaceTemerature: public ThermalElement {

  // sigma = 5.67037321×10−8 (J/s) m−2 K−4 
  // sigma = 8165.3×10−8 (J/day) m−2 K−4
  // sigma = 0.081653 (kJ/day) m−2 K−4
  const double sigma;	//< Stefan–Boltzmann constant (kJ/K4/m2/day)

  GroundSurfaceTemerature(): ThermalElement(),  sigma(0.081653) {};

  struct Parameters {	
    double alpha; 	//< Solar albedo
    double Cfc; 	//< Surface heat/moisture transfer coefficient for forced convection
    double Cnc;		//< Coefficient for natural convection
    double CSh;		//< Wind sheltering coefficient
    double eps;		//< Pavement emissivity
    double rhoCp;	//< Density specific heat pavement (J/m3/°C)
  };

  struct Asphalt: public Parameters {
    alpha = 0.12;
    Cfc = 0.0015;
    Cnc = 0.0015;
    CSh = 1.;
    eps = 0.94;
    rhoCp = 2.0e06;
  };

  struct Concrete: public Parameters {
    alpha = 0.20;
    Cfc = 0.0015;
    Cnc = 0.0015;
    CSh = 1.;
    eps = 0.94;
    rhoCp = 2.0e06;
  };

  struct BareSoil: public Parameters {
    alpha = 0.15;
    Cfc = 0.003;
    Cnc = 0.0015;
    CSh = 1.;
    eps = 0.94;
    rhoCp = 2.0e06;
  };

  struct OpJacobian {
  };


};

#endif //__GROUNDSURFACETEMERATURE_HPP

