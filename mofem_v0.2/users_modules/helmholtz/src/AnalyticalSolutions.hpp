/* \file best_approximation_incident_wave.cpp
 
  Analytical solutions
  
 */

/* This file is part of MoFEM.
 * MoFEM is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.

 * MoFEM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>. 
*/

/// List of analytical solution
enum AnalyticalSolutionTypes {
  INCIDENT_WAVE, PLANE_WAVE
};

/// Line command list of analytical solutions
const char *analytical_solution_types[2] = { 
  "incident_wave", 
  "plane_wave" 
};

/** Calculate the analytical solution of impinging wave on sphere

  See paper: 
  Exact solution of Impinging sphere from Acoustic isogeometric boundary element analysis by R.N. Simpson etc.

  */
struct IncidentWaveAnalyticalSolution {
  
   ublas::vector<double> result1;
   double wAvenumber;
   bool useReal;
   
   IncidentWaveAnalyticalSolution(double wavenumber,bool use_real):
     wAvenumber(wavenumber),useReal(use_real) {}
   ~IncidentWaveAnalyticalSolution() {}
   
  ublas::vector<double>& operator()(double x, double y, double z) {

    double R = sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0)); //radius
    //Incident wave in X direction.
    double cos_theta = y/R;

    const double k = wAvenumber;  //Wave number
    const double a = 0.5;         //radius of the sphere,wait to modify by user
    //const double a = 1.0;
    const double const1 = k * a;
    double const2 = k * R;
    
    const complex< double > i( 0.0, 1.0 );
    
    // magnitude of incident wave
    //const double phi_incident_mag = 1.0;
    
    const double tol = 1.0e-10;
    //double max = 0.0;
    //double min = 999999.0;
    
    complex< double > result = 0.0;
    complex< double > prev_result;
    
    double error = 100.0;
    unsigned int n = 0; //initialized the infinite series loop
    
    while( error > tol )  //finding the acoustic potential in one single point.
    {

      //The derivative of bessel function
      double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );
      //The derivative of Hankel function
      complex<double> hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
      double Pn = legendre_p( n, cos_theta );
      complex<double> hn = sph_hankel_1( n, const2 );  //S Hankel first kind function
      prev_result = result;
      result -= pow( i, n ) * ( 2.0 * n + 1.0 ) * jn_der / hn_der * Pn * hn;
      error = abs( abs( result ) - abs( prev_result ) );
      ++n;

    }
    
    if(useReal) {
      result1.resize(1);
      result1[0] = std::real(result);  
      return result1;
    } else {
      result1.resize(1);
      result1[0] = std::imag(result);  
      return result1;
    }
    
  }
  
};

/** \brief Calculate the analytical solution of plane wave guide propagating in direction theta
  
  Paper: ????

*/
struct PlaneWave {
  
   ublas::vector<double> result1;
   double wAvenumber;
   double tHeta;
   bool useReal;
   
   PlaneWave(double wavenumber,double theta,bool use_real):
     wAvenumber(wavenumber),tHeta(theta),useReal(use_real) {}
   ~PlaneWave() {}
   
  ublas::vector<double>& operator()(double x, double y, double z) {
    
    const double k = wAvenumber;  //Wave number

    const complex< double > i( 0.0, 1.0 );
    complex< double > result = 0.0;
    // magnitude of incident wave
    //const double phi_incident_mag = 3.0;
    
    //const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
    //const complex< double > total_field = inc_field + result;
    //ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
    
    /* cube */
    result = exp(i*(k*cos(tHeta)*x+k*sin(tHeta)*y));
    /* cube */
    
    if(useReal) {
      result1.resize(1);
      result1[0] = std::real(result);  
      return result1;
    } else {
      result1.resize(1);
      result1[0] = std::imag(result);  
      return result1;
    }
    
  }
  
};

