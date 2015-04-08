/* \file AnalyticalSolutions.hpp
 
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

/// Generic structure for analytical function
struct GenericAnalyticalSolution {

  enum VALUE_TYPE { REAL = 0, IMAG, LAST_VAL_TYPE };

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) = 0;

};

/// List of analytical solution
enum AnalyticalSolutionTypes {
  SPHERE_INCIDENT_WAVE, 
  PLANE_WAVE,
  CYLINDER_INCIDENT_WAVE
};

/// Line command list of analytical solutions
const char *analytical_solution_types[] = { 
  "sphere_incident_wave", 
  "plane_wave",
  "cylinder_incident_wave"
};

/** Calculate the analytical solution of impinging wave on sphere

  See paper: 
  Exact solution of Impinging sphere from Acoustic isogeometric boundary element analysis by R.N. Simpson etc.

  */
struct SphereIncidentWave: public GenericAnalyticalSolution {
  
   vector<ublas::vector<double> > rEsult;
   double wAvenumber;
   
   SphereIncidentWave(double wavenumber): wAvenumber(wavenumber) {}
   virtual ~SphereIncidentWave() {}
   
  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;

    double x2 = pow(x,2.0),y2 = pow(y,2.0);
    double R = sqrt(x2+y2+pow(z,2.0)); 
    double cos_theta = sqrt(x2+y2)/R;

    const double k = wAvenumber;  //Wave number
    const double a = 0.5;         //radius of the sphere,wait to modify by user
    //const double a = 1.0;
    const double const1 = k * a;
    double const2 = k * R;
    
    const complex< double > i( 0.0, 1.0 );
    
    // magnitude of incident wave
    //const double phi_incident_mag = 1.0;
    
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

    rEsult.resize(2);
    rEsult[REAL].resize(1);
    (rEsult[REAL])[0] = std::real(result);
    rEsult[IMAG].resize(1);
    (rEsult[IMAG])[0] = std::imag(result);

    return rEsult;
    
  }
  
};

/** \brief Calculate the analytical solution of plane wave guide propagating in direction theta
  
  FIXME: Paper: ????

*/
struct PlaneWave: public GenericAnalyticalSolution {
  
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double tHeta;
   
  PlaneWave(double wavenumber,double theta):
    wAvenumber(wavenumber),tHeta(theta) {}
  virtual ~PlaneWave() {}
   
  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {
    
    const double k = wAvenumber;  //Wave number

    const complex< double > i( 0.0, 1.0 );
    complex< double > result = 0.0;
    // magnitude of incident wave
    //const double phi_incident_mag = 3.0;
    
    //const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
    //const complex< double > total_field = inc_field + result;
    //ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
    
    result = exp(i*(k*cos(tHeta)*x+k*sin(tHeta)*y));
    
    rEsult.resize(2);
    rEsult[REAL].resize(1);
    (rEsult[REAL])[0] = std::real(result);
    rEsult[IMAG].resize(1);
    (rEsult[IMAG])[0] = std::imag(result);

    return rEsult;

  }
  
};

/** \brief Calculate the analytical solution of impinging wave on sphere

  Paper: 
    The generalized finite element method for Helmholtz equation: Theory, computation, and open problems
    Theofanis Strouboulis, Ivo Babuska, Realino Hidaja

*/
struct CylinderIncidentWave: public GenericAnalyticalSolution {
  
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
   
  CylinderIncidentWave(double wavenumber): wAvenumber(wavenumber) {}
  virtual ~CylinderIncidentWave() {}
   
  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;
    double x2 = pow(x,2.0),y2 = pow(y,2.0);
    double R = sqrt(x2+y2+pow(z,2.0));
    double theta = atan2(y,x);

    const double k = wAvenumber;  //Wave number
    const double a = 0.5;         //radius of the sphere,wait to modify by user
    //const double a = 1.0;
    const double const1 = k * a;
    double const2 = k * R;
    
    const complex< double > i( 0.0, 1.0 );
    
    // magnitude of incident wave
    //const double phi_incident_mag = 1.0;
    
    complex< double > result = 0.0;
    complex< double > prev_result;
    
    double error = 100.0;
    unsigned int n = 1; //initialized the infinite series loop
    
    double Jn_der_zero = ( - cyl_bessel_j( 1, const1 ));  
    complex< double > Hn_der_zero = ( - cyl_hankel_1( 1, const1 ));
    complex< double >Hn_zero = cyl_hankel_1( 0, const2 );  //S Hankel first kind function

    //n=0;
    result -= (Jn_der_zero * Hn_zero)/Hn_der_zero;
    
    
    while( error > tol )  //finding the acoustic potential in one single point.
    {

      prev_result = result;
      //The derivative of bessel function
      double Jn_der = (n / const1 * cyl_bessel_j( n, const1 ) - cyl_bessel_j( n + 1, const1 ));  
      //The derivative of Hankel function
      complex< double > Hn_der = (n / const1 * cyl_hankel_1( n, const1 ) - cyl_hankel_1( n + 1, const1 ));
      
      complex< double >Hn = cyl_hankel_1( n, const2 );  //S Hankel first kind function
      
      result -= 2.0 * pow( i, n ) * ( (Jn_der*Hn) / Hn_der ) * cos(n*theta);
      error = abs( abs( result ) - abs( prev_result ) );
      ++n;
    }
    
    //result *= phi_incident_mag;
    
    //const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
    //const complex< double > total_field = inc_field + result;
    //ofs << theta << "\t" << abs( result ) << "\t" << abs( inc_field ) << "\t" << abs( total_field ) <<  "\t" << R << endl; //write the file
    
    rEsult.resize(2);
    rEsult[REAL].resize(1);
    (rEsult[REAL])[0] = std::real(result);
    rEsult[IMAG].resize(1);
    (rEsult[IMAG])[0] = std::imag(result);
  
    return rEsult;
    
  }
  
};

