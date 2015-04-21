/** \file AnalyticalSolutions.hpp
  \ingroup mofem_helmholtz_elem

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

/** \brief Generic structure for analytical function
  \ingroup mofem_helmholtz_elem
*/
struct GenericAnalyticalSolution {

  enum VALUE_TYPE { REAL = 0, IMAG, LAST_VAL_TYPE };

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) = 0;

};

/** List of analytical solution
  \ingroup mofem_helmholtz_elem
*/
enum AnalyticalSolutionTypes {
  SOFT_SPHERE_SCATTER_WAVE, 
  PLANE_WAVE,
  CYLINDER_SCATTER_WAVE,
  INCIDENT_WAVE
};

/** Line command list of analytical solutions
  \ingroup mofem_helmholtz_elem
*/
const char *analytical_solution_types[] = { 
  "soft_sphere_incident_wave", 
  "plane_wave",
  "cylinder_scatter_wave"
  "incident_wave"
};

/** Incident wave 
  \ingroup mofem_helmholtz_elem


  \f[
  p_\textrm{inc} = \exp(ikz)
  \f]

  */
struct IncidentWave: public GenericAnalyticalSolution {

  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double pOwer;
 
  IncidentWave(double wavenumber,double power):
    wAvenumber(wavenumber),pOwer(power) {}
  ~IncidentWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const complex< double > i( 0.0, 1.0 );
    complex< double > result = 0.0;

    result = pOwer*exp(i*wAvenumber*z);

    rEsult.resize(2);
    rEsult[REAL].resize(1);
    (rEsult[REAL])[0] = std::real(result);
    rEsult[IMAG].resize(1);
    (rEsult[IMAG])[0] = std::imag(result);

    return rEsult;

  }	

};

/** Calculate the analytical solution of impinging wave on sphere
  \ingroup mofem_helmholtz_elem

  Equations from:
  <http://ansol.us/Products/Coustyx/Validation/MultiDomain/Scattering/PlaneWave/SoftSphere/Downloads/dataset_description.pdf>

  \f[
  p_\textrm{scattered} = \sum_0^N A_l h_l(kr)P_l(\cos(\phi))
  \f]

  where \f$h_l\f$ is the Hankel function of the first kind, \f$\phi\f$ is polar
  angle and \f$A_l\f$ is a constant. Constant is  should be caculated such that
  it satisfies both the Helmholtz wave equation and the Sommerfeld radiation
  condition.

  \f[
  A_l = -(2l+1)i^l \frac{j_l(ka)}{h_l(ka)}
  \f]
  where a is scatter sphere radius and \f$j_l\f$ Bessel function.

  */
struct SoftSphereScatterWave: public GenericAnalyticalSolution {
  
  vector<complex<double> > vecAl; ///< this is to calculate constant values of series only once
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double sphereRadius;
   
  SoftSphereScatterWave(double wavenumber,double sphere_radius = 1.): 
    wAvenumber(wavenumber),sphereRadius(sphere_radius) {}
  virtual ~SoftSphereScatterWave() {}
   
  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double R = sqrt(x2+y2+z2); 
    double cos_theta = z/R;

    const double k = wAvenumber;  	//Wave number
    const double a = sphereRadius;      //radius of the sphere,wait to modify by user

    const complex< double > i( 0.0, 1.0 );
    complex< double > Al;

    complex< double > result = 0.0;
    complex< double > prev_result;

    double error = 100.0;
    unsigned int n = 0; //initialized the infinite series loop
    
    while( error > tol )  //finding the acoustic potential in one single point.
    {

      if(vecAl.size()>n) {
	Al = vecAl[n];
      } else {
	// spherical Bessel function
	double jn = sph_bessel(n,k*a);
	// spherical Hankel function
	complex<double> hn = sph_hankel_1(n,k*a);
	//Constant term
	Al = -(2.0*n+1)*pow(i,n)*jn/hn;
	vecAl.push_back(Al);
      }
      
      prev_result = result;

      // Legendre function
      double Pn = legendre_p(n,cos_theta);
      result += Al*sph_hankel_1(n,k*R)*Pn;

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
  \ingroup mofem_helmholtz_elem

  WARRING: not tested
  
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
  \ingroup mofem_helmholtz_elem

  WARRING: not tested

  Paper: 
    The generalized finite element method for Helmholtz equation: Theory, computation, and open problems
    Theofanis Strouboulis, Ivo Babuska, Realino Hidaja

*/
struct CylinderScatterWave: public GenericAnalyticalSolution {
  
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double shereRadius;
   
  CylinderScatterWave(double wavenumber): wAvenumber(wavenumber) {}
  virtual ~CylinderScatterWave() {}
   
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

template <typename FUNEVAL>
PetscErrorCode solve_problem(FieldInterface& m_field,
  const string& problem_name,const string& fe_name,
  const string& re_field,const string &im_field,
  InsertMode mode,
  FUNEVAL &fun_evaluator,PetscBool is_partitioned) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  

  Mat A;
  ierr = m_field.MatCreateMPIAIJWithArrays(problem_name,&A); CHKERRQ(ierr);
  Vec D;

  vector<Vec> vec_F;
  vec_F.resize(2);

  ierr = m_field.VecCreateGhost(problem_name,ROW,&vec_F[0]); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost(problem_name,ROW,&vec_F[1]); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost(problem_name,COL,&D); CHKERRQ(ierr);

  FieldApproximationH1 field_approximation(m_field);
  // This increase rule for numerical intergaration. In case of 10 node
  // elements jacobian is varing lineary across element, that way to element
  // rule is added 1.
  field_approximation.addToRule = 1; 

  ierr = field_approximation.loopMatrixAndVector(
    problem_name,fe_name,re_field,A,vec_F,fun_evaluator); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  for(int ss = 0;ss<GenericAnalyticalSolution::LAST_VAL_TYPE;ss++) {

    // solve problem
    ierr = KSPSolve(solver,vec_F[ss],D); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    // save data on mesh
    if(ss == GenericAnalyticalSolution::REAL) {

      if(is_partitioned) {
	ierr = m_field.set_global_ghost_vector(problem_name,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
      } else {
	ierr = m_field.set_local_ghost_vector(problem_name,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
      }

      VecZeroEntries(D);
      ierr = VecGhostUpdateBegin(D,mode,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,mode,SCATTER_FORWARD); CHKERRQ(ierr);

    } else {
      if(is_partitioned) {
	ierr = m_field.set_other_local_ghost_vector(problem_name,re_field,im_field,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
      } else {
	ierr = m_field.set_other_global_ghost_vector(problem_name,re_field,im_field,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
      }
    }

  }

  // clean 
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&vec_F[GenericAnalyticalSolution::REAL]); CHKERRQ(ierr);
  ierr = VecDestroy(&vec_F[GenericAnalyticalSolution::IMAG]); CHKERRQ(ierr);

  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


