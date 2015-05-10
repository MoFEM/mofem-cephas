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
  \bug point source example not implemented.
*/
struct GenericAnalyticalSolution {

  enum VALUE_TYPE { REAL = 0, IMAG, LAST_VAL_TYPE };

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) = 0;
  virtual ~GenericAnalyticalSolution() {}

};

/** List of analytical solution
  \ingroup mofem_helmholtz_elem
*/
enum AnalyticalSolutionTypes {
  HARD_SPHERE_SCATTER_WAVE,
  SOFT_SPHERE_SCATTER_WAVE,
  PLANE_WAVE,
  HARD_CYLINDER_SCATTER_WAVE,
  SOFT_CYLINDER_SCATTER_WAVE,
  INCIDENT_WAVE,
  NO_ANALYTICAL_SOLUTION
};

/** Line command list of analytical solutions
  \ingroup mofem_helmholtz_elem
*/
const char *analytical_solution_types[] = {
  "hard_sphere_incident_wave",
  "soft_sphere_incident_wave",
  "plane_wave",
  "hard_cylinder_scatter_wave",
  "soft_cylinder_scatter_wave",
  "incident_wave",
  "no_analytical_solution"
};

// **** Analytical solutions ****

/** Incident wave
  \ingroup mofem_helmholtz_elem


  Equation from:
  Ihlenburg,Finite element analysis of acoustic scattering Springer Science & Business Media.

  Some details can be found here:
  <http://ansol.us/Products/Coustyx/Validation/MultiDomain/Scattering/PlaneWave/HardSphere/Downloads/dataset_description.pdf>

  \f[
  p_\textrm{inc} = \exp(ikd \cdot \mathbf{x})
  \f]


  */
struct IncidentWave: public GenericAnalyticalSolution {

  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  ublas::vector<double> dIrection;
  ublas::vector<double> cOordinate;
  double pOwerReal; ///< The real amplitude of the incident wave
  double pOwerImag;

  IncidentWave(double wavenumber,ublas::vector<double> d,double r_power = 1,double i_power = 0):
    wAvenumber(wavenumber),dIrection(d),pOwerReal(r_power),pOwerImag(i_power) {}
  ~IncidentWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const complex< double > i( 0.0, 1.0 );
    complex< double > result = 0.0;
    cOordinate.resize(3);
    cOordinate[0] = x;
    cOordinate[1] = y;
    cOordinate[2] = z;

    result = (pOwerReal+i*pOwerImag)*exp(i*wAvenumber*inner_prod(dIrection,cOordinate));

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


  equation from:
  Ihlenburg,Finite element analysis of acoustic scattering Springer Science & Business Media.

  \f[
  p_\textrm{scattered} = \sum_0^N A_l h_l(kr)P_l(\cos(\phi))
  \f]

  where \f$h_l\f$ is the Hankel function of the first kind, \f$\phi\f$ is polar
  angle and \f$A_l\f$ is a constant. Constant is  should be calculated such that
  it satisfies both the Helmholtz wave equation and the Sommerfeld radiation
  condition.

  \f[
  A_l = -(2l+1)i^l \frac{j_{l}'(ka)}{h_{l}'(ka)}
  \f]
  where a is scatter sphere radius and \f$j_l\f$ Spherical Bessel function.

  */
struct HardSphereScatterWave: public GenericAnalyticalSolution {

  vector<complex<double> > vecAl; ///< this is to calculate constant values of series only once
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double sphereRadius;

  HardSphereScatterWave(double wavenumber,double sphere_radius):
    wAvenumber(wavenumber),sphereRadius(sphere_radius) {}
  virtual ~HardSphereScatterWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double R = sqrt(x2+y2+z2);
    double cos_theta = z/R;

    const double k = wAvenumber;    //Wave number
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
        double const1 = k*a;
        double jn_der = n / const1 * sph_bessel( n, const1 ) - sph_bessel( n + 1, const1 );

        // spherical Hankel function
        complex< double > hn_der = n / const1 * sph_hankel_1( n, const1 ) - sph_hankel_1( n + 1, const1 );
        //Constant term
        Al = -(2.0*n+1)*pow(i,n)*jn_der/hn_der;
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

  where a is scatter sphere radius and \f$j_l\f$ Spherical Bessel function.


  */
struct SoftSphereScatterWave: public GenericAnalyticalSolution {

  vector<complex<double> > vecAl; ///< this is to calculate constant values of series only once
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  double sphereRadius;


  SoftSphereScatterWave(double wavenumber,double sphere_radius):
    wAvenumber(wavenumber),sphereRadius(sphere_radius) {}
  virtual ~SoftSphereScatterWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;

    double x2 = x*x;
    double y2 = y*y;
    double z2 = z*z;
    double R = sqrt(x2+y2+z2);
    double cos_theta = z/R;

    const double k = wAvenumber;    //Wave number
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


  \f[
  p_\textrm{scattered} = exp^{ik\mathbf{x}\Theta}
  \f]

  where:

  \f[
  \mathbf{x} = [x,y]
  \f]

  \f[
  \Theta = k[\cos(\theta),\sin(\theta)]
  \f]

  \theta is the wave propagating direction from range [0,2\pi]


   Paper: Gillman, A., Djellouli, R., & Amara, M. (2007).
   A mixed hybrid formulation based on oscillated finite element polynomials for solving Helmholtz problems.
   Journal of computational and applied mathematics

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


    //const complex< double > inc_field = exp( i * k * R * cos( theta ) );  //???? Incident wave
    //const complex< double > total_field = inc_field + result;

    result = exp(i*(k*cos(tHeta)*x+k*sin(tHeta)*y));

    rEsult.resize(2);
    rEsult[REAL].resize(1);
    (rEsult[REAL])[0] = std::real(result);
    rEsult[IMAG].resize(1);
    (rEsult[IMAG])[0] = std::imag(result);

    return rEsult;

  }

};

/** \brief Calculate the analytical solution of impinging wave on cylinder
  \ingroup mofem_helmholtz_elem

   \f[
  p_\textrm{scattered} = \sum_0^N A_l H_l(kr)\cos(l\theta)
  \f]

  where \f$H_l\f$ is the cylindrical Hankel function of the first kind, \f$\theta\f$
  is the polar angle in polar coordinate and \f$A_l\f$ is a constant. Constant is
  should be caculated such that it satisfies both the Helmholtz wave equation and the
  Sommerfeld radiation condition.

  \f[
  A_l = -\epsilon_{l} i^l \frac{J_{l}'(ka)}{H_{l}'(ka)}
  \f]

  where a is scatter sphere radius and \f$J_l\f$ Cylindrical Bessel function.

  \f[
  \epsilon_{l} = 1 \textrm{when}l=0
  \f]

   \f[
  \epsilon_{l} = 2 \textrm{when}l \neq 0
  \f]

  Paper:
    Kechroud, R., Soulaimani, A., & Antoine, X. (2009).
    A performance study of plane wave finite element methods with a Padé-type artificial boundary condition in acoustic scattering.
    Advances in Engineering Software, 40(8), 738-750.

*/

struct HardCylinderScatterWave: public GenericAnalyticalSolution {

  vector<complex<double> > vecAl; ///< this is to calculate constant values of series only once
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  //double shereRadius;
  double a;

  HardCylinderScatterWave(double wavenumber,double sphere_radius = 0.5): wAvenumber(wavenumber),a(sphere_radius) {}
  virtual ~HardCylinderScatterWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;
    double x2 = x*x,y2 = y*y;
    double R = sqrt(x2+y2);
    double theta = atan2(y,x)+2*M_PI;
  //double cos_theta = z/R;

    const double k = wAvenumber;  //Wave number
    const double const1 = k * a;
    double const2 = k * R;

    const complex< double > i( 0.0, 1.0 );
    complex< double > Al;
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
      if(vecAl.size()>n) {
        Al = vecAl[n-1];
      } else {
        // cylindrical Bessel function
        double Jn_der_ka = n / const1 * cyl_bessel_j( n, const1 ) - cyl_bessel_j( n + 1, const1 );
        // cylindrical Hankel function
        complex<double> Hn_der_ka = n / const1 * cyl_hankel_1( n, const1 ) - cyl_hankel_1( n + 1, const1 );
        //Constant term
        Al = -2.0*pow(i,n)*Jn_der_ka/Hn_der_ka;
        vecAl.push_back(Al);
      }

      prev_result = result;
      complex< double >Hn_kr = cyl_hankel_1( n, const2 );  //S Hankel first kind function

      result += Al * Hn_kr * cos(n*theta);
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


/** \brief Calculate the analytical solution of impinging wave on cylinder
  \ingroup mofem_helmholtz_elem

   \f[
  p_\textrm{scattered} = \sum_0^N A_l H_l(kr)\cos(l\theta)
  \f]

  where \f$H_l\f$ is the cylindrical Hankel function of the first kind, \f$\theta\f$
  is the polar angle in polar coordinate and \f$A_l\f$ is a constant. Constant is
  should be caculated such that it satisfies both the Helmholtz wave equation and the
  Sommerfeld radiation condition.

  \f[
  A_l = -\epsilon_{l} i^l \frac{J_{l}(ka)}{H_{l}(ka)}
  \f]

  where a is scatter sphere radius and \f$J_l\f$ Cylindrical Bessel function.

  \f[
  \epsilon_{l} = 1 \textrm{when}l=0
  \f]

   \f[
  \epsilon_{l} = 2 \textrm{when}l \neq 0
  \f]

  Paper:
    Kechroud, R., Soulaimani, A., & Antoine, X. (2009).
    A performance study of plane wave finite element methods with a Padé-type artificial boundary condition in acoustic scattering.
    Advances in Engineering Software, 40(8), 738-750.


*/
struct SoftCylinderScatterWave: public GenericAnalyticalSolution {

  vector<complex<double> > vecAl; ///< this is to calculate constant values of series only once
  vector<ublas::vector<double> > rEsult;
  double wAvenumber;
  //double shereRadius;
  double a;

  SoftCylinderScatterWave(double wavenumber,double sphere_radius = 0.5): wAvenumber(wavenumber),a(sphere_radius) {}
  virtual ~SoftCylinderScatterWave() {}

  virtual vector<ublas::vector<double> >& operator()(double x, double y, double z) {

    const double tol = 1.0e-10;
    double x2 = x*x,y2 = y*y;
    double R = sqrt(x2+y2);
    double theta = atan2(y,x)+2*M_PI;

    const double k = wAvenumber;  //Wave number
    const double const1 = k * a;
    double const2 = k * R;

    const complex< double > i( 0.0, 1.0 );
    complex< double > Al;
    // magnitude of incident wave
    //const double phi_incident_mag = 1.0;

    complex< double > result = 0.0;
    complex< double > prev_result;

    double error = 100.0;
    unsigned int n = 1; //initialized the infinite series loop

    double Jn_zero = cyl_bessel_j( 0, const1 );
    complex< double > Hn_zero_kr = cyl_hankel_1( 0, const2 );  //S Hankel first kind function
    complex< double > Hn_zero_ka = cyl_hankel_1( 0, const1 );  //S Hankel first kind function
    //n=0;
    result -= (Jn_zero * Hn_zero_kr)/Hn_zero_ka;

    while( error > tol )  //finding the acoustic potential in one single point.
    {


      if(vecAl.size()>n) {
        Al = vecAl[n-1];
      } else {
        // cylindrical Bessel function
        double Jn_ka = cyl_bessel_j( n, const1 );
        // cylindrical Hankel function
        complex<double> Hn_ka = cyl_hankel_1( n, const1 );
        //Constant term
        Al = -2.0*pow(i,n)*Jn_ka/Hn_ka;
        vecAl.push_back(Al);
      }


      prev_result = result;

      complex< double >Hn_kr = cyl_hankel_1( n, const2 );  //S Hankel first kind function

      result += Al * Hn_kr * cos(n*theta);
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

// **** Surface boundary conditions ( used to enforce BC on surface ) ****

// ** Dirichlet BC **

// **** Function to solve best approximation problem ****

template <typename FUNVAL>
PetscErrorCode calculate_matrix_and_vector(
  FieldInterface& m_field,const string &problem_name,const string &fe_name,const string &re_field,
  Mat A,vector<Vec> vec_F,FUNVAL &fun_evaluator
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  FieldApproximationH1 field_approximation(m_field);
  // This increase rule for numerical integration. In case of 10 node
  // elements jacobian is varying linearly across element, that way to element
  // rule is added 1.
  field_approximation.addToRule = 1;

  ierr = field_approximation.loopMatrixAndVector(
    problem_name,fe_name,re_field,A,vec_F,fun_evaluator
  ); CHKERRQ(ierr);

  if(A) {
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode save_data_on_mesh(
  FieldInterface& m_field,const string &problem_name,const string& re_field,const string &im_field,
  Vec D,int vt,InsertMode mode,PetscBool is_partitioned
) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  // save data on mesh
  if(vt == GenericAnalyticalSolution::REAL) {
    /* set data to field from solution vec */
    if(is_partitioned) {
      ierr = m_field.set_local_ghost_vector(problem_name,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
    } else {
      ierr = m_field.set_global_ghost_vector(problem_name,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
    }

  } else {

    if(is_partitioned) {
      ierr = m_field.set_other_local_ghost_vector(problem_name,re_field,im_field,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
    } else {
      ierr = m_field.set_other_global_ghost_vector(problem_name,re_field,im_field,COL,D,mode,SCATTER_REVERSE); CHKERRQ(ierr);
    }

  }

  PetscFunctionReturn(0);
}

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
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  Vec D;

  vector<Vec> vec_F;
  vec_F.resize(2);

  ierr = m_field.VecCreateGhost(problem_name,ROW,&vec_F[0]); CHKERRQ(ierr); /* real */
  ierr = m_field.VecCreateGhost(problem_name,ROW,&vec_F[1]); CHKERRQ(ierr); /* imag */
  ierr = m_field.VecCreateGhost(problem_name,COL,&D); CHKERRQ(ierr);

  for(int ss = 0;ss<2;ss++) {
    ierr = VecZeroEntries(vec_F[ss]); CHKERRQ(ierr);
    ierr = VecGhostUpdateBegin(vec_F[ss],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(vec_F[ss],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  }

  ierr = calculate_matrix_and_vector(m_field,problem_name,fe_name,re_field,A,vec_F,fun_evaluator); CHKERRQ(ierr);

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

    ierr = save_data_on_mesh(m_field,problem_name,re_field,im_field,D,ss,mode,is_partitioned); CHKERRQ(ierr);
  }

  // clean
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&vec_F[GenericAnalyticalSolution::REAL]); CHKERRQ(ierr);
  ierr = VecDestroy(&vec_F[GenericAnalyticalSolution::IMAG]); CHKERRQ(ierr);

  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
