//
//  main.cpp
//  Reliability_Analysis
//
//  Created by Xiaoyi Zhou on 22/05/2015.
//  Copyright (c) 2015 Xiaoyi Zhou. All rights reserved.
//
// This code adopts the open-source Matlab toolbox - FERUM 4.1 (Finite Element
//   Reliability Using Matlab) devloped by Jean-Marc BOURINET at the IFMA (Institut
//   Français de Mécanique Avancée) in Clermont-Ferrand, France, which is a new
//   version of FERUM originally developed and maintained by Terje Haukaas,
//   Armen Der Kiureghian and other contributors at the University of California
//   at Berkeley initiated in 1999.
//
// FERUM 4.1 is available at http://www.ifma.fr/Recherche/laboratoires_recherche/FERUM
// FERUM 3.0 is available at http://www.ce.berkeley.edu/projects/ferum/
//
// References:
//   [1] Merchers, R. (1999) Structural reliability analysis and prediction,
//       2nd edition, Wiley.
//   [2] Ditlevsen, O. (2005) Structural reliability methods, electronic version
//   [3] Choi, et al. (2007) Reliability-based structural design, Springer
//   [4] Bourinet, J.-M. (2009) FERUM 4.0 User's Guide
//

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//              PROCEDURE OF SECOND-ORDER RELIABILITY METHOD                  //
//                                                                            //
//  1. Read FORM results including:                                           //
//       design point in x and u spaces,                                      //
//       values of LSF and its first-order derivative;                        //
//       cosine direction or alpha                                            //
//  2. Compute the second-order derivative of LSF                             //
//  3. Construct an orthogonal matrix P based on alpha                        //
//  4. Compute the matrix A = PHP'/||G'|| and form A11 by deleting the last   //
//       row and column of A                                                  //
//  5. Compute the eigenvalues of A11 to get curvatures, k_i, i = 1, ..., n-1 //
//  6. Calculate SORM approximation of faiure probability by empirical formula//
//       Breitung                                                             //
//       Tvedt                                                                //
//       Curve fitting [Der Kiureghian and De Stefano]                        //
//       Point fitting [Der Kiureghian et al.]                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

struct SORM{
  
  //------------------------------------------------------------------------------
  // d2x/dudu = - u*(dx/du) - (df(x)/dx)/f*(dx/du)^2
  // dx/du = f/phi
  //
  
  virtual PetscErrorCode d2x_dudu(ublas::vector<double> x,
                                  ublas::vector<double> u,
                                  Stochastic_Model probdata,
                                  ublas::vector<double> &dxdz,
                                  ublas::matrix<double> &ddxddu) {
    PetscFunctionBegin;
    
    int nvars;
    nvars = probdata.num_vars;
    
    double imean, istd;
    double iloc, ishape, iscale;
    double ilambda;
    double ilower, iupper;
    
    int dist_type;
    
    ddxddu.resize(nvars,nvars); ddxddu.clear();
    dxdz.resize(nvars); dxdz.clear();
    
    ublas::vector<double> z; // y = Az
    z = prod(probdata.Lo, u);
    
    using boost::math::normal_distribution;
    normal_distribution<> snorm(0,1);
    
    for (int i=0; i<nvars; i++) {
      dist_type = probdata.marg(i,0);
      switch (dist_type) {
        case 1: { // Normal distribution
          imean = probdata.marg(i,4);
          istd  = probdata.marg(i,5);
          normal_distribution<> mynorm(imean,istd);
          dxdz(i) = pdf(snorm,z(i))/pdf(mynorm,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) + pow(dxdz(i),2)*(x(i)-imean)/pow(istd,2);
          break;
        }
        case 2: { // Lognormal distribution
          using boost::math::lognormal_distribution;
          iloc   = probdata.marg(i,4);
          iscale = probdata.marg(i,5);
          lognormal_distribution<> my_logn(iloc,iscale);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_logn,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) + pow(dxdz(i),2)*(pow(iscale,2) + log(x(i))-iloc)/pow(iscale,2)/x(i);
          break;
        }
        case 3: { // Exponential distribution
          using boost::math::exponential_distribution;
          ilambda = 1/probdata.marg(i,1);
          exponential_distribution<> my_exp(ilambda);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_exp,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) + ilambda*pow(dxdz(i),2);
          break;
        }
        case 4: { // Gumbel distribution
          using boost::math::extreme_value_distribution;
          iloc   = probdata.marg(i,4);
          iscale = probdata.marg(i,5);
          extreme_value_distribution<> my_EVI(iloc,iscale);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_EVI,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) - (exp(-(x(i)-iloc)/iscale) - 1)/iscale*pow(dxdz(i),2);
          break;
        }
        case 5: { // Weibull distribution
          using boost::math::weibull_distribution;
          ishape = probdata.marg(i,4);
          iscale = probdata.marg(i,5);
          weibull_distribution<> my_wbl(ishape,iscale);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_wbl,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) -((ishape -1)/x(i) - ishape*pow(x(i)/iscale,ishape-1)/iscale)*pow(dxdz(i),2);
          break;
        }
        case 6: { // Uniform distribution
          using boost::math::uniform_distribution;
          ilower = probdata.marg(i,4);
          iupper = probdata.marg(i,5);
          uniform_distribution<> my_unif(ilower,iupper);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_unif,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i);
          break;
        }
        case 7: { // Gamma distribution
          using boost::math::gamma_distribution;
          ishape = probdata.marg(i,4);
          iscale = probdata.marg(i,5);
          gamma_distribution<> my_gamma(ishape,iscale);
          dxdz(i) = pdf(snorm,z(i))/pdf(my_gamma,x(i));
          ddxddu(i,i) = - z(i)*dxdz(i) - (1/x(i) - 1/iscale)*pow(dxdz(i),2);
          break;
        }
        default:
          break;
      }
    }
    PetscFunctionReturn(0);
  }
  
  //----------------------------------------------------------------------------
  // Conduct Gram-Schmidt QR factorization of a matrix
  //// 1. classic
  //void gschmidt(ublas::matrix<double> A, ublas::matrix<double> &Q) {
  //
  //}
  
  // 2. modified
  //
  
  virtual PetscErrorCode gramschmidt(ublas::matrix<double> A,
                                     ublas::matrix<double> &Q) {
    
    PetscFunctionBegin;
    
    int m, n;
    m = A.size1(); // the number of rows
    n = A.size2(); // the number of columns
    
    //ublas::zero_matrix<double> R(m, n);
    ublas::matrix<double> R(m, n);
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
        R(i,j) = 0;
      }
    }
    
    Q = A;
    
    cout<<"Size row: "<<m<<"\t Size column: "<<n<<endl;
    
    ublas::vector<double> irow(n);
    ublas::matrix<double> TQ(n,m);
    for (int i=0; i<m; i++) {
      TQ = trans(Q);
      // Get i-th column
      ublas::matrix_row<ublas::matrix<double> > irow_Q(TQ, i);
      // Calculate normal of i-th column
      R(i, i) = norm_2(irow_Q);
      // Normalized the i-th column
      irow_Q = irow_Q/R(i,i);
      for (int j=0; j<n; j++) {
        Q(j, i) = irow_Q(j);
      }
      // --
      for (int k = i+1; k<n; k++) {
        ublas::matrix_row<ublas::matrix<double> > krow_Q(TQ, k);
        R(i,k) = inner_prod(irow_Q,krow_Q);
        krow_Q = krow_Q - R(i,k)*irow_Q;
        for (int j=0; j<n; j++) {
          Q(j,k) = krow_Q(j);
        }
      }
    }
    
    PetscFunctionReturn(0);
  }
  
  //----------------------------------------------------------------------------
  // Flip matrix left to right
  //
  
  virtual PetscErrorCode fliplr(ublas::matrix<double> &A) {
    
    PetscFunctionBegin;
    
    int nrow, ncol;
    nrow = A.size1(); ncol = A.size2();
    ublas::matrix<double> B;
    B.resize(nrow,ncol); B.clear();
    for (int i=0; i<nrow; i++) {
      for (int j=0; j<ncol; j++) {
        B(i,j) = A(i,ncol-j-1);
      }
    }
    A.clear();
    A = B;
    
    PetscFunctionReturn(0);
  }
  
  
  //------------------------------------------------------------------------------
  // Construct an orthogonal matrix to rotate vector
  //
  
  virtual PetscErrorCode orthonormal_matrix(ublas::vector<double> alpha,
                                            ublas::matrix<double> &Q) {
    
    PetscFunctionBegin;
    
    double nvars;
    
    nvars = alpha.size();
    ublas::identity_matrix<double> A1(nvars);
    ublas::matrix<double> A(nvars,nvars);
    A = A1; fliplr(A);
    for (int i=0; i<nvars; i++) {
      A(i,0) = alpha(i);
    }
    
    // conduct Gram-Schmidt orthonormalization
    gramschmidt(A,Q);
    fliplr(Q);
    Q = trans(Q);
    cout<<"\n"<<"Orthogonal matrix: "<<Q<<endl;
    
    PetscFunctionReturn(0);
  }
  
  //------------------------------------------------------------------------------
  // Construct Hessian matrix at U-space or standard normal space
  //
  
  virtual PetscErrorCode Hessian_Matrix(ublas::vector<double> dxdu,
                                        ublas::vector<double> grad_g,
                                        ublas::matrix<double> Hess_g,
                                        ublas::matrix<double> &Hess_G,
                                        ublas::matrix<double> Hess_x) {
    
    PetscFunctionBegin;
    
    //
    // ddG/dudu = ddg/dxdx*(dx/du)^2 + dg/dx*(ddx/dudu)
    //
    int nvars;
    nvars = dxdu.size();
    ublas::matrix<double> mat_grad_g;
    ublas::matrix<double> mat_dxdu;
    
    mat_grad_g.resize(nvars,nvars); mat_grad_g.clear();
    mat_dxdu.resize(nvars,nvars); mat_dxdu.clear();
    
    
    for (int i=0; i<nvars; i++) {
      mat_grad_g(i,i) = grad_g(i);
      mat_dxdu(i,i) = dxdu(i);
    }
    
    Hess_G.resize(nvars,nvars); Hess_G.clear();
    ublas::matrix<double> temp_Hess_G;
    temp_Hess_G = prod(Hess_g,mat_dxdu);cout<<"\nHess x: "<<Hess_x<<endl;
    Hess_G = prod(mat_dxdu,temp_Hess_G) + prod(mat_grad_g, Hess_x);
    
    PetscFunctionReturn(0);
  }
  
};