static char help[] ="Solves the simple ODE dy/dt=-10y, y(0)=1.";
/*
 Concepts: solving ordinary differential equations
 Processors: 1
 */
/* ------------------------------------------------------------------------
 This code demonstrates how one may solve an ordinary differential
 equation using the built-in ODE solvers.
 --------------------------------------------------------------------- */
#include "petscts.h"
#include <iostream>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;
using namespace boost::numeric;
/*
 Create an application context to contain data needed by the
 application-provided call-back routines, FormJacobian() and
 FormFunction().
 */
typedef struct {
  double G0;
  double beta;
  double Tg;

} AppCtx;

/*
 User-defined routines
 */
extern PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);
extern PetscErrorCode TSI_FormFunction (TS,PetscReal,Vec,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  TS             ts; /* timestepping context */
  Vec            y;  /* solution vectors */
  Mat            A;  /* Jacobian matrix */
  PetscErrorCode ierr;
  PetscReal      end_time, dt;
  PetscReal  ftime;
  PetscInt   num_time_steps;
  AppCtx         appctx;  /* user-defined work context */
  PetscScalar   *initial_condition;

  PetscInitialize(&argc,&argv,PETSC_NULL,help);
  //Set up the timestep (can be an option from command line)

  appctx.G0=3.76;    // shear modulus in dry state
  appctx.beta=-0.001682; // model parameters
  appctx.Tg=126+273.15;  // in K (glass transition temprature)


  dt = 0.1;  //time step
  end_time=200; //simulation time duration
  ierr = PetscOptionsGetReal(PETSC_NULL,"-dt",&dt,PETSC_NULL);CHKERRQ(ierr);
  num_time_steps = round(end_time/dt);
  /*
   Create vector to hold the solution
   */
  ierr = VecCreateSeq(PETSC_COMM_SELF,3,&y);CHKERRQ(ierr);

  /* Create matrix to hold Jacobian.  */
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,3,3,3,0,&A);CHKERRQ(ierr);

  /*
   Create timestepper context
   */
  ierr = TSCreate(PETSC_COMM_SELF,&ts);CHKERRQ(ierr);
  //ierr = TSSetProblemType(ts,TS_LINEAR);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,Monitor,&appctx,NULL);CHKERRQ(ierr);

  /*
   Set initial condition
   */
  ierr = VecGetArray(y,&initial_condition);CHKERRQ(ierr);

  initial_condition[0]=3.76;
  initial_condition[1]=60.0+273.15;  //temprature K
  initial_condition[2]=0.5;

  ierr = VecRestoreArray(y,&initial_condition);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,y);CHKERRQ(ierr);
  /*
   calculate resediual F(u,u',t)=0 in TSI_FormFunction
  */

  ierr = TSSetIFunction(ts,NULL,TSI_FormFunction,&appctx); CHKERRQ(ierr);

  /*
   Set the Jacobian matrix and the function used to compute Jacobians.
   */
  ierr = TSSetIJacobian(ts,A,A,FormJacobian,&appctx); CHKERRQ(ierr);

  /*
   This indicates that we are using Backward Euler’s method (implicit method).
   */
  ierr = TSSetType(ts,	TSBEULER);CHKERRQ(ierr);
  /*
   Set the initial time and the initial timestep given above.
//   */
  ierr = TSSetInitialTimeStep(ts,0.0,dt);CHKERRQ(ierr);
  /*
   Set a maximum number of timesteps and final simulation time.
   */
  ierr = TSSetDuration(ts,num_time_steps,end_time);
  /*
   Set any additional options from the options database. This
   includes all options for the nonlinear and linear solvers used
   internally the the timestepping routines.
   */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetUp(ts);CHKERRQ(ierr);

  /*
   Perform the solve. This is where the timestepping takes place.
   */
  ierr = TSSolve(ts, y);CHKERRQ(ierr);
  /*
   View information about the time-stepping method and the solution
   at the end time.
   */
//  TSView(ts, PETSC_VIEWER_STDOUT_SELF);
//  VecView(y,  PETSC_VIEWER_STDOUT_SELF);
  printf("\nThis is ftime: %f\n", ftime);
  /*
   Free the data structures constructed above
   */
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}



#undef __FUNCT__
#define __FUNCT__ "TSI_FormFunction"
PetscErrorCode TSI_FormFunction(TS ts,PetscReal t,Vec X,Vec X_t,Vec F,void *ctx)
{
  AppCtx   *appctx = (AppCtx*) ctx;   /* user-defined application context */
  PetscErrorCode ierr;
  PetscScalar    *y, *f, *y_t;
  double G;

  double G0=appctx->G0;
  double beta=appctx->beta;
  double Tg=appctx->Tg;

  ierr = VecGetArray(X,&y);CHKERRQ(ierr);
  ierr = VecGetArray(X_t,&y_t);CHKERRQ(ierr);
  ierr = VecGetArray(F,&f);CHKERRQ(ierr);

  double T=y[1];
  double c=y[2];

//  cout<<" G = "<< y[0] <<endl;
//  cout<<" T = "<< T <<endl;
//  cout<<" c = "<< c <<endl;
//  cout<<" beta = "<< beta <<endl;
//  cout<<" Tg = "<< Tg <<endl;

  G=G0*exp(-c*beta*log(1.0-T/Tg)*t);

  double dG_dT=(-c*beta*t*G)/(Tg*(T/Tg-1));
  double dG_dc=-beta*t*log(1-T/Tg)*G;
  double dG_dt=-c*beta*log(1-T/Tg)*G;

//  f[0]=y_t[0]-dG_dT*y_t[1]-dG_dc*y_t[2]-dG_dt;  //orignal
  f[0]=y_t[0]-dG_dt;  // after draping -dG_dT*y_t[1]-dG_dc*y_t[2] part
  f[1]=y_t[1]-20.0*M_PI*cos(2.0*M_PI*t);  //here M_PI is phi defined in math.h
  f[2]=y_t[2]+0.6*M_PI*cos(2.0*M_PI*t);

//  cout<<" t = "<< t <<endl;
//  cout<<" Gdot = "<< y_t[0] <<endl;
//  cout<<" Tdot = "<< y_t[1] <<endl;
//  cout<<" cdot = "<< y_t[2] <<endl;

//  cout<<" f[0] = "<< f[0] <<endl;
//  cout<<" f[1] = "<< f[1] <<endl;
//  cout<<" f[2] = "<< f[2] <<endl;

  ierr = VecRestoreArray(X,&y);CHKERRQ(ierr);
  ierr = VecRestoreArray(X_t,&y_t);CHKERRQ(ierr);
  ierr = VecRestoreArray(F,&f);CHKERRQ(ierr);
  return 0;
}

/* --------------------  Evaluate Jacobian F'(x) -------------------- */

#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
/*
 Calculate the Jacobian matrix J(X,t).

 Note: We put the Jacobian in the preconditioner storage B instead of J. This
 way we can give the -snes_mf_operator flag to check our work. This replaces
 J with a finite difference approximation, using our analytic Jacobian B for
 the preconditioner.
 */
PetscErrorCode FormJacobian(TS ts,PetscReal t,Vec X, Vec X_t, PetscReal a, Mat A, Mat B,void *ctx)
{
  AppCtx   *appctx = (AppCtx*) ctx;   /* user-defined application context */
  PetscErrorCode ierr;
  PetscScalar    *y, *y_t;
  double G;

  double G0=appctx->G0;
  double beta=appctx->beta;
  double Tg=appctx->Tg;

  ierr = VecGetArray(X,&y);CHKERRQ(ierr);
  ierr = VecGetArray(X_t,&y_t);CHKERRQ(ierr);

  double T=y[1];
  double c=y[2];

  double Tdot=y_t[1];
  double cdot=y_t[2];

  ublas::matrix<double> dF_du, dF_dudot, Amat;

  dF_du.resize(3,3);     dF_du.clear();
  dF_dudot.resize(3,3);  dF_dudot.clear();

  ublas::vector<int> row, col;
  row.resize(3);  col.resize(3);
  row[0]=0; row[1]=1; row[2]=2;
  col[0]=0; col[1]=1; col[2]=2;

//  cout<<"row = "<<row[0]<<" "<<row[1]<<" "<<row[2]<<endl;
//  cout<<"col = "<<col[0]<<" "<<col[1]<<" "<<col[2]<<endl;

  G=G0*exp(-c*beta*log(1-T/Tg)*t);

//  cout<<" G = "<< G <<endl;
//  cout<<" t = "<< t <<endl;
//  cout<<" T = "<< T <<endl;
//  cout<<" c = "<< c <<endl;

  dF_du(0,0)=(c*beta*t*Tdot)/(Tg*(T/Tg-1)) + beta*t*log(1-T/Tg)*cdot +c*beta*log(1-T/Tg);
  dF_du(0,1)=( (c*c*beta*t*t*Tdot)/(Tg*(T/Tg-1)) - beta*t*t*c*cdot*log(1-T/Tg)-c*c*beta*t*log(1-T/Tg)
              -c*t*Tdot/(Tg*(T/Tg-1))+t*cdot+c) * beta*G/(Tg*(T/Tg-1));
  dF_du(0,2)=( (-c*beta*t*t*Tdot)/(Tg*(T/Tg-1)) - beta*t*t*cdot*log(1-T/Tg)-c*beta*t*log(1-T/Tg)
              +(t*Tdot)/(Tg*(T/Tg-1)*log(1-T/Tg))+1  )*beta*log(1-T/Tg)*G;

//  cout<<" dF_du(0,0) = "<< dF_du(0,0) <<endl;
//  cout<<" dF_du(0,1) = "<< dF_du(0,1) <<endl;
//  cout<<" dF_du(0,2) = "<< dF_du(0,2) <<endl;

  dF_dudot(0,0)=1.0;  dF_dudot(1,1)=1.0;   dF_dudot(2,2)=1.0;
  Amat=dF_du + a*dF_dudot;

  ierr = MatSetValues(A,3,&(row)[0],3,&(col)[0],&Amat(0,0),INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}




/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u,void *ctx)
{
  PetscErrorCode ierr;
  PetscScalar  *y;

  ierr = VecGetArray(u,&y);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"time = %-6g Sol = %-11g\n",(double)time,(double)y[0]);  CHKERRQ(ierr);

  return 0;
}






















///*=================================================================================================
// This is only for constrant tempratrure
//=================================================================================================*/
//
//static char help[] ="Solves the simple ODE dy/dt=-10y, y(0)=1.";
///*
// Concepts: solving ordinary differential equations
// Processors: 1
// */
///* ------------------------------------------------------------------------
// This code demonstrates how one may solve an ordinary differential
// equation using the built-in ODE solvers.
// --------------------------------------------------------------------- */
//#include "petscts.h"
//#include <iostream>
//
//#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/numeric/ublas/vector.hpp>
//
//using namespace std;
//using namespace boost::numeric;
///*
// Create an application context to contain data needed by the
// application-provided call-back routines, FormJacobian() and
// FormFunction().
// */
//typedef struct {
//  double G0;
//  double beta;
//  double Tg;
//  double T;
//  double c;
//} AppCtx;
//
///*
// User-defined routines
// */
//extern PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);
//extern PetscErrorCode TSI_FormFunction (TS,PetscReal,Vec,Vec,Vec,void*);
//extern PetscErrorCode FormJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
//
//#undef __FUNCT__
//#define __FUNCT__ "main"
//int main(int argc,char **argv)
//{
//  TS             ts; /* timestepping context */
//  Vec            y;  /* solution vectors */
//  Mat            A;  /* Jacobian matrix */
//  PetscErrorCode ierr;
//  PetscReal      end_time, dt;
//  PetscReal  ftime;
//  PetscInt   num_time_steps;
//  AppCtx         appctx;  /* user-defined work context */
//  PetscScalar   *initial_condition;
//
//  PetscInitialize(&argc,&argv,PETSC_NULL,help);
//  //Set up the timestep (can be an option from command line)
//
//  appctx.G0=3.76;    // shear modulus in dry state
//  appctx.beta=-0.001682; // model parameters
//  appctx.Tg=126+273.15;  // in K (glass transition temprature)
//  appctx.T=80+273.15;  // in K (glass transition temprature)
//  appctx.c=1;  // in K (glass transition temprature)
//
//
//  dt = 0.1;  //time step
//  end_time=200; //simulation time duration
//  ierr = PetscOptionsGetReal(PETSC_NULL,"-dt",&dt,PETSC_NULL);CHKERRQ(ierr);
//  num_time_steps = round(end_time/dt);
//  /*
//   Create vector to hold the solution
//   */
//  ierr = VecCreateSeq(PETSC_COMM_SELF,1,&y);CHKERRQ(ierr);
//
//  /* Create matrix to hold Jacobian.  */
//  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,1,1,1,0,&A);CHKERRQ(ierr);
//
//  /*
//   Create timestepper context
//   */
//  ierr = TSCreate(PETSC_COMM_SELF,&ts);CHKERRQ(ierr);
//  //ierr = TSSetProblemType(ts,TS_LINEAR);CHKERRQ(ierr);
//  ierr = TSMonitorSet(ts,Monitor,&appctx,NULL);CHKERRQ(ierr);
//
//  /*
//   Set initial condition
//   */
//  ierr = VecGetArray(y,&initial_condition);CHKERRQ(ierr);
//
//  initial_condition[0]=3.76;
// 
//  ierr = VecRestoreArray(y,&initial_condition);CHKERRQ(ierr);
//  ierr = TSSetSolution(ts,y);CHKERRQ(ierr);
//  /*
//   calculate resediual F(u,u',t)=0 in TSI_FormFunction
//  */
//
//  ierr = TSSetIFunction(ts,NULL,TSI_FormFunction,&appctx); CHKERRQ(ierr);
//
//  /*
//   Set the Jacobian matrix and the function used to compute Jacobians.
//   */
//  ierr = TSSetIJacobian(ts,A,A,FormJacobian,&appctx); CHKERRQ(ierr);
//
//  /*
//   This indicates that we are using Backward Euler’s method (implicit method).
//   */
//  ierr = TSSetType(ts,	TSBEULER);CHKERRQ(ierr);
//  /*
//   Set the initial time and the initial timestep given above.
////   */
//  ierr = TSSetInitialTimeStep(ts,0.0,dt);CHKERRQ(ierr);
//  /*
//   Set a maximum number of timesteps and final simulation time.
//   */
//  ierr = TSSetDuration(ts,num_time_steps,end_time);
//  /*
//   Set any additional options from the options database. This
//   includes all options for the nonlinear and linear solvers used
//   internally the the timestepping routines.
//   */
//  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
//  ierr = TSSetUp(ts);CHKERRQ(ierr);
//
//  /*
//   Perform the solve. This is where the timestepping takes place.
//   */
//  ierr = TSSolve(ts, y);CHKERRQ(ierr);
//  /*
//   View information about the time-stepping method and the solution
//   at the end time.
//   */
////  TSView(ts, PETSC_VIEWER_STDOUT_SELF);
////  VecView(y,  PETSC_VIEWER_STDOUT_SELF);
//  printf("\nThis is ftime: %f\n", ftime);
//  /*
//   Free the data structures constructed above
//   */
//  ierr = VecDestroy(&y);CHKERRQ(ierr);
//  ierr = TSDestroy(&ts);CHKERRQ(ierr);
//  ierr = PetscFinalize();CHKERRQ(ierr);
//  return 0;
//}
//
//
//
//#undef __FUNCT__
//#define __FUNCT__ "TSI_FormFunction"
//PetscErrorCode TSI_FormFunction(TS ts,PetscReal t,Vec X,Vec X_t,Vec F,void *ctx)
//{
//  AppCtx   *appctx = (AppCtx*) ctx;   /* user-defined application context */
//  PetscErrorCode ierr;
//  PetscScalar    *y, *f, *y_t;
//  double G;
//
//  double G0=appctx->G0;
//  double beta=appctx->beta;
//  double Tg=appctx->Tg;
//  double T=appctx->T;
//  double c=appctx->c;
//
//  ierr = VecGetArray(X,&y);CHKERRQ(ierr);
//  ierr = VecGetArray(X_t,&y_t);CHKERRQ(ierr);
//  ierr = VecGetArray(F,&f);CHKERRQ(ierr);
//
////  cout<<" G = "<< y[0] <<endl;
////  cout<<" T = "<< T <<endl;
////  cout<<" c = "<< c <<endl;
////  cout<<" beta = "<< beta <<endl;
////  cout<<" Tg = "<< Tg <<endl;
//
//  G=G0*exp(-c*beta*log(1.0-T/Tg)*t);
//  double dG_dt=-c*beta*log(1-T/Tg)*G;
//
//  f[0]=y_t[0]-dG_dt;
// 
////  cout<<" t = "<< t <<endl;
////  cout<<" Gdot = "<< y_t[0] <<endl;
////  cout<<" Tdot = "<< y_t[1] <<endl;
////  cout<<" cdot = "<< y_t[2] <<endl;
////
////  cout<<" f[0] = "<< f[0] <<endl;
////  cout<<" f[1] = "<< f[1] <<endl;
////  cout<<" f[2] = "<< f[2] <<endl;
//
//  ierr = VecRestoreArray(X,&y);CHKERRQ(ierr);
//  ierr = VecRestoreArray(X_t,&y_t);CHKERRQ(ierr);
//  ierr = VecRestoreArray(F,&f);CHKERRQ(ierr);
//  return 0;
//}
//
///* --------------------  Evaluate Jacobian F'(x) -------------------- */
//
//#undef __FUNCT__
//#define __FUNCT__ "FormJacobian"
///*
// Calculate the Jacobian matrix J(X,t).
//
// Note: We put the Jacobian in the preconditioner storage B instead of J. This
// way we can give the -snes_mf_operator flag to check our work. This replaces
// J with a finite difference approximation, using our analytic Jacobian B for
// the preconditioner.
// */
//PetscErrorCode FormJacobian(TS ts,PetscReal t,Vec X, Vec X_t, PetscReal a, Mat A, Mat B,void *ctx)
//{
//  AppCtx   *appctx = (AppCtx*) ctx;   /* user-defined application context */
//  PetscErrorCode ierr;
//  PetscScalar    *y, *y_t;
//  double G;
//
//  double G0=appctx->G0;
//  double beta=appctx->beta;
//  double Tg=appctx->Tg;
//  double T=appctx->T;
//  double c=appctx->c;
//
//  ierr = VecGetArray(X,&y);CHKERRQ(ierr);
//  ierr = VecGetArray(X_t,&y_t);CHKERRQ(ierr);
//
//  double Tdot=y_t[1];
//  double cdot=y_t[2];
//
//  ublas::matrix<double> dF_du, dF_dudot, Amat;
//
//  dF_du.resize(1,1);     dF_du.clear();
//  dF_dudot.resize(1,1);  dF_dudot.clear();
//
//  ublas::vector<int> row, col;
//  row.resize(3);  col.resize(3);
//  row[0]=0;
//  col[0]=0;
//
////  cout<<"row = "<<row[0]<<" "<<row[1]<<" "<<row[2]<<endl;
////  cout<<"col = "<<col[0]<<" "<<col[1]<<" "<<col[2]<<endl;
//
//  G=G0*exp(-c*beta*log(1-T/Tg)*t);
//
////  cout<<" G = "<< G <<endl;
////  cout<<" t = "<< t <<endl;
////  cout<<" T = "<< T <<endl;
////  cout<<" c = "<< c <<endl;
//
//  dF_du(0,0)=c*beta*log(1-T/Tg);
//
////  cout<<" dF_du(0,0) = "<< dF_du(0,0) <<endl;
////  cout<<" dF_du(0,1) = "<< dF_du(0,1) <<endl;
////  cout<<" dF_du(0,2) = "<< dF_du(0,2) <<endl;
//
//  dF_dudot(0,0)=1.0;
//  Amat=dF_du + a*dF_dudot;
//
//  ierr = MatSetValues(A,1,&(row)[0],1,&(col)[0],&Amat(0,0),INSERT_VALUES);CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//
//  return 0;
//}
//
//
//
//
///* --------------------------------------------------------------------- */
//#undef __FUNCT__
//#define __FUNCT__ "Monitor"
//
//PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u,void *ctx)
//{
//  PetscErrorCode ierr;
//  PetscScalar  *y;
//
//  ierr = VecGetArray(u,&y);CHKERRQ(ierr);
//  ierr = PetscPrintf(PETSC_COMM_WORLD,"time = %-6g Sol = %-11g\n",(double)time,(double)y[0]);  CHKERRQ(ierr);
//
//  return 0;
//}

































