/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 *
 * Description: Implementation of thermal stress, i.e. right hand side as result of thermal stresses
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
 */

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

#ifndef __FE2_MACRO_SOLVER_HPP
#define __FE2_MACRO_SOLVER_HPP

#include "ElasticFE_RVELagrange_Periodic.hpp"
#include "ElasticFE_RVELagrange_RigidBodyTranslation.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"
#include "ElasticFE_RVELagrange_Periodic_Multi_Rhs.hpp"
#include "ElasticFE_RVELagrange_Homogenized_Stress_Periodic.hpp"

#include <FE2_ElasticFEMethod.hpp>
#include <FE2_Rhs_r_PSFEM.hpp>
#include <FE2_Rhs_rs_PSFEM.hpp>

namespace MoFEM {
  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                         MULTIPLE LAYERS LAMINATE                          *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct FE2_Macro_Solver_Laminate {
    
    ublas::vector<ublas::matrix<double> > RVE_Dmat_r;
    ublas::matrix<ublas::matrix<double> > RVE_Dmat_rs;
    
    ublas::vector<ublas::matrix<double> > Ply_1st_Dmat_r;
    ublas::vector<ublas::matrix<double> > Ply_2nd_Dmat_r;
    ublas::vector<ublas::matrix<double> > Ply_3rd_Dmat_r;
    ublas::vector<ublas::matrix<double> > Ply_4th_Dmat_r;
    
    ublas::matrix<ublas::matrix<double> > Ply_1st_Dmat_rs;
    ublas::matrix<ublas::matrix<double> > Ply_2nd_Dmat_rs;
    ublas::matrix<ublas::matrix<double> > Ply_3rd_Dmat_rs;
    ublas::matrix<ublas::matrix<double> > Ply_4th_Dmat_rs;
    
    //global variable Dmat
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r_Em;
    ublas::matrix<double> Dmat_r_NUm;
    ublas::matrix<double> Dmat_r_Ep;
    ublas::matrix<double> Dmat_r_Ez;
    ublas::matrix<double> Dmat_r_NUp;
    ublas::matrix<double> Dmat_r_NUpz;
    ublas::matrix<double> Dmat_r_Gzp;
    ublas::matrix<double> Dmat_r_Ef;
    ublas::matrix<double> Dmat_r_NUf;
    
    ublas::matrix<double> Dmat_rs_EmEm;
    ublas::matrix<double> Dmat_rs_NUmNUm;
    ublas::matrix<double> Dmat_rs_EpEp;
    ublas::matrix<double> Dmat_rs_EzEz;
    ublas::matrix<double> Dmat_rs_NUpNUp;
    ublas::matrix<double> Dmat_rs_NUpzNUpz;
    ublas::matrix<double> Dmat_rs_GzpGzp;
    ublas::matrix<double> Dmat_rs_EfEf;
    ublas::matrix<double> Dmat_rs_NUfNUf;
    
    ublas::matrix<double> Dmat_1st_Ply;
    ublas::matrix<double> Dmat_1st_Ply_r_Em;
    ublas::matrix<double> Dmat_1st_Ply_r_NUm;
    ublas::matrix<double> Dmat_1st_Ply_r_Ep;
    ublas::matrix<double> Dmat_1st_Ply_r_Ez;
    ublas::matrix<double> Dmat_1st_Ply_r_NUp;
    ublas::matrix<double> Dmat_1st_Ply_r_NUpz;
    ublas::matrix<double> Dmat_1st_Ply_r_Gzp;
    ublas::matrix<double> Dmat_1st_Ply_r_Ef;
    ublas::matrix<double> Dmat_1st_Ply_r_NUf;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_1st_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_2nd_Ply;
    ublas::matrix<double> Dmat_2nd_Ply_r_Em;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUm;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ep;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ez;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUp;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_2nd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_2nd_Ply_r_Ef;
    ublas::matrix<double> Dmat_2nd_Ply_r_NUf;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_2nd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_3rd_Ply;
    ublas::matrix<double> Dmat_3rd_Ply_r_Em;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUm;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ep;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ez;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUp;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUpz;
    ublas::matrix<double> Dmat_3rd_Ply_r_Gzp;
    ublas::matrix<double> Dmat_3rd_Ply_r_Ef;
    ublas::matrix<double> Dmat_3rd_Ply_r_NUf;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_3rd_Ply_r_Theta_4;
    
    ublas::matrix<double> Dmat_4th_Ply;
    ublas::matrix<double> Dmat_4th_Ply_r_Em;
    ublas::matrix<double> Dmat_4th_Ply_r_NUm;
    ublas::matrix<double> Dmat_4th_Ply_r_Ep;
    ublas::matrix<double> Dmat_4th_Ply_r_Ez;
    ublas::matrix<double> Dmat_4th_Ply_r_NUp;
    ublas::matrix<double> Dmat_4th_Ply_r_NUpz;
    ublas::matrix<double> Dmat_4th_Ply_r_Gzp;
    ublas::matrix<double> Dmat_4th_Ply_r_Ef;
    ublas::matrix<double> Dmat_4th_Ply_r_NUf;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_1;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_2;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_3;
    ublas::matrix<double> Dmat_4th_Ply_r_Theta_4;
    
    //------------------------------------------------------------------------------
    // To transform constitutive matrix
    //
    
    virtual PetscErrorCode Dmat_Transformation(double theta,
                                               ublas::matrix<FieldData> Dmat_123,
                                               ublas::matrix<FieldData> &Dmat_xyz) {
      PetscFunctionBegin;
      
      double l1, l2, l3, m1, m2, m3, n1, n2, n3;
      l1=cos(theta);
      m1=sin(theta);
      n1=0.0;
      l2=-sin(theta);
      m2=cos(theta);
      n2=0.0;
      l3=0.0;
      m3=0.0;
      n3=1.0;
      
      ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6); T_strain.clear();
      T_strain(0,0)=l1*l1;
      T_strain(0,1)=m1*m1;
      T_strain(0,2)=n1*n1;
      T_strain(0,3)=l1*m1;
      T_strain(0,4)=m1*n1;
      T_strain(0,5)=l1*n1;
      
      T_strain(1,0)=l2*l2;
      T_strain(1,1)=m2*m2;
      T_strain(1,2)=n2*n2;
      T_strain(1,3)=l2*m2;
      T_strain(1,4)=m2*n2;
      T_strain(1,5)=l2*n2;
      
      T_strain(2,0)=l3*l3;
      T_strain(2,1)=m3*m3;
      T_strain(2,2)=n3*n3;
      T_strain(2,3)=l3*m3;
      T_strain(2,4)=m3*n3;
      T_strain(2,5)=l3*n3;
      
      T_strain(3,0)=2*l1*l2;
      T_strain(3,1)=2*m1*m2;
      T_strain(3,2)=2*n1*n2;
      T_strain(3,3)=l1*m2+m1*l2;
      T_strain(3,4)=m1*n2+n1*m2;
      T_strain(3,5)=l1*n2+n1*l2;
      
      T_strain(4,0)=2*l2*l3;
      T_strain(4,1)=2*m2*m3;
      T_strain(4,2)=2*n2*n3;
      T_strain(4,3)=l2*m3+m2*l3;
      T_strain(4,4)=m2*n3+n2*m3;
      T_strain(4,5)=l2*n3+n2*l3;
      
      T_strain(5,0)=2*l1*l3;
      T_strain(5,1)=2*m1*m3;
      T_strain(5,2)=2*n1*n3;
      T_strain(5,3)=l1*m3+m1*l3;
      T_strain(5,4)=m1*n3+n1*m3;
      T_strain(5,5)=l1*n3+n1*l3;
      
      //cout<<"\n\nT_strain = "<<T_strain<<endl;
      
      ublas::matrix<FieldData> Mat1=prod(Dmat_123,T_strain);
      Dmat_xyz = prod(trans(T_strain), Mat1);
      
      
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode Dmat_Transformation_r_Theta(double theta,
                                                       ublas::matrix<FieldData> Dmat_123,
                                                       ublas::matrix<FieldData> &Dmat_xyz) {
      PetscFunctionBegin;
      
      double l1, l2, l3, m1, m2, m3, n1, n2, n3;
      double l1_r, l2_r, l3_r, m1_r, m2_r, m3_r, n1_r, n2_r, n3_r;
      l1   =  cos(theta); l1_r = -sin(theta);
      m1   =  sin(theta); m1_r =  cos(theta);
      n1   =  0.0;        n1_r =  0.0;
      l2   = -sin(theta); l2_r = -cos(theta);
      m2   =  cos(theta); m2_r = -sin(theta);
      n2   =  0.0;        n2_r =  0.0;
      l3   =  0.0;        l3_r =  0.0;
      m3   =  0.0;        m3_r =  0.0;
      n3   =  1.0;        n3_r =  0.0;
      
      ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6); T_strain.clear();
      T_strain(0,0) = l1*l1;
      T_strain(0,1) = m1*m1;
      T_strain(0,2) = n1*n1;
      T_strain(0,3) = l1*m1;
      T_strain(0,4) = m1*n1;
      T_strain(0,5) = l1*n1;
      
      T_strain(1,0) = l2*l2;
      T_strain(1,1) = m2*m2;
      T_strain(1,2) = n2*n2;
      T_strain(1,3) = l2*m2;
      T_strain(1,4) = m2*n2;
      T_strain(1,5) = l2*n2;
      
      T_strain(2,0) = l3*l3;
      T_strain(2,1) = m3*m3;
      T_strain(2,2) = n3*n3;
      T_strain(2,3) = l3*m3;
      T_strain(2,4) = m3*n3;
      T_strain(2,5) = l3*n3;
      
      T_strain(3,0) = 2*l1*l2;
      T_strain(3,1) = 2*m1*m2;
      T_strain(3,2) = 2*n1*n2;
      T_strain(3,3) = l1*m2+m1*l2;
      T_strain(3,4) = m1*n2+n1*m2;
      T_strain(3,5) = l1*n2+n1*l2;
      
      T_strain(4,0) = 2*l2*l3;
      T_strain(4,1) = 2*m2*m3;
      T_strain(4,2) = 2*n2*n3;
      T_strain(4,3) = l2*m3+m2*l3;
      T_strain(4,4) = m2*n3+n2*m3;
      T_strain(4,5) = l2*n3+n2*l3;
      
      T_strain(5,0) = 2*l1*l3;
      T_strain(5,1) = 2*m1*m3;
      T_strain(5,2) = 2*n1*n3;
      T_strain(5,3) = l1*m3+m1*l3;
      T_strain(5,4) = m1*n3+n1*m3;
      T_strain(5,5) = l1*n3+n1*l3;
      
      ublas::matrix<FieldData> T_strain_r_Theta;    T_strain_r_Theta.resize(6,6); T_strain_r_Theta.clear();
      T_strain_r_Theta(0,0) = 2*l1*l1_r;
      T_strain_r_Theta(0,1) = 2*m1*m1_r;
      T_strain_r_Theta(0,2) = 2*n1*n1_r;
      T_strain_r_Theta(0,3) = l1_r*m1 + l1*m1_r;
      T_strain_r_Theta(0,4) = m1_r*n1 + m1*n1_r;
      T_strain_r_Theta(0,5) = l1_r*n1 + l1*n1_r;
      
      T_strain_r_Theta(1,0) = 2*l2*l2_r;
      T_strain_r_Theta(1,1) = 2*m2*m2_r;
      T_strain_r_Theta(1,2) = 2*n2*n2_r;
      T_strain_r_Theta(1,3) = l2_r*m2 + l2*m2_r;
      T_strain_r_Theta(1,4) = m2_r*n2 + m2*n2_r;
      T_strain_r_Theta(1,5) = l2_r*n2 + l2*n2_r;
      
      T_strain_r_Theta(2,0) = 2*l3*l3_r;
      T_strain_r_Theta(2,1) = 2*m3*m3_r;
      T_strain_r_Theta(2,2) = 2*n3*n3_r;
      T_strain_r_Theta(2,3) = l3_r*m3 + l3*m3_r;
      T_strain_r_Theta(2,4) = m3_r*n3 + m3*n3_r;
      T_strain_r_Theta(2,5) = l3_r*n3 + l3*n3_r;
      
      T_strain_r_Theta(3,0) = 2*(l1_r*l2 + l1*l2_r);
      T_strain_r_Theta(3,1) = 2*(m1_r*m2 + m1*m2_r);
      T_strain_r_Theta(3,2) = 2*(n1_r*n2 + n1*n2_r);
      T_strain_r_Theta(3,3) = l1_r*m2 + l1*m2_r + m1_r*l2 + m1*l2_r;
      T_strain_r_Theta(3,4) = m1_r*n2 + m1*n2_r + n1_r*m2 + n1*m2_r;
      T_strain_r_Theta(3,5) = l1_r*n2 + l1*n2_r + n1_r*l2 + n1*l2_r;
      
      T_strain_r_Theta(4,0) = 2*(l2_r*l3 + l2*l3_r);
      T_strain_r_Theta(4,1) = 2*(m2_r*m3 + m2*m3_r);
      T_strain_r_Theta(4,2) = 2*(n2_r*n3 + n2*n3_r);
      T_strain_r_Theta(4,3) = l2_r*m3 + l2*m3_r + m2_r*l3 + m2*l3_r;
      T_strain_r_Theta(4,4) = m2_r*n3 + m2*n3_r + n2_r*m3 + n2*m3_r;
      T_strain_r_Theta(4,5) = l2_r*n3 + l2*n3_r + n2_r*l3 + n2*l3_r;
      
      T_strain_r_Theta(5,0) = 2*(l1_r*l3 + l1*l3_r);
      T_strain_r_Theta(5,1) = 2*(m1_r*m3 + m1*m3_r);
      T_strain_r_Theta(5,2) = 2*(n1_r*n3 + n1*n3_r);
      T_strain_r_Theta(5,3) = l1_r*m3 + l1*m3_r + m1_r*l3 + m1*l3_r;
      T_strain_r_Theta(5,4) = m1_r*n3 + m1*n3_r + n1_r*m3 + n1*m3_r;
      T_strain_r_Theta(5,5) = l1_r*n3 + l1*n3_r + n1_r*l3 + n1*l3_r;
      
      T_strain_r_Theta = T_strain_r_Theta*(M_PI/180);
      
      //cout<<"\n\nT_strain = "<<T_strain<<endl;
      
      ublas::matrix<FieldData> Mat1 = prod(Dmat_123,T_strain);
      ublas::matrix<FieldData> Mat2 = prod(trans(T_strain_r_Theta), Mat1);
      
      ublas::matrix<FieldData> Mat3 = prod(Dmat_123,T_strain_r_Theta);
      ublas::matrix<FieldData> Mat4 = prod(trans(T_strain), Mat3);
      
      Dmat_xyz = Mat2 + Mat4;
         
      PetscFunctionReturn(0);
    }
    
    
    virtual PetscErrorCode Dmat_Transformation_rs_Theta(double theta,
                                                       ublas::matrix<FieldData> Dmat_123,
                                                       ublas::matrix<FieldData> &Dmat_xyz) {
      PetscFunctionBegin;
      
      double l1,    l2,    l3,    m1,    m2,    m3,    n1,    n2,    n3;
      double l1_r,  l2_r,  l3_r,  m1_r,  m2_r,  m3_r,  n1_r,  n2_r,  n3_r;
      double l1_rs, l2_rs, l3_rs, m1_rs, m2_rs, m3_rs, n1_rs, n2_rs, n3_rs;
      l1   =  cos(theta); l1_r = -sin(theta); l1_rs = -cos(theta);
      m1   =  sin(theta); m1_r =  cos(theta); m1_rs = -sin(theta);
      n1   =  0.0;        n1_r =  0.0;        n1_rs =  0.0;
      l2   = -sin(theta); l2_r = -cos(theta); l2_rs =  sin(theta);
      m2   =  cos(theta); m2_r = -sin(theta); m2_rs = -cos(theta);
      n2   =  0.0;        n2_r =  0.0;        n2_rs =  0.0;
      l3   =  0.0;        l3_r =  0.0;        l3_rs =  0.0;
      m3   =  0.0;        m3_r =  0.0;        m3_rs =  0.0;
      n3   =  1.0;        n3_r =  0.0;        n3_rs =  0.0;
      
      ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6); T_strain.clear();
      T_strain(0,0) = l1*l1;
      T_strain(0,1) = m1*m1;
      T_strain(0,2) = n1*n1;
      T_strain(0,3) = l1*m1;
      T_strain(0,4) = m1*n1;
      T_strain(0,5) = l1*n1;
      
      T_strain(1,0) = l2*l2;
      T_strain(1,1) = m2*m2;
      T_strain(1,2) = n2*n2;
      T_strain(1,3) = l2*m2;
      T_strain(1,4) = m2*n2;
      T_strain(1,5) = l2*n2;
      
      T_strain(2,0) = l3*l3;
      T_strain(2,1) = m3*m3;
      T_strain(2,2) = n3*n3;
      T_strain(2,3) = l3*m3;
      T_strain(2,4) = m3*n3;
      T_strain(2,5) = l3*n3;
      
      T_strain(3,0) = 2*l1*l2;
      T_strain(3,1) = 2*m1*m2;
      T_strain(3,2) = 2*n1*n2;
      T_strain(3,3) = l1*m2+m1*l2;
      T_strain(3,4) = m1*n2+n1*m2;
      T_strain(3,5) = l1*n2+n1*l2;
      
      T_strain(4,0) = 2*l2*l3;
      T_strain(4,1) = 2*m2*m3;
      T_strain(4,2) = 2*n2*n3;
      T_strain(4,3) = l2*m3+m2*l3;
      T_strain(4,4) = m2*n3+n2*m3;
      T_strain(4,5) = l2*n3+n2*l3;
      
      T_strain(5,0) = 2*l1*l3;
      T_strain(5,1) = 2*m1*m3;
      T_strain(5,2) = 2*n1*n3;
      T_strain(5,3) = l1*m3+m1*l3;
      T_strain(5,4) = m1*n3+n1*m3;
      T_strain(5,5) = l1*n3+n1*l3;
      
      ublas::matrix<FieldData> T_strain_r_Theta;    T_strain_r_Theta.resize(6,6); T_strain_r_Theta.clear();
      T_strain_r_Theta(0,0) = 2*l1*l1_r;
      T_strain_r_Theta(0,1) = 2*m1*m1_r;
      T_strain_r_Theta(0,2) = 2*n1*n1_r;
      T_strain_r_Theta(0,3) = l1_r*m1 + l1*m1_r;
      T_strain_r_Theta(0,4) = m1_r*n1 + m1*n1_r;
      T_strain_r_Theta(0,5) = l1_r*n1 + l1*n1_r;
      
      T_strain_r_Theta(1,0) = 2*l2*l2_r;
      T_strain_r_Theta(1,1) = 2*m2*m2_r;
      T_strain_r_Theta(1,2) = 2*n2*n2_r;
      T_strain_r_Theta(1,3) = l2_r*m2 + l2*m2_r;
      T_strain_r_Theta(1,4) = m2_r*n2 + m2*n2_r;
      T_strain_r_Theta(1,5) = l2_r*n2 + l2*n2_r;
      
      T_strain_r_Theta(2,0) = 2*l3*l3_r;
      T_strain_r_Theta(2,1) = 2*m3*m3_r;
      T_strain_r_Theta(2,2) = 2*n3*n3_r;
      T_strain_r_Theta(2,3) = l3_r*m3 + l3*m3_r;
      T_strain_r_Theta(2,4) = m3_r*n3 + m3*n3_r;
      T_strain_r_Theta(2,5) = l3_r*n3 + l3*n3_r;
      
      T_strain_r_Theta(3,0) = 2*(l1_r*l2 + l1*l2_r);
      T_strain_r_Theta(3,1) = 2*(m1_r*m2 + m1*m2_r);
      T_strain_r_Theta(3,2) = 2*(n1_r*n2 + n1*n2_r);
      T_strain_r_Theta(3,3) = l1_r*m2 + l1*m2_r + m1_r*l2 + m1*l2_r;
      T_strain_r_Theta(3,4) = m1_r*n2 + m1*n2_r + n1_r*m2 + n1*m2_r;
      T_strain_r_Theta(3,5) = l1_r*n2 + l1*n2_r + n1_r*l2 + n1*l2_r;
      
      T_strain_r_Theta(4,0) = 2*(l2_r*l3 + l2*l3_r);
      T_strain_r_Theta(4,1) = 2*(m2_r*m3 + m2*m3_r);
      T_strain_r_Theta(4,2) = 2*(n2_r*n3 + n2*n3_r);
      T_strain_r_Theta(4,3) = l2_r*m3 + l2*m3_r + m2_r*l3 + m2*l3_r;
      T_strain_r_Theta(4,4) = m2_r*n3 + m2*n3_r + n2_r*m3 + n2*m3_r;
      T_strain_r_Theta(4,5) = l2_r*n3 + l2*n3_r + n2_r*l3 + n2*l3_r;
      
      T_strain_r_Theta(5,0) = 2*(l1_r*l3 + l1*l3_r);
      T_strain_r_Theta(5,1) = 2*(m1_r*m3 + m1*m3_r);
      T_strain_r_Theta(5,2) = 2*(n1_r*n3 + n1*n3_r);
      T_strain_r_Theta(5,3) = l1_r*m3 + l1*m3_r + m1_r*l3 + m1*l3_r;
      T_strain_r_Theta(5,4) = m1_r*n3 + m1*n3_r + n1_r*m3 + n1*m3_r;
      T_strain_r_Theta(5,5) = l1_r*n3 + l1*n3_r + n1_r*l3 + n1*l3_r;
      
      T_strain_r_Theta = T_strain_r_Theta*(M_PI/180);
      
      // The second-order derivative of rotation matrix
      ublas::matrix<FieldData> T_strain_rs_Theta;    T_strain_rs_Theta.resize(6,6); T_strain_rs_Theta.clear();
      T_strain_rs_Theta(0,0) = 2*l1_r*l1_r + 2*l1*l1_rs;
      T_strain_rs_Theta(0,1) = 2*m1_r*m1_r + 2*m1*m1_rs;
      T_strain_rs_Theta(0,2) = 2*n1_r*n1_r + 2*n1*n1_rs;
      T_strain_rs_Theta(0,3) = l1_rs*m1 + 2*l1_r*m1_r +  l1*m1_rs;
      T_strain_rs_Theta(0,4) = m1_rs*n1 + 2*m1_r*n1_r + m1*n1_rs;
      T_strain_rs_Theta(0,5) = l1_rs*n1 + 2*l1_r*n1_r + l1*n1_rs;
      
      T_strain_rs_Theta(1,0) = 2*l2_r*l2_r + 2*l2*l2_rs;
      T_strain_rs_Theta(1,1) = 2*m2_r*m2_r + 2*m2*m2_rs;
      T_strain_rs_Theta(1,2) = 2*n2_r*n2_r + 2*n2*n2_rs;
      T_strain_rs_Theta(1,3) = l2_rs*m2 + l2_r*m2_r + l2_r*m2_r + l2*m2_rs;
      T_strain_rs_Theta(1,4) = m2_rs*n2 + m2_r*n2_r + m2_r*n2_r + m2*n2_rs;
      T_strain_rs_Theta(1,5) = l2_rs*n2 + l2_r*n2_r + l2_r*n2_r + l2*n2_rs;
      
      T_strain_rs_Theta(2,0) = 2*l3_r*l3_r + 2*l3*l3_rs;
      T_strain_rs_Theta(2,1) = 2*m3_r*m3_r + 2*m3*m3_rs;
      T_strain_rs_Theta(2,2) = 2*n3_r*n3_r + 2*n3*n3_rs;
      T_strain_rs_Theta(2,3) = l3_rs*m3 + l3_r*m3_r + l3_r*m3_r + l3*m3_rs;
      T_strain_rs_Theta(2,4) = m3_rs*n3 + m3_r*n3_r + m3_r*n3_r + m3*n3_rs;
      T_strain_rs_Theta(2,5) = l3_rs*n3 + l3_r*n3_r + l3_r*n3_r + l3*n3_rs;
      
      T_strain_rs_Theta(3,0) = 2*(l1_rs*l2 + l1_r*l2_r + l1_r*l2_r + l1*l2_rs);
      T_strain_rs_Theta(3,1) = 2*(m1_rs*m2 + m1_r*m2_r + m1_r*m2_r + m1*m2_rs);
      T_strain_rs_Theta(3,2) = 2*(n1_rs*n2 + n1_r*n2_r + n1_r*n2_r + n1*n2_rs);
      T_strain_rs_Theta(3,3) = l1_rs*m2 + l1_r*m2_r + l1_r*m2_r + l1*m2_rs + m1_rs*l2 + m1_r*l2_r + m1_r*l2_r + m1*l2_rs;
      T_strain_rs_Theta(3,4) = m1_rs*n2 + m1_r*n2_r + m1_r*n2_r + m1*n2_rs + n1_rs*m2 + n1_r*m2_r + n1_r*m2_r + n1*m2_rs;
      T_strain_rs_Theta(3,5) = l1_rs*n2 + l1_r*n2_r + l1_r*n2_r + l1*n2_rs + n1_rs*l2 + n1_r*l2_r + n1_r*l2_r + n1*l2_rs;
      
      T_strain_rs_Theta(4,0) = 2*(l2_rs*l3 + l2_r*l3_r + l2_r*l3_r + l2*l3_rs);
      T_strain_rs_Theta(4,1) = 2*(m2_rs*m3 + m2_r*m3_r + m2_r*m3_r + m2*m3_rs);
      T_strain_rs_Theta(4,2) = 2*(n2_rs*n3 + n2_r*n3_r + n2_r*n3_r + n2*n3_rs);
      T_strain_rs_Theta(4,3) = l2_rs*m3 + l2_r*m3_r + l2_r*m3_r + l2*m3_rs + m2_rs*l3 + m2_r*l3_r + m2_r*l3_r + m2*l3_rs;
      T_strain_rs_Theta(4,4) = m2_rs*n3 + m2_r*n3_r + m2_r*n3_r + m2*n3_rs + n2_rs*m3 + n2_r*m3_r + n2_r*m3_r + n2*m3_rs;
      T_strain_rs_Theta(4,5) = l2_rs*n3 + l2_r*n3_r + l2_r*n3_r + l2*n3_rs + n2_rs*l3 + n2_r*l3_r + n2_r*l3_r + n2*l3_rs;
      
      T_strain_rs_Theta(5,0) = 2*(l1_rs*l3 + l1_r*l3_r + l1_r*l3_r + l1*l3_rs);
      T_strain_rs_Theta(5,1) = 2*(m1_rs*m3 + m1_r*m3_r + m1_r*m3_r + m1*m3_rs);
      T_strain_rs_Theta(5,2) = 2*(n1_rs*n3 + n1_r*n3_r + n1_r*n3_r + n1*n3_rs);
      T_strain_rs_Theta(5,3) = l1_rs*m3 + l1_r*m3_r + l1_r*m3_r + l1*m3_rs + m1_rs*l3 + m1_r*l3_r + m1_r*l3_r + m1*l3_rs;
      T_strain_rs_Theta(5,4) = m1_rs*n3 + m1_r*n3_r + m1_r*n3_r + m1*n3_rs + n1_rs*m3 + n1_r*m3_r + n1_r*m3_r + n1*m3_rs;
      T_strain_rs_Theta(5,5) = l1_rs*n3 + l1_r*n3_r + l1_r*n3_r + l1*n3_rs + n1_rs*l3 + n1_r*l3_r + n1_r*l3_r + n1*l3_rs;
      
      T_strain_rs_Theta = T_strain_rs_Theta*pow(M_PI/180,2);
      //cout<<"\n\nT_strain = "<<T_strain<<endl;
      
      ublas::matrix<FieldData> Mat1 = prod(Dmat_123,T_strain);
      ublas::matrix<FieldData> Mat2 = prod(trans(T_strain_rs_Theta), Mat1);
      
      ublas::matrix<FieldData> Mat3 = prod(Dmat_123,T_strain_r_Theta);
      ublas::matrix<FieldData> Mat4 = prod(trans(T_strain_r_Theta), Mat3);
      
      ublas::matrix<FieldData> Mat5 = prod(Dmat_123,T_strain_rs_Theta);
      ublas::matrix<FieldData> Mat6 = prod(trans(T_strain), Mat5);
      
      Dmat_xyz = Mat2 + Mat4 + Mat6;
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
	  virtual PetscErrorCode RVE_Dmat_Periodic(FieldInterface &m_field_RVE,
                                             int &nvars, int &nders,
                                             vector<string> &stochastic_fields,
                                             ublas::vector<double> matprop,
                                             int num_rvars,
                                             vector<string> vars_name) {
      PetscFunctionBegin;
      cout <<"\nHi from RVE_Dmat_Periodic!"<<endl;
      
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);   Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      Dmat_rs_EfEf.resize(6,6);     Dmat_rs_EfEf.clear();
      Dmat_rs_NUfNUf.resize(6,6);   Dmat_rs_NUfNUf.clear();
      
      /*****************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ****************************************************************************/
      
      //create matrices
      Vec dF,ddF,D,dD,ddD;
      
      vector<Vec> F(6);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F[0]); CHKERRQ(ierr);
      for(int ii = 1;ii<6;ii++) {
        ierr = VecDuplicate(F[0],&F[ii]); CHKERRQ(ierr);
      }
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD); CHKERRQ(ierr);
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;
      double PoissonRatio = 0.3;
      double YoungModulus_Fibre, PoissonRatio_Fibre;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D,F[0], RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus_Fibre = mydata.data.Young;
          PoissonRatio_Fibre = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre of isotropic material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            //cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Poissonp  = matprop(2);
           mydata.data.Poissonpz = matprop(3);
           mydata.data.Youngp    = matprop(4);
           mydata.data.Youngz    = matprop(5);
           mydata.data.Shearzp   = matprop(6);
           */
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre of transversely isotropic material:\n" << mydata;
        }
      }
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      //applied_strain.resize(6); applied_strain.clear(); applied_strain(0) = 1.0;
      
      MyElasticFEMethod MyFE(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D,F[0],"DISP_RVE");
      
      //ElasticFE_RVELagrange_Periodic MyFE_RVELagrangePeriodic(m_field_RVE,Aij,D,F[0],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ElasticFE_RVELagrange_Periodic_Multi_Rhs MyFE_RVELagrangePeriodic(m_field_RVE,Aij,D,F[0],F[1],F[2],F[3],F[4],F[5],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ElasticFE_RVELagrange_RigidBodyTranslation MyFE_RVELagrangeRigidBodyTrans(m_field_RVE,Aij,D,F[0],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank,"Lagrange_mul_disp_rigid_trans");
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecZeroEntries(F[ii]); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrangePeriodic);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE_rigid_trans",MyFE_RVELagrangeRigidBodyTrans);  CHKERRQ(ierr);
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecGhostUpdateBegin(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
      }
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      
      //----------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //----------------------------------------------------------------------------
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo;
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
      Vec Stress_Homo_r;
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
      
      //ublas::matrix<double> Dmat(6,6); Dmat.clear();
      ublas::vector<ublas::matrix<double> > Dmat_r(nvars);
      ublas::vector<ublas::matrix<double> > Dmat_rs(nvars);
      for(int irv=0; irv<nvars; irv++) {
        Dmat_r(irv).resize(6,6);   Dmat_r(irv).clear();
        Dmat_rs(irv).resize(6,6); Dmat_rs(irv).clear();
      }
      
      string VariableName;
      for(int ics = 0; ics<6; ics++) {
        
        cout<<"===============================================================\n";
        switch (ics) {
          case 0:
            cout<<"        Applied strain [1 0 0 0 0 0]^T\n"; break;
          case 1:
            cout<<"        Applied strain [0 1 0 0 0 0]^T\n"; break;
          case 2:
            cout<<"        Applied strain [0 0 1 0 0 0]^T\n"; break;
          case 3:
            cout<<"        Applied strain [0 0 0 1 0 0]^T\n"; break;
          case 4:
            cout<<"        Applied strain [0 0 0 0 1 0]^T\n"; break;
          case 5:
            cout<<"        Applied strain [0 0 0 0 0 1]^T\n"; break;
        }
        cout<<"===============================================================\n";
        
        ierr = VecZeroEntries(D); CHKERRQ(ierr);
        ierr = KSPSolve(solver,F[ics],D); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        //----------------------------------------------------------------------------
        // 3.2 Calculating zeroth-order homogenized stress using volume averaging theorem
        //----------------------------------------------------------------------------
        
        ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
        ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic(m_field_RVE,Aij,D,F[ics],&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressPeriodic);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec;
          VecGetArray(Stress_Homo, &avec);
          //cout<< "\nStress_Homo =\n";
          for(int kk=0; kk<6; kk++) {
            Dmat(kk,ics) = *avec;
            //cout.precision(15);
            //cout <<*avec<<endl;
            avec++;
          }
          VecRestoreArray(Stress_Homo,&avec);
        }
        
        //cout<< "\n\n";
        
        /*****************************************************************************
         *  4. SOLVE THE FIRST-ORDER AND SECOND-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
         *     1st order-[K][U_r] = -[K_r][U}
         *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
         *
         ****************************************************************************/
        
        int mat_rv = 0, PSFE_order = 0, idx_sf =0;
        for(int irv=1; irv <= num_rvars; irv++) {
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          VariableName = vars_name[irv]; //cout<<"\n variable name: "<<VariableName<<endl;
          if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
            mat_rv = 1; PSFE_order = 1; idx_sf = 0;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio in p-direction of fibre
            mat_rv = 1; PSFE_order = 1; idx_sf = 1;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
            mat_rv = 1; PSFE_order = 1; idx_sf = 2;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
            mat_rv = 1; PSFE_order = 1; idx_sf = 3;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUz);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
            mat_rv = 1; PSFE_order = 1; idx_sf = 4;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
            mat_rv = 1; PSFE_order = 1; idx_sf = 5;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
          }
          else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
            mat_rv = 1; PSFE_order = 1; idx_sf = 6;
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
          }
          VariableName.clear();
          
          if (mat_rv == 1) {
            ostringstream ss_field;
            ss_field << "DISP_RVE" << stochastic_fields[idx_sf];
            if (PSFE_order == 1) { // solution for first-order problem
              ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if (PSFE_order == 2) { // solution for second-order problem
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
              ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            //----------------------------------------------------------------------------
            // b. Calculating first-order homogenized stress
            //----------------------------------------------------------------------------
            ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
            
            ElasticFE_RVELagrange_Homogenized_Stress_Periodic MyFE_RVEHomoStressPeriodic_r(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressPeriodic_r);  CHKERRQ(ierr);
            
            if(pcomm->rank()==0) {
              PetscScalar    *avec_r;
              VecGetArray(Stress_Homo_r, &avec_r);
              //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
              //cout<< "\n"<<ss_field<<" = \n\n";
              for(int kk=0; kk<6; kk++) {
                //cout.precision(15);
                //cout<<*avec_r<<endl;
                //Dmat_r(irv-1)(kk,ics) = *avec_r;
                avec_r++;
              }
              VecRestoreArray(Stress_Homo_r,&avec_r);
            }
          }
          mat_rv = 0; PSFE_order = 0;
        }
      }
      
      // Asign results to corresponding matrix
      for(int irv=1; irv <= num_rvars; irv++) {
        VariableName = vars_name[irv];
        if (VariableName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          Dmat_r_Em = Dmat_r(irv-1); //cout<<"\nDmat_r_Em: "<<Dmat_r_Em<<endl;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUm = Dmat_r(irv-1);
        }
        else if (VariableName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUp = Dmat_r(irv-1);
        }
        else if (VariableName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUpz = Dmat_r(irv-1);
        }
        else if (VariableName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          Dmat_r_Ep = Dmat_r(irv-1);
        }
        else if (VariableName.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
          Dmat_r_Ez = Dmat_r(irv-1); //cout<<"\nDmat_r_Ez: "<<Dmat_r_Ez<<endl;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          Dmat_r_Gzp = Dmat_r(irv-1);
        }
        VariableName.clear();
      }
      
      //cout<<"\nDmat for periodic BC: "<<Dmat<<endl;
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      for(int ii = 0;ii<6;ii++) {
        ierr = VecDestroy(&F[ii]); CHKERRQ(ierr);
      }
      
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo_r); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    
  }
    
    
    virtual PetscErrorCode IndentifyMaterialVariable(int num_vars,
                                                     vector<string> vars_name,
                                                     vector<string> &NameMatVars,
                                                     int &num_mat_vars) {
      PetscFunctionBegin;
      // Reset "num_mat_vars", "NameMatVars"
      num_mat_vars = 0;
      NameMatVars.clear();
      
      string VarName;
      for (int irv=1;irv<=num_vars;irv++) {
        VarName = vars_name[irv];
        if (VarName.compare(0,2,"Em") == 0) {       // Young's modulus of matrix - isotropic
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUm") == 0) { // Poisson's ratio of matrix - isotropic
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUp") == 0) { // Poisson's ratio in p-direction of fibre
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUz") == 0) { // Poisson's ratio in z-direction of fibre
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,2,"Ep") == 0) {  // Young's modulus in p-direction of fibre
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,2,"Ez") == 0) {  // Young's modulus in longitudinal direction
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        else if (VarName.compare(0,3,"Gzp") == 0) { // Shear modulus in z-direction of fibre
          num_mat_vars ++;
          NameMatVars.push_back(VarName);
        }
        VarName.clear();
      }
      
      PetscFunctionReturn(0);
    }
    
    virtual PetscErrorCode SeperateRandomVariables(int num_vars,
                                                   vector<string> vars_name,
                                                   vector<string> &rve_vars_name,
                                                   vector<string> &ply_vars_name,
                                                   int &num_rve_vars,
                                                   int &num_ply_vars,
                                                   vector<int> &rve_vars_pos,
                                                   vector<int> &ply_vars_pos) {
                                                   
      PetscFunctionBegin;
      // Reset "num_rve_vars", "rve_vars_name"
      num_rve_vars = 0; num_ply_vars = 0;
      rve_vars_name.clear();
      
      string VarName;
      for (int irv=1;irv<=num_vars;irv++) {
        VarName = vars_name[irv];
        if (VarName.compare(0,2,"Em") == 0) {       // Young's modulus of matrix - isotropic
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUm") == 0) { // Poisson's ratio of matrix - isotropic
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUp") == 0) { // Poisson's ratio in p-direction of fibre
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUz") == 0) { // Poisson's ratio in z-direction of fibre
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,2,"Ep") == 0) {  // Young's modulus in p-direction of fibre
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,2,"Ez") == 0) {  // Young's modulus in longitudinal direction
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,3,"Gzp") == 0) { // Shear modulus in z-direction of fibre
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,2,"Ef") == 0) { // Young's modulus of fibre of isotropic material
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,3,"NUf") == 0) { // Poisson's ratio of fibre of isotropic material
          num_rve_vars ++;
          num_ply_vars ++;
          rve_vars_pos.push_back(irv);
          ply_vars_pos.push_back(irv);
          rve_vars_name.push_back(VarName);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,5,"force") == 0) { // Applied force on structure
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,11,"orientation") == 0) { // Ply angle
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,6,"theta1") == 0) { // Ply angle for 1st ply
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,6,"theta2") == 0) { // Ply angle for 2nd ply
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,6,"theta3") == 0) { // Ply angle for 3rd ply
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        else if (VarName.compare(0,6,"theta4") == 0) { // Ply angle for 4th ply
          num_ply_vars ++;
          ply_vars_pos.push_back(irv);
          ply_vars_name.push_back(VarName);
        }
        VarName.clear();
      }
      PetscFunctionReturn(0);
    }
    
    
    virtual PetscErrorCode RVE_Dmat_Disp(FieldInterface &m_field_RVE,
                                         vector<string> stochastic_fields,
                                         ublas::vector<double> matprop,
                                         int num_vars,
                                         vector<string> vars_name,
                                         int PSFE_order) {
                                     
      PetscFunctionBegin;
      cout<<"\n"<<endl;
      cout<<"///////////////////////////////////////////////////////////"<<endl;
      cout<<"//                                                       //"<<endl;
      cout<<"//               Hi from RVE_Dmat_Disp!                  //"<<endl;
      cout<<"//                                                       //"<<endl;
      cout<<"///////////////////////////////////////////////////////////"<<endl;
      
      PetscErrorCode ierr;
      
      // Identify random variables for material properties
      vector<string> rve_vars_name;
      vector<string> ply_vars_name;
      int num_rve_vars, num_ply_vars;
      vector<int> rve_vars_pos;
      vector<int> ply_vars_pos;
      
      SeperateRandomVariables(num_vars, vars_name, rve_vars_name, ply_vars_name,
                              num_rve_vars, num_ply_vars, rve_vars_pos, ply_vars_pos);
      
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);   Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      Dmat_rs_EfEf.resize(6,6);     Dmat_rs_EfEf.clear();
      Dmat_rs_NUfNUf.resize(6,6);   Dmat_rs_NUfNUf.clear();
      
      /*****************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ****************************************************************************/
      
      //create matrices
      Vec dF,ddF,D,dD,ddD;
      
      vector<Vec> F(6);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F[0]); CHKERRQ(ierr);
      for(int ii = 1;ii<6;ii++) {
        ierr = VecDuplicate(F[0],&F[ii]); CHKERRQ(ierr);
      }
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD); CHKERRQ(ierr);
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;
      double PoissonRatio = 0.3;
      double YoungModulus_Fibre, PoissonRatio_Fibre;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D,F[0], RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_vars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_vars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus_Fibre = mydata.data.Young;
          PoissonRatio_Fibre = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre of isotropic material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_vars;i++) {
            ParameterName = vars_name[i];
            //cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Poissonp  = matprop(2);
           mydata.data.Poissonpz = matprop(3);
           mydata.data.Youngp    = matprop(4);
           mydata.data.Youngz    = matprop(5);
           mydata.data.Shearzp   = matprop(6);
           */
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre of transversely isotropic material:\n" << mydata;
        }
      }
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation but required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      //applied_strain.resize(6); applied_strain.clear(); applied_strain(0) = 1.0;
      
      MyElasticFEMethod MyFE(m_field_RVE,Aij,D,F[0],LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D,F[0],"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D,F[0],F[1],F[2],F[3],F[4],F[5],applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecZeroEntries(F[ii]); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      }
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      for(int ii = 0; ii<6; ii++) {
        ierr = VecGhostUpdateBegin(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F[ii],ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F[ii]); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F[ii]); CHKERRQ(ierr);
      }
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      
      
      /*****************************************************************************
       *
       *  3. SOLVE THE ZEROTH-ORDER FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ****************************************************************************/
      
      //----------------------------------------------------------------------------
      // 3.1 Prepapration for solving FE equations
      //----------------------------------------------------------------------------
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo;
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
      Vec Stress_Homo_r;
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo_r);  CHKERRQ(ierr);
      
      //ublas::matrix<double> Dmat(6,6); Dmat.clear();
//      ublas::vector<ublas::matrix<double> > Dmat_r(num_mat_vars);
//      ublas::matrix<ublas::matrix<double> > Dmat_rs(num_mat_vars,num_mat_vars);
//      for(int irv=0; irv<num_mat_vars; irv++) {
//        Dmat_r(irv).resize(6,6);   Dmat_r(irv).clear();
//        for(int jrv=0; jrv<num_mat_vars; jrv++) {
//          Dmat_rs(irv,jrv).resize(6,6); Dmat_rs(irv,jrv).clear();
//        }
//      }
      
      RVE_Dmat_r.resize(num_vars);
      RVE_Dmat_rs.resize(num_vars,num_vars);
      for(int irv=0; irv<num_vars; irv++) {
        RVE_Dmat_r(irv).resize(6,6);   RVE_Dmat_r(irv).clear();
        for(int jrv=0; jrv<num_vars; jrv++) {
          RVE_Dmat_rs(irv,jrv).resize(6,6); RVE_Dmat_rs(irv,jrv).clear();
        }
      }
      
      string VarName, VarName_1st, VarName_2nd;
      string first_var, second_var, first_field, material_type, material_function;
      
      for(int ics = 0; ics<6; ics++) {
        
        cout<<"===============================================================\n";
        switch (ics) {
          case 0:
            cout<<"        Applied strain [1 0 0 0 0 0]^T\n"; break;
          case 1:
            cout<<"        Applied strain [0 1 0 0 0 0]^T\n"; break;
          case 2:
            cout<<"        Applied strain [0 0 1 0 0 0]^T\n"; break;
          case 3:
            cout<<"        Applied strain [0 0 0 1 0 0]^T\n"; break;
          case 4:
            cout<<"        Applied strain [0 0 0 0 1 0]^T\n"; break;
          case 5:
            cout<<"        Applied strain [0 0 0 0 0 1]^T\n"; break;
        }
        cout<<"===============================================================\n";
        
        //----------------------------------------------------------------------
        // 3.1 Solving the zeroth-order equation
        //     [K][U] = [F]
        //----------------------------------------------------------------------
        
        ierr = VecZeroEntries(D); CHKERRQ(ierr);
        ierr = KSPSolve(solver,F[ics],D); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field_RVE,Aij,D,F[ics],&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec;
          VecGetArray(Stress_Homo, &avec);
          //cout<< "\nStress_Homo =\n";
          for(int kk=0; kk<6; kk++) {
            Dmat(kk,ics) = *avec;
            //cout.precision(15);
            //cout <<*avec<<endl;
            avec++;
          }
          VecRestoreArray(Stress_Homo,&avec);
        }

        //cout<< "\n\n";
        
        //----------------------------------------------------------------------
        // 3.2 Solving the first-order equation
        //     [K][U_r] = -[K_r][U]
        //----------------------------------------------------------------------
        
        int ix_mat_rv = 0, idx_sf =0;
        
        for(int irv = 0; irv < num_rve_vars; irv++) {
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          VarName = rve_vars_name[irv];
          
          if (VarName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D,dF,"DISP_RVE","Young","isotropic","matrix");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D,dF,"DISP_RVE","Poisson","isotropic","matrix");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUz(m_field_RVE,Aij,D,dF,"DISP_RVE","PoissonZ","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUz);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D,dF,"DISP_RVE","YoungZ","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
          }
          else if (VarName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
            Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D,dF,"DISP_RVE","ShearZP","transversely_isotropic","reinforcement");
            ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
          }
          VarName.clear();

          ostringstream ss_field;
          ss_field << "DISP_RVE" << stochastic_fields[irv];
          
          //------------------------------------------------------------------
          // a. Getting the 1st-order derivative of nodal displacement
          //------------------------------------------------------------------
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          
          //------------------------------------------------------------------
          // b. Calculating first-order homogenized stress
          //------------------------------------------------------------------
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int kk=0; kk<6; kk++) {
              //cout.precision(15);
              //cout<<*avec_r<<endl;
              RVE_Dmat_r(irv)(kk,ics) = *avec_r;
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r,&avec_r);
          }

          
        }

        //----------------------------------------------------------------------
        // 3.3 Solving the second-order equation
        //     [K][U_rs] = -[K_rs][U]-2[K_r][U_s]
        //----------------------------------------------------------------------
        //
        if (PSFE_order == 2) {
          cout<<"\n==========================="<<endl;
          cout<<"     Second order          "<<endl;
          cout<<"==========================="<<endl;
          for(int irv = 0; irv < num_rve_vars; irv++) {
            VarName_1st = rve_vars_name[irv];
            if (VarName_1st.compare(0,2,"Em") == 0) {       // due to Young's modulus of matrix - isotropic
              first_var         = "Young";
              first_field       = "DISP_RVE_r_Em";
              material_function = "matrix";
              material_type     = "isotropic";
              idx_sf            = 0;
            }
            else if (VarName_1st.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
              first_var         = "Poisson";
              first_field       = "DISP_RVE_r_NUm";
              material_function = "matrix";
              material_type     = "isotropic";
              idx_sf            = 1;
            }
            else if (VarName_1st.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
              first_var         = "PoissonP";
              first_field       = "DISP_RVE_r_NUp";
              material_function = "reinforcement";
              material_type     = "transversely_isotropic";
              idx_sf            = 2;
            }
            else if (VarName_1st.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
              first_var         = "PoissonZ";
              first_field       = "DISP_RVE_r_NUpz";
              material_function = "reinforcement";
              material_type     = "transversely_isotropic";
              idx_sf            = 3;
            }
            else if (VarName_1st.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
              first_var         = "YoungP";
              first_field       = "DISP_RVE_r_Ep";
              material_function = "reinforcement";
              material_type     = "transversely_isotropic";
              idx_sf            = 4;
            }
            else if (VarName_1st.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
              first_var         = "YoungZ";
              first_field       = "DISP_RVE_r_Ez";
              material_function = "reinforcement";
              material_type     = "transversely_isotropic";
              idx_sf            = 5;
            }
            else if (VarName_1st.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
              first_var         = "ShearZP";
              first_field       = "DISP_RVE_r_Gzp";
              material_function = "reinforcement";
              material_type     = "transversely_isotropic";
              idx_sf            = 6;
            }
            idx_sf = irv;
            first_field.clear();
            ostringstream first_field;
            first_field << "DISP_RVE" << stochastic_fields[idx_sf];
            //VarName_1st.clear();
            
            for(int jrv=irv; jrv < num_rve_vars; jrv++) {
              VarName_2nd = rve_vars_name[jrv]; cout<<"\n The 1st and 2nd variable names are: "<<VarName_1st<<"\t"<<VarName_2nd<<endl;
              if (VarName_2nd.compare(0,2,"Em") == 0) {       // due to Young's modulus of matrix - isotropic
                second_var = "Young";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,3,"NUm") == 0) { // due to Poisson's ratio of matrix - isotropic
                second_var = "Poisson";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
                second_var = "PoissonP";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
                second_var = "PoissonZ";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
                second_var = "YoungP";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
                second_var = "YoungZ";
                ix_mat_rv  = 1;
              }
              else if (VarName_2nd.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
                second_var = "ShearZP";
                ix_mat_rv  = 1;
              }
              VarName_2nd.clear();
              
              //if (irv == jrv){cout<<"The variable is "<<irv<<"\t"<<jrv<<endl;
              {
                //Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs(m_field_RVE,Aij,D,ddF,"DISP_RVE",first_field,first_var,second_var, material_type, material_function);
                Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs(m_field_RVE,Aij,D,ddF,"DISP_RVE",first_field.str().c_str(),first_var,second_var, material_type, material_function);
                if (material_type.compare(0,5,"trans") == 0) {
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs);  CHKERRQ(ierr);
                } else {
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs);  CHKERRQ(ierr);
                }
                second_var.clear();
                //first_field.clear(); first_var.clear();  material_type.clear(); material_function.clear();
                
                if (ix_mat_rv == 1) {
                  ostringstream ss_field;
                  ss_field << "DISP_RVE" << stochastic_fields[idx_sf];
                  
                  //--------------------------------------------------------------
                  // a. Getting the 2nd-order derivative of nodal displacement
                  //--------------------------------------------------------------
                  ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
                  ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
                  
                  ierr = KSPSolve(solver,ddF,ddD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
                  // ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
                  
                  //--------------------------------------------------------------
                  // b. Calculating first-order homogenized stress
                  //--------------------------------------------------------------
                  ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
                  
                  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r(m_field_RVE,Aij,dD,dF,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
                  ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r);  CHKERRQ(ierr);
                  
                  if(pcomm->rank()==0) {
                    PetscScalar    *avec_r;
                    VecGetArray(Stress_Homo_r, &avec_r);
                    //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
                    //cout<< "\n"<<ss_field<<" = \n\n";
                    for(int kk=0; kk<6; kk++) {
                      //cout.precision(15);
                      //cout<<*avec_r<<endl;
                      RVE_Dmat_rs(rve_vars_pos[irv]-1,rve_vars_pos[jrv]-1)(kk,ics) = *avec_r;
                      RVE_Dmat_rs(rve_vars_pos[jrv]-1,rve_vars_pos[irv]-1)(kk,ics) = *avec_r;
                      avec_r++;
                    }
                    VecRestoreArray(Stress_Homo_r,&avec_r);
                  }
                  
                }
                ix_mat_rv = 0;
              }
            }
            VarName_1st.clear();
            first_field.clear(); first_var.clear(); material_type.clear(); material_function.clear();
          }
        }
      }
      
      // Asign results to corresponding matrix
      for(int irv = 0; irv < num_rve_vars; irv++) {
        VarName = rve_vars_name[irv];
        if (VarName.compare(0,2,"Em") == 0) { // due to Young's modulus of matrix - isotropic
          Dmat_r_Em = RVE_Dmat_r(rve_vars_pos[irv]-1); //cout<<"\nDmat_r_Em: "<<Dmat_r_Em<<endl;
        }
        else if (VarName.compare(0,3,"NUm") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUm = RVE_Dmat_r(rve_vars_pos[irv]-1);
        }
        else if (VarName.compare(0,3,"NUp") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUp = RVE_Dmat_r(rve_vars_pos[irv]-1);
        }
        else if (VarName.compare(0,3,"NUz") == 0) { // due to Poisson's ratio in p-direction of fibre
          Dmat_r_NUpz = RVE_Dmat_r(rve_vars_pos[irv]-1);
        }
        else if (VarName.compare(0,2,"Ep") == 0) { // due to Young's modulus in p-direction of fibre
          Dmat_r_Ep = RVE_Dmat_r(rve_vars_pos[irv]-1);
        }
        else if (VarName.compare(0,2,"Ez") == 0) { // due to Young's modulus in longitudinal direction
          Dmat_r_Ez = RVE_Dmat_r(rve_vars_pos[irv]-1); //cout<<"\nDmat_r_Ez: "<<Dmat_r_Ez<<endl;
        }
        else if (VarName.compare(0,3,"Gzp") == 0) { // due to shear modulus in z-direction of fibre
          Dmat_r_Gzp = RVE_Dmat_r(rve_vars_pos[irv]-1);
        }
        
        // ==============
        if (PSFE_order == 2) {
          for (int jrv = 0; jrv < num_rve_vars; jrv++) {
            VarName_2nd = rve_vars_name[jrv];
            cout<<"\n The 1st and 2nd variable names are: "<<VarName<<"\t"<<VarName_2nd<<endl;
            cout<<RVE_Dmat_rs(rve_vars_pos[irv]-1,rve_vars_pos[jrv]-1)<<endl;
          }
          VarName_2nd.clear();
        }
        VarName.clear();
      }
      
      //cout<<"\nDmat for periodic BC: "<<Dmat<<endl;
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      for(int ii = 0;ii<6;ii++) {
        ierr = VecDestroy(&F[ii]); CHKERRQ(ierr);
      }
      
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);
      ierr = VecDestroy(&Stress_Homo_r); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
      
    }
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    virtual PetscErrorCode Micro_FE_Dmat(FieldInterface &m_field_RVE,
                                         int &nvars, int &nders,
                                         vector<string> &stochastic_fields,
                                         ublas::vector<double> matprop,
                                         int num_rvars,
                                         vector<string> vars_name) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      PetscErrorCode ierr;
     
      RVE_Dmat_r.resize(num_rvars);
      for (int i = 0; i<num_rvars; i++) {
        RVE_Dmat_r(i).resize(6,6); RVE_Dmat_r(i).clear();
      }
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      Dmat_r_Ef.resize(6,6);   Dmat_r_Ef.clear();
      Dmat_r_NUf.resize(6,6);  Dmat_r_NUf.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);   Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      Dmat_rs_EfEf.resize(6,6);     Dmat_rs_EfEf.clear();
      Dmat_rs_NUfNUf.resize(6,6);   Dmat_rs_NUfNUf.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;
      double PoissonRatio = 0.3;
      double YoungModulus_Fibre, PoissonRatio_Fibre;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        // ---------------------------------
        //
        // Modify matrix material properties
        //
        // ---------------------------------
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify matrix fibre properties
        // Case 1: Isotropic material
        //
        // ---------------------------------
        if (name.compare(0,19,"MAT_FIBRE_ISOTROPIC") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Ef") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {//cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus_Fibre = mydata.data.Young;
          PoissonRatio_Fibre = mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
        
        // ---------------------------------
        //
        // Modify fibre material properties
        // Case 2: Transversely isotropic material
        //
        // ---------------------------------
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ef") == 0) {
              mydata.data.Youngz = matprop(i-1);
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUf") == 0) {
              mydata.data.Poissonp  = matprop(i-1);
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if ((ParameterName.compare(0,3,"NUf") == 0) || (ParameterName.compare(0,2,"Ef") == 0)) {
              mydata.data.Shearzp = mydata.data.Youngz/(2*(1+mydata.data.Poissonp));
            }
            ParameterName.clear();
          }
          
          /*
           mydata.data.Poissonp  = matprop(2);
           mydata.data.Poissonpz = matprop(3);
           mydata.data.Youngp    = matprop(4);
           mydata.data.Youngz    = matprop(5);
           mydata.data.Shearzp   = matprop(6);
           */
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs "<<endl;
      
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
      VecGetArray(Stress_Homo, &avec);
      for (int ii=0; ii<6; ii++) {
        Dmat(ii,0)=*avec;
        avec++;
      }
      
      /*if (pcomm->rank()==0) {
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,0)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;//cout<<"\n\nThe random variable is Ef"<<endl;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;//cout<<"\n\nThe random variable is NUf"<<endl;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++){
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,0) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,0) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,0) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,0) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,0) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,0) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,0) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,0) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,0) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,0) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,0) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,0) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,0) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,0) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
        
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      // solve for F2 and D2
      ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,1)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,1)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,1) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,1) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,1) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,1) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,1) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,1) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,1) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,1) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,1) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,1) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,1) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,1) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,1) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,1) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F3 and D3
      ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,2)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,2)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        VariableName.clear();
        
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,2) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,2) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,2) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,2) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,2) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,2) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,2) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,2) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,2) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,2) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,2) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,2) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,2) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,2) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F4 and D4
      ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,3)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,3)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,3) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,3) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,3) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,3) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,3) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,3) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,3) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,3) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,3) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,3) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,3) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,3) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,3) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,3) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F5 and D5
      ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,4)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,4)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,4) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,4) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,4) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,4) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,4) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,4) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,4) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,4) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,4) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,4) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,4) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,4) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,4) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,4) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F6 and D6
      ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,5)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,5)<<endl;
       }
       cout<< "\n\n";
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // due to Young's modulus of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ef(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ef);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // due to Poisson's ratio of fibre - isotropic material
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUf(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 14) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 15) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 16) { // 2nd order derivative due to Young's modulus of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EfEf(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ef","Young","Young", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EfEf);  CHKERRQ(ierr);
        }
        else if (idx_disp == 17) { // 2nd order derivative due to Poisson's ratio of fibre - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUfNUf(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUf","Poisson","Poisson", "isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUfNUf);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,5) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,5) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,5) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,5) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,5) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,5) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,5) = *avec_r; break;
                case 7:
                  Dmat_r_Ef(irow,5) = *avec_r; break;
                case 8:
                  Dmat_r_NUf(irow,5) = *avec_r; break;
                case 9:
                  Dmat_rs_EmEm(irow,5) = *avec_r; break;
                case 10:
                  Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
                case 11:
                  Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
                case 12:
                  Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
                case 13:
                  Dmat_rs_EpEp(irow,5) = *avec_r; break;
                case 14:
                  Dmat_rs_EzEz(irow,5) = *avec_r; break;
                case 15:
                  Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
                case 16:
                  Dmat_rs_EfEf(irow,5) = *avec_r; break;
                case 17:
                  Dmat_rs_NUfNUf(irow,5) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      ierr = VecDestroy(&F1); CHKERRQ(ierr);
      ierr = VecDestroy(&F2); CHKERRQ(ierr);
      ierr = VecDestroy(&F3); CHKERRQ(ierr);
      ierr = VecDestroy(&F4); CHKERRQ(ierr);
      ierr = VecDestroy(&F5); CHKERRQ(ierr);
      ierr = VecDestroy(&F6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    
    //================
    //
    //  This version of Micro_FE_REL intends to be feasible for both FORM and SORM.
    //
    //  History:
    //    2016.07.19: create
    //
    //================
    
    PetscErrorCode Macro_FE_REL_FSORM(FieldInterface &m_field_Macro,
                                      int num_rve_vars, int num_ply_vars,
                                      vector<string> rve_vars_name,
                                      vector<string> ply_vars_name,
                                      vector<string> stochastic_fields_ply,
                                      ublas::vector<double> TheVariables,
                                      ublas::vector<double> PlyAngle,
                                      PetscInt NO_Layers,
                                      int PSFE_order) {
      PetscFunctionBegin;
      cout<<"\n"<<endl;
      cout<<"///////////////////////////////////////////////////////\n//"<<endl;
      cout<<"//        Hi from macrolevel FE calculation!           \n//"<<endl;
      cout<<"/////////////////////////////////////////////////////////\n"<<endl;
      
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      Vec dF, dD;
      Vec ddF, ddD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&ddF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&ddD); CHKERRQ(ierr);
      
      /*****************************************************************************
       *
       *  1. Read the saved Dmat mechancial
       *     (from the computational homgenisaiton of the 0deg RVE)
       *
       ****************************************************************************/
      
      double theta;
      
      /*Dmat.clear();
       Dmat(0,0) = 1.2894e5; Dmat(0,1) = 0.0525e5; Dmat(0,2) = 0.0525e5;
       Dmat(1,0) = 0.0525e5; Dmat(1,1) = 0.1331e5; Dmat(1,2) = 0.0545e5;
       Dmat(2,0) = 0.0525e5; Dmat(2,1) = 0.0545e5; Dmat(2,2) = 0.1331e5;
       Dmat(3,3) = 0.0660e5; Dmat(4,4) = 0.0393e5; Dmat(5,5) = 0.0660e5;*/
      // First-layer
      theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
      Dmat_1st_Ply.resize(6,6);   Dmat_1st_Ply.clear();//cout<<"\nThe Dmat is: "<<Dmat<<endl;
      ierr = Dmat_Transformation(theta, Dmat, Dmat_1st_Ply); CHKERRQ(ierr);
      // Second-layer
      if (NO_Layers > 1) {
        theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_2nd_Ply.resize(6,6);   Dmat_2nd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_2nd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 2: \t"<<PlyAngle(1)<<"\n"<<Dmat_2nd_Ply<<endl;
      }
      // Third-layer
      if (NO_Layers > 2) {
        theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_3rd_Ply.resize(6,6);   Dmat_3rd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_3rd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 3: \t"<<PlyAngle(2)<<"\n"<<Dmat_3rd_Ply<<endl;
      }
      // Fourth-layer
      if (NO_Layers > 3) {
        theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_4th_Ply.resize(6,6);   Dmat_4th_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_4th_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 4: \t"<<PlyAngle(3)<<"\n"<<Dmat_4th_Ply<<endl;
      }
      
      /*************************************************************************
       *
       *  2. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
        MyElasticFEMethod_Macro(FieldInterface& _m_field_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
        FE2_ElasticFEMethod(_m_field_Macro,_A,_D,_F, _Dmat, _field_name) {};
        
        virtual PetscErrorCode RhsAndLhs() {
          PetscFunctionBegin;
          
          ierr = Lhs(); CHKERRQ(ierr);
          
          PetscFunctionReturn(0);
        }
      };
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      MyElasticFEMethod_Macro my_fe_1st_ply (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_2nd_ply (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_3rd_ply (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_4th_ply (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      // First layer
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_ply);  CHKERRQ(ierr);
      // Second layer
      if (NO_Layers > 1) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_ply);  CHKERRQ(ierr);
      }
      // Third layer
      if (NO_Layers > 2) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_ply);  CHKERRQ(ierr);
      }
      // Fourth layer
      if (NO_Layers > 3) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_ply);  CHKERRQ(ierr);
      }
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      
      //
      // Check whether force is considered as random variable or not
      //
      int idx_force = 100;
      for (int ii = 0; ii<num_ply_vars; ii++) {
        string VariableName;
        VariableName = ply_vars_name[ii];
        if (VariableName.compare(0,5,"force") == 0) {
          idx_force = ii;
        }
      }
      
      if (idx_force==100) { // Applied force is constant value or not a variable
        MetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      } else {
        MyMetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,
                                                          neumann_forces,F,
                                                          "DISP_MACRO",
                                                          "MESH_NODE_POSITIONS",
                                                          TheVariables(idx_force-1)); CHKERRQ(ierr);
      }
      
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      // elastic analys
      ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      
      cout<<"================================================\n";
      cout<<"  The zeroth-order FE equation has been solved. \n";
      cout<<"================================================\n";
      
      /*************************************************************************
       *
       *  4. SOLVE THE FIRST-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *
       ************************************************************************/
      
      /*************************************************************************
       *
       * 4.1 Case 1: Constituent maaterial properties are random variables
       *
       ************************************************************************/
      int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
      string VariableName;
      int var_pos;
      
      //int ply_mat_size = RVE_Dmat_r.size();
      
      Ply_1st_Dmat_r.resize(num_ply_vars);
      Ply_2nd_Dmat_r.resize(num_ply_vars);
      Ply_3rd_Dmat_r.resize(num_ply_vars);
      Ply_4th_Dmat_r.resize(num_ply_vars);
      
      for (int ivar=0; ivar<num_rve_vars; ivar++) {
        ierr = VecZeroEntries(dD); CHKERRQ(ierr);
        ierr = VecZeroEntries(dF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        // First-layer
        theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
        Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
        ierr = Dmat_Transformation(theta, RVE_Dmat_r(ivar), Ply_1st_Dmat_r(ivar)); CHKERRQ(ierr);
        // Second-layer
        if (NO_Layers > 1) {
          theta = PlyAngle(1)*(M_PI/180.0);  //rotation angle about the Z-axis
          Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation(theta, RVE_Dmat_r(ivar), Ply_2nd_Dmat_r(ivar)); CHKERRQ(ierr);
        }
        // Third-layer
        if (NO_Layers > 2) {
          theta = PlyAngle(2)*(M_PI/180.0);  //rotation angle about the Z-axis
          Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation(theta, RVE_Dmat_r(ivar), Ply_3rd_Dmat_r(ivar)); CHKERRQ(ierr);
        }
        // Fourth layer
        if (NO_Layers > 3) {
          theta = PlyAngle(3)*(M_PI/180.0);  //rotation angle about the Z-axis
          Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation(theta, RVE_Dmat_r(ivar), Ply_4th_Dmat_r(ivar)); CHKERRQ(ierr);
        }
        
        FE2_Rhs_r_PSFEM my_fe2_k_r_1st_Ply(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
        FE2_Rhs_r_PSFEM my_fe2_k_r_2nd_Ply(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
        FE2_Rhs_r_PSFEM my_fe2_k_r_3rd_Ply(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
        FE2_Rhs_r_PSFEM my_fe2_k_r_4th_Ply(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
        
        DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r(m_field_Macro,"DISP_MACRO",A,dD,dF);
        ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
        
        // First layer
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_1st_Ply);  CHKERRQ(ierr);
        // Second layer
        if (NO_Layers > 1) {
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_2nd_Ply);  CHKERRQ(ierr);
        }
        // Third layer
        if (NO_Layers > 2) {
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_3rd_Ply);  CHKERRQ(ierr);
        }
        // Fourth layer
        if (NO_Layers > 3) {
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_4th_Ply);  CHKERRQ(ierr);
        }
        
        ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r); CHKERRQ(ierr);
        
        //
        // post-processing
        ostringstream ss_field;
        ss_field.str(""); ss_field.clear();
        ss_field << "DISP_MACRO" << stochastic_fields_ply[ivar];
        
        ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
        
        //cout<<"First order derivative of dD"<<endl;
        //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        
        ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        
        //idx_mat_rv = 0;
        cout<<"===========================================================================\n";
        cout<<"  The first-order FE equation has been solved for "<<ss_field.str().c_str()<<"."<<endl;
        cout<<"===========================================================================\n";
        //cout<<"Solving the first-order equation "<<ss_field.str().c_str()<<" is finish. \n";
      }
      
      /*************************************************************************
       *
       * 4.2 Case 2: Ply geometrical parameters are random variables
       *
       ************************************************************************/
      
      for (int ivar=0; ivar<num_ply_vars; ivar++) {
        VariableName.clear(); VariableName = ply_vars_name[ivar];
        
        ostringstream first_field;
        first_field.str(""); first_field.clear();
        first_field << "DISP_MACRO" << stochastic_fields_ply[ivar];
        
        if (VariableName.compare(0,11,"orientation") == 0) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_1st_Dmat_r(ivar)); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0);  //rotation angle about the Z-axis
            Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_2nd_Dmat_r(ivar)); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0);  //rotation angle about the Z-axis
            Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_3rd_Dmat_r(ivar)); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0);  //rotation angle about the Z-axis
            Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_4th_Dmat_r(ivar)); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_1st_Ply(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_2nd_Ply(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_3rd_Ply(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_4th_Ply(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
          
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Theta_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Theta_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Theta_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Theta_4th_Ply);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for Theta."<<endl;
          cout<<"===========================================================================\n";
        }
        else if (VariableName.compare(0,6,"theta1") == 0) { // first layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_1st_Dmat_r(ivar)); CHKERRQ(ierr);
          //cout<<"\n\nThe first ply angle is "<<PlyAngle(0)<<endl;
          //cout<<"\n\n Dmat_1st_Ply_r_Theta_1 is "<<Dmat_1st_Ply_r_Theta_1<<endl;
          // Second-layer
          if (NO_Layers > 1) {
            Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          }
          // Third-layer
          if (NO_Layers > 2) {
            Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
          
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_1(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_1);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          // ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_1st_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for Theta_1."<<endl;
          cout<<"===========================================================================\n";
        }
        else if (VariableName.compare(0,6,"theta2") == 0) { // second layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          // Second-layer
          theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
          Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_2nd_Dmat_r(ivar)); CHKERRQ(ierr);
          // Third-layer
          if (NO_Layers > 2) {
            Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_2(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_2nd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for Theta_2."<<endl;
          cout<<"===========================================================================\n";
        }
        else if (VariableName.compare(0,6,"theta3") == 0) { // third layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          // Second-layer
          Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          // Third-layer
          theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
          Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_3rd_Dmat_r(ivar)); CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
          
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_3(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_3);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_3rd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for Theta_3."<<endl;
          cout<<"===========================================================================\n";
        }
        else if (VariableName.compare(0,6,"theta4") == 0) { // fourth layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          // Second-layer
          Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          // Third-layer
          Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          // Fourth layer
          theta = PlyAngle(3)*(M_PI/180.0);    //rotation angle about the Z-axis
          Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Ply_4th_Dmat_r(ivar)); CHKERRQ(ierr);
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Ply_1st_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Ply_2nd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Ply_3rd_Dmat_r(ivar),"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Ply_4th_Dmat_r(ivar),"DISP_MACRO");
          
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_4(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Fourth layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_4);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_4th_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for Theta_4."<<endl;
          cout<<"===========================================================================\n";
        }
        VariableName.clear();
      }
      
      
      /*************************************************************************
       *
       * 4.3 Case 3: Applied force is random variable
       *
       ************************************************************************/
      
      for (int ivar=0; ivar<num_ply_vars; ivar++) {

        idx_disp = 99;
        VariableName.clear(); VariableName = ply_vars_name[ivar];
        
        ostringstream first_field;
        first_field.str(""); first_field.clear();
        first_field << "DISP_MACRO" << stochastic_fields_ply[ivar];
        
        if (VariableName.compare(0,5,"force") == 0) {
          // Initiate the involved parameters to zero
          // ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // First-layer
          Ply_1st_Dmat_r(ivar).resize(6,6);   Ply_1st_Dmat_r(ivar).clear();
          // Second-layer
          Ply_2nd_Dmat_r(ivar).resize(6,6);   Ply_2nd_Dmat_r(ivar).clear();
          // Third-layer
          Ply_3rd_Dmat_r(ivar).resize(6,6);   Ply_3rd_Dmat_r(ivar).clear();
          // Fourth-layer
          Ply_4th_Dmat_r(ivar).resize(6,6);   Ply_4th_Dmat_r(ivar).clear();
          
          // Establish an object of elastic FE method
          MyElasticFEMethod_Macro my_fe_1st_Ply_r_F (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_2nd_Ply_r_F (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_3rd_Ply_r_F (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_4th_Ply_r_F (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
          
          //preproc
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          // Calculate applied forces and preassures on surface and
          // assemble force vector
          boost::ptr_map<string,NeummanForcesSurface> my_neumann_forces;
          MyMetaNeummanForces_r_PSFEM First_FE;
          ierr = First_FE.addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
          ierr = First_FE.setNeumannFiniteElementOperators(m_field_Macro,my_neumann_forces,dF,"DISP_MACRO"); CHKERRQ(ierr);
          boost::ptr_map<string,NeummanForcesSurface>::iterator mitt = my_neumann_forces.begin();
          for(;mitt!=my_neumann_forces.end();mitt++) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mitt->first,mitt->second->getLoopFe()); CHKERRQ(ierr);
          }
          
          // Assemble stiffness matrix
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_Ply_r_F);     CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_Ply_r_F);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
          // Insert value into the force vector
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          //cout<<"First order derivative of F"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          // Solve the FE equation
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",first_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_F",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          cout<<"===========================================================================\n";
          cout<<"  The first-order FE equation has been solved for the applied force, F."<<endl;
          cout<<"===========================================================================\n";
        }
        VariableName.clear();
      }
      
      
      /*************************************************************************
       *
       *  5. SOLVE THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *        [K][U_rs] = [F_rs] - [K_rs][U] - [K_r][U_s] - [K_s][U_r]  (r!=s)
       *     or
       *        [K][U_rs] = [F_rs] - [K_rs][U] - 2 [K_r][U_s]             (r==s)
       *
       ************************************************************************/
      
      Ply_1st_Dmat_rs.resize(num_ply_vars,num_ply_vars);
      Ply_2nd_Dmat_rs.resize(num_ply_vars,num_ply_vars);
      Ply_3rd_Dmat_rs.resize(num_ply_vars,num_ply_vars);
      Ply_4th_Dmat_rs.resize(num_ply_vars,num_ply_vars);
      
      for (int i=0; i<num_ply_vars; i++) {
        for (int j=0; j<num_ply_vars; j++) {
          Ply_1st_Dmat_rs(i,j).resize(6,6); Ply_1st_Dmat_rs(i,j).clear();
          Ply_2nd_Dmat_rs(i,j).resize(6,6); Ply_2nd_Dmat_rs(i,j).clear();
          Ply_3rd_Dmat_rs(i,j).resize(6,6); Ply_3rd_Dmat_rs(i,j).clear();
          Ply_4th_Dmat_rs(i,j).resize(6,6); Ply_4th_Dmat_rs(i,j).clear();
        }
      }
      
      // int PSFE_order = 2;
      if (PSFE_order == 2) {
        cout<<"\n"<<endl;
        cout<<"///////////////////////////////////////////////////////////"<<endl;
        cout<<"//"<<endl;
        cout<<"//  The solving of the 2nd macrolevel FE starts from here. "<<endl;
        cout<<"//"<<endl;
        cout<<"///////////////////////////////////////////////////////////\n"<<endl;
        // Second order
        ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
        ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        int sub_nvars = 0;
        for (int ivar = 0; ivar<num_ply_vars; ivar++) {          // num_ply_vars
          ostringstream first_field_r;
          first_field_r.str(""); first_field_r.clear();
          first_field_r << "DISP_MACRO" << stochastic_fields_ply[ivar];
          
          for (int jvar = ivar; jvar<num_ply_vars; jvar++) {     // num_ply_vars
            cout<<"The two variables are: "<<ply_vars_name[ivar]<<"\t"<<ply_vars_name[jvar]<<endl;
            
            ostringstream first_field_s;
            first_field_s.str(""); first_field_s.clear();
            first_field_s << "DISP_MACRO" << stochastic_fields_ply[jvar];
            
            // Get the second-order field
            if (jvar == ivar) {
              sub_nvars = sub_nvars + (jvar + 1);
              var_pos = sub_nvars;
            } else {
              var_pos = var_pos + jvar;
            }
            
            ostringstream second_field;
            second_field.str(""); second_field.clear();
            second_field << "DISP_MACRO" << stochastic_fields_ply[var_pos + num_ply_vars -1];
            
            if ((ivar<num_rve_vars) && (jvar<num_rve_vars)) {
              // ===============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //
              // ===============================================================
              // Construct element constitutive matrix
              // First-layer
              theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6); Ply_1st_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation(theta, RVE_Dmat_rs(ivar,jvar), Ply_1st_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Second-layer
              if (NO_Layers > 1) {
                theta = PlyAngle(1)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6); Ply_2nd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation(theta, RVE_Dmat_rs(ivar,jvar), Ply_2nd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Third-layer
              if (NO_Layers > 2) {
                theta = PlyAngle(2)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6); Ply_3rd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation(theta, RVE_Dmat_rs(ivar,jvar), Ply_3rd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                theta = PlyAngle(3)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_4th_Dmat_rs(ivar,jvar).resize(6,6); Ply_4th_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation(theta, RVE_Dmat_rs(ivar,jvar), Ply_4th_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              //ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_rs_EmEm",ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,11,"orientation")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6); Ply_1st_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_1st_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Second-layer
              if (NO_Layers > 1) {
                theta = PlyAngle(1)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6); Ply_2nd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_2nd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Third-layer
              if (NO_Layers > 2) {
                theta = PlyAngle(2)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6); Ply_3rd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_3rd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                theta = PlyAngle(3)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_4th_Dmat_rs(ivar,jvar).resize(6,6); Ply_4th_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_4th_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,6,"theta1")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_1st_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Second-layer
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              // Third-layer
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              // Fourth layer
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,6,"theta2")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              // Second-layer
              theta = PlyAngle(1)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_2nd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Third-layer
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              // Fourth layer
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,6,"theta3")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              // Second-layer
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              // Third-layer
              theta = PlyAngle(2)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_3rd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Fourth layer
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,6,"theta4")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              // Second-layer
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              // Third-layer
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              // Fourth layer
              theta = PlyAngle(3)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_r_Theta(theta, RVE_Dmat_r(ivar), Ply_4th_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ivar<num_rve_vars) && (ply_vars_name[jvar].compare(0,5,"force")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              // Second-layer
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              // Third-layer
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              // Fourth layer
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            else if ((ply_vars_name[ivar].compare(0,11,"orientation")==0) && (ply_vars_name[jvar].compare(0,11,"orientation")==0)) {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              theta = PlyAngle(0)*(M_PI/180.0);    //rotation angle about the Z-axis
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6); Ply_1st_Dmat_rs(ivar,jvar).clear();
              ierr = Dmat_Transformation_rs_Theta(theta, Dmat, Ply_1st_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              // Second-layer
              if (NO_Layers > 1) {
                theta = PlyAngle(1)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6); Ply_2nd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_rs_Theta(theta, Dmat, Ply_2nd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Third-layer
              if (NO_Layers > 2) {
                theta = PlyAngle(2)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6); Ply_3rd_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_rs_Theta(theta, Dmat, Ply_3rd_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                theta = PlyAngle(3)*(M_PI/180.0);  //rotation angle about the Z-axis
                Ply_4th_Dmat_rs(ivar,jvar).resize(6,6); Ply_4th_Dmat_rs(ivar,jvar).clear();
                ierr = Dmat_Transformation_rs_Theta(theta, Dmat, Ply_4th_Dmat_rs(ivar,jvar)); CHKERRQ(ierr);
              }
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            //else if ((ply_vars_name[ivar].compare(0,5,"force")==0) && (ply_vars_name[jvar].compare(0,5,"force")==0)) {
            else {
              // =============================================================
              //
              // The second-order derivaties w.r.t. RVE material properties
              //   and orientation
              //
              // =============================================================
              
              // First-layer
              Ply_1st_Dmat_rs(ivar,jvar).resize(6,6);   Ply_1st_Dmat_rs(ivar,jvar).clear();
              // Second-layer
              Ply_2nd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_2nd_Dmat_rs(ivar,jvar).clear();
              // Third-layer
              Ply_3rd_Dmat_rs(ivar,jvar).resize(6,6);   Ply_3rd_Dmat_rs(ivar,jvar).clear();
              // Fourth layer
              Ply_4th_Dmat_rs(ivar,jvar).resize(6,6);   Ply_4th_Dmat_rs(ivar,jvar).clear();
              
              // Assemble stiffness matrix
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_1st(m_field_Macro,A,D,ddF,Ply_1st_Dmat_r(ivar),Ply_1st_Dmat_r(jvar),"DISP_MACRO",Ply_1st_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_2nd(m_field_Macro,A,D,ddF,Ply_2nd_Dmat_r(ivar),Ply_2nd_Dmat_r(jvar),"DISP_MACRO",Ply_2nd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_3rd(m_field_Macro,A,D,ddF,Ply_3rd_Dmat_r(ivar),Ply_3rd_Dmat_r(jvar),"DISP_MACRO",Ply_3rd_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              FE2_Rhs_rr_PSFEM my_fe2_k_rs_ply_4th(m_field_Macro,A,D,ddF,Ply_4th_Dmat_r(ivar),Ply_4th_Dmat_r(jvar),"DISP_MACRO",Ply_4th_Dmat_rs(ivar,jvar),first_field_r.str().c_str(),first_field_s.str().c_str());
              
              DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs(m_field_Macro,"DISP_MACRO",A,ddD,ddF);
              
              ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Assemble stiffness matrix
              // First layer
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe2_k_rs_ply_1st);     CHKERRQ(ierr);
              // Third layer
              if (NO_Layers > 1) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe2_k_rs_ply_2nd);  CHKERRQ(ierr);
              }
              // Third layer
              if (NO_Layers > 2) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe2_k_rs_ply_3rd);  CHKERRQ(ierr);
              }
              // Fourth layer
              if (NO_Layers > 3) {
                ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe2_k_rs_ply_4th);  CHKERRQ(ierr);
              }
              ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs); CHKERRQ(ierr);
              
              // Solve finite element equation
              ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
              ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
              ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
              
              ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
              ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
              ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",second_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            }
            cout<<"===========================================================================\n";
            cout<<"  The second-order FE equation has been solved for the case of "<<ply_vars_name[ivar]<<"  "<<ply_vars_name[jvar]<<endl;
            cout<<"===========================================================================\n";
            
          } // for - jvar
        } //  for - ivar
      } // if (PSFE_order == 2)
        
      
      
      /***************************************************************************
       *
       *  6. FINISH
       *
       **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    //================
    //
    //  This version of Micro_FE_REL only works for the first-order reliability
    //    analysis.
    //
    //  History:
    //    2016.07.19: add note
    //================
    
    PetscErrorCode Macro_FE_REL(FieldInterface &m_field_Macro,
                                        int &nvars, int &nders,
                                        vector<string> &stochastic_fields,
                                        ublas::vector<double> TheVariables,
                                        int num_rvars,
                                        vector<string> vars_name,
                                        ublas::vector<double> PlyAngle,
                                        PetscInt NO_Layers) {
      PetscFunctionBegin;
      cout<<"\n Hi from macrolevel FE calculation!"<<endl;
      
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      Vec dF, dD;
      Vec ddF, ddD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&ddF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&ddD); CHKERRQ(ierr);
      
      /*****************************************************************************
       *
       *  1. Read the saved Dmat mechancial
       *     (from the computational homgenisaiton of the 0deg RVE)
       *
       ****************************************************************************/
      
      double theta;
      
      /*Dmat.clear();
       Dmat(0,0) = 1.2894e5; Dmat(0,1) = 0.0525e5; Dmat(0,2) = 0.0525e5;
       Dmat(1,0) = 0.0525e5; Dmat(1,1) = 0.1331e5; Dmat(1,2) = 0.0545e5;
       Dmat(2,0) = 0.0525e5; Dmat(2,1) = 0.0545e5; Dmat(2,2) = 0.1331e5;
       Dmat(3,3) = 0.0660e5; Dmat(4,4) = 0.0393e5; Dmat(5,5) = 0.0660e5;*/
      // First-layer
      theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
      Dmat_1st_Ply.resize(6,6);   Dmat_1st_Ply.clear();//cout<<"\nThe Dmat is: "<<Dmat<<endl;
      ierr = Dmat_Transformation(theta, Dmat, Dmat_1st_Ply); CHKERRQ(ierr);
      // Second-layer
      if (NO_Layers > 1) {
        theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_2nd_Ply.resize(6,6);   Dmat_2nd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_2nd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 2: \t"<<PlyAngle(1)<<"\n"<<Dmat_2nd_Ply<<endl;
      }
      // Third-layer
      if (NO_Layers > 2) {
        theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_3rd_Ply.resize(6,6);   Dmat_3rd_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_3rd_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 3: \t"<<PlyAngle(2)<<"\n"<<Dmat_3rd_Ply<<endl;
      }
      // Fourth-layer
      if (NO_Layers > 3) {
        theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
        Dmat_4th_Ply.resize(6,6);   Dmat_4th_Ply.clear();
        ierr = Dmat_Transformation(theta, Dmat, Dmat_4th_Ply); CHKERRQ(ierr);//cout<<"\n\nLayer 4: \t"<<PlyAngle(3)<<"\n"<<Dmat_4th_Ply<<endl;
      }
      
      /*************************************************************************
       *
       *  2. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
        MyElasticFEMethod_Macro(FieldInterface& _m_field_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
        FE2_ElasticFEMethod(_m_field_Macro,_A,_D,_F, _Dmat, _field_name) {};
        
        virtual PetscErrorCode RhsAndLhs() {
          PetscFunctionBegin;
          
          ierr = Lhs(); CHKERRQ(ierr);
          
          PetscFunctionReturn(0);
        }
      };
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);

      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      MyElasticFEMethod_Macro my_fe_1st_ply (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_2nd_ply (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_3rd_ply (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_4th_ply (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");

      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      // First layer
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_ply);  CHKERRQ(ierr);
      // Second layer
      if (NO_Layers > 1) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_ply);  CHKERRQ(ierr);
      }
      // Third layer
      if (NO_Layers > 2) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_ply);  CHKERRQ(ierr);
      }
      // Fourth layer
      if (NO_Layers > 3) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_ply);  CHKERRQ(ierr);
      }
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      
      // Check whether force is considered as random variable or not
      int idx_force = 100;
      for (int ii = 1; ii<=num_rvars; ii++) {
        string VariableName;
        VariableName = vars_name[ii];
        if (VariableName.compare(0,5,"force") == 0) {
          idx_force = ii;cout<<"The force is: "<<TheVariables(ii-1)<<endl;
        }
      }
      
      
      if (idx_force==100) {
        MetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      } else {
        MyMetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,
                                                          neumann_forces,F,
                                                          "DISP_MACRO",
                                                          "MESH_NODE_POSITIONS",
                                                          TheVariables(idx_force-1)); CHKERRQ(ierr);
      }
      
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      // elastic analys
      ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      
      cout<<"Solving the zeroth-order equation is finish. \n";
      
      
      /*************************************************************************
       *
       *  4. SOLVE THE FIRST-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *
       ************************************************************************/
      for (int ii=1; ii<=num_rvars; ii++) {
        
        /***********************************************************************
         *
         * 4.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,3,"NUz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,2,"Ef") == 0) {
          idx_disp = 7;
        }
        else if (VariableName.compare(0,3,"NUf") == 0) {
          idx_disp = 8;
        }
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 80;
        }
        else if (VariableName.compare(0,11,"orientation") == 0) {
          idx_disp = 90;
        }
        else if (VariableName.compare(0,6,"theta1") == 0) {
          idx_disp = 91;
        }
        else if (VariableName.compare(0,6,"theta2") == 0) {
          idx_disp = 92;
        }
        else if (VariableName.compare(0,6,"theta3") == 0) {
          idx_disp = 93;
        }
        else if (VariableName.compare(0,6,"theta4") == 0) {
          idx_disp = 94;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
          ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        
        switch (idx_disp) {
          case 0: {// due to Young's modulus of matrix (Em)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Em.resize(6,6);   Dmat_1st_Ply_r_Em.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_1st_Ply_r_Em); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Em.resize(6,6);   Dmat_2nd_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_2nd_Ply_r_Em); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Em.resize(6,6);   Dmat_3rd_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_3rd_Ply_r_Em); CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Em.resize(6,6);   Dmat_4th_Ply_r_Em.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_4th_Ply_r_Em); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Em,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Em,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Em_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Em_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Em_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Em_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { // due to Poisson's ratio of matrix (NUm)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUm.resize(6,6);   Dmat_1st_Ply_r_NUm.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_1st_Ply_r_NUm); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUm.resize(6,6);   Dmat_2nd_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_2nd_Ply_r_NUm); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUm.resize(6,6);   Dmat_3rd_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_3rd_Ply_r_NUm); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUm.resize(6,6);   Dmat_4th_Ply_r_NUm.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_4th_Ply_r_NUm); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUm,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUm,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUm_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUm_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUm_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUm_4th_Ply);  CHKERRQ(ierr);
            }
            
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: {// due to transversal Poisson's ratio of fibre (NUp)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUp.resize(6,6);   Dmat_1st_Ply_r_NUp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_1st_Ply_r_NUp); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUp.resize(6,6);   Dmat_2nd_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_2nd_Ply_r_NUp); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUp.resize(6,6);   Dmat_3rd_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_3rd_Ply_r_NUp); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUp.resize(6,6);   Dmat_4th_Ply_r_NUp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_4th_Ply_r_NUp); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUp_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUp_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUp_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUp_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: {// due to axial Poisson's ratio of fibre (NUpz)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUpz.resize(6,6);   Dmat_1st_Ply_r_NUpz.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_1st_Ply_r_NUpz); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUpz.resize(6,6);   Dmat_2nd_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_2nd_Ply_r_NUpz); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUpz.resize(6,6);   Dmat_3rd_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_3rd_Ply_r_NUpz); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUpz.resize(6,6);   Dmat_4th_Ply_r_NUpz.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_4th_Ply_r_NUpz); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUpz,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUpz,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUpz_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUpz_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUpz_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUpz_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 4: {// due to transversal modulus of fibre (Ep)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ep.resize(6,6);   Dmat_1st_Ply_r_Ep.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_1st_Ply_r_Ep); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ep.resize(6,6);   Dmat_2nd_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_2nd_Ply_r_Ep); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ep.resize(6,6);   Dmat_3rd_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_3rd_Ply_r_Ep); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ep.resize(6,6);   Dmat_4th_Ply_r_Ep.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_4th_Ply_r_Ep); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ep,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ep,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ep_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ep_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ep_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ep_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            break;
          }
          case 5: {// due to axial modulus of fibre (Ez)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ez.resize(6,6);   Dmat_1st_Ply_r_Ez.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_1st_Ply_r_Ez); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ez.resize(6,6);   Dmat_2nd_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_2nd_Ply_r_Ez); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ez.resize(6,6);   Dmat_3rd_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_3rd_Ply_r_Ez); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ez.resize(6,6);   Dmat_4th_Ply_r_Ez.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_4th_Ply_r_Ez); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ez,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ez,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ez_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ez_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ez_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ez_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: {// due to shear modulus of fibre (Gzp)
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Gzp.resize(6,6);   Dmat_1st_Ply_r_Gzp.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_1st_Ply_r_Gzp); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Gzp.resize(6,6);   Dmat_2nd_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_2nd_Ply_r_Gzp); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Gzp.resize(6,6);   Dmat_3rd_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_3rd_Ply_r_Gzp); CHKERRQ(ierr);
            }
            // Fourth-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Gzp.resize(6,6);   Dmat_4th_Ply_r_Gzp.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_4th_Ply_r_Gzp); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Gzp,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Gzp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Gzp_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Gzp_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Gzp_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Gzp_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
          case 7: {// due to Young's modulus of fibre (Ef) - isotropic
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_Ef.resize(6,6);   Dmat_1st_Ply_r_Ef.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_1st_Ply_r_Ef); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_Ef.resize(6,6);   Dmat_2nd_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_2nd_Ply_r_Ef); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_Ef.resize(6,6);   Dmat_3rd_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_3rd_Ply_r_Ef); CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_Ef.resize(6,6);   Dmat_4th_Ply_r_Ef.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_Ef, Dmat_4th_Ply_r_Ef); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Ef,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ef_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Ef,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ef(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Ef_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Ef_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Ef_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Ef_4th_Ply);  CHKERRQ(ierr);
            }
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ef); CHKERRQ(ierr);
            break;
          }
          case 8: { // due to Poisson's ratio of fibre (NUf) - isotropic
            // First-layer
            theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_1st_Ply_r_NUf.resize(6,6);   Dmat_1st_Ply_r_NUf.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_1st_Ply_r_NUf); CHKERRQ(ierr);
            // Second-layer
            if (NO_Layers > 1) {
              theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_2nd_Ply_r_NUf.resize(6,6);   Dmat_2nd_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_2nd_Ply_r_NUf); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 2) {
              theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_3rd_Ply_r_NUf.resize(6,6);   Dmat_3rd_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_3rd_Ply_r_NUf); CHKERRQ(ierr);
            }
            // Third-layer
            if (NO_Layers > 3) {
              theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
              Dmat_4th_Ply_r_NUf.resize(6,6);   Dmat_4th_Ply_r_NUf.clear();
              ierr = Dmat_Transformation(theta, Dmat_r_NUf, Dmat_4th_Ply_r_NUf); CHKERRQ(ierr);
            }
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_NUf,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUf_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_NUf,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUf(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
            // First layer
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_NUf_1st_Ply);  CHKERRQ(ierr);
            // Second layer
            if (NO_Layers > 1) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_NUf_2nd_Ply);  CHKERRQ(ierr);
            }
            // Third layer
            if (NO_Layers > 2) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_NUf_3rd_Ply);  CHKERRQ(ierr);
            }
            // Fourth layer
            if (NO_Layers > 3) {
              ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_NUf_4th_Ply);  CHKERRQ(ierr);
            }
            
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUf); CHKERRQ(ierr);
            break;
          }
        }
        // post-processing
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ss_field.str(""); ss_field.clear();
          ss_field << "DISP_MACRO" << stochastic_fields[idx_disp];
          if (idx_disp<nvars) {
            
            ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
            
            //cout<<"First order derivative of dD"<<endl;
            //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the first-order equation "<<ss_field.str().c_str()<<" is finish. \n";
            
            //cout<<"First order derivative of F"<<endl;
            //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
          }
          else {
            ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the second-order equation "<<ss_field.str().c_str()<<" is finish. \n";
          }
          //ierr = KSPReset(solver_Macro); CHKERRQ(ierr);
        }
        
        /***********************************************************************
         *
         * 4.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        
        if (idx_disp == 80) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // Establish an object of elastic FE method
          MyElasticFEMethod_Macro my_fe_1st_Ply_r_F (m_field_Macro,A,D,F,Dmat_1st_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_2nd_Ply_r_F (m_field_Macro,A,D,F,Dmat_2nd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_3rd_Ply_r_F (m_field_Macro,A,D,F,Dmat_3rd_Ply,"DISP_MACRO");
          MyElasticFEMethod_Macro my_fe_4th_Ply_r_F (m_field_Macro,A,D,F,Dmat_4th_Ply,"DISP_MACRO");
          
          //preproc
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          // Calculate applied forces and preassures on surface and
          // assemble force vector
          boost::ptr_map<string,NeummanForcesSurface> my_neumann_forces;
          MyMetaNeummanForces_r_PSFEM First_FE;
          ierr = First_FE.addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
          ierr = First_FE.setNeumannFiniteElementOperators(m_field_Macro,my_neumann_forces,dF,"DISP_MACRO"); CHKERRQ(ierr);
          boost::ptr_map<string,NeummanForcesSurface>::iterator mitt = my_neumann_forces.begin();
          for(;mitt!=my_neumann_forces.end();mitt++) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mitt->first,mitt->second->getLoopFe()); CHKERRQ(ierr);
          }
          
          // Assemble stiffness matrix
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply", my_fe_1st_Ply_r_F);     CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply", my_fe_2nd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply", my_fe_3rd_Ply_r_F);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply", my_fe_4th_Ply_r_F);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
          // Insert value into the force vector
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          //cout<<"First order derivative of F"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          // Solve the FE equation
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_F",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        }
        
        /***********************************************************************
         *
         * 4.3. Case 3: Ply orientation angle is treated as random variable
         *
         **********************************************************************/
        if (idx_disp == 90) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta.resize(6,6);   Dmat_1st_Ply_r_Theta.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta); CHKERRQ(ierr);
          // Second-layer
          if (NO_Layers > 1) {
            theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_2nd_Ply_r_Theta.resize(6,6);   Dmat_2nd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Third-layer
          if (NO_Layers > 2) {
            theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_3rd_Ply_r_Theta.resize(6,6);   Dmat_3rd_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta); CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
            Dmat_4th_Ply_r_Theta.resize(6,6);   Dmat_4th_Ply_r_Theta.clear();
            ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta); CHKERRQ(ierr);
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_1st_Ply(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_2nd_Ply(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_3rd_Ply(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_r_Theta_4th_Ply(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_r_Theta_1st_Ply);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_r_Theta_2nd_Ply);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_r_Theta_3rd_Ply);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_r_Theta_4th_Ply);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if (idx_disp == 91) { // first layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          //cout<<"\n\nThe angle case Theta_1 starts from here!!!"<<endl;
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          theta = PlyAngle(0)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_1st_Ply_r_Theta_1.resize(6,6);   Dmat_1st_Ply_r_Theta_1.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_1st_Ply_r_Theta_1); CHKERRQ(ierr);
          //cout<<"\n\nThe first ply angle is "<<PlyAngle(0)<<endl;
          //cout<<"\n\n Dmat_1st_Ply_r_Theta_1 is "<<Dmat_1st_Ply_r_Theta_1<<endl;
          // Second-layer
          if (NO_Layers > 1) {
            Dmat_2nd_Ply_r_Theta_1.resize(6,6);   Dmat_2nd_Ply_r_Theta_1.clear();
          }
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_1.resize(6,6);   Dmat_3rd_Ply_r_Theta_1.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_1.resize(6,6);   Dmat_4th_Ply_r_Theta_1.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_1,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_1(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_1,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_1(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_1);  CHKERRQ(ierr);
          // Second layer
          if (NO_Layers > 1) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_1);  CHKERRQ(ierr);
          }
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_1); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_1st_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_1"<<endl;
        }
        else if (idx_disp == 92) { // second layer angle is a random variable
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          //cout<<"\n\nThe angle case Theta_2 starts from here!!!"<<endl;
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_2.resize(6,6);   Dmat_1st_Ply_r_Theta_2.clear();
          // Second-layer
          theta = PlyAngle(1)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_2nd_Ply_r_Theta_2.resize(6,6);   Dmat_2nd_Ply_r_Theta_2.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_2nd_Ply_r_Theta_2); CHKERRQ(ierr);
          // Third-layer
          if (NO_Layers > 2) {
            Dmat_3rd_Ply_r_Theta_2.resize(6,6);   Dmat_3rd_Ply_r_Theta_2.clear();
          }
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_2.resize(6,6);   Dmat_4th_Ply_r_Theta_2.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_2,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_2(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_2,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_2(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_2);  CHKERRQ(ierr);
          // Third layer
          if (NO_Layers > 2) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_2);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_2); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_2nd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_2"<<endl;
        }
        else if (idx_disp == 93) { // third layer angle is a random variable
          // Initiate the involved parameters to zero
          //cout<<"\n\nThe angle case Theta_3 starts from here!!!"<<endl;
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_3.resize(6,6);   Dmat_1st_Ply_r_Theta_3.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_3.resize(6,6);   Dmat_2nd_Ply_r_Theta_3.clear();
          // Third-layer
          theta = PlyAngle(2)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_3rd_Ply_r_Theta_3.resize(6,6);   Dmat_3rd_Ply_r_Theta_3.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_3rd_Ply_r_Theta_3); CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            Dmat_4th_Ply_r_Theta_3.resize(6,6);   Dmat_4th_Ply_r_Theta_3.clear();
          }
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_3,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_3(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_3,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_3(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_3);  CHKERRQ(ierr);
          // Fourth layer
          if (NO_Layers > 3) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_3);  CHKERRQ(ierr);
          }
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_3); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_3rd_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_3"<<endl;
        }
        else if (idx_disp == 94) { // fourth layer angle is a random variable
          // Initiate the involved parameters to zero
          // cout<<"\n\nThe angle case Theta_4 starts from here!!!"<<endl;
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          // First-layer
          Dmat_1st_Ply_r_Theta_4.resize(6,6);   Dmat_1st_Ply_r_Theta_4.clear();
          // Second-layer
          Dmat_2nd_Ply_r_Theta_4.resize(6,6);   Dmat_2nd_Ply_r_Theta_4.clear();
          // Third-layer
          Dmat_3rd_Ply_r_Theta_4.resize(6,6);   Dmat_3rd_Ply_r_Theta_4.clear();
          // Fourth layer
          theta = PlyAngle(3)*(M_PI/180.0); //rotation angle about the Z-axis
          Dmat_4th_Ply_r_Theta_4.resize(6,6);   Dmat_4th_Ply_r_Theta_4.clear();
          ierr = Dmat_Transformation_r_Theta(theta, Dmat, Dmat_4th_Ply_r_Theta_4); CHKERRQ(ierr);
          
          FE2_Rhs_r_PSFEM my_fe2_k_1st_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_1st_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_2nd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_2nd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_3rd_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_3rd_Ply_r_Theta_4,"DISP_MACRO");
          FE2_Rhs_r_PSFEM my_fe2_k_4th_Ply_r_Theta_4(m_field_Macro,A,dD,dF,Dmat_4th_Ply_r_Theta_4,"DISP_MACRO");
          DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Theta_4(m_field_Macro,"DISP_MACRO",A,dD,dF);
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          // First layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_1st_Ply",my_fe2_k_1st_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Second layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_2nd_Ply",my_fe2_k_2nd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Third layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_3rd_Ply",my_fe2_k_3rd_Ply_r_Theta_4);  CHKERRQ(ierr);
          // Fourth layer
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_4th_Ply",my_fe2_k_4th_Ply_r_Theta_4);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Theta_4); CHKERRQ(ierr);
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_Theta_4th_Ply",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first order derivative of Theta_4"<<endl;
        }
        
      }
      
      /***************************************************************************
       *
       *  5. FINISH
       *
       **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation for
    //
    //        Two-Layer laminate
    //
    // =========================================================================
    
    virtual PetscErrorCode Macro_FE_Solver_Laminate(FieldInterface &m_field_Macro,
                                                    int &nvars, int &nders,
                                                    vector<string> &stochastic_fields) {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      Vec dF, dD;
      Vec ddF, ddD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&ddF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&ddD); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      //Matrix View
      //MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      //std::string wait;
      //std::cin >> wait;
      
      
      struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
        MyElasticFEMethod_Macro(FieldInterface& _m_field_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
        FE2_ElasticFEMethod(_m_field_Macro,_A,_D,_F, _Dmat, _field_name) {};
        
        virtual PetscErrorCode RhsAndLhs() {
          PetscFunctionBegin;
          
          ierr = Lhs(); CHKERRQ(ierr);
          
          PetscFunctionReturn(0);
        }
      };
      
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      
      
      /*****************************************************************************
       *
       * Read the saved Dmat mechancial (from the computational homgenisaiton of the 0deg RVE)
       *
       ****************************************************************************/
      
      //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
      double theta; theta=90*(M_PI/180.0); //rotation angle about the Z-axis
      ublas::matrix<FieldData> Dmat_90deg;  Dmat_90deg.resize(6,6);   Dmat_90deg.clear();
      ierr = Dmat_Transformation(theta, Dmat, Dmat_90deg); CHKERRQ(ierr);
      cout<<"\n\nDmat_90deg = "<<Dmat_90deg<<endl;
      
      //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
      theta=25*(M_PI/180.0); //rotation angle about the Z-axis
      ublas::matrix<FieldData> Dmat_pos25deg;  Dmat_pos25deg.resize(6,6);   Dmat_pos25deg.clear();
      ierr = Dmat_Transformation(theta, Dmat, Dmat_pos25deg); CHKERRQ(ierr);
      cout<<"\n\n Dmat_pos25deg = "<<Dmat_pos25deg<<endl;
      
      //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
      theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
      ublas::matrix<FieldData> Dmat_neg25deg;  Dmat_neg25deg.resize(6,6);   Dmat_neg25deg.clear();
      ierr = Dmat_Transformation(theta, Dmat, Dmat_neg25deg); CHKERRQ(ierr);
      cout<<"\n\n Dmat_neg25deg = "<<Dmat_neg25deg<<endl;
      
      
      MyElasticFEMethod_Macro my_fe_90deg    (m_field_Macro,A,D,F,Dmat_90deg,   "DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_pos25deg (m_field_Macro,A,D,F,Dmat_pos25deg,"DISP_MACRO");
      MyElasticFEMethod_Macro my_fe_neg25deg (m_field_Macro,A,D,F,Dmat_neg25deg,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg"   , my_fe_90deg);     CHKERRQ(ierr);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg", my_fe_pos25deg);  CHKERRQ(ierr);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg", my_fe_neg25deg);  CHKERRQ(ierr);
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      
      // Check whether force is considered as random variable or not
      
      
      MetaNeummanForces Zeroth_FE;
      ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      
      //VecCopy(ElemForce,F);
      
      
      /*************************************************************************
       *
       *  2. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      
      // elastic analys
      ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      
      cout<<"Solving the zeroth-order equation is finish. \n";
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ************************************************************************/
      for (int ii=0; ii<nvars; ii++) {
        
        if (ii < nvars) {
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((ii >= nvars) && (ii < nders)) {
          ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
          ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        
        switch (ii) {
          case 0: {// due to Young's modulus of matrix (Em)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_90deg;  Dmat_r_Em_90deg.resize(6,6);   Dmat_r_Em_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Em_90deg = "<<Dmat_r_Em_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_pos25deg;  Dmat_r_Em_pos25deg.resize(6,6);   Dmat_r_Em_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Em_pos25deg = "<<Dmat_r_Em_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_neg25deg;  Dmat_r_Em_neg25deg.resize(6,6);   Dmat_r_Em_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Em_neg25deg = "<<Dmat_r_Em_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_90deg(m_field_Macro,A,dD,dF,Dmat_r_Em_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_pos25deg(m_field_Macro,A,dD,dF,Dmat_r_Em_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em_neg25deg(m_field_Macro,A,dD,dF,Dmat_r_Em_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_Em_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_Em_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_Em_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { // due to Poisson's ratio of matrix (NUm)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_90deg;  Dmat_r_NUm_90deg.resize(6,6);   Dmat_r_NUm_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUm_90deg = "<<Dmat_r_NUm_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_pos25deg;  Dmat_r_NUm_pos25deg.resize(6,6);   Dmat_r_NUm_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUm_pos25deg = "<<Dmat_r_NUm_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_neg25deg;  Dmat_r_NUm_neg25deg.resize(6,6);   Dmat_r_NUm_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUm_neg25deg = "<<Dmat_r_NUm_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_90deg(m_field_Macro,A,D,dF,Dmat_r_NUm_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_pos25deg(m_field_Macro,A,D,dF,Dmat_r_NUm_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm_neg25deg(m_field_Macro,A,D,dF,Dmat_r_NUm_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_NUm_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_NUm_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_NUm_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: {// due to transversal Poisson's ratio of fibre (NUp)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_90deg;  Dmat_r_NUp_90deg.resize(6,6);   Dmat_r_NUp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUp_90deg = "<<Dmat_r_NUp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_pos25deg;  Dmat_r_NUp_pos25deg.resize(6,6);   Dmat_r_NUp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUp_pos25deg = "<<Dmat_r_NUp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_neg25deg;  Dmat_r_NUp_neg25deg.resize(6,6);   Dmat_r_NUp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUp_neg25deg = "<<Dmat_r_NUp_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_90deg(m_field_Macro,A,D,dF,Dmat_r_NUp_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_pos25deg(m_field_Macro,A,D,dF,Dmat_r_NUp_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp_neg25deg(m_field_Macro,A,D,dF,Dmat_r_NUp_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_NUp_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_NUp_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_NUp_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: {// due to axial Poisson's ratio of fibre (NUpz)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_90deg;  Dmat_r_NUpz_90deg.resize(6,6);   Dmat_r_NUpz_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUpz_90deg = "<<Dmat_r_NUpz_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_pos25deg;  Dmat_r_NUpz_pos25deg.resize(6,6);   Dmat_r_NUpz_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_pos25deg = "<<Dmat_r_NUpz_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_neg25deg;  Dmat_r_NUpz_neg25deg.resize(6,6);   Dmat_r_NUpz_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_neg25deg = "<<Dmat_r_NUpz_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_90deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_pos25deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_neg25deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_NUpz_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_NUpz_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_NUpz_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 4: {// due to transversal modulus of fibre (Ep)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_90deg;  Dmat_r_NUpz_90deg.resize(6,6);   Dmat_r_NUpz_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUpz_90deg = "<<Dmat_r_NUpz_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_pos25deg;  Dmat_r_NUpz_pos25deg.resize(6,6);   Dmat_r_NUpz_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_pos25deg = "<<Dmat_r_NUpz_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_neg25deg;  Dmat_r_NUpz_neg25deg.resize(6,6);   Dmat_r_NUpz_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_neg25deg = "<<Dmat_r_NUpz_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_90deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_pos25deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz_neg25deg(m_field_Macro,A,D,dF,Dmat_r_NUpz_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_NUpz_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_NUpz_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_NUpz_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 5: {// due to axial modulus of fibre (Ez)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_90deg;  Dmat_r_Ez_90deg.resize(6,6);   Dmat_r_Ez_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Ez_90deg = "<<Dmat_r_Ez_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_pos25deg;  Dmat_r_Ez_pos25deg.resize(6,6);   Dmat_r_Ez_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ez_pos25deg = "<<Dmat_r_Ez_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_neg25deg;  Dmat_r_Ez_neg25deg.resize(6,6);   Dmat_r_Ez_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ez_neg25deg = "<<Dmat_r_Ez_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_90deg(m_field_Macro,A,dD,dF,Dmat_r_Ez_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_pos25deg(m_field_Macro,A,dD,dF,Dmat_r_Ez_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez_neg25deg(m_field_Macro,A,dD,dF,Dmat_r_Ez_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_Ez_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_Ez_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_Ez_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: {// due to shear modulus of fibre (Gzp)
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_90deg;  Dmat_r_Gzp_90deg.resize(6,6);   Dmat_r_Gzp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Gzp_90deg = "<<Dmat_r_Gzp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_pos25deg;  Dmat_r_Gzp_pos25deg.resize(6,6);   Dmat_r_Gzp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Gzp_pos25deg = "<<Dmat_r_Gzp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_neg25deg;  Dmat_r_Gzp_neg25deg.resize(6,6);   Dmat_r_Gzp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Gzp_neg25deg = "<<Dmat_r_Gzp_neg25deg<<endl;
            
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_90deg(m_field_Macro,A,D,dF,Dmat_r_Gzp_90deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_pos25deg(m_field_Macro,A,D,dF,Dmat_r_Gzp_pos25deg,"DISP_MACRO");
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp_neg25deg(m_field_Macro,A,D,dF,Dmat_r_Gzp_neg25deg,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_r_Gzp_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_r_Gzp_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_r_Gzp_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
          case 7: {// 2nd order due to Em & Em
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta = 90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_90deg;  Dmat_r_Em_90deg.resize(6,6);   Dmat_r_Em_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Em_90deg = "<<Dmat_r_Em_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta = 25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_pos25deg;  Dmat_r_Em_pos25deg.resize(6,6);   Dmat_r_Em_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Em_pos25deg = "<<Dmat_r_Em_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta = -25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Em_neg25deg;  Dmat_r_Em_neg25deg.resize(6,6);   Dmat_r_Em_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Em, Dmat_r_Em_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Em_neg25deg = "<<Dmat_r_Em_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta = 90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EmEm_90deg;  Dmat_rs_EmEm_90deg.resize(6,6);   Dmat_rs_EmEm_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EmEm, Dmat_rs_EmEm_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_EmEm_90deg = "<<Dmat_rs_EmEm_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta = 25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EmEm_pos25deg;  Dmat_rs_EmEm_pos25deg.resize(6,6);   Dmat_rs_EmEm_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EmEm, Dmat_rs_EmEm_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EmEm_pos25deg = "<<Dmat_rs_EmEm_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta = -25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EmEm_neg25deg;  Dmat_rs_EmEm_neg25deg.resize(6,6);   Dmat_rs_EmEm_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EmEm, Dmat_rs_EmEm_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EmEm_neg25deg = "<<Dmat_rs_EmEm_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EmEm_90deg(m_field_Macro,A,D,ddF,Dmat_r_Em_90deg,"DISP_MACRO",Dmat_rs_EmEm_90deg,"DISP_MACRO_r_Em");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EmEm_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_Em_pos25deg,"DISP_MACRO",Dmat_rs_EmEm_pos25deg,"DISP_MACRO_r_Em");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EmEm_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_Em_neg25deg,"DISP_MACRO",Dmat_rs_EmEm_neg25deg,"DISP_MACRO_r_Em");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_EmEm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            cout<<"\nPostion 1"<<endl;
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_EmEm); CHKERRQ(ierr);
            cout<<"\nPostion 2"<<endl;
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_EmEm_90deg);  CHKERRQ(ierr);
            cout<<"\nPostion 3"<<endl;
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_EmEm_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_EmEm_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_EmEm); CHKERRQ(ierr);
            cout<<"\nTo 1\n"<<endl;
            break;
          }
          case 8: {// 2nd order due to NUm & NUm
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_90deg;  Dmat_r_NUm_90deg.resize(6,6);   Dmat_r_NUm_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUm_90deg = "<<Dmat_r_NUm_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_pos25deg;  Dmat_r_NUm_pos25deg.resize(6,6);   Dmat_r_NUm_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUm_pos25deg = "<<Dmat_r_NUm_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUm_neg25deg;  Dmat_r_NUm_neg25deg.resize(6,6);   Dmat_r_NUm_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUm, Dmat_r_NUm_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUm_neg25deg = "<<Dmat_r_NUm_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUmNUm_90deg;  Dmat_rs_NUmNUm_90deg.resize(6,6);   Dmat_rs_NUmNUm_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUmNUm, Dmat_rs_NUmNUm_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_NUmNUm_90deg = "<<Dmat_rs_NUmNUm_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUmNUm_pos25deg;  Dmat_rs_NUmNUm_pos25deg.resize(6,6);   Dmat_rs_NUmNUm_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUmNUm, Dmat_rs_NUmNUm_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUmNUm_pos25deg = "<<Dmat_rs_NUmNUm_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUmNUm_neg25deg;  Dmat_rs_NUmNUm_neg25deg.resize(6,6);   Dmat_rs_NUmNUm_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUmNUm, Dmat_rs_NUmNUm_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUmNUm_neg25deg = "<<Dmat_rs_NUmNUm_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUmNUm_90deg(m_field_Macro,A,D,ddF,Dmat_r_NUm_90deg,"DISP_MACRO",Dmat_rs_NUmNUm_90deg,"DISP_MACRO_r_NUm");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUmNUm_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_NUm_pos25deg,"DISP_MACRO",Dmat_rs_NUmNUm_pos25deg,"DISP_MACRO_r_NUm");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUmNUm_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_NUm_neg25deg,"DISP_MACRO",Dmat_rs_NUmNUm_neg25deg,"DISP_MACRO_r_NUm");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_NUmNUm_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_NUmNUm_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_NUmNUm_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            break;
          }
          case 9: {// 2nd order due to NUp & NUp
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_90deg;  Dmat_r_NUp_90deg.resize(6,6);   Dmat_r_NUp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUp_90deg = "<<Dmat_r_NUp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_pos25deg;  Dmat_r_NUp_pos25deg.resize(6,6);   Dmat_r_NUp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUp_pos25deg = "<<Dmat_r_NUp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUp_neg25deg;  Dmat_r_NUp_neg25deg.resize(6,6);   Dmat_r_NUp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUp, Dmat_r_NUp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUp_neg25deg = "<<Dmat_r_NUp_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpNUp_90deg;  Dmat_rs_NUpNUp_90deg.resize(6,6);   Dmat_rs_NUpNUp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpNUp, Dmat_rs_NUpNUp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_NUpNUp_90deg = "<<Dmat_rs_NUpNUp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpNUp_pos25deg;  Dmat_rs_NUpNUp_pos25deg.resize(6,6);   Dmat_rs_NUpNUp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpNUp, Dmat_rs_NUpNUp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUpNUp_pos25deg = "<<Dmat_rs_NUpNUp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpNUp_neg25deg;  Dmat_rs_NUpNUp_neg25deg.resize(6,6);   Dmat_rs_NUpNUp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpNUp, Dmat_rs_NUpNUp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUpNUp_neg25deg = "<<Dmat_rs_NUpNUp_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpNUp_90deg(m_field_Macro,A,D,ddF,Dmat_r_NUp_90deg,"DISP_MACRO",Dmat_rs_NUpNUp_90deg,"DISP_MACRO_r_NUp");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpNUp_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_NUp_pos25deg,"DISP_MACRO",Dmat_rs_NUpNUp_pos25deg,"DISP_MACRO_r_NUp");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpNUp_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_NUp_neg25deg,"DISP_MACRO",Dmat_rs_NUpNUp_neg25deg,"DISP_MACRO_r_NUp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_NUpNUp_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_NUpNUp_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_NUpNUp_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            break;
          }
          case 10: {// 2nd order due to NUpz & NUpz
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_90deg;  Dmat_r_NUpz_90deg.resize(6,6);   Dmat_r_NUpz_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_NUpz_90deg = "<<Dmat_r_NUpz_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_pos25deg;  Dmat_r_NUpz_pos25deg.resize(6,6);   Dmat_r_NUpz_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_pos25deg = "<<Dmat_r_NUpz_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_NUpz_neg25deg;  Dmat_r_NUpz_neg25deg.resize(6,6);   Dmat_r_NUpz_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_NUpz, Dmat_r_NUpz_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_NUpz_neg25deg = "<<Dmat_r_NUpz_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpzNUpz_90deg;  Dmat_rs_NUpzNUpz_90deg.resize(6,6);   Dmat_rs_NUpzNUpz_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpzNUpz, Dmat_rs_NUpzNUpz_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_NUpzNUpz_90deg = "<<Dmat_rs_NUpzNUpz_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpzNUpz_pos25deg;  Dmat_rs_NUpzNUpz_pos25deg.resize(6,6);   Dmat_rs_NUpzNUpz_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpzNUpz, Dmat_rs_NUpzNUpz_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUpzNUpz_pos25deg = "<<Dmat_rs_NUpzNUpz_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_NUpzNUpz_neg25deg;  Dmat_rs_NUpzNUpz_neg25deg.resize(6,6);   Dmat_rs_NUpzNUpz_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_NUpzNUpz, Dmat_rs_NUpzNUpz_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_NUpzNUpz_neg25deg = "<<Dmat_rs_NUpzNUpz_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpzNUpz_90deg(m_field_Macro,A,D,ddF,Dmat_r_NUpz_90deg,"DISP_MACRO",Dmat_rs_NUpzNUpz_90deg,"DISP_MACRO_r_NUpz");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpzNUpz_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_NUpz_pos25deg,"DISP_MACRO",Dmat_rs_NUpzNUpz_pos25deg,"DISP_MACRO_r_NUpz");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpzNUpz_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_NUpz_neg25deg,"DISP_MACRO",Dmat_rs_NUpzNUpz_neg25deg,"DISP_MACRO_r_NUpz");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_NUpzNUpz_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_NUpzNUpz_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_NUpzNUpz_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            break;
          }
          case 11: {// 2nd order due to Ep & Ep
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ep_90deg;  Dmat_r_Ep_90deg.resize(6,6);   Dmat_r_Ep_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_r_Ep_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Ep_90deg = "<<Dmat_r_Ep_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ep_pos25deg;  Dmat_r_Ep_pos25deg.resize(6,6);   Dmat_r_Ep_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_r_Ep_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ep_pos25deg = "<<Dmat_r_Ep_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ep_neg25deg;  Dmat_r_Ep_neg25deg.resize(6,6);   Dmat_r_Ep_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ep, Dmat_r_Ep_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ep_neg25deg = "<<Dmat_r_Ep_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EpEp_90deg;  Dmat_rs_EpEp_90deg.resize(6,6);   Dmat_rs_EpEp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EpEp, Dmat_rs_EpEp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_EpEp_90deg = "<<Dmat_rs_EpEp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EpEp_pos25deg;  Dmat_rs_EpEp_pos25deg.resize(6,6);   Dmat_rs_EpEp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EpEp, Dmat_rs_EpEp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EpEp_pos25deg = "<<Dmat_rs_EpEp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EpEp_neg25deg;  Dmat_rs_EpEp_neg25deg.resize(6,6);   Dmat_rs_EpEp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EpEp, Dmat_rs_EpEp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EpEp_neg25deg = "<<Dmat_rs_EpEp_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EpEp_90deg(m_field_Macro,A,D,ddF,Dmat_r_Ep_90deg,"DISP_MACRO",Dmat_rs_EpEp_90deg,"DISP_MACRO_r_Ep");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EpEp_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_Ep_pos25deg,"DISP_MACRO",Dmat_rs_EpEp_pos25deg,"DISP_MACRO_r_Ep");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EpEp_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_Ep_neg25deg,"DISP_MACRO",Dmat_rs_EpEp_neg25deg,"DISP_MACRO_r_Ep");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_EpEp_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_EpEp_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_EpEp_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            break;
          }
          case 12: {// 2nd order due to Ez & Ez
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_90deg;  Dmat_r_Ez_90deg.resize(6,6);   Dmat_r_Ez_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Ez_90deg = "<<Dmat_r_Ez_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_pos25deg;  Dmat_r_Ez_pos25deg.resize(6,6);   Dmat_r_Ez_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ez_pos25deg = "<<Dmat_r_Ez_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Ez_neg25deg;  Dmat_r_Ez_neg25deg.resize(6,6);   Dmat_r_Ez_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Ez, Dmat_r_Ez_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Ez_neg25deg = "<<Dmat_r_Ez_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EzEz_90deg;  Dmat_rs_EzEz_90deg.resize(6,6);   Dmat_rs_EzEz_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EzEz, Dmat_rs_EzEz_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_EzEz_90deg = "<<Dmat_rs_EzEz_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EzEz_pos25deg;  Dmat_rs_EzEz_pos25deg.resize(6,6);   Dmat_rs_EzEz_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EzEz, Dmat_rs_EzEz_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EzEz_pos25deg = "<<Dmat_rs_EzEz_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_EzEz_neg25deg;  Dmat_rs_EzEz_neg25deg.resize(6,6);   Dmat_rs_EzEz_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_EzEz, Dmat_rs_EzEz_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_EzEz_neg25deg = "<<Dmat_rs_EzEz_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EzEz_90deg(m_field_Macro,A,D,ddF,Dmat_r_Ez_90deg,"DISP_MACRO",Dmat_rs_EzEz_90deg,"DISP_MACRO_r_Ez");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EzEz_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_Ez_pos25deg,"DISP_MACRO",Dmat_rs_EzEz_pos25deg,"DISP_MACRO_r_Ez");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EzEz_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_Ez_neg25deg,"DISP_MACRO",Dmat_rs_EzEz_neg25deg,"DISP_MACRO_r_Ez");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_EzEz_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_EzEz_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_EzEz_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            break;
          }
          case 13: {// 2nd order due to Gzp & Gzp
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_90deg;  Dmat_r_Gzp_90deg.resize(6,6);   Dmat_r_Gzp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_r_Gzp_90deg = "<<Dmat_r_Gzp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_pos25deg;  Dmat_r_Gzp_pos25deg.resize(6,6);   Dmat_r_Gzp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Gzp_pos25deg = "<<Dmat_r_Gzp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_r_Gzp_neg25deg;  Dmat_r_Gzp_neg25deg.resize(6,6);   Dmat_r_Gzp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_r_Gzp, Dmat_r_Gzp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_r_Gzp_neg25deg = "<<Dmat_r_Gzp_neg25deg<<endl;
            
            //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
            theta=90*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_GzpGzp_90deg;  Dmat_rs_GzpGzp_90deg.resize(6,6);   Dmat_rs_GzpGzp_90deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_GzpGzp, Dmat_rs_GzpGzp_90deg); CHKERRQ(ierr);
            cout<<"\n\nDmat_rs_GzpGzp_90deg = "<<Dmat_rs_GzpGzp_90deg<<endl;
            
            //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_GzpGzp_pos25deg;  Dmat_rs_GzpGzp_pos25deg.resize(6,6);   Dmat_rs_GzpGzp_pos25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_GzpGzp, Dmat_rs_GzpGzp_pos25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_GzpGzp_pos25deg = "<<Dmat_rs_GzpGzp_pos25deg<<endl;
            
            //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
            theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
            ublas::matrix<FieldData> Dmat_rs_GzpGzp_neg25deg;  Dmat_rs_GzpGzp_neg25deg.resize(6,6);   Dmat_rs_GzpGzp_neg25deg.clear();
            ierr = Dmat_Transformation(theta, Dmat_rs_GzpGzp, Dmat_rs_GzpGzp_neg25deg); CHKERRQ(ierr);
            cout<<"\n\n Dmat_rs_GzpGzp_neg25deg = "<<Dmat_rs_GzpGzp_neg25deg<<endl;
            
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_GzpGzp_90deg(m_field_Macro,A,D,ddF,Dmat_r_Gzp_90deg,"DISP_MACRO",Dmat_rs_GzpGzp_90deg,"DISP_MACRO_r_Gzp");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_GzpGzp_pos25deg(m_field_Macro,A,D,ddF,Dmat_r_Gzp_pos25deg,"DISP_MACRO",Dmat_rs_GzpGzp_pos25deg,"DISP_MACRO_r_Gzp");
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_GzpGzp_neg25deg(m_field_Macro,A,D,ddF,Dmat_r_Gzp_neg25deg,"DISP_MACRO",Dmat_rs_GzpGzp_neg25deg,"DISP_MACRO_r_Gzp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_90deg",my_fe2_k_rs_GzpGzp_90deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_pos25deg",my_fe2_k_rs_GzpGzp_pos25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_neg25deg",my_fe2_k_rs_GzpGzp_neg25deg);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            break;
          }
        }
        // post-processing
        ostringstream ss_field;
        ss_field.str(""); ss_field.clear();
        ss_field << "DISP_MACRO" << stochastic_fields[ii];
        cout<<"The field is "<<ss_field.str().c_str()<<endl;
        if (ii<nvars) {
          
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          
          //cout<<"First order derivative of dD"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the first-order equation "<<ss_field.str().c_str()<<" is finish. \n";
          
          //cout<<"First order derivative of F"<<endl;
          //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
        }
        else {
          ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          cout<<"Solving the second-order equation "<<ss_field.str().c_str()<<" is finish. \n";
        }
        
      }
      
      /***************************************************************************
       *
       *  4. FINISH
       *
       **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    
    virtual PetscErrorCode Calculate_RVEDmat_PSFE(FieldInterface &m_field_RVE,
                                                  int &nvars, int &nders,
                                                  vector<string> &stochastic_fields) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);   Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);     Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;cout<<"\n Young modulus: "<<YoungModulus<<endl;
      double PoissonRatio = 0.3;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume = 0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE" ,MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
      VecGetArray(Stress_Homo, &avec);
      for (int ii=0; ii<6; ii++) {
        Dmat(ii,0)=*avec;
        avec++;
      }
      
      /*if (pcomm->rank()==0) {
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,0)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++){
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,0) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,0) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,0) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,0) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,0) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,0) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,0) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,0) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,0) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,0) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
        
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      // solve for F2 and D2
      ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,1)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,1)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,1) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,1) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,1) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,1) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,1) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,1) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,1) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,1) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,1) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,1) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F3 and D3
      ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,2)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,2)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,2) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,2) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,2) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,2) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,2) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,2) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,2) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,2) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,2) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,2) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F4 and D4
      ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,3)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,3)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,3) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,3) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,3) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,3) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,3) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,3) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,3) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,3) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,3) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,3) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F5 and D5
      ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,4)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,4)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,4) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,4) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,4) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,4) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,4) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,4) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,4) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,4) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,4) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,4) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F6 and D6
      ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,5)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,5)<<endl;
       }
       cout<< "\n\n";
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,5) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,5) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,5) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,5) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,5) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,5) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,5) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,5) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,5) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,5) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      ierr = VecDestroy(&F1); CHKERRQ(ierr);
      ierr = VecDestroy(&F2); CHKERRQ(ierr);
      ierr = VecDestroy(&F3); CHKERRQ(ierr);
      ierr = VecDestroy(&F4); CHKERRQ(ierr);
      ierr = VecDestroy(&F5); CHKERRQ(ierr);
      ierr = VecDestroy(&F6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
  };
  
  /*****************************************************************************
   *                                                                           *
   *                                                                           *
   *                         SINGLE LAYER LAMINA                               *
   *                                                                           *
   *                                                                           *
   ****************************************************************************/
  
  struct FE2_Macro_Solver {
    
    //global variable Dmat
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r_Em;
    ublas::matrix<double> Dmat_r_NUm;
    ublas::matrix<double> Dmat_r_Ep;
    ublas::matrix<double> Dmat_r_Ez;
    ublas::matrix<double> Dmat_r_NUp;
    ublas::matrix<double> Dmat_r_NUpz;
    ublas::matrix<double> Dmat_r_Gzp;
    
    ublas::matrix<double> Dmat_rs_EmEm;
    ublas::matrix<double> Dmat_rs_NUmNUm;
    ublas::matrix<double> Dmat_rs_EpEp;
    ublas::matrix<double> Dmat_rs_EzEz;
    ublas::matrix<double> Dmat_rs_NUpNUp;
    ublas::matrix<double> Dmat_rs_NUpzNUpz;
    ublas::matrix<double> Dmat_rs_GzpGzp;
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    
    virtual PetscErrorCode Calculate_RVEDmat(FieldInterface &m_field_RVE,
                                     int &nvars, int &nders,
                                     vector<string> &stochastic_fields,
                                     ublas::vector<double> matprop,
                                     int num_rvars,
                                     vector<string> vars_name) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);  Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);   Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;
      double PoissonRatio = 0.3;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field_RVE,BLOCKSET,it)) {
        cout << endl << *it << endl;
        
        //Get block name
        string name = it->get_name();
        
        if (name.compare(0,13,"MAT_ELASTIC_1") == 0) {
          Mat_Elastic mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int ii=1;ii<=num_rvars;ii++) {
            ParameterName = vars_name[ii];
            
            if (ParameterName.compare(0,2,"Em") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Young = matprop(ii-1);
            }
            else if (ParameterName.compare(0,3,"NUm") == 0) {cout<<"the variable name is "<<ParameterName<<endl;
              mydata.data.Poisson = matprop(ii-1);
            }
            ParameterName.clear();
          }
        
          /*
           mydata.data.Young   = matprop(0);
           mydata.data.Poisson = matprop(1);
           */
          YoungModulus=mydata.data.Young;
          PoissonRatio=mydata.data.Poisson;
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Matrix material:\n" << mydata;
        }
        
        if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0) {
          Mat_Elastic_TransIso mydata;
          ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
          
          string ParameterName;
          for (int i=1;i<=num_rvars;i++) {
            ParameterName = vars_name[i];
            cout<<"the variable name is "<<vars_name[i]<<endl;
            if (ParameterName.compare(0,2,"Ez") == 0) {
              mydata.data.Youngz = matprop(i-1);
            }
            else if (ParameterName.compare(0,2,"Ep") == 0) {
              mydata.data.Youngp = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"NUp") == 0) {
              mydata.data.Poissonp = matprop(i-1);
            }
            else if (ParameterName.compare(0,4,"NUpz") == 0) {
              mydata.data.Poissonpz = matprop(i-1);
            }
            else if (ParameterName.compare(0,3,"Gzp") == 0) {
              mydata.data.Shearzp = matprop(i-1);
            }
            ParameterName.clear();
          }

          /*
          mydata.data.Poissonp  = matprop(2);
          mydata.data.Poissonpz = matprop(3);
          mydata.data.Youngp    = matprop(4);
          mydata.data.Youngz    = matprop(5);
          mydata.data.Shearzp   = matprop(6);
           */
          ierr = it->set_attribute_data_structure(mydata); CHKERRQ(ierr);
          cout << "Fibre material:\n" << mydata;
        }
      }
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
      VecGetArray(Stress_Homo, &avec);
      for (int ii=0; ii<6; ii++) {
        Dmat(ii,0)=*avec;
        avec++;
      }
      
      /*if (pcomm->rank()==0) {
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,0)<<endl;
        }
      }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }

        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
                 
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++){
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,0) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,0) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,0) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,0) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,0) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,0) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,0) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,0) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,0) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,0) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      // solve for F2 and D2
      ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,1)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,1)<<endl;
        }
      }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }

        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,1) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,1) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,1) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,1) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,1) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,1) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,1) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,1) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,1) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,1) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F3 and D3
      ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,2)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,2)<<endl;
        }
      }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
          
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,2) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,2) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,2) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,2) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,2) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,2) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,2) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,2) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,2) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,2) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F4 and D4
      ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,3)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,3)<<endl;
        }
      }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,3) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,3) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,3) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,3) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,3) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,3) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,3) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,3) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,3) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,3) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F5 and D5
      ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,4)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,4)<<endl;
        }
      }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,4) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,4) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,4) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,4) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,4) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,4) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,4) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,4) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,4) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,4) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F6 and D6
      ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,5)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
        cout<< "\nStress_Homo = \n\n";
        for(int ii=0; ii<6; ii++){
          cout <<Dmat(ii,5)<<endl;
        }
        cout<< "\n\n";
      }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=1; ii <= num_rvars; ii++) {
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (idx_disp == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (idx_disp == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (idx_disp == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (idx_disp == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (idx_disp == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (idx_disp == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (idx_disp == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) { // solution for first-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) { // solution for second-order problem
          ss_field << "DISP_RVE" << stochastic_fields[idx_disp];
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        if ((idx_disp>=0) && (idx_disp<nders)) {
          ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
          
          ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
          
          if(pcomm->rank()==0) {
            PetscScalar    *avec_r;
            VecGetArray(Stress_Homo_r, &avec_r);
            
            //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
            //cout<< "\n"<<ss_field<<" = \n\n";
            for(int irow=0; irow<6; irow++) {
              cout.precision(15);
              //cout<<*avec_r<<endl;
              switch (idx_disp) {
                case 0:
                  Dmat_r_Em(irow,5) = *avec_r; break;
                case 1:
                  Dmat_r_NUm(irow,5) = *avec_r; break;
                case 2:
                  Dmat_r_NUp(irow,5) = *avec_r; break;
                case 3:
                  Dmat_r_NUpz(irow,5) = *avec_r; break;
                case 4:
                  Dmat_r_Ep(irow,5) = *avec_r; break;
                case 5:
                  Dmat_r_Ez(irow,5) = *avec_r; break;
                case 6:
                  Dmat_r_Gzp(irow,5) = *avec_r; break;
                case 7:
                  Dmat_rs_EmEm(irow,5) = *avec_r; break;
                case 8:
                  Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
                case 9:
                  Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
                case 10:
                  Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
                case 11:
                  Dmat_rs_EpEp(irow,5) = *avec_r; break;
                case 12:
                  Dmat_rs_EzEz(irow,5) = *avec_r; break;
                case 13:
                  Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
              }
              
              // write result to output file
              //TheFile<<setprecision(15)<<*avec_r<<'\n';
              avec_r++;
            }
            VecRestoreArray(Stress_Homo_r, &avec_r);
          }
        }
        //cout<< "\n\n";
      }
      
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      ierr = VecDestroy(&F1); CHKERRQ(ierr);
      ierr = VecDestroy(&F2); CHKERRQ(ierr);
      ierr = VecDestroy(&F3); CHKERRQ(ierr);
      ierr = VecDestroy(&F4); CHKERRQ(ierr);
      ierr = VecDestroy(&F5); CHKERRQ(ierr);
      ierr = VecDestroy(&F6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    // =========================================================================
    //
    //  B.VI. SOLUTION PHASE:
    //        Solve Macroscale FE equation
    //
    // =========================================================================
    
    virtual PetscErrorCode Macro_FE_Solver(FieldInterface &m_field_Macro,
                                           int &nvars, int &nders,
                                           vector<string> &stochastic_fields,
                                           ublas::vector<double> TheVariables,
                                           int num_rvars,
                                           vector<string> vars_name) {
      PetscFunctionBegin;
      
      PetscErrorCode ierr;
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      //create matrices
      
      Vec F, D;
      Vec dF, dD;
      Vec ddF, ddD;
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&F); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&D); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&dF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&dD); CHKERRQ(ierr);
      
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,&ddF); CHKERRQ(ierr);
      ierr = m_field_Macro.VecCreateGhost("ELASTIC_PROBLEM_MACRO",COL,&ddD); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat A;
      ierr = m_field_Macro.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_MACRO",&A); CHKERRQ(ierr);
      
      //Matrix View
      //MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      //std::string wait;
      //std::cin >> wait;
      
      
      struct MyElasticFEMethod_Macro: public FE2_ElasticFEMethod {
        MyElasticFEMethod_Macro(FieldInterface& _m_field_Macro,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
        FE2_ElasticFEMethod(_m_field_Macro,_A,_D,_F, _Dmat, _field_name) {};
        
        virtual PetscErrorCode RhsAndLhs() {
          PetscFunctionBegin;
          
          ierr = Lhs(); CHKERRQ(ierr);
          
          PetscFunctionReturn(0);
        }
      };
      
      
      Projection10NodeCoordsOnField ent_method_material_Macro(m_field_Macro,"MESH_NODE_POSITIONS");
      ierr = m_field_Macro.loop_dofs("MESH_NODE_POSITIONS",ent_method_material_Macro); CHKERRQ(ierr);
      
      //Assemble F and A
      DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field_Macro,"DISP_MACRO",A,D,F);
      MyElasticFEMethod_Macro my_fe_Macro(m_field_Macro,A,D,F,Dmat,"DISP_MACRO");
      
      ierr = VecZeroEntries(F); CHKERRQ(ierr);
      ierr = VecZeroEntries(D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = MatZeroEntries(A); CHKERRQ(ierr);
      
      ierr = m_field_Macro.set_global_VecCreateGhost("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //preproc
      ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      //loop elems
      //PetscBarrier(PETSC_NULL);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe_Macro);  CHKERRQ(ierr);
      ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe_Macro);  CHKERRQ(ierr);
      
      //forces and preassures on surface
      boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
      // Check whether force is considered as random variable or not
      int idx_force = 100;
      for (int ii = 1; ii<=num_rvars; ii++) {
        string VariableName;
        VariableName = vars_name[ii];
        if (VariableName.compare(0,5,"force") == 0) {
          idx_force = ii;cout<<"The force is: "<<TheVariables(ii-1)<<endl;
        }
      }
      
      
      if (idx_force==100) {
        MetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,neumann_forces,F,"DISP_MACRO"); CHKERRQ(ierr);
      } else {
        MyMetaNeummanForces Zeroth_FE;
        ierr = Zeroth_FE.setNeumannFiniteElementOperators(m_field_Macro,
                                                                    neumann_forces,F,
                                                                    "DISP_MACRO",
                                                                    "MESH_NODE_POSITIONS",
                                                                    TheVariables(idx_force-1)); CHKERRQ(ierr);
      }
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
      }
      
      //postproc
      ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
      
      //set matrix possitives define and symetric for cholesky and icc preceonditionser
      ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
      
      ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
      
      //VecCopy(ElemForce,F);
      
      
      /*************************************************************************
       *
       *  2. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      //Solver
      KSP solver_Macro;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver_Macro,A,A); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver_Macro); CHKERRQ(ierr);
      ierr = KSPSetUp(solver_Macro); CHKERRQ(ierr);
      
      //MatView(A,PETSC_VIEWER_STDOUT_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
      
      // elastic analys
      ierr = KSPSolve(solver_Macro,F,D); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      
      //Save data on mesh
      ierr = m_field_Macro.set_global_ghost_vector("ELASTIC_PROBLEM_MACRO",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      //VecView(D,PETSC_VIEWER_STDOUT_WORLD);
      //VecView(F,PETSC_VIEWER_STDOUT_WORLD);
      
      cout<<"Solving the zeroth-order equation is finish. \n";
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FIRST-ORDER AND THE SECOND-ORDER FE EQUILIBRIUM EQUATION
       *     1st order-[K][U_r] = -[K_r][U}
       *     2nd order-[K][U_rs] = -[K_rs][U]-2[K_r][U_s]
       *
       ************************************************************************/
      for (int ii=1; ii<=num_rvars; ii++) {
        
        /***********************************************************************
         *
         * 3.1. Case 1: Material properties are  treated as random variables
         *
         **********************************************************************/
        
        int idx_disp = 99; // index of displacement field in the field name vector <stochastic_fields>
        string VariableName = vars_name[ii];
        if (VariableName.compare(0,2,"Em") == 0) {
          idx_disp = 0;
        }
        else if (VariableName.compare(0,3,"NUm") == 0) {
          idx_disp = 1;
        }
        else if (VariableName.compare(0,3,"NUp") == 0) {
          idx_disp = 2;
        }
        else if (VariableName.compare(0,4,"NUpz") == 0) {
          idx_disp = 3;
        }
        else if (VariableName.compare(0,2,"Ep") == 0) {
          idx_disp = 4;
        }
        else if (VariableName.compare(0,2,"Ez") == 0) {
          idx_disp = 5;
        }
        else if (VariableName.compare(0,3,"Gzp") == 0) {
          idx_disp = 6;
        }
        else if (VariableName.compare(0,5,"force") == 0) {
          idx_disp = 80;
        }
        VariableName.clear();
        
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        else if ((idx_disp >= nvars) && (idx_disp < nders)) {
          ierr = VecZeroEntries(ddD); CHKERRQ(ierr);
          ierr = VecZeroEntries(ddF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        
        switch (idx_disp) {
          case 0: {// due to Young's modulus of matrix (Em)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Em(m_field_Macro,A,dD,dF,Dmat_r_Em,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Em);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Em);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Em); CHKERRQ(ierr);
            break;
          }
          case 1: { // due to Poisson's ratio of matrix (NUm)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUm(m_field_Macro,A,D,dF,Dmat_r_NUm,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUm); CHKERRQ(ierr);
            break;
          }
          case 2: {// due to transversal Poisson's ratio of fibre (NUp)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUp(m_field_Macro,A,D,dF,Dmat_r_NUp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUp); CHKERRQ(ierr);
            break;
          }
          case 3: {// due to axial Poisson's ratio of fibre (NUpz)
            FE2_Rhs_r_PSFEM my_fe2_k_r_NUpz(m_field_Macro,A,D,dF,Dmat_r_NUpz,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_NUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_NUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_NUpz); CHKERRQ(ierr);
            break;
          }
          case 4: {// due to transversal modulus of fibre (Ep)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ep(m_field_Macro,A,D,dF,Dmat_r_Ep,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Ep);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Ep);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ep); CHKERRQ(ierr);
            break;
          }
          case 5: {// due to axial modulus of fibre (Ez)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Ez(m_field_Macro,A,dD,dF,Dmat_r_Ez,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Ez);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Ez);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Ez); CHKERRQ(ierr);
            break;
          }
          case 6: {// due to shear modulus of fibre (Gzp)
            FE2_Rhs_r_PSFEM my_fe2_k_r_Gzp(m_field_Macro,A,D,dF,Dmat_r_Gzp,"DISP_MACRO");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_r_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_r_Gzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_r_Gzp); CHKERRQ(ierr);
            break;
          }
          case 7: {// 2nd order due to Em & Em
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EmEm(m_field_Macro,A,D,ddF,Dmat_r_Em,"DISP_MACRO",Dmat_rs_EmEm,"DISP_MACRO_r_Em");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Em(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Em); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EmEm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EmEm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Em); CHKERRQ(ierr);
            break;
          }
          case 8: {// 2nd order due to NUm & NUm
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUmNUm(m_field_Macro,A,D,ddF,Dmat_r_NUm,"DISP_MACRO",Dmat_rs_NUmNUm,"DISP_MACRO_r_NUm");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUm(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUmNUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUmNUm);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUm); CHKERRQ(ierr);
            break;
          }
          case 9: {// 2nd order due to NUp & NUp
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpNUp(m_field_Macro,A,D,ddF,Dmat_r_NUp,"DISP_MACRO",Dmat_rs_NUpNUp,"DISP_MACRO_r_NUp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUpNUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUpNUp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUp); CHKERRQ(ierr);
            break;
          }
          case 10: {// 2nd order due to NUpz & NUpz
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_NUpzNUpz(m_field_Macro,A,D,ddF,Dmat_r_NUpz,"DISP_MACRO",Dmat_rs_NUpzNUpz,"DISP_MACRO_r_NUpz");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_NUpz(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_NUpzNUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_NUpzNUpz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_NUpz); CHKERRQ(ierr);
            break;
          }
          case 11: {// 2nd order due to Ep & Ep
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EpEp(m_field_Macro,A,D,ddF,Dmat_r_Ep,"DISP_MACRO",Dmat_rs_EpEp,"DISP_MACRO_r_Ep");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ep(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EpEp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EpEp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ep); CHKERRQ(ierr);
            break;
          }
          case 12: {// 2nd order due to Ez & Ez
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_EzEz(m_field_Macro,A,D,ddF,Dmat_r_Ez,"DISP_MACRO",Dmat_rs_EzEz,"DISP_MACRO_r_Ez");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Ez(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_EzEz);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_EzEz);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Ez); CHKERRQ(ierr);
            break;
          }
          case 13: {// 2nd order due to Gzp & Gzp
            FE2_Rhs_rs_PSFEM my_fe2_k_rs_GzpGzp(m_field_Macro,A,D,ddF,Dmat_r_Gzp,"DISP_MACRO",Dmat_rs_GzpGzp,"DISP_MACRO_r_Gzp");
            DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc_rs_Gzp(m_field_Macro,"DISP_MACRO",A,dD,dF);
            ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe2_k_rs_GzpGzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe2_k_rs_GzpGzp);  CHKERRQ(ierr);
            ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc_rs_Gzp); CHKERRQ(ierr);
            break;
          }
        }
        // post-processing
        ostringstream ss_field;
        if ((idx_disp >=0) && (idx_disp < nvars)) {
          ss_field.str(""); ss_field.clear();
          ss_field << "DISP_MACRO" << stochastic_fields[idx_disp];
          if (idx_disp<nvars) {
            
            ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
            
            //cout<<"First order derivative of dD"<<endl;
            //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);//ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the first-order equation "<<ss_field.str().c_str()<<" is finish. \n";
            
            //cout<<"First order derivative of F"<<endl;
            //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
            
          }
          else {
            ierr = VecGhostUpdateBegin(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            ierr = VecAssemblyBegin(ddF); CHKERRQ(ierr);
            ierr = VecAssemblyEnd(ddF); CHKERRQ(ierr);
            
            ierr = KSPSolve(solver_Macro,ddF,ddD); CHKERRQ(ierr);
            ierr = VecGhostUpdateBegin(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGhostUpdateEnd(ddD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO",ss_field.str().c_str(),ROW,ddD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
            cout<<"Solving the second-order equation "<<ss_field.str().c_str()<<" is finish. \n";
          }
          //ierr = KSPReset(solver_Macro); CHKERRQ(ierr);
        }
        
        /***********************************************************************
         *
         * 3.2. Case 2: Applied forces are treated as random variables
         *
         **********************************************************************/
        
        if (idx_disp == 80) {
          // Initiate the involved parameters to zero
          //ierr = VecZeroEntries(D); CHKERRQ(ierr);
          ierr = VecZeroEntries(dF); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = VecZeroEntries(dD); CHKERRQ(ierr);
          ierr = MatZeroEntries(A); CHKERRQ(ierr);
          
          // Establish an object of elastic FE method
          MyElasticFEMethod_Macro my_fe_Macro_r_F(m_field_Macro,A,dD,dF,Dmat,"DISP_MACRO");
          
          //preproc
          ierr = m_field_Macro.problem_basic_method_preProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          // Calculate applied forces and preassures on surface and
          // assemble force vector
          boost::ptr_map<string,NeummanForcesSurface> my_neumann_forces;
          MyMetaNeummanForces_r_PSFEM First_FE;
          ierr = First_FE.addNeumannBCElements(m_field_Macro,"DISP_MACRO"); CHKERRQ(ierr);
          ierr = First_FE.setNeumannFiniteElementOperators(m_field_Macro,my_neumann_forces,dF,"DISP_MACRO"); CHKERRQ(ierr);
          boost::ptr_map<string,NeummanForcesSurface>::iterator mitt = my_neumann_forces.begin();
          for(;mitt!=my_neumann_forces.end();mitt++) {
            ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO",mitt->first,mitt->second->getLoopFe()); CHKERRQ(ierr);
          }
          
          // Assemble stiffness matrix
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO",my_fe_Macro_r_F);  CHKERRQ(ierr);
          ierr = m_field_Macro.loop_finite_elements("ELASTIC_PROBLEM_MACRO","ELASTIC_FE_MACRO_REL",my_fe_Macro_r_F);  CHKERRQ(ierr);
          
          ierr = m_field_Macro.problem_basic_method_postProcess("ELASTIC_PROBLEM_MACRO",my_dirichlet_bc); CHKERRQ(ierr);
          
          ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
          
          // Insert value into the force vector
          ierr = VecGhostUpdateBegin(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF); CHKERRQ(ierr);
          //cout<<"First order derivative of F"<<endl;
          //ierr = VecView(dF,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
          
          // Solve the FE equation
          ierr = KSPSolve(solver_Macro,dF,dD); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_Macro.set_other_global_ghost_vector("ELASTIC_PROBLEM_MACRO","DISP_MACRO","DISP_MACRO_r_F",ROW,dD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          //ierr = VecView(dD,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        }
        
      }
      
    /***************************************************************************
     *
     *  4. FINISH
     *
     **************************************************************************/
      
      //Destroy matrices
      ierr = VecDestroy(&F); CHKERRQ(ierr);
      ierr = VecDestroy(&dF); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF); CHKERRQ(ierr);
      ierr = VecDestroy(&D); CHKERRQ(ierr);
      ierr = VecDestroy(&dD); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD); CHKERRQ(ierr);
      ierr = MatDestroy(&A); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver_Macro); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
    // =========================================================================
    //
    //  A.VI. SOLUTION PHASE:
    //        Caculate RVE constitutive matrix Dmat
    //
    // =========================================================================
    
    virtual PetscErrorCode Calculate_RVEDmat_PSFE(FieldInterface &m_field_RVE,
                                             int &nvars, int &nders,
                                             vector<string> &stochastic_fields) {
      PetscFunctionBegin;
      cout <<"Hi from Calculate_RVEDmat"<<endl;
      
      PetscErrorCode ierr;
      
      Dmat.resize(6,6); Dmat.clear();
      
      Dmat_r_Em.resize(6,6);   Dmat_r_Em.clear();
      Dmat_r_NUm.resize(6,6);   Dmat_r_NUm.clear();
      Dmat_r_Ep.resize(6,6);   Dmat_r_Ep.clear();
      Dmat_r_Ez.resize(6,6);   Dmat_r_Ez.clear();
      Dmat_r_NUp.resize(6,6);  Dmat_r_NUp.clear();
      Dmat_r_NUpz.resize(6,6); Dmat_r_NUpz.clear();
      Dmat_r_Gzp.resize(6,6);  Dmat_r_Gzp.clear();
      
      Dmat_rs_EmEm.resize(6,6);     Dmat_rs_EmEm.clear();
      Dmat_rs_NUmNUm.resize(6,6);     Dmat_rs_NUmNUm.clear();
      Dmat_rs_EpEp.resize(6,6);     Dmat_rs_EpEp.clear();
      Dmat_rs_EzEz.resize(6,6);     Dmat_rs_EzEz.clear();
      Dmat_rs_NUpNUp.resize(6,6);   Dmat_rs_NUpNUp.clear();
      Dmat_rs_NUpzNUpz.resize(6,6); Dmat_rs_NUpzNUpz.clear();
      Dmat_rs_GzpGzp.resize(6,6);   Dmat_rs_GzpGzp.clear();
      
      
      /*************************************************************************
       *
       *  0. PREPARATION FOR PROCESSING SOLVE
       *
       ************************************************************************/
      Vec F1,F2,F3,F4,F5,F6,D1,D2,D3,D4,D5,D6;
      Vec dF1,dF2,dF3,dF4,dF5,dF6;
      Vec ddF1,ddF2,ddF3,ddF4,ddF5,ddF6;
      Vec dD1,dD2,dD3,dD4,dD5,dD6;
      Vec ddD1,ddD2,ddD3,ddD4,ddD5,ddD6;
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&F6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&dF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,&ddF6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&dD6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD1); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD2); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD3); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD4); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD5); CHKERRQ(ierr);
      ierr = m_field_RVE.VecCreateGhost("ELASTIC_PROBLEM_RVE",COL,&ddD6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  1. Assembling global stiffness matrix K
       *     and external force vector F
       ************************************************************************/
      Mat Aij;
      ierr = m_field_RVE.MatCreateMPIAIJWithArrays("ELASTIC_PROBLEM_RVE",&Aij); CHKERRQ(ierr);
      
      struct MyElasticFEMethod: public ElasticFEMethod {
        MyElasticFEMethod(FieldInterface& _m_field,
                          Mat& _Aij,Vec& _D,Vec& _F,double _lambda,double _mu, string _field_name = "DISPLACEMENT"):
        ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu,_field_name) {};
        
        PetscErrorCode Fint(Vec F_int) {
          PetscFunctionBegin;
          ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
          for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
          }
          PetscFunctionReturn(0);
        }
      };
      
      //Assemble F and Aij
      double YoungModulus = 3500;cout<<"\n Young modulus: "<<YoungModulus<<endl;
      double PoissonRatio = 0.3;
      double alpha;
      int field_rank=3;
      
      /*************************************************************************
       *
       *  2. Get the volume of RVE
       *
       ************************************************************************/
      double RVE_volume;    RVE_volume = 0.0;  //RVE volume for full RVE We need this for stress calculation
      Vec RVE_volume_Vec;
      ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field_RVE.get_moab(),MYPCOMM_INDEX);
      ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
      ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);
      
      RVEVolume MyRVEVol(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
      RVEVolumeTrans MyRVEVolTrans(m_field_RVE,Aij,D1,F1, RVE_volume_Vec);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE" ,MyRVEVol);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyRVEVolTrans);  CHKERRQ(ierr);
      //    ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
      ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
      cout<<"Final RVE_volume = "<< RVE_volume <<endl;
      
      
      ublas::vector<FieldData> applied_strain;  //it is not used in the calculation, it is required by ElasticFE_RVELagrange_Disp as input
      applied_strain.resize(1.5*field_rank+1.5); applied_strain.clear();
      
      MyElasticFEMethod my_fe(m_field_RVE,Aij,D1,F1,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),"DISP_RVE");
      TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field_RVE,Aij,D1,F1,"DISP_RVE");
      ElasticFE_RVELagrange_Disp_Multi_Rhs MyFE_RVELagrange(m_field_RVE,Aij,D1,F1,F2,F3,F4,F5,F6,applied_strain,"DISP_RVE","Lagrange_mul_disp",field_rank);
      
      cout<<"After ElasticFE_RVELagrange_Disp_Multi_Rhs = "<<endl;
      ierr = VecZeroEntries(F1); CHKERRQ(ierr);
      ierr = VecZeroEntries(F2); CHKERRQ(ierr);
      ierr = VecZeroEntries(F3); CHKERRQ(ierr);
      ierr = VecZeroEntries(F4); CHKERRQ(ierr);
      ierr = VecZeroEntries(F5); CHKERRQ(ierr);
      ierr = VecZeroEntries(F6); CHKERRQ(ierr);
      
      ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
      
      ierr = VecZeroEntries(D1); CHKERRQ(ierr);
      ierr = VecZeroEntries(D2); CHKERRQ(ierr);
      ierr = VecZeroEntries(D3); CHKERRQ(ierr);
      ierr = VecZeroEntries(D4); CHKERRQ(ierr);
      ierr = VecZeroEntries(D5); CHKERRQ(ierr);
      ierr = VecZeroEntries(D6); CHKERRQ(ierr);
      
      ierr = m_field_RVE.set_global_VecCreateGhost("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",MyTIsotFE);  CHKERRQ(ierr);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVELagrange);  CHKERRQ(ierr);
      
      
      ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F1); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F1); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F2); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F2); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F3); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F3); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F4); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F4); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F5); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F5); CHKERRQ(ierr);
      
      ierr = VecAssemblyBegin(F6); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(F6); CHKERRQ(ierr);
      
      
      /*************************************************************************
       *
       *  3. SOLVE THE FINITE ELEMENT EQUILIBRIUM EQUATION
       *     [K][U] = [F]
       *
       ************************************************************************/
      
      //Solver
      KSP solver;
      ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
      ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
      ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
      ierr = KSPSetUp(solver); CHKERRQ(ierr);
      
      //create a vector for 6 components of homogenized stress
      Vec Stress_Homo, Stress_Homo_r, Stress_Homo_rs;
      PetscScalar *avec;
      
      if(pcomm->rank()==0) {
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,6,6,0,PETSC_NULL,&Stress_Homo_rs);
      } else {
        int ghost[] = {0,1,2,3,4,5};
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_r);
        VecCreateGhost(PETSC_COMM_WORLD,0,6,6,ghost,&Stress_Homo_rs);
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 1: applied macro strain: [1 0 0 0 0 0]^T
      //------------------------------------------------------------------------
      
      cout<<"===============================================================\n";
      cout<<"        Applied strain [1 0 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      //solve for F1 and D1
      ierr = KSPSolve(solver,F1,D1); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // calculate homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_1(m_field_RVE,Aij,D1,F1,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_1);  CHKERRQ(ierr);
      VecGetArray(Stress_Homo, &avec);
      for (int ii=0; ii<6; ii++) {
        Dmat(ii,0)=*avec;
        avec++;
      }
      
      /*if (pcomm->rank()==0) {
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,0)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D1,dF1,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D1,dF1,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D1,dF1,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D1,dF1,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D1,ddF1,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF1,dD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF1,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF1); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF1); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF1,ddD1); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD1,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD1,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r1(m_field_RVE,Aij,dD1,dF1,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r1);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++){
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,0) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,0) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,0) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,0) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,0) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,0) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,0) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,0) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,0) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,0) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,0) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,0) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,0) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,0) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
        
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 2: applied macro strain: [0 1 0 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 1 0 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      
      // solve for F2 and D2
      ierr = KSPSolve(solver,F2,D2); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_2(m_field_RVE,Aij,D2,F2,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_2);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,1)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,1)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo,&avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D2,dF2,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D2,dF2,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D2,dF2,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D2,dF2,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D2,ddF2,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF2,dD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF2,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF2); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF2); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF2,ddD2); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD2,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD2,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r2(m_field_RVE,Aij,dD2,dF2,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r2);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,1) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,1) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,1) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,1) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,1) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,1) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,1) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,1) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,1) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,1) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,1) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,1) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,1) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,1) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 3: applied macro strain: [0 0 1 0 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 1 0 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F3 and D3
      ierr = KSPSolve(solver,F3,D3); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_3(m_field_RVE,Aij,D3,F3,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_3);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,2)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,2)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D3,dF3,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D3,dF3,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D3,dF3,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D3,dF3,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D3,ddF3,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF3,dD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF3,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF3); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF3); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF3,ddD3); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD3,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD3,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r3(m_field_RVE,Aij,dD3,dF3,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r3);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,2) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,2) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,2) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,2) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,2) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,2) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,2) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,2) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,2) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,2) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,2) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,2) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,2) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,2) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 4: applied macro strain: [0 0 0 1 0 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 1 0 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F4 and D4
      ierr = KSPSolve(solver,F4,D4); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      
      // Extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_4(m_field_RVE,Aij,D4,F4,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_4);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,3)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,3)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D4,dF4,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D4,dF4,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D4,dF4,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D4,dF4,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D4,ddF4,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF4,dD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF4,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF4); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF4); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF4,ddD4); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD4,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD4,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r4(m_field_RVE,Aij,dD4,dF4,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r4);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,3) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,3) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,3) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,3) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,3) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,3) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,3) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,3) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,3) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,3) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,3) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,3) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,3) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,3) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 5: applied macro strain: [0 0 0 0 1 0]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"        Applied strain [0 0 0 0 1 0]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F5 and D5
      ierr = KSPSolve(solver,F5,D5); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_5(m_field_RVE,Aij,D5,F5,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_5);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,4)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,4)<<endl;
       }
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D5,dF5,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D5,dF5,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D5,dF5,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D5,dF5,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D5,ddF5,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF5,dD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF5,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF5); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF5); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF5,ddD5); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD5,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD5,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r5(m_field_RVE,Aij,dD5,dF5,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r5);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,4) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,4) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,4) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,4) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,4) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,4) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,4) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,4) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,4) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,4) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,4) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,4) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,4) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,4) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      //------------------------------------------------------------------------
      // 3.1 Solving the equation to get nodal displacement
      //     case 6: applied macro strain: [0 0 0 0 0 1]^T
      //------------------------------------------------------------------------
      cout<<"===============================================================\n";
      cout<<"         Applied strain [0 0 0 0 0 1]^T\n";
      cout<<"===============================================================\n";
      
      //-------------
      // ZEROTH-ORDER
      //-------------
      //solve for F6 and D6
      ierr = KSPSolve(solver,F6,D6); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(D6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = m_field_RVE.set_global_ghost_vector("ELASTIC_PROBLEM_RVE",ROW,D6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      
      // Calculate and extract homogenized stress
      ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);
      ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_6(m_field_RVE,Aij,D6,F6,&RVE_volume, applied_strain, Stress_Homo,"DISP_RVE","Lagrange_mul_disp",field_rank);
      ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_6);  CHKERRQ(ierr);
      
      VecGetArray(Stress_Homo, &avec);
      for(int ii=0; ii<6; ii++){
        Dmat(ii,5)=*avec;
        avec++;
      }
      
      /*if(pcomm->rank()==0){
       cout<< "\nStress_Homo = \n\n";
       for(int ii=0; ii<6; ii++){
       cout <<Dmat(ii,5)<<endl;
       }
       cout<< "\n\n";
       }*/
      VecRestoreArray(Stress_Homo, &avec);
      
      
      //-------------
      // FIRST & SECOND-ORDER
      //-------------
      for(int ii=0; ii < nders; ii++) {
        if (ii < nvars) {
          ierr = VecZeroEntries(dF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        } else {
          ierr = VecZeroEntries(ddF6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        
        if (ii == 0) { // due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Em(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Young","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_Em);  CHKERRQ(ierr);
        }
        else if (ii == 1) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUm(m_field_RVE,Aij,D6,dF6,"DISP_RVE","Poisson","isotropic","matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_r_NUm);  CHKERRQ(ierr);
        }
        else if (ii == 2) { // due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonP","transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUp);  CHKERRQ(ierr);
        }
        else if (ii == 3) { // due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_NUpz(m_field_RVE,Aij,D6,dF6,"DISP_RVE","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_NUpz);  CHKERRQ(ierr);
        }
        else if (ii == 4) { // due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ep(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ep);  CHKERRQ(ierr);
        }
        else if (ii == 5) { // due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Ez(m_field_RVE,Aij,D6,dF6,"DISP_RVE","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Ez);  CHKERRQ(ierr);
        }
        else if (ii == 6) { // due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_r_PSFEM my_fe_k_r_Gzp(m_field_RVE,Aij,D6,dF6,"DISP_RVE","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_r_Gzp);  CHKERRQ(ierr);
        }
        else if (ii == 7) { // 2nd order derivative due to Young's modulus of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EmEm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Em","Young","Young", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_EmEm);  CHKERRQ(ierr);
        }
        else if (ii == 8) { // 2nd order derivative due to Poisson's ratio of matrix - isotropic
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUmNUm(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUm","Poisson","Poisson", "isotropic", "matrix");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","ELASTIC_FE_RVE",my_fe_k_rs_NUmNUm);  CHKERRQ(ierr);
        }
        else if (ii == 9) { // 2nd order derivative due to Poisson's ratio in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpNUp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUp","PoissonP","PoissonP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_NUpNUp);  CHKERRQ(ierr);
        }
        else if (ii == 10) { // 2nd order derivative due to Poisson's ratio in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_NUpzNUpz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_NUpz","PoissonZ","PoissonZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE", "TRAN_ISO_FE_RVE", my_fe_k_rs_NUpzNUpz);  CHKERRQ(ierr);
        }
        else if (ii == 11) { // 2nd order derivative due to Young's modulus in p-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EpEp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ep","YoungP","YoungP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EpEp);  CHKERRQ(ierr);
        }
        else if (ii == 12) { // 2nd order derivative due to Young's modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_EzEz(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Ez","YoungZ","YoungZ", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_EzEz);  CHKERRQ(ierr);
        }
        else if (ii == 13) { // 2nd order derivative due to shear modulus in z-direction of fibre
          Trans_Iso_Rhs_rs_PSFEM my_fe_k_rs_GzpGzp(m_field_RVE,Aij,D6,ddF6,"DISP_RVE","DISP_RVE_r_Gzp","ShearZP","ShearZP", "transversely_isotropic", "reinforcement");
          ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","TRAN_ISO_FE_RVE",my_fe_k_rs_GzpGzp);  CHKERRQ(ierr);
        }
        
        //
        ostringstream ss_field;
        ss_field << "DISP_RVE" << stochastic_fields[ii];
        if (ii < nvars){ // solution for first-order problem
          //------------
          // a. Solving
          //------------
          ierr = VecGhostUpdateBegin(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,dF6,dD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(dD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,dD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        else { // solution for second-order problem
          ierr = VecGhostUpdateBegin(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddF6,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(ddF6); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(ddF6); CHKERRQ(ierr);
          
          ierr = KSPSolve(solver,ddF6,ddD6); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(ddD6,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","DISP_RVE",ss_field.str().c_str(),ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = m_field_RVE.set_other_global_ghost_vector("ELASTIC_PROBLEM_RVE","Lagrange_mul_disp","Lagrange_mul_disp",ROW,ddD6,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
        }
        
        //------------
        // b. Calculating first-order homogenized stress
        //------------
        ierr = VecZeroEntries(Stress_Homo_r); CHKERRQ(ierr);
        
        ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp_r6(m_field_RVE,Aij,dD6,dF6,&RVE_volume, applied_strain, Stress_Homo_r,"DISP_RVE","Lagrange_mul_disp",field_rank);
        ierr = m_field_RVE.loop_finite_elements("ELASTIC_PROBLEM_RVE","Lagrange_FE",MyFE_RVEHomoStressDisp_r6);  CHKERRQ(ierr);
        
        if(pcomm->rank()==0) {
          PetscScalar    *avec_r;
          VecGetArray(Stress_Homo_r, &avec_r);
          
          //cout<<"\n"<<ss_field.str().c_str()<<" =\n"<<endl;
          //cout<< "\n"<<ss_field<<" = \n\n";
          for(int irow=0; irow<6; irow++) {
            cout.precision(15);
            //cout<<*avec_r<<endl;
            switch (ii) {
              case 0:
                Dmat_r_Em(irow,5) = *avec_r; break;
              case 1:
                Dmat_r_NUm(irow,5) = *avec_r; break;
              case 2:
                Dmat_r_NUp(irow,5) = *avec_r; break;
              case 3:
                Dmat_r_NUpz(irow,5) = *avec_r; break;
              case 4:
                Dmat_r_Ep(irow,5) = *avec_r; break;
              case 5:
                Dmat_r_Ez(irow,5) = *avec_r; break;
              case 6:
                Dmat_r_Gzp(irow,5) = *avec_r; break;
              case 7:
                Dmat_rs_EmEm(irow,5) = *avec_r; break;
              case 8:
                Dmat_rs_NUmNUm(irow,5) = *avec_r; break;
              case 9:
                Dmat_rs_NUpNUp(irow,5) = *avec_r; break;
              case 10:
                Dmat_rs_NUpzNUpz(irow,5) = *avec_r; break;
              case 11:
                Dmat_rs_EpEp(irow,5) = *avec_r; break;
              case 12:
                Dmat_rs_EzEz(irow,5) = *avec_r; break;
              case 13:
                Dmat_rs_GzpGzp(irow,5) = *avec_r; break;
            }
            
            // write result to output file
            //TheFile<<setprecision(15)<<*avec_r<<'\n';
            avec_r++;
          }
          VecRestoreArray(Stress_Homo_r, &avec_r);
        }
        cout<< "\n\n";
      }
      
      
      /*************************************************************************
       *
       *  4. FINISH
       *
       ************************************************************************/
      
      //Destroy matrices and vectors
      ierr = VecDestroy(&F1); CHKERRQ(ierr);
      ierr = VecDestroy(&F2); CHKERRQ(ierr);
      ierr = VecDestroy(&F3); CHKERRQ(ierr);
      ierr = VecDestroy(&F4); CHKERRQ(ierr);
      ierr = VecDestroy(&F5); CHKERRQ(ierr);
      ierr = VecDestroy(&F6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&dF1); CHKERRQ(ierr);
      ierr = VecDestroy(&dF2); CHKERRQ(ierr);
      ierr = VecDestroy(&dF3); CHKERRQ(ierr);
      ierr = VecDestroy(&dF4); CHKERRQ(ierr);
      ierr = VecDestroy(&dF5); CHKERRQ(ierr);
      ierr = VecDestroy(&dF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddF1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddF6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&D1); CHKERRQ(ierr);
      ierr = VecDestroy(&D2); CHKERRQ(ierr);
      ierr = VecDestroy(&D3); CHKERRQ(ierr);
      ierr = VecDestroy(&D4); CHKERRQ(ierr);
      ierr = VecDestroy(&D5); CHKERRQ(ierr);
      ierr = VecDestroy(&D6); CHKERRQ(ierr);
      
      
      ierr = VecDestroy(&dD1); CHKERRQ(ierr);
      ierr = VecDestroy(&dD2); CHKERRQ(ierr);
      ierr = VecDestroy(&dD3); CHKERRQ(ierr);
      ierr = VecDestroy(&dD4); CHKERRQ(ierr);
      ierr = VecDestroy(&dD5); CHKERRQ(ierr);
      ierr = VecDestroy(&dD6); CHKERRQ(ierr);
      
      ierr = VecDestroy(&ddD1); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD2); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD3); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD4); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD5); CHKERRQ(ierr);
      ierr = VecDestroy(&ddD6); CHKERRQ(ierr);
      
      ierr = MatDestroy(&Aij); CHKERRQ(ierr);
      ierr = KSPDestroy(&solver); CHKERRQ(ierr);
      
      PetscFunctionReturn(0);
    }
    
    
  };
  
}

#endif //__FE2_MACRO_SOLVER_HPP