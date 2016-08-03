/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 * FIXME: DESCRIPTION
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

#ifndef __FE2_POSTPROCSTRESSFORRELIABILITY_HPP__
#define __FE2_POSTPROCSTRESSFORRELIABILITY_HPP__

namespace ObosleteUsersModules {
  
  struct FE2_PostProcStressForReliability_Zeroth: public PostProcDisplacemenysAndStarinOnRefMesh {
    
    FieldInterface& mField;
    
    ublas::matrix<double> Dmat;
    ublas::matrix<double> StressGP;
    ublas::matrix<double> StrainGP;
    
    Tag th_stress,th_prin_stress_vect1,th_prin_stress_vect2,th_prin_stress_vect3,th_prin_stress_vals;
    
    FE2_PostProcStressForReliability_Zeroth(FieldInterface& _mField,
                                            string _field_name,
                                            ublas::matrix<double> _Dmat):
    PostProcDisplacemenysAndStarinOnRefMesh(_mField.get_moab(),_field_name),mField(_mField),Dmat(_Dmat){
      double def_VAL2[3] = { 0.0, 0.0, 0.0 };
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT1",3,MB_TYPE_DOUBLE,th_prin_stress_vect1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT2",3,MB_TYPE_DOUBLE,th_prin_stress_vect2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT3",3,MB_TYPE_DOUBLE,th_prin_stress_vect3,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VALUES",3,MB_TYPE_DOUBLE,th_prin_stress_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
      rval = moab_post_proc.tag_get_handle("STRESS_VAL",9,MB_TYPE_DOUBLE,th_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      
    }
    
    
    vector< ublas::matrix<FieldData> > invH;
    vector< FieldData > detH;
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      cout<<"Post-processing for zeroth order"<<endl;
      
      try {//EntityHandle fe_ent = VolumeElementForcesAndSourcesCore(mField)::UserDataOperator()::getMoFEMFEPtr()->get_ent();
        
        //Loop over elements
        ierr = do_operator(); CHKERRQ(ierr);
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Strains to Nodes in PostProc Mesh: create vector containing matrices
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        //du/dx
        ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
        
        int gg = 0;
        for(;viit!=GradU_at_GaussPt.end();viit++,mit++,gg++) {
          
          
          ublas::matrix< FieldData > GradU = *viit;
          if(!invH.empty()) {
            //GradU =
            //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
            //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
            //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
            //H =
            //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
            //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
            //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
            //invH =
            //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
            //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
            //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
            //GradU =
            //[ dU/dX1 dU/dX2 dU/dX3 ]
            //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
            //[ dW/dX1 dW/dX2 dW/dX3 ]
            GradU = prod( GradU, invH[gg] );
          }
          ublas::matrix< FieldData > Strain; Strain.resize(3,3); Strain = 0.5*( GradU + trans(GradU) );
          rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
          //caluate stress and save it into tag
          ublas::vector<FieldData> Strain_VectorNotation(6);
          Strain_VectorNotation[0] = Strain(0,0);
          Strain_VectorNotation[1] = Strain(1,1);
          Strain_VectorNotation[2] = Strain(2,2);
          Strain_VectorNotation[3] = 2*Strain(0,1);
          Strain_VectorNotation[4] = 2*Strain(1,2);
          Strain_VectorNotation[5] = 2*Strain(2,0);//cout<<"Strain: "<<Strain_VectorNotation<<endl;
          
          //double w = V*G_TET_W[gg];
          ublas::vector< FieldData > Stress_VectorNotation = prod( Dmat, Strain_VectorNotation );
          
          ublas::matrix< FieldData > Stress = ublas::zero_matrix<FieldData>(3,3);//cout<<"Stress: "<<Stress_VectorNotation<<endl;
          Stress(0,0) = Stress_VectorNotation[0];
          Stress(1,1) = Stress_VectorNotation[1];
          Stress(2,2) = Stress_VectorNotation[2];
          Stress(0,1) = Stress(1,0) = Stress_VectorNotation[3];
          Stress(1,2) = Stress(2,1) = Stress_VectorNotation[4];
          Stress(2,0) = Stress(0,2) = Stress_VectorNotation[5];
          
          rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(Stress.data()[0])); CHKERR_PETSC(rval);
          
          //Principal stress vectors
          //Calculated from finding the eigenvalues and eigenvectors
          ublas::matrix< FieldData > eigen_vectors = Stress;
          ublas::vector<double> eigen_values(3);
          
          //LAPACK - eigenvalues and vectors. Applied twice for initial creates memory space
          int n = 3, lda = 3, info, lwork = -1;
          double wkopt;
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),&wkopt,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          lwork = (int)wkopt;
          double work[lwork];
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),work,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          
          //Combine eigenvalues and vectors to create principal stress vector
          ublas::vector<double> prin_stress_vect1(3);
          ublas::vector<double> prin_stress_vect2(3);
          ublas::vector<double> prin_stress_vect3(3);
          ublas::vector<double> prin_vals_vect(3);
          
          //eigen_vectors = trans(eigen_vectors);
          for (int ii=0; ii < 3; ii++) {
            prin_vals_vect[0] = eigen_values[0];
            prin_vals_vect[1] = eigen_values[1];
            prin_vals_vect[2] = eigen_values[2];
            prin_stress_vect1[ii] = eigen_vectors.data()[ii+3*0];
            prin_stress_vect2[ii] = eigen_vectors.data()[ii+3*1];
            prin_stress_vect3[ii] = eigen_vectors.data()[ii+3*2];
          }
          //cout<<"D matrix: "<<Dmat<<endl;
          //cout<<"Cauchy stress at GP "<<gg<<" : \t"<<Stress<<endl;
          //cout<<"Principal stress at GP "<<gg<<" : \t"<<prin_vals_vect<<endl;
          if (gg==0) {
            cout.precision(15);
            // Output stresses at Gauss point
            StressGP.resize(3,3); StressGP.clear();
            StressGP = Stress;
            // Output straines at Gauss point
            StrainGP.resize(3,3); StrainGP.clear();
            StrainGP = Strain;
            //cout<<"Cauchy stress at GP "<<gg<<" : \t"<<Stress<<endl;
          }
          //Tag principle stress vectors 1, 2, 3
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect1,&mit->second,1,&prin_stress_vect1[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect2,&mit->second,1,&prin_stress_vect2[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect3,&mit->second,1,&prin_stress_vect3[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vals,&mit->second,1,&prin_vals_vect[0]); CHKERR_PETSC(rval);
          
        }
        
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      
      PetscFunctionReturn(0);
    }
    
  };
  
  // ===========================================================================
  //
  // Calculationg the first order derivative of stress
  //
  // ===========================================================================
  
  struct FE2_PostProcStressForReliability_First: public PostProcDisplacemenysAndStarinOnRefMesh {
    
    FieldInterface& mField;
    
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r;
    ublas::matrix<double> StressGP_r;
    ublas::matrix<double> StrainGP_r;
    string zeroth_field;
    string first_field;
    
    Tag th_stress,th_prin_stress_vect1,th_prin_stress_vect2,th_prin_stress_vect3,th_prin_stress_vals;
    
    FE2_PostProcStressForReliability_First(FieldInterface& _mField,
                                     string _zeroth_field,
                                     string _first_field,
                                     ublas::matrix<double> _Dmat,
                                     ublas::matrix<double> _Dmat_r):
    PostProcDisplacemenysAndStarinOnRefMesh(_mField.get_moab(),_zeroth_field),
                                            mField(_mField),
                                            zeroth_field(_zeroth_field),
                                            first_field(_first_field),
                                            Dmat(_Dmat), Dmat_r(_Dmat_r){
      double def_VAL2[3] = { 0.0, 0.0, 0.0 };
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT1",3,MB_TYPE_DOUBLE,th_prin_stress_vect1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT2",3,MB_TYPE_DOUBLE,th_prin_stress_vect2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT3",3,MB_TYPE_DOUBLE,th_prin_stress_vect3,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VALUES",3,MB_TYPE_DOUBLE,th_prin_stress_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
      rval = moab_post_proc.tag_get_handle("STRESS_VAL",9,MB_TYPE_DOUBLE,th_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      
    }
    
    
    vector< ublas::matrix<FieldData> > invH;
    vector< FieldData > detH;
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      try {
        
        //Loop over elements
        ierr = do_operator(); CHKERRQ(ierr);
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Strains to Nodes in PostProc Mesh: create vector containing matrices
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        vector< ublas::matrix< FieldData > > Grad_dU_at_GaussPt;
        //du/dx
        ierr = GetGaussDiffDataVector(zeroth_field,GradU_at_GaussPt); CHKERRQ(ierr);
        ierr = GetGaussDiffDataVector(first_field,Grad_dU_at_GaussPt); CHKERRQ(ierr);
        
        vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
        vector< ublas::matrix< FieldData > >::iterator viit_du = Grad_dU_at_GaussPt.begin();
        
        map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
        int gg = 0;
        
        // calculate zeroth-order displacement iduced strain/stress
        for(;viit!=GradU_at_GaussPt.end();viit++,viit_du++,mit++,gg++) {
          
          
          ublas::matrix< FieldData > GradU = *viit;
          ublas::matrix< FieldData > Grad_dU = *viit_du;
          
          if(!invH.empty()) {
            //GradU =
            //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
            //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
            //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
            //H =
            //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
            //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
            //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
            //invH =
            //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
            //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
            //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
            //GradU =
            //[ dU/dX1 dU/dX2 dU/dX3 ]
            //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
            //[ dW/dX1 dW/dX2 dW/dX3 ]
            GradU   = prod( GradU,   invH[gg] );
            Grad_dU = prod( Grad_dU, invH[gg] );
          }
          ublas::matrix< FieldData > Strain; Strain.resize(3,3); Strain = 0.5*( GradU   + trans(GradU)   );
          ublas::matrix< FieldData > Strain_du; Strain_du.resize(3,3); Strain_du = 0.5*( Grad_dU + trans(Grad_dU) );
          
          //rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
          //rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
          
          //caluate stress and save it into tag
          ublas::vector<FieldData> Strain_VectorNotation(6);
          Strain_VectorNotation[0] = Strain(0,0);
          Strain_VectorNotation[1] = Strain(1,1);
          Strain_VectorNotation[2] = Strain(2,2);
          Strain_VectorNotation[3] = 2*Strain(0,1);
          Strain_VectorNotation[4] = 2*Strain(1,2);
          Strain_VectorNotation[5] = 2*Strain(2,0);
          
          ublas::vector<FieldData> Strain_du_VectorNotation(6);
          Strain_du_VectorNotation[0] = Strain_du(0,0);
          Strain_du_VectorNotation[1] = Strain_du(1,1);
          Strain_du_VectorNotation[2] = Strain_du(2,2);
          Strain_du_VectorNotation[3] = 2*Strain_du(0,1);
          Strain_du_VectorNotation[4] = 2*Strain_du(1,2);
          Strain_du_VectorNotation[5] = 2*Strain_du(2,0);
          
          // dS/dx = [dD/dx] u + D [du/dx]
          ublas::vector< FieldData > dStress_dx_VectorNotation(6); dStress_dx_VectorNotation.clear();
          dStress_dx_VectorNotation = prod( Dmat, Strain_du_VectorNotation ) + prod( Dmat_r, Strain_VectorNotation );
          
          ublas::matrix< FieldData > dStress_dx = ublas::zero_matrix<FieldData>(3,3);
          dStress_dx(0,0) = dStress_dx_VectorNotation[0];
          dStress_dx(1,1) = dStress_dx_VectorNotation[1];
          dStress_dx(2,2) = dStress_dx_VectorNotation[2];
          dStress_dx(0,1) = dStress_dx(1,0) = dStress_dx_VectorNotation[3];
          dStress_dx(1,2) = dStress_dx(2,1) = dStress_dx_VectorNotation[4];
          dStress_dx(2,0) = dStress_dx(0,2) = dStress_dx_VectorNotation[5];
          //cout<<"The d stress: "<<dStress_dx<<"\n"<<prod( Dmat, Strain_du_VectorNotation ) << "\n"<< prod( Dmat_r, Strain_VectorNotation )<<endl;
          
          rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(dStress_dx.data()[0])); CHKERR_PETSC(rval);
          
          //Principal stress vectors
          //Calculated from finding the eigenvalues and eigenvectors
          ublas::matrix< FieldData > eigen_vectors = dStress_dx;
          ublas::vector<double> eigen_values(3);
          
          //LAPACK - eigenvalues and vectors. Applied twice for initial creates memory space
          int n = 3, lda = 3, info, lwork = -1;
          double wkopt;
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),&wkopt,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          lwork = (int)wkopt;
          double work[lwork];
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),work,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          
          //Combine eigenvalues and vectors to create principal stress vector
          ublas::vector<double> dprin_stress_dx_vect1(3);
          ublas::vector<double> dprin_stress_dx_vect2(3);
          ublas::vector<double> dprin_stress_dx_vect3(3);
          ublas::vector<double> dprin_vals_vect(3);
          
          //eigen_vectors = trans(eigen_vectors);
          for (int ii=0; ii < 3; ii++) {
            dprin_vals_vect[0] = eigen_values[0];
            dprin_vals_vect[1] = eigen_values[1];
            dprin_vals_vect[2] = eigen_values[2];
            dprin_stress_dx_vect1[ii] = eigen_vectors.data()[ii+3*0];
            dprin_stress_dx_vect2[ii] = eigen_vectors.data()[ii+3*1];
            dprin_stress_dx_vect3[ii] = eigen_vectors.data()[ii+3*2];
          }
          //cout.precision(15);
          //cout<<"1st order derivative of Cauchy stress at GP "<<gg<<" : \t"<<dStress_dx<<endl;
          //cout<<"1st order derivative of principa stress at GP "<<gg<<" :\t"<<dprin_vals_vect<<endl;
          if (gg==0) {
            // Output 1st-order derivative of stress at Gauss point
            StressGP_r.resize(3,3); StressGP_r.clear(); StressGP_r = dStress_dx;
            // Output 1st-order derivative of strain at Gauss point
            StrainGP_r.resize(3,3); StrainGP_r.clear(); StrainGP_r = Strain_du;
          }
          //Tag principle stress vectors 1, 2, 3
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect1,&mit->second,1,&dprin_stress_dx_vect1[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect2,&mit->second,1,&dprin_stress_dx_vect2[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect3,&mit->second,1,&dprin_stress_dx_vect3[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vals,&mit->second,1,&dprin_vals_vect[0]); CHKERR_PETSC(rval);
          
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }
    
  };
  
  // ===========================================================================
  //
  // Calculationg the second order derivative of stress
  //
  // ===========================================================================
  
  struct FE2_PostProcStressForReliability_Second: public PostProcDisplacemenysAndStarinOnRefMesh {
    
    FieldInterface& mField;
    
    ublas::matrix<double> Dmat;
    ublas::matrix<double> Dmat_r;
    ublas::matrix<double> Dmat_rs;
    ublas::matrix<double> StressGP_r;
    ublas::matrix<double> StressGP_rs;
    ublas::matrix<double> StrainGP_r;
    ublas::matrix<double> StrainGP_rs;
    string zeroth_field;
    string first_field;
    string second_field;
    
    Tag th_stress,th_prin_stress_vect1,th_prin_stress_vect2,th_prin_stress_vect3,th_prin_stress_vals;
    
    FE2_PostProcStressForReliability_Second(FieldInterface& _mField,
                                            string _zeroth_field,
                                            string _first_field,
                                            string _second_field,
                                            ublas::matrix<double> _Dmat,
                                            ublas::matrix<double> _Dmat_r,
                                            ublas::matrix<double> _Dmat_rs):
    PostProcDisplacemenysAndStarinOnRefMesh(_mField.get_moab(),_zeroth_field),
    mField(_mField),
    zeroth_field(_zeroth_field),
    first_field(_first_field),
    second_field(_second_field),
    Dmat(_Dmat), Dmat_r(_Dmat_r), Dmat_rs(_Dmat_rs){
      double def_VAL2[3] = { 0.0, 0.0, 0.0 };
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT1",3,MB_TYPE_DOUBLE,th_prin_stress_vect1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT2",3,MB_TYPE_DOUBLE,th_prin_stress_vect2,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VECT3",3,MB_TYPE_DOUBLE,th_prin_stress_vect3,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("PRIN_STRESS_VALUES",3,MB_TYPE_DOUBLE,th_prin_stress_vals,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0 };
      rval = moab_post_proc.tag_get_handle("STRESS_VAL",9,MB_TYPE_DOUBLE,th_stress,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      
    }
    
    
    vector< ublas::matrix<FieldData> > invH;
    vector< FieldData > detH;
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      
      try {
        
        //Loop over elements
        ierr = do_operator(); CHKERRQ(ierr);
        
        //Higher order approximation of geometry
        ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
        
        //Strains to Nodes in PostProc Mesh: create vector containing matrices
        vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
        vector< ublas::matrix< FieldData > > Grad_dU_at_GaussPt;
        vector< ublas::matrix< FieldData > > Grad_ddU_at_GaussPt;
        
        //du/dx
        ierr = GetGaussDiffDataVector(zeroth_field,GradU_at_GaussPt); CHKERRQ(ierr);
        ierr = GetGaussDiffDataVector(first_field,Grad_dU_at_GaussPt); CHKERRQ(ierr);
        ierr = GetGaussDiffDataVector(second_field,Grad_ddU_at_GaussPt); CHKERRQ(ierr);
        
        vector< ublas::matrix< FieldData > >::iterator viit     = GradU_at_GaussPt.begin();
        vector< ublas::matrix< FieldData > >::iterator viit_du  = Grad_dU_at_GaussPt.begin();
        vector< ublas::matrix< FieldData > >::iterator viit_ddu = Grad_ddU_at_GaussPt.begin();
        
        map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
        int gg = 0;
        
        // calculate zeroth-order displacement iduced strain/stress
        for(;viit!=GradU_at_GaussPt.end();viit++,viit_du++,viit_ddu++,mit++,gg++) {
          
          
          ublas::matrix< FieldData > GradU = *viit;
          ublas::matrix< FieldData > Grad_dU = *viit_du;
          ublas::matrix< FieldData > Grad_ddU = *viit_ddu;
          
          if(!invH.empty()) {
            //GradU =
            //[ dU/dChi1 dU/dChi2 dU/dChi3 ]
            //[ dV/dChi1 dV/dChi2 dU/dChi3 ]
            //[ dW/dChi1 dW/dChi2 dW/dChi3 ]
            //H =
            //[ dX1/dChi1 dX1/dChi2 dX1/dChi3 ]
            //[ dX2/dChi1 dX2/dChi2 dX2/dChi3 ]
            //[ dX3/dChi1 dX3/dChi2 dX3/dChi3 ]
            //invH =
            //[ dChi1/dX1 dChi1/dX2 dChi1/dX3 ]
            //[ dChi2/dX1 dChi2/dX2 dChi2/dX3 ]
            //[ dChi3/dX1 dChi3/dX2 dChi3/dX3 ]
            //GradU =
            //[ dU/dX1 dU/dX2 dU/dX3 ]
            //[ dV/dX1 dV/dX2 dV/dX3 ] = GradU * invH
            //[ dW/dX1 dW/dX2 dW/dX3 ]
            GradU    = prod( GradU,    invH[gg] );
            Grad_dU  = prod( Grad_dU,  invH[gg] );
            Grad_ddU = prod( Grad_ddU, invH[gg] );
          }
          ublas::matrix< FieldData > Strain;     Strain.resize(3,3);     Strain     = 0.5*( GradU    + trans(GradU)   );
          ublas::matrix< FieldData > Strain_du;  Strain_du.resize(3,3);  Strain_du  = 0.5*( Grad_dU  + trans(Grad_dU) );
          ublas::matrix< FieldData > Strain_ddu; Strain_ddu.resize(3,3); Strain_ddu = 0.5*( Grad_ddU + trans(Grad_ddU) );
          
          //rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
          //rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
          
          //caluate stress and save it into tag
          ublas::vector<FieldData> Strain_VectorNotation(6);
          Strain_VectorNotation[0] = Strain(0,0);
          Strain_VectorNotation[1] = Strain(1,1);
          Strain_VectorNotation[2] = Strain(2,2);
          Strain_VectorNotation[3] = 2*Strain(0,1);
          Strain_VectorNotation[4] = 2*Strain(1,2);
          Strain_VectorNotation[5] = 2*Strain(2,0);
          
          ublas::vector<FieldData> Strain_du_VectorNotation(6);
          Strain_du_VectorNotation[0] = Strain_du(0,0);
          Strain_du_VectorNotation[1] = Strain_du(1,1);
          Strain_du_VectorNotation[2] = Strain_du(2,2);
          Strain_du_VectorNotation[3] = 2*Strain_du(0,1);
          Strain_du_VectorNotation[4] = 2*Strain_du(1,2);
          Strain_du_VectorNotation[5] = 2*Strain_du(2,0);
          
          ublas::vector<FieldData> Strain_ddu_VectorNotation(6);
          Strain_ddu_VectorNotation[0] = Strain_ddu(0,0);
          Strain_ddu_VectorNotation[1] = Strain_ddu(1,1);
          Strain_ddu_VectorNotation[2] = Strain_ddu(2,2);
          Strain_ddu_VectorNotation[3] = 2*Strain_ddu(0,1);
          Strain_ddu_VectorNotation[4] = 2*Strain_ddu(1,2);
          Strain_ddu_VectorNotation[5] = 2*Strain_ddu(2,0);
          
          // [ddSigma/ddX] = [ddD/ddX][Epsilon] + 2* [dD/dX][dEpsilon/dX] + [D] [ddEpsilon/ddX]
          ublas::vector< FieldData > ddStress_ddx_VectorNotation(6); ddStress_ddx_VectorNotation.clear();
          ddStress_ddx_VectorNotation = prod( Dmat_rs, Strain_VectorNotation )
                                      + prod( Dmat_r,  Strain_du_VectorNotation )
                                      + prod( Dmat,    Strain_ddu_VectorNotation );
          
          ublas::matrix< FieldData > ddStress_ddx = ublas::zero_matrix<FieldData>(3,3);
          ddStress_ddx(0,0) = ddStress_ddx_VectorNotation[0];
          ddStress_ddx(1,1) = ddStress_ddx_VectorNotation[1];
          ddStress_ddx(2,2) = ddStress_ddx_VectorNotation[2];
          ddStress_ddx(0,1) = ddStress_ddx(1,0) = ddStress_ddx_VectorNotation[3];
          ddStress_ddx(1,2) = ddStress_ddx(2,1) = ddStress_ddx_VectorNotation[4];
          ddStress_ddx(2,0) = ddStress_ddx(0,2) = ddStress_ddx_VectorNotation[5];
          //cout<<"The d stress: "<<dStress_dx<<"\n"<<prod( Dmat, Strain_du_VectorNotation ) << "\n"<< prod( Dmat_r, Strain_VectorNotation )<<endl;
          
          rval = moab_post_proc.tag_set_data(th_stress,&mit->second,1,&(ddStress_ddx.data()[0])); CHKERR_PETSC(rval);
          
          //Principal stress vectors
          //Calculated from finding the eigenvalues and eigenvectors
          ublas::matrix< FieldData > eigen_vectors = ddStress_ddx;
          ublas::vector<double> eigen_values(3);
          
          //LAPACK - eigenvalues and vectors. Applied twice for initial creates memory space
          int n = 3, lda = 3, info, lwork = -1;
          double wkopt;
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),&wkopt,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          lwork = (int)wkopt;
          double work[lwork];
          info = lapack_dsyev('V','U',n,&(eigen_vectors.data()[0]),lda,&(eigen_values.data()[0]),work,lwork);
          if(info != 0) SETERRQ1(PETSC_COMM_SELF,1,"is something wrong with lapack_dsyev info = %d",info);
          
          //Combine eigenvalues and vectors to create principal stress vector
          ublas::vector<double> ddprin_stress_ddx_vect1(3);
          ublas::vector<double> ddprin_stress_ddx_vect2(3);
          ublas::vector<double> ddprin_stress_ddx_vect3(3);
          ublas::vector<double> ddprin_vals_vect(3);
          
          //eigen_vectors = trans(eigen_vectors);
          for (int ii=0; ii < 3; ii++) {
            ddprin_vals_vect[0] = eigen_values[0];
            ddprin_vals_vect[1] = eigen_values[1];
            ddprin_vals_vect[2] = eigen_values[2];
            ddprin_stress_ddx_vect1[ii] = eigen_vectors.data()[ii+3*0];
            ddprin_stress_ddx_vect2[ii] = eigen_vectors.data()[ii+3*1];
            ddprin_stress_ddx_vect3[ii] = eigen_vectors.data()[ii+3*2];
          }
          //cout.precision(15);
          //cout<<"1st order derivative of Cauchy stress at GP "<<gg<<" : \t"<<dStress_dx<<endl;
          //cout<<"1st order derivative of principa stress at GP "<<gg<<" :\t"<<dprin_vals_vect<<endl;
          if (gg==0) {
            // Output 2nd-order derivative of stress at Gauss point
            StressGP_rs.resize(3,3); StressGP_rs.clear(); StressGP_rs = ddStress_ddx;
            // Output 2nd-order derivative of strain at Gauss point
            StrainGP_rs.resize(3,3); StrainGP_rs.clear(); StrainGP_rs = Strain_ddu;
          }
          //Tag principle stress vectors 1, 2, 3
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect1,&mit->second,1,&ddprin_stress_ddx_vect1[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect2,&mit->second,1,&ddprin_stress_ddx_vect2[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vect3,&mit->second,1,&ddprin_stress_ddx_vect3[0]); CHKERR_PETSC(rval);
          rval = moab_post_proc.tag_set_data(th_prin_stress_vals,&mit->second,1,&ddprin_vals_vect[0]); CHKERR_PETSC(rval);
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << "throw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }
      
      PetscFunctionReturn(0);
    }
    
  };
}

#endif //__FE2_POSTPROCSTRESSFORRELIABILITY_HPP__
