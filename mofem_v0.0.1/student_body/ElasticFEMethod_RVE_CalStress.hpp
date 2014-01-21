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

#ifndef __ElasticFEMethod_RVE_CalStress_HPP__
#define __ElasticFEMethod_RVE_CalStress_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethod.hpp"
#include "FieldInterface.hpp"
#include "FieldCore.hpp"


namespace MoFEM {


struct ElasticFEMethod_RVE_CalStress: public ElasticFEMethod{
    
    Vec F_stress;
    Interface& moab;
    string field_name;
    
    
    ElasticFEMethod_RVE_CalStress(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
          double _lambda,double _mu, Interface& _moab, string _field_name = "DISPLACEMENT" ):
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), moab(_moab), field_name(_field_name){};
    
  
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        g_NTET.resize(4*45);
        ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
        G_W_TET = G_TET_W45;
        g_NTRI.resize(3*13);
        ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13);
        G_W_TRI = G_TRI_W13;
        
        // See FEAP - - A Finite Element Analysis Program
        D_lambda.resize(6,6);
        D_lambda.clear();
        for(int rr = 0;rr<3;rr++) {
            for(int cc = 0;cc<3;cc++) {
                D_lambda(rr,cc) = 1;
            }
        }
        D_mu.resize(6,6);
        D_mu.clear();
        for(int rr = 0;rr<6;rr++) {
            D_mu(rr,rr) = rr<3 ? 2 : 1;
        }
        
        ierr = VecDuplicate(F,&F_stress); CHKERRQ(ierr);
        ierr = VecZeroEntries(F_stress); CHKERRQ(ierr);
        ierr = VecGhostUpdateBegin(F_stress,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        ierr = VecGhostUpdateEnd(F_stress,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = VecAssemblyBegin(F_stress); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F_stress); CHKERRQ(ierr);
        
//        cout<<"\nFinal F_stress vector\n";
//        ierr = VecView(F_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

        
        ublas::matrix<double> F_stress_vec, coord_stress, Sigma_bar;
        PetscScalar *F_stress_array;
        
        Range nodes;
        double coords[3], disp_applied[3];
        EntityHandle node = 0;
        rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
        
        F_stress_vec.resize(3,nodes.size());  F_stress_vec.clear();
        coord_stress.resize(nodes.size(),3);  coord_stress.clear();
        Sigma_bar.resize(3,3);  Sigma_bar.clear();
        
        ierr = VecGetArray(F_stress,&F_stress_array); CHKERRQ(ierr);
        
        for(int ii=0; ii<nodes.size(); ii++){
            F_stress_vec(0,ii)=F_stress_array[3*ii+0];
            F_stress_vec(1,ii)=F_stress_array[3*ii+1];
            F_stress_vec(2,ii)=F_stress_array[3*ii+2];
        }
        
        for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
                EntityHandle ent = dof->get_ent();
                if(node!=ent) {
                    rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
                    //cout<<"\n\n coord = " << coords[0]<<" "<< coords[1]<<" " << coords[2] << "\n\n";
                    int nn=dof->get_petsc_gloabl_dof_idx(); nn=nn/3;
                    //cout<< " \n n ="<< nn << endl;
                    coord_stress(nn,0)=coords[0];  coord_stress(nn,1)=coords[1];  coord_stress(nn,2)=coords[2];
                    node = ent;
                }
            }
        }
//        cout<<"Sigma_bar.size1() "<<Sigma_bar.size1()<<endl;
//        cout<<"Sigma_bar.size2() "<<Sigma_bar.size2()<<endl;
//        cout<<"F_stress_vec.size2() "<<F_stress_vec.size2()<<endl;
//        cout<<"F_stress_vec.size2() "<<*F_stress_vec.data().begin()<<endl;
        cout<<"Sigma_bar "<<Sigma_bar<<endl;
        
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                    Sigma_bar.size1(),Sigma_bar.size2(),F_stress_vec.size2(),
                    1.,&*F_stress_vec.data().begin(),F_stress_vec.size2(),
                    &*coord_stress.data().begin(),coord_stress.size2(),
                    0.,&*Sigma_bar.data().begin(),Sigma_bar.size2());
        
        cout<<"Sigma_bar "<<Sigma_bar<<endl;
        
        //cout<<"coord_stress "<<coord_stress<<endl;
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        PetscFunctionReturn(0);
    }

    
    vector<ublas::vector<FieldData> > f_int;
    virtual PetscErrorCode Fint() {
        PetscFunctionBegin;
            double _lambda,_mu;
            ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
            ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);
            
            //Gradient at Gauss points;
            vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
            ierr = GetGaussDiffDataVector(field_name,GradU_at_GaussPt); CHKERRQ(ierr);
            
            //cout<<"GradU_at_GaussPt "<< GradU_at_GaussPt.size() << endl;
            //cout<<"GradU_at_GaussPt "<< GradU_at_GaussPt[0] << endl;
        
            unsigned int g_dim = g_NTET.size()/4;
            assert(GradU_at_GaussPt.size() == g_dim);
            NOT_USED(g_dim);
            vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
            int gg = 0;
            for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
                    ublas::matrix< FieldData > GradU = *viit;
                    ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
                    ublas::vector< FieldData > VoightStrain(6);
                    VoightStrain[0] = Strain(0,0);
                    VoightStrain[1] = Strain(1,1);
                    VoightStrain[2] = Strain(2,2);
                    VoightStrain[3] = 2*Strain(0,1);
                    VoightStrain[4] = 2*Strain(1,2);
                    VoightStrain[5] = 2*Strain(2,0);
                    double w = V*G_W_TET[gg];
                    ublas::vector<FieldData> VoightStress = prod(w*D,VoightStrain);
                    //BT * VoigtStress
                    for(int rr = 0;rr<row_mat;rr++) {
                        f_int.resize(row_mat);
                        ublas::matrix<FieldData> &B = (rowBMatrices[rr])[gg];
                        if(gg == 0) {
                            f_int[rr] = prod( trans(B), VoightStress );
                        } else {
                            f_int[rr] += prod( trans(B), VoightStress );
                        }
                    }
            }
        
        PetscFunctionReturn(0);
    }

    
    virtual PetscErrorCode Fint(Vec F_stress) {
        PetscFunctionBegin;
            ierr = Fint(); CHKERRQ(ierr);
            for(int rr = 0;rr<row_mat;rr++) {
                if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
                if(RowGlob[rr].size()==0) continue;
                ierr = VecSetValues(F_stress,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
            }
        PetscFunctionReturn(0);
    }
  

    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        
        //Assembly F_stress
        //cout<<"\nBefore\n";
        //ierr = VecView(F_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        ierr = Fint(F_stress); CHKERRQ(ierr);
        //cout<<"\nAfter\n";
        //ierr = VecView(F_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
 };

    
}

#endif //ElasticFEMethod_RVE_CalStress
