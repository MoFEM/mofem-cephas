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

#ifndef __ElasticFEMethod_RVE_HomoStress_HPP__
#define __ElasticFEMethod_RVE_HomoStress_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include "ElasticFEMethod.hpp"
#include "FieldInterface.hpp"
#include "FieldCore.hpp"


namespace MoFEM {


struct ElasticFEMethod_RVE_HomoStress: public ElasticFEMethod{
    
    Vec F_stress, coord_stress;
    Interface& moab;
    string field_name;
    double *RVE_volume_Elastic;
    

    ElasticFEMethod_RVE_HomoStress(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
          double _lambda,double _mu, Interface& _moab, Vec&_F_stress, Vec&_coord_stress, double *_RVE_volume_Elastic, string _field_name = "DISPLACEMENT"):
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), moab(_moab), F_stress(_F_stress), coord_stress(_coord_stress), RVE_volume_Elastic(_RVE_volume_Elastic), field_name(_field_name){};
    
  
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = VecAssemblyBegin(F_stress); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F_stress); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(coord_stress); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(coord_stress); CHKERRQ(ierr);
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
  

    
    virtual PetscErrorCode store_coord(Vec coord_stress) {
        PetscFunctionBegin;
        
        double coords[3];
        int coord_pos[3];
        EntityHandle fe_handle = fe_ptr->get_ent();
        Range tetNodes;
        EntityHandle node = 0;
        rval = moab.get_connectivity(&fe_handle,1,tetNodes); CHKERR_THROW(rval);
        
        
        for(Range::iterator nit = tetNodes.begin(); nit!=tetNodes.end(); nit++) {
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
                EntityHandle ent = dof->get_ent();
                if(node!=ent) {
                    //cout<<"coord "<<endl;
                    rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
                    //cout<<"\n\n coord = " << coords[0]<<" "<< coords[1]<<" " << coords[2] << "\n\n";
                    int nn=dof->get_petsc_gloabl_dof_idx(); nn=nn/3;
                    //cout<< " \n n ="<< nn << endl;
                    coord_pos[0]=3*nn;  coord_pos[1]=3*nn+1;  coord_pos[2]=3*nn+2;
                    ierr = VecSetValues(coord_stress,3, coord_pos ,coords,INSERT_VALUES); CHKERRQ(ierr);
                    node = ent;
                }
            }
        }
        PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        
        //volume of RVE
        *RVE_volume_Elastic+=V;
//        cout<< "V ="<<V<<endl;
//        cout<< "RVE_volume_Elastic ="<<*RVE_volume_Elastic<<endl;
        
        //Assembly F_stress
        ierr = Fint(F_stress); CHKERRQ(ierr);
        ierr = store_coord(coord_stress);
        
        //cout<<"\nAfter\n";
        //ierr = VecView(F_stress,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
 };

    
}

#endif //ElasticFEMethod_RVE_HomoStress
