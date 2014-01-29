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

#ifndef __ELASTICFEMETHOD_RVE_TRY_HPP__
#define __ELASTICFEMETHOD_RVE_TRY_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethod.hpp"

namespace MoFEM {

struct ElasticFEMethod_RVE_Try: public ElasticFEMethod {
    
    Vec D_star;
    ublas::matrix<double> strain_app;
    ElasticFEMethod_RVE_Try(
                    FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
                    double _lambda,double _mu, Vec& _D_star, ublas::matrix<double> _strain_app):
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), D_star(_D_star), strain_app(_strain_app){};
    
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        Range nodes;
        double coords[3], disp_applied[3];
        EntityHandle node = 0;
        int index_D_star[3];
        rval = moab.get_entities_by_type(0,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
        for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
            //cout<<"\n New postProcess \n";
            int jj=0;
            for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
                EntityHandle ent = dof->get_ent();
                if(node!=ent) {
                    rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
                    //cout<<"\n\n coord = " << coords[0]<<" "<< coords[1]<<" " << coords[2] << "\n\n";
                    for(int ii=0; ii<3; ii++) disp_applied[ii]=strain_app(ii,0)*coords[0] + strain_app(ii,1)*coords[1] + strain_app(ii,2)*coords[2];
                    //cout<<"\n\n strain_app "<<strain_app<<"\n\n";
                    //cout<<"\n\n disp_applied = " << disp_applied[0]<<" "<< disp_applied[1]<<" " << disp_applied[2] << "\n\n";
                    node = ent;
                }
                index_D_star[jj]=dof->get_petsc_gloabl_dof_idx(); jj++;
            }
            ierr = VecSetValues(D_star,3,&(index_D_star[0]),&(disp_applied[0]),INSERT_VALUES); CHKERRQ(ierr);
        }

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
        ierr = VecZeroEntries(Data); CHKERRQ(ierr);
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_FieldData(this,Data); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_MatrixDiagonal(this,Aij); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,F); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        
        ierr = VecAssemblyBegin(D_star); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(D_star); CHKERRQ(ierr);
        
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode Fint(Vec F_int) {
        PetscFunctionBegin;
        ierr = Fint(); CHKERRQ(ierr);
        for(int rr = 0;rr<row_mat;rr++) {
            if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
            if(RowGlob[rr].size()==0) continue;
            f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
            ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
    }


    virtual PetscErrorCode Fint() {
        PetscFunctionBegin;
        
         try {
            
            double _lambda,_mu;
            ierr = GetMatParameters(&_lambda,&_mu); CHKERRQ(ierr);
            ierr = calulateD(_lambda,_mu); CHKERRQ(ierr);
            
            //Gradient at Gauss points;
            vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
            ierr = GetGaussDiffDataVector("DISPLACEMENT_FROM_APP_STRAIN",GradU_at_GaussPt); CHKERRQ(ierr);
            
            //cout<<"GradU_at_GaussPt "<< GradU_at_GaussPt.size() << endl;
            
            unsigned int g_dim = g_NTET.size()/4;
            assert(GradU_at_GaussPt.size() == g_dim);
            NOT_USED(g_dim);
            vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
            int gg = 0;
            for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
                try {
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
                } catch (const std::exception& ex) {
                    ostringstream ss;
                    ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
                    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
                } 
            }
            
        } catch (const std::exception& ex) {
            ostringstream ss;
            ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
            SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
        }
        
        PetscFunctionReturn(0);
    }
   
};

    
}

#endif //ELASTICFEMETHOD_RVE_TRY_HPP
