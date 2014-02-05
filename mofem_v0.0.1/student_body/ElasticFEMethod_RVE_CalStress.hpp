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
#include <boost/numeric/ublas/storage.hpp>
#include "ElasticFEMethod.hpp"
#include "FieldInterface.hpp"
#include "FieldCore.hpp"


namespace MoFEM {


struct ElasticFEMethod_RVE_CalStress: public ElasticFEMethod{
    
    Vec F_stress, coord_stress;
    Interface& moab;
    string field_name;
    double RVE_volume;
    ublas::matrix<double> Sigma_homo;

    
    ElasticFEMethod_RVE_CalStress(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
          double _lambda,double _mu, Interface& _moab, Vec&_F_stress, Vec&_coord_stress, double _RVE_volume, ublas::matrix<double> _Sigma_homo, string _field_name = "DISPLACEMENT" ):
    ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu), moab(_moab), F_stress(_F_stress), coord_stress(_coord_stress), RVE_volume(_RVE_volume), Sigma_homo(_Sigma_homo), field_name(_field_name){};
    
  
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
        
        RVE_volume=0.0;
        PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = VecAssemblyBegin(F_stress); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F_stress); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(coord_stress); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(coord_stress); CHKERRQ(ierr);
        
        ublas::matrix< double > F_stress_mat, coord_stress_mat, Sigma_homo;
        int array_indices[3];
        double array_output[3];
        
        for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|UnknownCubitName,it)) {
            std::size_t BoundNodesFound=it->get_Cubit_name().find("BoundNodes");
            if (BoundNodesFound==std::string::npos) continue;
            Range tris, nodes;
            rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
            rval = moab.get_connectivity(tris,nodes,true); CHKERR_PETSC(rval);
            
            F_stress_mat.resize(3,nodes.size()); coord_stress_mat.resize(nodes.size(),3);
            F_stress_mat.clear();  coord_stress_mat.clear();
            
            int ii=0;
            for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
                int jj=0;
                for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problem_ptr,*nit,dof)) {
                    array_indices[jj]=dof->get_petsc_gloabl_dof_idx(); jj++;
                }
                ierr = VecGetValues(F_stress,3,array_indices,array_output); CHKERRQ(ierr);
                F_stress_mat(0,ii)=array_output[0]; F_stress_mat(1,ii)=array_output[1],  F_stress_mat(2,ii)=array_output[2];
                ierr = VecGetValues(coord_stress,3,array_indices,array_output); CHKERRQ(ierr);
                coord_stress_mat(ii,0)=array_output[0]; coord_stress_mat(ii,1)=array_output[1], coord_stress_mat(ii,2)=array_output[2];
                ii++;
            }
        }
        
//        cout<<"\n\n\n\n\n\n\n F_stress_mat \n\n\n\n\n\n\\n\n\n\n\n";
//        cout<<F_stress_mat;
//        cout<<"\n\n\n\n\n\n\n coord_stress_mat \n\n\n\n\n\n\\n\n\n\n\n";
//        cout<<coord_stress_mat;
//        cout<<"\n\n\n\n\n\n\n RVE_volume \n\n\n\n\n\n\\n\n\n\n\n";
        cout<<RVE_volume;
        Sigma_homo.resize(3,3);  Sigma_homo.clear();
        
        Sigma_homo=prod(F_stress_mat, coord_stress_mat);
        Sigma_homo=(1.0/ RVE_volume)*(Sigma_homo);
        cout<<"\n\n\nSigma_bar "<<Sigma_homo<<endl<<endl<<endl;
      
        
        

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
        
            unsigned int g_dim = g_NTET.size()/4;
            assert(GradU_at_GaussPt.size() == g_dim);
            NOT_USED(g_dim);
            vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
            int gg = 0;
            for(;viit!=GradU_at_GaussPt.end();viit++,gg++) {
                    ublas::matrix< FieldData > GradU = *viit;
                    ublas::matrix< FieldData > Strain = 0.5*( GradU + trans(GradU) );
                
//                    cout<<"Strain "<< Strain << endl;

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
        RVE_volume+=V;
        
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

#endif //ElasticFEMethod_RVE_CalStress
