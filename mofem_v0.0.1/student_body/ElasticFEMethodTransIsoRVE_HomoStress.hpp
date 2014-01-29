/* Copyright (C) 2013, Michel Cortis <mikecortis at gmail.com>
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

#ifndef __ELASTICFEMETHODTRANSISORVE_HomoStress_HPP__
#define __ELASTICFEMETHODTRANSISORVE_HomoStress_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethodTransIso.hpp"


namespace MoFEM {
    
    struct ElasticFEMethodTransIsoRVE_HomoStress: public TranIsotropicFibreDirRotElasticFEMethod {
        
        Vec F_stress, coord_stress;
        Interface& moab;
        string field_name;
        double *RVE_volume_TransIso;

        ElasticFEMethodTransIsoRVE_HomoStress(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _D,Vec& _F, double _lambda,double _mu,double _E_p,double _E_z, double _nu_p,double _nu_pz, double _G_zp, Interface& _moab, Vec&_F_stress, Vec&_coord_stress, double *_RVE_volume_TransIso, string _field_name = "DISPLACEMENT"):
        
        TranIsotropicFibreDirRotElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu,_E_p, _E_z, _nu_p, _nu_pz, _G_zp), moab(_moab), F_stress(_F_stress), coord_stress(_coord_stress), RVE_volume_TransIso(_RVE_volume_TransIso), field_name(_field_name)  {
            
            double def_VAL2[3] = {0,0,0};
            // create a new tag th_fibre_dir
            rval = moab.tag_get_handle( "POT_FLOW_FIBRE_DIR",3,MB_TYPE_DOUBLE,th_fibre_dir,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL2); CHKERR_THROW(rval);
        };
        
        
        vector< ublas::matrix<FieldData> > D_At_GaussPoint;
         PetscErrorCode Fint() {
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
                ublas::vector<FieldData> VoightStress = prod(w*D_At_GaussPoint[gg],VoightStrain);
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
        
        PetscErrorCode Fint(Vec F_stress) {
            PetscFunctionBegin;
            ierr = Fint(); CHKERRQ(ierr);
            for(int rr = 0;rr<row_mat;rr++) {
                if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
                if(RowGlob[rr].size()==0) continue;
                f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
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
            
            ///Get Stiffness Matrix
            ublas::symmetric_matrix<FieldData,ublas::upper> StiffnessMatrix;
            StiffnessMatrix.resize(6);
            StiffnessMatrix.clear();
            TransverseIsotropicStiffnessMatrix TranIsoMat(nu_p, nu_pz, E_p, E_z, G_zp);
            StiffnessMatrix=TranIsoMat.StiffnessMatrix;
            //        IsotropicStiffnessMatrix IsoMat(lambda, mu);
            //        StiffnessMatrix=IsoMat.StiffnessMatrix;
            
            ///Rotating the Stiffness matrix according a set of axes of rotations and their respective angle
            
            D_At_GaussPoint.resize(coords_at_Gauss_nodes.size());
            
            for(unsigned int gg=0;gg<coords_at_Gauss_nodes.size();gg++){
                
                double fVec[3];
                ierr = ComputeFibreDirection(&fVec[0]); CHKERRQ(ierr);
                
                int noOfRotations = 1; //Number of Rotations
                
                double zVec[3]={0.0,0.0,1.0};
                double AxVector[3]={fVec[1]*zVec[2]-fVec[2]*zVec[1] , fVec[2]*zVec[0]-fVec[0]*zVec[2] , fVec[0]*zVec[1]-fVec[1]*zVec[0]};
                double AxAngle[1]= {asin((sqrt(pow(AxVector[0],2)+pow(AxVector[1],2)+pow(AxVector[2],2)))/(sqrt(pow(fVec[0],2)+pow(fVec[1],2)+pow(fVec[2],2)))*(sqrt(pow(zVec[0],2)+pow(zVec[1],2)+pow(zVec[2],2))))};
                
                double negAxAngle[noOfRotations];
                for (int aa=0; aa<noOfRotations; aa++) negAxAngle[aa]=-AxAngle[aa];
                
                ublas::matrix<double> DummyMatrix,DummyMatrix2;
                DummyMatrix = ublas::zero_matrix<FieldData>(6,6);
                DummyMatrix = StiffnessMatrix;
                
                ///Rotating Stiffness over a number of axis/angle rotations
                for (int aa=0; aa<noOfRotations; aa++) {
                    
                    StressTransformation StressRotMat(&AxVector[3*aa], AxAngle[aa]);
                    StrainTransformation invStrainRotMat(&AxVector[3*aa], negAxAngle[aa]);
                    
                    ublas::matrix<double> TrpMatrixStress;
                    TrpMatrixStress = ublas::zero_matrix<FieldData>(6,6);
                    TrpMatrixStress=StressRotMat.StressRotMat;
                    
                    ublas::matrix<double> TrpMatrixInvStrain;
                    TrpMatrixInvStrain = ublas::zero_matrix<FieldData>(6,6); 
                    TrpMatrixInvStrain=invStrainRotMat.StrainRotMat;
                    
                    DummyMatrix2 = ublas::zero_matrix<FieldData>(6,6); 
                    ublas::matrix< FieldData > dummyA = prod( DummyMatrix , TrpMatrixInvStrain );
                    DummyMatrix2 = prod(TrpMatrixStress,dummyA);
                    DummyMatrix = ublas::zero_matrix<FieldData>(6,6); 
                    DummyMatrix = DummyMatrix2;
                }
                
                D_At_GaussPoint[gg].resize(6,6);
                D_At_GaussPoint[gg].clear();
                D_At_GaussPoint[gg] = DummyMatrix;
//                cout<<"\n\n\n D_At_GaussPoint "<<D_At_GaussPoint[gg] << endl;
            }

            //Assembly Aij and F
            ierr = Fint(F_stress); CHKERRQ(ierr);
            ierr = store_coord(coord_stress);
            *RVE_volume_TransIso+=V;
            

            ierr = OpStudentEnd(); CHKERRQ(ierr);
            PetscFunctionReturn(0); }
        
    };

}
#endif //__ELASTICFEMETHODTRANSISORVE_HomoStress_HPP__

