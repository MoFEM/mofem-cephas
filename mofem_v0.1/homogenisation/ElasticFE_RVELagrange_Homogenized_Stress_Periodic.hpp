/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ElasticFE_RVELagrange_Homogenized_Stress_Periodic_HPP__
#define __ElasticFE_RVELagrange_Homogenized_Stress_Periodic__


#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Periodic.hpp"

namespace MoFEM {

struct ElasticFE_RVELagrange_Homogenized_Stress_Periodic: public ElasticFE_RVELagrange_Periodic {

    Vec DVec;
    Vec Stress_Homo;
    double *RVE_volume;

    ElasticFE_RVELagrange_Homogenized_Stress_Periodic(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,double *_RVE_volume):
    ElasticFE_RVELagrange_Periodic(_mField, _dirihlet_ptr,_Aij, _D, _F), DVec(_D),RVE_volume(_RVE_volume){};
    
    
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        
        //create a vector for 6 components of homogenized stress
        ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
        ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);

        
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = VecAssemblyBegin(Stress_Homo); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(Stress_Homo); CHKERRQ(ierr);
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        
        ierr = VecScale(Stress_Homo, 1.0/(*RVE_volume)); CHKERRQ(ierr);
        if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
        ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }
    
    
    ublas::vector<ublas::vector<FieldData> > Lamda;
    
    virtual PetscErrorCode Calculate_Homo_Stress() {
        PetscFunctionBegin;
        X_mat.resize(3,6);    X_mat.clear();
        nodes_coord.resize(3,3);
        gauss_coord.resize(3,g_TRI_dim);
        D_mat.resize(row_mat);
        Lamda.resize(row_mat);
        
        for(int ff=0; ff<2; ff++){
            ublas::vector<FieldData>  Stress_Homo_elem; //define and clear here
            Stress_Homo_elem.resize(6);   Stress_Homo_elem.clear();   //homogenised stress for one element (triangle)

            D_mat.resize(row_mat);
            //used to calculate the coordinates of a Gauss points
            nodes_coord(0,0)=coords_face[ff][0]; nodes_coord(0,1)=coords_face[ff][3]; nodes_coord(0,2)=coords_face[ff][6];
            nodes_coord(1,0)=coords_face[ff][1]; nodes_coord(1,1)=coords_face[ff][4]; nodes_coord(1,2)=coords_face[ff][7];
            nodes_coord(2,0)=coords_face[ff][2]; nodes_coord(2,1)=coords_face[ff][5]; nodes_coord(2,2)=coords_face[ff][8];
            
            //coordinates for all gauss points
            gauss_coord=prod(nodes_coord, g_NTRI_mat);
            
//            cout<<"g_NTRI_mat "<<g_NTRI_mat<<endl<<endl;
//            cout<<"nodes_coord "<<nodes_coord<<endl<<endl;
//            cout<<"gauss_coord "<<gauss_coord<<endl<<endl;
//            std::string wait;
//            std::cin >> wait;
//              cout<<"row_mat "<<row_mat<<endl;
            
            for(int rr=0; rr<row_mat; rr++){
                for(int gg = 0;gg<g_TRI_dim;gg++) {
                    double w = area*G_W_TRI[gg];
                    X_mat(0,0)=2.0*gauss_coord(0,gg);  X_mat(0,3)=gauss_coord(1,gg);  X_mat(0,4)=gauss_coord(2,gg);
                    X_mat(1,1)=2.0*gauss_coord(1,gg);  X_mat(1,3)=gauss_coord(0,gg);  X_mat(1,5)=gauss_coord(2,gg);
                    X_mat(2,2)=2.0*gauss_coord(2,gg);  X_mat(2,4)=gauss_coord(0,gg);  X_mat(2,5)=gauss_coord(1,gg);
                    X_mat=0.5*X_mat;
                    
                    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
                    ublas::matrix<FieldData> &col_Mat = X_mat;
                    
                    ublas::matrix<FieldData>  D_mat1;    //Dmat1=NT*X_mat
                    D_mat1.resize(row_Mat.size2(),col_Mat.size2());
                    //Integrate D_mat
                    if(gg == 0) {
                        D_mat[rr].resize(H_mat[ff][rr].size1(),D_mat1.size2());
                        
                        //calculate (D_mat1= w * NT * X_mat)
//                        cout<<"\n row_Mat "<<row_Mat<<endl;
//                        cout<<"\n col_Mat "<<rr<<col_Mat;
//                        cout<<"\n w "<<w<<endl;
                        
                        D_mat1=prod(w*trans(row_Mat), col_Mat);
                        //calculate (D_mat = H_mat * D_mat1)
//                        cout<<"\n rr "<<rr<<endl;
//                        cout<<"\n D_mat[rr] "<<D_mat[rr]<<endl;
                        D_mat[rr]=prod(H_mat[ff][rr], D_mat1);
                    }else{
                        //calculate (D_mat1= w * NT * X_mat)
                        D_mat1=prod(w*trans(row_Mat), col_Mat);
                        
                        //calculate (D_mat = H_mat * D_mat1)
                        D_mat[rr]+=prod(H_mat[ff][rr], D_mat1);
                    }
                }
                
                //To get Lambda for stress calculation
                Lamda[rr].resize(RowGlob[ff][rr].size());  Lamda[rr].clear();
                
                //Canonical numbering of two triangles (actual canonical number is [3, 4], here we use only 3 as we have lagrange multipliers only for one
                int ff_arr[]={3, 3};
                const EntityHandle*  conn_Prism;  int num_nodes_Prism;
                prism_periodic=fe_ptr->get_ent();
                rval = moab.get_connectivity(prism_periodic,conn_Prism,num_nodes_Prism,true); CHKERR_PETSC(rval);
//                cout<<"num_nodes_Prism  "<<num_nodes_Prism<<endl;
                EntityHandle Tri_Prism;
                rval = moab.side_element(prism_periodic,2,ff_arr[ff],Tri_Prism); CHKERR_PETSC(rval);
                switch(rr) {
                    case 0:  //for nodes
//                        cout<<"For nodes"<<endl;
                        for(int nn = 0;nn<num_nodes_Prism/2; nn++) {
                            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"Lagrange_mul_disp",conn_Prism[nn],iit)) {
                                Lamda[rr][3*nn+iit->get_dof_rank()]=iit->get_FieldData();
                            }
                        }
//                        for(int ii=0; ii<Lamda[rr].size(); ii++) cout<<Lamda[rr][ii]<<" ";
//                        cout<<endl;
                        break;
                    case 1:  case 2:  case 3:  //For edges
//                    cout<<"For Edges"<<endl;
                        for(int ee=0; ee<3; ee++) {
                            EntityHandle edge;
                            rval = moab.side_element(Tri_Prism,1,rr-1,edge); CHKERR_PETSC(rval);
                            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"Lagrange_mul_disp",edge,iit)) {
                                Lamda[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
                            }
                        }
//                        for(int ii=0; ii<Lamda[rr].size(); ii++) cout<<Lamda[rr][ii]<<" ";
//                        cout<<endl;
                        break;
                    case 4: //for face
//                    cout<<"For Face"<<endl;
                        for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"Lagrange_mul_disp",Tri_Prism,iit)) {
                            Lamda[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
                        }
//                        for(int ii=0; ii<Lamda[rr].size(); ii++) cout<<Lamda[rr][ii]<<" ";
//                        cout<<endl;
                        break;
                }
                
//                cout<< " Stress_Homo_elem before =  "<< Stress_Homo_elem << endl;
                Stress_Homo_elem+=prod(trans(D_mat[rr]), -1*Lamda[rr]);   //Lamda is reaction force (so multiply for -1 to get the force)
            }
            int Indices[6]={0, 1, 2, 3, 4, 5};
            ierr = VecSetValues(Stress_Homo,6,Indices,&(Stress_Homo_elem.data())[0],ADD_VALUES); CHKERRQ(ierr);
        }
        PetscFunctionReturn(0);
    }
    

     PetscErrorCode operator()() {
        PetscFunctionBegin;
//        cout<<"Hi from class ElasticFE_RVELagrange_Homogenized_Stress_Periodic"<<endl;
         
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
        ierr = Get_H_mat();   //It will be used from the Class ElasticFE_RVELagrange_Disp
        ierr = Calculate_Homo_Stress(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    
};

    
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
