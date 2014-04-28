/* Copyright (C) 2013, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#ifndef __ElasticFE_RVELagrange_Homogenized_Stress_HPP__
#define __ElasticFE_RVELagrange_Homogenized_Stress__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"

namespace MoFEM {

struct ElasticFE_RVELagrange_Homogenized_Stress_Disp: public ElasticFE_RVELagrange_Disp {

    Vec DVec;

    ElasticFE_RVELagrange_Homogenized_Stress_Disp(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F):
    ElasticFE_RVELagrange_Disp(_mField, _dirihlet_ptr,_Aij, _D, _F), DVec(_D){};
    
    
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscTime(&v1); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
        
        Stress_Homo.resize(6);   Stress_Homo.clear();
        
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        
        cout<< " Stress_Homo =  "<<endl;
        for(int ii=0; ii<6; ii++) cout<<Stress_Homo(ii)<<endl; 
        
        
        PetscFunctionReturn(0);
    }
 
    
    
     //Indices for rows only (used to extract lagrange multipliers from D vector to calculate Homogenized stress)
    double coords_face[9];
    virtual PetscErrorCode GetN_and_Indices() {
        PetscFunctionBegin;
        
        //For nodes
        row_mat = 0;
        RowGlob[row_mat].resize(9);
        typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
        const EntityHandle* conn_face;
        int num_nodes;
        EntityHandle face_tri;  face_tri=fe_ptr->get_ent();
        rval = moab.get_connectivity(face_tri,conn_face,num_nodes,true); CHKERR_PETSC(rval);
//        cout<< "num_nodes ="<<num_nodes << endl;
//        cout<< "conn_face ="<<conn_face << endl;
        rval = moab.get_coords(conn_face,num_nodes,coords_face); CHKERR_PETSC(rval);
//        for(int ii=0; ii<9; ii++) cout<<"coord "<<coords_face[ii]<<endl;
//        cout<<endl<<endl;
        ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
        ublas::vector<FieldData,ublas::bounded_array<double,3> > normal(3);
        ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
        area = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;   // area of each face of triangle
//        cout<<" area = "<<area<<endl;

        
        int nn = 0;
        for(;nn<3;nn++) {
            dofs_iterator niit,hi_niit;   //for rows
            
            //minimum and maximum row and column indices for each node on the surface
            niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
            hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
            
            for(;niit!=hi_niit;niit++) {
                RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
            }
        }
//        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<RowGlob[row_mat].size(); ii++) cout<<RowGlob[row_mat][ii]<<" ";
//        cout<<"\n\n\n"<<endl;
        
        
        // Find row indices for Edges
        vector<int> FaceEdgeSense;
        vector<int> FaceEdgeOrder;
        vector<vector<double> > N_edge_data;
        vector<vector<double> > diffN_edge_data;
        double* N_edge[3];
        double* diffN_edge[3];
        
        FaceEdgeSense.resize(3);
        FaceEdgeOrder.resize(3);
        N_edge_data.resize(3);
        diffN_edge_data.resize(3);
        int ee = 0; row_mat++;
        for(;ee<3;ee++) {
            EntityHandle edge;
            rval = moab.side_element(face_tri,1,ee,edge); CHKERR_PETSC(rval);
            int side_number,offset;
            rval = moab.side_number(face_tri,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
            dofs_iterator eiit,hi_eiit;
//            cout<<"side_number "<<side_number<<endl;
//            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
//            cout<<"edge "<<edge<<endl;
            eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            hi_eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            
            if(eiit!=hi_eiit) {
                FaceEdgeOrder[ee] = eiit->get_max_order();
                if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
                    assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
                    RowGlob[row_mat].resize(distance(eiit,hi_eiit));
//                    cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
                    for(;eiit!=hi_eiit;eiit++) {
                        RowGlob[row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
                    }
                    
                    N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    N_edge[ee] = &(N_edge_data[ee][0]);
                    diffN_edge[ee] = &(diffN_edge_data[ee][0]);
                    row_mat++;
                }
            }else {
                FaceEdgeOrder[ee] = 0;
                N_edge[ee] = NULL;
                diffN_edge[ee] = NULL;
            }
        }
//        cout<<"\nFor Edges "<<endl;
//        cout<<"row_mat  =  "<<row_mat<<endl;
//        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[1].size()<<endl;
//        for(int jj=0; jj<3; jj++) for(int ii=0; ii<RowGlob[1].size(); ii++) cout<<RowGlob[jj][ii]<<" ";
//        cout<<"\n\n\n";
        //Find the shape function N at each gauss point for all the edges and then re-araange in the form as mentioned for nodes
        if(N_edge[0] != NULL){
            ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
            //            cout<<"Hi from insidie"<<endl<<endl;
            ee = 0; int row_mat1=1;
            for (;ee<3;ee++,row_mat1++){
                rowNMatrices[row_mat1].resize(g_TRI_dim);
                unsigned int gg = 0; unsigned int gg1=0;
                int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
//                cout<<"nodes_edge  "<<nodes_edge<<endl;
                for(;gg<g_TRI_dim;gg++) {
                    rowNMatrices[row_mat1][gg].resize(3,3*nodes_edge);
                    rowNMatrices[row_mat1][gg].clear();
                    int kk=0;
                    for(int ii=0; ii<nodes_edge; ii++){
                        for (int jj=0; jj<3; jj++){
//                              cout<<"jj  "<<jj<<endl;
//                              cout<<"kk  "<<kk<<endl;
//                              cout<<"gg1  "<<gg1<<endl;
//                              cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
                            rowNMatrices[row_mat1][gg](jj,kk)=N_edge[ee][gg1]; kk++;
                        }
                        gg1++;
                    }
                }
            }
        }

        
        
        //Find the rows and column indices for face of the triangle
        dofs_iterator fiit,hi_fiit;
        fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",face_tri));
        hi_fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",face_tri));
        
        //cout<<"fiit!=hi_fiit  "<<(fiit!=hi_fiit) << endl;
        if(fiit!=hi_fiit) {
            RowGlob[row_mat].resize(distance(fiit,hi_fiit));
//            cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
            int face_order = fiit->get_max_order();
            assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
            if(NBFACE_H1(face_order)>0) {
                for(;fiit!=hi_fiit;fiit++) {
                    RowGlob[row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
                }
            }
            
//        cout<<"\nFor Faces "<<endl;
//        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[4].size()<<endl;
//        for(int ii=0; ii<RowGlob[4].size(); ii++) cout<<RowGlob[4][ii]<<" ";
//        cout<<endl<<endl;
            
            //Find out the shape funcition and rearrange in the form as shown for the nodes
//            cout<<"NBFACE_H1(face_order) "<<NBFACE_H1(face_order)<<endl;
            vector<double> N_face, diffN_face;
            N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
            diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
            int face_nodes[] = { 0,1,2 };
            ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
            
            rowNMatrices[row_mat].resize(g_TRI_dim);
            unsigned int gg = 0;  unsigned int gg1=0;
            int nodes_face=NBFACE_H1(face_order);
//            for(int ii=0; ii<N_face.size(); ii++) cout<<"N_face  "<<N_face[ii]<<endl;
            for(;gg<g_TRI_dim;gg++) {
                rowNMatrices[row_mat][gg].resize(3,3*nodes_face);
                rowNMatrices[row_mat][gg].clear();
                int kk=0;
                for(int ii=0; ii<nodes_face; ii++){
                    for (int jj=0; jj<3; jj++){
//                          cout<<"jj  "<<jj<<endl;
//                          cout<<"kk  "<<kk<<endl;
//                          cout<<"gg1  "<<gg1<<endl<<endl;
//                          cout<<"N_face[gg1]  "<<N_face[gg1]<<endl<<endl;
                          rowNMatrices[row_mat][gg](jj,kk)=N_face[gg1]; kk++;
                    }
                    gg1++;
                }
//                cout<<"rowNMatrices[row_mat][0](jj,kk)  "<<rowNMatrices[row_mat][gg]<<endl;
//                std::string wait;
//                std::cin >> wait;
            }
            row_mat++;
        }
        
        PetscFunctionReturn(0);
    }

    
    
    ublas::matrix<FieldData> X_mat, nodes_coord, gauss_coord;
    ublas::vector<ublas::matrix<FieldData> > D_mat;
    ublas::vector<FieldData>  Stress_Homo;
    ublas::vector<ublas::vector<FieldData> > Lamda;
    
    
    virtual PetscErrorCode Calculate_Homo_Stress() {
        PetscFunctionBegin;
        X_mat.resize(3,6);    X_mat.clear();
        nodes_coord.resize(3,3);
        gauss_coord.resize(3,g_TRI_dim);
        D_mat.resize(row_mat);
        Lamda.resize(row_mat);

        //used to calculate the coordinates of a Gauss points
        nodes_coord(0,0)=coords_face[0]; nodes_coord(0,1)=coords_face[3]; nodes_coord(0,2)=coords_face[6];
        nodes_coord(1,0)=coords_face[1]; nodes_coord(1,1)=coords_face[4]; nodes_coord(1,2)=coords_face[7];
        nodes_coord(2,0)=coords_face[2]; nodes_coord(2,1)=coords_face[5]; nodes_coord(2,2)=coords_face[8];
        
        //coordinates for all gauss points
        gauss_coord=prod(nodes_coord, g_NTRI_mat);
        
//        cout<<"g_NTRI_mat "<<g_NTRI_mat<<endl<<endl;
//        cout<<"nodes_coord "<<nodes_coord<<endl<<endl;
//        cout<<"gauss_coord "<<gauss_coord<<endl<<endl;
//        std::string wait;
//        std::cin >> wait;

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
                    D_mat[rr].resize(H_mat[rr].size1(),D_mat1.size2());
//                    cout<<"\n row_Mat "<<row_Mat<<endl;
//                    cout<<"\n col_Mat "<<rr<<col_Mat;
//                    cout<<"\n w "<<w<<col_Mat;
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    D_mat[rr]=prod(H_mat[rr], D_mat1);
                } else {
                    //calculate (D_mat1= w * NT * X_mat)
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    //calculate (D_mat = H_mat * D_mat1)
                    D_mat[rr]+=prod(H_mat[rr], D_mat1);
                }
            }
//            cout<<"row_mat  =  "<<row_mat<<endl;
//            cout<< " D_mat[rr] =  "<<D_mat[rr]<<endl;
            Lamda[rr].resize(RowGlob[rr].size());
            ierr = VecGetValues(DVec,RowGlob[rr].size(),&(RowGlob[rr])[0],&(Lamda[rr].data())[0]); CHKERRQ(ierr);
//            cout<< " Lamda[rr] =  "<<Lamda[rr]<<endl;
            
//            cout<< " Stress_Homo before =  "<<Stress_Homo<<endl;
//            cout<< "  New part          =  "<<prod(trans(D_mat[rr]), Lamda[rr])<<endl;
            Stress_Homo+=prod(trans(D_mat[rr]), -1*Lamda[rr]);   //Lamda is reaction force (so multiply for -1 to get the force)
//            cout<< " Stress_Homo after  =   "<<Stress_Homo<<endl;
        }
        
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
//        cout<<"Hi from class ElasticFE_RVELagrange_Homogenized_Stress"<<endl;
        
        
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
        ierr = Get_H_mat();   //It will be used from the Class ElasticFE_RVELagrange_Disp
        ierr = Calculate_Homo_Stress(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    
};

    
}

#endif //__ElasticFE_RVELagrange_Periodic_RigidBodyMotion
