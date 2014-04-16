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

#ifndef __ElasticFE_RVELagrange_Periodic_HPP__
#define __ElasticFE_RVELagrange_Periodic_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange.hpp"

namespace MoFEM {

struct ElasticFE_RVELagrange_Periodic: public ElasticFE_RVELagrange {

    ElasticFE_RVELagrange_Periodic(
                                   FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F):
    ElasticFE_RVELagrange(_mField, _dirihlet_ptr,_Aij, _D, _F ){};
    

    
    
    
    virtual PetscErrorCode GetN_and_Indices() {
        PetscFunctionBegin;
        
        
        //Find out indices for row and column for nodes on the surface, i.e. triangles
        row_mat = 0;
        RowGlob[row_mat].resize(9);
        ColGlob[row_mat].resize(9);
        typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
        const EntityHandle* conn_face;
        int num_nodes;
        EntityHandle face_tri;  face_tri=fe_ptr->get_ent();
        rval = moab.get_connectivity(face_tri,conn_face,num_nodes,true); CHKERR_PETSC(rval);
        cout<< "num_nodes ="<<num_nodes << endl;
        cout<< "conn_face ="<<conn_face << endl;
        //Stop code
        std::string wait;
        std::cin >> wait;

//        int nn = 0;
//        for(;nn<3;nn++) {
//            dofs_iterator niit,hi_niit;   //for rows
//            dofs_iterator col_niit,hi_col_niit;  // for columns
//            string field_name;
//            
//            //minimum and maximum row and column indices for each node on the surface
//            niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
//            hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
//            col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
//            hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
//            
//            
//            // two different loops, i.e. one for row and one for column (may be need it for multiphysics problems)
//            for(;niit!=hi_niit;niit++) {
//                RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
//            }
//            for(;col_niit!=hi_col_niit;col_niit++) {
//                ColGlob[row_mat][nn*col_niit->get_max_rank()+col_niit->get_dof_rank()] = col_niit->get_petsc_gloabl_dof_idx();
//            }
//        }
//        
//        
//        cout<<"\nFor nodes "<<endl;
//        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<RowGlob[row_mat].size(); ii++) cout<<RowGlob[row_mat][ii]<<" ";
//        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[row_mat][ii]<<" ";
//        cout<<"\n\n\n"<<endl;
        
        PetscFunctionReturn(0);
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
//    // Find colum indices for Edges (row are always 6)
//    vector<int> FaceEdgeSense;
//    vector<int> FaceEdgeOrder;
//    vector<vector<double> > N_edge_data;
//    vector<vector<double> > diffN_edge_data;
//    double* N_edge[3];
//    double* diffN_edge[3];
//    
//    FaceEdgeSense.resize(3);
//    FaceEdgeOrder.resize(3);
//    N_edge_data.resize(3);
//    diffN_edge_data.resize(3);
//    int ee = 0; row_mat++;
//    for(;ee<3;ee++) {
//        EntityHandle edge;
//        rval = moab.side_element(face_tri,1,ee,edge); CHKERR_PETSC(rval);
//        int side_number,offset;
//        rval = moab.side_number(face_tri,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
//        dofs_iterator col_eiit,col_hi_eiit;
//        
//        //            cout<<"side_number "<<side_number<<endl;
//        //            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
//        //            cout<<"edge "<<edge<<endl<<endl;
//        
//        col_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",edge));
//        col_hi_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",edge));
//        
//        if(col_eiit!=col_hi_eiit) {
//            //                cout<<"Hello "<<endl;
//            FaceEdgeOrder[ee] = col_eiit->get_max_order();
//            if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
//                assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(col_eiit,col_hi_eiit));
//                ColGlob[row_mat].resize(distance(col_eiit,col_hi_eiit));
//                //                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//                for(;col_eiit!=col_hi_eiit;col_eiit++) {
//                    ColGlob[row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
//                }
//                
//                //shape functions for edges
//                N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
//                diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
//                N_edge[ee] = &(N_edge_data[ee][0]);
//                diffN_edge[ee] = &(diffN_edge_data[ee][0]);
//                row_mat++;
//            }
//        }else {
//            FaceEdgeOrder[ee] = 0;
//            N_edge[ee] = NULL;
//            diffN_edge[ee] = NULL;
//        }
//    }
//    
//    //        cout<<"row_mat  =  "<<row_mat<<endl;
//    //        cout<<"\nFor Edges "<<endl;
//    //        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[1].size()<<endl;
//    //        for(int jj=1; jj<4; jj++) for(int ii=0; ii<ColGlob[1].size(); ii++) cout<<ColGlob[jj][ii]<<" ";
//    //        cout<<"\n\n\n";
//    
//    //Find the shape function N at each gauss point for all the edges and then re-araange in the form as mentioned for nodes
//    if(N_edge[0] != NULL){
//        ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
//        //            cout<<"Hi from insidie"<<endl<<endl;
//        ee = 0; int row_mat1=1;
//        for (;ee<3;ee++,row_mat1++){
//            rowNMatrices[row_mat1].resize(g_TRI_dim);
//            unsigned int gg = 0; unsigned int gg1=0;
//            int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
//            //                cout<<"nodes_edge  "<<nodes_edge<<endl;
//            for(;gg<g_TRI_dim;gg++) {
//                rowNMatrices[row_mat1][gg].resize(3,3*nodes_edge);
//                rowNMatrices[row_mat1][gg].clear();
//                int kk=0;
//                for(int ii=0; ii<nodes_edge; ii++){
//                    for (int jj=0; jj<3; jj++){
//                        //                              cout<<"jj  "<<jj<<endl;
//                        //                              cout<<"kk  "<<kk<<endl;
//                        //                              cout<<"gg1  "<<gg1<<endl;
//                        //                              cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
//                        rowNMatrices[row_mat1][gg](jj,kk)=N_edge[ee][gg1]; kk++;
//                    }
//                    gg1++;
//                }
//                //                    cout<<"rowNMatrices "<<rowNMatrices[row_mat1][gg]<<endl;
//            }
//        }
//    }
//    
//    //Find the column indices for face of the triangle (indices for the rows are the same as for nodes)
//    dofs_iterator col_fiit, col_hi_fiit;
//    col_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",face_tri));
//    col_hi_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",face_tri));
//    
//    if(col_fiit!=col_hi_fiit) {
//        ColGlob[row_mat].resize(distance(col_fiit,col_hi_fiit));
//        //            cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//        //            cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//        int face_order = col_fiit->get_max_order();
//        assert((unsigned int)3*NBFACE_H1(face_order)==distance(col_fiit,col_hi_fiit));
//        if(NBFACE_H1(face_order)>0) {
//            for(;col_fiit!=col_hi_fiit;col_fiit++) {
//                ColGlob[row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
//            }
//        }
//        //        cout<<"row_mat  =  "<<row_mat<<endl;
//        //        cout<<"\nFor Faces "<<endl;
//        //        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//        //        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[row_mat][ii]<<" ";
//        //        cout<<"\n\n\n";
//        //Find out the shape funcition and rearrange in the form as shown for the nodes
//        //            cout<<"NBFACE_H1(face_order) "<<NBFACE_H1(face_order)<<endl;
//        vector<double> N_face, diffN_face;
//        N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
//        diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
//        int face_nodes[] = { 0,1,2 };
//        ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
//        
//        rowNMatrices[row_mat].resize(g_TRI_dim);
//        unsigned int gg = 0;  unsigned int gg1=0;
//        int nodes_face=NBFACE_H1(face_order);
//        //            for(int ii=0; ii<N_face.size(); ii++) cout<<"N_face  "<<N_face[ii]<<endl;
//        for(;gg<g_TRI_dim;gg++) {
//            rowNMatrices[row_mat][gg].resize(3,3*nodes_face);
//            rowNMatrices[row_mat][gg].clear();
//            int kk=0;
//            for(int ii=0; ii<nodes_face; ii++){
//                for (int jj=0; jj<3; jj++){
//                    //                          cout<<"jj  "<<jj<<endl;
//                    //                          cout<<"kk  "<<kk<<endl;
//                    //                          cout<<"gg1  "<<gg1<<endl<<endl;
//                    //                          cout<<"N_face[gg1]  "<<N_face[gg1]<<endl<<endl;
//                    rowNMatrices[row_mat][gg](jj,kk)=N_face[gg1]; kk++;
//                }
//                gg1++;
//            }
//            //                cout<<"rowNMatrices[row_mat][0](jj,kk)  "<<rowNMatrices[row_mat][gg]<<endl;
//            //                std::string wait;
//            //                std::cin >> wait;
//        }
//        row_mat++;
//    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      PetscErrorCode operator()() {
      PetscFunctionBegin;
        cout<<"Hi from class"<<endl;
        
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
//        ierr = Get_H_mat();
//        ierr = Lhs(); CHKERRQ(ierr);
//        ierr = Rhs(); CHKERRQ(ierr);
          
        PetscFunctionReturn(0);
    }

};

    
}

#endif //__ElasticFE_RVELagrange_Periodic__
