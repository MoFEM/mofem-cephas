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
        
        RowGlob.clear();  rowNMatrices.clear();  ColGlob.clear();
        RowGlob.resize(1+6+2);    // 1-node, 6-edges   2-face (Two opposite triangles for prisms)
        rowNMatrices.resize(1+6+2);
        ColGlob.resize(1+6+2);
        
        
        //Indices for row and column for nodes for Prisms element
        row_mat = 0;  RowGlob[row_mat].clear();   ColGlob[row_mat].clear();
        RowGlob[row_mat].resize(18);
        ColGlob[row_mat].resize(18);
        typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
        const EntityHandle* conn_prism;
        int num_nodes;
        EntityHandle prism_periodic;  prism_periodic=fe_ptr->get_ent();
        rval = moab.get_connectivity(prism_periodic,conn_prism,num_nodes,true); CHKERR_PETSC(rval);
//        cout<< "num_nodes ="<<num_nodes << endl;
        
        int nn = 0;
        for(;nn<num_nodes;nn++) {
            dofs_iterator niit,hi_niit;   //for rows
            dofs_iterator col_niit,hi_col_niit;  // for columns
            string field_name;
            
            //minimum and maximum row and column indices for each node on the surface
            niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",conn_prism[nn]));
            hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",conn_prism[nn]));
            
            col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_prism[nn]));
            hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_prism[nn]));
        
            for(;niit!=hi_niit;niit++) {
                RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
            }
            for(;col_niit!=hi_col_niit;col_niit++) {
                ColGlob[row_mat][nn*col_niit->get_max_rank()+col_niit->get_dof_rank()] = col_niit->get_petsc_gloabl_dof_idx();
            }
        }
//        cout<<"\nFor nodes "<<endl;
//        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<RowGlob[row_mat].size(); ii++) cout<<RowGlob[row_mat][ii]<<" ";
//        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[row_mat][ii]<<" ";
//        cout<<"\n\n\n"<<endl;
        
        // Find row and colum indices for Edges
        vector<int> FaceEdgeSense;
        vector<int> FaceEdgeOrder;
        vector<vector<double> > N_edge_data;
        vector<vector<double> > diffN_edge_data;
        double* N_edge[6];
        double* diffN_edge[6];
        
        FaceEdgeSense.resize(6);
        FaceEdgeOrder.resize(6);
        N_edge_data.resize(6);
        diffN_edge_data.resize(6);

        int ee_arr[]={0, 1, 2, 6, 7, 8};  //side numbers (Canonical numbering) of edges belong to triangles
        int ee = 0; row_mat++;
        for(;ee<6;ee++) {
            EntityHandle edge;
            rval = moab.side_element(prism_periodic,1,ee_arr[ee],edge); CHKERR_PETSC(rval);
            int side_number,offset;
            rval = moab.side_number(prism_periodic,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
            dofs_iterator eiit,hi_eiit,col_eiit,col_hi_eiit;
//            cout<<"side_number "<<side_number<<endl;
//            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
//            cout<<"edge "<<edge<<endl;
            
            eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            hi_eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            col_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",edge));
            col_hi_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",edge));
            
            //For rows (we will get indices only for -ve edges, i.e. eiit!=hi_eiit===True)
            RowGlob[row_mat].clear();  //without clear it will get the previous value
            if(eiit!=hi_eiit) {
                FaceEdgeOrder[ee] = eiit->get_max_order();
                if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
                    assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
                    RowGlob[row_mat].resize(distance(eiit,hi_eiit));
//                    cout<<"row_mat "<<row_mat<<endl;
//                    cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
                    for(;eiit!=hi_eiit;eiit++) {
                        RowGlob[row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
                    }
                }
            }

            ColGlob[row_mat].clear(); //without clear it will get the previous value
            if(col_eiit!=col_hi_eiit) {
                FaceEdgeOrder[ee] = col_eiit->get_max_order();
                if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
                    assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(col_eiit,col_hi_eiit));
                    ColGlob[row_mat].resize(distance(col_eiit,col_hi_eiit));
//                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
                    for(;col_eiit!=col_hi_eiit;col_eiit++) {
                        ColGlob[row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
                    }
                    N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    N_edge[ee] = &(N_edge_data[ee][0]);
                    diffN_edge[ee] = &(diffN_edge_data[ee][0]);
                    row_mat++;   //row_mat++ here as column indices exists for every edge (this increment will create 0 entries in row indices)
                }
            }
        }

//        cout<<"\nFor Edges "<<endl;
//        cout<<"\nRowGlob[ii].size() "<<endl;
//        for(int ii=1; ii<row_mat; ii++) cout<<RowGlob[ii].size() <<"  ";
//        cout<<"\nColGlob[ii].size() "<<endl;
//        for(int ii=1; ii<row_mat; ii++) cout<<ColGlob[ii].size() <<"  ";
//        cout<<"\nCol indices "<<endl;
//        for(int jj=1; jj<row_mat; jj++) for(int ii=0; ii<ColGlob[jj].size(); ii++) cout<<ColGlob[jj][ii]<<" "; cout<<endl;
//        cout<<"Row indices "<<endl;
//        for(int jj=1; jj<row_mat; jj++) for(int ii=0; ii<RowGlob[jj].size(); ii++) cout<<RowGlob[jj][ii]<<" "; cout<<endl;
//        cout<<"\n\n\n";
        
//        //Find the shape function N at each gauss point for all the edges and then re-araange in the form as mentioned for nodes
//        if(N_edge[0] != NULL){
//            ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
//            //            cout<<"Hi from insidie"<<endl<<endl;
//            ee = 0; int row_mat1=1;
//            for (;ee<3;ee++,row_mat1++){
//                rowNMatrices[row_mat1].resize(g_TRI_dim);
//                unsigned int gg = 0; unsigned int gg1=0;
//                int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
//                //                cout<<"nodes_edge  "<<nodes_edge<<endl;
//                for(;gg<g_TRI_dim;gg++) {
//                    rowNMatrices[row_mat1][gg].resize(3,3*nodes_edge);
//                    rowNMatrices[row_mat1][gg].clear();
//                    int kk=0;
//                    for(int ii=0; ii<nodes_edge; ii++){
//                        for (int jj=0; jj<3; jj++){
//                            //                              cout<<"jj  "<<jj<<endl;
//                            //                              cout<<"kk  "<<kk<<endl;
//                            //                              cout<<"gg1  "<<gg1<<endl;
//                            //                              cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
//                            rowNMatrices[row_mat1][gg](jj,kk)=N_edge[ee][gg1]; kk++;
//                        }
//                        gg1++;
//                    }
//                }
//            }
//        }
        
        
        
        
        vector<int> FaceOrder;
        FaceOrder.resize(2);

        int ff_arr[]={3, 4};  //Canonical numbering of two triangles only
        int ff = 0; row_mat++;
        for(;ff<2;ff++) {
            EntityHandle Tri_Prism;
            rval = moab.side_element(prism_periodic,2,ff_arr[ff],Tri_Prism); CHKERR_PETSC(rval);
            
            dofs_iterator fiit,hi_fiit, col_fiit, col_hi_fiit;
            fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",Tri_Prism));
            hi_fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",Tri_Prism));
            
            col_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",Tri_Prism));
            col_hi_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",Tri_Prism));

             //For rows (we will get indices only for -ve edges, i.e. eiit!=hi_eiit===True)
            RowGlob[row_mat].clear();  //without clear it will get the previous value
            if(fiit!=hi_fiit) {
                FaceOrder[ff] = fiit->get_max_order();
//                cout<<"FaceOrder[ff] "<<FaceOrder[ff]<<endl;
                assert((unsigned int)3*NBFACE_H1(FaceOrder[ff])==distance(fiit,hi_fiit));
                if(NBFACE_H1(FaceOrder[ff])>0) {
                    RowGlob[row_mat].resize(distance(fiit,hi_fiit));
//                    cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
                    for(;fiit!=hi_fiit;fiit++) {
                        RowGlob[row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
                    }
                }
            }
            
            ColGlob[row_mat].clear(); //without clear it will get the previous value
            if(col_fiit!=col_hi_fiit) {
                FaceOrder[ff] = col_fiit->get_max_order();
//                cout<<"FaceOrder[ff] "<<FaceOrder[ff]<<endl;
                assert((unsigned int)3*NBFACE_H1(FaceOrder[ff])==distance(col_fiit,col_hi_fiit));
                if(NBFACE_H1(FaceOrder[ff])>0) {
                    ColGlob[row_mat].resize(distance(col_fiit,col_hi_fiit));
//                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
                    for(;col_fiit!=col_hi_fiit;col_fiit++) {
                        ColGlob[row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
                    }
                }
            }
            

            
            cout<<"\n\n";
        }
        
        
        
        
        
        
        
//        //Stop code
//        std::string wait;
//        std::cin >> wait;
//
        
        PetscFunctionReturn(0);
    }
    
    
    
    
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
