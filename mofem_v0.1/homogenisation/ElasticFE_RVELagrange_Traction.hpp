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

#ifndef __ElasticFE_RVELagrange_Traction_HPP__
#define __ElasticFE_RVELagrange_Traction_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"


namespace MoFEM {

struct ElasticFE_RVELagrange_Traction: public ElasticFE_RVELagrange_Disp {

    ElasticFE_RVELagrange_Traction(
                          FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,ublas::vector<FieldData> _applied_strain):
    ElasticFE_RVELagrange_Disp(_mField, _dirihlet_ptr,_Aij, _D, _F, _applied_strain){};
    
    
    vector<DofIdx> DirihletBC;
    
    
    double coords_face[9];
    double area;
     virtual PetscErrorCode GetN_and_Indices() {
        PetscFunctionBegin;
        
        
        //Find out indices for row and column for nodes on the surface, i.e. triangles
        row_mat = 0;
        RowGlob[row_mat].resize(6);
        ColGlob[row_mat].resize(9);
        
        typedef FENumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator row_dofs_iterator;
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
         //cout<<" area = "<<area<<endl;
         
         
        //minimum and maximum rows indices for each node on the surface
        row_dofs_iterator niit,hi_niit;   //iterator for rows
        niit = row_multiIndex->get<FieldName_mi_tag>().lower_bound("Lagrange_mul_disp");
        hi_niit = row_multiIndex->get<FieldName_mi_tag>().upper_bound("Lagrange_mul_disp");
        int nn = 0;
        for(;niit!=hi_niit;niit++) {
            RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
        }

        
        nn = 0;
        for(;nn<3;nn++) {
            dofs_iterator col_niit,hi_col_niit;  // iterator for columns
            string field_name;
            
            //minimum and maximum row and column indices for each node on the surface
            col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
            hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
            
            // two different loops, i.e. one for row and one for column (may be need it for multiphysics problems)
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
//        //Stop code
//        std::string wait;
//        std::cin >> wait;

        
        // Find colum indices for Edges (row are always 6)
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
            dofs_iterator col_eiit,col_hi_eiit;
            
//            cout<<"side_number "<<side_number<<endl;
//            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
//            cout<<"edge "<<edge<<endl<<endl;
            
            col_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",edge));
            col_hi_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",edge));
            
            if(col_eiit!=col_hi_eiit) {
//                cout<<"Hello "<<endl;
                FaceEdgeOrder[ee] = col_eiit->get_max_order();
                if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
                    assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(col_eiit,col_hi_eiit));
                    ColGlob[row_mat].resize(distance(col_eiit,col_hi_eiit));
//                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
                    for(;col_eiit!=col_hi_eiit;col_eiit++) {
                        ColGlob[row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
                    }
                    
                    //shape functions for edges
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
        
//        cout<<"row_mat  =  "<<row_mat<<endl;
//        cout<<"\nFor Edges "<<endl;
//        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[1].size()<<endl;
//        for(int jj=1; jj<4; jj++) for(int ii=0; ii<ColGlob[1].size(); ii++) cout<<ColGlob[jj][ii]<<" ";
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
//                    cout<<"rowNMatrices "<<rowNMatrices[row_mat1][gg]<<endl;
                }
            }
        }
        
        //Find the column indices for face of the triangle (indices for the rows are the same as for nodes)
        dofs_iterator col_fiit, col_hi_fiit;
        col_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",face_tri));
        col_hi_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",face_tri));
        
        if(col_fiit!=col_hi_fiit) {
            ColGlob[row_mat].resize(distance(col_fiit,col_hi_fiit));
//            cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//            cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
            int face_order = col_fiit->get_max_order();
            assert((unsigned int)3*NBFACE_H1(face_order)==distance(col_fiit,col_hi_fiit));
            if(NBFACE_H1(face_order)>0) {
                for(;col_fiit!=col_hi_fiit;col_fiit++) {
                    ColGlob[row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
                }
            }
//        cout<<"row_mat  =  "<<row_mat<<endl;
//        cout<<"\nFor Faces "<<endl;
//        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
//        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[row_mat][ii]<<" ";
//        cout<<"\n\n\n";
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

    
    ublas::vector<ublas::matrix<FieldData> > H_mat;
    
    virtual PetscErrorCode Get_H_mat() {
        PetscFunctionBegin;
        
        double coords_face[9];
        //To find the side of the element (+X, -X, +Y, -Y, +Z, -Z)
        const EntityHandle *conn_face;
        int num_nodes_face;
        rval = moab.get_connectivity(fe_ptr->get_ent(),conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
//        cout<<"num_nodes_face "<<num_nodes_face<<endl<<endl;
        rval = moab.get_coords(conn_face,num_nodes_face,coords_face); CHKERR_PETSC(rval);
//        for(int ii=0; ii<9; ii++) cout<<"coord "<<coords_face[ii]<<"   ";
//        cout<<endl<<endl;
        ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
        ublas::vector<FieldData,ublas::bounded_array<double,3> > normal(3);
        ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
        
        
//        cout<<"normal(0) "<<normal(0)<<endl;
//        cout<<"normal(1) "<<normal(1)<<endl;
//        cout<<"normal(2) "<<normal(2)<<endl;
//        std::string wait;
//        std::cin >> wait;
        
        ublas::matrix<FieldData> H_mat_1Node;
        H_mat_1Node.resize(6,3);  H_mat_1Node.clear();
        
        if(normal(0)>0){     //+X face of the RVE
            H_mat_1Node(0,0)=1.0;  H_mat_1Node(3,1)=1.0;  H_mat_1Node(4,2)=1.0;
        }

        if(normal(0)<0){    //-X face of the RVE
            H_mat_1Node(0,0)=-1.0;  H_mat_1Node(3,1)=-1.0;  H_mat_1Node(4,2)=-1.0;
        }
        
        if(normal(1)>0){     //+Y face of the RVE
            H_mat_1Node(1,1)=1.0;  H_mat_1Node(3,0)=1.0;  H_mat_1Node(5,2)=1.0;
        }
        
        if(normal(1)<0){    //-Y face of the RVE
            H_mat_1Node(1,1)=-1.0;  H_mat_1Node(3,0)=-1.0;  H_mat_1Node(5,2)=-1.0;
        }
        
        if(normal(2)>0){    //+Z face of the RVE
            H_mat_1Node(2,2)=1.0;  H_mat_1Node(4,0)=1.0;  H_mat_1Node(5,1)=1.0;
        }
        
        if(normal(2)<0){    //-Z face of the RVE
            H_mat_1Node(2,2)=-1.0;  H_mat_1Node(4,0)=-1.0;  H_mat_1Node(5,1)=-1.0;
        }
        
        H_mat.resize(row_mat);
//        cout<<"normal "<<normal<<endl;
        for(int rrr=0; rrr<row_mat; rrr++){ //for row_mat
//            cout<<"(rowNMatrices[rr])[0] "<<(rowNMatrices[rrr])[0]<<endl;
            int num_col=(rowNMatrices[rrr])[0].size2();
            H_mat[rrr].resize(6,num_col);
            int cc1=0;
            for(int bb = 0; bb<num_col/3; bb++) {  //blocks of 6x3
//                cout<<"cc1 "<<cc1<<endl;
                for(int rr = 0; rr<6; rr++) {
//                     cout<<"rr "<<rr<<endl;
                    for(int cc = 0; cc<3; cc++) {
//                        cout<<"cc "<<cc<<endl;
//                        cout<<"H_mat_1Node(rr,cc)"<<H_mat_1Node(rr,cc)<<endl<<endl;
                        H_mat[rrr](rr,(cc+cc1))=H_mat_1Node(rr,cc);
                    }
                }
                cc1+=3;
            }
//            cout<<"H_mat[rrr] "<<H_mat[rrr]<<endl;
        }
//        std::string wait;
//        std::cin >> wait;
        PetscFunctionReturn(0);
    }

    
    
    
    //Calculate and assemble NT x N matrix
//    ublas::matrix<ublas::matrix<FieldData> > NTN;
    virtual PetscErrorCode Stiffness() {
        PetscFunctionBegin;
        //        cout<<" row_mat; "<<row_mat<<endl;
        NTN.resize(1,row_mat);
        const EntityHandle *conn_face;
        int num_nodes_face;
        rval = moab.get_connectivity(fe_ptr->get_ent(),conn_face,num_nodes_face,true); CHKERR_PETSC(rval);
        //        cout<<"num_nodes_face "<<num_nodes_face<<endl<<endl;
        
        //Calculate C Matrix, i.e. (NT x N)
        for(int rr = 0;rr<1;rr++) {   //only for nodes (i.e. only 6 rows we have no row indices for higher order DOFs)
            //            cout<<" rr = "<<rr<<endl;
            for(int gg = 0;gg<g_TRI_dim;gg++) {
                ublas::matrix<double> &row_Mat = (rowNMatrices[rr])[gg];
//                cout<<" row_Mat; "<<row_Mat<<endl;
                double w = area*G_W_TRI[gg];
                for(int cc = 0;cc<row_mat;cc++) {
                    ublas::matrix<FieldData> &col_Mat = (rowNMatrices[cc])[gg];
                    ublas::matrix<FieldData>  NTN1;
                    NTN1.resize(row_Mat.size2(),col_Mat.size2());

                    if(gg == 0) {
//                        cout<<" row_Mat; "<<row_Mat<<endl;
//                        cout<<" col_Mat; "<<col_Mat<<endl;
                        
                        //calculate NTN1=(w* NT * N)
                        NTN1 = prod( w*trans(row_Mat), col_Mat);
                        //calculate NTN=(H_mat* w* NT * N)
                        NTN(rr,cc).resize(H_mat[rr].size1(),NTN1.size2());
                        NTN(rr,cc) = prod(H_mat[rr], NTN1);
//                        cout<<" NTN= "<<NTN(rr,cc)<<endl;
                    }else {
                        NTN1 = prod( w*trans(row_Mat), col_Mat);
                        NTN(rr,cc)+=prod(H_mat[rr], NTN1);
                     }
                }
            }
        }
        PetscFunctionReturn(0);
    }
   
    
    
    virtual PetscErrorCode Lhs() {
        PetscFunctionBegin;
        ierr = Stiffness(); CHKERRQ(ierr);
        
        //Assembly C with size (6 x 3N), where M is number of nodes on boundary and N are total nodes
        for(int rr = 0;rr<1;rr++) { //only for nodes (i.e. only 6 rows, we have no row indices for higher order DOFs)
            for(int cc = 0;cc<row_mat;cc++) {
//                cout<<"\n RowGlob[rr].size() "<<RowGlob[rr].size()<<endl;
//                for(int ii=0; ii<RowGlob[rr].size(); ii++) cout<<RowGlob[rr][ii]<<" ";
//                cout<<"\n ColGlob[cc].size() "<<ColGlob[cc].size()<<endl;
//                for(int ii=0; ii<ColGlob[cc].size(); ii++) cout<<ColGlob[cc][ii]<<" ";
//                cout<<"\n\n NTN(rr,cc).size() "<<NTN(rr,cc)<<endl;
//                cout<<"\n\n\n"<<endl;
                if(ColGlob[cc].size()==0) continue;
                if(RowGlob[rr].size()!=NTN(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
                if(ColGlob[cc].size()!=NTN(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

                ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(NTN(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
            }
        }


        //Assembly trans(C) with size (3Nx6), where M is number of nodes on boundary and N are total nodes
        for(int rr = 0;rr<1;rr++) {
            for(int cc = 0;cc<row_mat;cc++) {
                
//                cout<<"\n\nNTN(rr,cc) "<<NTN(rr,cc)<<endl;
                NTN(rr,cc) = trans(NTN(rr,cc));
//                cout<<"\n\nNTN(rr,cc) "<<NTN(rr,cc)<<endl;
//                std::string wait;
//                std::cin >> wait;

                ierr = MatSetValues(Aij,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(NTN(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
            }
        }
        PetscFunctionReturn(0);
    }

    
    
    
    ublas::matrix<FieldData> X_mat, nodes_coord, gauss_coord;
    ublas::vector<ublas::matrix<FieldData> > D_mat;
    ublas::vector<FieldData> f; //f.resize(9);
    
    //Calculate the right hand side vector, i.e. f=D_mat * applied_strain and assemble it into the global force vector F
    virtual PetscErrorCode Rhs() {
        PetscFunctionBegin;
        X_mat.resize(3,6);    X_mat.clear();
        nodes_coord.resize(3,3);
        gauss_coord.resize(3,g_TRI_dim);
        D_mat.resize(row_mat);
        
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
        //cout<<"area "<<area << endl;
        
        
        for(int rr=0; rr<1; rr++){
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
                    
                    //calculate (D_mat1= w * NT * X_mat)
//                    cout<<"\n row_Mat "<<row_Mat<<endl;
//                    cout<<"\n col_Mat "<<rr<<col_Mat;
//                    cout<<"\n w "<<w<<col_Mat;
                    
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    //calculate (D_mat = H_mat * D_mat1)
//                    cout<<"\n rr "<<rr<<endl;
//                    cout<<"\n D_mat[rr] "<<D_mat[rr]<<endl;
                    D_mat[rr]=prod(H_mat[rr], D_mat1);
                } else {
                    //calculate (D_mat1= w * NT * X_mat)
                    D_mat1=prod(w*trans(row_Mat), col_Mat);
                    
                    //calculate (D_mat = H_mat * D_mat1)
                    D_mat[rr]+=prod(H_mat[rr], D_mat1);
                }
            }
//            cout<<"\n rr "<<rr<<endl;
//            cout<<"\n RowGlob[0].size() "<<RowGlob[0].size()<<endl;
//            for(int ii=0; ii<RowGlob[0].size(); ii++) cout<<RowGlob[0][ii]<<" ";
//            cout<< "\nD_mat[rr] =  "<<D_mat[rr]<<endl<<endl;
            //Assemble D_mat into global force vector F
//            cout<<"\n D_mat[rr] "<<D_mat[rr]<<endl;
            f=prod(D_mat[rr], applied_strain);
//            cout<<"\n f "<<f<<endl;
            ierr = VecSetValues(F,RowGlob[0].size(),&(RowGlob[0])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
        }
        
        PetscFunctionReturn(0);
    }

    
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
//        cout<<"Hi from class"<<endl;
        
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
        ierr = Get_H_mat();
        
//        cout<<"Before BCs "<<endl;
//        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndiciesCol(this,ColGlob,DirihletBC); CHKERRQ(ierr);
//        cout<<"After BCs "<<endl;
        
        ierr = Lhs(); CHKERRQ(ierr);
        ierr = Rhs(); CHKERRQ(ierr);
//        std::string wait;
//        std::cin >> wait;
        PetscFunctionReturn(0);
    }

   

};

    
}

#endif //__ElasticFE_RVELagrange_Traction_HPP__
