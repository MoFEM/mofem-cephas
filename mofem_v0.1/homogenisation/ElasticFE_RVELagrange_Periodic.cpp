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

#include "ElasticFE_RVELagrange_Periodic.hpp"

namespace MoFEM {
  
  PetscErrorCode ElasticFE_RVELagrange_Periodic::GetN_and_Indices() {
    PetscFunctionBegin;
    
    RowGlob.resize(2); ColGlob.resize(2);
    row_mat=0;
    RowGlob[0].resize(1+3+1);    // 1-node, 3-edges   1-face (-ve faces of prisms)
    ColGlob[0].resize(1+3+1);    // 1-node, 3-edges   1-face (-ve faces of prisms)
    RowGlob[1].resize(1+3+1);    // 1-node, 3-edges   1-face (+ve faces of prisms)
    ColGlob[1].resize(1+3+1);    // 1-node, 3-edges   1-face (+ve faces of prisms)
    rowNMatrices.resize(1+3+1);  // shape functions are the same for +ve and negative triangles (so only calculate for one face)
    
    typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
    int /*num_nodes,*/ num_nodes1;
    dofs_iterator niit,hi_niit;   //for rows
    dofs_iterator col_niit,hi_col_niit;  // for columns
    prism_periodic=fePtr->get_ent();
    //Indices for row and column for nodes for Prisms element
    
    //        cout<<"row_mat "<<row_mat<<endl;
    RowGlob[0][0].clear();     ColGlob[0][0].clear();    RowGlob[1][0].clear();     ColGlob[1][0].clear();
    RowGlob[0][0].resize(9);   ColGlob[0][0].resize(9);  RowGlob[1][0].resize(9);   ColGlob[1][0].resize(9);
    
    rval = moab.get_connectivity(prism_periodic,conn_Prism,num_nodes1,true); CHKERR_PETSC(rval);
    rval = moab.get_coords(conn_Prism,num_nodes1,coords_prism); CHKERR_PETSC(rval);
    
    //        cout<<"num_nodes1 "<<num_nodes1<<endl;
    //        cout<<"coord   ";
    //        for(int ii=0; ii<18; ii++) cout<<coords_prism[ii]<<"  ";
    //        cout<<endl<<endl;
    
    int count=0;
    for(int ff=0; ff<2; ff++) {
      for(int nn = 0;nn<num_nodes1/2;nn++) {
        coords_face[ff][3*nn+0]=coords_prism[3*count+0];    coords_face[ff][3*nn+1]=coords_prism[3*count+1];    coords_face[ff][3*nn+2]=coords_prism[3*count+2];
        ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
        ublas::vector<FieldData,ublas::bounded_array<double,3> > normal(3);
        ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face[ff],&*normal.data().begin()); CHKERRQ(ierr);
        area = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;   // area of each face of triangle
        
        //minimum and maximum row and column indices for each node on the prism triangle
        niit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",conn_Prism[count]));
        hi_niit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",conn_Prism[count]));
        
        col_niit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_Prism[count]));
        hi_col_niit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_Prism[count]));
        
        if(ff==0){
          for(;niit!=hi_niit;niit++) {
            RowGlob[0][0][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
            RowGlob[1][0][nn*niit->get_max_rank()+niit->get_dof_rank()]  = niit->get_petsc_gloabl_dof_idx();
          }}
        for(;col_niit!=hi_col_niit;col_niit++) {
          ColGlob[ff][0][nn*col_niit->get_max_rank()+col_niit->get_dof_rank()] = col_niit->get_petsc_gloabl_dof_idx();
        }
        count++;
      }
    }
    row_mat++;
    
    
    //        cout<<"\nFor Nodes "<<endl;
    //        cout<<"\nRowGlob[ii].size() "<<endl;
    //        for(int ii=0; ii<2; ii++) cout<<RowGlob[ii][0].size() <<"  ";
    //        cout<<"\nColGlob[ii].size() "<<endl;
    //        for(int ii=0; ii<2; ii++) cout<<ColGlob[ii][0].size() <<"  ";
    //        cout<<endl;
    //        for(int jj=0; jj<2; jj++) for(int ii=0; ii<RowGlob[jj][0].size(); ii++) cout<<RowGlob[jj][0][ii]<<" ";
    //        cout<<endl;
    //        for(int jj=0; jj<2; jj++) for(int ii=0; ii<ColGlob[jj][0].size(); ii++) cout<<ColGlob[jj][0][ii]<<" ";
    //        cout<<"\n\n\n"<<endl;
    
    //        cout<<"coord   ";
    //        for(int ff=0; ff<2; ff++)for(int ii=0; ii<9; ii++) cout<<coords_face[ff][ii]<<"  ";
    //        cout<<endl<<endl;
    //        //Stop code
    //        std::string wait;
    //        std::cin >> wait;
    
    // Find row and colum indices for Edges
    vector<int> FaceEdgeSense;
    vector<int> FaceEdgeOrder;
    vector<vector<double> > N_edge_data;
    vector<vector<double> > diffN_edge_data;
    
    
    double* N_edge[3];  //Only for one side of prism (Triangle)
    double* diffN_edge[3];
    FaceEdgeSense.resize(3);
    FaceEdgeOrder.resize(3);
    N_edge_data.resize(3);
    diffN_edge_data.resize(3);
    
    count=0;
    int ee_arr[]={0, 1, 2, 6, 7, 8};  //side numbers (Canonical numbering) of edges belong to triangles
    
    for(int ff=0; ff<2; ff++) { //face
      row_mat=1;
      for(int ee=0; ee<3; ee++){ //edge
        EntityHandle edge;
        rval = moab.side_element(prism_periodic,1,ee_arr[count],edge); CHKERR_PETSC(rval);
        int side_number,offset;
        rval = moab.side_number(prism_periodic,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
        
        //                cout<<"side_number "<<side_number<<endl;
        //                cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
        //                cout<<"edge "<<edge<<endl;
        //                std::string wait;
        //                std::cin >> wait;
        
        dofs_iterator eiit,hi_eiit,col_eiit,col_hi_eiit;
        eiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",edge));
        hi_eiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",edge));
        col_eiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",edge));
        col_hi_eiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",edge));
        
        //For rows (we will get indices only for -ve edges, i.e. eiit!=hi_eiit===True)
        if(ff==0){
          RowGlob[0][row_mat].clear();  RowGlob[1][row_mat].clear();//without clear it will get the previous value
          if(eiit!=hi_eiit) {
            FaceEdgeOrder[ee] = eiit->get_max_order();
            if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
              assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
              
              RowGlob[0][row_mat].resize(distance(eiit,hi_eiit));
              RowGlob[1][row_mat].resize(distance(eiit,hi_eiit));
              
              //                            cout<<"row_mat "<<row_mat<<endl;
              //                            cout<<"RowGlob[row_mat].size() "<<RowGlob[0][row_mat].size()<<endl;
              for(;eiit!=hi_eiit;eiit++) {
                RowGlob[0][row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
                RowGlob[1][row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
              }
              
            }
          }
        }
        
        ColGlob[ff][row_mat].clear(); //without clear it will get the previous value
        if(col_eiit!=col_hi_eiit) {
          FaceEdgeOrder[ee] = col_eiit->get_max_order();
          if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
            assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(col_eiit,col_hi_eiit));
            ColGlob[ff][row_mat].resize(distance(col_eiit,col_hi_eiit));
            //                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
            for(;col_eiit!=col_hi_eiit;col_eiit++) {
              ColGlob[ff][row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
            }
            N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
            diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
            N_edge[ee] = &(N_edge_data[ee][0]);
            diffN_edge[ee] = &(diffN_edge_data[ee][0]);
          }
          count++;  row_mat++;
          //                    cout<<"row_mat  "<<row_mat<<endl;
        }
      }
    }
    
    //        cout<<"\nFor Edges "<<endl;
    //        cout<<"\nrow_mat "<<row_mat << endl;
    //        cout<<"\nRowGlob[ii].size() "<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=1; jj<row_mat; jj++) cout<<RowGlob[ff][jj].size() <<"  ";
    //        cout<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=1; jj<row_mat; jj++) for(int ii=0; ii<RowGlob[ff][jj].size(); ii++) cout<<RowGlob[ff][jj][ii]<<" ";
    //
    //        cout<<"\nFor Edges "<<endl;
    //        cout<<"\nrow_mat "<<row_mat << endl;
    //        cout<<"\ColGlob[ii].size() "<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=1; jj<row_mat; jj++) cout<<ColGlob[ff][jj].size() <<"  ";
    //        cout<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=1; jj<row_mat; jj++) for(int ii=0; ii<ColGlob[ff][jj].size(); ii++) cout<<ColGlob[ff][jj][ii]<<" ";
    //
    //        cout<<endl;
    //        std::string wait;
    //        std::cin >> wait;
    
    
    //Find the shape function N at each gauss point for all the edges and then re-araange in the form as mentioned for nodes
    if(N_edge[0] != NULL){
      ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
      //            cout<<"Hi from insidie"<<endl<<endl;
      int ee = 0; int row_mat1=1;
      for (;ee<3;ee++,row_mat1++){
        rowNMatrices[row_mat1].resize(g_TRI_dim);
        int gg = 0;
        int gg1=0;
        int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
        //                cout<<"nodes_edge  "<<nodes_edge<<endl;
        for(;gg<g_TRI_dim;gg++) {
          rowNMatrices[row_mat1][gg].resize(3,3*nodes_edge);
          rowNMatrices[row_mat1][gg].clear();
          int kk=0;
          for(int ii=0; ii<nodes_edge; ii++){
            for (int jj=0; jj<3; jj++){
              //                            cout<<"jj  "<<jj<<endl;
              //                            cout<<"kk  "<<kk<<endl;
              //                            cout<<"gg1  "<<gg1<<endl;
              //                            cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
              rowNMatrices[row_mat1][gg](jj,kk)=N_edge[ee][gg1]; kk++;
            }
            gg1++;
          }
        }
      }
    }
    
    
    
    
    // Indices and shape funcitons for the faces
    int ff_arr[]={3, 4};  //Canonical numbering of two triangles only
    for(int ff=0;ff<2;ff++) {
      
      EntityHandle Tri_Prism;
      rval = moab.side_element(prism_periodic,2,ff_arr[ff],Tri_Prism); CHKERR_PETSC(rval);
      
      dofs_iterator fiit,hi_fiit, col_fiit, col_hi_fiit;
      fiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",Tri_Prism));
      hi_fiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",Tri_Prism));
      
      col_fiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",Tri_Prism));
      col_hi_fiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",Tri_Prism));
      
      if(col_fiit != col_hi_fiit)  row_mat=4;
      //            cout<<"row_mat hi "<<row_mat<<endl;
      if(ff==0){
        RowGlob[0][row_mat].clear();  RowGlob[1][row_mat].clear();//without clear it will get the previous value
        if(fiit!=hi_fiit) {
          int face_order = fiit->get_max_order();
          //                    cout<<"face_order "<<face_order<<endl;
          assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
          if(NBFACE_H1(face_order)>0) {
            RowGlob[0][row_mat].resize(distance(fiit,hi_fiit));  RowGlob[1][row_mat].resize(distance(fiit,hi_fiit));
            //                        cout<<"RowGlob[0][row_mat].size() "<<RowGlob[0][row_mat].size()<<endl;
            //                        cout<<"RowGlob[1][row_mat].size() "<<RowGlob[1][row_mat].size()<<endl;
            for(;fiit!=hi_fiit;fiit++) {
              RowGlob[0][row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
              RowGlob[1][row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
            }
          }
          
          
          //Find out the shape funcition and rearrange in the form as shown for the nodes (This should be inside fiit!=hi_fiit) and we do it only once for only one side of Prism and use the same to assemble
          vector<double> N_face, diffN_face;
          N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
          diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
          int face_nodes[] = { 0,1,2 };
          ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
          
          rowNMatrices[row_mat].resize(g_TRI_dim);
          int gg = 0;  int gg1=0;
          int nodes_face=NBFACE_H1(face_order);
          
          for(;gg<g_TRI_dim;gg++) {
            rowNMatrices[row_mat][gg].resize(3,3*nodes_face);
            rowNMatrices[row_mat][gg].clear();
            int kk=0;
            for(int ii=0; ii<nodes_face; ii++){
              for (int jj=0; jj<3; jj++){
                //                                  cout<<"jj  "<<jj<<endl;
                //                                  cout<<"kk  "<<kk<<endl;
                //                                  cout<<"gg1  "<<gg1<<endl<<endl;
                //                                  cout<<"N_face[gg1]  "<<N_face[gg1]<<endl<<endl;
                rowNMatrices[row_mat][gg](jj,kk)=N_face[gg1]; kk++;
              }
              gg1++;
            }
            //                        cout<<"rowNMatrices[row_mat][0](jj,kk)  "<<rowNMatrices[row_mat][gg]<<endl;
            //                        std::string wait;
            //                        std::cin >> wait;
          }
        }
      }
      
      
      ColGlob[ff][row_mat].clear(); //without clear it will get the previous value
      if(col_fiit!=col_hi_fiit) {
        int face_order = col_fiit->get_max_order();
        //                cout<<"FaceOrder[ff] "<<face_order<<endl;
        assert((unsigned int)3*NBFACE_H1(face_order)==distance(col_fiit,col_hi_fiit));
        if(NBFACE_H1(face_order)>0) {
          ColGlob[ff][row_mat].resize(distance(col_fiit,col_hi_fiit));
          //                    cout<<"ColGlob[row_mat].size() "<<ColGlob[ff][row_mat].size()<<endl;
          for(;col_fiit!=col_hi_fiit;col_fiit++) {
            ColGlob[ff][row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
          }
        }
        row_mat++;
      }
      
    }
    
    ////        cout<<"\nFor Faces "<<endl;
    //        cout<<"\nRowGlob[ii].size() "<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=4; jj<row_mat; jj++) cout<<RowGlob[ff][jj].size() <<"  ";
    //        cout<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=4; jj<row_mat; jj++) for(int ii=0; ii<RowGlob[ff][jj].size(); ii++) cout<<RowGlob[ff][jj][ii]<<" ";
    //        cout<<"\nColGlob[ii].size() "<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=4; jj<row_mat; jj++) cout<<ColGlob[ff][jj].size() <<"  ";
    //        cout<<endl;
    //        for(int ff=0; ff<2; ff++) for(int jj=4; jj<row_mat; jj++) for(int ii=0; ii<ColGlob[ff][jj].size(); ii++) cout<<ColGlob[ff][jj][ii]<<" ";
    //        cout<<endl;
    //        cout<<"row_mat  "<<row_mat<<endl;
    //        std::string wait;
    //        std::cin >> wait;
    
    PetscFunctionReturn(0);
  }
  
  
  PetscErrorCode ElasticFE_RVELagrange_Periodic::Get_H_mat() {
    PetscFunctionBegin;
    H_mat.resize(2);  //one for -ve triangle and one for +ve triangle of the prism
    for(int ff=0; ff<2; ff++){
      H_mat(ff).resize(row_mat);
      for(int rr=0; rr<row_mat; rr++){
        //                cout<<"rr "<<rr<<endl;
        //                cout<<"(rowNMatrices[rr])[0].size2() "<<(rowNMatrices[rr])[0].size2()<<endl;
        int num_col=(rowNMatrices[rr])[0].size2();
        H_mat[ff][rr].resize(num_col,num_col);
        H_mat[ff][rr].clear();
        for(int ii = 0; ii<num_col; ii++) {
          if(ff==0) H_mat[ff][rr](ii,ii) = -1.0;
          else H_mat[ff][rr](ii,ii) = +1.0;
        }
        //                cout<<"H_mat "<<H_mat[ff][rr]<<endl;
      }
      //            cout<<"ff "<<ff<<endl;
    }
    PetscFunctionReturn(0);
  }
  
  
  //Calculate and assemble NT x N matrix
  PetscErrorCode ElasticFE_RVELagrange_Periodic::Stiffness() {
    PetscFunctionBegin;
    NTN.resize(2);
    for(int ff=0; ff<2; ff++){
      NTN(ff).resize(row_mat,row_mat);
      
      //Calculate C Matrix, i.e. (NT x N)
      for(int rr = 0;rr<row_mat;rr++) {
        //                cout<<" rr = "<<rr<<endl;
        for(int gg = 0;gg<g_TRI_dim;gg++) {
          ublas::matrix<double> &row_Mat = (rowNMatrices[rr])[gg];
          //                    cout<<" row_Mat; "<<row_Mat<<endl;
          double w = area*G_W_TRI[gg];
          for(int cc = 0;cc<row_mat;cc++) {
            ublas::matrix<FieldData> &col_Mat = (rowNMatrices[cc])[gg];
            ublas::matrix<FieldData>  NTN1;
            NTN1.resize(row_Mat.size2(),col_Mat.size2());
            
            //                        cout<<" G_W_TRI[gg] "<<G_W_TRI[gg]<<endl;
            //                        cout<<" row_Mat; "<<row_Mat<<endl;
            //                        cout<<" col_Mat; "<<col_Mat<<endl;
            if(gg == 0) {
              //                            cout<<" row_Mat; "<<row_Mat<<endl;
              //                            cout<<" col_Mat; "<<col_Mat<<endl;
              //calculate NTN1=(w* NT * N)
              NTN1 = prod( w*trans(row_Mat), col_Mat);
              //calculate NTN=(H_mat* w* NT * N)
              NTN(ff)(rr,cc).resize(H_mat[ff][rr].size1(),NTN1.size2());
              NTN(ff)(rr,cc) = prod(H_mat[ff][rr], NTN1);
              //                            cout<<" NTN= "<<NTN(ff)(rr,cc)<<endl;
            }else {
              NTN1 = prod( w*trans(row_Mat), col_Mat);
              NTN(ff)(rr,cc)+=prod(H_mat[ff][rr], NTN1);
            }
          }
        }
      }
      //            cout<<endl<<endl<<endl;
    }
    PetscFunctionReturn(0);
  }
  
  
  PetscErrorCode ElasticFE_RVELagrange_Periodic::Lhs() {
    PetscFunctionBegin;
    ierr = Stiffness(); CHKERRQ(ierr);
    
    //Assembly C
    for(int ff=0; ff<2; ff++){
      for(int rr = 0;rr<row_mat;rr++) {
        for(int cc = 0;cc<row_mat;cc++) {
          //                    cout<<"HI from assembly "<<endl;
          ierr = MatSetValues(Aij,RowGlob[ff][rr].size(),&(RowGlob[ff][rr])[0],ColGlob[ff][cc].size(),&(ColGlob[ff][cc])[0],&(NTN[ff](rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
        }
      }
    }
    
    //Assembly trans(C)
    for(int ff=0; ff<2; ff++){
      for(int rr = 0;rr<row_mat;rr++) {
        for(int cc = 0;cc<row_mat;cc++) {
          NTN[ff](rr,cc) = trans(NTN[ff](rr,cc));
          ierr = MatSetValues(Aij,ColGlob[ff][cc].size(),&(ColGlob[ff][cc])[0],RowGlob[ff][rr].size(),&(RowGlob[ff][rr])[0],&(NTN[ff](rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
        }
      }
    }
    PetscFunctionReturn(0);
  }
  
  
  
  //Calculate the right hand side vector, i.e. f=D_mat * applied_strain and assemble it into the global force vector F
  PetscErrorCode ElasticFE_RVELagrange_Periodic::Rhs() {
    PetscFunctionBegin;
    X_mat.resize(3,6);    X_mat.clear();
    nodes_coord.resize(3,3);
    gauss_coord.resize(3,g_TRI_dim);
    
    for(int ff=0; ff<2; ff++){
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
      //            cout<<"area "<<area << endl;
      
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
            
            //                        //calculate (D_mat1= w * NT * X_mat)
            //                        cout<<"\n row_Mat "<<row_Mat<<endl;
            //                        cout<<"\n col_Mat "<<rr<<col_Mat;
            //                        cout<<"\n w "<<w<<endl;
            
            D_mat1=prod(w*trans(row_Mat), col_Mat);
            //calculate (D_mat = H_mat * D_mat1)
            //                    cout<<"\n rr "<<rr<<endl;
            //                    cout<<"\n D_mat[rr] "<<D_mat[rr]<<endl;
            D_mat[rr]=prod(H_mat[ff][rr], D_mat1);
          }else{
            //calculate (D_mat1= w * NT * X_mat)
            D_mat1=prod(w*trans(row_Mat), col_Mat);
            
            //calculate (D_mat = H_mat * D_mat1)
            D_mat[rr]+=prod(H_mat[ff][rr], D_mat1);
          }
        }
        //            cout<<"\n rr "<<rr<<endl;
        //            cout<<"\n RowGlob[ff][0].size() "<<RowGlob[ff][0].size()<<endl;
        //            for(int ii=0; ii<RowGlob[ff][0].size(); ii++) cout<<RowGlob[ff][0][ii]<<" ";
        //            cout<< "\nD_mat[rr] =  "<<D_mat[rr]<<endl<<endl;
        //Assemble D_mat into global force vector F
        f=prod(D_mat[rr], applied_strain);
        //            cout<<"\n f "<<f<<endl;
        ierr = VecSetValues(F,RowGlob[ff][rr].size(),&(RowGlob[ff][rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }
  
  
  
  PetscErrorCode ElasticFE_RVELagrange_Periodic::operator()() {
    PetscFunctionBegin;
    //        cout<<"Hi from class ElasticFE_RVELagrange_Periodic"<<endl;
    ierr = GetN_and_Indices(); CHKERRQ(ierr);
    ierr = Get_H_mat();
    ierr = Lhs(); CHKERRQ(ierr);
    ierr = Rhs(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
}
