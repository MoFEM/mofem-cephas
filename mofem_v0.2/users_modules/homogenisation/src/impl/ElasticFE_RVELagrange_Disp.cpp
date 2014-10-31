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


#include <MoFEM.hpp>
using namespace MoFEM;

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <ElasticFEMethod.hpp>
#include "ElasticFE_RVELagrange_Disp.hpp"

using namespace ObosleteUsersModules;

namespace MoFEM {
  
  ElasticFE_RVELagrange_Disp::ElasticFE_RVELagrange_Disp(
                                                         FieldInterface& _mField,Mat &_Aij,Vec &_D,Vec& _F,ublas::vector<FieldData> _applied_strain,
                                                         const string& _field_main, const string& _field_lagrange, int _rank_field):
  FEMethod_UpLevelStudent(_mField.get_moab(),1), mField(_mField),
  Aij(_Aij),F(_F), applied_strain(_applied_strain),field_main(_field_main), field_lagrange(_field_lagrange), rank_field(_rank_field){
    pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    
    snes_B = Aij;
    snes_x = _D;
    snes_f = F;
    
    RowGlob.resize(1+3+1);    // 1-node, 3-edges 1-face
    rowNMatrices.resize(1+3+1);
    ColGlob.resize(1+3+1);
    
    //      g_TRI_dim = 13;
    //      g_NTRI.resize(3*g_TRI_dim);
    //      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,g_TRI_dim);
    //      G_W_TRI = G_TRI_W13;
    
    g_TRI_dim = 28;
    g_NTRI.resize(3*g_TRI_dim);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X28,G_TRI_Y28,g_TRI_dim);
    G_W_TRI = G_TRI_W28;
    
    
    row_mat=0;  //row_mat=0 for nodes   [1,2,3] for edges, 4 for face (for triangle)
    /*
     Changing the shape function matrix for nodes to somthing like
     N=[N1 0  0  N2 0  0  N3 0  0 ]
     0  N1 0  0  N2 0  0  N3 0
     0  0  N1 0  0  N2 0  0  N3
     */
    
    rowNMatrices[row_mat].resize(g_TRI_dim);
    int gg = 0;
    int gg1=0;
//    cout<<"rank_field "<<rank_field<<endl;
    for(;gg<g_TRI_dim;gg++) {
      rowNMatrices[row_mat][gg].resize(rank_field,3*rank_field);   rowNMatrices[row_mat][gg].clear();
      int kk=0;
      for(int ii=0; ii<3; ii++){
        for (int jj=0; jj<rank_field; jj++){
          rowNMatrices[row_mat][gg](jj,kk)=g_NTRI[gg1]; kk++;
        }
        gg1++;
      }
    }
    
    //shape functions matrix (g_NTRI_mat) to find the gauss point coordinates
    g_NTRI_mat.resize(3,g_TRI_dim);
    int kk=0;
    for (int ii=0; ii<g_TRI_dim; ii++){
      for(int jj=0; jj<3; jj++)
      {
        g_NTRI_mat(jj,ii)=g_NTRI[kk]; kk++;
      }
    }
  };
  
  
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::preProcess() {
    PetscFunctionBegin;
//    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
//    ierr = PetscTime(&v1); CHKERRQ(ierr);
//    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::postProcess() {
    PetscFunctionBegin;
    
    switch(snes_ctx) {
      case CTX_SNESNONE: {
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
        break;
      case CTX_SNESSETFUNCTION: {
        ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
        break;
      case CTX_SNESSETJACOBIAN: {
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
        break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::GetN_and_Indices() {
    PetscFunctionBegin;
    
    //Find out indices for row and column for nodes on the surface, i.e. triangles
    row_mat = 0;
    RowGlob[row_mat].resize(3*rank_field);
    ColGlob[row_mat].resize(3*rank_field);
    typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
    const EntityHandle* conn_face;
    int num_nodes;
    EntityHandle face_tri;  face_tri=fePtr->get_ent();
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
    
    int nn = 0;
    for(;nn<3;nn++) {
      dofs_iterator niit,hi_niit;   //for rows
      dofs_iterator col_niit,hi_col_niit;  // for columns
      string field_name;
      
      //cout<<"field_main  ="<<field_main<<endl<<endl;
      //cout<<"field_lagrange  ="<<field_lagrange<<endl<<endl;
      
      //-------field_main = "DISPLACEMENT"
      //-------field_lagrange = "Lagrange_mul_disp"
      
      //minimum and maximum row and column indices for each node on the surface
      niit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_lagrange, conn_face[nn]));
      hi_niit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_lagrange, conn_face[nn]));
      col_niit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_main, conn_face[nn]));
      hi_col_niit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_main, conn_face[nn]));
      
      
      // two different loops, i.e. one for row and one for column (may be need it for multiphysics problems)
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
    //        std::string wait;
    //        std::cin >> wait;
    
    // Find row and colum indices for Edges
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
      dofs_iterator eiit,hi_eiit,col_eiit,col_hi_eiit;
      //            cout<<"side_number "<<side_number<<endl;
      //            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
      //            cout<<"edge "<<edge<<endl;
      eiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_lagrange, edge));
      hi_eiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_lagrange,edge));
      
      col_eiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_main,edge));
      col_hi_eiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_main,edge));
      
      if(eiit!=hi_eiit) {
        FaceEdgeOrder[ee] = eiit->get_max_order();
        if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
          assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
          RowGlob[row_mat].resize(distance(eiit,hi_eiit));
          ColGlob[row_mat].resize(distance(col_eiit,col_hi_eiit));
          //                    cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
          //                    cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
          for(;eiit!=hi_eiit;eiit++) {
            RowGlob[row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
          }
          
          for(;col_eiit!=col_hi_eiit;col_eiit++) {
            ColGlob[row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
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
    
    //        cout<<"row_mat  =  "<<row_mat<<endl;
    //        cout<<"\nFor Edges "<<endl;
    //        cout<<"\n RowGlob[row_mat].size() "<<RowGlob[1].size()<<endl;
    //        for(int jj=0; jj<3; jj++) for(int ii=0; ii<RowGlob[1].size(); ii++) cout<<RowGlob[jj][ii]<<" ";
    //        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[1].size()<<endl;
    //        for(int jj=0; jj<3; jj++) for(int ii=0; ii<ColGlob[1].size(); ii++) cout<<ColGlob[jj][ii]<<" ";
    //        cout<<"\n\n\n";
    ////        cout<<"\n ColGlob[row_mat].size() "<<ColGlob[1].size()<<endl;
    ////        for(int ii=0; ii<ColGlob[row_mat].size(); ii++) cout<<ColGlob[1][ii]<<" ";
    //        cout<<"\n\n\n"<<endl;
    
    
    //Find the shape function N at each gauss point for all the edges and then re-araange in the form as mentioned for nodes
    if(N_edge[0] != NULL){
      ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
      //            cout<<"Hi from insidie"<<endl<<endl;
      ee = 0; int row_mat1=1;
      for (;ee<3;ee++,row_mat1++){
        rowNMatrices[row_mat1].resize(g_TRI_dim);
        int gg = 0; unsigned int gg1=0;
        int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
        //                cout<<"nodes_edge  "<<nodes_edge<<endl;
        for(;gg<g_TRI_dim;gg++) {
          rowNMatrices[row_mat1][gg].resize(rank_field,rank_field*nodes_edge);   rowNMatrices[row_mat1][gg].clear(); // resize(rank_field,rank_field*nodes_edge) = resize(3,9) for rank_field=3
          int kk=0;
          for(int ii=0; ii<nodes_edge; ii++){
            for (int jj=0; jj<rank_field; jj++){
              //                              cout<<"jj  "<<jj<<endl;
              //                              cout<<"kk  "<<kk<<endl;
              //                              cout<<"gg1  "<<gg1<<endl;
              //                              cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
              rowNMatrices[row_mat1][gg](jj,kk)=N_edge[ee][gg1]; kk++;
            }
            gg1++;
          }
        }
        //          for(int jj=0; jj<13; jj++){
        //            cout<<rowNMatrices[1][jj]<<endl<<endl;
        //          }
        //          std::string wait;
        //          std::cin >> wait;
      }
    }
    
    
    //Find the rows and column indices for face of the triangle
    dofs_iterator fiit,hi_fiit, col_fiit, col_hi_fiit;
    fiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_lagrange,face_tri));
    hi_fiit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_lagrange,face_tri));
    
    col_fiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(field_main,face_tri));
    col_hi_fiit = colPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(field_main,face_tri));
    
    //cout<<"fiit!=hi_fiit  "<<(fiit!=hi_fiit) << endl;
    if(fiit!=hi_fiit) {
      RowGlob[row_mat].resize(distance(fiit,hi_fiit));
      ColGlob[row_mat].resize(distance(col_fiit,col_hi_fiit));
      //            cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
      //            cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
      int face_order = fiit->get_max_order();
      assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
      if(NBFACE_H1(face_order)>0) {
        for(;fiit!=hi_fiit;fiit++) {
          RowGlob[row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
        }
        
        for(;col_fiit!=col_hi_fiit;col_fiit++) {
          ColGlob[row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
        }
      }
      
      //Find out the shape funcition and rearrange in the form as shown for the nodes
      //            cout<<"NBFACE_H1(face_order) "<<NBFACE_H1(face_order)<<endl;
      vector<double> N_face, diffN_face;
      N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
      diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
      int face_nodes[] = { 0,1,2 };
      ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
      
      rowNMatrices[row_mat].resize(g_TRI_dim);
      int gg = 0;
      int gg1=0;
      int nodes_face=NBFACE_H1(face_order);
      
      
      //            for(int ii=0; ii<N_face.size(); ii++) cout<<"N_face  "<<N_face[ii]<<endl;
      
      for(;gg<g_TRI_dim;gg++) {
        rowNMatrices[row_mat][gg].resize(rank_field,rank_field*nodes_face);  rowNMatrices[row_mat][gg].clear();
        int kk=0;
        for(int ii=0; ii<nodes_face; ii++){
          for (int jj=0; jj<rank_field; jj++){
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
  
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::Get_H_mat() {
    PetscFunctionBegin;
    H_mat.resize(row_mat);
    for(int rr=0; rr<row_mat; rr++){
      //            cout<<"rr "<<rr<<endl;
      //            cout<<"(rowNMatrices[rr])[0].size2() "<<(rowNMatrices[rr])[0].size2()<<endl;
      int num_col=(rowNMatrices[rr])[0].size2();
      H_mat[rr].resize(num_col,num_col);
      H_mat[rr].clear();
      for(int ii = 0; ii<num_col; ii++) {
        H_mat[rr](ii,ii) = 1.0;
      }
      //        cout<<"H_mat "<<H_mat[rr]<<endl;
    }
    PetscFunctionReturn(0);
  }
  
  
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::Stiffness() {
    PetscFunctionBegin;
    //        cout<<" row_mat; "<<row_mat<<endl;
    NTN.resize(row_mat,row_mat);
    
    //Calculate C Matrix, i.e. (NT x N)
    //        cout<<" row_mat; "<<row_mat<<endl;
    //        for(int ii=0; ii<39; ii++) cout<<" g_NTRI; "<<g_NTRI[ii]<<endl;
    //        for(int ii=0; ii<13; ii++) cout<<" rowNMatrices; "<<rowNMatrices[0][ii]<<endl;
    for(int rr = 0;rr<row_mat;rr++) {
      //            cout<<" rr = "<<rr<<endl;
      for(int gg = 0;gg<g_TRI_dim;gg++) {
        ublas::matrix<double> &row_Mat = (rowNMatrices[rr])[gg];
        double w = area*G_W_TRI[gg];
        for(int cc = 0;cc<row_mat;cc++) {
          
          ublas::matrix<FieldData> &col_Mat = (rowNMatrices[cc])[gg];
          ublas::matrix<FieldData>  NTN1;
          NTN1.resize(row_Mat.size2(),col_Mat.size2());
          
          if(gg == 0) {
            //                        cout<<" w; "<<w<<endl;
            //                        cout<<" row_Mat; "<<row_Mat<<endl;
            //                        cout<<" col_Mat; "<<col_Mat<<endl;
            
            NTN(rr,cc).resize(row_Mat.size2(),col_Mat.size2());
            //calculate (w* NT * N)
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                        w,&*row_Mat.data().begin(),row_Mat.size2(),
                        &*col_Mat.data().begin(),col_Mat.size2(),
                        0.,&*NTN1.data().begin(),NTN1.size2());
            //                        cout<<" NTN1; "<<NTN1<<endl;
            //calculate (H_mat * w* NT * N)
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        H_mat[rr].size2(),NTN1.size2(),H_mat[rr].size1(),
                        1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                        &*NTN1.data().begin(),NTN1.size2(),
                        0.,&*NTN(rr,cc).data().begin(),NTN(rr,cc).size2());
            //                        cout<<" NTN; "<<NTN<<endl;
          } else {
            //                        cout<<" NTN1.size1() gg; "<<NTN1.size1()<<endl;
            //                        cout<<" NTN1.size2() gg; "<<NTN1.size2()<<endl;
            
            //calculate (w* NT * N)
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                        row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                        w,&*row_Mat.data().begin(),row_Mat.size2(),
                        &*col_Mat.data().begin(),col_Mat.size2(),
                        0.,&*NTN1.data().begin(),NTN1.size2());
            
            //calculate (H_mat * w* NT * N)
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                        H_mat[rr].size2(),NTN1.size2(),H_mat[rr].size1(),
                        1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                        &*NTN1.data().begin(),NTN1.size2(),
                        1.,&*NTN(rr,cc).data().begin(),NTN(rr,cc).size2());
          }
        }
      }
    }
    //        cout<<" NTN "<<NTN<<endl;
    PetscFunctionReturn(0);
  }
  
  //********************************************************************************
  PetscErrorCode ElasticFE_RVELagrange_Disp::Lhs() {
    PetscFunctionBegin;
    ierr = Stiffness(); CHKERRQ(ierr);
    
    //Assembly C with size (3M x 3N), where M is number of nodes on boundary and N are total nodes
    for(int rr = 0;rr<row_mat;rr++) {
      for(int cc = 0;cc<row_mat;cc++) {
        if(ColGlob[cc].size()==0) continue;
        if(RowGlob[rr].size()!=NTN(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(ColGlob[cc].size()!=NTN(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        ierr = MatSetValues(snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(NTN(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
    }
    
    //Assembly trans(C) with size (3Nx3M), where M is number of nodes on boundary and N are total nodes
    for(int rr = 0;rr<row_mat;rr++) {
      for(int cc = 0;cc<row_mat;cc++) {
        ublas::matrix<FieldData> trans_NTN = trans(NTN(rr,cc));
        ierr = MatSetValues(snes_B,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(trans_NTN.data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }
  
  
  //********************************************************************************
  //Calculate the right hand side vector, i.e. f=D_max * applied_strain and assemble it into the global force vector F
  PetscErrorCode ElasticFE_RVELagrange_Disp::Rhs() {
    PetscFunctionBegin;
    X_mat.resize(rank_field,1.5*rank_field+1.5);    X_mat.clear();
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
    
    //cout<<"area "<<area << endl;
    for(int rr=0; rr<row_mat; rr++){
      for(int gg = 0;gg<g_TRI_dim;gg++) {
        double w = area*G_W_TRI[gg];
        
        switch(rank_field) {
          case 3:
            X_mat(0,0)=2.0*gauss_coord(0,gg);  X_mat(0,3)=gauss_coord(1,gg);  X_mat(0,4)=gauss_coord(2,gg);
            X_mat(1,1)=2.0*gauss_coord(1,gg);  X_mat(1,3)=gauss_coord(0,gg);  X_mat(1,5)=gauss_coord(2,gg);
            X_mat(2,2)=2.0*gauss_coord(2,gg);  X_mat(2,4)=gauss_coord(0,gg);  X_mat(2,5)=gauss_coord(1,gg);
            X_mat=0.5*X_mat;
            break;
          case 1:
            X_mat(0,0)=gauss_coord(0,gg);  X_mat(0,1)=gauss_coord(1,gg);  X_mat(0,2)=gauss_coord(2,gg);
            break;
          default:
            SETERRQ(PETSC_COMM_SELF,1,"not implemented");
        }
        
        
        ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
        ublas::matrix<FieldData> &col_Mat = X_mat;
        
        ublas::matrix<FieldData>  D_mat1;    //Dmat1=NT*X_mat
        D_mat1.resize(row_Mat.size2(),col_Mat.size2());
        
        //Integrate D_mat
        if(gg == 0) {
          D_mat[rr].resize(row_Mat.size2(),col_Mat.size2());
          
          //calculate (D_mat1= w * NT * X_mat)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                      row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                      w,&*row_Mat.data().begin(),row_Mat.size2(),
                      &*col_Mat.data().begin(),col_Mat.size2(),
                      0.,&*D_mat1.data().begin(),D_mat1.size2());
          
          //calculate (D_mat = H_mat * D_mat1)
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                      H_mat[rr].size2(),D_mat1.size2(),H_mat[rr].size1(),
                      1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                      &*D_mat1.data().begin(),D_mat1.size2(),
                      0.,&*D_mat[rr].data().begin(),D_mat[rr].size2());
          
        } else {
          
          //calculate (D_mat1= w * NT * X_mat)
          cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
                      row_Mat.size2(),col_Mat.size2(),row_Mat.size1(),
                      w,&*row_Mat.data().begin(),row_Mat.size2(),
                      &*col_Mat.data().begin(),col_Mat.size2(),
                      0.,&*D_mat1.data().begin(),D_mat1.size2());
          
          //calculate (D_mat = H_mat * D_mat1)
          cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
                      H_mat[rr].size2(),D_mat1.size2(),H_mat[rr].size1(),
                      1.0,&*H_mat[rr].data().begin(),H_mat[rr].size2(),
                      &*D_mat1.data().begin(),D_mat1.size2(),
                      1.,&*D_mat[rr].data().begin(),D_mat[rr].size2());
        }
      }
      //cout<< " D_mat[rr] =  "<<D_mat[rr]<<endl<<endl;
      
      f=prod(D_mat[rr], applied_strain);
//      cout<<"f "<<f<<endl;
      
//      if(rank_field==1){  //RHS=D_mat*applied_strain + initial_macro_concentraion
//        ublas::vector<FieldData> f_unit; f_unit.resize(f.size());
//        for(int ii=0; ii<f.size(); ii++){
//          f_unit(ii)=1.0;
//        }
//        f=0*f+0.2*f_unit;   //adding initial moisture concentraiton
//        cout<<"f "<<f<<endl<<endl;
//      }
      
      
      if (snes_ctx==CTX_SNESSETFUNCTION) {f*=-1;}
      
      //Assemble D_mat into global force vector F
      ierr = VecSetValues(snes_f,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
  }
  
  
  //********************************************************************************
  vector<ublas::vector<FieldData> > f_ext;
  vector<ublas::vector<FieldData> > f_ext_trans;
  
  PetscErrorCode ElasticFE_RVELagrange_Disp::Rhs_fext() {
    PetscFunctionBegin;
    ierr = Stiffness(); CHKERRQ(ierr);
    
    ublas::vector<ublas::vector<FieldData> > Data_elm_main;  //u for each element
    ublas::vector<ublas::vector<FieldData> > Data_elm_Lamda; //lambda (Lagrange multiplieres) for each element

//    cout<<"row_mat = "<< row_mat << endl;
    Data_elm_main.resize(row_mat);
    Data_elm_Lamda.resize(row_mat);
    
    for(int rr=0; rr<row_mat; rr++){
        Data_elm_main[rr].resize(RowGlob[rr].size());
        Data_elm_Lamda[rr].resize(RowGlob[rr].size());

        switch(rr) {
          case 0:  //for nodes
//            cout<<"For nodes"<<endl;
            const EntityHandle* conn;
            int num_nodes;
            rval = mField.get_moab().get_connectivity(fePtr->get_ent(),conn,num_nodes,true); CHKERR_PETSC(rval);
            //                    cout<<"num_nodes  =  "<<num_nodes<<endl;
            for(int nn = 0;nn<num_nodes; nn++) {
              for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_main,conn[nn],iit)) {
                Data_elm_main[rr][rank_field*nn+iit->get_dof_rank()]=iit->get_FieldData();
              }
              for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_lagrange,conn[nn],iit)) {
                Data_elm_Lamda[rr][rank_field*nn+iit->get_dof_rank()]=iit->get_FieldData();
              }
            }
//            cout<<"field "<<field_main<<endl;
//            for(int ii=0; ii<Data_elm_main[rr].size(); ii++) cout<<Data_elm_main[rr][ii]<<" ";
//            cout<<endl;
//            cout<<"field "<<field_lagrange<<endl;
//            for(int ii=0; ii<Data_elm_Lamda[rr].size(); ii++) cout<<Data_elm_Lamda[rr][ii]<<" ";
//            cout<<endl;
            break;
          case 1:  case 2:  case 3: { //For edges
//            cout<<"For Edges"<<endl;
            EntityHandle edge;
            rval = moab.side_element(fePtr->get_ent(),1,rr-1,edge); CHKERR_PETSC(rval);
            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_main,edge,iit)) {
              Data_elm_main[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
            }
            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_lagrange,edge,iit)) {
              Data_elm_Lamda[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
            }
//            cout<<"field "<<field_main<<endl;
//            for(int ii=0; ii<Data_elm_main[rr].size(); ii++) cout<<Data_elm_main[rr][ii]<<" ";
//            cout<<endl;
//            cout<<"field "<<field_lagrange<<endl;
//            for(int ii=0; ii<Data_elm_Lamda[rr].size(); ii++) cout<<Data_elm_Lamda[rr][ii]<<" ";
//            cout<<endl;
            break;
          }
            
          case 4: //for face
//            cout<<"For Face"<<endl;
            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_main,fePtr->get_ent(),iit)) {
              Data_elm_main[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
            }
            for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,field_lagrange,fePtr->get_ent(),iit)) {
              Data_elm_Lamda[rr][iit->get_EntDofIdx()]=iit->get_FieldData();
            }
//            cout<<"field "<<field_main<<endl;
//            for(int ii=0; ii<Data_elm_main[rr].size(); ii++) cout<<Data_elm_main[rr][ii]<<" ";
//            cout<<endl;
//            cout<<"field "<<field_lagrange<<endl;
//            for(int ii=0; ii<Data_elm_Lamda[rr].size(); ii++) cout<<Data_elm_Lamda[rr][ii]<<" ";
//            cout<<endl;
            break;
        }
      }

    f_ext.resize(row_mat);    //cu * u
    f_ext_trans.resize(row_mat);  //trans(c)*lamda
    
    for(int rr = 0;rr<row_mat;rr++) {
      if(RowGlob[rr].size()==0) continue;
      int rr_start=0;
      for(int cc = 0;cc<row_mat;cc++) {
        if(ColGlob[cc].size()==0) continue;
//        cout<<"rr "<<rr<<endl;
//        cout<<"cc "<<cc<<endl;
//        cout<<"NTN(rr,cc) "<<NTN(rr,cc)<<endl;
//        cout<<"Data_elm_main[cc] "<<Data_elm_main[cc]<<endl;
//        cout<<"Data_elm_Lamda[cc] "<<Data_elm_Lamda[cc]<<endl;
//        cout<<"f_ext[rr] "<<f_ext[rr]<<endl;
//        cout<<"f_ext_trans[rr] "<<f_ext_trans[rr]<<endl;
        if(rr_start == 0) {
          f_ext[rr]       =  prod(NTN(rr,cc),Data_elm_main[cc]);
          f_ext_trans[rr] =  prod(Data_elm_Lamda[cc], trans(NTN(rr,cc)));
          rr_start++;
        } else {
          f_ext[rr]       +=  prod(NTN(rr,cc),Data_elm_main[cc]);
          f_ext_trans[rr] +=  prod(Data_elm_Lamda[cc], trans(NTN(rr,cc)));
        }
//        cout<<"f_ext[rr] "<<f_ext[rr]<<endl;
//        cout<<"f_ext_trans[rr] "<<f_ext_trans[rr]<<endl;
      }
    }
    
    for(int rr = 0;rr<row_mat;rr++) {
      if(RowGlob[rr].size()==0) continue;
      if(RowGlob[rr].size()!=f_ext[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//      cout<<"ColGlob[rr] "<<ColGlob[rr].size()<<endl;
//      cout<<"f_ext_trans[rr] "<<f_ext_trans[rr]<<endl;
      if(ColGlob[rr].size()!=f_ext_trans[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      ierr = VecSetValues(snes_f,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_ext[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      ierr = VecSetValues(snes_f,ColGlob[rr].size(),&(ColGlob[rr])[0],&(f_ext_trans[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
    }
//    std::string wait;
//    std::cin >> wait;
    PetscFunctionReturn(0);
  }
//********************************************************************************
  
  PetscErrorCode ElasticFE_RVELagrange_Disp::operator()() {
    PetscFunctionBegin;
    //        cout<<"Hi from class ElasticFE_RVELagrange_Disp"<<endl;
    ierr = GetN_and_Indices(); CHKERRQ(ierr);
    ierr = Get_H_mat();
    
    switch(snes_ctx) {
      case CTX_SNESNONE: {
        ierr = Lhs(); CHKERRQ(ierr);
        ierr = Rhs(); CHKERRQ(ierr);
        ierr = Rhs_fext(); CHKERRQ(ierr);
      }
        break;
      case CTX_SNESSETFUNCTION: {
        ierr = Rhs(); CHKERRQ(ierr);
        ierr = Rhs_fext(); CHKERRQ(ierr);
      }
        break;
      case CTX_SNESSETJACOBIAN: {
        ierr = Lhs(); CHKERRQ(ierr);
      }
        break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }
  
}

