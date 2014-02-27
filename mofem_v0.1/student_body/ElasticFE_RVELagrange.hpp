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

#ifndef __ElasticFE_RVELagrange_HPP__
#define __ElasticFE_RVELagrange_HPP__

#include <boost/numeric/ublas/symmetric.hpp>

namespace MoFEM {

struct ElasticFE_RVELagrange: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat Aij;
    Vec F;

    bool propeties_from_BlockSet_Mat_ElasticSet;
    ElasticFE_RVELagrange(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F):
      FEMethod_UpLevelStudent(_mField.get_moab(),_dirihlet_ptr,1), mField(_mField),
      Aij(_Aij),F(_F){
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);


//      if(F!=PETSC_NULL) {
//	//VEC & MAT Options
//	//If index is set to -1 ingonre its assembly
//	VecSetOption(F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
//      }
          
      RowGlob.resize(1+3+1);    // 1-node, 3-edges 1-face
      rowNMatrices.resize(1+3+1);
      ColGlob.resize(1+3+1);
    
      g_TRI_dim = 13;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13);
      G_W_TRI = G_TRI_W13;
      
      row_mat=0;
      rowNMatrices[row_mat].resize(g_TRI_dim);
      unsigned int gg = 0;
      unsigned int gg1=0;
    for(;gg<g_TRI_dim;gg++) {
         rowNMatrices[row_mat][gg].resize(3,9);
         rowNMatrices[row_mat][gg].clear();
          int kk=0;
          for(int ii=0; ii<3; ii++){
              for (int jj=0; jj<3; jj++){
//                  cout<<"jj  "<<jj<<endl;
//                  cout<<"kk  "<<kk<<endl;
//                  cout<<"gg1  "<<gg1<<endl<<endl;
                  rowNMatrices[row_mat][gg](jj,kk)=g_NTRI[gg1]; kk++;
              }
              gg1++;
            }
//        cout<<endl<<endl<<endl;
      }
    };

    ErrorCode rval;
    
    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    vector<double> g_NTRI;
    const double* G_W_TRI;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
      PetscSynchronizedFlush(PETSC_COMM_WORLD); 
      ierr = PetscTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
    
    int row_mat,col_mat,g_TRI_dim;
    vector<vector<DofIdx> > RowGlob;
    vector<vector<DofIdx> > ColGlob;
    vector<vector<ublas::matrix<double> > > rowNMatrices;

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
      PetscFunctionReturn(0);
    }

    virtual PetscErrorCode GetN_and_Indices() {
        PetscFunctionBegin;
        
        row_mat = 0;
        RowGlob[row_mat].resize(9);
        ColGlob[row_mat].resize(9);
        typedef FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dofs_iterator;
        const EntityHandle* conn_face;
        int num_nodes;
        EntityHandle face_tri;  face_tri=fe_ptr->get_ent();
        rval = moab.get_connectivity(face_tri,conn_face,num_nodes,true); CHKERR_PETSC(rval);
//        cout<< "num_nodes ="<<num_nodes << endl;
//        cout<< "conn_face ="<<conn_face << endl;
        int nn = 0;
        for(;nn<3;nn++) {
            dofs_iterator niit,hi_niit;
            dofs_iterator col_niit,hi_col_niit;
            string field_name;
            niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
            hi_niit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",conn_face[nn]));
            col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
            hi_col_niit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",conn_face[nn]));
            
            for(;niit!=hi_niit;niit++,col_niit++) {
                RowGlob[row_mat][nn*niit->get_max_rank()+niit->get_dof_rank()] = niit->get_petsc_gloabl_dof_idx();
                ColGlob[row_mat][nn*col_niit->get_max_rank()+col_niit->get_dof_rank()] = col_niit->get_petsc_gloabl_dof_idx();
            }
        }
//        cout<<"RowGlob[row_mat] "<<RowGlob[row_mat][0]<<endl;
//        cout<<"ColGlob[row_mat] "<<ColGlob[row_mat][0]<<endl;
//        cout<<"\n\n\n\n\n\n\n\n\n\n\n";
        
        
        row_mat=4;
        dofs_iterator fiit,hi_fiit, col_fiit, col_hi_fiit;
        fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",face_tri));
        hi_fiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",face_tri));
        
        col_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",face_tri));
        col_hi_fiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",face_tri));
        
        //cout<<"fiit!=hi_fiit  "<<(fiit!=hi_fiit) << endl;
        if(fiit!=hi_fiit) {
            RowGlob[row_mat].resize(distance(fiit,hi_fiit));
            ColGlob[row_mat].resize(distance(col_fiit,col_hi_fiit));
//            cout<<"RowGlob[row_mat].size() "<<RowGlob[row_mat].size()<<endl;
//            cout<<"ColGlob[row_mat].size() "<<ColGlob[row_mat].size()<<endl;
            int face_order = fiit->get_max_order();
            assert((unsigned int)3*NBFACE_H1(face_order)==distance(fiit,hi_fiit));
            if(NBFACE_H1(face_order)>0) {

                for(;fiit!=hi_fiit;fiit++,col_fiit++) {
                      RowGlob[row_mat][fiit->get_EntDofIdx()]=fiit->get_petsc_gloabl_dof_idx();
                      ColGlob[row_mat][col_fiit->get_EntDofIdx()]=col_fiit->get_petsc_gloabl_dof_idx();
                }
            }
//            cout<<"RowGlob[row_mat] "<<RowGlob[row_mat][0]<<endl;
//            cout<<"ColGlob[row_mat] "<<ColGlob[row_mat][0]<<endl;
            vector<double> N_face, diffN_face;

            N_face.resize(g_TRI_dim*NBFACE_H1(face_order));
            diffN_face.resize(2*g_TRI_dim*NBFACE_H1(face_order));
            int face_nodes[] = { 0,1,2 };
            ierr = H1_FaceShapeFunctions_MBTRI(face_nodes,face_order,&g_NTRI[0],&diffNTRI[0],&N_face[0],&diffN_face[0],g_TRI_dim); CHKERRQ(ierr);
        
            rowNMatrices[row_mat].resize(g_TRI_dim);
            unsigned int gg = 0;
            unsigned int gg1=0;
            int nodes_face=NBFACE_H1(face_order);
            
            for(;gg<g_TRI_dim;gg++) {
                rowNMatrices[row_mat][gg].resize(3,3*nodes_face);
                rowNMatrices[row_mat][gg].clear();
                int kk=0;
                for(int ii=0; ii<nodes_face; ii++){
                    for (int jj=0; jj<3; jj++){
//                          cout<<"jj  "<<jj<<endl;
//                          cout<<"kk  "<<kk<<endl;
//                          cout<<"gg1  "<<gg1<<endl<<endl;
//                            cout<<"N_face[gg1]  "<<N_face[gg1]<<endl<<endl;
                        rowNMatrices[row_mat][gg](jj,kk)=N_face[gg1]; kk++;
                    }
                    gg1++;
                }
//                cout<<endl<<endl<<endl;
            }
        }
        
        
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
        
//        cout << "list dofs\n";
//        for(_IT_GET_FEDATA_BY_TYPE_DOFS_FOR_LOOP_(this,"DISPLACEMENT",MBTRI,dof)) {
//            cout << *dof << endl;
//        }
//        cout << "list<-end\n";
        
        int ee = 0; row_mat=1;
        for(;ee<3;ee++,row_mat++) {
            EntityHandle edge;
            rval = moab.side_element(face_tri,1,ee,edge); CHKERR_PETSC(rval);
            int side_number,offset;
            rval = moab.side_number(face_tri,edge,side_number,FaceEdgeSense[ee],offset); CHKERR_PETSC(rval);
            dofs_iterator eiit,hi_eiit,col_eiit,col_hi_eiit;
//            cout<<"side_number "<<side_number<<endl;
//            cout<<"FaceEdgeSense[ee] "<<FaceEdgeSense[ee]<<endl;
//            cout<<"edge "<<edge<<endl;
            eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            hi_eiit = row_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("Lagrange_mul_disp",edge));
            
            col_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",edge));
            col_hi_eiit = col_multiIndex->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",edge));

            if(eiit!=hi_eiit) {
                FaceEdgeOrder[ee] = eiit->get_max_order();
                if(NBEDGE_H1(FaceEdgeOrder[ee])>0) {
                    assert(3*NBEDGE_H1(FaceEdgeOrder[ee]) == distance(eiit,hi_eiit));
                    RowGlob[row_mat].resize(distance(eiit,hi_eiit));
                    ColGlob[row_mat].resize(distance(col_eiit,col_hi_eiit));
                    for(;eiit!=hi_eiit;eiit++,col_hi_eiit++) {
                        RowGlob[row_mat][eiit->get_EntDofIdx()]=eiit->get_petsc_gloabl_dof_idx();
                        ColGlob[row_mat][col_eiit->get_EntDofIdx()]=col_eiit->get_petsc_gloabl_dof_idx();
                    }
//                    cout<<"FaceEdgeOrder  "<<FaceEdgeOrder[ee]<<endl;
//                    cout<<"NBEDGE_H1(FaceEdgeOrder[ee])  "<<NBEDGE_H1(FaceEdgeOrder[ee])<<endl;
                    N_edge_data[ee].resize(g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    diffN_edge_data[ee].resize(2*g_TRI_dim*NBEDGE_H1(FaceEdgeOrder[ee]));
                    N_edge[ee] = &(N_edge_data[ee][0]);
                    diffN_edge[ee] = &(diffN_edge_data[ee][0]);
                }
            }else {
                FaceEdgeOrder[ee] = 0;
                N_edge[ee] = NULL;
                diffN_edge[ee] = NULL;
            }
        }
        
        if(N_edge[0] != NULL){
            ierr = H1_EdgeShapeFunctions_MBTRI(&FaceEdgeSense[0],&FaceEdgeOrder[0],&g_NTRI[0],&diffNTRI[0],N_edge,diffN_edge,g_TRI_dim); CHKERRQ(ierr);
//            cout<<"Hi from insidie"<<endl<<endl;

            ee = 0; row_mat=1;
            for (;ee<3;ee++,row_mat++){
                rowNMatrices[row_mat].resize(g_TRI_dim);
                unsigned int gg = 0; unsigned int gg1=0;
                int nodes_edge=NBEDGE_H1(FaceEdgeOrder[ee]);
//                cout<<"nodes_edge  "<<nodes_edge<<endl;
                
                
                for(;gg<g_TRI_dim;gg++) {
                    rowNMatrices[row_mat][gg].resize(3,3*nodes_edge);
                    rowNMatrices[row_mat][gg].clear();
                    int kk=0;
                    for(int ii=0; ii<nodes_edge; ii++){
                        for (int jj=0; jj<3; jj++){
//                              cout<<"jj  "<<jj<<endl;
//                              cout<<"kk  "<<kk<<endl;
//                              cout<<"gg1  "<<gg1<<endl;
//                              cout<<"N_face[gg1]  "<<N_edge[ee][gg1]<<endl<<endl;
                              rowNMatrices[row_mat][gg](jj,kk)=N_edge[ee][gg1]; kk++;
                        }
                        gg1++;
                    }
//                    cout<<endl;
                }
            }
        }

        PetscFunctionReturn(0);
    }
    
    
    
    
    
    
    
    
    
    ublas::matrix<ublas::matrix<FieldData> > NTN;
    
    virtual PetscErrorCode Stiffness() {
        PetscFunctionBegin;
        row_mat=4;
        NTN.resize(row_mat,row_mat);
        double coords_face[9];
        rval = moab.get_coords(conn_face,num_nodes_face,coords_face); CHKERR_PETSC(rval);
//        ierr = ShapeDiffMBTRI(diffNTRI); CHKERRQ(ierr);
//        ublas::vector<FieldData,ublas::bounded_array<double,3> > normal(3);
//        ierr = ShapeFaceNormalMBTRI(diffNTRI,coords_face,&*normal.data().begin()); CHKERRQ(ierr);
//        double area = cblas_dnrm2(3,&*normal.data().begin(),1)*0.5;
        
        
        for(int rr = 0;rr<row_mat;rr++) {
            for(int gg = 0;gg<g_TRI_dim;gg++) {
//                ublas::matrix<double> &row_Mat = (rowBMatrices[rr])[gg];
                
//                double w,area_at_Gauss_pt;
//                area_at_Gauss_pt = cblas_dnrm2(3,Normals_at_Gauss_pts[gg].data().begin(),1)*0.5;
//                w = area_at_Gauss_pt*G_W_TRI[gg];
                
//                double w = V*G_W_TET[gg];
//                if(detH.size()>0) {
//                    w *= detH[gg];
//                }
//                BD.resize(6,row_Mat.size2());
//                //ublas::noalias(BD) = prod( w*D,row_Mat );
//                cblas_dsymm(CblasRowMajor,CblasLeft,CblasUpper,
//                            BD.size1(),BD.size2(),
//                            w,&*D.data().begin(),D.size2(),
//                            &*row_Mat.data().begin(),row_Mat.size2(),
//                            0.,&*BD.data().begin(),BD.size2());
//                for(int cc = rr;cc<col_mat;cc++) {
//                    ublas::matrix<FieldData> &col_Mat = (colBMatrices[cc])[gg];
//                    if(gg == 0) {
//                        K(rr,cc).resize(BD.size2(),col_Mat.size2());
//                        //ublas::noalias(K(rr,cc)) = prod(trans(BD) , col_Mat ); // int BT*D*B
//                        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                                    BD.size2(),col_Mat.size2(),BD.size1(),
//                                    1.,&*BD.data().begin(),BD.size2(),
//                                    &*col_Mat.data().begin(),col_Mat.size2(),
//                                    0.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
//                    } else {
//                        //ublas::noalias(K(rr,cc)) += prod(trans(BTD) , col_Mat ); // int BT*D*B
//                        cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
//                                    BD.size2(),col_Mat.size2(),BD.size1(),
//                                    1.,&*BD.data().begin(),BD.size2(),
//                                    &*col_Mat.data().begin(),col_Mat.size2(),
//                                    1.,&*K(rr,cc).data().begin(),K(rr,cc).size2());
//                    }
//                }
            }
        }
        PetscFunctionReturn(0);
    }

    
    
    
    
    
    
    
     virtual PetscErrorCode Lhs() {
        PetscFunctionBegin;
        ierr = Stiffness(); CHKERRQ(ierr);
//        for(int rr = 0;rr<row_mat;rr++) {
//            for(int cc = rr;cc<col_mat;cc++) {
//                if(ColGlob[cc].size()==0) continue;
//                if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//                if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
//                ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
//                if(rr!=cc) {
//                    K(cc,rr) = trans(K(rr,cc));
//                    ierr = MatSetValues(Aij,ColGlob[cc].size(),&(ColGlob[cc])[0],RowGlob[rr].size(),&(RowGlob[rr])[0],&(K(cc,rr).data())[0],ADD_VALUES); CHKERRQ(ierr);
//                }
//            }
//        }
        PetscFunctionReturn(0);
    }

    
    
    
    
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
        cout<<"Hi from class"<<endl;
        
        ierr = GetN_and_Indices(); CHKERRQ(ierr);
        ierr = Lhs(); CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
    }

};

    
}

#endif //__ELASTICFEMETHOD_HPP__
