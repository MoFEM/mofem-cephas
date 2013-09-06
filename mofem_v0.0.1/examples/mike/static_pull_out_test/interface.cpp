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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "ElasticFEMethodForInterface_mike.hpp"

using namespace MoFEM;

/*struct MyElasticFEMethod: public InterfaceElasticFEMethod,InterfaceFEMethod {
 
 //    ErrorCode rval;
 //    PetscErrorCode ierr;
 //    Range dummy;
 
 Range& SideSet3, SideSet4;
 
 MyElasticFEMethod(Interface& _moab, double _YoungModulus): InterfaceElasticFEMethod(_moab), InterfaceFEMethod(_moab,_YoungModulus), SideSet3(dummy), SideSet4(dummy){
 DispData.resize(1+6+2);
 };
 
 MyElasticFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,
 double _lambda,double _mu,double _YoungModulus,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,Range &_SideSet4): 
 InterfaceElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), 
 InterfaceFEMethod(_moab,_Aij,_F,_YoungModulus,_SideSet1,_SideSet2), 
 SideSet3(_SideSet3), SideSet4(_SideSet4) {
 DispData.resize(1+6+2);
 };
 
 PetscErrorCode ApplyDirihletBC() {
 PetscFunctionBegin;
 
 //        Tag th_bcs;
 //        double def_VAL=0;
 //        rval = moab.tag_get_handle("BC",1,MB_TYPE_DOUBLE,th_bcs,MB_TAG_CREAT|MB_TAG_DENSE,&def_VAL); CHKERR(rval);
 
 Range SideSet3Edges, SideSet3Nodes, SideSet4Edges, SideSet4Nodes;
 Range SideSet3_, SideSet4_;
 
 rval = moab.get_adjacencies(SideSet3,1,false,SideSet3Edges,Interface::UNION); CHKERR_THROW(rval);
 rval = moab.get_connectivity(SideSet3,SideSet3Nodes,true); CHKERR_THROW(rval);
 rval = moab.get_adjacencies(SideSet4,1,false,SideSet4Edges,Interface::UNION); CHKERR_THROW(rval);
 rval = moab.get_connectivity(SideSet4,SideSet4Nodes,true); CHKERR_THROW(rval);
 
 SideSet3_.insert(SideSet3.begin(),SideSet3.end());
 SideSet3_.insert(SideSet3Edges.begin(),SideSet3Edges.end());
 SideSet3_.insert(SideSet3Nodes.begin(),SideSet3Nodes.end());
 
 SideSet4_.insert(SideSet4.begin(),SideSet4.end());
 SideSet4_.insert(SideSet4Edges.begin(),SideSet4Edges.end());
 SideSet4_.insert(SideSet4Nodes.begin(),SideSet4Nodes.end());
 
 //Set Dirihlet for SideSet1
 DirihletBC.resize(0);
 Range::iterator siit1 = SideSet1_.begin();
 for(;siit1!=SideSet1_.end();siit1++) {
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
 for(;riit!=hi_riit;riit++) {
 //                if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
 DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
 for(int cc = 0;cc<col_mat;cc++) {
 vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
 }
 for(int rr = 0;rr<row_mat;rr++) {
 vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
 }
 //                } //apply fix Y for all entities to sideset 3
 }
 }
 
 //Set Dirihlet for SideSet3
 DirihletBC.resize(0);
 Range::iterator siit2 = SideSet3_.begin();
 for(;siit2!=SideSet3_.end();siit2++) {
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit2);
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit2);
 for(;riit!=hi_riit;riit++) {
 if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
 DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
 for(int cc = 0;cc<col_mat;cc++) {
 vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
 }
 for(int rr = 0;rr<row_mat;rr++) {
 vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
 }
 } //apply fix Y for all entities to sideset 3
 }
 }
 
 //Set Dirihlet for SideSet4
 DirihletBC.resize(0);
 Range::iterator siit3 = SideSet4_.begin();
 for(;siit3!=SideSet4_.end();siit3++) {
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit3);
 FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit3);
 for(;riit!=hi_riit;riit++) {
 if(riit->get_dof_rank()==2){//0-X, 1-Y, 2-Z
 DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
 for(int cc = 0;cc<col_mat;cc++) {
 vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
 }
 for(int rr = 0;rr<row_mat;rr++) {
 vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
 if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
 }
 } //apply fix Z for all entities to sideset 4
 }
 }
 
 PetscFunctionReturn(0);
 }
 
 PetscErrorCode operator()() {
 PetscFunctionBegin;
 ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
 ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);
 ierr = GetMatrices(); CHKERRQ(ierr);
 //Dirihlet Boundary Condition
 ApplyDirihletBC();
 if(Diagonal!=PETSC_NULL) {
 if(DirihletBC.size()>0) {
 DirihletBCDiagVal.resize(DirihletBC.size());
 fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
 ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
 }
 }
 
 //Assembly Aij and F
 ierr = RhsAndLhs(); CHKERRQ(ierr);
 
 //Neumann Boundary Conditions
 ierr = NeumannBC(); CHKERRQ(ierr);
 
 ierr = OpStudentEnd(); CHKERRQ(ierr);
 PetscFunctionReturn(0);
 }
 
 };*/

struct MyElasticFEMethod: public InterfaceElasticFEMethod {
    
    Range& SideSet3, SideSet4;
    
    MyElasticFEMethod(Interface& _moab): InterfaceElasticFEMethod(_moab), SideSet3(dummy), SideSet4(dummy) {};
    
    MyElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,
                      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,Range &_SideSet4): 
    InterfaceElasticFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), SideSet3(_SideSet3), SideSet4(_SideSet4) {};
        
    PetscErrorCode ApplyDirihletBC() {
        PetscFunctionBegin;
        
        Range SideSet1Edges, SideSet1Nodes, SideSet3Edges, SideSet3Nodes, SideSet4Edges, SideSet4Nodes;
        Range SideSet1_,SideSet3_, SideSet4_;

        rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet1_,SideSet1Nodes,true); CHKERR_THROW(rval);
        rval = moab.get_adjacencies(SideSet3,1,false,SideSet3Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet3,SideSet3Nodes,true); CHKERR_THROW(rval);
        rval = moab.get_adjacencies(SideSet4,1,false,SideSet4Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet4,SideSet4Nodes,true); CHKERR_THROW(rval);
        
        SideSet1_.insert(SideSet1.begin(),SideSet1.end());
        SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
        SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
        
        SideSet3_.insert(SideSet3.begin(),SideSet3.end());
        SideSet3_.insert(SideSet3Edges.begin(),SideSet3Edges.end());
        SideSet3_.insert(SideSet3Nodes.begin(),SideSet3Nodes.end());
        
        SideSet4_.insert(SideSet4.begin(),SideSet4.end());
        SideSet4_.insert(SideSet4Edges.begin(),SideSet4Edges.end());
        SideSet4_.insert(SideSet4Nodes.begin(),SideSet4Nodes.end());
        
        //Set Dirihlet for SideSet1
        DirihletBC.resize(0);
        Range::iterator siit1 = SideSet1_.begin();
        for(;siit1!=SideSet1_.end();siit1++) {
            FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
            FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
            for(;riit!=hi_riit;riit++) {
                //                if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
                DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
                for(int cc = 0;cc<col_mat;cc++) {
                    vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
                    if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
                }
                for(int rr = 0;rr<row_mat;rr++) {
                    vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
                    if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
                }
                //                } //apply fix Y for all entities to sideset 3
            }
        }
        
        /* //Set Dirihlet for SideSet3
         DirihletBC.resize(0);
         Range::iterator siit2 = SideSet3_.begin();
         for(;siit2!=SideSet3_.end();siit2++) {
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit2);
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit2);
         for(;riit!=hi_riit;riit++) {
         if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
         DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
         for(int cc = 0;cc<col_mat;cc++) {
         vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
         }
         for(int rr = 0;rr<row_mat;rr++) {
         vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
         }
         } //apply fix Y for all entities to sideset 3
         }
         }
         
         //Set Dirihlet for SideSet4
         DirihletBC.resize(0);
         Range::iterator siit3 = SideSet4_.begin();
         for(;siit3!=SideSet4_.end();siit3++) {
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit3);
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit3);
         for(;riit!=hi_riit;riit++) {
         if(riit->get_dof_rank()==2){//0-X, 1-Y, 2-Z
         DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
         for(int cc = 0;cc<col_mat;cc++) {
         vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
         }
         for(int rr = 0;rr<row_mat;rr++) {
         vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
         }
         } //apply fix Z for all entities to sideset 4
         }
         }*/
        
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        //Dirihlet Boundary Condition
        ApplyDirihletBC();
        if(Diagonal!=PETSC_NULL) {
            if(DirihletBC.size()>0) {
                DirihletBCDiagVal.resize(DirihletBC.size());
                fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
                ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
            }
        }
        
        //Assembly Aij and F
        ierr = RhsAndLhs(); CHKERRQ(ierr);
        
        //Neumann Boundary Conditions
        ierr = NeumannBC(); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
    
};

struct MyInterfaceFEMethod: public InterfaceFEMethod {
    
    Range& SideSet3, SideSet4;
    
    MyInterfaceFEMethod(Interface& _moab,double _YoungModulus): InterfaceFEMethod(_moab,_YoungModulus), SideSet3(dummy), SideSet4(dummy) {
        DispData.resize(1+6+2);
    };
    
    MyInterfaceFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec& _F,double _YoungModulus,
                        Range &_SideSet1,Range &_SideSet2,Range &_SideSet3,Range &_SideSet4):
    InterfaceFEMethod(_moab,_dirihlet_ptr,_Aij,_F,_YoungModulus,_SideSet1,_SideSet2), SideSet3(_SideSet3), SideSet4(_SideSet4) {
        DispData.resize(1+6+2);
    };
    
    PetscErrorCode ApplyDirihletBC() {
        PetscFunctionBegin;
        
        //        Tag th_bcs;
        //        double def_VAL=0;
        //        rval = moab.tag_get_handle("BC",1,MB_TYPE_DOUBLE,th_bcs,MB_TAG_CREAT|MB_TAG_DENSE,&def_VAL); CHKERR(rval);
        
        Range SideSet1Edges, SideSet1Nodes,SideSet3Edges, SideSet3Nodes, SideSet4Edges, SideSet4Nodes;
        Range SideSet1_ ,SideSet3_, SideSet4_;

        rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
        rval = moab.get_adjacencies(SideSet3,1,false,SideSet3Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet3,SideSet3Nodes,true); CHKERR_THROW(rval);
        rval = moab.get_adjacencies(SideSet4,1,false,SideSet4Edges,Interface::UNION); CHKERR_THROW(rval);
        rval = moab.get_connectivity(SideSet4,SideSet4Nodes,true); CHKERR_THROW(rval);
        
        SideSet1_.insert(SideSet1.begin(),SideSet1.end());
        SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
        SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
        
        SideSet3_.insert(SideSet3.begin(),SideSet3.end());
        SideSet3_.insert(SideSet3Edges.begin(),SideSet3Edges.end());
        SideSet3_.insert(SideSet3Nodes.begin(),SideSet3Nodes.end());
        
        SideSet4_.insert(SideSet4.begin(),SideSet4.end());
        SideSet4_.insert(SideSet4Edges.begin(),SideSet4Edges.end());
        SideSet4_.insert(SideSet4Nodes.begin(),SideSet4Nodes.end());
        
        //Set Dirihlet for SideSet1
        /* DirihletBC.resize(0);
         Range::iterator siit1 = SideSet1_.begin();
         for(;siit1!=SideSet1_.end();siit1++) {
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
         for(;riit!=hi_riit;riit++) {
         //                if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
         DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
         for(int cc = 0;cc<col_mat;cc++) {
         vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
         }
         for(int rr = 0;rr<row_mat;rr++) {
         vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
         }
         //                } //apply fully fixex for all entities to sideset 1
         }
         }*/
        
        /* //Set Dirihlet for SideSet3
         DirihletBC.resize(0);
         Range::iterator siit2 = SideSet3_.begin();
         for(;siit2!=SideSet3_.end();siit2++) {
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit2);
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit2);
         for(;riit!=hi_riit;riit++) {
         if(riit->get_dof_rank()==1){//0-X, 1-Y, 2-Z
         DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
         for(int cc = 0;cc<col_mat;cc++) {
         vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
         }
         for(int rr = 0;rr<row_mat;rr++) {
         vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
         }
         } //apply fix Y for all entities to sideset 3
         }
         }
         
         //Set Dirihlet for SideSet4
         DirihletBC.resize(0);
         Range::iterator siit3 = SideSet4_.begin();
         for(;siit3!=SideSet4_.end();siit3++) {
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit3);
         FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit3);
         for(;riit!=hi_riit;riit++) {
         if(riit->get_dof_rank()==2){//0-X, 1-Y, 2-Z
         DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
         for(int cc = 0;cc<col_mat;cc++) {
         vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
         }
         for(int rr = 0;rr<row_mat;rr++) {
         vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
         if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
         }
         } //apply fix Z for all entities to sideset 4
         }
         }*/
        
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);
        
        ierr = RhsAndLhs(); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
};

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  moabField_Core core(moab);
  moabField& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  //Interface
  EntityHandle meshset_interface;
  ierr = mField.get_msId_meshset(5,SideSet,meshset_interface); CHKERRQ(ierr);
  ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
  // stl::bitset see for more details
  BitRefLevel bit_level_interface;
  bit_level_interface.set(0);
  ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,true,true); CHKERRQ(ierr);
  EntityHandle meshset_level_interface;
  rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);

  //update BC for refined (with interface) mesh
  EntityHandle meshset_SideSet1; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(1,SideSet,meshset_SideSet1); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet1,bit_level_interface,meshset_SideSet1,MBTRI,true); CHKERRQ(ierr);
  EntityHandle meshset_SideSet2; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(2,SideSet,meshset_SideSet2); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet2,bit_level_interface,meshset_SideSet2,MBTRI,true); CHKERRQ(ierr);
  EntityHandle meshset_SideSet3; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(3,SideSet,meshset_SideSet3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet3,bit_level_interface,meshset_SideSet3,MBTRI,true); CHKERRQ(ierr);
    EntityHandle meshset_SideSet4; //Dirihlet BC is there
    ierr = mField.get_msId_meshset(4,SideSet,meshset_SideSet4); CHKERRQ(ierr);
    ierr = mField.refine_get_childern(meshset_SideSet4,bit_level_interface,meshset_SideSet4,MBTRI,true); CHKERRQ(ierr);
    
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(1);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(meshset_level_interface,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  /*BitRefLevel bit_level1;
  bit_level1.set(2);
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_PRISM(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet1,bit_level1,meshset_SideSet1,MBTRI,true,3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet2,bit_level1,meshset_SideSet2,MBTRI,true,3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet3,bit_level1,meshset_SideSet3,MBTRI,true,3); CHKERRQ(ierr);*/

  BitRefLevel problem_bit_level = bit_level0;

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  //FE Interface
  ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  /****/
  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2,SideSet3, SideSet4;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(4,SideSet,2,SideSet4,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 4 : %u\n",SideSet4.size());

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  const double alpha = 0.05;
  ExampleDiriheltBC myDirihletBC(moab,SideSet1);
  MyElasticFEMethod MyFE(moab,&myDirihletBC,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2,SideSet3,SideSet4);
  MyInterfaceFEMethod IntMyFE(moab,&myDirihletBC,Aij,F,YoungModulus*alpha,SideSet1,SideSet2,SideSet3,SideSet4);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);


  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PostProcCohesiveForces fe_post_proc_prisms(moab,YoungModulus*alpha);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
  }


  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);


  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

