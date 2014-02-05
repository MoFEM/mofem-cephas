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

#ifndef __POSTPROCVERTEXMETHOD_HPP__
#define __POSTPROCVERTEXMETHOD_HPP__

#include "FieldInterface.hpp"
#include "FieldCore.hpp"


using namespace MoFEM;

// Write Displacements DOFS on Vertices
struct PostProcVertexMethodTemp: public FieldInterface::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;
    
    Tag th_disp;
    string field_name;
    Vec V;
    string tag_name;
    PostProcVertexMethodTemp(Interface& _moab,
                         string _field_name = "TEMPERATURE",Vec _V = PETSC_NULL,string _tag_name = "__NotSet__"):
    EntMethod(),moab(_moab),field_name(_field_name),V(_V),tag_name(_tag_name) {
    }
    
    vector<double> vals;
    Range nodes;
    
    VecScatter ctx;
    Vec V_glob;
    PetscScalar *V_glob_array;
    
    PetscErrorCode preProcess() {
        PetscFunctionBegin;
        PetscPrintf(PETSC_COMM_WORLD,"Start postprocess\n");
        
        MoFEMField_multiIndex::index<FieldName_mi_tag>::type::iterator field_it = moabfields->get<FieldName_mi_tag>().find(field_name);
        if(field_it==moabfields->get<FieldName_mi_tag>().end()) SETERRQ1(PETSC_COMM_SELF,1,"field < %s > not found (top tip: check spelling)",field_name.c_str());
        
        double def_VAL[field_it->get_max_rank()];
        bzero(def_VAL,field_it->get_max_rank()*sizeof(double));
        // create TAG
        if(tag_name == "__NotSet__") {
            tag_name = field_name+"_VAL";
        }
        rval = moab.tag_get_handle(tag_name.c_str(),field_it->get_max_rank(),MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL);
        if(rval==MB_ALREADY_ALLOCATED) rval = MB_SUCCESS;
        CHKERR(rval);
        
        rval = moab.get_entities_by_type(field_it->get_meshset(),MBVERTEX,nodes,true); CHKERR_PETSC(rval);
        vals.resize(nodes.size()*field_it->get_max_rank());
        if(V!=PETSC_NULL) {
            ierr = VecScatterCreateToAll(V,&ctx,&V_glob); CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx,V,V_glob,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
            ierr = VecGetArray(V_glob,&V_glob_array); CHKERRQ(ierr);
        } else {
            V_glob_array = NULL;
        }
        PetscFunctionReturn(0);
    }
    
    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        
        if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
        if(dof_ptr->get_name() != field_name) PetscFunctionReturn(0);
        
        EntityHandle ent = dof_ptr->get_ent();
        int dof_rank = dof_ptr->get_dof_rank();
        double fval;
        if(V_glob_array == NULL) {
            fval = dof_ptr->get_FieldData();
        } else {
            fval = V_glob_array[dof_ptr->get_petsc_gloabl_dof_idx()];
        }
        Range::iterator nit = find(nodes.begin(),nodes.end(),ent);
        if(nit==nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        unsigned int pos = std::distance(nodes.begin(),nit);
        pos = dof_ptr->get_max_rank()*pos+dof_rank;
        if(pos>vals.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        vals[pos] = fval;
        //cerr << pos << " --> " << fval << " ent " << ent << endl;
        PetscFunctionReturn(0);
    }
    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        ierr = moab.tag_set_data(th_disp,nodes,&vals[0]); CHKERRQ(ierr);
        if(V!=PETSC_NULL) {
            ierr = VecRestoreArray(V_glob,&V_glob_array); CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx); CHKERRQ(ierr);
            ierr = VecDestroy(&V_glob); CHKERRQ(ierr);
        }
        PetscPrintf(PETSC_COMM_WORLD,"End postprocess\n");
        PetscFunctionReturn(0);
    }
    
};

#endif // __POSTPROCVERTEXMETHOD_HPP__
