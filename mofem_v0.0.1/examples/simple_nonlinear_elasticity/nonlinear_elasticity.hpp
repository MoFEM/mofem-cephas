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

#ifndef __NONLINEAR_ELASTICITY_HPP__
#define __NONLINEAR_ELASTICITY_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"
#include "moabFEMethod_DriverComplexForLazy.hpp"

#include "complex_for_lazy.h"

namespace MoFEM {

struct SetPositionsEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;
    Interface& moab;

    EntityHandle node;
    double coords[3];

    SetPositionsEntMethod(Interface& _moab): EntMethod(),moab(_moab),node(no_handle) {}
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Set Positions\n");
      PetscFunctionReturn(0);
    } 
     
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"End Set Positions\n");
      PetscFunctionReturn(0);
    }

};

struct ExampleDiriheltBC: public BaseDirihletBC {
    Range SideSet1_;

    ExampleDiriheltBC(Interface &moab,Range& SideSet1) {
	ErrorCode rval;
	Range SideSet1Edges,SideSet1Nodes;
	rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
	rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
	SideSet1_.insert(SideSet1.begin(),SideSet1.end());
	SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
	SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
    }

    PetscErrorCode ApplyDirihletBC(
      moabField::FEMethod *fe_method_ptr,
      vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = SideSet1_.begin();
      for(;siit1!=SideSet1_.end();siit1++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;riit!=hi_riit;riit++) {
	    if(riit->get_name()!="DISPLACEMENT") continue;
	    // all fixed
	    // if some ranks are selected then we could apply BC in particular direction
	    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
	      vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	    }
	  }
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;ciit!=hi_ciit;ciit++) {
	    if(ciit->get_name()!="DISPLACEMENT") continue;
	    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
	      vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	      if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	    }
	  }
      }
      PetscFunctionReturn(0);
    }
};

struct ElasticFEMethod: public FEMethod_DriverComplexForLazy {

  Range& SideSet1;
  Range& SideSet2;
  Range SideSet1_;

  ElasticFEMethod(Interface& _moab,BaseDirihletBC *_dirihlet_bc_method_ptr,double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,int _verbose = 0): 
      FEMethod_DriverComplexForLazy(_moab,_dirihlet_bc_method_ptr,_lambda,_mu,_verbose), SideSet1(_SideSet1),SideSet2(_SideSet2)  {

    set_PhysicalEquationNumber(neohookean);

    Range SideSet1Edges,SideSet1Nodes;
    rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
    rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
    SideSet1_.insert(SideSet1.begin(),SideSet1.end());
    SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
    SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());

  }


  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = FEMethod_DriverComplexForLazy::operator()(SideSet1_,SideSet2); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

}

#endif //__NONLINEAR_ELASTICITY_HPP__
