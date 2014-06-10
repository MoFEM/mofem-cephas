/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Test for linar elastic dynamics.
 *
 * This is not exactly procedure for linear elatic dynamics, since jacobian is
 * evaluated at every time step and snes procedure is involved. However it is
 * implemented like that, to test methodology for general nonlinear problem.
 *
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


#ifndef __NEUMANM_FORCES_HPP
#define __NEUMANN_FORCES_HPP

#include "ForcesAndSurcesCore.hpp"

namespace MoFEM {

struct NeummanForcesSurface {

  FieldInterface &mField;

  struct MyTriangleFE: public TriElementForcesAndSurcesCore {
    MyTriangleFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return ceil(order/2); };
  };

  MyTriangleFE fe;
  MyTriangleFE& getLoopFe() { return fe; }

  NeummanForcesSurface(
    FieldInterface &m_field):
    mField(m_field),fe(m_field) {}

  
  struct OpNeumannForce: public TriElementForcesAndSurcesCore::UserDataOperator {

    Vec &F;
    force_cubit_bc_data &dAta;

    OpNeumannForce(const string field_name,Vec &_F,force_cubit_bc_data &data):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getIndices().size()/rank;
      
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = getArea()*getGaussPts()(2,gg);
	for(int rr = 0;rr<rank;rr++) {

	  double force;
	  if(rr == 0) {
	    force = dAta.data.value3;
	  } else if(rr == 1) {
	    force = dAta.data.value4;
	  } else if(rr == 2) {
	    force = dAta.data.value5;
	  } else {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  force *= dAta.data.value1;
	  cblas_daxpy(nb_row_dofs,val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

	}

      }
    
      //cerr << Nf << endl;

      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpNeumannPreassure: public TriElementForcesAndSurcesCore::UserDataOperator {

    Vec &F;
    pressure_cubit_bc_data &dAta;
    bool ho_geometry;

    OpNeumannPreassure(const string field_name,Vec &_F,
      pressure_cubit_bc_data &data,bool _ho_geometry = false):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),ho_geometry(_ho_geometry) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getIndices().size()/rank;
      
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      //cerr << getNormal() << endl;
      //cerr << getNormals_at_GaussPt() << endl;

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = getGaussPts()(2,gg);
	for(int rr = 0;rr<rank;rr++) {

	  double force;
	  if(ho_geometry) {
	    force = dAta.data.value1*getNormals_at_GaussPt()(gg,rr);
	  } else {
	    force = dAta.data.value1*getNormal()[rr];
	  }
	  cblas_daxpy(nb_row_dofs,val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

	}

      }
    
      /*cerr << "VecSetValues\n";
      cerr << Nf << endl;
      cerr << data.getIndices() << endl;*/
      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode addForce(const string field_name,Vec &F,int ms_id) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_Cubit_msId(ms_id,NodeSet,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_cubit_bc_data_structure(mapForce[ms_id]); CHKERRQ(ierr);
    fe.get_op_to_do_Rhs().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id]));
    PetscFunctionReturn(0);
  }

   PetscErrorCode addPreassure(const string field_name,Vec &F,int ms_id,bool ho_geometry = false) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_Cubit_msId(ms_id,SideSet,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_cubit_bc_data_structure(mapPreassure[ms_id]); CHKERRQ(ierr);
    fe.get_op_to_do_Rhs().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],ho_geometry));
    PetscFunctionReturn(0);
  }

  private:

  map<int,force_cubit_bc_data> mapForce;
  map<int,pressure_cubit_bc_data> mapPreassure;

};

struct MetaNeummanForces {

  static PetscErrorCode addNeumannBCElements(
    FieldInterface &mField,
    const string problem_name,
    const string field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|ForceSet,it)) {
      ostringstream fe_name;
      fe_name << "FORCE_FE_" << it->get_msId();
      ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(fe_name.str(),field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(fe_name.str(),field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name.str(),field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
	ierr = mField.modify_finite_element_add_field_data(fe_name.str(),mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,fe_name.str()); CHKERRQ(ierr);
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,fe_name.str()); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {
      ostringstream fe_name;
      fe_name << "PRESSURE_FE_" << it->get_msId();
      ierr = mField.add_finite_element(fe_name.str()); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row(fe_name.str(),field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col(fe_name.str(),field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data(fe_name.str(),field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
	ierr = mField.modify_finite_element_add_field_data(fe_name.str(),mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,fe_name.str()); CHKERRQ(ierr);
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,fe_name.str()); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  static PetscErrorCode setFiniteElementOperators( 
    FieldInterface &mField,
    boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
    Vec &F,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NodeSet|ForceSet,it)) {
      ostringstream fe_name;
      fe_name << "FORCE_FE_" << it->get_msId();
      string fe_name_str = fe_name.str();
      neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
      ierr = neumann_forces.at(fe_name_str).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
      /*force_cubit_bc_data data;
      ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|PressureSet,it)) {
      ostringstream fe_name;
      fe_name << "PRESSURE_FE_" << it->get_msId();
      string fe_name_str = fe_name.str();
      neumann_forces.insert(fe_name_str,new NeummanForcesSurface(mField));
      bool ho_geometry = mField.check_field(mesh_nodals_positions);
      neumann_forces.at(fe_name_str).addPreassure(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
      /*pressure_cubit_bc_data data;
      ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    PetscFunctionReturn(0);
  }


};

}

#endif //__NEUMAN_FORCES_HPP

