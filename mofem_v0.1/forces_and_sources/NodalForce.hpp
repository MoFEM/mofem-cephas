/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Description: FIXME
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

#ifndef __NODAL_FORCES_HPP
#define __NODAL_FORCES_HPP

namespace MoFEM {

/** \brief Force applied to nodes
  * \ingroup mofem_static_boundary_conditions
  */
struct NodalForce {

  FieldInterface &mField;
  NodalForce(FieldInterface &m_field): mField(m_field),fe(m_field) {}

  struct MyFE: public VertexElementForcesAndSurcesCore {
    MyFE(FieldInterface &_m_field): VertexElementForcesAndSurcesCore(_m_field) {}
  };

  MyFE fe;
  MyFE& getLoopFe() { return fe; }

  struct bCForce {
    ForceCubitBcData data;
    Range nOdes;
  };
  map<int,bCForce> mapForce;

  boost::ptr_vector<MethodsForOp> methodsOp;

  struct OpNodalForce: public VertexElementForcesAndSurcesCore::UserDataOperator {

    Vec &F;
    bool useSnesF;
    bCForce &dAta;
    boost::ptr_vector<MethodsForOp> &methodsOp;

    OpNodalForce(const string field_name,Vec &_F,bCForce &data,
      boost::ptr_vector<MethodsForOp> &methods_op,bool use_snes_f = false):
      VertexElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),useSnesF(use_snes_f),dAta(data),methodsOp(methods_op) {}

    ublas::vector<FieldData> Nf;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      EntityHandle ent = getMoFEMFEPtr()->get_ent();
      if(dAta.nOdes.find(ent)==dAta.nOdes.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();
      if(rank != 3) {
	SETERRQ(PETSC_COMM_SELF,1,"wrong field rank");
      }
      if(data.getIndices().size()!=(unsigned int)rank) {
	SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }

      Nf.resize(3);
      for(int rr = 0;rr<rank;rr++) {
	if(rr == 0) {
	  Nf[0] = dAta.data.data.value3;
	} else if(rr == 1) {
	  Nf[1] = dAta.data.data.value4;
	} else if(rr == 2) {
	  Nf[2] = dAta.data.data.value5;
	} else {
	  SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	}
      }

      ierr = MethodsForOp::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
      Vec myF = F;
      if(useSnesF) {
	myF = getFEMethod()->snes_f;
      } 
      ierr = VecSetValues(myF,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  PetscErrorCode addForce(const string field_name,Vec &F,int ms_id,bool use_snes_f = false) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_Cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_cubit_bc_data_structure(mapForce[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBVERTEX,mapForce[ms_id].nOdes,true); CHKERR_PETSC(rval);
    fe.get_op_to_do_Rhs().push_back(new OpNodalForce(field_name,F,mapForce[ms_id],methodsOp,use_snes_f));
    PetscFunctionReturn(0);
  }

};

struct MetaNodalForces {

  //nodal forces
  struct TagForceScale: public MethodsForOp {
    FieldInterface &mField;
    double *sCale;
    Tag thScale;
    TagForceScale(FieldInterface &m_field): mField(m_field) {
      ErrorCode rval;
      double def_scale = 1.;
      const EntityHandle root_meshset = mField.get_moab().get_root_set();
      rval = mField.get_moab().tag_get_handle("_LoadFactor_Scale_",1,MB_TYPE_DOUBLE,thScale,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_scale); 
      if(rval == MB_ALREADY_ALLOCATED) {
	rval = mField.get_moab().tag_get_by_ptr(thScale,&root_meshset,1,(const void**)&sCale); CHKERR_THROW(rval);
      } else {
	CHKERR_THROW(rval);
	rval = mField.get_moab().tag_set_data(thScale,&root_meshset,1,&def_scale); CHKERR_THROW(rval);
	rval = mField.get_moab().tag_get_by_ptr(thScale,&root_meshset,1,(const void**)&sCale); CHKERR_THROW(rval);
      }
    }
    PetscErrorCode scaleNf(const FieldInterface::FEMethod *fe,ublas::vector<FieldData> &Nf) {
      PetscFunctionBegin;
      Nf *= *sCale;
      PetscFunctionReturn(0);
    }
  };

  static PetscErrorCode addNodalForceElement (
    FieldInterface &mField,
    const string problem_name,
    const string field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    ierr = mField.add_finite_element("FORCE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_problem_add_finite_element(problem_name,"FORCE_FE"); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      Range edges;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBEDGE,edges,true); CHKERR_PETSC(rval);
      Range tris_nodes;
      rval = mField.get_moab().get_connectivity(tris,tris_nodes); CHKERR_PETSC(rval);
      Range edges_nodes;
      rval = mField.get_moab().get_connectivity(edges,edges_nodes); CHKERR_PETSC(rval);
      Range nodes;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
      nodes = subtract(nodes,tris_nodes);
      nodes = subtract(nodes,edges_nodes);
      ierr = mField.add_ents_to_finite_element_by_VERTICEs(nodes,"FORCE_FE"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  static PetscErrorCode setNodalForceElementOperators( 
    FieldInterface &mField,
    boost::ptr_map<string,NodalForce> &nodal_forces,
    Vec &F,const string field_name) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    string fe_name;
    fe_name = "FORCE_FE";
    nodal_forces.insert(fe_name,new NodalForce(mField));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
      ierr = nodal_forces.at(fe_name).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
      /*ForceCubitBcData data;
      ierr = it->get_cubit_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    PetscFunctionReturn(0);
  }

};

}

#endif //__NODAL_FORCES_HPP
