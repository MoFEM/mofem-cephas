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


#ifndef __SURFACE_PERSSURE_HPP
#define __SURFACE_PERSSURE_HPP

namespace MoFEM {

struct MethodsForOp {

  virtual PetscErrorCode scaleNf(const FEMethod *fe,ublas::vector<FieldData> &Nf) = 0;

  static PetscErrorCode applyScale(
    const FEMethod *fe,
    boost::ptr_vector<MethodsForOp> &methodsOp,ublas::vector<FieldData> &Nf) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
    boost::ptr_vector<MethodsForOp>::iterator vit = methodsOp.begin();
    for(;vit!=methodsOp.end();vit++) {
      ierr = vit->scaleNf(fe,Nf); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  
  virtual ~MethodsForOp() {}

};

/** \brief Force and pressures applied to surfaces
  * \ingroup mofem_static_boundary_conditions
  */
struct NeummanForcesSurface {

  FieldInterface &mField;

  struct MyTriangleFE: public FaceElementForcesAndSourcesCore {
    MyTriangleFE(FieldInterface &_mField): FaceElementForcesAndSourcesCore(_mField) {}
    int getRule(int order) { return order; };
  };

  MyTriangleFE fe;
  MyTriangleFE& getLoopFe() { return fe; }

  NeummanForcesSurface(
    FieldInterface &m_field):
    mField(m_field),fe(m_field) {}

  struct bCForce {
    ForceCubitBcData data;
    Range tRis;
  };
  map<int,bCForce> mapForce;
  struct bCPreassure {
    PressureCubitBcData data;
    Range tRis;
  };
  map<int,bCPreassure> mapPreassure;

  boost::ptr_vector<MethodsForOp> methodsOp;

  struct OpNeumannForce: public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCForce &dAta;
    boost::ptr_vector<MethodsForOp> &methodsOp;

    OpNeumannForce(const string field_name,Vec &_F,bCForce &data,
      boost::ptr_vector<MethodsForOp> &methods_op):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),methodsOp(methods_op) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      EntityHandle ent = getMoFEMFEPtr()->get_ent();
      if(dAta.tRis.find(ent)==dAta.tRis.end()) PetscFunctionReturn(0);

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
	    force = dAta.data.data.value3;
	  } else if(rr == 1) {
	    force = dAta.data.data.value4;
	  } else if(rr == 2) {
	    force = dAta.data.data.value5;
	  } else {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  force *= dAta.data.data.value1;
	  cblas_daxpy(nb_row_dofs,val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

	}

      }

      ierr = MethodsForOp::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpNeumannPreassure:public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCPreassure &dAta;
    boost::ptr_vector<MethodsForOp> &methodsOp;
    bool ho_geometry;

    OpNeumannPreassure(const string field_name,Vec &_F,
      bCPreassure &data,boost::ptr_vector<MethodsForOp> &methods_op,
      bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),methodsOp(methods_op),ho_geometry(_ho_geometry) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

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
	    force = dAta.data.data.value1*getNormals_at_GaussPt()(gg,rr);
	  } else {
	    force = dAta.data.data.value1*getNormal()[rr];
	  }
	  cblas_daxpy(nb_row_dofs,0.5*val*force,&data.getN()(gg,0),1,&Nf[rr],rank);

	}

      }
    
      /*cerr << "VecSetValues\n";
      cerr << Nf << endl;
      cerr << data.getIndices() << endl;*/
      ierr = MethodsForOp::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };

  struct OpNeumannPreassureFlux:public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCPreassure &dAta;
    boost::ptr_vector<MethodsForOp> &methodsOp;
    bool ho_geometry;

    OpNeumannPreassureFlux(const string field_name,Vec &_F,
      bCPreassure &data,boost::ptr_vector<MethodsForOp> &methods_op,
      bool _ho_geometry = false):
      FaceElementForcesAndSourcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),methodsOp(methods_op),ho_geometry(_ho_geometry) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(dAta.tRis.find(getMoFEMFEPtr()->get_ent())==dAta.tRis.end()) PetscFunctionReturn(0);

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
	double flux;
	if(ho_geometry) {
	  double area = cblas_dnrm2(3,&getNormals_at_GaussPt()(gg,0),1);
	  flux = dAta.data.data.value1*area;
	} else {
	  flux = dAta.data.data.value1*getArea();
	}
	cblas_daxpy(nb_row_dofs,val*flux,&data.getN()(gg,0),1,&*Nf.data().begin(),1);

      }
    
      //cerr << "VecSetValues\n";
      //cerr << Nf << endl;
      //cerr << data.getIndices() << endl;
      ierr = MethodsForOp::applyScale(getFEMethod(),methodsOp,Nf); CHKERRQ(ierr);
      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  };


  PetscErrorCode addForce(const string field_name,Vec &F,int ms_id) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,NODESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_bc_data_structure(mapForce[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapForce[ms_id].tRis,true); CHKERR_PETSC(rval);
    fe.getRowOpPtrVector().push_back(new OpNeumannForce(field_name,F,mapForce[ms_id],methodsOp));
    PetscFunctionReturn(0);
  }

   PetscErrorCode addPreassure(const string field_name,Vec &F,int ms_id,bool ho_geometry = false) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERR_PETSC(rval);
    fe.getRowOpPtrVector().push_back(new OpNeumannPreassure(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
    PetscFunctionReturn(0);
  }

  PetscErrorCode addFlux(const string field_name,Vec &F,int ms_id,bool ho_geometry = false) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_cubit_msId(ms_id,SIDESET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_bc_data_structure(mapPreassure[ms_id].data); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(cubit_meshset_ptr->meshset,MBTRI,mapPreassure[ms_id].tRis,true); CHKERR_PETSC(rval);
    fe.getRowOpPtrVector().push_back(new OpNeumannPreassureFlux(field_name,F,mapPreassure[ms_id],methodsOp,ho_geometry));
    PetscFunctionReturn(0);
  }

  

};

struct MetaNeummanForces {

  static PetscErrorCode addNeumannBCElements(
    FieldInterface &mField,
    const string field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("FORCE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("FORCE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("FORCE_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("FORCE_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"FORCE_FE"); CHKERRQ(ierr);
    }

    ierr = mField.add_finite_element("PRESSURE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("PRESSURE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("PRESSURE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("PRESSURE_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("PRESSURE_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"PRESSURE_FE"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  static PetscErrorCode setNeumannFiniteElementOperators( 
    FieldInterface &mField,
    boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
    Vec &F,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    string fe_name;
    fe_name = "FORCE_FE";
    neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
      ierr = neumann_forces.at(fe_name).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
      /*ForceCubitBcData data;
      ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    fe_name = "PRESSURE_FE";
    neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
      bool ho_geometry = mField.check_field(mesh_nodals_positions);
      ierr =  neumann_forces.at(fe_name).addPreassure(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
      /*PressureCubitBcData data;
      ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    PetscFunctionReturn(0);
  }

  static PetscErrorCode addNeumannFluxBCElements(
    FieldInterface &mField,
    const string field_name,
    const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    ierr = mField.add_finite_element("FLUX_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("FLUX_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("FLUX_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }

    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
      Range tris;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"FLUX_FE"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  static PetscErrorCode setNeumannFluxFiniteElementOperators( 
    FieldInterface &mField,
    boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
    Vec &F,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    string fe_name;
    fe_name = "FLUX_FE";
    neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
      bool ho_geometry = mField.check_field(mesh_nodals_positions);
      ierr = neumann_forces.at(fe_name).addFlux(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
      /*PressureCubitBcData data;
      ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
      my_split << *it << endl;
      my_split << data << endl;*/
    }
    PetscFunctionReturn(0);
  }

};

}

#endif //__SURFACE_PERSSURE_HPP


/***************************************************************************//**
 * \defgroup mofem_static_boundary_conditions Pressure and force boundary conditions
 * \ingroup mofem_forces_and_sources 
 ******************************************************************************/


