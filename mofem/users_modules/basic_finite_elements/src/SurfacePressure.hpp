/* \fiele SurfacePressure.hpp
  \brief Implementation of pressure and forces on triangles surface

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

#ifndef __SURFACE_PERSSURE_HPP__
#define __SURFACE_PERSSURE_HPP__

/** \brief Force and pressures applied to surfaces
  * \ingroup mofem_static_boundary_conditions
  */
struct NeummanForcesSurface {

  FieldInterface &mField;

  struct MyTriangleFE: public FaceElementForcesAndSourcesCore {
    MyTriangleFE(FieldInterface &m_field);
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

  boost::ptr_vector<MethodForForceScaling> methodsOp;

  /// Operator for force element
  struct OpNeumannForce: public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCForce &dAta;
    boost::ptr_vector<MethodForForceScaling> &methodsOp;

    OpNeumannForce(
      const string field_name,Vec &_F,bCForce &data,
      boost::ptr_vector<MethodForForceScaling> &methods_op);

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

  };

  /// Operator for pressure element
  struct OpNeumannPreassure:public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCPreassure &dAta;
    boost::ptr_vector<MethodForForceScaling> &methodsOp;
    bool hoGeometry;

    OpNeumannPreassure(
      const string field_name,Vec &_F,
      bCPreassure &data,boost::ptr_vector<MethodForForceScaling> &methods_op,
      bool ho_geometry = false);

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

  };

  /// Operator for flux element
  struct OpNeumannFlux:public FaceElementForcesAndSourcesCore::UserDataOperator {

    Vec &F;
    bCPreassure &dAta;
    boost::ptr_vector<MethodForForceScaling> &methodsOp;
    bool hoGeometry;

    OpNeumannFlux(
      const string field_name,Vec &_F,
      bCPreassure &data,boost::ptr_vector<MethodForForceScaling> &methods_op,
      bool ho_geometry);

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data);

  };

  DEPRECATED typedef OpNeumannFlux OpNeumannPreassureFlux;

  /// Add force element operator  (integration on face)
  PetscErrorCode addForce(const string field_name,Vec &F,int ms_id);

  /// Add pressure element operator (integration on face)
  PetscErrorCode addPreassure(const string field_name,Vec &F,int ms_id,bool ho_geometry = false);

  /// Add flux element operator (integration on face)
  PetscErrorCode addFlux(const string field_name,Vec &F,int ms_id,bool ho_geometry = false);

};

/// Meta functions to add elements from blocksets
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

#endif //__SURFACE_PERSSURE_HPP__


/***************************************************************************//**
 * \defgroup mofem_static_boundary_conditions Pressure and force boundary conditions
 * \ingroup user_modules
 ******************************************************************************/
