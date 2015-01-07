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

#include <MoFEM.hpp>
#include <ThermalElement.hpp>

using namespace boost::numeric;

namespace MoFEM {

PetscErrorCode ThermalElement::addThermalElements(
  const string problem_name,const string field_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  ierr = mField.add_finite_element("THERMAL_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_FE",field_name); CHKERRQ(ierr);
  if(mField.check_field(mesh_nodals_positions)) {
    ierr = mField.modify_finite_element_add_field_data("THERMAL_FE",mesh_nodals_positions); CHKERRQ(ierr);
  }
  ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_FE"); CHKERRQ(ierr);

  //takes skin of block of entities
  //Skinner skin(&mField.get_moab());
  // loop over all blocksets and get data which name is FluidPressure
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_THERMALSET,it)) {

    Mat_Thermal temp_data;
    ierr = it->get_attribute_data_structure(temp_data); CHKERRQ(ierr);
    setOfBlocks[it->get_msId()].cOnductivity = temp_data.data.Conductivity;
    setOfBlocks[it->get_msId()].cApacity = temp_data.data.HeatCapacity;
    rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"THERMAL_FE"); CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);
}

PetscErrorCode ThermalElement::addThermalFluxElement(
  const string problem_name,const string field_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  ierr = mField.add_finite_element("THERMAL_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_FLUX_FE",field_name); CHKERRQ(ierr);
  if(mField.check_field(mesh_nodals_positions)) {
    ierr = mField.modify_finite_element_add_field_data("THERMAL_FLUX_FE",mesh_nodals_positions); CHKERRQ(ierr);
  }
  ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_FLUX_FE"); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|HEATFLUXSET,it)) {
    ierr = it->get_cubit_bc_data_structure(setOfFluxes[it->get_msId()].dAta); CHKERRQ(ierr);
    rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"THERMAL_FLUX_FE"); CHKERRQ(ierr);
  }

  //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
  //not elegant, but good enough
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
    if(it->get_Cubit_name().compare(0,9,"HEAT_FLUX") == 0) {
      vector<double> data;
      ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
      if(data.size()!=1) {
        SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
      }
      strcpy(setOfFluxes[it->get_msId()].dAta.data.name,"HeatFlu");
      setOfFluxes[it->get_msId()].dAta.data.flag1 = 1;
      setOfFluxes[it->get_msId()].dAta.data.value1 = data[0];
      //cerr << setOfFluxes[it->get_msId()].dAta << endl;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"THERMAL_FLUX_FE"); CHKERRQ(ierr);

    }
  }


  PetscFunctionReturn(0);
}



PetscErrorCode ThermalElement::addThermalConvectionElement(
  const string problem_name,const string field_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  ierr = mField.add_finite_element("THERMAL_CONVECTION_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL_CONVECTION_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_CONVECTION_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_CONVECTION_FE",field_name); CHKERRQ(ierr);
  if(mField.check_field(mesh_nodals_positions)) {
    ierr = mField.modify_finite_element_add_field_data("THERMAL_CONVECTION_FE",mesh_nodals_positions); CHKERRQ(ierr);
  }
  ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_CONVECTION_FE"); CHKERRQ(ierr);

  //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
  //not elegant, but good enough
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
    if(it->get_Cubit_name().compare(0,10,"CONVECTION") == 0) {
      
      vector<double> data;
      ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
      if(data.size()!=2) {
        SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
      }
      setOfConvection[it->get_msId()].cOnvection = data[0];
      setOfConvection[it->get_msId()].tEmperature = data[1];
      //cerr << setOfFluxes[it->get_msId()].dAta << endl;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfConvection[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(setOfConvection[it->get_msId()].tRis,"THERMAL_CONVECTION_FE"); CHKERRQ(ierr);

    }
  }

  PetscFunctionReturn(0);
}


PetscErrorCode ThermalElement::addThermalRadiationElement(
  const string problem_name,const string field_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;
  
  PetscErrorCode ierr;
  ErrorCode rval;
  
  ierr = mField.add_finite_element("THERMAL_RADIATION_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("THERMAL_RADIATION_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("THERMAL_RADIATION_FE",field_name); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("THERMAL_RADIATION_FE",field_name); CHKERRQ(ierr);
  if(mField.check_field(mesh_nodals_positions)) {
    ierr = mField.modify_finite_element_add_field_data("THERMAL_RADIATION_FE",mesh_nodals_positions); CHKERRQ(ierr);
  }
  ierr = mField.modify_problem_add_finite_element(problem_name,"THERMAL_RADIATION_FE"); CHKERRQ(ierr);
  
  //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
  //not elegant, but good enough
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
    if(it->get_Cubit_name().compare(0,9,"RADIATION") == 0) {
      vector<double> data;
      ierr = it->get_Cubit_attributes(data); CHKERRQ(ierr);
      if(data.size()!=3) {
        SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
      }
      setOfRadiation[it->get_msId()].sIgma = data[0];
      setOfRadiation[it->get_msId()].eMissivity = data[1];
      //setOfRadiation[it->get_msId()].aBsorption = data[2];
      setOfRadiation[it->get_msId()].aMbienttEmp = data[2];
      //cerr << setOfFluxes[it->get_msId()].dAta << endl;
      rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfRadiation[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TRIs(setOfRadiation[it->get_msId()].tRis,"THERMAL_RADIATION_FE"); CHKERRQ(ierr);
      
    }
  }
  
  
  PetscFunctionReturn(0);
}


PetscErrorCode ThermalElement::setThermalFiniteElementRhsOperators(string field_name,Vec &F) {
  PetscFunctionBegin;
  map<int,BlockData>::iterator sit = setOfBlocks.begin();
  feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
  for(;sit!=setOfBlocks.end();sit++) {
    //add finite element
    feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,F,sit->second,commonData));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ThermalElement::setThermalFiniteElementLhsOperators(string field_name,Mat A) {
  PetscFunctionBegin;
  map<int,BlockData>::iterator sit = setOfBlocks.begin();
  for(;sit!=setOfBlocks.end();sit++) {
    //add finite elemen
    feLhs.get_op_to_do_Lhs().push_back(new OpThermalLhs(field_name,A,sit->second,commonData));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ThermalElement::setThermalFluxFiniteElementRhsOperators(string field_name,Vec &F,const string mesh_nodals_positions) {
  PetscFunctionBegin;
  bool ho_geometry = false;
  if(mField.check_field(mesh_nodals_positions)) {
    ho_geometry = true;
  }
  map<int,FluxData>::iterator sit = setOfFluxes.begin();
  for(;sit!=setOfFluxes.end();sit++) {
    //add finite element
    feFlux.get_op_to_do_Rhs().push_back(new OpHeatFlux(field_name,F,sit->second,ho_geometry));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ThermalElement::setThermalConvectionFiniteElementRhsOperators(string field_name,Vec &F,const string mesh_nodals_positions) {
  PetscFunctionBegin;
  bool ho_geometry = false;
  if(mField.check_field(mesh_nodals_positions)) {
    ho_geometry = true;
  }
  map<int,ConvectionData>::iterator sit = setOfConvection.begin();
  for(;sit!=setOfConvection.end();sit++) {
    //add finite element
    feConvectionRhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData));
    feConvectionRhs.get_op_to_do_Rhs().push_back(new OpConvectionRhs(field_name,F,sit->second,commonData,ho_geometry));
  }
  PetscFunctionReturn(0);
} 

PetscErrorCode ThermalElement::setThermalConvectionFiniteElementLhsOperators(string field_name,Mat A,const string mesh_nodals_positions) {
  PetscFunctionBegin;
  bool ho_geometry = false;
  if(mField.check_field(mesh_nodals_positions)) {
    ho_geometry = true;
  }
  map<int,ConvectionData>::iterator sit = setOfConvection.begin();
  for(;sit!=setOfConvection.end();sit++) {
    //add finite element
    feConvectionLhs.get_op_to_do_Lhs().push_back(new OpConvectionLhs(field_name,A,sit->second,ho_geometry));
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ThermalElement::setTimeSteppingProblem(string field_name,string rate_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;
 
  bool ho_geometry = false;
  if(mField.check_field(mesh_nodals_positions)) {
    ho_geometry = true;
  }

  {
    map<int,BlockData>::iterator sit = setOfBlocks.begin();
    for(;sit!=setOfBlocks.end();sit++) {
      //add finite element
      //those methods are to calulate matrices on Lhs
      //  feLhs.get_op_to_do_Lhs().push_back(new OpGetTetTemperatureAtGaussPts(field_name,commonData));
      feLhs.get_op_to_do_Lhs().push_back(new OpThermalLhs(field_name,sit->second,commonData));
      feLhs.get_op_to_do_Lhs().push_back(new OpHeatCapacityLsh(field_name,sit->second,commonData));
      //those methods are to calulate vectors on Rhs
      feRhs.get_op_to_do_Rhs().push_back(new OpGetTetTemperatureAtGaussPts(field_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpGetTetRateAtGaussPts(rate_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpGetGradAtGaussPts(field_name,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpThermalRhs(field_name,sit->second,commonData));
      feRhs.get_op_to_do_Rhs().push_back(new OpHeatCapacityRhs(field_name,sit->second,commonData));
    }
  }

  //Flux
  {
    map<int,FluxData>::iterator sit = setOfFluxes.begin();
    for(;sit!=setOfFluxes.end();sit++) {
      //add finite element
      feFlux.get_op_to_do_Rhs().push_back(new OpHeatFlux(field_name,sit->second,ho_geometry));
    }
  }
  
  
  // Convection
  {
    map<int,ConvectionData>::iterator sit = setOfConvection.begin();
    for(;sit!=setOfConvection.end();sit++) {
      //add finite element
	feConvectionRhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData));
      feConvectionRhs.get_op_to_do_Rhs().push_back(new OpConvectionRhs(field_name,sit->second,commonData,ho_geometry));
    }
  }
  {
    map<int,ConvectionData>::iterator sit = setOfConvection.begin();
    for(;sit!=setOfConvection.end();sit++) {
	//add finite element
	feConvectionLhs.get_op_to_do_Lhs().push_back(new OpConvectionLhs(field_name,sit->second,ho_geometry));
    }
  }

  //Radiation
  {
    map<int,RadiationData>::iterator sit = setOfRadiation.begin();
    for(;sit!=setOfRadiation.end();sit++) {
	//add finite element
	feRadiationRhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData));
      feRadiationRhs.get_op_to_do_Rhs().push_back(new OpRadiationRhs(field_name,sit->second,commonData,ho_geometry));
    }
  }
  {
    map<int,RadiationData>::iterator sit = setOfRadiation.begin();
    for(;sit!=setOfRadiation.end();sit++) {
      //add finite element
	feRadiationLhs.get_op_to_do_Rhs().push_back(new OpGetTriTemperatureAtGaussPts(field_name,commonData));
      feRadiationLhs.get_op_to_do_Lhs().push_back(new OpRadiationLhs(field_name,sit->second,commonData,ho_geometry));
    }
  }


  PetscFunctionReturn(0);
}


PetscErrorCode ThermalElement::setTimeSteppingProblem(TsCtx &ts_ctx,string field_name,string rate_name,const string mesh_nodals_positions) {
  PetscFunctionBegin;
  
  PetscErrorCode ierr;
  ierr = setTimeSteppingProblem(field_name,rate_name,mesh_nodals_positions); CHKERRQ(ierr);

  //rhs
  TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_FE",&feRhs));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_FLUX_FE",&feFlux));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_CONVECTION_FE",&feConvectionRhs));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("THERMAL_RADIATION_FE",&feRadiationRhs));
  
  //lhs
  TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("THERMAL_FE",&feLhs));    
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("THERMAL_CONVECTION_FE",&feConvectionLhs));
  loops_to_do_Mat.push_back(TsCtx::loop_pair_type("THERMAL_RADIATION_FE",&feRadiationLhs));    
  //monitor
  //TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();

  PetscFunctionReturn(0);
}

}

