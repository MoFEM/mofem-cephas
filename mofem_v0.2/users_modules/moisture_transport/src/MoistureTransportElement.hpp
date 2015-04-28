/** \file MoistureTransportElement.hpp
 * \brief Operators and data structures for moisture analys
 *
 * Implementation of moisture transport element for unsteady and steady case.
 *
 */

/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 *
 * This file is part of MoFEM.
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

#ifndef __MOISTURE_TRANSPORT_ELEMENT_HPP
#define __MOISTURE_TRANSPORT_ELEMENT_HPP

/** \brief struture grouping operators and data used for moisture problems
 * \ingroup mofem_moisture_transport_elem
 *
 * In order to assemble matrices and right hand vectors, the loops over
 * elements, enetities over that elememnts and finally loop over intergration
 * points are executed.
 *
 * Following implementation separte those three cegories of loops and to eeach
 * loop attach operator.
 *
 */

struct MoistureTransportElement: public ThermalElement {
  MoistureTransportElement(FieldInterface &m_field): ThermalElement(m_field) {}
  
  /** \brief add diffusion element on tets
   * \infroup mofem_moisture_transport_elem
   *
   * It get data from block set and define elemenet in moab
   *
   * \param problem name
   * \param field name
   * \param name of mesh nodal positions (if not defined nodal coordinates are used)
   */
  PetscErrorCode addDiffusionElement(const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    //      cout<<"insides the addDiffusionElements = "<<endl;
    ierr = mField.add_finite_element("DIFFUSION_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("DIFFUSION_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("DIFFUSION_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"DIFFUSION_FE"); CHKERRQ(ierr);
    
    // loop over all blocksets
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
      if(it->get_name().compare(0,12,"MAT_MOISTURE") == 0){
        Mat_Moisture diffusion_data;
        ierr = it->get_attribute_data_structure(diffusion_data); CHKERRQ(ierr);
        //          cout<<"diffusion_data.data.Diffusivity = "<<diffusion_data.data.Diffusivity<<endl;
        //          cout<<"it->get_msId() = "<<it->get_msId()<<endl;
        //It is moisture conductivity acting the same as heat conductivity in thermal problem
        
        setOfBlocks[it->get_msId()].cOnductivity_mat.resize(3,3); //(3X3) conductivity matrix
        setOfBlocks[it->get_msId()].cOnductivity_mat.clear();
        setOfBlocks[it->get_msId()].cOnductivity_mat(0,0)=diffusion_data.data.Diffusivity;
        setOfBlocks[it->get_msId()].cOnductivity_mat(1,1)=diffusion_data.data.Diffusivity;
        setOfBlocks[it->get_msId()].cOnductivity_mat(2,2)=diffusion_data.data.Diffusivity;

        //setOfBlocks[it->get_msId()].cOnductivity = diffusion_data.data.Diffusivity;
        setOfBlocks[it->get_msId()].cApacity = 1.0; //moisture capicity is 1

        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        //          cout<<"setOfBlocks[it->get_msIdx()].tEts.size() = "<<setOfBlocks[it->get_msId()].tEts.size()<<endl;
        ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"DIFFUSION_FE"); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }
  
  
  /** \brief add darceys flow element on tets
   * \infroup mofem_moisture_transport_elem
   *
   * It get data from block set and define elemenet in moab
   *
   * \param problem name
   * \param field name
   * \param name of mesh nodal positions (if not defined nodal coordinates are used)
   */
  PetscErrorCode addDarceysFlowElement(const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;
    //      cout<<"insides the addDiffusionElements = "<<endl;
    ierr = mField.add_finite_element("DARCEYS_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("DARCEYS_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("DARCEYS_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("DARCEYS_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("DARCEYS_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"DARCEYS_FE"); CHKERRQ(ierr);
    
    // loop over all blocksets
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BLOCKSET|MAT_MOISTURESET,it)) {
      if(it->get_name().compare(0,19,"MAT_MOISTURE") == 0){
        Mat_Moisture darceys_data;
        ierr = it->get_attribute_data_structure(darceys_data); CHKERRQ(ierr);
        
//        cout<<"darceys_data.data.Viscosity = "<<darceys_data.data.Viscosity<<endl;
//        cout<<"darceys_data.data.Permeability = "<<darceys_data.data.Permeability<<endl;

//        cout<<"it->get_msId() = "<<it->get_msId()<<endl;
        //It is moisture conductivity acting the same as heat conductivity in thermal problem
        setOfBlocks[it->get_msId()].cOnductivity_mat.resize(3,3); //(3X3) conductivity matrix
        setOfBlocks[it->get_msId()].cOnductivity_mat.clear();
        setOfBlocks[it->get_msId()].cOnductivity_mat(0,0)=darceys_data.data.Permeability/darceys_data.data.Viscosity;
        setOfBlocks[it->get_msId()].cOnductivity_mat(1,1)=darceys_data.data.Permeability/darceys_data.data.Viscosity;
        setOfBlocks[it->get_msId()].cOnductivity_mat(2,2)=darceys_data.data.Permeability/darceys_data.data.Viscosity;

        //setOfBlocks[it->get_msId()].cOnductivity = darceys_data.data.Permeability/darceys_data.data.Viscosity;
        setOfBlocks[it->get_msId()].cApacity = 1.0; //moisture capicity is 1 (will see this for the Darceys unsteady flow ????)
        
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
        //          cout<<"setOfBlocks[it->get_msIdx()].tEts.size() = "<<setOfBlocks[it->get_msId()].tEts.size()<<endl;
        ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"DARCEYS_FE"); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  /** \brief add diffusion flux element
   * \infroup mofem_moisture_transport_elem
   *
   * It get data from het flux set and define elemenet in moab. Aletrantively
   * uses block set with name HET_FLUX.
   *
   * \param problem name
   * \param field name
   * \param name of mesh nodal positions (if not defined nodal coordinates are used)
   */
  PetscErrorCode addDiffusionFluxElement(const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    ErrorCode rval;
    
    ierr = mField.add_finite_element("DIFFUSION_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FLUX_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FLUX_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"DIFFUSION_FLUX_FE"); CHKERRQ(ierr);
    
    //this is alternative method for setting boundary conditions, to bypass bu in cubit file reader.
    //not elegant, but good enough
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
      if(it->get_name().compare(0,9,"MASS_FLUX") == 0) {
        vector<double> data;
        ierr = it->get_attributes(data); CHKERRQ(ierr);
        if(data.size()!=1) {
          SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
        }
//          cout<<"data[0]   "<<data[0]<<endl;
//          std::string wait;
//          std::cin >> wait;
        setOfFluxes[it->get_msId()].dAta.data.flag1 = 1;
        setOfFluxes[it->get_msId()].dAta.data.value1 = data[0];
        //cerr << setOfFluxes[it->get_msId()].dAta << endl;
        rval = mField.get_moab().get_entities_by_type(it->meshset,MBTRI,setOfFluxes[it->get_msId()].tRis,true); CHKERR_PETSC(rval);
        ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluxes[it->get_msId()].tRis,"DIFFUSION_FLUX_FE"); CHKERRQ(ierr);
        
      }
    }
    
    PetscFunctionReturn(0);
  }


  /** \brief set up operators for unsedy moisture transport problem
   * \infroup mofem_moisture_transport_elem
   */
  PetscErrorCode setTimeSteppingProblem(TsCtx &ts_ctx,string field_name,string rate_name,const string mesh_nodals_positions= "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    ierr = ThermalElement::setTimeSteppingProblem(field_name,rate_name,mesh_nodals_positions); CHKERRQ(ierr);
    
    //rhs
    TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("DIFFUSION_FE",&feRhs));
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("DIFFUSION_FLUX_FE",&feFlux));
    
    //lhs
    TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("DIFFUSION_FE",&feLhs));
    //monitor
    //TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();
    
    PetscFunctionReturn(0);
  }

};


#endif //__MOISTURE_TRANSPORT_ELEMENT_HPP

/***************************************************************************//**
* \defgroup mofem_moisture_elem Moisture element
* \ingroup mofem_forces_and_sources
******************************************************************************/

