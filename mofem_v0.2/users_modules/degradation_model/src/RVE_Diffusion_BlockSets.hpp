/* Copyright (C) 2015, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
 * --------------------------------------------------------------
 *
 * Description: Reading homogenised thermal conductivity matrix from a binary file already stored otherwise read from MAT_THERMAL blocksets
 *
 * This can be used, if a composite block is placed between two steel blocks, so the thermal propeties of steel blocks will be read
 * from MAT_THERMAL and the themal conductivity matirx will be read from already homogenised and stored binary files. Thermal_capacity
 * for the composite is not homogenised and can be read from the input cubit RVE lockset MAT_RVE_THERMAL.
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

#ifndef __RVE_DIFFUSION_BLOCKSETS_HPP
#define __RVE_DIFFUSION_BLOCKSETS_HPP


namespace MoFEM {
  
  
  struct RVE_Diffusion_BlockSets {

    map<int,ThermalElement::BlockData> &setOfBlocks; ///< maps block set id with appropiate BlockData
    FieldInterface &mField;

    RVE_Diffusion_BlockSets(FieldInterface &m_field, map<int,ThermalElement::BlockData> &set_of_blocks): mField(m_field), setOfBlocks(set_of_blocks) {}
    
    PetscErrorCode addDiffusionElement(const string problem_name, const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;

      PetscErrorCode ierr;
      ErrorCode rval;
      ierr = mField.add_finite_element("DIFFUSION_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_row("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_col("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",field_name); CHKERRQ(ierr);
      if(mField.check_field(mesh_nodals_positions)) {
        ierr = mField.modify_finite_element_add_field_data("DIFFUSION_FE",mesh_nodals_positions); CHKERRQ(ierr);
      }
      ierr = mField.modify_problem_add_finite_element(problem_name,"DIFFUSION_FE"); CHKERRQ(ierr);
      
      
      //There can be situation that composite is attached to steel, so MAT_THERMAL for steel and MAT_RVE_THERMAL for RVE
      for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)) {
//        cout<<"_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_  "<<endl;
        if(it->get_name().compare(0,12,"MAT_MOISTURE") == 0) {
//          cout<<"MAT_THERMAL  "<<endl;
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
        else if(it->get_name().compare(0,17,"MAT_RVE_DIFFUSION") == 0){
          cout<<"MAT_RVE_THERMAL  "<<endl;
          setOfBlocks[it->get_msId()].cOnductivity_mat.resize(3,3); //(3X3) conductivity matrix
          setOfBlocks[it->get_msId()].cOnductivity_mat.clear();
          cout<< "setOfBlocks[it->get_msId()].cOnductivity_mat Before Reading= "<<setOfBlocks[it->get_msId()].cOnductivity_mat<<endl;
          int fd;
          PetscViewer view_in;
          PetscViewerBinaryOpen(PETSC_COMM_WORLD,"diffusion_Dmat.out",FILE_MODE_READ,&view_in);
          PetscViewerBinaryGetDescriptor(view_in,&fd);
          PetscBinaryRead(fd,&setOfBlocks[it->get_msId()].cOnductivity_mat(0,0),9,PETSC_DOUBLE);
          PetscViewerDestroy(&view_in);
          cout<< "setOfBlocks[it->get_msId()].cOnductivity_mat After Reading= "<<setOfBlocks[it->get_msId()].cOnductivity_mat<<endl;
          
          vector<double> RVE_diffusion_data;
          ierr = it->get_Cubit_attributes(RVE_diffusion_data); CHKERRQ(ierr);
          setOfBlocks[it->get_msId()].cApacity = 1.0;
          cout<< "setOfBlocks[it->get_msId()].cApacity= "<<setOfBlocks[it->get_msId()].cApacity<<endl;

          rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,setOfBlocks[it->get_msId()].tEts,true); CHKERR_PETSC(rval);
          ierr = mField.add_ents_to_finite_element_by_TETs(setOfBlocks[it->get_msId()].tEts,"DIFFUSION_FE"); CHKERRQ(ierr);

        }
      }

      PetscFunctionReturn(0);
    }
    
  };
}

#endif //__RVE_DIFFUSION_BLOCKSETS_HPP




