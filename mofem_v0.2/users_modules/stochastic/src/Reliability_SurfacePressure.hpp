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


#ifndef __RELIABILITY_SURFACEPRESSURE_HPP
#define __RELIABILITY_SURFACEPRESSURE_HPP

namespace MoFEM {

  struct ReliabilityForceScale: public MethodsForOp {
    double ScaleFactor;
    ReliabilityForceScale(double _force_val):MethodsForOp(),ScaleFactor(_force_val) {}
    
    /*PetscErrorCode set_ScaleFactor(double force_val) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      ScaleFactor = force_val;cout<<"Scale factor is "<<ScaleFactor<<endl;
      PetscFunctionReturn(0);
    }*/
    virtual PetscErrorCode scaleNf(const FEMethod *fe,ublas::vector<FieldData> &Nf) {
      PetscErrorCode ierr;
      PetscFunctionBegin;
      //ierr = set_ScaleFactor(); CHKERRQ(ierr);
      Nf *= ScaleFactor;
      
      //cout<<"\nHello from scale Nf: "<<ScaleFactor<<endl;
      
      PetscFunctionReturn(0);
    }
    
    virtual ~ReliabilityForceScale() {}
  };
  
  
  struct MyMetaNeummanForces: public MetaNeummanForces {
    
    static PetscErrorCode setNeumannFiniteElementOperators(FieldInterface &mField,
                                                           boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
                                                           Vec &F,
                                                           const string field_name,
                                                           const string mesh_nodals_positions = "MESH_NODE_POSITIONS",
                                                           double Updated_Force_Val = 0.0) {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      string fe_name;
      fe_name = "FORCE_FE";
      double force_val;
      double theScale;
      neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
        ierr = neumann_forces.at(fe_name).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
        ForceCubitBcData data;
        ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
        cout<<"Applied UDF A: "<<data.data.value1<<endl;
        force_val = data.data.value1;
         /*my_split << *it << endl;
         my_split << data << endl;*/
      }
      theScale = Updated_Force_Val/force_val;
      neumann_forces.at(fe_name).methodsOp.push_back(new ReliabilityForceScale(theScale));
      fe_name = "PRESSURE_FE";
      neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
        bool ho_geometry = mField.check_field(mesh_nodals_positions);
        ierr = neumann_forces.at(fe_name).addPreassure(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
        /*PressureCubitBcData data;
         ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
         my_split << *it << endl;
         my_split << data << endl;*/
      }
      neumann_forces.at(fe_name).methodsOp.push_back(new ReliabilityForceScale(1));
      //cout<<"\n\nHello from new meta \n\n";
      PetscFunctionReturn(0);
    }
    
  };
  
  
  struct MyMetaNeummanForces_r_PSFEM: public MetaNeummanForces {
    
    static PetscErrorCode setNeumannFiniteElementOperators(
                                                           FieldInterface &mField,
                                                           boost::ptr_map<string,NeummanForcesSurface> &neumann_forces,
                                                           Vec &F,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
      PetscFunctionBegin;
      PetscErrorCode ierr;
      string fe_name;
      fe_name = "FORCE_FE";
      double force_val;
      neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
        ierr = neumann_forces.at(fe_name).addForce(field_name,F,it->get_msId());  CHKERRQ(ierr);
        ForceCubitBcData data;
        ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
        cout<<"Applied UDF B: "<<data.data.value1<<endl;
        force_val = data.data.value1;
        /*my_split << *it << endl;
         my_split << data << endl;*/
      }
      neumann_forces.at(fe_name).methodsOp.push_back(new ReliabilityForceScale(1/force_val));
      fe_name = "PRESSURE_FE";
      neumann_forces.insert(fe_name,new NeummanForcesSurface(mField));
      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
        bool ho_geometry = mField.check_field(mesh_nodals_positions);
        ierr = neumann_forces.at(fe_name).addPreassure(field_name,F,it->get_msId(),ho_geometry); CHKERRQ(ierr);
        /*PressureCubitBcData data;
         ierr = it->get_bc_data_structure(data); CHKERRQ(ierr);
         my_split << *it << endl;
         my_split << data << endl;*/
      }
      neumann_forces.at(fe_name).methodsOp.push_back(new ReliabilityForceScale(1));
      //cout<<"\n\nHello from new meta for the firsr order derivative of force \n\n";
      PetscFunctionReturn(0);
    }
    
  };

}

#endif //__RELIABILITY_SURFACE_PERSSURE_HPP


/***************************************************************************//**
 * \defgroup mofem_static_boundary_conditions Pressure and force boundary conditions
 * \ingroup mofem_forces_and_sources
 ******************************************************************************/
