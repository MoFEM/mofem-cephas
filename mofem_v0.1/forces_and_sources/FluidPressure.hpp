/* Copyright (C) 2014, Lukasz Kaczmarczyk (likask AT wp.pl)
 * --------------------------------------------------------------
 *
 * Description: Implementation of fluid pressure element
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


#ifndef __FLUID_PRESSURE_HPP
#define __FLUID_PRESSURE_HPP

#include "ForcesAndSurcesCore.hpp"

namespace MoFEM {

struct FluidPressure {

  FieldInterface &mField;
  struct MyTriangleFE: public TriElementForcesAndSurcesCore {
    MyTriangleFE(FieldInterface &_mField): TriElementForcesAndSurcesCore(_mField) {}
    int getRule(int order) { return ceil(order/2); };
  };
  MyTriangleFE fe;
  MyTriangleFE& getLoopFe() { return fe; }

  FluidPressure(FieldInterface &m_field): mField(m_field),fe(mField) {}

  typedef int MeshSetId;
  struct FluidData {
    double dEnsity; ///< fluid density [kg/m^2] or any consistent unit
    ublas::vector<double> aCCeleration; ///< acceleration [m/s^2]
    ublas::vector<double> zEroPressure; ///< fluid level of reference zero pressure.
    Range tRis; ///< range of surface elemennt to which fluid pressure is applied
    friend ostream& operator<<(ostream& os,const FluidPressure::FluidData &e);
  };
  map<MeshSetId,FluidData> setOfFluids;

  PetscErrorCode ierr;
  ErrorCode rval;

  struct OpCalulatePressure: public TriElementForcesAndSurcesCore::UserDataOperator {
    Vec F;
    FluidData &dAta;
    bool allowNegativePressure; ///< allows for negative pressures
    bool hoGeometry;
    OpCalulatePressure(const string field_name,Vec _F,FluidData &data,
      bool allow_negative_pressure,bool ho_geometry):
      TriElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data),allowNegativePressure(allow_negative_pressure),hoGeometry(ho_geometry) {}
    ublas::vector<FieldData> Nf;
    PetscErrorCode ierr;
    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;
      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      EntityHandle ent = getMoFEMFEPtr()->get_ent();
      if(dAta.tRis.find(ent)==dAta.tRis.end()) PetscFunctionReturn(0);

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();
      int nb_row_dofs = data.getIndices().size()/rank;
      
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	ublas::vector<double> dist;
	dist = ublas::matrix_row<ublas::matrix<double> >(getCoordsAtGaussPts(),gg);
	dist -= dAta.zEroPressure;
	double dot = cblas_ddot(3,&dist[0],1,&dAta.aCCeleration[0],1);
	if(!allowNegativePressure) dot = -fmax(0,-dot);
	double pressure = dot*dAta.dEnsity;

	for(int rr = 0;rr<rank;rr++) {
	  double force;
	  if(hoGeometry) {
	    force = pressure*getNormals_at_GaussPt()(gg,rr);
	  } else {
	    force = pressure*getNormal()[rr];
	  }
	  cblas_daxpy(nb_row_dofs,getGaussPts()(2,gg)*force,&data.getN()(gg,0),1,&Nf[rr],rank);
	}

      
      }

      //cerr << Nf << endl;
      //cerr << data.getIndices() << endl;
      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);


      PetscFunctionReturn(0);
    }
  };

  PetscErrorCode addNeumannFluidPressureBCElements(
    const string problem_name,const string field_name,const string mesh_nodals_positions = "MESH_NODE_POSITIONS") {
    PetscFunctionBegin;

    ierr = mField.add_finite_element("FLUID_PRESSURE_FE",MF_ZERO); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_row("FLUID_PRESSURE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_col("FLUID_PRESSURE_FE",field_name); CHKERRQ(ierr);
    ierr = mField.modify_finite_element_add_field_data("FLUID_PRESSURE_FE",field_name); CHKERRQ(ierr);
    if(mField.check_field(mesh_nodals_positions)) {
      ierr = mField.modify_finite_element_add_field_data("FLUID_PRESSURE_FE",mesh_nodals_positions); CHKERRQ(ierr);
    }
    ierr = mField.modify_problem_add_finite_element(problem_name,"FLUID_PRESSURE_FE"); CHKERRQ(ierr);

    //takes skin of block of entities
    Skinner skin(&mField.get_moab());
    // loop over all blocksets and get data which name is FluidPressure
    for(_IT_CUBITMESHSETS_BY_NAME_FOR_LOOP_(mField,"FLUIDPRESSURE",bit)) {

      //get block attributes
      vector<double> attributes;
      ierr = bit->get_Cubit_attributes(attributes); CHKERRQ(ierr);
      if(attributes.size()<7) {
	SETERRQ(PETSC_COMM_SELF,1,"not enough block attributes to deffine fluid pressure element");
      }
      setOfFluids[bit->get_msId()].dEnsity = attributes[0];
      setOfFluids[bit->get_msId()].aCCeleration.resize(3);
      setOfFluids[bit->get_msId()].aCCeleration[0] = attributes[1];
      setOfFluids[bit->get_msId()].aCCeleration[1] = attributes[2];
      setOfFluids[bit->get_msId()].aCCeleration[2] = attributes[3];
      setOfFluids[bit->get_msId()].zEroPressure.resize(3);
      setOfFluids[bit->get_msId()].zEroPressure[0] = attributes[4];
      setOfFluids[bit->get_msId()].zEroPressure[1] = attributes[5];
      setOfFluids[bit->get_msId()].zEroPressure[2] = attributes[6];
      //get blok tetrahedrals and triangles
      Range tets;
      rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
      Range tris;
      rval = mField.get_moab().get_entities_by_type(bit->meshset,MBTRI,setOfFluids[bit->get_msId()].tRis,true); CHKERR_PETSC(rval);
      //this get triangles only on block surfaces
      Range tets_skin_tris;
      rval = skin.find_skin(tets,false,tets_skin_tris); CHKERR(rval);
      setOfFluids[bit->get_msId()].tRis.merge(tets_skin_tris);
      ostringstream ss;
      ss << setOfFluids[bit->get_msId()] << endl;
      PetscPrintf(PETSC_COMM_WORLD,ss.str().c_str());

      ierr = mField.add_ents_to_finite_element_by_TRIs(setOfFluids[bit->get_msId()].tRis,"FLUID_PRESSURE_FE"); CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);
  }


  PetscErrorCode setNeumannFluidPressureFiniteElementOperators(string field_name,Vec F,
    bool allow_negative_pressure = true,bool ho_geometry = false) {
    PetscFunctionBegin;
    map<MeshSetId,FluidData>::iterator sit = setOfFluids.begin();
    for(;sit!=setOfFluids.end();sit++) {
      //add finite element
      fe.get_op_to_do_Rhs().push_back(new OpCalulatePressure(field_name,F,sit->second,allow_negative_pressure,ho_geometry));
    }
    PetscFunctionReturn(0);
  }

};

ostream& operator<<(ostream& os,const FluidPressure::FluidData &e) {
  os << "dEnsity " << e.dEnsity << endl;
  os << "aCCeleration " << e.aCCeleration << endl;
  os << "zEroPressure " << e.zEroPressure << endl;
  return os;
}

}

#endif //__FLUID_PRESSSURE_HPP

