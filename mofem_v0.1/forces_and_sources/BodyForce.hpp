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


#ifndef __BODY_FORCE_HPP
#define __BODY_FORCE_HPP

#include "ForcesAndSurcesCore.hpp"

namespace MoFEM {

struct BodyFroceConstantField {

  FieldInterface &mField;
  VolumeH1H1ElementForcesAndSurcesCore fe;

  VolumeH1H1ElementForcesAndSurcesCore& getLoopFe() { return fe; }

  BodyFroceConstantField(
    FieldInterface &m_field):
    mField(m_field),fe(m_field) {}

  struct OpBodyForce: public VolumeH1H1ElementForcesAndSurcesCore::UserDataOperator {

    Vec F;
    Block_BodyForces &dAta;
    OpBodyForce(const string field_name,Vec _F,Block_BodyForces &data):
      VolumeH1H1ElementForcesAndSurcesCore::UserDataOperator(field_name),
      F(_F),dAta(data) {}

    ublas::vector<FieldData> Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofMoFEMEntity *dof_ptr;
      ierr = ptrFE->fe_ptr->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);
      int rank = dof_ptr->get_max_rank();

      int nb_row_dofs = data.getIndices().size()/rank;
      
      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {

	double val = ptrFE->vOlume*ptrFE->gaussPts(3,gg);
	//cerr <<  ptrFE->vOlume << " " << ptrFE->gaussPts(3,gg) << " " << data.getN().size1() << endl << endl;

	for(int rr = 0;rr<rank;rr++) {

	  double acc;
	  if(rr == 0) {
	    acc = dAta.data.acceleration_x;
	  } else if(rr == 1) {
	    acc = dAta.data.acceleration_y;
	  } else if(rr == 2) {
	    acc = dAta.data.acceleration_z;
	  } else {
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  acc *= dAta.data.density;
	  cblas_daxpy(nb_row_dofs,val*acc,&data.getN()(gg,0),1,&Nf[rr],rank);

	}

      }
    
      //cerr << Nf << endl;

      ierr = VecSetValues(F,data.getIndices().size(),
	&data.getIndices()[0],&Nf[0],ADD_VALUES); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  PetscErrorCode addBlock(const string field_name,Vec &F,int ms_id) {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    const CubitMeshSets *cubit_meshset_ptr;
    ierr = mField.get_Cubit_msId(ms_id,BlockSet,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->get_attribute_data_structure(mapData[ms_id]); CHKERRQ(ierr);     
    fe.get_op_to_do_Rhs().push_back(new OpBodyForce(field_name,F,mapData[ms_id]));
    PetscFunctionReturn(0);
  } 


  private:

  map<int,Block_BodyForces> mapData;

};

}

#endif //__BODY_FORCE_HPP

