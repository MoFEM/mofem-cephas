/** \file BodyForce.hpp
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

/** \brief Body forces elements
  * \ingroup mofem_body_forces
  */
struct BodyFroceConstantField {

  MoFEM::Interface &mField;

  struct MyVolumeFE: public MoFEM::VolumeElementForcesAndSourcesCore {
    MyVolumeFE(MoFEM::Interface &m_field): MoFEM::VolumeElementForcesAndSourcesCore(m_field) {}
    int getRule(int order) { return order; };
  };

  MyVolumeFE fe;
  MyVolumeFE& getLoopFe() { return fe; }

  BodyFroceConstantField(
    MoFEM::Interface &m_field):
    mField(m_field),fe(m_field) {}

  struct OpBodyForce: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

    Vec F;
    Block_BodyForces &dAta;
    Range blockTets;
    OpBodyForce(const std::string field_name,Vec _F,Block_BodyForces &data,Range block_tets):
      MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,UserDataOperator::OPROW),
      F(_F),dAta(data),blockTets(block_tets) {}

    VectorDouble Nf;

    PetscErrorCode doWork(
      int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      PetscFunctionBegin;

      if(data.getIndices().size()==0) PetscFunctionReturn(0);
      if(blockTets.find(getNumeredEntFiniteElementPtr()->getEnt())==blockTets.end()) PetscFunctionReturn(0);

      PetscErrorCode ierr;

      const FENumeredDofEntity *dof_ptr;
      ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(
        data.getIndices()[0],&dof_ptr
      ); CHKERRQ(ierr);
      int rank = dof_ptr->getNbOfCoeffs();

      int nb_row_dofs = data.getIndices().size()/rank;

      Nf.resize(data.getIndices().size());
      bzero(&*Nf.data().begin(),data.getIndices().size()*sizeof(FieldData));

      for(unsigned int gg = 0;gg<data.getN().size1();gg++) {
        double val = getVolume()*getGaussPts()(3,gg);
        if(getHoGaussPtsDetJac().size()>0) {
          val *= getHoGaussPtsDetJac()[gg]; ///< higher order geometry
        }
        for(int rr = 0;rr<rank;rr++) {

          double acc;
          if(rr == 0) {
            acc = -dAta.data.acceleration_x;
          } else if(rr == 1) {
            acc = -dAta.data.acceleration_y;
          } else if(rr == 2) {
            acc = -dAta.data.acceleration_z;
          } else {
            SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
          }
          acc *= dAta.data.density;
          cblas_daxpy(nb_row_dofs,val*acc,&data.getN()(gg,0),1,&Nf[rr],rank);

        }
      }
      // std::cerr << dAta.data.acceleration_x << std::endl;
      // std::cerr << dAta.data.acceleration_y << std::endl;
      // std::cerr << dAta.data.acceleration_z << std::endl;
      // std::cerr << dAta.data.density << std::endl;
      // std::cerr << Nf << std::endl;

      ierr = VecSetValues(F,data.getIndices().size(),
        &data.getIndices()[0],&Nf[0],ADD_VALUES
      ); CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


  };

  PetscErrorCode addBlock(const std::string field_name,Vec &F,int ms_id) {
    PetscErrorCode ierr;
    ErrorCode rval;
    const CubitMeshSets *cubit_meshset_ptr;
    MeshsetsManager *mmanager_ptr;
    PetscFunctionBegin;
    ierr = mField.query_interface(mmanager_ptr); CHKERRQ(ierr);
    ierr = mmanager_ptr->getCubitMeshsetPtr(ms_id,BLOCKSET,&cubit_meshset_ptr); CHKERRQ(ierr);
    ierr = cubit_meshset_ptr->getAttributeDataStructure(mapData[ms_id]); CHKERRQ(ierr);
    EntityHandle meshset = cubit_meshset_ptr->getMeshset();
    Range tets;
    rval = mField.get_moab().get_entities_by_type(meshset,MBTET,tets,true); CHKERRQ_MOAB(rval);
    fe.getOpPtrVector().push_back(new OpBodyForce(field_name,F,mapData[ms_id],tets));
    PetscFunctionReturn(0);
  }


  private:

  std::map<int,Block_BodyForces> mapData;

};

#endif //__BODY_FORCE_HPP

/***************************************************************************//**
 * \defgroup mofem_body_forces Body forces elements
 * \ingroup user_modules
 ******************************************************************************/
