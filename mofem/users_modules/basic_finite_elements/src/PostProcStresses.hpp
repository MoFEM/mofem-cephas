/** \file PostProcStresses.hpp
 * \brief Post-processing stresses for non-linear analysis
 * \ingroup nonlinear_elastic_elem
 *
 * Implementation of method for post-processing stresses.
 */

/*
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

#ifndef __POSTPROCSTRESSES_HPP__
#define __POSTPROCSTRESSES_HPP__

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif

struct PostPorcStress: public MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator {

  moab::Interface &postProcMesh;
  std::vector<EntityHandle> &mapGaussPts;

  NonlinearElasticElement::BlockData &dAta;
  PostProcVolumeOnRefinedMesh::CommonData &commonData;
  bool fieldDisp;

  PostPorcStress(
    moab::Interface &post_proc_mesh,
    std::vector<EntityHandle> &map_gauss_pts,
    const std::string field_name,
    NonlinearElasticElement::BlockData &data,
    PostProcVolumeOnRefinedMesh::CommonData &common_data,
    bool field_disp = false
  ):
  MoFEM::VolumeElementForcesAndSourcesCore::UserDataOperator(field_name,ForcesAndSurcesCore::UserDataOperator::OPROW),
  postProcMesh(post_proc_mesh),
  mapGaussPts(map_gauss_pts),
  dAta(data),
  commonData(common_data),
  fieldDisp(field_disp) {}

  NonlinearElasticElement::CommonData nonLinearElementCommonData;

  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;

    if(type != MBVERTEX) PetscFunctionReturn(0);
    if(data.getIndices().size()==0) PetscFunctionReturn(0);
    if(dAta.tEts.find(getNumeredEntFiniteElementPtr()->getEnt()) == dAta.tEts.end()) {
      PetscFunctionReturn(0);
    }

    ErrorCode rval;
    PetscErrorCode ierr;

    const FENumeredDofEntity *dof_ptr;
    ierr = getNumeredEntFiniteElementPtr()->getRowDofsByPetscGlobalDofIdx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);

    int id  = dAta.iD;

    Tag th_id;
    int def_block_id = -1;
    rval = postProcMesh.tag_get_handle(
      "BLOCK_ID",1,MB_TYPE_INTEGER,th_id,MB_TAG_CREAT|MB_TAG_SPARSE,&def_block_id); CHKERRQ_MOAB(rval);
    Range::iterator tit = commonData.tEts.begin();
    for(;tit!=commonData.tEts.end();tit++) {
      rval = postProcMesh.tag_set_data(th_id,&*tit,1,&id);  CHKERRQ_MOAB(rval);
    }

    string tag_name_piola1 = dof_ptr->getName()+"_PIOLA1_STRESS";
    string tag_name_energy = dof_ptr->getName()+"_ENERGY_DENSITY";

    int tag_length = 9;
    double def_VAL[tag_length];
    bzero(def_VAL,tag_length*sizeof(double));
    Tag th_piola1,th_energy;
    rval = postProcMesh.tag_get_handle(
      tag_name_piola1.c_str(),tag_length,MB_TYPE_DOUBLE,th_piola1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);
    rval = postProcMesh.tag_get_handle(
      tag_name_energy.c_str(),1,MB_TYPE_DOUBLE,th_energy,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERRQ_MOAB(rval);

    int nb_gauss_pts = data.getN().size1();
    if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
      SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency");
    }
    if(commonData.gradMap[rowFieldName].size()!=(unsigned int)nb_gauss_pts) {
      SETERRQ1(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCY,"data inconsistency, filed <%s> not found",rowFieldName.c_str());
    }

    ublas::matrix<double> H,invH;
    double detH;

    dAta.materialDoublePtr->commonDataPtr = &nonLinearElementCommonData;
    dAta.materialDoublePtr->opPtr = this;
    ierr = dAta.materialDoublePtr->getDataOnPostProcessor(
      commonData.fieldMap,commonData.gradMap
    ); CHKERRQ(ierr);

    nonLinearElementCommonData.dataAtGaussPts = commonData.fieldMap;
    nonLinearElementCommonData.gradAtGaussPts = commonData.gradMap;

    for(int gg = 0;gg<nb_gauss_pts;gg++) {

      dAta.materialDoublePtr->gG = gg;
      dAta.materialDoublePtr->F.resize(3,3);
      noalias(dAta.materialDoublePtr->F) = (commonData.gradMap[rowFieldName])[gg];
      if(fieldDisp) {
        for(int dd = 0;dd<3;dd++) {
          dAta.materialDoublePtr->F(dd,dd) += 1;
        }
      }
      if(commonData.gradMap["MESH_NODE_POSITIONS"].size()==(unsigned int)nb_gauss_pts) {
        H.resize(3,3);
        invH.resize(3,3);
        noalias(H) = (commonData.gradMap["MESH_NODE_POSITIONS"])[gg];
        ierr = dAta.materialDoublePtr->dEterminatnt(H,detH);  CHKERRQ(ierr);
        ierr = dAta.materialDoublePtr->iNvert(detH,H,invH); CHKERRQ(ierr);
        noalias(dAta.materialDoublePtr->F) = prod(dAta.materialDoublePtr->F,invH);
      }

      ierr = dAta.materialDoublePtr->calculateP_PiolaKirchhoffI(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);
      rval = postProcMesh.tag_set_data(th_piola1,&mapGaussPts[gg],1,&dAta.materialDoublePtr->P(0,0)); CHKERRQ_MOAB(rval);
      dAta.materialDoublePtr->calculateElasticEnergy(dAta,getNumeredEntFiniteElementPtr()); CHKERRQ(ierr);
      rval = postProcMesh.tag_set_data(th_energy,&mapGaussPts[gg],1,&dAta.materialDoublePtr->eNergy); CHKERRQ_MOAB(rval);

    }

    PetscFunctionReturn(0);

  }

};

#endif //__POSTPROCSTRESSES_HPP__
