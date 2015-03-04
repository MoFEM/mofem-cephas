/**  
 * \brief Postprocessing stresses for nonolinear analyis
 *
 * Implementation of method for postprocessing stresses.
 *
 */

/* Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
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

#ifndef __POSTPROCSTRESSES_HPP
#define __POSTPROCSTRESSES_HPP

#ifndef WITH_ADOL_C
  #error "MoFEM need to be compiled with ADOL-C"
#endif 

struct PostPorcStress: public TetElementForcesAndSourcesCore::UserDataOperator {

  Interface &postProcMesh;
  vector<EntityHandle> &mapGaussPts;

  NonlinearElasticElement::BlockData &dAta;
  PostPocOnRefinedMesh::CommonData &commonData;
  NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<double> &fUn;

  PostPorcStress(
    Interface &post_proc_mesh,
    vector<EntityHandle> &map_gauss_pts,
    const string field_name,
    NonlinearElasticElement::BlockData &data,
    PostPocOnRefinedMesh::CommonData &common_data,
    NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<double> &fun):
    TetElementForcesAndSourcesCore::UserDataOperator(field_name),
    postProcMesh(post_proc_mesh),mapGaussPts(map_gauss_pts),
    dAta(data),commonData(common_data),fUn(fun) {}


  PetscErrorCode doWork(
    int side,
    EntityType type,
    DataForcesAndSurcesCore::EntData &data) {
    PetscFunctionBegin;

    if(type != MBVERTEX) PetscFunctionReturn(0);
    if(data.getIndices().size()==0) PetscFunctionReturn(0);
    if(dAta.tEts.find(getMoFEMFEPtr()->get_ent()) == dAta.tEts.end()) {
	PetscFunctionReturn(0);
    }

    ErrorCode rval;
    PetscErrorCode ierr;
     
    const FENumeredDofMoFEMEntity *dof_ptr;
    ierr = getMoFEMFEPtr()->get_row_dofs_by_petsc_gloabl_dof_idx(data.getIndices()[0],&dof_ptr); CHKERRQ(ierr);

    string tag_name_piola1 = dof_ptr->get_name()+"_PIOLA1_STRESS";

    int tag_length = 9;
    double def_VAL[tag_length];
    bzero(def_VAL,tag_length*sizeof(double));
    Tag th_piola1;
    rval = postProcMesh.tag_get_handle(
	tag_name_piola1.c_str(),tag_length,MB_TYPE_DOUBLE,th_piola1,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_PETSC(rval);

    int nb_gauss_pts = data.getN().size1();
    if(mapGaussPts.size()!=(unsigned int)nb_gauss_pts) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }
    if(commonData.gradMap[row_field_name].size()!=(unsigned int)nb_gauss_pts) {
	SETERRQ(PETSC_COMM_SELF,MOFEM_DATA_INCONSISTENCT,"data inconsistency");
    }

    ublas::matrix<double> H,invH;
    double detH;

    for(int gg = 0;gg<nb_gauss_pts;gg++) {

	fUn.F.resize(3,3);
	noalias(fUn.F) = (commonData.gradMap[row_field_name])[gg];
	if(commonData.gradMap["MESH_NODE_POSITIONS"].size()==(unsigned int)nb_gauss_pts) {
	  H.resize(3,3);
	  invH.resize(3,3);
	  noalias(H) = (commonData.gradMap["MESH_NODE_POSITIONS"])[gg];
	  ierr = fUn.dEterminatnt(H,detH);  CHKERRQ(ierr);
	  ierr = fUn.iNvert(detH,H,invH); CHKERRQ(ierr);
	  noalias(fUn.F) = prod(fUn.F,invH);  
	}

	ierr = fUn.CalualteP_PiolaKirchhoffI(dAta,getMoFEMFEPtr()); CHKERRQ(ierr);
	rval = postProcMesh.tag_set_data(th_piola1,&mapGaussPts[gg],1,&fUn.P(0,0)); CHKERR_PETSC(rval);

    }

    PetscFunctionReturn(0);

  }

};

#endif //__POSTPROCSTRESSES_HPP



