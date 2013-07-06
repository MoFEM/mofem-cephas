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

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"

#include "complex_for_lazy.h"

/* 
 * From: Demkowicz
 * To discuss basic a priori error estimates, we return now to basic mathemat-
 * ical issues. We have learned from Cea’s lemma that, for the coercive case,
 * the actual approximation error can always be bounded by the product of a
 * mesh-independent constant and the best approximation error. Thus, to es-
 * timate the approximation error, it is sufficient to estimate the best
 * approxi- mation error. By definition, the best approximation error is always
 * bounded by the norm of the difference between the exact solution and any
 * particular choice of a function that “lives” in the FE space. */


using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_BitFieldId("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_MoFEMFE("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_MoFEMFE_row_add_bit("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_col_add_bit("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_data_add_bit("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_BitProblemId("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_MoFEMFE_add_bit("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_MoFEMFE_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",3); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  /*ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = mField.build_fields(); CHKERRQ(ierr);*/

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  struct SetPositionsEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;

    EntityHandle node;
    double coords[3];

    SetPositionsEntMethod(Interface& _moab): EntMethod(_moab),node(no_handle) {}
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start Set Positions\n");
      PetscFunctionReturn(0);
    } 
     
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"End Set Positions\n");
      PetscFunctionReturn(0);
    }

  };

  SetPositionsEntMethod set_positions(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,set_positions); CHKERRQ(ierr);

  struct ElasticFEMethod: public FEMethod_ComplexForLazy {

    Range& SideSet1;
    Range& SideSet2;
    Range SideSet1_;
    vector<DofIdx> DirihletBC;
    vector<FieldData> DirihletBCDiagVal;
    Vec Diagonal;

    Tag th_res;

    ElasticFEMethod(Interface& _moab,
      double _lambda,double _mu, Range &_SideSet1,Range &_SideSet2, int _verbose = 0): 
      FEMethod_ComplexForLazy(_moab,spatail_analysis,_lambda,_mu,_verbose),
      SideSet1(_SideSet1),SideSet2(_SideSet2) {
      
      Range SideSet1Edges,SideSet1Nodes;
      rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
      rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
      SideSet1_.insert(SideSet1.begin(),SideSet1.end());
      SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
      SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());

      double def_VAL[3] = {0,0,0};
      // create TAG
      rval = moab.tag_get_handle("RESIDUAL_VAL",3,MB_TYPE_DOUBLE,th_res,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);

    };

    PetscErrorCode ApplyDirihletBC() {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = SideSet1_.begin();
      for(;siit1!=SideSet1_.end();siit1++) {
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	for(;riit!=hi_riit;riit++) {
	  if(riit->get_name()!="SPATIAL_POSITION") continue;
	  // all fixed
	  // if some ranks are selected then we could apply BC in particular direction
	  DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	  for(int cc = 0;cc<(1+6+4+1);cc++) {
	    vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	  }
	  for(int rr = 0;rr<(1+6+4+1);rr++) {
	    vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode ApplyDirihletBCFace() {
      PetscFunctionBegin;
      vector<DofIdx>::iterator dit = DirihletBC.begin();
      for(;dit!=DirihletBC.end();dit++) {
	vector<DofIdx>::iterator it = find(FaceNodeIndices.begin(),FaceNodeIndices.end(),*dit);
	if(it!=FaceNodeIndices.end()) *it = -1; // of idx is set -1 row is not assembled
	for(int ee = 0;ee<6;ee++) {
	  it = find(FaceEdgeIndices_data[ee].begin(),FaceEdgeIndices_data[ee].end(),*dit);
	  if(it!=FaceEdgeIndices_data[ee].end()) *it = -1; // of idx is set -1 row is not assembled
	}
      }
      PetscFunctionReturn(0);
    }

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    //set_PhysicalEquationNumber(neohookean);
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	Diagonal = PETSC_NULL;
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatZeroEntries(*snes_A); CHKERRQ(ierr);
	ierr = VecDuplicate(snes_f,&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()() {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(); CHKERRQ(ierr);

    ierr = ApplyDirihletBC(); CHKERRQ(ierr);
    if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
    }

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    double t[] = { 0,0,1e-2, 0,0,1e-2, 0,0,1e-2 };

    switch(ctx) {
      case ctx_SNESSetFunction: { 
        ierr = GetFint(); CHKERRQ(ierr);
	VecSetOption(snes_f, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(snes_f,RowGlob[i_nodes].size(),&(RowGlob[i_nodes][0]),&(Fint_h.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  if(RowGlob[1+ee].size()>0) {
	    ierr = VecSetValues(snes_f,RowGlob[1+ee].size(),&(RowGlob[1+ee][0]),&(Fint_h_edge_data[ee].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	for(int ff = 0;ff<4;ff++) {
	  if(RowGlob[1+6+ff].size()>0) {
	    ierr = VecSetValues(snes_f,RowGlob[1+6+ff].size(),&(RowGlob[1+6+ff][0]),&(Fint_h_face_data[ff].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	for(;siit!=hi_siit;siit++) {
	  Range::iterator fit = find(SideSet2.begin(),SideSet2.end(),siit->ent);
	  if(fit==SideSet2.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetFExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  //cerr << "FExt " << FExt << endl;
	  //cerr << "FaceNodeIndices.size() " << FaceNodeIndices.size() << endl;
	  ierr = ApplyDirihletBCFace(); CHKERRQ(ierr);
	  ierr = VecSetValues(snes_f,FaceNodeIndices.size(),&(FaceNodeIndices[0]),&*FExt.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceEdgeIndices_data[ee].size()>0) {
	      ierr = VecSetValues(snes_f,FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),&*FExt_edge_data[ee].data().begin(),ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(FaceIndices.size()>0) {
	    ierr = VecSetValues(snes_f,FaceIndices.size(),&(FaceIndices[0]),&*FExt_face.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = GetTangent(); CHKERRQ(ierr);
	//cerr << "Khh " << Khh << endl;
	ierr = MatSetValues(*snes_A,
	  RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	  ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	  &*(Khh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	//
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(*snes_A,
	    RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	    ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	    &*(Kedgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(*snes_A,
	    RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	    ColGlob[1+ee].size(),&*(ColGlob[1+ee].begin()),
	    &*(Khedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(*snes_A,
	      RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	      ColGlob[1+eee].size(),&*(ColGlob[1+eee].begin()),
	      &*(Khh_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(*snes_A,
	      RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	      ColGlob[1+6+fff].size(),&*(ColGlob[1+6+fff].begin()),
	      &*(Khh_edgeface_data(ee,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(*snes_A,
	    RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	    ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	    &*(Kfaceh_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(*snes_A,
	    RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	    ColGlob[1+6+ff].size(),&*(ColGlob[1+6+ff].begin()),
	    &*(Khedge_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(*snes_A,
	      RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	      ColGlob[1+eee].size(),&*(ColGlob[1+eee].begin()),
	      &*(Khh_faceedge_data(ff,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(*snes_A,
	      RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	      ColGlob[1+6+fff].size(),&*(ColGlob[1+6+fff].begin()),
	      &*(Khh_faceface_data(ff,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	//
	/*for(;siit!=hi_siit;siit++) {
	  Range::iterator fit = find(SideSet2.begin(),SideSet2.end(),siit->ent);
	  if(fit==SideSet2.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetTangentExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  ierr = ApplyDirihletBCFace(); CHKERRQ(ierr);
	  ierr = MatSetValues(*snes_A,
	    FaceNodeIndices.size(),&(FaceNodeIndices[0]),FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	    &*(KExt_hh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceNodeIndices.size()==0) continue;
	    ierr = MatSetValues(*snes_A,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      &*(KExt_edgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(*snes_A,
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_hedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    for(int eee = 0;eee<3;eee++) {
	      ierr = MatSetValues(*snes_A,
		FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
		FaceEdgeIndices_data[eee].size(),&(FaceEdgeIndices_data[eee][0]),
		&*(KExt_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    }
	    ierr = MatSetValues(*snes_A,
	      FaceIndices.size(),&(FaceIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_faceedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(*snes_A,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceIndices.size(),&(FaceIndices[0]),
	      &*(KExt_edgeface_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(*snes_A,
	    FaceIndices.size(),&(FaceIndices[0]),FaceIndices.size(),&(FaceIndices[0]),
	    &*(KExt_faceface.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}*/
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    /*EntityHandle tet = fe_ent_ptr->get_ent();
    Range faces;
    rval = moab.get_adjacencies(&tet,1,2,false,faces); CHKERR_PETSC(rval);
    for(Range::iterator fit = faces.begin();fit!=faces.end();fit++) {
      ierr = GetFaceIndicesAndData(*fit); CHKERRQ(ierr);
      ublas::vector<double> t(9,0);
      ierr = GetFExt(*fit,&t.data()[0],NULL,NULL); CHKERRQ(ierr);
      ierr = GetTangentExt(*fit,&t.data()[0],NULL,NULL); CHKERRQ(ierr);
    }*/

    PetscFunctionReturn(0);
  }
  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	NumeredDofMoFEMEntity_multiIndex& numered_dofs_rows = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_rows);
	NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator niit,hi_niit;
	niit = numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
	hi_niit= numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(problem_ptr->get_nb_local_dofs_row());
	PetscScalar *array;
	ierr = VecGetArray(snes_f,&array); CHKERRQ(ierr);
	Range local_ents;
	for(;niit!=hi_niit;niit++) {
	  local_ents.insert(niit->get_ent());
	}
	rval = moab.tag_set_data(th_res,local_ents,array); CHKERR_PETSC(rval);
	ierr = VecRestoreArray(snes_f,&array); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(*snes_A,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*snes_A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	//Matrix View
	//MatView(*snes_A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	//std::string wait;
	//std::cin >> wait;
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");

    }

    ierr = PetscGetTime(&v2); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscFunctionReturn(0);
  }

  };

  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet2.size());

  const double YoungModulus = 1;
  const double PoissonRatio = 0.25;
  ElasticFEMethod MyFE(moab,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2);

  moabSnesCtx SnesCtx(mField,"ELASTIC_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("ELASTIC",&MyFE));

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  //ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  //PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //PetscFinalize();
  //return 0;

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  /*KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  */

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // Write Displacements DOFS on Vertices
  struct MyEntMethod: public moabField::EntMethod {
    ErrorCode rval;
    PetscErrorCode ierr;

    Tag th_disp;
    MyEntMethod(Interface& _moab): EntMethod(_moab) {
      double def_VAL[3] = {0,0,0};
      // create TAG
      rval = moab.tag_get_handle("SPATIAL_POSITIONS_VAL",3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
    }

    vector<double> vals;
    Range nodes;
    
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscPrintf(PETSC_COMM_WORLD,"Start postporcess\n");
      rval = moab.get_entities_by_type(0,MBVERTEX,nodes); CHKERR_PETSC(rval);
      vals.resize(nodes.size()*3);
      PetscFunctionReturn(0);
    }
    
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      if(dof_ptr->get_ent_type()!=MBVERTEX) PetscFunctionReturn(0);
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double fval = dof_ptr->get_FieldData();
      Range::iterator nit = find(nodes.begin(),nodes.end(),ent);
      if(nit==nodes.end()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      unsigned int pos = std::distance(nodes.begin(),nit);
      pos = 3*pos+dof_rank;
      if(pos>vals.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      vals[pos] = fval;
      //cerr << pos << " --> " << fval << " ent " << ent << endl;
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = moab.tag_set_data(th_disp,nodes,&vals[0]); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,"End postporcess\n");
      PetscFunctionReturn(0);
    }

  };

  MyEntMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  struct PostProcFEMethod: public FEMethod_UpLevelStudent {

    Interface& moab_post_proc;
    Core mb_instance_post_proc;
    Interface& moab_ref;
    Core mb_instance_ref;

    const int max_level;
    vector<EntityHandle> meshset_level;
    bool init_ref;

    ParallelComm* pcomm;
    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;

    PostProcFEMethod(Interface& _moab): FEMethod_UpLevelStudent(_moab,1),
      moab_post_proc(mb_instance_post_proc),moab_ref(mb_instance_ref),
      max_level(2),init_ref(false) {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

      meshset_level.resize(max_level+1);

    }

    Tag th_disp,th_strain;
    vector<double> g_NTET;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start PostProc\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

      if(init_ref) PetscFunctionReturn(0);
      
      double base_coords[] = {
	0,0,0,
	1,0,0,
	0,1,0,
	0,0,1 };
      
      EntityHandle nodes[4];
      for(int nn = 0;nn<4;nn++) {
	rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
      }
      EntityHandle tet;
      rval = moab_ref.create_element(MBTET,nodes,4,tet); CHKERR_PETSC(rval);

      moabField_Core core_ref(moab_ref);
      moabField& mField_ref = core_ref;
      ierr = mField_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

      for(int ll = 0;ll<max_level;ll++) {
	PetscPrintf(PETSC_COMM_WORLD,"Refine Level %d\n",ll);
	rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = mField_ref.refine_get_ents(BitRefLevel().set(ll),meshset_level[ll]); CHKERRQ(ierr);
	ierr = mField_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = mField_ref.refine_TET(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = mField_ref.refine_get_ents(BitRefLevel().set(max_level),meshset_level[max_level]); CHKERRQ(ierr);

      std::vector<double> ref_coords;
      rval = moab_ref.get_vertex_coordinates(ref_coords); CHKERR_PETSC(rval);
      g_NTET.resize(4*ref_coords.size()/3);
      ShapeMBTET(&g_NTET[0],&ref_coords[0],&ref_coords[ref_coords.size()/3],&ref_coords[2*ref_coords.size()/3],ref_coords.size()/3);

      double def_VAL[9] = {0,0,0, 0,0,0, 0,0,0};
      // create TAG
      rval = moab_post_proc.tag_get_handle("SPATIAL_POSITIONS_VAL",3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
      rval = moab_post_proc.tag_get_handle("STRAIN_VAL",9,MB_TYPE_DOUBLE,th_strain,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);

      init_ref = true;

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      map<EntityHandle,EntityHandle> node_map;
      Range::iterator nit = ref_nodes.begin();
      for(int nn = 0;nit!=ref_nodes.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      Range ref_tets;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBTET,ref_tets); CHKERR_PETSC(rval);
      Range::iterator tit = ref_tets.begin();
      for(;tit!=ref_tets.end();tit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*tit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_tet;
	rval = moab_post_proc.create_element(MBTET,conn_post_proc,4,ref_tet); CHKERR_PETSC(rval);
      }

      Data_at_Gauss_pt::iterator diit = data_at_gauss_pt.find("SPATIAL_POSITION");
      if(diit==data_at_gauss_pt.end()) SETERRQ(PETSC_COMM_SELF,1,"no SPATIAL_POSITION !!!");
      vector< ublas::vector<FieldData> > &data = diit->second;
      vector< ublas::vector<FieldData> >::iterator vit = data.begin();
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;vit!=data.end();vit++,mit++) {
	rval = moab_post_proc.tag_set_data(th_disp,&mit->second,1,&vit->data()[0]); CHKERR_PETSC(rval);
      }

      //Strains to Noades in PostProc Mesh
      vector< ublas::matrix< FieldData > > GradU_at_GaussPt;
      ierr = GetGaussDiffDataVector("SPATIAL_POSITION",GradU_at_GaussPt); CHKERRQ(ierr);
      vector< ublas::matrix< FieldData > >::iterator viit = GradU_at_GaussPt.begin();
      mit = node_map.begin();
      for(;viit!=GradU_at_GaussPt.end();viit++,mit++) {
	ublas::matrix< FieldData > Strain = 0.5*( (*viit) + trans(*viit) );
	rval = moab_post_proc.tag_set_data(th_strain,&mit->second,1,&(Strain.data()[0])); CHKERR_PETSC(rval);
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End PostProc: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      ParallelComm* pcomm_post_proc = ParallelComm::get_pcomm(&moab_post_proc,MYPCOMM_INDEX);
      if(pcomm_post_proc == NULL) pcomm_post_proc =  new ParallelComm(&moab_post_proc,PETSC_COMM_WORLD);
      for(unsigned int rr = 1; rr<pcomm_post_proc->size();rr++) {
	Range tets;
	rval = moab_post_proc.get_entities_by_type(0,MBTET,tets); CHKERR_PETSC(rval);
	rval = pcomm_post_proc->broadcast_entities(rr,tets); CHKERR(rval);
      }
      PetscFunctionReturn(0);
    }

  };

  PostProcFEMethod fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  //ierr = KSPDestroy(&solver); CHKERRQ(ierr);


  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

