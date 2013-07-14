/** \file moabField_Core.hpp
 * \brief Core moabField::FEMethod class for user interface
 * 
 * Low level data structures not used directly by user
 *
 * Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl) <br>
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef __MOABFEMETHOD_DRIVERCOMPLEXFORLAZY_HPP__
#define __MOABFEMETHOD_DRIVERCOMPLEXFORLAZY_HPP__

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "moabSnes.hpp"
#include "moabFEMethod_ComplexForLazy.hpp"

#include "complex_for_lazy.h"

using namespace boost::numeric;

namespace MoFEM {

struct FEMethod_DriverComplexForLazy: public FEMethod_ComplexForLazy {

    enum bc_type { fixed_x = 1,fixed_y = 1<<1, fixed_z = 1<<2 };

    double t_val;
    PetscErrorCode set_t_val(double t_val_) {
      PetscFunctionBegin;
      t_val = t_val_;
      PetscFunctionReturn(0);
    }

    FEMethod_DriverComplexForLazy(Interface& _moab,
      double _lambda,double _mu,int _verbose = 0): 
      FEMethod_ComplexForLazy(_moab,FEMethod_ComplexForLazy::spatail_analysis,_lambda,_mu,_verbose) { };

    vector<DofIdx> DirihletBC;
    PetscErrorCode ApplyDirihletBC(Range &SideSet,unsigned int bc = fixed_x|fixed_y|fixed_z,bool zero_bc = true) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      if(zero_bc) DirihletBC.resize(0);
      Range::iterator siit1 = SideSet.begin();
      for(;siit1!=SideSet.end();siit1++) {
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	for(;riit!=hi_riit;riit++) {
	  if(riit->get_name()!=field_name) continue;
    	  unsigned int my_bc = 0;
	  switch(riit->get_dof_rank()) {
	    case 0: my_bc = fixed_x; break;
	    case 1: my_bc = fixed_y; break;
	    case 2: my_bc = fixed_z; break;
	    default:
	      SETERRQ(PETSC_COMM_SELF,1,"not implemented");
	  }
	  if((bc&my_bc) == 0) continue;
	  // all fixed
	  // if some ranks are selected then we could apply BC in particular direction
	  DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	  for(int cc = 0;cc<(i_last);cc++) {
	    vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),riit->get_petsc_gloabl_dof_idx());
	    if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	  }
	  for(int rr = 0;rr<(i_last);rr++) {
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
	for(int ee = 0;ee<3;ee++) {
	  it = find(FaceEdgeIndices_data[ee].begin(),FaceEdgeIndices_data[ee].end(),*dit);
	  if(it!=FaceEdgeIndices_data[ee].end()) *it = -1; // of idx is set -1 row is not assembled
	}
	it = find(FaceIndices.begin(),FaceIndices.end(),*dit);
	if(it!=FaceIndices.end()) *it = -1; // of idx is set -1 row is not assembled
      }
      PetscFunctionReturn(0);
    }

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;

  Vec Diagonal;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n");
    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	Diagonal = PETSC_NULL;
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatZeroEntries(*snes_B); CHKERRQ(ierr);
	ierr = VecDuplicate(snes_f,&Diagonal); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  vector<FieldData> DirihletBCDiagVal;
  
  PetscErrorCode CalculateFint(Vec f) {
    PetscFunctionBegin;

    switch(ctx) {
      case ctx_SNESSetFunction: { 
        ierr = GetFint(); CHKERRQ(ierr);
	VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(f,RowGlob[i_nodes].size(),&(RowGlob[i_nodes][0]),&(Fint_h.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	for(int ee = 0;ee<6;ee++) {
	  if(RowGlob[1+ee].size()>0) {
	    ierr = VecSetValues(f,RowGlob[1+ee].size(),&(RowGlob[1+ee][0]),&(Fint_h_edge_data[ee].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	for(int ff = 0;ff<4;ff++) {
	  if(RowGlob[1+6+ff].size()>0) {
	    ierr = VecSetValues(f,RowGlob[1+6+ff].size(),&(RowGlob[1+6+ff][0]),&(Fint_h_face_data[ff].data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	if(RowGlob[i_volume].size()>0) {
	  ierr = VecSetValues(f,RowGlob[i_volume].size(),&(RowGlob[i_volume][0]),&(Fint_h_volume.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode CalculateTangent(Mat B) {
    PetscFunctionBegin;

    switch(ctx) {
      case ctx_SNESSetJacobian:
	ierr = GetTangent(); CHKERRQ(ierr);
	//cerr << "Khh " << Khh << endl;
	ierr = MatSetValues(B,
	  RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	  ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	  &*(Khh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	//
	for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	    ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	    &*(Kedgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	    ColGlob[1+ee].size(),&*(ColGlob[1+ee].begin()),
	    &*(Khedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(B,
	      RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	      ColGlob[1+eee].size(),&*(ColGlob[1+eee].begin()),
	      &*(Khh_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(B,
	      RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	      ColGlob[1+6+fff].size(),&*(ColGlob[1+6+fff].begin()),
	      &*(Khh_edgeface_data(ee,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    RowGlob[1+ee].size(),&*(RowGlob[1+ee].begin()),
	    ColGlob[i_volume].size(),&*(ColGlob[i_volume].begin()),
	    &*(Khh_edgevolume_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlob[i_volume].size(),&*(RowGlob[i_volume].begin()),
	    ColGlob[1+ee].size(),&*(ColGlob[1+ee].begin()),
	    &*(Khh_volumeedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	    ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	    &*(Kfaceh_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	    ColGlob[1+6+ff].size(),&*(ColGlob[1+6+ff].begin()),
	    &*(Khface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int eee = 0;eee<6;eee++) {
	    ierr = MatSetValues(B,
	      RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	      ColGlob[1+eee].size(),&*(ColGlob[1+eee].begin()),
	      &*(Khh_faceedge_data(ff,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  for(int fff = 0;fff<4;fff++) {
	    ierr = MatSetValues(B,
	      RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	      ColGlob[1+6+fff].size(),&*(ColGlob[1+6+fff].begin()),
	      &*(Khh_faceface_data(ff,fff).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    RowGlob[1+6+ff].size(),&*(RowGlob[1+6+ff].begin()),
	    ColGlob[i_volume].size(),&*(ColGlob[i_volume].begin()),
	    &*(Khh_facevolume_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    RowGlob[i_volume].size(),&*(RowGlob[i_volume].begin()),
	    ColGlob[1+6+ff].size(),&*(ColGlob[1+6+ff].begin()),
	    &*(Khh_volumeface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	ierr = MatSetValues(B,
	  RowGlob[i_volume].size(),&*(RowGlob[i_volume].begin()),
	  ColGlob[i_nodes].size(),&*(ColGlob[i_nodes].begin()),
	  &*(Kvolumeh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlob[i_nodes].size(),&*(RowGlob[i_nodes].begin()),
	  ColGlob[i_volume].size(),&*(ColGlob[i_volume].begin()),
	  &*(Khvolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(B,
	  RowGlob[i_volume].size(),&*(RowGlob[i_volume].begin()),
	  ColGlob[i_volume].size(),&*(ColGlob[i_volume].begin()),
	  &*(Khh_volumevolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);

  }

  PetscErrorCode CaluclateFext(Vec f,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    switch(ctx) {
      case ctx_SNESSetFunction: { 
	for(;siit!=hi_siit;siit++) {
	  VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetFExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  //cerr << "FExt " << FExt << endl;
	  //cerr << "FaceNodeIndices.size() " << FaceNodeIndices.size() << endl;
	  ierr = ApplyDirihletBCFace(); CHKERRQ(ierr);
	  ierr = VecSetValues(f,FaceNodeIndices.size(),&(FaceNodeIndices[0]),&*FExt.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceEdgeIndices_data[ee].size()>0) {
	      ierr = VecSetValues(f,FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),&*FExt_edge_data[ee].data().begin(),ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(FaceIndices.size()>0) {
	    ierr = VecSetValues(f,FaceIndices.size(),&(FaceIndices[0]),&*FExt_face.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode CalculateTangentExt(Mat B,double *t,Range& NeumannSideSet) {
    PetscFunctionBegin;

    SideNumber_multiIndex& side_table = const_cast<SideNumber_multiIndex&>(fe_ent_ptr->get_side_number_table());
    SideNumber_multiIndex::nth_index<1>::type::iterator siit = side_table.get<1>().lower_bound(boost::make_tuple(MBTRI,0));
    SideNumber_multiIndex::nth_index<1>::type::iterator hi_siit = side_table.get<1>().upper_bound(boost::make_tuple(MBTRI,4));

    switch(ctx) {
      case ctx_SNESSetJacobian:
	for(;siit!=hi_siit;siit++) {
	  Range::iterator fit = find(NeumannSideSet.begin(),NeumannSideSet.end(),siit->ent);
	  if(fit==NeumannSideSet.end()) continue;
	  ierr = GetFaceIndicesAndData(siit->ent); CHKERRQ(ierr);
	  ierr = GetTangentExt(siit->ent,t,NULL,NULL); CHKERRQ(ierr);
	  ierr = ApplyDirihletBCFace(); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceNodeIndices.size(),&(FaceNodeIndices[0]),FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	    &*(KExt_hh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  for(int ee = 0;ee<3;ee++) {
	    if(FaceNodeIndices.size()==0) continue;
	    ierr = MatSetValues(B,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      &*(KExt_edgeh_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(B,
	      FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_hedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    for(int eee = 0;eee<3;eee++) {
	      ierr = MatSetValues(B,
		FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
		FaceEdgeIndices_data[eee].size(),&(FaceEdgeIndices_data[eee][0]),
		&*(KExt_edgeedge_data(ee,eee).data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    }
	    ierr = MatSetValues(B,
	      FaceIndices.size(),&(FaceIndices[0]),
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      &*(KExt_faceedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValues(B,
	      FaceEdgeIndices_data[ee].size(),&(FaceEdgeIndices_data[ee][0]),
	      FaceIndices.size(),&(FaceIndices[0]),
	      &*(KExt_edgeface_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  }
	  ierr = MatSetValues(B,
	    FaceNodeIndices.size(),&(FaceNodeIndices[0]),FaceIndices.size(),&(FaceIndices[0]),
	    &*(KExt_hface.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceIndices.size(),&(FaceIndices[0]),FaceNodeIndices.size(),&(FaceNodeIndices[0]),
	    &*(KExt_faceh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	  ierr = MatSetValues(B,
	    FaceIndices.size(),&(FaceIndices[0]),FaceIndices.size(),&(FaceIndices[0]),
	    &*(KExt_faceface.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }


    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()(Range& DirihletSideSet,Range& NeumannSideSet) {
    PetscFunctionBegin;

    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndices(); CHKERRQ(ierr);

    ierr = ApplyDirihletBC(DirihletSideSet); CHKERRQ(ierr);
    if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
    }

    double t[] = { 0,0,t_val, 0,0,t_val, 0,0,t_val };

    switch(ctx) {
      case ctx_SNESSetFunction: { 
	ierr = CalculateFint(snes_f); CHKERRQ(ierr);
	ierr = CaluclateFext(snes_f,t,NeumannSideSet); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian:
	ierr = CalculateTangent(*snes_B); CHKERRQ(ierr);
	ierr = CalculateTangentExt(*snes_B,t,NeumannSideSet); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

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
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	ierr = MatDiagonalSet(*snes_B,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	//Matrix View
	//MatView(*snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
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

}

#endif

