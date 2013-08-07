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

#ifndef __LINEAR_DYNAMICS_ELASTICITY_HPP
#define __LINEAR_DYNAMICS_ELASTICITY_HPP

#include "moabField.hpp"
#include "moabField_Core.hpp"
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "moabTs.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

  /*
   * M*u'' + K*u' - F = 0
   *
   * F( t, [ dot_u, u], [ dot_u', u'] ) = [ 0 -1 ][ dot_u' ] + [ 1 0 ][ dot_u ] + [ 0    ] = [ 0 ]
   *                                      [ M  0 ][ u'     ]   [ 0 K ][ u     ]   [ F(t) ]   [ 0 ]
   * 
   * Fu  = [ 1 0 ][ dot_u ]
   *       [ 0 K ][ u     ] 
   * 
   * Fu' = [ 0 -1 ][ dot_u' ]
   *       [ M  0 ][ u'     ]
   *
   * G = Fu + a*Fu'
   *
   * dG/(du,ddot_u) = [ -1 0 ] + a*[ 0 1 ]
   *                  [  0 K ] +   [ M 0 ]
   *
   */

  struct DynamicElasticFEMethod: public ElasticFEMethod {

    struct BC {
      const bool Dirihlet_BC_on_SideSet2;
      BC(): Dirihlet_BC_on_SideSet2(false),final_time(0) {}
      const double final_time;
      BC(const double _final_time): Dirihlet_BC_on_SideSet2(true), final_time(_final_time) {}

      virtual PetscErrorCode f_CalcTraction(double ts_t,ublas::vector<double,ublas::bounded_array<double,3> >& traction) const = 0;
      virtual PetscErrorCode f_CalcDisp(double ts_t,
	 ublas::vector<double,ublas::bounded_array<double,3> >& disp,
	 ublas::vector<double,ublas::bounded_array<double,3> >& vel) const = 0;

    };


    moabField& mField;
    PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method;
    double rho;
    const int debug;
    BC *bc;

    DynamicElasticFEMethod(Interface& _moab,moabField& _mField,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,double _rho,Range &_SideSet1,Range &_SideSet2,BC *_bc): 
      ElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,
      _SideSet1,_SideSet2),mField(_mField),fe_post_proc_method(moab,_lambda,_mu),rho(_rho),debug(1),bc(_bc) {
      rval = moab.get_connectivity(SideSet2,SideSet2Nodes,true); CHKERR_THROW(rval);

      PetscInt ghosts[1] = { 0 };
      if(pcomm->rank() == 0) {
	VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostU);
	VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&GhostK);
      } else {
	VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostU);
	VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&GhostK);
      }

      mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&u_by_row);

      rval = moab.get_adjacencies(SideSet2,1,false,SideSet2Edges,Interface::UNION); CHKERR_THROW(rval);
      SideSet2_.insert(SideSet2.begin(),SideSet2.end());
      SideSet2_.insert(SideSet2Edges.begin(),SideSet2Edges.end());
      SideSet2_.insert(SideSet2Nodes.begin(),SideSet2Nodes.end());

    };
  
    ~DynamicElasticFEMethod() {
      //VecDestroy(&GhostU);
      //VecDestroy(&GhostK);
      //VecDestroy(&u_by_row);
    }

    Range SideSet2Nodes,SideSet2Edges,SideSet2_;
    Vec GhostU,GhostK;
    Vec u_by_row;
    set<UId> DirihletBC_SideSet2;

    PetscErrorCode ApplyDirihletBC() {
      PetscFunctionBegin;
      ierr = ElasticFEMethod::ApplyDirihletBC(); CHKERRQ(ierr);
      
      DirihletBC_SideSet2.clear();
      if(bc->Dirihlet_BC_on_SideSet2) {

	if(ts_t > bc->final_time) PetscFunctionReturn(0);

	Range::iterator siit2 = SideSet2_.begin();
	for(;siit2!=SideSet2_.end();siit2++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit2);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit2);
	  for(;riit!=hi_riit;riit++) {
	    if(riit->get_name()!="DISPLACEMENT") continue;
	    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	    for(int rr = 0;rr<row_mat;rr++) {
	      vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlob[rr].end() ) {
		if( it == RowGlob[rr].end() ) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
		DirihletBC_SideSet2.insert(riit->get_unique_id());
		*it = -1; // of idx is set -1 row is not assembled
	      }
	    }
	  }
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator ciit = col_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit2);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_ciit = col_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit2);
	  for(;ciit!=hi_ciit;ciit++) {
	    if(ciit->get_name()!="DISPLACEMENT") continue;
	    for(int cc = 0;cc<col_mat;cc++) {
	      vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	      if( it != ColGlob[cc].end() ) {
		*it = -1; // of idx is set -1 column is not assembled
	      }
	    }
	  }
	}

      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode DirihletBCDiagonalSet() {
      PetscFunctionBegin;
      Range DirihletSet;
      DirihletSet.insert(SideSet1_.begin(),SideSet1_.end());

      if(bc->Dirihlet_BC_on_SideSet2) {
	if(ts_t < bc->final_time) {
	  DirihletSet.insert(SideSet2_.begin(),SideSet2_.end());
	}
      }

      NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
      dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
      hi_dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(problem_ptr->get_nb_local_dofs_row());
      for(;dit!=hi_dit;dit++) {
	if(find(DirihletSet.begin(),DirihletSet.end(),dit->get_ent()) == DirihletSet.end()) continue;
	ierr = MatSetValue(*ts_B,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),ts_a,INSERT_VALUES); CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode DirihletBCSet() {
      PetscFunctionBegin;

      if(bc->Dirihlet_BC_on_SideSet2) {

	if(ts_t > bc->final_time) PetscFunctionReturn(0);

	ublas::vector<double,ublas::bounded_array<double,3> > disp(3);
	ublas::vector<double,ublas::bounded_array<double,3> > vel(3);
	ierr = bc->f_CalcDisp(ts_t,disp,vel); CHKERRQ(ierr);

	NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator dit,hi_dit;
	dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().lower_bound(0);
	hi_dit = problem_ptr->numered_dofs_rows.get<PetscLocalIdx_mi_tag>().upper_bound(problem_ptr->get_nb_local_dofs_row());
	double *array,*array2,*array3;
	ierr = VecGetArray(ts_F,&array); CHKERRQ(ierr);
	ierr = VecGetArray(ts_u_t,&array2); CHKERRQ(ierr);
	ierr = VecGetArray(ts_u,&array3); CHKERRQ(ierr);

	for(;dit!=hi_dit;dit++) {
	  if(find(SideSet2_.begin(),SideSet2_.end(),dit->get_ent()) == SideSet2_.end()) continue;
	  for(int dd = 0;dd<3;dd++) {
	    if(dit->get_dof_rank()==dd) {	
	      array[dit->get_petsc_local_dof_idx()] = -(vel[dd]) + array2[dit->get_petsc_local_dof_idx()];
	    }
	  } 
	}

	dit = problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().lower_bound(0);
	hi_dit = problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().upper_bound(
	  problem_ptr->get_nb_local_dofs_col()+problem_ptr->get_nb_ghost_dofs_col());
	for(;dit!=hi_dit;dit++) {
	  if(find(SideSet2_.begin(),SideSet2_.end(),dit->get_ent()) == SideSet2_.end()) continue;
	  for(int dd = 0;dd<3;dd++) {
	    if(dit->get_dof_rank()==dd) {
	      dit->get_FieldData() = disp[dd];
	      array2[dit->get_petsc_local_dof_idx()] = vel[dd];
	    }
	  }
	}

	ierr = VecRestoreArray(ts_F,&array); CHKERRQ(ierr);
	ierr = VecRestoreArray(ts_u_t,&array2); CHKERRQ(ierr);


      }


      PetscFunctionReturn(0);

    }

    /// Set Neumann Boundary Conditions on SideSet2
    PetscErrorCode NeumannBC(Vec F_ext) {
      PetscFunctionBegin;
      if(bc->Dirihlet_BC_on_SideSet2) PetscFunctionReturn(0);

      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
      ierr = bc->f_CalcTraction(ts_t,traction); CHKERRQ(ierr);

      //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
      ierr = ElasticFEMethod::NeumannBC(F_ext,traction,SideSet2); CHKERRQ(ierr);
    
      PetscFunctionReturn(0);
    }

    //This is for L2 space 
    vector<DofIdx> VelRowGlob;
    vector<DofIdx> VelRowLocal;
    vector<DofIdx> VelColGlob;
    vector<DofIdx> VelColLocal;
    vector<ublas::matrix<FieldData> > VelRowNMatrix;
    vector<ublas::matrix<FieldData> > VelColNMatrix;
 
    PetscErrorCode GetMatricesVelocities() {
      PetscFunctionBegin;
      VelRowGlob.clear();
      VelRowLocal.clear();
      VelColGlob.clear();
      VelColLocal.clear();
      VelRowNMatrix.clear();
      VelColNMatrix.clear();
      //ROWs
      ierr = GetRowGlobalIndices("VELOCITIES",MBTET,VelRowGlob); CHKERRQ(ierr);
      ierr = GetRowLocalIndices("VELOCITIES",MBTET,VelRowLocal); CHKERRQ(ierr);
      ierr = GetGaussRowNMatrix("VELOCITIES",MBTET,VelRowNMatrix); CHKERRQ(ierr);
      //COLs
      ierr = GetColGlobalIndices("VELOCITIES",MBTET,VelColGlob); CHKERRQ(ierr);
      ierr = GetColLocalIndices("VELOCITIES",MBTET,VelColLocal); CHKERRQ(ierr);
      ierr = GetGaussColNMatrix("VELOCITIES",MBTET,VelColNMatrix); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    ublas::vector<FieldData> accelerations;
    ublas::vector<FieldData> dot_velocities;
    vector<ublas::vector<FieldData> > displacements;
    vector<ublas::vector<FieldData> > velocities;
    PetscErrorCode GetVelocities_Form_TS_u_t() {
      PetscFunctionBegin;
      accelerations.clear();
      dot_velocities.clear();
      displacements.clear();
      velocities.clear();
      Vec u_local;
      ierr = VecGhostGetLocalForm(ts_u,&u_local); CHKERRQ(ierr);
      Vec u_t_local;
      ierr = VecGhostGetLocalForm(ts_u_t,&u_t_local); CHKERRQ(ierr);
      int local_size;
      ierr = VecGetLocalSize(u_t_local,&local_size); CHKERRQ(ierr);
      double *array,*array2;
      ierr = VecGetArray(u_t_local,&array); CHKERRQ(ierr);
      ierr = VecGetArray(u_local,&array2); CHKERRQ(ierr);
      accelerations.resize(VelColLocal.size());
      dot_velocities.resize(VelColLocal.size());
      int ii = 0;
      vector<DofIdx>::iterator it = VelColLocal.begin();
      for(;it!=VelColLocal.end();it++,ii++) {
	if(*it < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	if(*it >= local_size) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	accelerations[ii] = array[*it];
	dot_velocities[ii] = array2[*it];
      }
      int cc = 0;
      velocities.resize(col_mat);
      displacements.resize(col_mat);
      for(;cc<col_mat;cc++) {
	if(ColLocal[cc].size()==0) continue;
	velocities[cc].resize(ColLocal[cc].size());
	displacements[cc].resize(ColLocal[cc].size());
	vector<DofIdx>::iterator iit = ColLocal[cc].begin();
	for(int iii = 0;iit!=ColLocal[cc].end();iit++,iii++) {
	  if(*iit < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(*iit >= local_size) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  (velocities[cc])[iii] = array[*iit];
	  NumeredDofMoFEMEntity_multiIndex::index<PetscLocalIdx_mi_tag>::type::iterator ciit,hi_ciit;
	  ciit = problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().find(*iit);
	  hi_ciit = problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().end();
	  if(ciit==hi_ciit) {
	    for(ciit =  problem_ptr->numered_dofs_cols.get<PetscLocalIdx_mi_tag>().begin();
	      ciit!=hi_ciit;ciit++) {
	      cerr << *it << " " << ColLocal.size() << endl;
	      cerr << *ciit << endl;
	    }
	    SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  }
	  (displacements[cc])[iii] = ciit->get_FieldData();
	}
      }
      ierr = VecRestoreArray(u_local,&array2); CHKERRQ(ierr);
      ierr = VecRestoreArray(u_t_local,&array); CHKERRQ(ierr);
      ierr = VecGhostRestoreLocalForm(ts_u,&u_local); CHKERRQ(ierr);
      ierr = VecGhostRestoreLocalForm(ts_u_t,&u_t_local); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    ublas::vector<ublas::matrix<FieldData> > Mass;
    PetscErrorCode MassLhs() {
      PetscFunctionBegin;
      Mass.resize(row_mat);
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	  ublas::matrix<FieldData> &col_Mat = VelColNMatrix[gg];
	  double w = rho*V*G_TET_W45[gg];
	  if(gg == 0) {
	    Mass[rr] = w*prod( trans(row_Mat), col_Mat );
	  } else {
	    Mass[rr] += w*prod( trans(row_Mat), col_Mat );
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    ublas::matrix<FieldData> VV;
    PetscErrorCode VVLhs() {
      PetscFunctionBegin;
      int g_dim = g_NTET.size()/4;
      for(int gg = 0;gg<g_dim;gg++) {
	ublas::matrix<FieldData> &row_Mat = VelRowNMatrix[gg];
	ublas::matrix<FieldData> &col_Mat = VelColNMatrix[gg];
	double w = rho*V*G_TET_W45[gg];
	if(gg == 0) {
	  VV = w*prod( trans(row_Mat), col_Mat );
	} else {
	  VV += w*prod( trans(row_Mat), col_Mat );
	}
      }
      PetscFunctionReturn(0);
    }

    ublas::vector<ublas::matrix<FieldData> > VU;
    PetscErrorCode VULhs() {
      PetscFunctionBegin;
      VU.resize(col_mat);
      int g_dim = g_NTET.size()/4;
      for(int cc = 0;cc<col_mat;cc++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = VelRowNMatrix[gg];
	  ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	  double w = rho*V*G_TET_W45[gg];
	  if(gg == 0) {
	    VU[cc] = -w*prod( trans(row_Mat), col_Mat );
	  } else {
	    VU[cc] -= w*prod( trans(row_Mat), col_Mat );
	  }
	}
      }
      PetscFunctionReturn(0);
    }

    // STIFFNESS - run first
    // COPUPLING_VU - last
    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      g_NTET.resize(4*45);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 

      if(fe_name=="STIFFNESS") {
	// See FEAP - - A Finite Element Analysis Program
	D_lambda = ublas::zero_matrix<FieldData>(6,6);
	for(int rr = 0;rr<3;rr++) {
	  ublas::matrix_row<ublas::matrix<FieldData> > row_D_lambda(D_lambda,rr);
	  for(int cc = 0;cc<3;cc++) {
	    row_D_lambda[cc] = 1;
	  }
	}
	D_mu = ublas::zero_matrix<FieldData>(6,6);
	for(int rr = 0;rr<6;rr++) {
	  D_mu(rr,rr) = rr<3 ? 2 : 1;
	}
	D = lambda*D_lambda + mu*D_mu;
	
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    PetscReal ftime;
	    ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
	    PetscInt steps,snesfails,rejects,nonlinits,linits;
	    ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
	    ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
	    ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
	    ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
	    ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,
	      "\tsteps %D (%D rejected, %D SNES fails), ftime %G, nonlinits %D, linits %D\n",steps,rejects,snesfails,ftime,nonlinits,linits);
	    ierr = mField.set_global_VecCreateGhost(problem_ptr->get_name(),Col,ts_u,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = mField.set_local_VecCreateGhost(problem_ptr->get_name(),Row,u_by_row,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    NumeredDofMoFEMEntity_multiIndex &numered_dofs_cols = const_cast<NumeredDofMoFEMEntity_multiIndex&>(problem_ptr->numered_dofs_cols);
	    NumeredDofMoFEMEntity_multiIndex::index<FieldName_mi_tag>::type::iterator lit;
	    Range::iterator nit = SideSet2Nodes.begin();
	    for(;nit!=SideSet2Nodes.end();nit++) {
	      NumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator dit,hi_dit;
	      dit = numered_dofs_cols.get<MoABEnt_mi_tag>().lower_bound(*nit);
	      hi_dit = numered_dofs_cols.get<MoABEnt_mi_tag>().upper_bound(*nit);
	      double _coords_[3];
	      rval = moab.get_coords(&*nit,1,_coords_);  CHKERR_THROW(rval);
	      for(;dit!=hi_dit;dit++) {
		PetscPrintf(PETSC_COMM_WORLD,"%s [ %d ] %6.4e ",dit->get_name().c_str(),dit->get_dof_rank(),dit->get_FieldData());
		PetscPrintf(PETSC_COMM_WORLD,"-> %3.4f %3.4f %3.4f ",_coords_[0],_coords_[1],_coords_[2]);
		PetscPrintf(PETSC_COMM_WORLD,"-> time %6.4e\n",ftime);
	      }
	    }
	    if(steps%10==0) { 
	      rval = fe_post_proc_method.moab_post_proc.delete_mesh(); CHKERR_PETSC(rval);
	      ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","STIFFNESS",fe_post_proc_method);  CHKERRQ(ierr);
	      if(pcomm->rank()==0) {
		ostringstream sss;
		sss << (int)(ftime*1e3) << "_out.vtk";
		rval = fe_post_proc_method.moab_post_proc.write_file(sss.str().c_str(),"VTK",""); CHKERR_PETSC(rval);
	      }
	    }
	    ierr = VecZeroEntries(GhostU); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(GhostU,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostU,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    } break;
	  case ctx_TSSetRHSFunction:
	  case ctx_TSSetIFunction:
	    ierr = VecZeroEntries(ts_F); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(ts_F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(ts_F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    //
	    ierr = DirihletBCSet(); CHKERRQ(ierr);
	    //
	    break;
	  case ctx_TSSetRHSJacobian:
	  case ctx_TSSetIJacobian:
	    ierr = MatZeroEntries(*ts_B); CHKERRQ(ierr);
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	
	PetscFunctionReturn(0);
      }
      if(fe_name=="MASS") {
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VV") {
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    ierr = VecZeroEntries(GhostK); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(GhostK,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostK,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  } break;
	  default: {
	  } break;
	}
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VU") {
	PetscFunctionReturn(0);
      }

      SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      if(fe_name=="STIFFNESS") {
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    ierr = VecGhostUpdateBegin(GhostU,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostU,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(GhostU); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(GhostU); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(GhostU,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostU,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    } break;
	  default: {
	    } break;
	}
	PetscFunctionReturn(0);
      }
      if(fe_name=="MASS") {
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VV") {
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    ierr = VecGhostUpdateBegin(GhostK,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostK,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(GhostK); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(GhostK); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(GhostK,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(GhostK,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    double *array,*array2;
	    ierr = VecGetArray(GhostU,&array); CHKERRQ(ierr);
	    ierr = VecGetArray(GhostK,&array2); CHKERRQ(ierr);
	    PetscReal ftime;
	    ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
	    PetscPrintf(PETSC_COMM_WORLD,
	      "Time  = %6.4e Potential Energy = %6.4e Kinematic Energy = %6.4e Total = %6.4e\n",
	      ftime,array[0],array2[0],array[0]+array2[0]);
	    ierr = VecRestoreArray(GhostU,&array); CHKERRQ(ierr);
	    ierr = VecRestoreArray(GhostK,&array2); CHKERRQ(ierr);
	    } break;
	  default: {
	    } break;
	}
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VU") {
	switch (ts_ctx) {
	  case ctx_TSSetRHSFunction:
	  case ctx_TSSetIFunction:
	    ierr = VecGhostUpdateBegin(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
	    break;
	  case ctx_TSSetRHSJacobian:
	  case ctx_TSSetIJacobian: {
	    ierr = MatAssemblyBegin(*ts_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*ts_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	    //
	    ierr = DirihletBCDiagonalSet(); CHKERRQ(ierr);
	    //
	    ierr = MatAssemblyBegin(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    *ts_flag = SAME_NONZERO_PATTERN; 
	    //Matrix View
	    //MatView(*ts_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	    //std::string wait;
	    //std::cin >> wait;
	    } break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}

	PetscFunctionReturn(0);
      }

      SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
     PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      if(ts_F!=PETSC_NULL) {
	//If index is set to -1 ingonre its assembly
	VecSetOption(ts_F, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); 
      }

      if(fe_name=="STIFFNESS") {
	ierr = GetMatrices(); CHKERRQ(ierr);
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    ierr = Fint(); CHKERRQ(ierr);
	    Vec u_local;
	    ierr = VecGhostGetLocalForm(u_by_row,&u_local); CHKERRQ(ierr);
	    double *array2;
	    ierr = VecGetArray(u_local,&array2); CHKERRQ(ierr);
	    double *array_U;
	    ierr = VecGetArray(GhostU,&array_U); CHKERRQ(ierr);
	    for(int rr = 0;rr<row_mat;rr++) {	
	      ublas::vector<FieldData> displacements;
	      displacements.resize(RowLocal[rr].size());
	      vector<DofIdx>::iterator iit = RowLocal[rr].begin();
	      for(int iii = 0;iit!=RowLocal[rr].end();iit++,iii++) {
		displacements[iii] = array2[*iit];
	      }
	      array_U[0] += 0.5*inner_prod(f_int[rr],displacements);
	    }
	    ierr = VecRestoreArray(u_local,&array2); CHKERRQ(ierr);
	    ierr = VecRestoreArray(GhostU,&array_U); CHKERRQ(ierr);
	    ierr = VecGhostRestoreLocalForm(ts_u,&u_local); CHKERRQ(ierr);
	    } break;
	  case ctx_TSSetRHSFunction: {
	    } break;
	  case ctx_TSSetRHSJacobian: {
	    } break;
	  case ctx_TSSetIFunction: {
	    ierr = NeumannBC(ts_F); CHKERRQ(ierr);
	    ierr = Fint(ts_F); CHKERRQ(ierr);
	  } break;
	  case ctx_TSSetIJacobian: {
	    ierr = Stiffness(); CHKERRQ(ierr);
	    for(int rr = 0;rr<row_mat;rr++) {
	      if(RowGlob[rr].size()==0) continue;
	      for(int cc = 0;cc<col_mat;cc++) {
		if(ColGlob[cc].size()==0) continue;
		if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
		if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
		ierr = MatSetValues(*ts_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	      }
	    }
	    } break;
  	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	ierr = OpStudentEnd(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
      }
      if(fe_name=="MASS") {
	for(int cc = 0;cc<col_mat;cc++) {
	  ColGlob[cc].resize(0);
	  ColLocal[cc].resize(0);
	}
	ierr = GetMatricesRows(); CHKERRQ(ierr);
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	ierr = MassLhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetRHSFunction: {
	    } break;
	  case ctx_TSSetRHSJacobian: {
	    } break;
	  case ctx_TSSetIFunction: {
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    for(int rr = 0;rr<row_mat;rr++) {
	      if(RowGlob[rr].size()==0) continue;
	      ublas::vector<FieldData> Mu_t = prod(Mass[rr],accelerations);
	      if(RowGlob[rr].size()!=Mu_t.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(ts_F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(Mu_t.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
	  case ctx_TSSetIJacobian: {
	    for(int rr = 0;rr<row_mat;rr++) {
	      if(RowGlob[rr].size()==0) continue;
	      if(RowGlob[rr].size()!=Mass[rr].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(VelColGlob.size()!=Mass[rr].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      Mass[rr] *= ts_a;
	      ierr = MatSetValues(*ts_B,RowGlob[rr].size(),&(RowGlob[rr])[0],VelColGlob.size(),&VelColGlob[0],&(Mass[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
  	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	ierr = OpStudentEnd(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VV") {
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = VVLhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSTSMonitorSet: {
	    Vec u_local;
	    ierr = VecGhostGetLocalForm(ts_u,&u_local); CHKERRQ(ierr);
	    double *array;
	    ierr = VecGetArray(u_local,&array); CHKERRQ(ierr);
	    ublas::vector<FieldData> velocities_row,velocities_col;
	    velocities_col.resize(VelColLocal.size());
	    vector<DofIdx>::iterator it = VelColLocal.begin();
	    for(int ii = 0;it!=VelColLocal.end();it++,ii++) {
	      velocities_col[ii] = array[*it];
	    }
	    velocities_row.resize(VelRowLocal.size());
	    it = VelRowLocal.begin();
	    for(int ii = 0;it!=VelRowLocal.end();it++,ii++) {
	      velocities_row[ii] = array[*it];
	    }
	    ierr = VecRestoreArray(u_local,&array); CHKERRQ(ierr);
	    ierr = VecGhostRestoreLocalForm(ts_u,&u_local); CHKERRQ(ierr);
	    ublas::vector<FieldData> VVvelocities_col = prod(VV,velocities_col);  
	    double *array_K;
	    ierr = VecGetArray(GhostK,&array_K); CHKERRQ(ierr);
	    array_K[0] += 0.5*inner_prod(velocities_row,VVvelocities_col);
	    ierr = VecRestoreArray(GhostK,&array_K); CHKERRQ(ierr);
	    } break;
	  case ctx_TSSetRHSFunction: {
	    } break;
	  case ctx_TSSetRHSJacobian: {
	    } break;
	  case ctx_TSSetIFunction: {
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    if(VV.size2()!=dot_velocities.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ublas::vector<FieldData> VVu = prod(VV,dot_velocities);
	    if(VelRowGlob.size()!=VVu.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    //cerr << "VVu " << VVu << endl;
	    ierr = VecSetValues(ts_F,VelRowGlob.size(),&VelRowGlob[0],&(VVu.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	  } break;
	  case ctx_TSSetIJacobian: {
	    if(VelRowGlob.size()!=VV.size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    if(VelColGlob.size()!=VV.size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = MatSetValues(*ts_B,VelRowGlob.size(),&(VelRowGlob)[0],VelColGlob.size(),&VelColGlob[0],&(VV.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  } break;
  	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	ierr = OpStudentEnd(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
      }
      if(fe_name=="COPUPLING_VU") {
	ierr = GetMatricesCols(); CHKERRQ(ierr);
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	ierr = VULhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetRHSFunction: {
	    } break;
	  case ctx_TSSetRHSJacobian: {
	    } break;
	  case ctx_TSSetIFunction: {
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    for(int cc = 0;cc<col_mat;cc++) {
	      if(VU[cc].size2()!=velocities[cc].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ublas::vector<FieldData> VUu_t = prod(VU[cc],velocities[cc]);
	      if(VelRowGlob.size()!=VUu_t.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(ts_F,VelRowGlob.size(),&VelRowGlob[0],&(VUu_t.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
	  case ctx_TSSetIJacobian: {
	    for(int cc = 0;cc<col_mat;cc++) {
	      if(VelRowGlob.size()!=VU[cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(ColGlob[cc].size()!=VU[cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      VU[cc] *= ts_a;
	      ierr = MatSetValues(*ts_B,VelRowGlob.size(),&(VelRowGlob)[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(VU[cc].data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
  	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	ierr = OpStudentEnd(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
      }
      SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
      PetscFunctionReturn(0);
    }

  };

#endif //__LINEAR_DYNAMICS_ELASTICITY_HPP


