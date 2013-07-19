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

#include "ElasticFEMethod.hpp"
#include "PostProcDisplacementOnMesh.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
#include "ElasticFEMethodForInterface.hpp"
#include "ArcLeghtTools.hpp"

#ifdef __cplusplus
extern "C" {
#endif
#include <petsc-private/snesimpl.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

struct ArcInterfaceElasticFEMethod: public InterfaceElasticFEMethod {

  ArcInterfaceElasticFEMethod(Interface& _moab): InterfaceElasticFEMethod(_moab) {};

  ArcInterfaceElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3): 
      InterfaceElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2,_SideSet3) {};


  PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
      ierr = GetMatrices(); CHKERRQ(ierr);

      //Dirihlet Boundary Condition
      ApplyDirihletBC();
      if(Diagonal!=PETSC_NULL) {
	if(DirihletBC.size()>0) {
	  DirihletBCDiagVal.resize(DirihletBC.size());
	  fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	  ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	}
      }

      //Assembly Aij and F
      ierr = RhsAndLhs(); CHKERRQ(ierr);

      //Neumann Boundary Conditions
      ierr = NeumannBC(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);

    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
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
    ierr = VecDuplicate(F,&Diagonal); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
    ierr = MatDiagonalSet(Aij,Diagonal,ADD_VALUES); CHKERRQ(ierr);
    ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
    // Note MAT_FLUSH_ASSEMBLY
    ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = PetscGetTime(&v2); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscFunctionReturn(0);
  }

};





}
