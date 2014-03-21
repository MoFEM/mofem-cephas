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

#ifndef __ElasticFEMethod_NewBCs__
#define __ElasticFEMethod_NewBCs__

#include <boost/numeric/ublas/symmetric.hpp>
#include "ElasticFEMethod.hpp"

namespace MoFEM {
    
    
struct ElasticFEMethod_NewBCs: public ElasticFEMethod{
    
    ElasticFEMethod_NewBCs(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,
double _lambda,double _mu): ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu){};
    
    PetscErrorCode postProcess() {
        PetscFunctionBegin;
        // Note MAT_FLUSH_ASSEMBLY
        ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_ZerosRowsColumns(this,Aij); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_RHS(this,F); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
        ierr = PetscTime(&v2); CHKERRQ(ierr);
        ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        PetscFunctionReturn(0);
    }

    
    PetscErrorCode operator()() {
        PetscFunctionBegin;
        ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
        ierr = GetMatrices(); CHKERRQ(ierr);
        
        //Assembly Aij and F
        ierr = RhsAndLhs(); CHKERRQ(ierr);
        
        //Neumann Boundary Conditions
        ierr = NeumannBC(F); CHKERRQ(ierr);
        
        ierr = OpStudentEnd(); CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
};
    
}

#endif //__ElasticFEMethod_NewBCs__
