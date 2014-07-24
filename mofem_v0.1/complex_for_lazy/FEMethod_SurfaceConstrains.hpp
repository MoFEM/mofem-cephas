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

#ifndef __MOABSURFACECONSTRAINS_HPP__
#define __MOABSURFACECONSTRAINS_HPP__

namespace MoFEM {

struct ConstrainSurfacGeometry:public FieldInterface::FEMethod {
  ErrorCode rval;
  PetscErrorCode ierr;
  FieldInterface& mField;
  Mat C;
  string lambdaFieldName;

  vector<double> diffNTRI;
  vector<double> g_NTRI;
  const double *G_TRI_W;

  Tag thProjection;
  Range crackFrontEdgesNodes;
  bool useProjectionFromCrackFront;

  ConstrainSurfacGeometry(FieldInterface& _mField,Mat _C,string _lambdaFieldName,int _verbose = 0);
  ConstrainSurfacGeometry(FieldInterface& _mField,Mat _C,int _verbose = 0);
  void runInConstructor();

  PetscErrorCode preProcess();

  ublas::bounded_matrix<double,3,9 > C_MAT_ELEM,iC_MAT_ELEM;
  ublas::bounded_matrix<double,9,3 > CT_MAT_ELEM;
  ublas::bounded_matrix<double,3,9 > dC_MAT_ELEM;
  ublas::bounded_matrix<double,9,9 > dCT_MAT_ELEM;
  ublas::vector<double,ublas::bounded_array<double,9> > ig_VEC_ELEM;
  ublas::vector<double,ublas::bounded_array<double,9> > if_VEC_ELEM;

  ublas::vector<DofIdx,ublas::bounded_array<DofIdx,3> > lambdaGlobalRowIndices,lambdaGlobalColIndices;
  vector<vector<DofIdx> > dofGlobalRowIndices,dofGlobalColIndices;
  ublas::vector<double,ublas::bounded_array<double,3> > entLambdaData;
  ublas::vector<double,ublas::bounded_array<double,9> > entDofsData,ent_idofs_data;
  ublas::vector<double,ublas::bounded_array<double,9> > cOords;

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = operator()(false,false); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  };

  EntityHandle fAce;
  PetscErrorCode cOnstrain(double *dofs_iX,double *C,double *iC,double *g,double *ig);
  virtual PetscErrorCode iNtegrate(bool transpose,bool nonlinear);
  virtual PetscErrorCode aSsemble(bool transpose,bool nonlinear);
  virtual PetscErrorCode operator()(bool transpose,bool nonlinear);
  PetscErrorCode postProcess();

};

struct ConstraunSurfaceGeometryRhs: public ConstrainSurfacGeometry {
  Vec g;
  ConstraunSurfaceGeometryRhs(FieldInterface& _mField,Vec _g,string _lambdaFieldName,int _verbose = 0); 
  ConstraunSurfaceGeometryRhs(FieldInterface& _mField,Vec _g,int _verbose = 0); 

  ublas::vector<double,ublas::bounded_array<double,3> > g_VEC_ELEM;
  ublas::vector<double,ublas::bounded_array<double,9> > f_VEC_ELEM;
  PetscErrorCode iNtegrate(bool transpose,bool nonlinear);
  PetscErrorCode aSsemble(bool transpose,bool nonlinear);
};

struct SnesConstrainSurfacGeometryTools: public FieldInterface::FEMethod {

  PetscErrorCode setElemData(FieldInterface::FEMethod &e) {
    PetscFunctionBegin;
    //copy all elem that to other FEMethod class
    e.problemPtr = problemPtr;
    e.fieldsPtr = fieldsPtr;
    e.entitiesPtr = entitiesPtr;
    e.dofsPtr = dofsPtr;
    e.finiteElementsPtr = finiteElementsPtr;
    e.finiteElementsEntitiesPtr = finiteElementsEntitiesPtr;
    e.adjacenciesPtr = adjacenciesPtr;
    e.feName = feName;
    e.fePtr = fePtr;
    e.dataPtr = dataPtr;
    e.rowPtr = rowPtr;
    e.colPtr = colPtr;
    PetscFunctionReturn(0);
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct SnesConstrainSurfacGeometry: public SnesConstrainSurfacGeometryTools {

  FieldInterface& mField;
  ConstrainSurfacGeometry matMethod;
  ConstraunSurfaceGeometryRhs vecMethod;
  bool nonlinear;
  bool useProjectionFromCrackFront;


  SnesConstrainSurfacGeometry(FieldInterface& _mField,int _verbose = 0):
   mField(_mField),
    matMethod(_mField,PETSC_NULL),
    vecMethod(_mField,snes_f),
    nonlinear(false),useProjectionFromCrackFront(false) {}

  SnesConstrainSurfacGeometry(FieldInterface& _mField,string _lambdaFieldName,int _verbose = 0):
   mField(_mField),
    matMethod(_mField,PETSC_NULL,_lambdaFieldName),
    vecMethod(_mField,snes_f,_lambdaFieldName),
    nonlinear(false),useProjectionFromCrackFront(false) {}

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    matMethod.useProjectionFromCrackFront = useProjectionFromCrackFront;
    vecMethod.useProjectionFromCrackFront = useProjectionFromCrackFront; 
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	vecMethod.g = snes_f;
	ierr = setElemData(vecMethod); CHKERRQ(ierr);
	ierr = vecMethod.operator()(true,nonlinear); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	matMethod.C = *snes_B;
	ierr = setElemData(matMethod); CHKERRQ(ierr);
	ierr = matMethod.operator()(true,nonlinear); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN: {
	ierr = matMethod.postProcess(); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }


};

}

#endif //__MOABSURFACECONSTRAINS_HPP__
