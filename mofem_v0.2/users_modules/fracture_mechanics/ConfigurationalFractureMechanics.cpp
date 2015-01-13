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

#ifdef WITH_TETGEN

#include <tetgen.h>
#ifdef REAL
  #undef REAL
#endif

#endif

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>
#include <FEMethod_SurfaceConstrains.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

extern "C" {
  #include <complex_for_lazy.h>
}

#include <ArcLengthTools.hpp>
#include <MatShellConstrainsByMarkAinsworth.hpp>
#include <FEMethod_ComplexForLazy.hpp>
#include <FEMethod_DriverComplexForLazy.hpp>

#include <SurfacePressureComplexForLazy.hpp>
#include <PostProcNonLinearElasticityStresseOnRefindedMesh.hpp>

#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>

#include <petsctime.h>

#include <FEMethod_ComplexConstArea.hpp>

#include <FaceSplittingTool.hpp>
#include <ConfigurationalFractureMechanics.hpp>

using namespace ObosleteUsersModules;

using namespace MoFEM;

//phisical_equation_volume eq_solid = hooke; /*stvenant_kirchhoff;*/

struct CrackFrontData {

  Range crackFrontEdgeNodes; 

  PetscErrorCode initCrackFrontData(FieldInterface& m_field) {
    PetscFunctionBegin;
    ErrorCode rval;
    PetscErrorCode ierr;
    Range crack_corners_edges;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_corners_edges,true); CHKERRQ(ierr);
    rval = m_field.get_moab().get_connectivity(crack_corners_edges,crackFrontEdgeNodes,true); CHKERR_PETSC(rval);
    PetscFunctionReturn(0);
  }

  PetscErrorCode setCrackFrontIndices(FEMethod *fe_ptr,string &material_field_name,vector<DofIdx>& GlobIndices,bool not_at_crack_front) {
    PetscFunctionBegin;
    if(!crackFrontEdgeNodes.empty()) {
    for(_IT_GET_FEROW_BY_NAME_DOFS_FOR_LOOP_(fe_ptr,material_field_name,dof)) {
      Range::iterator nit = find(crackFrontEdgeNodes.begin(),crackFrontEdgeNodes.end(),dof->get_ent());
      if(not_at_crack_front) {
	//if nit is not a part of crack front set
	if(nit != crackFrontEdgeNodes.end()) continue;
      } else {
	//if nit is part of crack front set
	if(nit == crackFrontEdgeNodes.end()) continue;
      }
      vector<DofIdx>::iterator it = find(GlobIndices.begin(),GlobIndices.end(),dof->get_petsc_gloabl_dof_idx());
      if(it != GlobIndices.end()) {
	*it = -1;
      }
    }}
    PetscFunctionReturn(0);
  }

};

struct MyNonLinearSpatialElasticFEMthod: public NonLinearSpatialElasticFEMthod,CrackFrontData {

  MyNonLinearSpatialElasticFEMthod(FieldInterface& _m_field,double _lambda,double _mu,int _verbose = 0): 
    FEMethod_ComplexForLazy_Data(_m_field,_verbose),
    NonLinearSpatialElasticFEMthod(_m_field,_lambda,_mu,0,_verbose) {}

  PetscErrorCode AssembleSpatialCoupledTangent(Mat B) {
    PetscFunctionBegin;
    vector<DofIdx> frontRowGlobMaterial = RowGlobMaterial[0];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial,true); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN:
	if(KHh.size1()!=frontRowGlobMaterial.size()) {
	  SETERRQ(PETSC_COMM_SELF,1,"KHh.size()!=frontRowGlobMaterial.size()");
	}
        ierr = MatSetValues(B,
	  frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	  ColGlobSpatial[i_nodes].size(),&*(ColGlobSpatial[i_nodes].begin()),
	  &*(KHh.data().begin()),ADD_VALUES); CHKERRQ(ierr);
        for(int ee = 0;ee<6;ee++) {
	  ierr = MatSetValues(B,
	    frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	    ColGlobSpatial[1+ee].size(),&*(ColGlobSpatial[1+ee].begin()),
	    &*(KHedge_data[ee].data().begin()),ADD_VALUES); CHKERRQ(ierr);
        }
        for(int ff = 0;ff<4;ff++) {
	  ierr = MatSetValues(B,
	    frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	    ColGlobSpatial[1+6+ff].size(),&*(ColGlobSpatial[1+6+ff].begin()),
	    &*(KHface_data[ff].data().begin()),ADD_VALUES); CHKERRQ(ierr);
        }
        ierr = MatSetValues(B,
	  frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	  ColGlobSpatial[i_volume].size(),&*(ColGlobSpatial[i_volume].begin()),
	  &*(KHvolume.data().begin()),ADD_VALUES); CHKERRQ(ierr);
        break;
      default:
        SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpComplexForLazyStart(); CHKERRQ(ierr);
    ierr = GetIndicesSpatial(); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION: { 
        ierr = GetFint(); CHKERRQ(ierr);
	ierr = AssembleSpatialFint(snes_f); CHKERRQ(ierr);
      }
      break;
      case CTX_SNESSETJACOBIAN:
	ierr = GetTangent(); CHKERRQ(ierr);
	ierr = AssembleSpatialTangent(snes_B); CHKERRQ(ierr);
	if(isCoupledProblem) {
	  ierr = GetIndicesRow(RowGlobMaterial,material_field_name); CHKERRQ(ierr);
	  ierr = AssembleSpatialCoupledTangent(snes_B); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};


struct MyEshelbyFEMethod: public EshelbyFEMethod,CrackFrontData {

  MyEshelbyFEMethod(FieldInterface& _m_field,double _lambda,double _mu,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_m_field,_verbose), 
    EshelbyFEMethod(_m_field,_lambda,_mu,_verbose) {
    type_of_analysis = material_analysis;
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch (ts_ctx) {
      case CTX_TSSETIFUNCTION: {
	snes_ctx = CTX_SNESSETFUNCTION;
	snes_f = ts_F;
	break;
      }
      case CTX_TSSETIJACOBIAN: {
	snes_ctx = CTX_SNESSETJACOBIAN;
	snes_B = ts_B;
	break;
      }
      default:
      break;
    }
    ierr = EshelbyFEMethod::preProcess(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  PetscErrorCode AssembleMaterialTangent(Mat B) {
    PetscFunctionBegin;

    vector<DofIdx> frontRowGlobMaterial = RowGlobMaterial[0];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial,true); CHKERRQ(ierr);

    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION:
      case CTX_SNESSETJACOBIAN:
	ierr = MatSetValues(B,
	  frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	  ColGlobMaterial[0].size(),&*(ColGlobMaterial[0].begin()),
	  &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode AssembleMaterialFint(Vec f) {
    PetscFunctionBegin;
    vector<DofIdx> frontRowGlobMaterial = RowGlobMaterial[0];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial,true); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESNONE:
      case CTX_SNESSETFUNCTION: { 
	ierr = VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
	ierr = VecSetValues(f,frontRowGlobMaterial.size(),&(frontRowGlobMaterial[0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct MyMeshSmoothingFEMethod: public MeshSmoothingFEMethod,CrackFrontData {

  Vec frontF;
  Vec tangentFrontF;

  MyMeshSmoothingFEMethod(FieldInterface& _m_field,int _verbose = 0):
    FEMethod_ComplexForLazy_Data(_m_field,_verbose), 
    MeshSmoothingFEMethod(_m_field,_verbose),
      frontF(PETSC_NULL),tangentFrontF(PETSC_NULL),
      stabilise(false) {
      type_of_analysis = mesh_quality_analysis;

      g_NTET.resize(4*1);
      ShapeMBTET(&g_NTET[0],G_TET_X1,G_TET_Y1,G_TET_Z1,1);
      g_TET_W = G_TET_W1;

    }

  ~MyMeshSmoothingFEMethod() {
    if(frontF!=PETSC_NULL) {
      ierr = VecDestroy(&frontF); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      frontF = PETSC_NULL;
    }
    if(tangentFrontF!=PETSC_NULL) {
      ierr = VecDestroy(&tangentFrontF); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      tangentFrontF = PETSC_NULL;
    }
  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;

    switch (ts_ctx) {
      case CTX_TSSETIFUNCTION: {
	snes_ctx = CTX_SNESSETFUNCTION;
	snes_f = ts_F;
	break;
      }
      case CTX_TSSETIJACOBIAN: {
	snes_ctx = CTX_SNESSETJACOBIAN;
	snes_B = ts_B;
	break;
      }
      default:
      break;
    }

    if(crackFrontEdgeNodes.size()>0) {
      if(frontF == PETSC_NULL) {
	ierr = VecDuplicate(snes_f,&frontF); CHKERRQ(ierr);
	ierr = VecSetOption(frontF,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
	ierr = VecDuplicate(snes_f,&tangentFrontF); CHKERRQ(ierr);
	ierr = VecSetOption(tangentFrontF,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE); CHKERRQ(ierr);
      }
      ierr = VecZeroEntries(frontF); CHKERRQ(ierr);
      ierr = VecGhostUpdateBegin(frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      ierr = VecGhostUpdateEnd(frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	if(!crackFrontEdgeNodes.empty()) {
	  ierr = VecAssemblyBegin(frontF); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(frontF); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(frontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(frontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecGhostUpdateBegin(frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	  ierr = VecGhostUpdateEnd(frontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	}
	break;
	default:
	break;
      }
    }
    PetscFunctionReturn(0);
  }

  bool stabilise;
  PetscErrorCode AssembleMeshSmoothingTangent(Mat B) {
    PetscFunctionBegin;
    vector<DofIdx> frontRowGlobMaterial = RowGlobMaterial[i_nodes];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial,false); CHKERRQ(ierr);
    vector<DofIdx> frontRowGlobMaterial_front_only = RowGlobMaterial[i_nodes];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial_front_only,true); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESSETJACOBIAN:
	ierr = MatSetValues(B,
	  frontRowGlobMaterial.size(),&*(frontRowGlobMaterial.begin()),
	  ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	  &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	if(stabilise) {
	  ierr = MatSetValues(B,
	    frontRowGlobMaterial_front_only.size(),&*(frontRowGlobMaterial_front_only.begin()),
	    ColGlobMaterial[i_nodes].size(),&*(ColGlobMaterial[i_nodes].begin()),
	    &*(KHH.data().begin()),ADD_VALUES); CHKERRQ(ierr);
	}  
	if(!crackFrontEdgeNodes.empty()) {
	  double *f_tangent_front_mesh_array;
	  if(tangentFrontF==PETSC_NULL) SETERRQ(PETSC_COMM_SELF,1,"vector for crack front not created");
	  ierr = VecGetArray(tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);
	  for(int nn = 0;nn<4;nn++) {
	    FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator dit,hi_dit;
	    dit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple("LAMBDA_CRACK_TANGENT_CONSTRAIN",conn[nn]));
	    hi_dit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple("LAMBDA_CRACK_TANGENT_CONSTRAIN",conn[nn]));
	    if(distance(dit,hi_dit)>0) {
	      FENumeredDofMoFEMEntity_multiIndex::index<Composite_Name_And_Ent_mi_tag>::type::iterator diit,hi_diit;
	      diit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().lower_bound(boost::make_tuple(material_field_name,conn[nn]));
	      hi_diit = rowPtr->get<Composite_Name_And_Ent_mi_tag>().upper_bound(boost::make_tuple(material_field_name,conn[nn]));
	      for(;diit!=hi_diit;diit++) {
		for(unsigned int ddd = 0;ddd<ColGlobMaterial[i_nodes].size();ddd++) {
		  if(frontRowGlobMaterial_front_only[3*nn+diit->get_dof_rank()]!=diit->get_petsc_gloabl_dof_idx()) {
		    SETERRQ2(PETSC_COMM_SELF,1,"data inconsistency %d != %d",
		      3*nn+diit->get_dof_rank(),diit->get_petsc_gloabl_dof_idx());
		  }
		  if(diit->get_petsc_local_dof_idx()==-1) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
		  double g = f_tangent_front_mesh_array[diit->get_petsc_local_dof_idx()]*KHH(3*nn+diit->get_dof_rank(),ddd);
		  DofIdx lambda_idx = dit->get_petsc_gloabl_dof_idx();
		  ierr = MatSetValues(B,1,&lambda_idx,1,&ColGlobMaterial[i_nodes][ddd],&g,ADD_VALUES); CHKERRQ(ierr);
		}
	      }
	    }
	  }
	  ierr = VecRestoreArray(tangentFrontF,&f_tangent_front_mesh_array); CHKERRQ(ierr);
	}
	break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode AssembleMeshSmoothingFint(Vec f) {
    PetscFunctionBegin;
    vector<DofIdx> frontRowGlobMaterial = RowGlobMaterial[i_nodes];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial,false); CHKERRQ(ierr);
    vector<DofIdx> frontRowGlobMaterial_front_only = RowGlobMaterial[i_nodes];
    ierr = setCrackFrontIndices(this,material_field_name,frontRowGlobMaterial_front_only,true); CHKERRQ(ierr);
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecSetOption(f,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
	//cerr << "Fint_h " << Fint_h << endl;
	ierr = VecSetValues(f,frontRowGlobMaterial.size(),&(frontRowGlobMaterial[0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	if(stabilise) {
	  ierr = VecSetValues(
	    f,frontRowGlobMaterial_front_only.size(),&(frontRowGlobMaterial_front_only[0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	}  
	if(!crackFrontEdgeNodes.empty()) {
	  if(frontF==PETSC_NULL) SETERRQ(PETSC_COMM_SELF,1,"vector for crack front not created");
	  ierr = VecSetValues(frontF,frontRowGlobMaterial_front_only.size(),&(frontRowGlobMaterial_front_only[0]),&(Fint_H.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	}
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

};

struct TangentWithMeshSmoothingFrontConstrain_FEMethod: public C_CONSTANT_AREA_FEMethod {

  MyMeshSmoothingFEMethod *meshFEPtr;
  const double eps;

  TangentWithMeshSmoothingFrontConstrain_FEMethod(FieldInterface& _mField,
    MyMeshSmoothingFEMethod *_mesh_fe_ptr,
    string _lambda_field_name,int _verbose = 0):
    C_CONSTANT_AREA_FEMethod(_mField,PETSC_NULL,PETSC_NULL,_lambda_field_name,_verbose),
    meshFEPtr(_mesh_fe_ptr),eps(1e-10) {}

  Tag thFrontTangent;
  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    switch (ts_ctx) {
      case CTX_TSSETIFUNCTION: {
	snes_ctx = CTX_SNESSETFUNCTION;
	snes_f = ts_F;
	break;
      }
      case CTX_TSSETIJACOBIAN: {
	snes_ctx = CTX_SNESSETJACOBIAN;
	snes_B = ts_B;
	break;
      }
      default:
      break;
    }

    ierr = C_CONSTANT_AREA_FEMethod::preProcess(); CHKERRQ(ierr);
    /*//TAG  - only for one proc analysis
    double def[] = {0,0,0};
    rval = mField.get_moab().tag_get_handle("FRONT_TANGENT",3,MB_TYPE_DOUBLE,
      thFrontTangent,MB_TAG_CREAT|MB_TAG_SPARSE,&def); CHKERR_THROW(rval);*/
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = VecZeroEntries(meshFEPtr->tangentFrontF); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(meshFEPtr->tangentFrontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(meshFEPtr->tangentFrontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(meshFEPtr->tangentFrontF); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(meshFEPtr->tangentFrontF); CHKERRQ(ierr);
	/*//resent tags - only for one proc analysis
	for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,lambda_field_name,dit)) {
	  EntityHandle ent = dit->get_ent();
	  rval = mField.get_moab().tag_set_data(thFrontTangent,&ent,1,def); CHKERR_PETSC(rval);
	}*/
      } break;
      case CTX_SNESSETJACOBIAN: {
	ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      } break;
      default: {
      } break;
    }
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    EntityHandle face = fePtr->get_ent();
    Range tet;
    rval = mField.get_moab().get_adjacencies(&face,1,3,false,tet); CHKERR_PETSC(rval);
    try {
      ierr = getData(true,true); CHKERRQ(ierr);
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    }
    try {
      ublas::vector<double,ublas::bounded_array<double,9> > ELEM_CONSTRAIN1(9);
      ierr = calcDirevatives(
	&*diffNTRI.data().begin(),&*dofs_X.data().begin(),
	NULL,NULL,NULL,&*ELEM_CONSTRAIN1.data().begin(),NULL); CHKERRQ(ierr);
      //take in account face orientation in respect crack surface
      Tag th_interface_side;
      rval = moab.tag_get_handle("INTERFACE_SIDE",th_interface_side); CHKERR_PETSC(rval);
      int side;
      rval = moab.tag_get_data(th_interface_side,&face,1,&side); CHKERR_PETSC(rval);
      if(side == 1) {
	ELEM_CONSTRAIN1 *= -1;
      }
      //quality 
      ublas::vector<double,ublas::bounded_array<double,9> > F_FRONT_MESH_SMOOTHING(9);
      ublas::noalias(F_FRONT_MESH_SMOOTHING) = ublas::zero_vector<double>(9);
      double *f_front_mesh_array;
      ierr = VecGetArray(meshFEPtr->frontF,&f_front_mesh_array); CHKERRQ(ierr);
      for(int nn = 0;nn<3;nn++) {
	for(int dd = 0;dd<3;dd++) {
	  if(local_disp_dofs_row_idx[3*nn+dd]==-1) continue;
	  F_FRONT_MESH_SMOOTHING[3*nn+dd] = f_front_mesh_array[local_disp_dofs_row_idx[3*nn+dd]];
	}
      }
      ierr = VecRestoreArray(meshFEPtr->frontF,&f_front_mesh_array); CHKERRQ(ierr);
      //tangent
      if(snes_ctx == CTX_SNESSETJACOBIAN) {
	double center[3]; 
	tricircumcenter3d_tp(&coords.data()[0],&coords.data()[3],&coords.data()[6],center,NULL,NULL);
	cblas_daxpy(3,-1,&coords.data()[0],1,center,1);
	double r = cblas_dnrm2(3,center,1);
	for(int nn = 0;nn<3;nn++) {
	  for(int dd = 0;dd<3;dd++) {
	    // ---> calculate tangent starts here
	    ublas::vector<double,ublas::bounded_array<double,9> > idofs_X(9,0);
	    idofs_X[nn*3+dd] = r*eps;
	    ublas::vector<double,ublas::bounded_array<double,9> > dELEM_CONSTRAIN1(9);
	    ierr = calcDirevatives(
	      &*diffNTRI.data().begin(),
	      &*dofs_X.data().begin(),
	      &*idofs_X.data().begin(),
	      NULL,NULL,NULL,&*dELEM_CONSTRAIN1.data().begin()); CHKERRQ(ierr);
	    if(side == 1) {
	      dELEM_CONSTRAIN1 /= -r*eps;
	    } else {
	      dELEM_CONSTRAIN1 /= +r*eps;
	    }
	    //dg -> C*q_quality
	    double g[3] = {0,0,0};
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      for(int ddd = 0;ddd<3;ddd++) {
		g[nnn] += dELEM_CONSTRAIN1[nnn*3+ddd]*F_FRONT_MESH_SMOOTHING[3*nnn+ddd];
	      }
	    }
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      ierr = MatSetValues(snes_B,
		1,&lambda_dofs_row_indx[nnn],1,&disp_dofs_col_idx[3*nn+dd],
		&g[nnn],ADD_VALUES); CHKERRQ(ierr);
	    }
	    //dCT_lambda
	    for(int nnn = 0;nnn<3;nnn++) {
	      for(int ddd = 0;ddd<3;ddd++) {
		dELEM_CONSTRAIN1[nnn*3+ddd] *= lambda[nnn];
	      }
	    }
	    for(int nnn = 0;nnn<3;nnn++) {
	      if(lambda_dofs_row_indx[nnn] == -1) continue;
	      ierr = MatSetValues(snes_B,
		3,&disp_dofs_row_idx[3*nnn],
		1,&disp_dofs_col_idx[3*nn+dd],
		&dELEM_CONSTRAIN1[3*nnn],ADD_VALUES); CHKERRQ(ierr);
	    }
	    // ---> calculate tangent end here
	  }
	}
      }
      switch(snes_ctx) {
	case CTX_SNESSETFUNCTION: { 
	  ublas::vector<double,ublas::bounded_array<double,3> > g(3);
	  for(int nn = 0;nn<3;nn++) {
	    g[nn] = cblas_ddot(3,&ELEM_CONSTRAIN1[3*nn],1,&F_FRONT_MESH_SMOOTHING[3*nn],1);
	  }
	  //cerr << "g : " << g << endl;
	  ierr = VecSetValues(snes_f,3,&*lambda_dofs_row_indx.data().begin(),&*g.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  ierr = VecSetValues(meshFEPtr->tangentFrontF,9,&disp_dofs_row_idx[0],&*ELEM_CONSTRAIN1.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  ublas::vector<double,ublas::bounded_array<double,9> > f(9);
	  for(int nn = 0;nn<3;nn++) {
	    for(int dd = 0;dd<3;dd++) {
	      f[nn*3+dd] = lambda[nn]*ELEM_CONSTRAIN1[3*nn+dd];
	    }
	  }
	  //cerr << "f : " << f << endl;
	  ierr = VecSetValues(snes_f,9,&disp_dofs_row_idx[0],&*f.data().begin(),ADD_VALUES); CHKERRQ(ierr);
	  /*//TAG - only for one proc analysis
	  for(int nn = 0;nn<3;nn++) { 
	    EntityHandle ent = lambda_dofs_row_ents[nn];
	    if(ent == no_handle) continue;
	    double *t;
	    rval = mField.get_moab().tag_get_by_ptr(thFrontTangent,&ent,1,(const void **)&t); CHKERR_PETSC(rval);
	    cblas_daxpy(3,+1,&ELEM_CONSTRAIN1[3*nn],1,t,1);
	  }*/
	} break;
	case CTX_SNESSETJACOBIAN: {
	  for(int nn = 0;nn<3;nn++) {
	    int lambda_dof_idx = lambda_dofs_col_indx[nn];
	    ierr = MatSetValues(snes_B,3,&disp_dofs_row_idx[3*nn],1,&lambda_dof_idx,&ELEM_CONSTRAIN1[3*nn],ADD_VALUES); CHKERRQ(ierr);
	  }
	} break;
	default:
	  break;
      }
    } catch (const std::exception& ex) {
      ostringstream ss;
      ss << "thorw in method: " << ex.what() << endl;
      SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
    } 
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(meshFEPtr->tangentFrontF); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(meshFEPtr->tangentFrontF); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(meshFEPtr->tangentFrontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(meshFEPtr->tangentFrontF,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(meshFEPtr->tangentFrontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(meshFEPtr->tangentFrontF,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
      } break;
      case CTX_SNESSETJACOBIAN: {
	ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      } break;
      default: {
      } break;
    }
    PetscFunctionReturn(0);
  }
  
};

struct BothSurfaceConstrains: public FEMethod {

  FieldInterface& mField;
  BothSurfaceConstrains(FieldInterface& m_field): mField(m_field) {} 

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  snes_ctx = CTX_SNESSETFUNCTION;
	  snes_f = ts_F;
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  snes_ctx = CTX_SNESSETJACOBIAN;
	  snes_B = ts_B;
	  break;
	}
	default:
	break;
    }
    switch(snes_ctx) {
	case CTX_SNESSETFUNCTION: { 
	  ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	}
	break;
	case CTX_SNESSETJACOBIAN: 
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	break;
	default:
	break;
    }
    PetscFunctionReturn(0);
  }
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    //cerr << "AAAAAAAA\n";
    PetscErrorCode ierr;
    vector<int> lambda_dofs(9,-1);
    vector<double> lambda_vals(9,0);
    for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"LAMBDA_BOTH_SIDES",MBVERTEX,it)) {
	lambda_dofs[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] = it->get_petsc_gloabl_dof_idx();
	lambda_vals[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] = it->get_FieldData();
	//cerr << "l " 
	  //<< 3*it->side_number_ptr->side_number+it->get_EntDofIdx() << " " 
	  //<< lambda_dofs[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] << " " 
	  //<< lambda_vals[3*it->side_number_ptr->side_number+it->get_EntDofIdx()] 
	  //<< endl;
    }
    vector<int> positions_dofs(18,-1);
    vector<double> positions_vals(18,0);
    for(_IT_GET_FEROW_BY_TYPE_DOFS_FOR_LOOP_(this,"MESH_NODE_POSITIONS",MBVERTEX,it)) {
	int dd = 3*it->side_number_ptr->side_number+it->get_EntDofIdx();
	if(lambda_dofs[dd>8 ? dd -9 : dd] == -1) continue;
	positions_dofs[dd] = it->get_petsc_gloabl_dof_idx();
	positions_vals[dd] = it->get_FieldData();
	//cerr << "p " << dd << " " << positions_dofs[dd] << " " << positions_vals[dd] << endl;
    }
    //cerr << endl;
    const double alpha = 0;
    const double betha = 1e2;
    switch(snes_ctx) {
	case CTX_SNESSETFUNCTION:  
	  for(int ii = 0;ii<9;ii++) {
	    if(lambda_dofs[ii] == -1) continue;
	    double val1 = betha*(positions_vals[0+ii] - positions_vals[9+ii]) + alpha*lambda_vals[ii];
	    ierr = VecSetValue(snes_f,lambda_dofs[ii],val1,INSERT_VALUES); CHKERRQ(ierr);
	    double val2 = betha*lambda_vals[ii];
	    ierr = VecSetValue(snes_f,positions_dofs[0+ii],+val2,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = VecSetValue(snes_f,positions_dofs[9+ii],-val2,INSERT_VALUES); CHKERRQ(ierr);
	  }
	break;
	case CTX_SNESSETJACOBIAN: 
	  for(int ii = 0;ii<9;ii++) {
	    if(lambda_dofs[ii] == -1) continue;
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],positions_dofs[0+ii],+1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],positions_dofs[9+ii],-1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,positions_dofs[0+ii],lambda_dofs[ii],+1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,positions_dofs[9+ii],lambda_dofs[ii],-1*betha,INSERT_VALUES); CHKERRQ(ierr);
	    ierr = MatSetValue(snes_B,lambda_dofs[ii],lambda_dofs[ii],alpha,INSERT_VALUES); CHKERRQ(ierr);
	  }
	break;
	default:
	break;
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
	case CTX_SNESSETJACOBIAN: 
	  ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	  ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	break;
	default:
	break;
    }
    PetscFunctionReturn(0);
  }

};

 
PetscErrorCode ConfigurationalFractureMechanics::set_material_fire_wall(FieldInterface& m_field) {
  PetscFunctionBegin;
  ErrorCode rval;
  BitRefLevel def_bit_level = 0;
  rval = m_field.get_moab().tag_get_handle("_Materiar_FireWall",sizeof(Material_FirelWall_def),MB_TYPE_OPAQUE,
    th_MaterialFireWall,MB_TAG_CREAT|MB_TAG_SPARSE|MB_TAG_BYTES,&def_bit_level); 
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  rval = m_field.get_moab().tag_get_by_ptr(th_MaterialFireWall,&root_meshset,1,(const void**)&material_FirelWall); CHKERR_PETSC(rval);
  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::thermal_field(FieldInterface& m_field) {
  PetscFunctionBegin;
  if(material_FirelWall->operator[](FW_thermal_field)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_thermal_field);

  PetscErrorCode ierr;

  ierr = m_field.add_field("TEMPERATURE",H1,1,MF_ZERO); CHKERRQ(ierr);

  Range level_tets;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"TEMPERATURE"); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::spatial_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_spatial_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_spatial_problem_definition);

  PetscErrorCode ierr;
  ErrorCode rval;

  //Fields
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("SPATIAL_DISPLACEMENT",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("ELASTIC",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","TEMPERATURE"); CHKERRQ(ierr);
  }

  //add neumman finite elemnets to add static boundary conditions
  ierr = m_field.add_finite_element("NEUAMNN_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUAMNN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  Range level_tris;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    Range tris;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    tris = intersect(level_tris,tris);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    Range tris;
    rval = m_field.get_moab().get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    tris = intersect(level_tris,tris);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);

  //add nodal force element
  ierr = MetaNodalForces::addNodalForceElement(m_field,"ELASTIC_MECHANICS","SPATIAL_POSITION"); CHKERRQ(ierr);

  Range level_tets;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(level_tets,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //ierr = m_field.add_ents_to_field_by_TETs(level_tets,"MESH_NODE_DISPLACEMENT"); CHKERRQ(ierr);
  Range blocked_tets;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBTET,blocked_tets); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(blocked_tets,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //ierr = m_field.add_ents_to_field_by_TETs(blocked_tets,"MESH_NODE_DISPLACEMENT"); CHKERRQ(ierr);

  PetscInt order;
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    Tag th_set_order;
    rval = m_field.get_moab().tag_get_handle("_SET_ORDER",th_set_order); CHKERR_PETSC(rval);
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    rval = m_field.get_moab().tag_get_data(th_set_order,&root_meshset,1,&order); CHKERR_PETSC(rval);
  }

  //set app. order
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_DISPLACEMENT",1); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  //ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_DISPLACEMENT",1); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::material_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_material_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_material_problem_definition);

  PetscErrorCode ierr;

  //Fields
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("MATERIAL",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data 
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  //
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_off_field_row("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("MATERIAL","TEMPERATURE"); CHKERRQ(ierr);
  }

  //define problems
  ierr = m_field.add_problem("MATERIAL_MECHANICS",MF_ZERO); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::coupled_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_coupled_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_coupled_problem_definition);

  PetscErrorCode ierr;

  //Fields
  //ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("MATERIAL_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  //FE
  ierr = m_field.add_finite_element("ELASTIC_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("MATERIAL_COUPLED",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("MESH_SMOOTHER",MF_ZERO); CHKERRQ(ierr);

  //fes definitions
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

  //
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","SPATIAL_POSITION"); CHKERRQ(ierr);
  if(m_field.check_field("TEMPERATURE")) {
    ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","TEMPERATURE"); CHKERRQ(ierr);
  }

  //
  ierr = m_field.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("COUPLED_PROBLEM",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","ELASTIC_COUPLED"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","NEUAMNN_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","MESH_SMOOTHER"); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","BOTH_SIDE_OF_CRACK"); CHKERRQ(ierr);

  bool cs = true;
  if(cs) {
    ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM",ss2.str()); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::arclength_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;

  ErrorCode rval;

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_arc_lenhghat_definition);

  PetscErrorCode ierr;
  
  ierr = m_field.add_field("LAMBDA",NOFIELD,1); CHKERRQ(ierr);

  ierr = m_field.add_finite_element("ARC_LENGTH"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);
  //elem data
  ierr = m_field.modify_finite_element_add_field_data("ARC_LENGTH","LAMBDA"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MATERIAL_COUPLED","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("MESH_SMOOTHER","LAMBDA"); CHKERRQ(ierr);

  ierr = m_field.modify_problem_add_finite_element("COUPLED_PROBLEM","ARC_LENGTH"); CHKERRQ(ierr);

  //this entity will carray data for this finite element
  EntityHandle meshset_FE_ARC_LENGTH;
  rval = m_field.get_moab().create_meshset(MESHSET_SET,meshset_FE_ARC_LENGTH); CHKERR_PETSC(rval);
  //get LAMBDA field meshset
  EntityHandle meshset_field_LAMBDA = m_field.get_field_meshset("LAMBDA");
  //add LAMBDA field meshset to finite element ARC_LENGTH
  rval = m_field.get_moab().add_entities(meshset_FE_ARC_LENGTH,&meshset_field_LAMBDA,1); CHKERR_PETSC(rval);
  //add finite element ARC_LENGTH meshset to refinment database (all ref bit leveles)
  ierr = m_field.seed_ref_level_MESHSET(meshset_FE_ARC_LENGTH,BitRefLevel().set()); CHKERRQ(ierr);
  //finally add created meshset to the ARC_LENGTH finite element
  ierr = m_field.add_ents_to_finite_element_by_MESHSET(meshset_FE_ARC_LENGTH,"ARC_LENGTH"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::constrains_problem_definition(FieldInterface& m_field) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_constrains_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_constrains_problem_definition);

  ErrorCode rval;
  PetscErrorCode ierr;

  bool cs = true;

  //Fields
  ierr = m_field.add_field("LAMBDA_SURFACE",H1,1,MF_ZERO); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ierr = m_field.add_field(ss.str(),H1,1,MF_ZERO); CHKERRQ(ierr);
  }
  //CRACK
  if(cs) {
    ierr = m_field.add_field("LAMBDA_CRACK_SURFACE",H1,1,MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.add_field("LAMBDA_BOTH_SIDES",H1,3,MF_ZERO); CHKERRQ(ierr);
  }

  //FE
  ierr = m_field.add_finite_element("C_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("CTC_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("CandCT_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.add_finite_element(ss0.str(),MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.add_finite_element(ss1.str(),MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.add_finite_element(ss2.str(),MF_ZERO); CHKERRQ(ierr);
  }

  //CRACK
  if(cs) {
    //ierr = m_field.add_finite_element("C_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    //ierr = m_field.add_finite_element("CTC_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("CandCT_CRACK_SURFACE_ELEM",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.add_finite_element("BOTH_SIDE_OF_CRACK",MF_ZERO); CHKERRQ(ierr);
  }

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("CandCT_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("CandCT_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_finite_element_add_field_row(ss0.str(),ss.str()); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(ss0.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(ss0.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(ss0.str(),ss.str()); CHKERRQ(ierr);
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_finite_element_add_field_row(ss1.str(),ss.str()); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(ss1.str(),ss.str()); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(ss1.str(),ss.str()); CHKERRQ(ierr);
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_finite_element_add_field_row(ss2.str(),ss.str()); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col(ss2.str(),ss.str()); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(ss2.str(),"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data(ss2.str(),ss.str()); CHKERRQ(ierr);
  }

  //CRACK
  if(cs) {
    /*ierr = m_field.modify_finite_element_add_field_row("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("C_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = m_field.modify_finite_element_add_field_row("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("CTC_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);*/

    ierr = m_field.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("CandCT_CRACK_SURFACE_ELEM","LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);

    ierr = m_field.modify_finite_element_add_field_row("BOTH_SIDE_OF_CRACK","LAMBDA_BOTH_SIDES"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("BOTH_SIDE_OF_CRACK","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("BOTH_SIDE_OF_CRACK","LAMBDA_BOTH_SIDES"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("BOTH_SIDE_OF_CRACK","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("BOTH_SIDE_OF_CRACK","LAMBDA_BOTH_SIDES"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("BOTH_SIDE_OF_CRACK","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  }

  //define problems
  ierr = m_field.add_problem("CCT_ALL_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_problem("C_ALL_MATRIX",MF_ZERO); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ierr = m_field.modify_problem_add_finite_element("C_ALL_MATRIX",ss0.str()); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("CCT_ALL_MATRIX",ss1.str()); CHKERRQ(ierr);
  }

  //CRACK //this for testing only
  /*if(cs) {
    ierr = m_field.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }*/

  Range level_tets,level_tris,level_edges,level_nodes;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

  Range corners_edges,corners_nodes,surfaces_faces;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  corners_edges = intersect(corners_edges,level_edges);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  corners_nodes = intersect(corners_nodes,level_nodes);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces_faces,true); CHKERRQ(ierr);
  Skinner skin(&m_field.get_moab());
  Range skin_faces; 
  rval = skin.find_skin(0,level_tets,false,skin_faces); CHKERR(rval);
  surfaces_faces = intersect(surfaces_faces,skin_faces);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 100 = %d\n",corners_edges.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of NODESET 101 = %d\n",corners_nodes.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 102 = %d\n",surfaces_faces.size()); CHKERRQ(ierr);

  Range crack_surfaces_faces,crack_front_edges;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surfaces_faces,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  crack_surfaces_faces = intersect(crack_surfaces_faces,level_tris);
  crack_front_edges = intersect(crack_front_edges,level_edges);

  Range crack_front_nodes;
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  Range corners_edges_nodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());

  //FE Couple nodes on two sides of crack
  Range one_side_crack_surface_faces;
  vector<int> sides(crack_surfaces_faces.size());
  Tag th_interface_side;
  const int def_side[] = {0};
  rval = m_field.get_moab().tag_get_handle("INTERFACE_SIDE",1,MB_TYPE_INTEGER,
    th_interface_side,MB_TAG_CREAT|MB_TAG_SPARSE,def_side); CHKERR_PETSC(rval);
  rval = m_field.get_moab().tag_get_data(th_interface_side,crack_surfaces_faces,&*sides.begin()); CHKERR_PETSC(rval);
  int nn = 0;
  for(Range::iterator fit = crack_surfaces_faces.begin();fit!=crack_surfaces_faces.end();fit++,nn++) {
    if(sides[nn] == 0) {
      one_side_crack_surface_faces.insert(*fit);
    } 
  }
  Range one_side_crack_surface_faces_nodes;
  rval = m_field.get_moab().get_connectivity(one_side_crack_surface_faces,one_side_crack_surface_faces_nodes,true); CHKERR_PETSC(rval);
  one_side_crack_surface_faces_nodes = subtract(one_side_crack_surface_faces_nodes,corners_nodes);
  Range one_side_crack_surface_faces_nodes_without_crack_front;
  one_side_crack_surface_faces_nodes_without_crack_front = subtract(one_side_crack_surface_faces_nodes,crack_front_nodes);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Surfaces Faces = %d\n",crack_surfaces_faces.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Surfaces Faces On One Side = %d\n",one_side_crack_surface_faces.size()); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Front Edges = %d\n",crack_front_edges.size()); CHKERRQ(ierr);

  Interface& moab = m_field.get_moab();

  ierr = m_field.seed_finite_elements(surfaces_faces); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces,"C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces,"CandCT_SURFACE_ELEM"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << msId;
    ostringstream ss1;
    ss1 << "CTC_SURFACE_ELEM_msId_" << msId;
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << msId;
    Range surfaces_faces_msId;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(msId,SIDESET,2,surfaces_faces_msId,true); CHKERRQ(ierr);
    surfaces_faces_msId = intersect(surfaces_faces_msId,skin_faces);
    ierr = m_field.seed_finite_elements(surfaces_faces_msId); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET ( mdId =  %d ) = %d\n",msId,surfaces_faces_msId.size()); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss0.str()); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss1.str()); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(surfaces_faces_msId,ss2.str()); CHKERRQ(ierr);
  }

  if(cs) {
    ierr = m_field.seed_finite_elements(crack_surfaces_faces); CHKERRQ(ierr);
    //ierr = m_field.add_ents_to_finite_element_by_TRIs(one_side_crack_surface_faces,"C_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    //ierr = m_field.add_ents_to_finite_element_by_TRIs(one_side_crack_surface_faces,"CTC_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(one_side_crack_surface_faces,"CandCT_CRACK_SURFACE_ELEM"); CHKERRQ(ierr);
  }

  Range surfaces_nodes;
  rval = moab.get_connectivity(surfaces_faces,surfaces_nodes,true); CHKERR_PETSC(rval);
  surfaces_nodes = subtract(surfaces_nodes,corners_nodes);
  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_VERTICEs(surfaces_nodes,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  //NOTE: always order should be 1
  ierr = m_field.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    Range surfaces_faces_msId;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(msId,SIDESET,2,surfaces_faces_msId,true); CHKERRQ(ierr);
    surfaces_faces_msId = intersect(surfaces_faces_msId,level_tris);
    Range surfaces_nodes_msId;
    rval = m_field.get_moab().get_connectivity(surfaces_faces_msId,surfaces_nodes_msId,true); CHKERR_PETSC(rval);
    surfaces_nodes_msId = subtract(surfaces_nodes_msId,corners_nodes);
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    ierr = m_field.add_ents_to_field_by_VERTICEs(surfaces_nodes_msId,ss.str()); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,ss.str(),1); CHKERRQ(ierr);
  }

  //CRCAK
  if(cs) {

    ierr = m_field.remove_ents_from_field("LAMBDA_SURFACE",one_side_crack_surface_faces_nodes_without_crack_front); CHKERRQ(ierr);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      int msId = it->get_msId();
      if((msId < 10200)||(msId >= 10300)) continue;
      ostringstream ss;
      ss << "LAMBDA_SURFACE_msId_" << msId;
      ierr = m_field.remove_ents_from_field(ss.str(),one_side_crack_surface_faces_nodes_without_crack_front); CHKERRQ(ierr);
    }

    ierr = m_field.add_ents_to_field_by_VERTICEs(one_side_crack_surface_faces_nodes_without_crack_front,"LAMBDA_CRACK_SURFACE"); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_SURFACE",1); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_VERTICEs(one_side_crack_surface_faces_nodes_without_crack_front,"LAMBDA_BOTH_SIDES"); CHKERRQ(ierr);
    ierr = m_field.set_field_order(0,MBVERTEX,"LAMBDA_BOTH_SIDES",1); CHKERRQ(ierr);

    Tag th_my_ref_level;
    rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); 
    const EntityHandle root_meshset = m_field.get_moab().get_root_set();
    BitRefLevel *ptr_bit_level0;
    rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
    BitRefLevel& bit_level0 = *ptr_bit_level0;

    Range prisms;
    ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBPRISM,prisms); CHKERRQ(ierr);
    ierr = m_field.seed_ref_level(prisms,bit_level0); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_PRISMs(prisms,"BOTH_SIDE_OF_CRACK"); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_crack_front_problem_definition(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_constrains_crack_front_problem_definition)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_constrains_crack_front_problem_definition);

  PetscErrorCode ierr;
  ErrorCode rval;

  Range level_tets,level_tris,level_edges,level_nodes;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTET,level_tets); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);

  Range body_surface_nodes;
  {
    Range surfaces;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,surfaces,true); CHKERRQ(ierr);
    surfaces = intersect(surfaces,level_tris);
    for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
      int msId = it->get_msId();
      if((msId < 10200)||(msId >= 10300)) continue;
      Range surf;
      rval = m_field.get_moab().get_entities_by_type(it->get_meshset(),MBTRI,surf,true); CHKERR_PETSC(rval);
      surf = intersect(surf,level_tris);
      surfaces.merge(surf);
    }
    rval = m_field.get_moab().get_connectivity(surfaces,body_surface_nodes,true); CHKERR_PETSC(rval);
  }


  //Fields
  ierr = m_field.add_field("LAMBDA_CRACKFRONT_AREA",H1,1,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("LAMBDA_CRACK_TANGENT_CONSTRAIN",H1,1,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("GRIFFITH_FORCE",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_field("GRIFFITH_FORCE_TANGENT",H1,3,MF_ZERO); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("C_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("CTC_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("dCT_CRACKFRONT_AREA_ELEM",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("C_TANGENT_ELEM",MF_ZERO); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("C_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("C_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("CTC_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("CTC_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("dCT_CRACKFRONT_AREA_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("dCT_CRACKFRONT_AREA_ELEM","LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);

  ierr = m_field.modify_finite_element_add_field_row("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_TANGENT_ELEM","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("C_TANGENT_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("C_CRACKFRONT_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_problem("CTC_CRACKFRONT_MATRIX",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("CTC_CRACKFRONT_MATRIX","CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  {
    Range crack_surfaces_faces,crack_front_edges;
    ierr = m_field.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,crack_surfaces_faces,true); CHKERRQ(ierr);
    ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
    crack_surfaces_faces = intersect(crack_surfaces_faces,level_tris);
    crack_front_edges = intersect(crack_front_edges,level_edges);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Faces = %d\n",crack_surfaces_faces.size()); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Crack Front = %d\n",crack_front_edges.size()); CHKERRQ(ierr);
    Range crack_front_nodes;
    rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
    crack_front_nodes = intersect(crack_front_nodes,level_nodes);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Nodes = %d\n",crack_front_nodes.size()); CHKERRQ(ierr);
    Range tangent_crack_front_nodes = subtract(crack_front_nodes,body_surface_nodes);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Tangent Crack Front Nodes = %d\n",tangent_crack_front_nodes.size()); CHKERRQ(ierr);

    Range crack_surfaces_edge_faces;
    rval = m_field.get_moab().get_adjacencies(
      crack_front_nodes,2,false,crack_surfaces_edge_faces,Interface::UNION); CHKERR_PETSC(rval);
    crack_surfaces_edge_faces = crack_surfaces_edge_faces.subset_by_type(MBTRI);
    crack_surfaces_edge_faces = intersect(crack_surfaces_edge_faces,crack_surfaces_faces);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Front Faces = %d\n",crack_surfaces_edge_faces.size()); CHKERRQ(ierr);

    Range level_tris;
    ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
    crack_surfaces_edge_faces = intersect(crack_surfaces_edge_faces,level_tris);

    ierr = m_field.seed_finite_elements(crack_surfaces_edge_faces); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"C_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"CTC_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(crack_surfaces_edge_faces,"C_TANGENT_ELEM"); CHKERRQ(ierr);

    Range level_edges;
    ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBEDGE,level_edges); CHKERRQ(ierr);
    crack_front_edges = intersect(crack_front_edges,level_edges);
    ierr = m_field.seed_finite_elements(crack_front_edges); CHKERRQ(ierr);

    //add entitities (by tets) to the field
    ierr = m_field.add_ents_to_field_by_VERTICEs(crack_front_nodes,"LAMBDA_CRACKFRONT_AREA"); CHKERRQ(ierr);
    ierr = m_field.add_ents_to_field_by_VERTICEs(tangent_crack_front_nodes,"LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);

  }

  //NOTE: always order should be 1
  ierr = m_field.set_field_order(0,MBVERTEX,"LAMBDA_CRACKFRONT_AREA",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"LAMBDA_CRACK_TANGENT_CONSTRAIN",1); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element(problem,"dCT_CRACKFRONT_AREA_ELEM"); CHKERRQ(ierr);

  if(problem == "COUPLED_PROBLEM" || problem == "COUPLED_DYNAMIC") {
    ierr = m_field.modify_problem_add_finite_element(problem,"C_TANGENT_ELEM"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("MESH_SMOOTHER","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("MESH_SMOOTHER","LAMBDA_CRACK_TANGENT_CONSTRAIN"); CHKERRQ(ierr);
  }



  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::spatial_partition_problems(FieldInterface& m_field) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = m_field.partition_problem("ELASTIC_MECHANICS",1); CHKERRQ(ierr);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
  ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS",false,0,pcomm->size()); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::material_partition_problems(FieldInterface& m_field) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition MATERIAL_MECHANICS
  ierr = m_field.partition_problem("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("MATERIAL_MECHANICS"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::coupled_partition_problems(FieldInterface& m_field) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = m_field.partition_problem("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::constrains_partition_problems(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = m_field.simple_partition_problem("CCT_ALL_MATRIX",0); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = m_field.compose_problem("C_ALL_MATRIX","CCT_ALL_MATRIX",false,problem,true); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("C_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("C_ALL_MATRIX"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::crackfront_partition_problems(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  //partition
  ierr = m_field.simple_partition_problem("CTC_CRACKFRONT_MATRIX",0); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("CTC_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = m_field.compose_problem("C_CRACKFRONT_MATRIX","CTC_CRACKFRONT_MATRIX",false,problem,true); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("C_CRACKFRONT_MATRIX"); CHKERRQ(ierr);
  ierr = m_field.partition_ghost_dofs("C_CRACKFRONT_MATRIX"); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode ConfigurationalFractureMechanics::set_spatial_positions(FieldInterface& m_field) {
  PetscFunctionBegin;
  if(material_FirelWall->operator[](FW_set_spatial_positions)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_set_spatial_positions);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"SPATIAL_POSITION",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) {
      double &fval = dof_ptr->get_FieldData();
      fval = 0;
      continue;
    }
    EntityHandle ent = dof_ptr->get_ent();
    int dof_rank = dof_ptr->get_dof_rank();
    double &fval = dof_ptr->get_FieldData();
    if(node!=ent) {
      rval = m_field.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
      node = ent;
    }
    fval = coords[dof_rank];
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::set_material_positions(FieldInterface& m_field) {
  PetscFunctionBegin;

  if(material_FirelWall->operator[](FW_set_material_positions)) PetscFunctionReturn(0);
  material_FirelWall->set(FW_set_material_positions);

  ErrorCode rval;

  EntityHandle node = 0;
  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) {
      double &fval = dof_ptr->get_FieldData();
      fval = 0;
      continue;
    }
    EntityHandle ent = dof_ptr->get_ent();
    int dof_rank = dof_ptr->get_dof_rank();
    double &fval = dof_ptr->get_FieldData();
    if(node!=ent) {
	rval = m_field.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
    }
    fval = coords[dof_rank];
  }

  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::set_coordinates_from_material_solution(FieldInterface& m_field,bool only_crack_front) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ErrorCode rval;

  Range crack_front_edges,crack_front_edges_nodes;
  if(only_crack_front) {
    ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
    rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  }

  double coords[3];
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    if(only_crack_front) {
      if(
	crack_front_edges_nodes.find(dof_ptr->get_ent())
	==crack_front_edges_nodes.end()) {
	continue;
      }
    }
    EntityHandle ent = dof_ptr->get_ent();
    int dof_rank = dof_ptr->get_dof_rank();
    double fval = dof_ptr->get_FieldData();
    rval = m_field.get_moab().get_coords(&ent,1,coords); CHKERR_PETSC(rval);
    coords[dof_rank] = fval;
    rval = m_field.get_moab().set_coords(&ent,1,coords); CHKERR_PETSC(rval);
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"LAMBDA_SURFACE",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"LAMBDA_CRACK_SURFACE",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"LAMBDA_CRACK_TANGENT_CONSTRAIN",dof_ptr)) {
    dof_ptr->get_FieldData() = 0;
  }
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,ss.str(),dof_ptr)) {
      dof_ptr->get_FieldData() = 0;
    }
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::solve_spatial_problem(FieldInterface& m_field,SNES *snes,bool postproc) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  set_PhysicalEquationNumber(hooke);

  //create matrices
  Vec F;
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  PetscBool flg;
  PetscReal my_tol;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal atol,rtol,stol;
    PetscInt maxit,maxf;
    ierr = SNESGetTolerances(*snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);
    atol = my_tol;
    rtol = atol*1e2;
    ierr = SNESSetTolerances(*snes,atol,rtol,stol,maxit,maxf); CHKERRQ(ierr);
  }

  SnesCtx snes_ctx(m_field,"ELASTIC_MECHANICS");
  
  ierr = SNESSetApplicationContext(*snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);

  const double young_modulus = 1;
  const double poisson_ratio = 0.25;
  MyNonLinearSpatialElasticFEMthod my_fe(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_thermal_expansion",&my_fe.thermal_expansion,&flg); CHKERRQ(ierr);

  NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,Aij,F);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_forces = neumann_forces.getLoopSpatialFe();
  //fe_forces.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = fe_forces.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = fe_forces.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->get_msId(),true);  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new MetaNodalForces::TagForceScale(m_field));
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //postproc
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

  ierr = m_field.set_local_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = m_field.set_global_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  if(postproc) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags with nodal postions for spatial problem\n"); CHKERRQ(ierr);
    PostProcVertexMethod ent_method(m_field.get_moab(),"SPATIAL_POSITION");
    ierr = m_field.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",COL,ent_method); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags with residuals for spatial problem\n"); CHKERRQ(ierr);
    PostProcVertexMethod ent_method_res(m_field.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
    ierr = m_field.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",COL,ent_method_res); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Save tags on refined mesh with stresses\n"); CHKERRQ(ierr);
    ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);
    if(pcomm->rank()==0) {
      if(fe_post_proc_stresses_method!=NULL) delete fe_post_proc_stresses_method;
      fe_post_proc_stresses_method = new PostProcStressNonLinearElasticity(m_field.get_moab(),my_fe);
      fe_post_proc_stresses_method->do_broadcast = false;
      ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",*fe_post_proc_stresses_method,0,pcomm->size());  CHKERRQ(ierr);
    }
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Spatial problem solved\n"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::surface_projection_data(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  ierr = delete_surface_projection_data(m_field); CHKERRQ(ierr);

  //C_ALL_MATRIX is problem used to constrain surface constrain matrices
  if(projSurfaceCtx==NULL) {
    projSurfaceCtx = new matPROJ_ctx(m_field,problem,"C_ALL_MATRIX");
    ierr = m_field.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&projSurfaceCtx->C); CHKERRQ(ierr);
  }

  Range corners_edges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range corners_edges_nodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.merge(corners_edges_nodes);
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  //Loops over body surface (CFE_SURFACE) and crack surface (CFE_CRACK_SURFACE)
  ConstrainSurfacGeometry CFE_SURFACE(m_field,projSurfaceCtx->C);
  ConstrainSurfacGeometry CFE_CRACK_SURFACE(m_field,projSurfaceCtx->C,"LAMBDA_CRACK_SURFACE");

  map<int,ConstrainSurfacGeometry*> CFE_SURFACE_msId_ptr;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    CFE_SURFACE_msId_ptr[msId] = new ConstrainSurfacGeometry(m_field,projSurfaceCtx->C,ss.str());
  }

  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projSurfaceCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projSurfaceCtx->C); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("C_ALL_MATRIX","C_SURFACE_ELEM",CFE_SURFACE);  CHKERRQ(ierr);
  //ierr = m_field.loop_finite_elements("C_ALL_MATRIX","C_CRACK_SURFACE_ELEM",CFE_CRACK_SURFACE);  CHKERRQ(ierr);
  for(map<int,ConstrainSurfacGeometry*>::iterator mit = CFE_SURFACE_msId_ptr.begin();
    mit!=CFE_SURFACE_msId_ptr.end();mit++) {
    ostringstream ss0;
    ss0 << "C_SURFACE_ELEM_msId_" << mit->first;
    ierr = m_field.loop_finite_elements("C_ALL_MATRIX",ss0.str(),*(mit->second));  CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projSurfaceCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  for(map<int,ConstrainSurfacGeometry*>::iterator mit = CFE_SURFACE_msId_ptr.begin();
    mit!=CFE_SURFACE_msId_ptr.end();mit++) {
    delete mit->second;
  }

  /*{
    //Matrix View
    MatView(projSurfaceCtx->C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::delete_surface_projection_data(FieldInterface& m_field) {
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if(projSurfaceCtx!=NULL) {
    ierr = projSurfaceCtx->DestroyQorP(); CHKERRQ(ierr);
    ierr = projSurfaceCtx->DestroyQTKQ(); CHKERRQ(ierr);
    ierr = MatDestroy(&projSurfaceCtx->C); CHKERRQ(ierr);
    delete projSurfaceCtx;
    projSurfaceCtx = NULL;
  }

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::project_force_vector(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  //ErrorCode rval;

  Vec F_Material;
  ierr = m_field.VecCreateGhost(problem,ROW,&F_Material); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",ROW,F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  int M,m;
  ierr = VecGetSize(F_Material,&M); CHKERRQ(ierr);
  ierr = VecGetLocalSize(F_Material,&m); CHKERRQ(ierr);
  Mat Q;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTF_Material;
  ierr = VecDuplicate(F_Material,&QTF_Material); CHKERRQ(ierr);
  ierr = MatMult(Q,F_Material,QTF_Material); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(m_field.get_moab(),"MESH_NODE_POSITIONS",QTF_Material,"MATERIAL_FORCE_PROJECTED");
  ierr = m_field.loop_dofs(problem,"MESH_NODE_POSITIONS",ROW,ent_method_material); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",ROW,QTF_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  double nrm2_material;
  ierr = VecNorm(QTF_Material,NORM_2,&nrm2_material);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_QTF_Material = %6.4e\n",nrm2_material); CHKERRQ(ierr);

  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = VecDestroy(&QTF_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::project_form_th_projection_tag(FieldInterface& m_field,string problem,bool do_not_project) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;

  Vec D,dD;
  ierr = m_field.VecCreateGhost(problem,COL,&D); CHKERRQ(ierr);
  ierr = VecDuplicate(D,&dD); CHKERRQ(ierr);
  ierr = VecZeroEntries(dD); CHKERRQ(ierr);
  ierr = VecSetOption(dD,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);

  ierr = m_field.set_local_VecCreateGhost(problem,COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.set_global_VecCreateGhost(problem,COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  Tag th_projection;
  rval = m_field.get_moab().tag_get_handle("PROJECTION_CRACK_SURFACE",th_projection); CHKERR_PETSC(rval);

  const MoFEMProblem *problemPtr;
  ierr = m_field.get_problem(problem,&problemPtr); CHKERRQ(ierr);

  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);

  Range crack_front_edges;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  Range crack_front_edges_nodes;
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_edges_nodes,true); CHKERR_PETSC(rval);
  Range::iterator nit = crack_front_edges_nodes.begin();
  for(;nit!=crack_front_edges_nodes.end();nit++) {

    NumeredDofMoFEMEntity_multiIndex::index<Composite_Name_Ent_And_Part_mi_tag>::type::iterator dit,hi_dit;
    dit = problemPtr->numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().lower_bound(
      boost::make_tuple("MESH_NODE_POSITIONS",*nit,pcomm->rank()));
    hi_dit = problemPtr->numered_dofs_cols.get<Composite_Name_Ent_And_Part_mi_tag>().upper_bound(
      boost::make_tuple("MESH_NODE_POSITIONS",*nit,pcomm->rank()));
    if(distance(dit,hi_dit)==0) continue;
    ublas::vector<PetscInt,ublas::bounded_array<PetscInt,3> > indices;
    indices.resize(3);
    fill(indices.begin(),indices.end(),-1);
    ublas::vector<double,ublas::bounded_array<double,3> > new_coords;
    new_coords.resize(3);
    rval = m_field.get_moab().tag_get_data(th_projection,&*nit,1,&*new_coords.data().begin()); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,3> > delta;
    delta.resize(3);
    for(;dit!=hi_dit;dit++) {
      indices[dit->get_dof_rank()] = dit->get_petsc_gloabl_dof_idx();
      delta[dit->get_dof_rank()] = new_coords[dit->get_dof_rank()] - dit->get_FieldData();
    }
    ierr = VecSetValues(dD,3,&*indices.data().begin(),&*delta.data().begin(),INSERT_VALUES); CHKERRQ(ierr);

  }

  ierr = VecAssemblyBegin(dD); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(dD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  int M,m;
  ierr = VecGetSize(dD,&M); CHKERRQ(ierr);
  ierr = VecGetLocalSize(dD,&m); CHKERRQ(ierr);
  Mat Q;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
  ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);

  Vec QTdD;
  ierr = VecDuplicate(D,&QTdD); CHKERRQ(ierr);
  ierr = MatMult(Q,dD,QTdD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(QTdD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTdD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecAXPY(D,1.,QTdD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
 
  //not project nodes on crack surface 
  if(do_not_project) {
    ierr = m_field.set_local_VecCreateGhost(problem,COL,QTdD,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  }
  ierr = m_field.set_global_VecCreateGhost(problem,COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(m_field.get_moab(),"MESH_NODE_POSITIONS",D,"PROJECTION_CRACK_SURFACE");
  ierr = m_field.loop_dofs(problem,"MESH_NODE_POSITIONS",COL,ent_method); CHKERRQ(ierr);

  if(do_not_project) {
    ierr = m_field.set_global_VecCreateGhost(problem,COL,QTdD,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  }

  VecDestroy(&D);
  VecDestroy(&dD);
  VecDestroy(&QTdD);
  MatDestroy(&Q);

  PetscFunctionReturn(0);
}


PetscErrorCode ConfigurationalFractureMechanics::front_projection_data(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  
  ierr = delete_front_projection_data(m_field); CHKERRQ(ierr);

  if(projFrontCtx==NULL) {
    projFrontCtx = new matPROJ_ctx(m_field,problem,"C_CRACKFRONT_MATRIX");
    ierr = m_field.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&projFrontCtx->C); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::delete_front_projection_data(FieldInterface& m_field) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  if(projFrontCtx!=NULL) {
    ierr = projFrontCtx->DestroyQorP(); CHKERRQ(ierr);
    ierr = projFrontCtx->DestroyQTKQ(); CHKERRQ(ierr);
    ierr = MatDestroy(&projFrontCtx->C); CHKERRQ(ierr);
    delete projFrontCtx;
    projFrontCtx = NULL;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::griffith_force_vector(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  Vec LambdaVec,GriffithForceVec;
  ierr = m_field.VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,&LambdaVec); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("C_CRACKFRONT_MATRIX",COL,&GriffithForceVec); CHKERRQ(ierr);

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  //if(flg != PETSC_TRUE) {
    //SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  //}
  if(flg == PETSC_TRUE) {
    ierr = VecSet(LambdaVec,gc); CHKERRQ(ierr);
  } else {
    ierr = VecSet(LambdaVec,1); CHKERRQ(ierr);
  }

  Mat Q;
  {
    int M,m;
    ierr = VecGetSize(GriffithForceVec,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(GriffithForceVec,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
  }

  Range corners_edges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range corners_edges_nodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.merge(corners_edges_nodes);
  Range blocked_nodes;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,blocked_nodes); CHKERRQ(ierr);
  corners_nodes.merge(blocked_nodes);

  C_CONSTANT_AREA_FEMethod C_AREA_ELEM(m_field,projFrontCtx->C,Q,"LAMBDA_CRACKFRONT_AREA");

  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projFrontCtx->C); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_AREA_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatMultTranspose(projFrontCtx->C,LambdaVec,GriffithForceVec); CHKERRQ(ierr);

  Vec QTGriffithForceVec;
  ierr = VecDuplicate(GriffithForceVec,&QTGriffithForceVec); CHKERRQ(ierr);
  ierr = MatMult(Q,GriffithForceVec,QTGriffithForceVec); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE",ROW,QTGriffithForceVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_area(m_field.get_moab(),"MESH_NODE_POSITIONS",QTGriffithForceVec,"GRIFFITH_FORCE");
  ierr = m_field.loop_dofs(problem,"MESH_NODE_POSITIONS",COL,ent_method_area); CHKERRQ(ierr);

  double nrm2_griffith_force;
  ierr = VecNorm(QTGriffithForceVec,NORM_2,&nrm2_griffith_force);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_QTGriffithForceVec = %6.4e\n",nrm2_griffith_force); CHKERRQ(ierr);

  //tangent front froce
  matPROJ_ctx projFrontCtx_tangent(m_field,problem,"C_CRACKFRONT_MATRIX");
  ierr = m_field.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&projFrontCtx_tangent.C); CHKERRQ(ierr);
  C_FRONT_TANGENT_FEMethod C_TANGENT_ELEM(m_field,projFrontCtx_tangent.C,PETSC_NULL,"LAMBDA_CRACKFRONT_AREA");
  ierr = MatZeroEntries(projFrontCtx_tangent.C); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_TANGENT_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx_tangent.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx_tangent.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatMultTranspose(projFrontCtx_tangent.C,LambdaVec,GriffithForceVec); CHKERRQ(ierr);
  //ierr = MatMult(Q,GriffithForceVec,QTGriffithForceVec); CHKERRQ(ierr);
  //ierr = VecGhostUpdateBegin(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  //ierr = VecGhostUpdateEnd(QTGriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(GriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(GriffithForceVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE_TANGENT",ROW,GriffithForceVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method_tangent(
    m_field.get_moab(),"MESH_NODE_POSITIONS",GriffithForceVec,"GRIFFITH_TANGENT_FORCE");
  ierr = m_field.loop_dofs(problem,"MESH_NODE_POSITIONS",COL,ent_method_tangent); CHKERRQ(ierr);
  ierr = MatDestroy(&projFrontCtx_tangent.C); CHKERRQ(ierr);

  //cleaning
  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);
  ierr = VecDestroy(&GriffithForceVec); CHKERRQ(ierr);
  ierr = VecDestroy(&QTGriffithForceVec); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::griffith_g(FieldInterface& m_field,string problem) {
  PetscFunctionBegin;
  
  PetscErrorCode ierr;
  //ErrorCode rval;

  Vec F_Material;
  ierr = m_field.VecCreateGhost(problem,ROW,&F_Material); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",ROW,F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  Vec F_Griffith;
  ierr = m_field.VecCreateGhost(problem,ROW,&F_Griffith); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE",ROW,F_Griffith,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  Vec F_Griffith_Tangent;
  ierr = m_field.VecCreateGhost(problem,ROW,&F_Griffith_Tangent); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","GRIFFITH_FORCE_TANGENT",ROW,F_Griffith_Tangent,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  Vec LambdaVec,LambdaVec_Tangent;
  ierr = m_field.VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,&LambdaVec); CHKERRQ(ierr);
  ierr = VecDuplicate(LambdaVec,&LambdaVec_Tangent); CHKERRQ(ierr);

  Mat Q;
  {
    int M,m;
    ierr = VecGetSize(F_Material,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(F_Material,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,m,M,M,projSurfaceCtx,&Q); CHKERRQ(ierr);
    ierr = MatShellSetOperation(Q,MATOP_MULT,(void(*)(void))matQ_mult_shell); CHKERRQ(ierr);
  }

  ErrorCode rval;
  Range corners_edges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edges_nodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edges_nodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edges_nodes.begin(),corners_edges_nodes.end());
  Range blocked_nodes;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,blocked_nodes); CHKERRQ(ierr);
  corners_nodes.merge(blocked_nodes);

  C_CONSTANT_AREA_FEMethod C_AREA_ELEM(m_field,projFrontCtx->C,Q,"LAMBDA_CRACKFRONT_AREA");

  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
  ierr = MatSetOption(projFrontCtx->C,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);

  ierr = MatZeroEntries(projFrontCtx->C); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_AREA_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx->C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /*{
    //Matrix View
    MatView(projFrontCtx->C,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
    std::string wait;
    std::cin >> wait;
  }*/

  // R = CT*(CC^T)^(-1) [ unit m/m^(-2) = 1/m ] [ 0.5/0.25 = 2 ]
  // R^T = (CC^T)^(-T)C [ unit m/m^(-2) = 1/m ] 
  Mat RT;
  {
    int N,n;
    int M,m;
    ierr = VecGetSize(F_Material,&N); CHKERRQ(ierr);
    ierr = VecGetLocalSize(F_Material,&n); CHKERRQ(ierr);
    ierr = VecGetSize(LambdaVec,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(LambdaVec,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,projFrontCtx,&RT); CHKERRQ(ierr);
    ierr = MatShellSetOperation(RT,MATOP_MULT,(void(*)(void))matRT_mult_shell); CHKERRQ(ierr);
  }
  
  //calculate tangent griffith force
  matPROJ_ctx projFrontCtx_tangent(m_field,problem,"C_CRACKFRONT_MATRIX");
  ierr = m_field.MatCreateMPIAIJWithArrays("C_CRACKFRONT_MATRIX",&projFrontCtx_tangent.C); CHKERRQ(ierr);
  C_FRONT_TANGENT_FEMethod C_TANGENT_ELEM(m_field,projFrontCtx_tangent.C,Q,"LAMBDA_CRACKFRONT_AREA");
  ierr = MatZeroEntries(projFrontCtx_tangent.C); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("C_CRACKFRONT_MATRIX","C_CRACKFRONT_AREA_ELEM",C_TANGENT_ELEM);  CHKERRQ(ierr);
  ierr = MatAssemblyBegin(projFrontCtx_tangent.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(projFrontCtx_tangent.C,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  Mat RT_Tangent;
  {
    int N,n;
    int M,m;
    ierr = VecGetSize(F_Material,&N); CHKERRQ(ierr);
    ierr = VecGetLocalSize(F_Material,&n); CHKERRQ(ierr);
    ierr = VecGetSize(LambdaVec_Tangent,&M); CHKERRQ(ierr);
    ierr = VecGetLocalSize(LambdaVec_Tangent,&m); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&projFrontCtx_tangent,&RT_Tangent); CHKERRQ(ierr);
    ierr = MatShellSetOperation(RT_Tangent,MATOP_MULT,(void(*)(void))matRT_mult_shell); CHKERRQ(ierr);
  }

  //clualte griffith forces
  double gc = 1;
  PetscBool flg_gc;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg_gc); CHKERRQ(ierr);

  ierr = projFrontCtx->InitQorP(F_Material); CHKERRQ(ierr);
  ierr = projFrontCtx_tangent.InitQorP(F_Material); CHKERRQ(ierr);

  // unit of LambdaVec [ N * 1/m = N*m/m^2 = J/m^2 ]
  ierr = VecScale(F_Material,-1./gc); CHKERRQ(ierr);
  ierr = MatMult(RT,F_Material,LambdaVec); CHKERRQ(ierr);  
  ierr = MatMult(RT_Tangent,F_Material,LambdaVec_Tangent); CHKERRQ(ierr);  
  ierr = VecScale(F_Material,gc); CHKERRQ(ierr);
  ierr = VecScale(LambdaVec,gc); CHKERRQ(ierr);
  ierr = VecScale(LambdaVec_Tangent,gc); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(LambdaVec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(LambdaVec_Tangent,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(LambdaVec_Tangent,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //get pointer to problem
  const MoFEMProblem *problemPtr;
  ierr = m_field.get_problem("C_CRACKFRONT_MATRIX",&problemPtr); CHKERRQ(ierr);
  //make maps to save griffith forces
  ierr = m_field.set_global_VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,LambdaVec,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  map_ent_g.clear();
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problemPtr,"LAMBDA_CRACKFRONT_AREA",diit)) {
    map_ent_g[diit->get_ent()] = diit->get_FieldData();
  }
  ierr = m_field.set_global_VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,LambdaVec_Tangent,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  map<EntityHandle,double> g_tangent_map;
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problemPtr,"LAMBDA_CRACKFRONT_AREA",diit)) {
    g_tangent_map[diit->get_ent()] = diit->get_FieldData();
  }

  //vector of matrerial force magnitudes
  Vec JVec;
  ierr = VecDuplicate(LambdaVec,&JVec); CHKERRQ(ierr);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);

  ublas::vector<double> coords;
  coords.resize(3);
  PetscPrintf(PETSC_COMM_WORLD,"\n\ngriffith force\n\n");

  Tag th_g,th_g_tangent;
  const double def_val = 0;
  rval = m_field.get_moab().tag_get_handle(
    "G1",1,MB_TYPE_DOUBLE,th_g,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_THROW(rval);
  rval = m_field.get_moab().tag_get_handle(
    "G3",1,MB_TYPE_DOUBLE,th_g_tangent,MB_TAG_CREAT|MB_TAG_SPARSE,&def_val); CHKERR_THROW(rval);

  map_ent_j.clear();
  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_NAME_FOR_LOOP_(problemPtr,"LAMBDA_CRACKFRONT_AREA",diit)) {
    EntityHandle ent = diit->get_ent();
    double g_val = map_ent_g[ent];
    double g_tangent_val = fabs(g_tangent_map[ent]);
    rval = m_field.get_moab().tag_set_data(th_g,&ent,1,&g_val); CHKERR_PETSC(rval);
    rval = m_field.get_moab().tag_set_data(th_g_tangent,&ent,1,&g_tangent_val); CHKERR_PETSC(rval);
    rval = m_field.get_moab().get_coords(&ent,1,&*coords.data().begin()); CHKERR_PETSC(rval);
    int dd = 0;
    ublas::vector<double,ublas::bounded_array<double,9> > material_force(3);
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(m_field,"MATERIAL_FORCE",diit->get_ent(),diiit)) {
      material_force[diiit->get_dof_rank()] = diiit->get_FieldData();
      dd++;
    }
    if(dd != 3) SETERRQ1(PETSC_COMM_SELF,1,"can not find material force vector at node %ld",diit->get_ent());
    ublas::vector<double,ublas::bounded_array<double,9> > griffith_force(3);
    for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(m_field,"GRIFFITH_FORCE",diit->get_ent(),diiit)) {
      griffith_force[diiit->get_dof_rank()] = diiit->get_FieldData();
      dd++;
    }
    if(dd != 6) SETERRQ1(PETSC_COMM_SELF,1,"can not find griffith force vector at node %ld",diit->get_ent());
    double j = norm_2(material_force)/(norm_2(griffith_force)/gc);

    j = copysign(j,g_val);
    map_ent_j[ent] = j;
    if(diit->get_part()==pcomm->rank()) {
      ierr = VecSetValue(JVec,diit->get_petsc_gloabl_dof_idx(),j,INSERT_VALUES); CHKERRQ(ierr);
    }
    ostringstream ss;
    ss << "griffith force at ent ";
    ss << setw(5) << ent;
    //coords
    ss << "\tcoords";
    ss << " " << setw(10) << setprecision(4) << coords[0];
    ss << " " << setw(10) << setprecision(4) << coords[1];
    ss << " " << setw(10) << setprecision(4) << coords[2];
    //g 
    ss << "\t\tg1 " << scientific << setprecision(4) << g_val;
    ss << " / " << scientific << setprecision(4) << j;
    ss << " ( " << scientific << setprecision(4) << g_val/j << " )";
    //g tangent
    ss << "\t\tg3 " << scientific << setprecision(4) << g_tangent_val;
    ss << " ( " << scientific << setprecision(4) << g_tangent_val/j << " )";
    if(flg_gc == PETSC_TRUE) {
      ss << "\t relative error (gc-g)/gc ";
      ss << scientific << setprecision(4) << (gc-g_val)/gc;
      ss << " / " << scientific << setprecision(4) << (gc-j)/gc;
    }
    ss << endl; 
    PetscPrintf(PETSC_COMM_WORLD,"%s",ss.str().c_str());
  }

  ierr = VecSum(LambdaVec,&ave_g); CHKERRQ(ierr);
  ierr = VecMin(LambdaVec,PETSC_NULL,&min_g); CHKERRQ(ierr);
  ierr = VecMax(LambdaVec,PETSC_NULL,&max_g); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JVec); CHKERRQ(ierr);
  ierr = VecSum(JVec,&ave_j); CHKERRQ(ierr);
  ierr = VecMin(JVec,PETSC_NULL,&min_j); CHKERRQ(ierr);
  ierr = VecMax(JVec,PETSC_NULL,&max_j); CHKERRQ(ierr);
  ierr = VecDestroy(&JVec); CHKERRQ(ierr);

  Range crackSurfacesFaces;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,crackSurfacesFaces,true); CHKERRQ(ierr);
  Range level_tris;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
  crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);

  double aRea = 0;
  ublas::vector<double,ublas::bounded_array<double,6> > diffNTRI(6);
  ShapeDiffMBTRI(&diffNTRI.data()[0]);
  for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = m_field.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,9> > dofsX(9);
    ublas::vector<double,ublas::bounded_array<double,3> > normal(3);
    for(int nn = 0;nn<num_nodes; nn++) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(m_field,"MESH_NODE_POSITIONS",conn[nn],dit)) {
	dofsX[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
      }
    }
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    //crack surface area is a half of crack top and bottom body surface
    aRea += norm_2(normal)*0.25;
  }

  {
    int N;
    ierr = VecGetSize(LambdaVec,&N); CHKERRQ(ierr);
    ave_g /= N;
    ave_j /= N;
  }
  PetscPrintf(PETSC_COMM_WORLD,"\naverage griffith force %6.4e / %6.4e Crack surface area %6.4e\n",ave_g,ave_j,aRea);
  PetscPrintf(PETSC_COMM_WORLD,"\n\n");

  PostProcVertexMethod ent_method(m_field.get_moab(),"LAMBDA_CRACKFRONT_AREA");
  ierr = m_field.loop_dofs("C_CRACKFRONT_MATRIX","LAMBDA_CRACKFRONT_AREA",ROW,ent_method); CHKERRQ(ierr);

  ierr = MatDestroy(&Q); CHKERRQ(ierr);
  ierr = MatDestroy(&RT); CHKERRQ(ierr);
  ierr = MatDestroy(&RT_Tangent); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Griffith); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Griffith_Tangent); CHKERRQ(ierr);
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec); CHKERRQ(ierr);
  ierr = VecDestroy(&LambdaVec_Tangent); CHKERRQ(ierr);

  ierr = projFrontCtx_tangent.DestroyQorP(); CHKERRQ(ierr);
  ierr = projFrontCtx_tangent.DestroyQTKQ(); CHKERRQ(ierr);
  ierr = MatDestroy(&projFrontCtx_tangent.C); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode MySnesConvernceTest(SNES snes,int it,double xnorm,double gnorm,double fnorm,SNESConvergedReason *reason,void *void_ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = SNESConvergedDefault(snes,it,xnorm,gnorm,fnorm,reason,PETSC_NULL); CHKERRQ(ierr);
  const PetscReal div = 10e3;
  if(fnorm > div) {
    *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode MySnesConvernceTest_SNESLINESEARCHBT(SNES snes,int it,double xnorm,double gnorm,double fnorm,SNESConvergedReason *reason,void *void_ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = SNESConvergedDefault(snes,it,xnorm,gnorm,fnorm,reason,PETSC_NULL); CHKERRQ(ierr);
  const PetscReal div = 10e3;
  if(fnorm > div) {
    *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
  }
  if(it>0) {
    SNESLineSearch linesearch;
    ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBT); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode MySnesConvernceTest_SNESLINESEARCHL2(SNES snes,int it,double xnorm,double gnorm,double fnorm,SNESConvergedReason *reason,void *void_ctx) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = SNESConvergedDefault(snes,it,xnorm,gnorm,fnorm,reason,PETSC_NULL); CHKERRQ(ierr);
  const PetscReal div = 10e3;
  if(fnorm > div) {
    *reason = SNES_DIVERGED_FUNCTION_DOMAIN;
  }
  if(it>0) {
    SNESLineSearch linesearch;
    ierr = SNESGetLineSearch(snes,&linesearch); CHKERRQ(ierr);
    ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::fix_all_but_one(FieldInterface& m_field,Range &fix_nodes,const double fraction_treshold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  ErrorCode rval;
  Tag th_freez;
  const int def_order = 0;
  rval = m_field.get_moab().tag_get_handle("FROZEN_NODE",1,MB_TYPE_INTEGER,th_freez,MB_TAG_CREAT|MB_TAG_SPARSE,&def_order); CHKERR_PETSC(rval);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"freeze front nodes:\n");
  double max_g_j;
  EntityHandle max_g_j_ent = 0;
  for(
    map<EntityHandle,double>::iterator mit = map_ent_j.begin();
    mit!=map_ent_j.end();mit++) {
    double fraction = (max_j-mit->second)/max_j;
    double g_j = map_ent_g[mit->first]/mit->second;
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "front node = %ld max_j = %6.4e j = %6.4e (%6.4e) g/j = %4.3f",
      mit->first,max_j,mit->second,fraction,g_j); CHKERRQ(ierr);
    bool freez_or_not_to_freez;
    if(fraction > fraction_treshold) {
      freez_or_not_to_freez = true;
    } else {
      freez_or_not_to_freez = false;
      if(max_g_j_ent==0) {
	max_g_j = g_j;
	max_g_j_ent = mit->first;
      } else {
	if(g_j < max_g_j) {
	  max_g_j = g_j;
	  max_g_j_ent = mit->first;
	}
      }
    }
    if(freeze_all_but_one) {
      freez_or_not_to_freez = freeze_all_but_one;
    } 
    if(freez_or_not_to_freez) {
      ierr = PetscPrintf(PETSC_COMM_WORLD," freeze\n");
      fix_nodes.insert(mit->first);
      int freez = 1;
      EntityHandle node = mit->first;
      rval = m_field.get_moab().tag_set_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
    } else {
      int freez;
      EntityHandle node = mit->first;
      rval = m_field.get_moab().tag_get_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
      if(freez == 0) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
      } else if(freez == 1) {
	ierr = PetscPrintf(PETSC_COMM_WORLD," unfreeze\n");
	int freez = 0;
	rval = m_field.get_moab().tag_set_data(th_freez,&node,1,&freez); CHKERR_PETSC(rval);
      } 
    }
  }
  if(freeze_all_but_one) {
    if(max_g_j_ent == 0) {
      SETERRQ(PETSC_COMM_SELF,1,"imosible case, at least one node shold not break");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nunfreez %ld g/j = %4.3f\n",max_g_j_ent,max_g_j);
    int freez = 0;
    rval = m_field.get_moab().tag_set_data(th_freez,&max_g_j_ent,1,&freez); CHKERR_PETSC(rval);
    fix_nodes.erase(max_g_j_ent);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::solve_coupled_problem(FieldInterface& m_field,SNES *snes,double da,const double fraction_treshold) {
  PetscFunctionBegin;

  PetscErrorCode ierr;
  ErrorCode rval;
  PetscBool flg;

  set_PhysicalEquationNumber(hooke);

  ierr = front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  //create matrices
  Mat K;
  ierr = m_field.MatCreateMPIAIJWithArrays("COUPLED_PROBLEM",&K); CHKERRQ(ierr);
  //create vectors
  Vec F;
  ierr = m_field.VecCreateGhost("COUPLED_PROBLEM",ROW,&F); CHKERRQ(ierr);
  Vec D;
  ierr = m_field.VecCreateGhost("COUPLED_PROBLEM",COL,&D); CHKERRQ(ierr);

  if(material_FirelWall->operator[](FW_arc_lenhghat_definition)) {
  } else {
    SETERRQ(PETSC_COMM_SELF,1,"arc length not initialised)");
  }

  Range corners_edges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,corners_edges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  Range corners_edgesNodes;
  rval = m_field.get_moab().get_connectivity(corners_edges,corners_edgesNodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(corners_edgesNodes.begin(),corners_edgesNodes.end());
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  Range crack_front_edges,crack_front_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRQ(ierr);
  rval = m_field.get_moab().get_connectivity(crack_front_edges,crack_front_nodes,true); CHKERR_PETSC(rval);
  Range level_nodes;
  ierr = m_field.get_entities_by_type_and_ref_level(*ptr_bit_level0,BitRefLevel().set(),MBVERTEX,level_nodes); CHKERRQ(ierr);
  crack_front_nodes = intersect(crack_front_nodes,level_nodes);
  //corners_nodes.merge(subtract(level_nodes,crack_front_nodes));

  Range fix_nodes;
  ierr = fix_all_but_one(m_field,fix_nodes,fraction_treshold); CHKERRQ(ierr);
  da *= fabs(crack_front_nodes.size()-fix_nodes.size())/(double)crack_front_nodes.size();
  corners_nodes.merge(fix_nodes);

  struct MyPrePostProcessFEMethod: public FEMethod {
    
    FieldInterface& m_field;
    ArcLengthCtx *arc_ptr;

    MyPrePostProcessFEMethod(FieldInterface& _m_field,
      ArcLengthCtx *_arc_ptr): 
      m_field(_m_field),arc_ptr(_arc_ptr) {}
  
    PetscErrorCode ierr;
      
    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  snes_ctx = CTX_SNESSETFUNCTION;
	  snes_f = ts_F;
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  snes_ctx = CTX_SNESSETJACOBIAN;
	  snes_B = ts_B;
	  break;
	}
	default:
	break;
      }
        
      //PetscAttachDebugger();
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
          ierr = VecZeroEntries(snes_f); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(snes_f,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecZeroEntries(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
        
      PetscFunctionReturn(0);
    }
      
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
	  //snes_f
          ierr = VecGhostUpdateBegin(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(snes_f,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
        }
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }
      PetscFunctionReturn(0);
    }

  };

  struct AssembleLambdaFEMethod: public FEMethod {
    
    FieldInterface& m_field;
    ArcLengthCtx *arc_ptr;

    SpatialPositionsBCFEMethodPreAndPostProc *bC;

    AssembleLambdaFEMethod(FieldInterface& _m_field,
      ArcLengthCtx *_arc_ptr,SpatialPositionsBCFEMethodPreAndPostProc *bc): 
      m_field(_m_field),arc_ptr(_arc_ptr),bC(bc) {}
  
    PetscErrorCode ierr;
      
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      switch (ts_ctx) {
	case CTX_TSSETIFUNCTION: {
	  snes_ctx = CTX_SNESSETFUNCTION;
	  snes_f = ts_F;
	  break;
	}
	case CTX_TSSETIJACOBIAN: {
	  snes_ctx = CTX_SNESSETJACOBIAN;
	  snes_B = ts_B;
	  break;
	}
	default:
	break;
      }
      PetscFunctionReturn(0);
    }
    PetscErrorCode operator()() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }
    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      switch(snes_ctx) {
        case CTX_SNESSETFUNCTION: {
	  //F_lambda
          ierr = VecGhostUpdateBegin(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecGhostUpdateEnd(arc_ptr->F_lambda,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
          ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
          ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  for(vector<int>::iterator vit = bC->dofsIndices.begin();vit!=bC->dofsIndices.end();vit++) {
	    ierr = VecSetValue(arc_ptr->F_lambda,*vit,0,INSERT_VALUES); CHKERRQ(ierr);
	  }
	  ierr = VecAssemblyBegin(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(arc_ptr->F_lambda); CHKERRQ(ierr);
	  ierr = VecDot(arc_ptr->F_lambda,arc_ptr->F_lambda,&arc_ptr->F_lambda2); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tFlambda2 = %6.4e\n",arc_ptr->F_lambda2);
	  //add F_lambda
	  ierr = VecAXPY(snes_f,arc_ptr->get_FieldData(),arc_ptr->F_lambda); CHKERRQ(ierr);
	  PetscPrintf(PETSC_COMM_WORLD,"\tlambda = %6.4e\n",arc_ptr->get_FieldData());  
	  double fnorm;
	  ierr = VecNorm(snes_f,NORM_2,&fnorm); CHKERRQ(ierr);	
	  PetscPrintf(PETSC_COMM_WORLD,"\tfnorm = %6.4e\n",fnorm);  
	}
        break;
        default:
          SETERRQ(PETSC_COMM_SELF,1,"not implemented");
      }  
      PetscFunctionReturn(0);
    }

  };

  //arc elem
  ArcLengthCtx arc_ctx(m_field,"COUPLED_PROBLEM");
  ArcLengthSnesCtx arc_snes_ctx(m_field,"COUPLED_PROBLEM",&arc_ctx);
  ArcLengthMatShell arc_mat_ctx(m_field,K,&arc_ctx,"COUPLED_PROBLEM");
  PCShellCtx pc_ctx(K,K,&arc_ctx);
  ArcLengthElemFEMethod arc_elem(m_field,this,&arc_ctx);
  //spatial and material forces
  const double young_modulus = 1;
  const double poisson_ratio = 0.;
  MyNonLinearSpatialElasticFEMthod fe_spatial(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_spatial.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_spatial.isCoupledProblem = true;
  MyEshelbyFEMethod fe_material(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  ierr = fe_material.initCrackFrontData(m_field); CHKERRQ(ierr);
  fe_material.isCoupledProblem = true;
  //meshs moothing
  MyMeshSmoothingFEMethod smoother(m_field);
  ierr = smoother.initCrackFrontData(m_field); CHKERRQ(ierr);
  set_qual_ver(3);
  //constrains
  SnesConstrainSurfacGeometry constrain_body_surface(m_field,"LAMBDA_SURFACE");
  constrain_body_surface.nonlinear = true;
  SnesConstrainSurfacGeometry constrain_crack_surface(m_field,"LAMBDA_CRACK_SURFACE");
  constrain_crack_surface.nonlinear = true;
  Snes_CTgc_CONSTANT_AREA_FEMethod ct_gc(m_field,projFrontCtx,"COUPLED_PROBLEM","LAMBDA_CRACKFRONT_AREA");
  Snes_dCTgc_CONSTANT_AREA_FEMethod dct_gc(m_field,K,"LAMBDA_CRACKFRONT_AREA");
  TangentWithMeshSmoothingFrontConstrain_FEMethod tangent_constrain(m_field,&smoother,"LAMBDA_CRACK_TANGENT_CONSTRAIN");
  map<int,SnesConstrainSurfacGeometry*> other_body_surface_constrains;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    other_body_surface_constrains[msId] = new SnesConstrainSurfacGeometry(m_field,ss.str());
    other_body_surface_constrains[msId]->nonlinear = true;
  }
  //bothsieds constrains
  BothSurfaceConstrains  both_sides_constrains(m_field);
  //dirichlet constrains
  FixBcAtEntities fix_material_pts(m_field,"MESH_NODE_POSITIONS",corners_nodes);
  fix_material_pts.fieldNames.push_back("LAMBDA_SURFACE");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_SURFACE");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT");
  fix_material_pts.fieldNames.push_back("LAMBDA_CRACK_TANGENT_CONSTRAIN");
  fix_material_pts.fieldNames.push_back("LAMBDA_BOTH_SIDES");
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
    int msId = it->get_msId();
    if((msId < 10200)||(msId >= 10300)) continue;
    ostringstream ss;
    ss << "LAMBDA_SURFACE_msId_" << msId;
    fix_material_pts.fieldNames.push_back(ss.str());
  }
  SpatialPositionsBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"SPATIAL_POSITION",K,D,F);
  my_dirichlet_bc.fixFields.push_back("MESH_NODE_POSITIONS");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_SURFACE");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_CRACK_SURFACE_WITH_CRACK_FRONT");
  my_dirichlet_bc.fixFields.push_back("LAMBDA_BOTH_SIDES");
  //boundary conditions
  Tag th_scale;
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_scale); CHKERR_PETSC(rval);
  double *force_scale;
  rval = m_field.get_moab().tag_get_by_ptr(th_scale,&root_meshset,1,(const void**)&force_scale); CHKERR_PETSC(rval);
  arc_ctx.get_FieldData() = (*force_scale);
  double scaled_reference_load = 1;
  double *scale_lhs = &(arc_ctx.get_FieldData());
  double *scale_rhs = &(scaled_reference_load);
  NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,K,arc_ctx.F_lambda,scale_lhs,scale_rhs);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_forces = neumann_forces.getLoopSpatialFe();
  //fe_forces.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
  fe_forces.uSeF = true; 
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = fe_forces.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    ierr = fe_forces.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }
  //portsproc
  MyPrePostProcessFEMethod pre_post_method(m_field,&arc_ctx);
  AssembleLambdaFEMethod assemble_F_lambda(m_field,&arc_ctx,&my_dirichlet_bc);

  //rhs
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&pre_post_method);
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&fix_material_pts);
  arc_snes_ctx.get_preProcess_to_do_Rhs().push_back(&ct_gc);
  SnesCtx::loops_to_do_type& loops_to_do_Rhs = arc_snes_ctx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",arc_ctx.F_lambda,it->get_msId());  CHKERRQ(ierr);
  }
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  //arc length
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NoNE",&assemble_F_lambda));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_elem));

  arc_snes_ctx.get_postProcess_to_do_Rhs().push_back(&fix_material_pts);
  arc_snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
  arc_snes_ctx.get_postProcess_to_do_Rhs().push_back(&pre_post_method);

  //lhs
  arc_snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
  arc_snes_ctx.get_preProcess_to_do_Mat().push_back(&fix_material_pts);
  SnesCtx::loops_to_do_type& loops_to_do_Mat = arc_snes_ctx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("BOTH_SIDE_OF_CRACK",&both_sides_constrains));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("dCT_CRACKFRONT_AREA_ELEM",&dct_gc));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_SURFACE_ELEM",&constrain_body_surface));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("CandCT_CRACK_SURFACE_ELEM",&constrain_crack_surface));
  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    ostringstream ss2;
    ss2 << "CandCT_SURFACE_ELEM_msId_" << mit->first;
    loops_to_do_Mat.push_back(SnesCtx::loop_pair_type(ss2.str(),mit->second));
  }
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MESH_SMOOTHER",&smoother));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("C_TANGENT_ELEM",&tangent_constrain));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC_COUPLED",&fe_spatial));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("MATERIAL_COUPLED",&fe_material));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_forces));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ARC_LENGTH",&arc_elem));
  arc_snes_ctx.get_postProcess_to_do_Mat().push_back(&fix_material_pts);
  arc_snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

  PetscInt M,N;
  ierr = MatGetSize(K,&M,&N); CHKERRQ(ierr);
  PetscInt m,n;
  MatGetLocalSize(K,&m,&n);
  Mat ShellK;
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,(void*)&arc_mat_ctx,&ShellK); CHKERRQ(ierr);
  ierr = MatShellSetOperation(ShellK,MATOP_MULT,(void(*)(void))arc_length_mult_shell); CHKERRQ(ierr);

  ierr = SNESSetApplicationContext(*snes,&arc_snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,SnesRhs,&arc_snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,ShellK,K,SnesMat,&arc_snes_ctx); CHKERRQ(ierr);
  PetscReal my_tol;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_tol",&my_tol,&flg); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal atol,rtol,stol;
    PetscInt maxit,maxf;
    ierr = SNESGetTolerances(*snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);
    atol = my_tol;
    rtol = atol*1e2;
    ierr = SNESSetTolerances(*snes,atol,rtol,stol,maxit,maxf); CHKERRQ(ierr);
  }

  KSP ksp;
  ierr = SNESGetKSP(*snes,&ksp); CHKERRQ(ierr);
  if(flg == PETSC_TRUE) {
    PetscReal rtol,atol,dtol;
    PetscInt maxits;
    ierr = KSPGetTolerances(ksp,&rtol,&atol,&dtol,&maxits); CHKERRQ(ierr);
    atol = my_tol*1e-2;
    rtol = atol*1e-2;
    ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxits); CHKERRQ(ierr);
  }
  PC pc;
  ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
  ierr = PCShellSetContext(pc,&pc_ctx); CHKERRQ(ierr);
  ierr = PCShellSetApply(pc,pc_apply_arc_length); CHKERRQ(ierr);
  ierr = PCShellSetSetUp(pc,pc_setup_arc_length); CHKERRQ(ierr);

  ierr = m_field.set_local_VecCreateGhost("COUPLED_PROBLEM",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = m_field.get_problem("COUPLED_PROBLEM",&(arc_elem.problemPtr)); CHKERRQ(ierr);
  ierr = arc_elem.set_dlambda_to_x(D,0); CHKERRQ(ierr);
  ierr = VecCopy(D,arc_ctx.x0); CHKERRQ(ierr);
  ierr = arc_elem.get_dlambda(D); CHKERRQ(ierr);
  //calculate rhs and F_lambda
  ierr = arc_ctx.set_alpha_and_beta(0,1); CHKERRQ(ierr);
  ierr = SnesRhs(*snes,D,F,&arc_snes_ctx); CHKERRQ(ierr);
  //set s
  ierr = arc_ctx.set_alpha_and_beta(1,0); CHKERRQ(ierr);
  ierr = arc_elem.calculate_lambda_int(); CHKERRQ(ierr);
  ierr = arc_ctx.set_s(arc_elem.lambda_int+da/arc_elem.aRea0); CHKERRQ(ierr);

  SNESLineSearch linesearch;
  ierr = SNESGetLineSearch(*snes,&linesearch); CHKERRQ(ierr);
  PetscReal atol,rtol,stol;
  PetscInt maxit,maxf;
  ierr = SNESGetTolerances(*snes,&atol,&rtol,&stol,&maxit,&maxf); CHKERRQ(ierr);

  Vec D0;
  ierr = m_field.VecCreateGhost("COUPLED_PROBLEM",COL,&D0); CHKERRQ(ierr);
  ierr = m_field.set_local_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  SNESConvergedReason reason;

  ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(*snes,MySnesConvernceTest_SNESLINESEARCHBT,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
  ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
  total_its = its;
  ierr = SNESGetConvergedReason(*snes,&reason); CHKERRQ(ierr);
  if(reason < 0) {
    ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
    ierr = SNESSetConvergenceTest(*snes,MySnesConvernceTest_SNESLINESEARCHL2,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
    total_its += its;
    ierr = SNESGetConvergedReason(*snes,&reason); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
    if(reason < 0) {
      ierr = m_field.set_global_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
      ierr = VecCopy(D0,D); CHKERRQ(ierr);
      smoother.stabilise = true;
      for(int ii = 1;ii<=3;ii++) {
	ierr = VecCopy(D,D0); CHKERRQ(ierr);
	smoother.alpha22 = exp(-ii);
	ierr = SNESLineSearchSetType(linesearch,SNESLINESEARCHL2); CHKERRQ(ierr);
	ierr = SNESSetConvergenceTest(*snes,MySnesConvernceTest_SNESLINESEARCHL2,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
	ierr = SNESSolve(*snes,PETSC_NULL,D); CHKERRQ(ierr);
	ierr = SNESGetIterationNumber(*snes,&its); CHKERRQ(ierr);
	total_its += its;
	ierr = SNESGetConvergedReason(*snes,&reason); CHKERRQ(ierr);
	if(reason < 0) {
	  ierr = m_field.set_global_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = VecCopy(D0,D); CHKERRQ(ierr);
	  break;
	} else {
	  ierr = m_field.set_global_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	  ierr = set_coordinates_from_material_solution(m_field,false); CHKERRQ(ierr);
	}
      }
      smoother.stabilise = false;
    }
  }
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  *force_scale = (arc_ctx.get_FieldData());
  ierr = VecDestroy(&D0); CHKERRQ(ierr);

  //Save data on mesh
  ierr = m_field.set_global_VecCreateGhost("COUPLED_PROBLEM",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  aRea = arc_elem.aRea;

  ierr = m_field.set_field(0,MBVERTEX,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.set_field(0,MBEDGE,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.set_field(0,MBTRI,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.set_field(0,MBTET,"SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.field_axpy(+1.,"SPATIAL_POSITION","SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.field_axpy(-1.,"MESH_NODE_POSITIONS","SPATIAL_DISPLACEMENT"); CHKERRQ(ierr);

  Vec DISP;
  ierr = VecDuplicate(D,&DISP); CHKERRQ(ierr);
  ierr = m_field.set_other_global_VecCreateGhost(
    "COUPLED_PROBLEM","SPATIAL_POSITION","SPATIAL_DISPLACEMENT",COL,DISP,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecDot(DISP,arc_ctx.F_lambda,&energy); CHKERRQ(ierr);
  lambda = arc_ctx.get_FieldData();
  energy = 0.5*fabs(lambda)*fabs(energy);

  int verb = 1;
  if(reason>=0) {

    //Save field on mesh
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_POSITION on tags\n"); CHKERRQ(ierr);
    }
    PostProcVertexMethod ent_method_spatial(m_field.get_moab(),"SPATIAL_POSITION");
    ierr = m_field.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",COL,ent_method_spatial); CHKERRQ(ierr);
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save MESH_NODE_POSITIONS on tags\n"); CHKERRQ(ierr);
    }
    PostProcVertexMethod ent_method_material(m_field.get_moab(),"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",COL,ent_method_material); CHKERRQ(ierr);
    if(verb>2) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save MATERIAL_RESIDUAL on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_res_mat(m_field.get_moab(),"MESH_NODE_POSITIONS",F,"MATERIAL_RESIDUAL");
      ierr = m_field.loop_dofs("COUPLED_PROBLEM","MESH_NODE_POSITIONS",COL,ent_method_res_mat); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_RESIDUAL on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_res_spat(m_field.get_moab(),"SPATIAL_POSITION",F,"SPATIAL_RESIDUAL");
      ierr = m_field.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",COL,ent_method_res_spat); CHKERRQ(ierr);
    }
    if(verb>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Save SPATIAL_DISPLACEMENT on tags\n"); CHKERRQ(ierr);
      PostProcVertexMethod ent_method_disp_spat(m_field.get_moab(),"SPATIAL_POSITION",DISP,"SPATIAL_DISPLACEMENT");
      ierr = m_field.loop_dofs("COUPLED_PROBLEM","SPATIAL_POSITION",COL,ent_method_disp_spat); CHKERRQ(ierr);
    }
 
  }

  ierr = VecDestroy(&DISP); CHKERRQ(ierr);
  ierr = MatDestroy(&ShellK); CHKERRQ(ierr);
  ierr = MatDestroy(&K); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);

  for(map<int,SnesConstrainSurfacGeometry*>::iterator mit = other_body_surface_constrains.begin();
    mit!=other_body_surface_constrains.end();mit++) {
    delete mit->second;
  }

  ierr = delete_front_projection_data(m_field); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode ConfigurationalFractureMechanics::calculate_material_forces(FieldInterface& m_field,string problem,string fe) {
  PetscFunctionBegin;

  PetscErrorCode ierr;

  set_PhysicalEquationNumber(hooke);

  Range CornersEdges,corners_nodes;
  ierr = m_field.get_Cubit_msId_entities_by_dimension(100,SIDESET,1,CornersEdges,true); CHKERRQ(ierr);
  ierr = m_field.get_Cubit_msId_entities_by_dimension(101,NODESET,0,corners_nodes,true); CHKERRQ(ierr);
  ErrorCode rval;
  Range CornersEdgesNodes;
  rval = m_field.get_moab().get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
  corners_nodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
  Range nodes_to_block;
  BitRefLevel bit_to_block = BitRefLevel().set(BITREFLEVEL_SIZE-2);
  ierr = m_field.get_entities_by_type_and_ref_level(bit_to_block,BitRefLevel().set(),MBVERTEX,nodes_to_block); CHKERRQ(ierr);
  corners_nodes.merge(nodes_to_block);

  Vec F_Material;
  ierr = m_field.VecCreateGhost(problem,COL,&F_Material); CHKERRQ(ierr);

  const double young_modulus = 1;
  const double poisson_ratio = 0;
  MyEshelbyFEMethod material_fe(m_field,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  PetscBool flg = PETSC_TRUE;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_thermal_expansion",&material_fe.thermal_expansion,&flg); CHKERRQ(ierr);
  material_fe.snes_f = F_Material;
  material_fe.set_snes_ctx(SnesMethod::CTX_SNESSETFUNCTION);
  FixBcAtEntities fix_material_pts(m_field,"MESH_NODE_POSITIONS",corners_nodes);
  fix_material_pts.snes_x = PETSC_NULL;
  fix_material_pts.snes_f = F_Material;
  fix_material_pts.set_snes_ctx(SnesMethod::CTX_SNESSETFUNCTION);

  ierr = VecZeroEntries(F_Material); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = m_field.problem_basic_method_preProcess(problem,fix_material_pts); CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements(problem,fe,material_fe);  CHKERRQ(ierr);
  ierr = m_field.problem_basic_method_postProcess(problem,fix_material_pts); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F_Material,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F_Material); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F_Material); CHKERRQ(ierr);

  PostProcVertexMethod ent_method_material(m_field.get_moab(),"MESH_NODE_POSITIONS",F_Material,"MATERIAL_FORCE");
  ierr = m_field.loop_dofs(problem,"MESH_NODE_POSITIONS",COL,ent_method_material); CHKERRQ(ierr);

  double nrm2_material;
  ierr = VecNorm(F_Material,NORM_2,&nrm2_material);   CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"nrm2_F_Material = %6.4e\n",nrm2_material); CHKERRQ(ierr);

  //Fields
  ierr = m_field.set_other_global_VecCreateGhost(
    problem,"MESH_NODE_POSITIONS","MATERIAL_FORCE",
    ROW,F_Material,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  
  //detroy matrices
  ierr = VecDestroy(&F_Material); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
ConfigurationalFractureMechanics::ArcLengthElemFEMethod::ArcLengthElemFEMethod(
  FieldInterface& _mField,ConfigurationalFractureMechanics *_conf_prob,ArcLengthCtx *_arc_ptr): 
    mField(_mField),conf_prob(_conf_prob),arc_ptr(_arc_ptr) {
    PetscErrorCode ierr;
    ErrorCode rval;

    //ghost dofs for diagonal value
    PetscInt ghosts[1] = { 0 };
    Interface &moab = mField.get_moab();
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm->rank() == 0) {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,1,1,0,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      ierr = VecCreateGhost(PETSC_COMM_WORLD,0,1,1,ghosts,&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    }

    //get crack surcace
    ierr = mField.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,crackSurfacesFaces,true); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    Range level_tris;
    ierr = mField.get_entities_by_type_and_ref_level(*conf_prob->ptr_bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    crackSurfacesFaces = intersect(crackSurfacesFaces,level_tris);
    //get crack surface faces adjacenet to crack front
    Range crack_front_edges;
    ierr = mField.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,crack_front_edges,true); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    mField.get_moab().get_connectivity(crack_front_edges,crackFrontNodes,true);

    //get pointer to coupled problem
    ierr = mField.get_problem("COUPLED_PROBLEM",&problemPtr); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    set<DofIdx> set_idx; //set of global dofs indexes for crack surface
    for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
      const EntityHandle* conn; 
      int num_nodes; 
      rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); 
      if (MB_SUCCESS != rval) { 
	PetscTraceBackErrorHandler(
	  PETSC_COMM_WORLD,
	  __LINE__,PETSC_FUNCTION_NAME,__FILE__,
	  MOFEM_DATA_INSONSISTENCY,PETSC_ERROR_INITIAL,"can not get connectibility",PETSC_NULL);
	CHKERRABORT(PETSC_COMM_SELF,rval);
      }
      for(int nn = 0;nn<num_nodes; nn++) {
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problemPtr,conn[nn],dit)) {
	  if(dit->get_name()!="MESH_NODE_POSITIONS") continue;
	  set_idx.insert(dit->get_petsc_gloabl_dof_idx());
	}
      }
    }

    //vector surfaceDofs keeps material positions for crack surface dofs
    ierr = PetscMalloc(set_idx.size()*sizeof(PetscInt),&isIdx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    copy(set_idx.begin(),set_idx.end(),isIdx);
    ierr = ISCreateGeneral(PETSC_COMM_SELF,set_idx.size(),isIdx,PETSC_USE_POINTER,&isSurface); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecCreateSeq(PETSC_COMM_SELF,set_idx.size(),&surfaceDofs); CHKERRABORT(PETSC_COMM_WORLD,ierr);

    ierr = mField.VecCreateGhost("C_CRACKFRONT_MATRIX",ROW,&lambdaVec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    ierr = VecSet(lambdaVec,1); CHKERRABORT(PETSC_COMM_WORLD,ierr);

}
ConfigurationalFractureMechanics::ArcLengthElemFEMethod::~ArcLengthElemFEMethod() {
  PetscErrorCode ierr;
  ierr = VecScatterDestroy(&surfaceScatter); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&ghostDiag); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&surfaceDofs); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = ISDestroy(&isSurface); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = PetscFree(isIdx); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  ierr = VecDestroy(&lambdaVec); CHKERRABORT(PETSC_COMM_WORLD,ierr);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calculate_area() {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  //scatter from snes_x to surfaceDofs
  ierr = VecScatterCreate(snes_x,isSurface,surfaceDofs,PETSC_NULL,&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(surfaceScatter,snes_x,surfaceDofs,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&surfaceScatter); CHKERRABORT(PETSC_COMM_SELF,ierr);

  //get crack front nodal positions form surfaceDofs
  double *array;
  ierr = VecGetArray(surfaceDofs,&array); CHKERRQ(ierr);  
  PetscInt size;
  ISGetSize(isSurface,&size);
  for(int ii = 0;ii<size;ii++) {
    problemPtr->numered_dofs_rows.get<
      PetscGlobalIdx_mi_tag>().find(isIdx[ii])->get_FieldData() = array[ii];
  }
  ierr = VecRestoreArray(surfaceDofs,&array); CHKERRABORT(PETSC_COMM_SELF,ierr);

  //calculate crack surface area
  aRea = 0;
  aRea0 = 0;
  ublas::vector<double,ublas::bounded_array<double,6> > diffNTRI(6);
  ShapeDiffMBTRI(&diffNTRI.data()[0]);
  for(Range::iterator fit = crackSurfacesFaces.begin();fit!=crackSurfacesFaces.end();fit++) {
    const EntityHandle* conn; 
    int num_nodes; 
    rval = mField.get_moab().get_connectivity(*fit,conn,num_nodes,true); CHKERR_PETSC(rval);
    ublas::vector<double,ublas::bounded_array<double,9> > dofsX(9);
    ierr = mField.get_moab().get_coords(conn,num_nodes,&*dofsX.data().begin()); CHKERRQ(ierr);
    ublas::vector<double,ublas::bounded_array<double,3> > normal(3);
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    aRea0 += norm_2(normal)*0.25;
    for(int nn = 0;nn<num_nodes; nn++) {
      if(find(crackFrontNodes.begin(),crackFrontNodes.end(),conn[nn])!=crackFrontNodes.end()) {
      for(_IT_GET_DOFS_FIELD_BY_NAME_AND_ENT_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",conn[nn],dit)) {
	dofsX[nn*3+dit->get_dof_rank()] = dit->get_FieldData();
      }}
    }
    ierr = ShapeFaceNormalMBTRI(&diffNTRI[0],&dofsX.data()[0],&normal.data()[0]); CHKERRQ(ierr);
    //crack surface area is a half of crack top and bottom body surface
    aRea += norm_2(normal)*0.25;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calculate_lambda_int() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  ierr = calculate_area(); CHKERRQ(ierr);
  lambda_int = arc_ptr->alpha*aRea/aRea0 + arc_ptr->dlambda*arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
  //PetscPrintf(PETSC_COMM_WORLD,"\tsurface crack area = %6.4e lambda_int = %6.4e\n",aRea,lambda_int);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::calculate_db() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  if(arc_ptr->alpha!=0) {
    ierr = MatMultTranspose(conf_prob->projFrontCtx->C,lambdaVec,arc_ptr->db); CHKERRQ(ierr);
    ierr = VecScale(arc_ptr->db,1./aRea0); CHKERRQ(ierr);
    if(arc_ptr->alpha!=1) {
      ierr = VecScale(arc_ptr->db,arc_ptr->alpha); CHKERRQ(ierr);
    }
  } else {
    ierr = VecZeroEntries(arc_ptr->db); CHKERRQ(ierr);
  }   
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::set_dlambda_to_x(Vec x,double dlambda) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //check if locl dof idx is non zero, i.e. that lambda is acessible from this processor
  if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
    double *array;
    ierr = VecGetArray(x,&array); CHKERRQ(ierr);
    double lambda_old = array[arc_ptr->get_petsc_local_dof_idx()];
    if(!(dlambda == dlambda)) {
      ostringstream sss;
      sss << "s " << arc_ptr->s << " " << arc_ptr->beta << " " << arc_ptr->F_lambda2;
      SETERRQ(PETSC_COMM_SELF,1,sss.str().c_str());
    }
    array[arc_ptr->get_petsc_local_dof_idx()] = lambda_old + dlambda;
    PetscPrintf(PETSC_COMM_WORLD,"\tload factor lambda_old,lambda_new/dlambda = %6.4e, %6.4e (%6.4e)\n",lambda_old, array[arc_ptr->get_petsc_local_dof_idx()], dlambda);
    ierr = VecRestoreArray(x,&array); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::get_dlambda(Vec x) {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //set dx
  ierr = VecCopy(x,arc_ptr->dx); CHKERRQ(ierr);
  ierr = VecAXPY(arc_ptr->dx,-1,arc_ptr->x0); CHKERRQ(ierr);
  //if LAMBDA dof is on this partition
  if(arc_ptr->get_petsc_local_dof_idx()!=-1) {
    double *array;
    ierr = VecGetArray(arc_ptr->dx,&array); CHKERRQ(ierr);
    arc_ptr->dlambda = array[arc_ptr->get_petsc_local_dof_idx()];
    array[arc_ptr->get_petsc_local_dof_idx()] = 0;
    ierr = VecRestoreArray(arc_ptr->dx,&array); CHKERRQ(ierr);
  }
  //brodcast dlambda
  unsigned int part = arc_ptr->get_part();
  //MPI_Bcast(&(arc_ptr->dlambda),1,MPI_DOUBLE,part,PETSC_COMM_WORLD);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
  Vec lambda_ghost;
  if(pcomm->rank()==part) {
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD,1,1,0,PETSC_NULL,&(arc_ptr->dlambda),&lambda_ghost); CHKERRQ(ierr);
  } else {
    int one[] = {0};
    ierr = VecCreateGhostWithArray(PETSC_COMM_WORLD,0,1,1,one,&(arc_ptr->dlambda),&lambda_ghost); CHKERRQ(ierr);
  }
  ierr = VecGhostUpdateBegin(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(lambda_ghost,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecDestroy(&lambda_ghost); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"\tload factor increment dlambda = %6.4e\n",arc_ptr->dlambda);
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::preProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch (ts_ctx) {
    case CTX_TSSETIFUNCTION: {
      snes_ctx = CTX_SNESSETFUNCTION;
      snes_f = ts_F;
      break;
    }
    case CTX_TSSETIJACOBIAN: {
      snes_ctx = CTX_SNESSETJACOBIAN;
      snes_B = ts_B;
      break;
    }
    default:
    break;
  }
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: { 
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = get_dlambda(snes_x); CHKERRQ(ierr);
	ierr = calculate_lambda_int(); CHKERRQ(ierr);
      }
      break;
    case CTX_SNESSETJACOBIAN: 
	ierr = MatAssemblyBegin(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = calculate_db(); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::operator()() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  //get dlambda dof 
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: {
      //calculate residual for arc length row
      arc_ptr->res_lambda = lambda_int - arc_ptr->s;
      ierr = VecSetValue(snes_f,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->res_lambda,INSERT_VALUES); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_SELF,"\n");
      PetscPrintf(PETSC_COMM_SELF,"\t** Arc-Length residual:\n");
      PetscPrintf(PETSC_COMM_SELF,"\t  residual of arc-length control res_lambda = %6.4e crack area/f_lambda_int = %6.4e ( %6.4e )\n"
	,arc_ptr->res_lambda,aRea,aRea/aRea0);
    } break; 
    case CTX_SNESSETJACOBIAN: {
      //calculate diagonal therm
      double diag = arc_ptr->beta*sqrt(arc_ptr->F_lambda2);
      ierr = VecSetValue(ghostDiag,0,diag,INSERT_VALUES); CHKERRQ(ierr);
      ierr = MatSetValue(snes_B,arc_ptr->get_petsc_gloabl_dof_idx(),arc_ptr->get_petsc_gloabl_dof_idx(),1,INSERT_VALUES); CHKERRQ(ierr);
    } break;
    default:
    break;
  }	
  PetscFunctionReturn(0);
}
PetscErrorCode ConfigurationalFractureMechanics::ArcLengthElemFEMethod::postProcess() {
  PetscFunctionBegin;
  PetscErrorCode ierr;
  switch(snes_ctx) {
    case CTX_SNESSETFUNCTION: { 
	ParallelComm* pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
	double res_nrm2[6];
	Vec res_nrm2_vec;
	if(pcomm->rank()==0) {
	  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,6,6,res_nrm2,&res_nrm2_vec); CHKERRQ(ierr);
	} else {
	  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,1,0,6,res_nrm2,&res_nrm2_vec); CHKERRQ(ierr);
	}
	//
	Range CrackFrontEdges;
	ierr = mField.get_Cubit_msId_entities_by_dimension(201,SIDESET,1,CrackFrontEdges,true); CHKERRQ(ierr);
	Range CrackFrontNodes;
	ErrorCode rval;
	rval = mField.get_moab().get_connectivity(CrackFrontEdges,CrackFrontNodes,true); CHKERR_PETSC(rval);
	//
	double *array;
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
	ierr = VecGetArray(snes_f,&array); CHKERRQ(ierr);
	ierr = VecZeroEntries(res_nrm2_vec); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(res_nrm2_vec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(res_nrm2_vec,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(problemPtr,dof)) {
	  if(dof->get_part()!=pcomm->rank()) continue;
	  double val = pow(array[dof->get_petsc_local_dof_idx()],2);
	  if(dof->get_name() == "SPATIAL_POSITION") {
	    ierr = VecSetValue(res_nrm2_vec,0,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "MESH_NODE_POSITIONS") {
	    if(find(CrackFrontNodes.begin(),CrackFrontNodes.end(),dof->get_ent())!=CrackFrontNodes.end()) {
	      ierr = VecSetValue(res_nrm2_vec,1,val,ADD_VALUES); CHKERRQ(ierr);
	    } else {
	      ierr = VecSetValue(res_nrm2_vec,2,val,ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	  if(dof->get_name().compare(0,14,"LAMBDA_SURFACE") == 0) {
	    ierr = VecSetValue(res_nrm2_vec,3,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "LAMBDA_CRACK_SURFACE") {
	    ierr = VecSetValue(res_nrm2_vec,4,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	  if(dof->get_name() == "LAMBDA_CRACK_TANGENT_CONSTRAIN") {
	    ierr = VecSetValue(res_nrm2_vec,5,val,ADD_VALUES); CHKERRQ(ierr);
	  }
	}
	ierr = VecAssemblyBegin(res_nrm2_vec); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(res_nrm2_vec); CHKERRQ(ierr);
	resSpatialNrm2 = sqrt(res_nrm2[0]);
	resCrackFrontNrm2 = sqrt(res_nrm2[1]);
	//if(pcomm->rank()==0) {
	  PetscPrintf(PETSC_COMM_WORLD,"\t** Equilibrium Residuals:\n");
	  PetscPrintf(PETSC_COMM_WORLD,"\t  equilibrium of physical forces res_spatial_nrm2 = %6.4e\n",resSpatialNrm2);
	  PetscPrintf(PETSC_COMM_WORLD,"\t  equilibrium of material forces res_crack_front_nrm2 = %6.4e\n",resCrackFrontNrm2);
	  PetscPrintf(PETSC_COMM_WORLD,"\t** Mesh Quality Residuals:\n");
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of mesh smoother res_mesh_smoother_nrm2 = %6.4e\n",sqrt(res_nrm2[2]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of surface constraint res_surface_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[3]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of crack surface constraint res_crack_surface_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[4]));
	  PetscPrintf(PETSC_COMM_WORLD,"\t  residual of crack front mesh smoother tangential constraint res_crack_front_tangent_constrain_nrm2 = %6.4e\n",sqrt(res_nrm2[5]));
	  PetscPrintf(PETSC_COMM_WORLD,"\n");
	//}
	ierr = VecRestoreArray(snes_f,&array); CHKERRQ(ierr);
	ierr = VecDestroy(&res_nrm2_vec); CHKERRQ(ierr);
    } break;
    case CTX_SNESSETJACOBIAN: {
	ierr = VecAssemblyBegin(ghostDiag); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(ghostDiag); CHKERRQ(ierr);
	ierr = VecGhostUpdateBegin(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	ierr = VecGhostUpdateEnd(ghostDiag,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	double *diag;
	ierr = VecGetArray(ghostDiag,&diag); CHKERRQ(ierr);
	arc_ptr->diag = *diag;
	ierr = VecRestoreArray(ghostDiag,&diag); CHKERRQ(ierr);
	//PetscPrintf(PETSC_COMM_WORLD,"\tdiag = %6.4e\n",arc_ptr->diag);
	/*ierr = MatAssemblyBegin(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(snes_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
 	//Matrix View
	MatView(snes_B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
	std::string wait;
	std::cin >> wait;*/
    } break;
    default:
	break;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode main_spatial_solution(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); 
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  EntityHandle meshset_level0;
  rval = m_field.get_moab().create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);

  //solve problem
  ierr = conf_prob.set_spatial_positions(m_field); CHKERRQ(ierr);
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  ierr = conf_prob.solve_spatial_problem(m_field,&snes); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}
PetscErrorCode main_material_forces(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr; 

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);

  //caculate material forces
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(m_field,"MATERIAL_MECHANICS","MATERIAL"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_force_vector(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_g(m_field,"MATERIAL_MECHANICS"); CHKERRQ(ierr);
  ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
  ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_setup(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;
  //PetscBool flg = PETSC_TRUE;

  ierr = conf_prob.set_material_fire_wall(m_field); CHKERRQ(ierr);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  ierr = conf_prob.spatial_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.coupled_problem_definition(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_crack_front_problem_definition(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.arclength_problem_definition(m_field); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL_COUPLED",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MATERIAL",MBTET); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_set_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("MATERIAL_MECHANICS",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("COUPLED_PROBLEM",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("C_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_ref_level_set_bit("CTC_CRACKFRONT_MATRIX",bit_level0); CHKERRQ(ierr);

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = conf_prob.spatial_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.material_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.coupled_partition_problems(m_field); CHKERRQ(ierr);
  ierr = conf_prob.constrains_partition_problems(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.crackfront_partition_problems(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_restart(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  Tag th_griffith_force;
  rval = m_field.get_moab().tag_get_handle("GRIFFITH_FORCE",th_griffith_force); CHKERR_PETSC(rval);
  rval = m_field.get_moab().tag_delete(th_griffith_force); CHKERR_PETSC(rval);
  Tag th_freez;
  rval = m_field.get_moab().tag_get_handle("FROZEN_NODE",th_freez); CHKERR_PETSC(rval);
  rval = m_field.get_moab().tag_delete(th_freez); CHKERR_PETSC(rval);

  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_spatial_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_material_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_coupled_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_constrains_crack_front_problem_definition) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
  conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;

  ierr = main_arc_length_setup(m_field,conf_prob); CHKERRQ(ierr);
  ierr = m_field.check_number_of_ents_in_ents_field(); CHKERRQ(ierr);
  ierr = m_field.check_number_of_ents_in_ents_finite_element(); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_rescale_load_factor(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  double gc;
  PetscBool flg;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
  }

  const EntityHandle root_meshset = m_field.get_moab().get_root_set();

  Tag th_t_val;
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  double max_j = conf_prob.max_j;
  double a = fabs(max_j)/pow(load_factor,2);
  double new_load_factor = copysign(sqrt(gc/a),load_factor);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\ncoefficient a = %6.4e\n",a); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"new load factor value = %6.4e\n\n",new_load_factor); CHKERRQ(ierr);
  load_factor = new_load_factor;
  SNES snes;
  //solve spatial problem
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  ierr = conf_prob.solve_spatial_problem(m_field,&snes,false); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
  ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
  ierr = conf_prob.griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode main_arc_length_solve(FieldInterface& m_field,ConfigurationalFractureMechanics& conf_prob) {
  PetscFunctionBegin;

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscBool flg = PETSC_TRUE;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&m_field.get_moab(),MYPCOMM_INDEX);

  //ref meshset ref level 0
  Tag th_my_ref_level;
  rval = m_field.get_moab().tag_get_handle("_MY_REFINMENT_LEVEL",th_my_ref_level); CHKERR_PETSC(rval);
  const EntityHandle root_meshset = m_field.get_moab().get_root_set();
  BitRefLevel *ptr_bit_level0;
  rval = m_field.get_moab().tag_get_by_ptr(th_my_ref_level,&root_meshset,1,(const void**)&ptr_bit_level0); CHKERR_PETSC(rval);
  BitRefLevel& bit_level0 = *ptr_bit_level0;

  //caculate material forces
  ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);

  double da_0 = 0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-my_da",&da_0,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_da (what is the crack area increment ?)");
  }
  double da = da_0;

  int def_nb_load_steps = 0;
  Tag th_nb_load_steps;
  rval = m_field.get_moab().tag_get_handle("_NB_LOAD_STEPS",1,MB_TYPE_INTEGER,
    th_nb_load_steps,MB_TAG_CREAT|MB_TAG_SPARSE,&def_nb_load_steps); 
  int *ptr_nb_load_steps;
  rval = m_field.get_moab().tag_get_by_ptr(th_nb_load_steps,&root_meshset,1,(const void**)&ptr_nb_load_steps); CHKERR_PETSC(rval);
  int &step = *ptr_nb_load_steps;
  if(step!=0) {
    step++;
  }

  PetscInt nb_load_steps = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_load_steps",&nb_load_steps,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_load_steps (what is the number of load_steps ?)");
  }

  if(step == 0) {
    //ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
  }

  Tag th_t_val;
  rval = m_field.get_moab().tag_get_handle("_LoadFactor_Scale_",th_t_val); CHKERR_PETSC(rval);
  double *load_factor_ptr;
  rval = m_field.get_moab().tag_get_by_ptr(th_t_val,&root_meshset,1,(const void**)&load_factor_ptr); CHKERR_THROW(rval);
  double& load_factor = *load_factor_ptr;

  FaceSplittingTools face_splitting_tools(m_field);

  for(int aa = 0;step<nb_load_steps;step++,aa++) {

    PetscLogDouble t1,t2;
    PetscLogDouble v1,v2;
    ierr = PetscTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n** number of step = %D\n",step); CHKERRQ(ierr);

    if(aa == 0) {
      ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
      ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
    }

    double gc;
    PetscBool flg;
    ierr = PetscOptionsGetReal(PETSC_NULL,"-my_gc",&gc,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_gc (what is the fracture energy ?)");
    }

    //calculate initial load factor
    if(step == 0) {
      ierr = main_rescale_load_factor(m_field,conf_prob); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,"* da = %6.4e\n",da); CHKERRQ(ierr);
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    
    Vec D0;
    ierr = m_field.VecCreateGhost("COUPLED_PROBLEM",COL,&D0); CHKERRQ(ierr);
    ierr = m_field.set_local_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    double load_factor0 = load_factor;
    
    int nb_sub_steps;
    ierr = PetscOptionsGetInt("","-my_nb_sub_steps",&nb_sub_steps,&flg); CHKERRQ(ierr);
    if(flg != PETSC_TRUE) {
      nb_sub_steps = 20;
    }

    int its_d;
    ierr = PetscOptionsGetInt("","-my_its_d",&its_d,&flg); CHKERRQ(ierr);

    bool at_least_one_step_converged = false;
    conf_prob.freeze_all_but_one = false;
    double _da_ = (aa == 0) ? 0 : da;
    int ii = 0;
    for(;ii<nb_sub_steps;ii++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* number of substeps = %D _da_ = %6.4e\n\n",ii,_da_); CHKERRQ(ierr);
      ierr = conf_prob.solve_coupled_problem(m_field,&snes,_da_); CHKERRQ(ierr);
      if(conf_prob.total_its == 0) break;
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason); CHKERRQ(ierr);
      if(reason > 0) {
	if(da > 0) {
	  if(aa > 0 && ii == 0) {
	    if(flg != PETSC_TRUE) {
	      its_d = 5;
	    }
	    double gamma = 0.5,reduction = 1;
	    reduction = pow((double)its_d/(double)(conf_prob.total_its+1),gamma);
	    const double max_da_reduction = 10;
	    if(reduction<1 || da < max_da_reduction*da_0) {
	      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e\n\n",reduction); CHKERRQ(ierr);
	      da = fmin(da*reduction,max_da_reduction*da_0);
	    }
	  }
	}
	at_least_one_step_converged = true;
	conf_prob.freeze_all_but_one = false;
	ierr = m_field.set_local_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	load_factor0 = load_factor;
	ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
	ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
	_da_ = 0;
      } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* reset unknowns vector\n"); CHKERRQ(ierr);
	ierr = m_field.set_global_VecCreateGhost("COUPLED_PROBLEM",COL,D0,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	load_factor = load_factor0;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, recalculate spatial positions\n"); CHKERRQ(ierr);
	{
	  SNES snes_spatial;
	  ierr = SNESCreate(PETSC_COMM_WORLD,&snes_spatial); CHKERRQ(ierr);  
	  ierr = SNESSetFromOptions(snes_spatial); CHKERRQ(ierr);
	  ierr = conf_prob.solve_spatial_problem(m_field,&snes_spatial,false); CHKERRQ(ierr);
	  ierr = SNESDestroy(&snes_spatial); CHKERRQ(ierr);
	  ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
	  ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	  ierr = conf_prob.griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
	}
	if(!at_least_one_step_converged) {
	  if(da > 0 && aa > 0) {
	    da = 0.5*da;
	    _da_ = da;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* failed to converge, set da = %6.4e ( 0.5 )\n",_da_); CHKERRQ(ierr);
	  } 
	}
	///PetscAttachDebugger();

	if(_da_ == 0) {
	  if(conf_prob.freeze_all_but_one) {
	    SETERRQ(PETSC_COMM_SELF,1,"* unable to converge");
	  } else {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"* freez all but one\n"); CHKERRQ(ierr);
	    conf_prob.freeze_all_but_one = true;
	  }
	}
      }

    }
    ierr = VecDestroy(&D0); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    double reduction = pow(2./(fmax(ii+1,3)-1),0.5);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n* change of da = %6.4e (nb. sub-steps %d) \n\n",reduction,ii); CHKERRQ(ierr);
    da *= reduction;
    da = fmax(da,1e-1*da_0);

    ierr = conf_prob.set_coordinates_from_material_solution(m_field,false); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "load_path: %4D Area %6.4e Lambda %6.4e Energy %6.4e\n",
      step,conf_prob.aRea,conf_prob.lambda,conf_prob.energy); CHKERRQ(ierr);
    if(pcomm->rank()==0) {
      EntityHandle out_meshset;
      rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
      ierr = m_field.problem_get_FE("COUPLED_PROBLEM","MATERIAL_COUPLED",out_meshset); CHKERRQ(ierr);
      ostringstream ss;
      ss << "out_load_step_" << step << ".vtk";
      rval = m_field.get_moab().write_file(ss.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
      rval = m_field.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);

      {
	Range SurfacesFaces;
	ierr = m_field.get_Cubit_msId_entities_by_dimension(102,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
	Range CrackSurfacesFaces;
	ierr = m_field.get_Cubit_msId_entities_by_dimension(200,SIDESET,2,CrackSurfacesFaces,true); CHKERRQ(ierr);
	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,SIDESET,it)) {
	  int msId = it->get_msId();
	  if((msId < 10200)||(msId >= 10300)) continue;
	  Range SurfacesFaces_msId;
	  ierr = m_field.get_Cubit_msId_entities_by_dimension(msId,SIDESET,2,SurfacesFaces_msId,true); CHKERRQ(ierr);
	  SurfacesFaces.insert(SurfacesFaces_msId.begin(),SurfacesFaces_msId.end());
	}

	Range level_tris;
	ierr = m_field.get_entities_by_type_and_ref_level(bit_level0,BitRefLevel().set(),MBTRI,level_tris); CHKERRQ(ierr);
	SurfacesFaces = intersect(SurfacesFaces,level_tris);
	CrackSurfacesFaces = intersect(CrackSurfacesFaces,level_tris);

	EntityHandle out_meshset;
	rval = m_field.get_moab().create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
	rval = m_field.get_moab().add_entities(out_meshset,CrackSurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss1;
	ss1 << "out_crack_surface_" << step << ".vtk";
	rval = m_field.get_moab().write_file(ss1.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = m_field.get_moab().add_entities(out_meshset,SurfacesFaces); CHKERR_PETSC(rval);
	ostringstream ss2;
	ss2 << "out_surface_" << step << ".vtk";
	rval = m_field.get_moab().write_file(ss2.str().c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
	rval = m_field.get_moab().delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
      }

      ostringstream ss3;
      ss3 << "restart_" << step << ".h5m";
      rval = m_field.get_moab().write_file(ss3.str().c_str()); CHKERR_PETSC(rval);

      const MoFEMProblem *problemPtr;
      ierr = m_field.get_problem("COUPLED_PROBLEM",&problemPtr); CHKERRQ(ierr);

      for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|UNKNOWNCUBITNAME,it)) {
	if(it->get_Cubit_name() != "LoadPath") continue;

	Range nodes;
	rval = m_field.get_moab().get_entities_by_type(it->meshset,MBVERTEX,nodes,true); CHKERR_PETSC(rval);
	for(Range::iterator nit = nodes.begin();nit!=nodes.end();nit++) {
	  double coords[3];
	  rval = m_field.get_moab().get_coords(&*nit,1,coords); CHKERR_PETSC(rval);
	  for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_ENT_FOR_LOOP_(problemPtr,*nit,dof)) {
	    if(dof->get_name()!="SPATIAL_POSITION") continue;
	    ierr = PetscPrintf(PETSC_COMM_WORLD,
	      "load_path_disp ent %ld dim %d "
	      "coords ( %8.6f %8.6f %8.6f ) "
	      "val %6.4e Lambda %6.4e\n",
	      dof->get_ent(),dof->get_dof_rank(),
	      coords[0],coords[1],coords[2],
	      dof->get_FieldData()-coords[dof->get_dof_rank()],      
	      conf_prob.lambda); CHKERRQ(ierr);
	  }
	}

      }
    }

    ierr = PetscTime(&v2);CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Step Time Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

    #ifdef WITH_TETGEN

      {
	Range edges_to_cat;
	ierr = face_splitting_tools.getCornerEdges(edges_to_cat,0); CHKERRQ(ierr);
	if(edges_to_cat.size()>0) {
	  Range new_nodes;
	  ierr = face_splitting_tools.propagateBySplit(new_nodes,edges_to_cat,0); CHKERRQ(ierr);
	  ierr = face_splitting_tools.conerProblem(new_nodes,0); CHKERRQ(ierr);
	}
      }

      {
	face_splitting_tools.moabTetGenMap.clear();
	face_splitting_tools.tetGenMoabMap.clear();
	face_splitting_tools.tetGenData.clear();
	vector<string> switches1;
	if(pcomm->rank() == 0) {
	  switches1.push_back("rp178sqRS0JVV");
	  ierr = face_splitting_tools.rebuildMeshWithTetGen(switches1,0); CHKERRQ(ierr);	
	} else {
	  switches1.push_back("rp178sqRS0JQ");
	  ierr = face_splitting_tools.rebuildMeshWithTetGen(switches1,0); CHKERRQ(ierr);	
	}
      }

      bit_level0 = BitRefLevel().set(face_splitting_tools.meshIntefaceBitLevels.back());
      //retart analysis
      ierr = main_arc_length_restart(m_field,conf_prob); CHKERRQ(ierr);
      //project and set coords
      conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_spatial_positions) = 0;
      conf_prob.material_FirelWall->operator[](ConfigurationalFractureMechanics::FW_set_material_positions) = 0;
      ierr = conf_prob.set_spatial_positions(m_field); CHKERRQ(ierr);
      ierr = conf_prob.set_material_positions(m_field); CHKERRQ(ierr);
      //solve spatial problem
      ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);  
      ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
      ierr = conf_prob.solve_spatial_problem(m_field,&snes,false); CHKERRQ(ierr);
      ierr = SNESDestroy(&snes); CHKERRQ(ierr);
      //calculate Griffth forces
      ierr = conf_prob.calculate_material_forces(m_field,"COUPLED_PROBLEM","MATERIAL_COUPLED"); CHKERRQ(ierr);
      ierr = conf_prob.front_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.surface_projection_data(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.project_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_force_vector(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.griffith_g(m_field,"COUPLED_PROBLEM"); CHKERRQ(ierr);
      ierr = conf_prob.delete_surface_projection_data(m_field); CHKERRQ(ierr);
      ierr = conf_prob.delete_front_projection_data(m_field); CHKERRQ(ierr);
    #endif


  }

  PetscFunctionReturn(0);
}

#ifdef WITH_ADOL_C

#include <ConfigurationalFractureForDynamics.cpp>

#endif //WITH_ADOL_C

