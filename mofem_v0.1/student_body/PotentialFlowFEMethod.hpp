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


#ifndef __POTENTIALFLOWFEMETHOD_HPP__
#define __POTENTIALFLOWFEMETHOD_HPP__

#include <boost/numeric/ublas/symmetric.hpp>
#include "FEMethod_UpLevelStudent.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"
extern "C" {
#include <gm_rule.h>
}

namespace MoFEM {

struct LaplacianElem: public FEMethod_UpLevelStudent {

    FieldInterface& mField;
    Mat A;
    Vec F;
    LaplacianElem(FieldInterface& _mField,Mat _A,Vec _F): 
      FEMethod_UpLevelStudent(_mField.get_moab()),mField(_mField),A(_A),F(_F) {
    }; 

    const double *G_TET_W,*G_TRI_W;
    vector<double> g_NTET,g_NTRI;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_TET_W = G_TET_W45;
      g_NTRI.resize(3*13);
      ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 
      G_TRI_W = G_TRI_W13;
      ierr = VecSetOption(F,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);  CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    vector<vector<DofIdx> > RowGlobDofs,ColGlobDofs;
    vector<vector< ublas::matrix<FieldData> > > diffRowNMatrix,diffColNMatrix;

    vector< ublas::matrix<FieldData> > invH;
    vector< FieldData > detH;

    PetscErrorCode get_ShapeFunctionsAndIndices() {
      PetscFunctionBegin;

      //Higher order approximation of geometry
      ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);
      //
      RowGlobDofs.resize(1+6+4+1);
      diffRowNMatrix.resize(1+6+4+1);
      ierr = GetRowGlobalIndices("POTENTIAL_FIELD",RowGlobDofs[0]); CHKERRQ(ierr);
      ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",diffRowNMatrix[0]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffRowNMatrix[0]); CHKERRQ(ierr);
      //
      int ee = 0;
      for(;ee<6;ee++) {
	ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBEDGE,RowGlobDofs[1+ee],ee); CHKERRQ(ierr);
	if(RowGlobDofs[1+ee].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBEDGE,diffRowNMatrix[1+ee],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffRowNMatrix[1+ee]); CHKERRQ(ierr);
	  //
	} 
      }
      int ff = 0;
      for(;ff<4;ff++) {
	ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBTRI,RowGlobDofs[1+6+ff],ff); CHKERRQ(ierr);
	if(RowGlobDofs[1+6+ff].size()>0) {
	  ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBTRI,diffRowNMatrix[1+6+ff],ff); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffRowNMatrix[1+6+ff]); CHKERRQ(ierr);
	  //
	}
      }
      ierr = GetRowGlobalIndices("POTENTIAL_FIELD",MBTET,RowGlobDofs[1+6+4]); CHKERRQ(ierr);
      if(RowGlobDofs[1+6+4].size()>0) {
	ierr = GetGaussRowDiffNMatrix("POTENTIAL_FIELD",MBTET,diffRowNMatrix[1+6+4]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffRowNMatrix[1+6+4]); CHKERRQ(ierr);
	//
      }

      ColGlobDofs.resize(1+6+4+1);
      diffColNMatrix.resize(1+6+4+1);
      ierr = GetColGlobalIndices("POTENTIAL_FIELD",ColGlobDofs[0]); CHKERRQ(ierr);
      ierr = GetGaussColDiffNMatrix("POTENTIAL_FIELD",diffColNMatrix[0]); CHKERRQ(ierr);
      //HO gemometry
      ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffColNMatrix[0]); CHKERRQ(ierr);
      //
      ee = 0;
      for(;ee<6;ee++) {
	ierr = GetColGlobalIndices("POTENTIAL_FIELD",MBEDGE,ColGlobDofs[1+ee],ee); CHKERRQ(ierr);
	if(ColGlobDofs[1+ee].size()>0) {
	  ierr = GetGaussColDiffNMatrix("POTENTIAL_FIELD",MBEDGE,diffColNMatrix[1+ee],ee); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffColNMatrix[1+ee]); CHKERRQ(ierr);
	  //
	} 
      }
      ff = 0;
      for(;ff<4;ff++) {
	ierr = GetColGlobalIndices("POTENTIAL_FIELD",MBTRI,ColGlobDofs[1+6+ff],ff); CHKERRQ(ierr);
	if(ColGlobDofs[1+6+ff].size()>0) {
	  ierr = GetGaussColDiffNMatrix("POTENTIAL_FIELD",MBTRI,diffColNMatrix[1+6+ff],ff); CHKERRQ(ierr);
	  //HO gemometry
	  ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffColNMatrix[1+6+ff]); CHKERRQ(ierr);
	  //
	}
      }
      ierr = GetColGlobalIndices("POTENTIAL_FIELD",MBTET,ColGlobDofs[1+6+4]); CHKERRQ(ierr);
      if(ColGlobDofs[1+6+4].size()>0) {
	ierr = GetGaussColDiffNMatrix("POTENTIAL_FIELD",MBTET,diffColNMatrix[1+6+4]); CHKERRQ(ierr);
	//HO gemometry
	ierr = GetHierarchicalGeometryApproximation_ApplyToDiffShapeFunction(1,invH,diffColNMatrix[1+6+4]); CHKERRQ(ierr);
	//
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode compute_LHS(FieldData a = 1) {
      PetscFunctionBegin;

      ublas::matrix<FieldData> K;
      for(int rr = 0;rr<(1+6+4+1); rr++) {
	if(RowGlobDofs[rr].size()==0) continue;
	for(int cc = 0;cc<(1+6+4+1);cc++) {
	  if(RowGlobDofs[cc].size()==0) continue;
	  K.resize(RowGlobDofs[rr].size(),ColGlobDofs[cc].size());
	  for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) {
	    ublas::noalias(K) = prod( trans(diffRowNMatrix[rr][gg]),diffColNMatrix[cc][gg] ); 
	    if(detH.size()==0) {
	      K *= a*V*G_TET_W[gg];
	    } else {
	      K *= a*V*detH[gg]*G_TET_W[gg];
	    }
	    ierr = MatSetValues(A,
	      RowGlobDofs[rr].size(),&(RowGlobDofs[rr])[0],
	      ColGlobDofs[cc].size(),&(ColGlobDofs[cc])[0],
	      &(K.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	}
      }

      PetscFunctionReturn(0);
    }

    ublas::matrix<double> gaussPts;
    PetscErrorCode Get_g_NTET() {
      PetscFunctionBegin;

      int order = 1;
      for(_IT_GET_FEDATA_DOFS_FOR_LOOP_(this,"POTENTIAL_FIELD",dof)) {
	order = max(order,dof->get_max_order());
      }

      int rule = max(0,order-1);
      if( 2*rule + 1 < 2*(order-1) ) {
	SETERRQ2(PETSC_COMM_SELF,1,"wrong rule %d %d",order,rule);
      }
      int nb_gauss_pts = gm_rule_size(rule,3);
      if(gaussPts.size2() == (unsigned int)nb_gauss_pts) {
	PetscFunctionReturn(0);
      }
      gaussPts.resize(4,nb_gauss_pts);
      ierr = Grundmann_Moeller_integration_points_3D_TET(
	rule,&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),&gaussPts(3,0)); CHKERRQ(ierr);

      g_NTET.resize(4*nb_gauss_pts);
      ierr = ShapeMBTET(&g_NTET[0],&gaussPts(0,0),&gaussPts(1,0),&gaussPts(2,0),nb_gauss_pts); CHKERRQ(ierr);
      G_TET_W = &gaussPts(3,0);

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = Get_g_NTET(); CHKERRQ(ierr);
      ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

      ierr = get_ShapeFunctionsAndIndices(); CHKERRQ(ierr);
      ierr = compute_LHS(); CHKERRQ(ierr);

      ierr = OpStudentEnd(); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  };


struct PostProcPotentialFlowOnRefMesh: public PostProcDisplacemenysAndStarinOnRefMesh {

  Tag th_phi,th_u;
  PostProcPotentialFlowOnRefMesh(Interface& _moab): PostProcDisplacemenysAndStarinOnRefMesh(_moab,"DISPLACEMENT") {
    double def_VAL2[3] = { 0.0, 0.0, 0.0 };
    rval = moab_post_proc.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
    rval = moab_post_proc.tag_get_handle("U",3,MB_TYPE_DOUBLE,th_u,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL2); CHKERR_THROW(rval);
  }

  PetscErrorCode do_operator() {
    PetscFunctionBegin;
    ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);

    Range ref_nodes;
    rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
    if(4*ref_nodes.size()!=g_NTET.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    if(ref_nodes.size()!=coords_at_Gauss_nodes.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
    Range::iterator nit = ref_nodes.begin();
    node_map.clear();
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

    PetscFunctionReturn(0);
  }

  vector< ublas::matrix<FieldData> > invH;
  vector< FieldData > detH;

  PetscErrorCode operator()() {
    PetscFunctionBegin;

    //Loop over elements
    ierr = do_operator(); CHKERRQ(ierr);

    //Higher order approximation of geometry
    ierr = GetHierarchicalGeometryApproximation(invH,detH); CHKERRQ(ierr);

    //Strains to Nodes in PostProc Mesh: create vector containing matrices
    vector< ublas::matrix< FieldData > > negativeVelocities;
    vector< ublas::vector< FieldData > > phi;

    ierr = GetGaussDataVector("POTENTIAL_FIELD",phi); CHKERRQ(ierr);
    ierr = GetGaussDiffDataVector("POTENTIAL_FIELD",negativeVelocities); CHKERRQ(ierr);

    map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
    for(;mit!=node_map.end();mit++) {
      int gg = distance(node_map.begin(),mit);
      negativeVelocities[gg] = -prod( trans( invH[gg] ), trans(negativeVelocities[gg]) );
      rval = moab_post_proc.tag_set_data(th_u,&mit->second,1,&(negativeVelocities[gg].data()[0])); CHKERR_PETSC(rval);
      rval = moab_post_proc.tag_set_data(th_phi,&mit->second,1,&(phi[gg].data()[0])); CHKERR_PETSC(rval);
    }

    PetscFunctionReturn(0);

  }

};

};

#endif // __POTENTIALFLOWFEMETHOD_HPP__

