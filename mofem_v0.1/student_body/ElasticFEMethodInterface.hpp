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


#ifndef __ELASTICFEMETHODFORINTERFACE_HPP__
#define __ELASTICFEMETHODFORINTERFACE_HPP__


#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>
extern "C" {
#include <gm_rule.h>
}

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

namespace MoFEM {

struct ToolsInterfaceFEMethod {

  double youngModulus;
  ToolsInterfaceFEMethod(double _young_modulus): youngModulus(_young_modulus) {}

  ublas::matrix<double> R;
  ublas::matrix<double> Dglob;
  double tangent1[3],tangent2[3];

  PetscErrorCode CalcR(
    const double *diffNTRI,
    const double *coords_face3,
    const double *normal3,
    const double area3) {
    PetscFunctionBegin;
    bzero(tangent1,3*sizeof(double));
    bzero(tangent2,3*sizeof(double));
    int ii = 0;
    for(; ii<3; ii++) {
	tangent1[0] += coords_face3[3*ii + 0]*diffNTRI[2*ii+0];
	tangent1[1] += coords_face3[3*ii + 1]*diffNTRI[2*ii+0];
	tangent1[2] += coords_face3[3*ii + 2]*diffNTRI[2*ii+0];
	tangent2[0] += coords_face3[3*ii + 0]*diffNTRI[2*ii+1];
	tangent2[1] += coords_face3[3*ii + 1]*diffNTRI[2*ii+1];
	tangent2[2] += coords_face3[3*ii + 2]*diffNTRI[2*ii+1];
    }
    R = ublas::zero_matrix<double>(3,3);
    ublas::matrix_row<ublas::matrix<double> > R_normal(R,0);
    R_normal[0] = normal3[0];
    R_normal[1] = normal3[1];
    R_normal[2] = normal3[2];
    R_normal /= 2.*area3;
    ublas::matrix_row<ublas::matrix<double> > R_tangent1(R,1);
    R_tangent1[0] = tangent1[0];
    R_tangent1[1] = tangent1[1];
    R_tangent1[2] = tangent1[2];
    double nrm1 = cblas_dnrm2(3,tangent1,1);
    R_tangent1 /= nrm1;
    ublas::matrix_row<ublas::matrix<double> > R_tangent2(R,2);
    R_tangent2[0] = tangent2[0];
    R_tangent2[1] = tangent2[1];
    R_tangent2[2] = tangent2[2];
    double nrm2 = cblas_dnrm2(3,tangent2,1);
    R_tangent2 /= nrm2;
    //cerr << R << endl;
    PetscFunctionReturn(0);
  }

  PetscErrorCode CalcDglob() {
    PetscFunctionBegin;
    ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
    int ii = 0;
    for(;ii<3;ii++) Dloc(ii,ii) = youngModulus;
    //Dloc(0,0) = youngModulus;
    Dglob = prod( Dloc, R );
    Dglob = prod( trans(R), Dglob );
    PetscFunctionReturn(0);
  }



};

struct InterfaceFEMethod: public FEMethod_UpLevelStudent,ToolsInterfaceFEMethod {

  FieldInterface &mField;
  vector<ublas::vector<FieldData> > DispData;
  string field_name;

  InterfaceFEMethod(FieldInterface& _mField,Mat &_Aij,Vec _X,Vec _F,double _young_modulus,string _field_name = "DISPLACEMENT"):
    FEMethod_UpLevelStudent(_mField.get_moab(),1),ToolsInterfaceFEMethod(_young_modulus),
    mField(_mField),field_name(_field_name) {

    snes_B = &_Aij;
    snes_x = _X;
    snes_f = _F;

    DispData.resize(1+6+2);
  }

  PetscErrorCode CalcR() {
    PetscFunctionBegin;
    ierr = ToolsInterfaceFEMethod::CalcR(diffNTRI,coords_face3,normal3,area3); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

  vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap;
  vector<ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > > gap_loc;


  /* \brief calulate gap in global and local coorinates
    *
    * Function, make a loop for all gauss points, and calculate gap ( separation
    * of interface ). We have three types of shape functions, for nodes, edges and
    * face of interface itself.  Values of shepe functions, for each gauss pt, are
    * stored in matrixes, nodeNTRI, _H1edgeN_, _H1edgeN_, _H1faceN_, for nodes,
    * edges and faces, respectively. 
    *
  */
  virtual PetscErrorCode Calc_gap() {
    PetscFunctionBegin;
    SideNumber_multiIndex &side_table = fe_ptr->get_side_number_table();
    SideNumber_multiIndex::iterator sit = side_table.begin();
    map<EntityType,bitset<9> > ents_bits;
    for(;sit!=side_table.end();sit++) {
      if(sit->get_ent_type() == MBVERTEX) {
	ents_bits[MBVERTEX].set(sit->side_number);
      }	else if(sit->get_ent_type() == MBEDGE) {
	ents_bits[MBEDGE].set(sit->side_number);
      }
    }
    int g_dim = g_NTRI.size()/3;
    gap.resize(g_dim);
    gap_loc.resize(g_dim);
    for(int gg = 0;gg<g_dim;gg++) {
	gap[gg] = ublas::zero_vector<FieldData>(3);
	//nodes
	double *nodeNTRI = &g_NTRI[gg*3];
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX,0));
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  int side = dit->side_number_ptr->side_number;
	  if(ents_bits[MBVERTEX].test(side) && ents_bits[MBVERTEX].test(side+3)) {
	    FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator diit;
	    (gap[gg])[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	  }
	}
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBVERTEX,3));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  int side = dit->side_number_ptr->side_number;
	  if(ents_bits[MBVERTEX].test(side-3) && ents_bits[MBVERTEX].test(side)) {
	    (gap[gg])[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	  }
	}
	//edges
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBEDGE,2));
	for(;dit!=hi_dit;dit++) {
	  int side = dit->side_number_ptr->side_number;	
    	  if(ents_bits[MBEDGE].test(side) && ents_bits[MBEDGE].test(side+6)) {
	    assert(side >= 0);
	    assert(side <= 2);
	    int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[side]);
	    int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	    double *_H1edgeN_ = &*H1edgeN[side].begin();
	    double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	    (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  }
	} 
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  int side = dit->side_number_ptr->side_number;	
    	  if(ents_bits[MBEDGE].test(side-6) && ents_bits[MBEDGE].test(side)) {
	    double *_H1edgeN_ = &H1edgeN[side][0];
	    int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[side]);
	    int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	    double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	    (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  }
	} 
	//faces
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple(field_name,MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple(field_name,MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  (gap[gg])[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
	gap_loc[gg] = prod(R,gap[gg]);
    }
    PetscFunctionReturn(0);
  }

  int row_mat,col_mat;
  ublas::matrix<ublas::matrix<FieldData> > K;
  vector<vector<ublas::matrix<FieldData> > > rowNMatrices;
  vector<vector<ublas::matrix<FieldData> > > colNMatrices;
  vector<vector<DofIdx> > RowGlob;
  vector<vector<DofIdx> > ColGlob;

  vector<DofIdx> DirihletBC;

  virtual PetscErrorCode LhsInt() {
    PetscFunctionBegin;
		
    int g_dim = g_NTRI.size()/3;
    K.resize(row_mat,col_mat);
    for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	    ///K matrices
	    if(gg == 0) {
	      K(rr,cc) = ublas::zero_matrix<FieldData>(row_Mat.size2(),col_Mat.size2());
	    }
	    double w = area4*G_TRI_W[gg];
	    ublas::matrix<FieldData> NTD = prod( trans(row_Mat), w*Dglob );
	    K(rr,cc) += prod(NTD , col_Mat ); 
	  }
	}
	if(RowGlob[rr].size()==0) continue;
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(*snes_B,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode RhsInt() {
    PetscFunctionBegin;

    int g_dim = g_NTRI.size()/3;
    for(int gg = 0;gg<g_dim;gg++) {
      //Traction
      ublas::vector<FieldData,ublas::bounded_array<FieldData, 3> > traction;
      traction = prod(Dglob,gap[gg]);
      if(traction.size()!=3) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      double w = area3*G_TRI_W[gg];
      for(int rr = 0;rr<row_mat;rr++) {
	ublas::matrix<FieldData> &N = (rowNMatrices[rr])[gg];
	ublas::vector<FieldData> f_int = prod(trans(N),w*traction);
	if(RowGlob[rr].size()==0) continue;
	if(RowGlob[rr].size()!=f_int.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	ierr = VecSetValues(snes_f,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int.data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode Matrices() {
    PetscFunctionBegin;
    //rows
    RowGlob.resize(1+6+2);
    rowNMatrices.resize(1+6+2);
    row_mat = 0;
    ierr = GetRowGlobalIndices(field_name,RowGlob[row_mat]); CHKERRQ(ierr);
    ierr = GetGaussRowNMatrix(field_name,rowNMatrices[row_mat]); CHKERRQ(ierr);
    row_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowGlobalIndices(field_name,MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix(field_name,MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowGlobalIndices(field_name,MBEDGE,RowGlob[row_mat],ee+6); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix(field_name,MBEDGE,rowNMatrices[row_mat],ee+6); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    ierr = GetRowGlobalIndices(field_name,MBTRI,RowGlob[row_mat],3); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix(field_name,MBTRI,rowNMatrices[row_mat],3); CHKERRQ(ierr);
	row_mat++;
    }
    ierr = GetRowGlobalIndices(field_name,MBTRI,RowGlob[row_mat],4); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix(field_name,MBTRI,rowNMatrices[row_mat],4); CHKERRQ(ierr);
	row_mat++;
    }
    //cols
    ColGlob.resize(1+6+2);
    colNMatrices.resize(1+6+2);
    col_mat = 0;
    ierr = GetColGlobalIndices(field_name,ColGlob[col_mat]); CHKERRQ(ierr);
    ierr = GetGaussColNMatrix(field_name,colNMatrices[col_mat]); CHKERRQ(ierr);
    ierr = GetDataVector(field_name,DispData[col_mat]); CHKERRQ(ierr);
    col_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColGlobalIndices(field_name,MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix(field_name,MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetDataVector(field_name,MBEDGE,DispData[col_mat],ee); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColGlobalIndices(field_name,MBEDGE,ColGlob[col_mat],ee+6); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix(field_name,MBEDGE,colNMatrices[col_mat],ee+6); CHKERRQ(ierr);
	  ierr = GetDataVector(field_name,MBEDGE,DispData[col_mat],ee+6); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    ierr = GetColGlobalIndices(field_name,MBTRI,ColGlob[col_mat],3); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix(field_name,MBTRI,colNMatrices[col_mat],3); CHKERRQ(ierr);
	ierr = GetDataVector(field_name,MBTRI,DispData[col_mat],3); CHKERRQ(ierr);
	col_mat++;
    }
    ierr = GetColGlobalIndices(field_name,MBTRI,ColGlob[col_mat],4); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix(field_name,MBTRI,colNMatrices[col_mat],4); CHKERRQ(ierr);
	ierr = GetDataVector(field_name,MBTRI,DispData[col_mat],4); CHKERRQ(ierr);
	col_mat++;
    }
    PetscFunctionReturn(0);
  }

  vector<double> g_NTRI;
  const double *G_TRI_W;

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    g_NTRI.resize(3*28);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X28,G_TRI_Y28,28); 
    G_TRI_W = G_TRI_W28;

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    switch(snes_ctx) {
      case ctx_SNESNone: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: {
	ierr = MatAssemblyBegin(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(*snes_B,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetFunction:  {
	ierr = VecAssemblyBegin(snes_f); CHKERRQ(ierr);
	ierr = VecAssemblyEnd(snes_f); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }
    PetscFunctionReturn(0);
  }

  int iter;
  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

    //Rotation matrix
    ierr = CalcR(); CHKERRQ(ierr);
    ierr = CalcDglob(); CHKERRQ(ierr);

    //Calculate Matrices
    ierr = Matrices();    CHKERRQ(ierr);
    //Calcualte gap
    ierr = Calc_gap(); CHKERRQ(ierr);

    switch(snes_ctx) {
      case ctx_SNESNone: {
	ierr = RhsInt(); CHKERRQ(ierr);
	ierr = LhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetFunction: { 
	ierr = RhsInt(); CHKERRQ(ierr);
      }
      break;
      case ctx_SNESSetJacobian: { 
	ierr = LhsInt(); CHKERRQ(ierr);
      }
      break;
      default:
	SETERRQ(PETSC_COMM_SELF,1,"not implemented");
    }

    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }

};

struct PostProcCohesiveForces:public FEMethod_UpLevelStudent,PostProcOnRefMesh_Base,ToolsInterfaceFEMethod {
  
    FieldInterface &mField;
    ParallelComm* pcomm;

    PostProcCohesiveForces(FieldInterface& _mField,double _young_modulus): 
      FEMethod_UpLevelStudent(_mField.get_moab()), PostProcOnRefMesh_Base(),ToolsInterfaceFEMethod(_young_modulus),
      mField(_mField) {
      pcomm = ParallelComm::get_pcomm(&mField.get_moab(),MYPCOMM_INDEX);
    };

    vector<double> g_NTRI;
    vector<EntityHandle> nodes_on_face3;
    vector<EntityHandle> nodes_on_face4;
  
    Tag th_cohesive_force;
    Tag th_gap;
    Tag th_disp;

    PetscErrorCode preProcess() {
      PetscFunctionBegin;

      if(init_ref) PetscFunctionReturn(0);
      
      double def_VAL[3] = {0,0,0};
      // create TAG
      rval = moab_post_proc.tag_get_handle("COHESIVE_FORCE_VAL",3,MB_TYPE_DOUBLE,th_cohesive_force,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
      rval = moab_post_proc.tag_get_handle("GAP_VAL",3,MB_TYPE_DOUBLE,th_gap,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);
      rval = moab_post_proc.tag_get_handle("DISPLACEMENTS_VAL",3,MB_TYPE_DOUBLE,th_disp,MB_TAG_CREAT|MB_TAG_SPARSE,&def_VAL); CHKERR(rval);

      double base_coords[] = {
	//face3
	0,0,0,
	1,0,0,
	0,1,0,
	//face4
	0,0,1,
	1,0,1,
	0,1,1 
      };
      EntityHandle nodes[6];
      for(int nn = 0;nn<6;nn++) {
	rval = moab_ref.create_vertex(&base_coords[3*nn],nodes[nn]); CHKERR_PETSC(rval);
      }
      EntityHandle prism;
      rval = moab_ref.create_element(MBPRISM,nodes,6,prism); CHKERR_PETSC(rval);
      //create faces using get_adjacencies
      Range ref_faces;
      rval = moab_ref.get_adjacencies(&nodes[0],3,2,true,ref_faces,Interface::UNION); CHKERR_PETSC(rval);
      Range faces;
      rval = moab_ref.get_adjacencies(&nodes[3],3,2,true,ref_faces,Interface::UNION); CHKERR_PETSC(rval);
  
      //
      FieldCore core_ref(moab_ref);
      FieldInterface& mField_ref = core_ref;
      ierr = mField_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

      for(int ll = 0;ll<max_level;ll++) {
	PetscPrintf(PETSC_COMM_WORLD,"Refine Level %d\n",ll);
	rval = moab_ref.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = mField_ref.get_entities_by_ref_level(BitRefLevel().set(ll),BitRefLevel().set(),meshset_level[ll]); CHKERRQ(ierr);
	ierr = mField_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = mField_ref.refine_PRISM(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET|MESHSET_TRACK_OWNER,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = mField_ref.get_entities_by_ref_level(BitRefLevel().set(max_level),BitRefLevel().set(),meshset_level[max_level]); CHKERRQ(ierr);

      //if(pcomm->rank()==0) {
	//moab_ref.write_file("debug.vtk","VTK",""); CHKERR_PETSC(rval);
      //}

      //
      Range ref_nodes;
      rval = moab_ref.get_entities_by_type(0,MBVERTEX,ref_nodes); CHKERR_PETSC(rval);
      std::vector<double> ref_coords(3*ref_nodes.size());
      rval = moab_ref.get_coords(ref_nodes,&ref_coords[0]); CHKERR_PETSC(rval);
      std::vector<double> ref_coordsf3(2*ref_nodes.size()/2);
      std::vector<double> ref_coordsf4(2*ref_nodes.size()/2);
      int nnn = 0,mmm = 0;
      for(unsigned int nn = 0;nn<ref_nodes.size();nn++) {
	if(ref_coords[3*nn+2]>0) {
	  ref_coordsf4[mmm] = ref_coords[3*nn+0];
	  nodes_on_face4.push_back(ref_nodes[nn]);
	  mmm++;
	} else {
	  ref_coordsf3[nnn] = ref_coords[3*nn+0];
	  nodes_on_face3.push_back(ref_nodes[nn]);
	  nnn++;
	}
      }
      assert(nnn == mmm);
      nnn = 0; mmm = 0;
      for(unsigned int nn = 0;nn<ref_nodes.size();nn++) {
	if(ref_coords[3*nn+2]>0) {
	  ref_coordsf4[nodes_on_face3.size()+mmm] = ref_coords[3*nn+1];
	  mmm++;
	} else {
	  ref_coordsf3[nodes_on_face3.size()+nnn] = ref_coords[3*nn+1];
	  nnn++;
	}
      }
      assert(nnn == mmm);

      assert(nodes_on_face3.size() == nodes_on_face4.size());
      g_NTRI.resize(6*nodes_on_face3.size());
      ShapeMBTRI(&g_NTRI[0],&ref_coordsf3[0],&ref_coordsf3[nodes_on_face3.size()],nodes_on_face3.size());
      ShapeMBTRI(&g_NTRI[3*nodes_on_face3.size()],&ref_coordsf4[0],&ref_coordsf4[nodes_on_face4.size()],nodes_on_face4.size());

      //cerr << ref_coords_z0.size() << " " << ref_nodes.size() << " " << g_NTRI.size() << endl;

      init_ref = true;

      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;
      ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

      //cerr << nodes_on_face3.size() << " " << g_NTRI.size() << " " << coords_at_Gauss_nodes.size() << endl;
      if(6*nodes_on_face3.size()!=g_NTRI.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");

      coords_at_Gauss_nodes.resize(2*nodes_on_face3.size());
      for(unsigned int gg = 0;gg<nodes_on_face3.size();gg++) {
	coords_at_Gauss_nodes[gg].resize(3);
	coords_at_Gauss_nodes[nodes_on_face3.size()+gg].resize(3);
	for(int dd = 0;dd<3;dd++) {
	  (coords_at_Gauss_nodes[gg])[dd] = cblas_ddot(3,&coords_face3[dd],3,&g_NTRI[gg*3],1);
	  (coords_at_Gauss_nodes[nodes_on_face3.size()+gg])[dd] = cblas_ddot(3,&coords_face4[dd],3,&g_NTRI[3*nodes_on_face3.size()+gg*3],1);
	}
      }

      map<EntityHandle,EntityHandle> node_map;
      vector<EntityHandle>::iterator nit = nodes_on_face3.begin();
      for(int nn = 0;nit!=nodes_on_face3.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nn]).data()[0],node); CHKERR_PETSC(rval);
      }
      nit = nodes_on_face4.begin();
      for(int nn = 0;nit!=nodes_on_face4.end();nit++,nn++) {
	EntityHandle &node = node_map[*nit];
	rval = moab_post_proc.create_vertex(&(coords_at_Gauss_nodes[nodes_on_face3.size()+nn]).data()[0],node); CHKERR_PETSC(rval);
      }

      Range ref_prisms;
      rval = moab_ref.get_entities_by_type(meshset_level[max_level],MBPRISM,ref_prisms); CHKERR_PETSC(rval);
      Range::iterator pit = ref_prisms.begin();
      for(;pit!=ref_prisms.end();pit++) {
	const EntityHandle *conn_ref;
        int num_nodes;
	rval = moab_ref.get_connectivity(*pit,conn_ref,num_nodes,true); CHKERR_PETSC(rval);
	assert(num_nodes==6);
	EntityHandle conn_post_proc[num_nodes];
	for(int nn = 0;nn<num_nodes;nn++) {
	  map<EntityHandle,EntityHandle>::iterator mit = node_map.find(conn_ref[nn]);
	  assert(mit!=node_map.end());
	  NOT_USED(mit);
	  conn_post_proc[nn] = node_map[conn_ref[nn]];
	}
	EntityHandle ref_prism;
	rval = moab_post_proc.create_element(MBPRISM,conn_post_proc,6,ref_prism); CHKERR_PETSC(rval);
      }

      //Rotation matrix
      ierr = CalcR(diffNTRI,coords_face3,normal3,area3); CHKERRQ(ierr);

      //Dglob
      ierr = CalcDglob(); CHKERRQ(ierr);

      //cerr << Dglob << endl;

      //face3
      for(unsigned int gg = 0;gg<nodes_on_face3.size();gg++) {
	double *nodeNTRI = &g_NTRI[gg*3];
	EntityHandle node = node_map[nodes_on_face3[gg]];
	//gap
	double gap[] = {0,0,0};
	rval = moab_post_proc.tag_set_data(th_gap,&node,1,gap); CHKERR_PETSC(rval);
	double *gap_ptr;
	rval = moab_post_proc.tag_get_by_ptr(th_gap,&node,1,(const void **)&gap_ptr); CHKERR_PETSC(rval);
	//disp
	double disp[] = {0,0,0};
	rval = moab_post_proc.tag_set_data(th_disp,&node,1,disp); CHKERR_PETSC(rval);
	double *disp_ptr;
	rval = moab_post_proc.tag_get_by_ptr(th_disp,&node,1,(const void **)&disp_ptr); CHKERR_PETSC(rval);
	//nodes
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,0));
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  //cerr << *dit << endl;
	  disp_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	  gap_ptr[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	}
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,3));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  gap_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	}
	//edges
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,2));
	for(;dit!=hi_dit;dit++) {
	  int side_number = dit->side_number_ptr->side_number;	
	  assert(side_number >= 0);
	  assert(side_number <= 2);
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  /*cerr << "side_number " << side_number 
	    << " " << dit->get_EntDofIdx() 
	    << " " << approx_dof 
	    << " " << nb_dofs_H1edge 
	    << " " << H1edgeN[side_number].size()
	    << endl;*/
	  double *_H1edgeN_ = &*H1edgeN[side_number].begin();
	  double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	  disp_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	} 
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	} 
	//faces
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	} 
	
	ublas::vector<double,ublas::bounded_array<double, 3> > Gap(3);
	Gap[0] = gap_ptr[0];
 	Gap[1] = gap_ptr[1];
 	Gap[2] = gap_ptr[2];  
	ublas::vector<double,ublas::bounded_array<double, 3> > Trac(3);
	Trac = prod(Dglob,Gap);
	rval = moab_post_proc.tag_set_data(th_cohesive_force,&node,1,&Trac.data()[0]); CHKERR_PETSC(rval);

      }

      //face4
      for(unsigned int gg = 0;gg<nodes_on_face4.size();gg++) {
	double *nodeNTRI = &g_NTRI[3*nodes_on_face3.size()+gg*3];
	EntityHandle node = node_map[nodes_on_face4[gg]];
	//node
	double disp[] = {0,0,0};
	rval = moab_post_proc.tag_set_data(th_disp,&node,1,disp); CHKERR_PETSC(rval);
	double *disp_ptr;
	rval = moab_post_proc.tag_get_by_ptr(th_disp,&node,1,(const void **)&disp_ptr); CHKERR_PETSC(rval);
	//gap
	double gap[] = {0,0,0};
	rval = moab_post_proc.tag_set_data(th_gap,&node,1,gap); CHKERR_PETSC(rval);
	double *gap_ptr;
	rval = moab_post_proc.tag_get_by_ptr(th_gap,&node,1,(const void **)&gap_ptr); CHKERR_PETSC(rval);
	//nodes
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,3));
	FEDofMoFEMEntity_multiIndex::index<Composite_Name_Type_And_Side_Number_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  //cerr << *dit << " " << nodeNTRI[dit->side_number_ptr->side_number] << " " << dit->get_FieldData() << endl;
	  disp_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	  gap_ptr[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	}
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,0));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  gap_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	}
	//edges
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[nodes_on_face3.size()*nb_dofs_H1edge + gg*nb_dofs_H1edge + approx_dof];
	  disp_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus*/
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); //*minus*/
	}
 	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,2));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[nodes_on_face3.size()*nb_dofs_H1edge + gg*nb_dofs_H1edge + approx_dof];
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); //*minus*/
	}
	//faces
	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[nodes_on_face3.size()*nb_dofs_H1face + gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus/
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_Name_Type_And_Side_Number_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[nodes_on_face3.size()*nb_dofs_H1face + gg*nb_dofs_H1face + approx_dof];
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}

	ublas::vector<double,ublas::bounded_array<double, 3> > Gap(3);
	Gap[0] = gap_ptr[0];
 	Gap[1] = gap_ptr[1];
 	Gap[2] = gap_ptr[2];  
	ublas::vector<double,ublas::bounded_array<double, 3> > Trac(3);
	Trac = prod(Dglob,Gap);
	rval = moab_post_proc.tag_set_data(th_cohesive_force,&node,1,&Trac.data()[0]); CHKERR_PETSC(rval);

      }

      
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      ParallelComm* pcomm_post_proc = ParallelComm::get_pcomm(&moab_post_proc,MYPCOMM_INDEX);
      if(pcomm_post_proc == NULL) pcomm_post_proc =  new ParallelComm(&moab_post_proc,PETSC_COMM_WORLD);
      for(unsigned int rr = 1; rr<pcomm_post_proc->size();rr++) {
	Range prisms;
	rval = moab_post_proc.get_entities_by_type(0,MBPRISM,prisms); CHKERR_PETSC(rval);
	rval = pcomm_post_proc->broadcast_entities(rr,prisms); CHKERR(rval);
      }
      PetscFunctionReturn(0);
    }


};

} 

#endif //__ELASTICFEMETHODFORINTERFACE_HPP__
