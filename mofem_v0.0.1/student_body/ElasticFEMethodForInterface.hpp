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
#include "moabFEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

namespace MoFEM {

struct InterfaceFEMethod: public ElasticFEMethod {

  double YoungModulus; 
  ublas::matrix<double> R;
  ublas::matrix<double> Dglob;
  double tangent1[3],tangent2[3];

  vector<ublas::vector<FieldData> > DispData;

  InterfaceFEMethod(
      FieldInterface& _mField,double _YoungModulus): 
      ElasticFEMethod(_mField),YoungModulus(_YoungModulus) {
      DispData.resize(1+6+2);
    };

  InterfaceFEMethod(
      FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,Mat &_Aij,Vec &_D,Vec& _F,double _YoungModulus): 
	ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,0,0),YoungModulus(_YoungModulus) {
    DispData.resize(1+6+2);
    };

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start Assembly\n",pcomm->rank(),v2-v1,t2-t1);

    ierr = PetscGetTime(&v1); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);
    g_NTET.resize(4*45);
    ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
    g_NTRI.resize(3*13);
    ShapeMBTRI(&g_NTRI[0],G_TRI_X13,G_TRI_Y13,13); 

    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    // Note MAT_FLUSH_ASSEMBLY
    ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = PetscGetTime(&v2); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode CalcR() {
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
    R_normal /= area3;
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

  virtual PetscErrorCode CalcDglob() {
    PetscFunctionBegin;
    ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
    int ii = 0;
    for(;ii<3;ii++) Dloc(ii,ii) = YoungModulus;
    //Dloc(0,0) = YoungModulus;
    Dglob = prod( Dloc, R );
    Dglob = prod( trans(R), Dglob );
    PetscFunctionReturn(0);
  }

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
	    double w = area3*G_TRI_W13[gg];
	    ublas::matrix<FieldData> NTD = prod( trans(row_Mat), w*Dglob );
	    K(rr,cc) += prod(NTD , col_Mat ); 
	  }
	}
	if(RowGlob[rr].size()==0) continue;
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K(rr,cc).size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K(rr,cc).size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K(rr,cc).data())[0],ADD_VALUES); CHKERRQ(ierr);
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
    ierr = GetRowGlobalIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
    ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
    row_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowGlobalIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowGlobalIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee+6); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee+6); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    ierr = GetRowGlobalIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],3); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],3); CHKERRQ(ierr);
	row_mat++;
    }
    ierr = GetRowGlobalIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],4); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],4); CHKERRQ(ierr);
	row_mat++;
    }
    //cols
    ColGlob.resize(1+6+2);
    colNMatrices.resize(1+6+2);
    col_mat = 0;
    ierr = GetColGlobalIndices("DISPLACEMENT",ColGlob[col_mat]); CHKERRQ(ierr);
    ierr = GetGaussColNMatrix("DISPLACEMENT",colNMatrices[col_mat]); CHKERRQ(ierr);
    ierr = GetDataVector("DISPLACEMENT",DispData[col_mat]); CHKERRQ(ierr);
    col_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColGlobalIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  ierr = GetDataVector("DISPLACEMENT",MBEDGE,DispData[col_mat],ee); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColGlobalIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee+6); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee+6); CHKERRQ(ierr);
	  ierr = GetDataVector("DISPLACEMENT",MBEDGE,DispData[col_mat],ee+6); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    ierr = GetColGlobalIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],3); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],3); CHKERRQ(ierr);
	ierr = GetDataVector("DISPLACEMENT",MBTRI,DispData[col_mat],3); CHKERRQ(ierr);
	col_mat++;
    }
    ierr = GetColGlobalIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],4); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],4); CHKERRQ(ierr);
	ierr = GetDataVector("DISPLACEMENT",MBTRI,DispData[col_mat],4); CHKERRQ(ierr);
	col_mat++;
    }
    PetscFunctionReturn(0);
  }

  virtual PetscErrorCode RhsAndLhs() {
    PetscFunctionBegin;

    //Rotation matrix
    ierr = CalcR(); CHKERRQ(ierr);
    //Dglob
    ierr = CalcDglob(); CHKERRQ(ierr);
    //Calculate Matrices
    ierr = Matrices();    CHKERRQ(ierr);
    //Apply Dirihlet BC
    ierr = dirihlet_bc_method_ptr->SetDirihletBC_to_ElementIndicies(this,RowGlob,ColGlob,DirihletBC); CHKERRQ(ierr);

    //Assemble interface
    ierr = LhsInt(); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

    ierr = RhsAndLhs(); CHKERRQ(ierr);

    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

struct PostProcCohesiveForces: public InterfaceFEMethod,PostProcOnRefMesh_Base {
  
    PostProcCohesiveForces(FieldInterface& _mField,double _YoungModulus): 
      InterfaceFEMethod(_mField,_YoungModulus), PostProcOnRefMesh_Base() {
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
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Start PostProc\n",pcomm->rank(),v2-v1,t2-t1);
      ierr = PetscGetTime(&v1); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

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
	rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = mField_ref.refine_get_ents(BitRefLevel().set(ll),BitRefLevel().set(),meshset_level[ll]); CHKERRQ(ierr);
	ierr = mField_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = mField_ref.refine_PRISM(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = mField_ref.refine_get_ents(BitRefLevel().set(max_level),BitRefLevel().set(),meshset_level[max_level]); CHKERRQ(ierr);

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
      ierr = CalcR(); CHKERRQ(ierr);

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
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,0));
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  //cerr << *dit << endl;
	  disp_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	  gap_ptr[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	}
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  gap_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	}
	//edges
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,2));
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
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[gg*nb_dofs_H1edge + approx_dof];
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	} 
	//faces
	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
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
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,3));
	FEDofMoFEMEntity_multiIndex::index<Composite_mi_tag>::type::iterator 
	  hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,5));
	for(;dit!=hi_dit;dit++) {
	  //cerr << *dit << " " << nodeNTRI[dit->side_number_ptr->side_number] << " " << dit->get_FieldData() << endl;
	  disp_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	  gap_ptr[dit->get_dof_rank()] -= nodeNTRI[dit->side_number_ptr->side_number-3]*dit->get_FieldData();
	}
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,0));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBVERTEX,2));
	for(;dit!=hi_dit;dit++) {
	  gap_ptr[dit->get_dof_rank()] += nodeNTRI[dit->side_number_ptr->side_number]*dit->get_FieldData();
	}
	//edges
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,6));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,8));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[nodes_on_face3.size()*nb_dofs_H1edge + gg*nb_dofs_H1edge + approx_dof];
	  disp_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus*/
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); //*minus*/
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,0));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBEDGE,2));
	for(;dit!=hi_dit;dit++) {
	  double *_H1edgeN_ = &H1edgeN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1edge = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1edgeN_[nodes_on_face3.size()*nb_dofs_H1edge + gg*nb_dofs_H1edge + approx_dof];
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); //*minus*/
	}
	//faces
	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderFaceH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[nodes_on_face3.size()*nb_dofs_H1face + gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus/
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().lower_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
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
      ierr = PetscGetTime(&v2); CHKERRQ(ierr);
      ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);
      ParallelComm* pcomm_post_proc = ParallelComm::get_pcomm(&moab_post_proc,MYPCOMM_INDEX);
      if(pcomm_post_proc == NULL) pcomm_post_proc =  new ParallelComm(&moab_post_proc,PETSC_COMM_WORLD);
      for(unsigned int rr = 1; rr<pcomm_post_proc->size();rr++) {
	Range prisms;
	rval = moab_post_proc.get_entities_by_type(0,MBPRISM,prisms); CHKERR_PETSC(rval);
	rval = pcomm_post_proc->broadcast_entities(rr,prisms); CHKERR(rval);
      }
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End PostProc: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
      PetscFunctionReturn(0);
    }


};

} 

#endif //__ELASTICFEMETHODFORINTERFACE_HPP__
