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

using namespace MoFEM;

struct MyElasticFEMethod: public ElasticFEMethod {

  Range& SideSet3;

  MyElasticFEMethod(Interface& _moab): ElasticFEMethod(_moab),SideSet3(dummy) {};

  MyElasticFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3): 
      ElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,_SideSet1,_SideSet2), SideSet3(_SideSet3) {};

  PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction2(3);
      traction2[0] = 0;
      traction2[1] = +1;
      traction2[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(traction2,SideSet2); CHKERRQ(ierr);

      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction3(3);
      traction3[0] = 0;
      traction3[1] = -1;
      traction3[2] = 0;
      ierr = ElasticFEMethod::NeumannBC(traction3,SideSet3); CHKERRQ(ierr);

      PetscFunctionReturn(0);
  }


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

struct InterfaceFEMethod: public MyElasticFEMethod {

  double YoungModulus; 
  ublas::matrix<double> R;
  ublas::matrix<double> Dglob;
  double tangent1[3],tangent2[3];

  InterfaceFEMethod(
      Interface& _moab,double _YoungModulus): 
      MyElasticFEMethod(_moab),YoungModulus(_YoungModulus) {};

  InterfaceFEMethod(
      Interface& _moab,Mat &_Aij,Vec& _F,double _YoungModulus,Range &_SideSet1,Range &_SideSet2,Range &_SideSet3): 
	MyElasticFEMethod(_moab,_Aij,_F,0,0,_SideSet1,_SideSet2,_SideSet3),YoungModulus(_YoungModulus) {};

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
    ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
    // Note MAT_FLUSH_ASSEMBLY
    ierr = MatAssemblyBegin(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(Aij,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
    ierr = PetscGetTime(&v2); CHKERRQ(ierr);
    ierr = PetscGetCPUTime(&t2); CHKERRQ(ierr);

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"End Assembly: Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
    PetscFunctionReturn(0);
  }

  PetscErrorCode CalcR() {
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

  PetscErrorCode CalcDglob() {
    PetscFunctionBegin;
    ublas::matrix<double> Dloc = ublas::zero_matrix<double>(3,3);
    int ii = 0;
    for(;ii<3;ii++) Dloc(ii,ii) = YoungModulus;
    Dglob = prod( Dloc, R );
    Dglob = prod( trans(R), Dglob );
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    ierr = OpStudentStart_PRISM(g_NTRI); CHKERRQ(ierr);

    //Rotation matrix
    ierr = CalcR(); CHKERRQ(ierr);

    //Dglob
    ierr = CalcDglob(); CHKERRQ(ierr);

    //rows
    RowGlob.resize(1+6+2);
    rowNMatrices.resize(1+6+2);
    row_mat = 0;
    ierr = GetRowIndices("DISPLACEMENT",RowGlob[row_mat]); CHKERRQ(ierr);
    ierr = GetGaussRowNMatrix("DISPLACEMENT",rowNMatrices[row_mat]); CHKERRQ(ierr);
    row_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetRowIndices("DISPLACEMENT",MBEDGE,RowGlob[row_mat],ee+6); CHKERRQ(ierr);
	if(RowGlob[row_mat].size()!=0) {
	  ierr = GetGaussRowNMatrix("DISPLACEMENT",MBEDGE,rowNMatrices[row_mat],ee+6); CHKERRQ(ierr);
	  row_mat++;
	}
    }
    ierr = GetRowIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],3); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],3); CHKERRQ(ierr);
	row_mat++;
    }
    ierr = GetRowIndices("DISPLACEMENT",MBTRI,RowGlob[row_mat],4); CHKERRQ(ierr);
    if(RowGlob[row_mat].size()!=0) {
	ierr = GetGaussRowNMatrix("DISPLACEMENT",MBTRI,rowNMatrices[row_mat],4); CHKERRQ(ierr);
	row_mat++;
    }
    //cols
    ColGlob.resize(1+6+2);
    colNMatrices.resize(1+6+2);
    col_mat = 0;
    ierr = GetColIndices("DISPLACEMENT",ColGlob[col_mat]); CHKERRQ(ierr);
    ierr = GetGaussColNMatrix("DISPLACEMENT",colNMatrices[col_mat]); CHKERRQ(ierr);
    col_mat++;
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    for(int ee = 0;ee<3;ee++) { //edges matrices
	ierr = GetColIndices("DISPLACEMENT",MBEDGE,ColGlob[col_mat],ee+6); CHKERRQ(ierr);
	if(ColGlob[col_mat].size()!=0) {
	  ierr = GetGaussColNMatrix("DISPLACEMENT",MBEDGE,colNMatrices[col_mat],ee+6); CHKERRQ(ierr);
	  col_mat++;
	}
    }
    ierr = GetColIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],3); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],3); CHKERRQ(ierr);
	col_mat++;
    }
    ierr = GetColIndices("DISPLACEMENT",MBTRI,ColGlob[col_mat],4); CHKERRQ(ierr);
    if(ColGlob[col_mat].size()!=0) {
	ierr = GetGaussColNMatrix("DISPLACEMENT",MBTRI,colNMatrices[col_mat],4); CHKERRQ(ierr);
	col_mat++;
    }
   
    //Apply Dirihlet BC
    ApplyDirihletBC();

    //Assemble interface
    int g_dim = g_NTRI.size()/3;
    ublas::matrix<FieldData> K[row_mat][col_mat];
    for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	    ///K matrices
	    if(gg == 0) {
	      K[rr][cc] = ublas::zero_matrix<FieldData>(row_Mat.size2(),col_Mat.size2());
	    }
	    double w = area3*G_TRI_W13[gg];
	    ublas::matrix<FieldData> NTD = prod( trans(row_Mat), w*Dglob );
	    K[rr][cc] += prod(NTD , col_Mat ); 
	  }
	}
	if(RowGlob[rr].size()==0) continue;
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K[rr][cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K[rr][cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K[rr][cc].data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
    }

    ierr = OpStudentEnd(); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }


};

struct PostProcCohesiveForces: public InterfaceFEMethod,PostProcDisplacemenysAndStarinOnRefMesh_Base {
  
    PostProcCohesiveForces(Interface &_moab,double _YoungModulus): 
      InterfaceFEMethod(_moab,_YoungModulus), PostProcDisplacemenysAndStarinOnRefMesh_Base() {
      pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
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
      moabField_Core core_ref(moab_ref);
      moabField& mField_ref = core_ref;
      ierr = mField_ref.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);

      for(int ll = 0;ll<max_level;ll++) {
	PetscPrintf(PETSC_COMM_WORLD,"Refine Level %d\n",ll);
	rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[ll]); CHKERR_PETSC(rval);
	ierr = mField_ref.refine_get_ents(BitRefLevel().set(ll),meshset_level[ll]); CHKERRQ(ierr);
	ierr = mField_ref.add_verices_in_the_middel_of_edges(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
	ierr = mField_ref.refine_PRISM(meshset_level[ll],BitRefLevel().set(ll+1)); CHKERRQ(ierr);
      }
      rval = moab_ref.create_meshset(MESHSET_SET,meshset_level[max_level]); CHKERR_PETSC(rval);
      ierr = mField_ref.refine_get_ents(BitRefLevel().set(max_level),meshset_level[max_level]); CHKERRQ(ierr);

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
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	} 
	//faces
	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[gg*nb_dofs_H1face + approx_dof];
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
	  gap_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus*/
	}
	//faces
	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,4));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[nodes_on_face3.size()*nb_dofs_H1face + gg*nb_dofs_H1face + approx_dof];
	  disp_ptr[dit->get_dof_rank()] -= val*dit->get_FieldData(); //*minus/
	  gap_ptr[dit->get_dof_rank()] += val*dit->get_FieldData(); 
	}
 	dit = data_multiIndex->get<Composite_mi_tag>().find(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	hi_dit = data_multiIndex->get<Composite_mi_tag>().upper_bound(boost::make_tuple("DISPLACEMENT",MBTRI,3));
	for(;dit!=hi_dit;dit++) {
	  double *_H1faceN_ = &H1faceN[dit->side_number_ptr->side_number][0];
	  int nb_dofs_H1face = dit->get_order_nb_dofs(maxOrderEdgeH1[dit->side_number_ptr->side_number]);
	  int approx_dof = floor((double)dit->get_EntDofIdx()/(double)dit->get_max_rank());
	  double val = _H1faceN_[nodes_on_face3.size()*nb_dofs_H1face + gg*nb_dofs_H1face + approx_dof];
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

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //Reade parameters from line command
  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //We need that for code profiling
  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  moabField_Core core(moab);
  moabField& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  //Interface
  EntityHandle meshset_interface;
  ierr = mField.get_msId_meshset(4,SideSet,meshset_interface); CHKERRQ(ierr);
  ierr = mField.get_msId_3dENTS_sides(meshset_interface,true); CHKERRQ(ierr);
  // stl::bitset see for more details
  BitRefLevel bit_level_interface;
  bit_level_interface.set(0);
  ierr = mField.get_msId_3dENTS_split_sides(0,bit_level_interface,meshset_interface,true,true); CHKERRQ(ierr);
  EntityHandle meshset_level_interface;
  rval = moab.create_meshset(MESHSET_SET,meshset_level_interface); CHKERR_PETSC(rval);
  ierr = mField.refine_get_ents(bit_level_interface,meshset_level_interface); CHKERRQ(ierr);

  //update BC for refined (with interface) mesh
  EntityHandle meshset_SideSet1; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(1,SideSet,meshset_SideSet1); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet1,bit_level_interface,meshset_SideSet1,MBTRI,true,3); CHKERRQ(ierr);
  EntityHandle meshset_SideSet2; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(2,SideSet,meshset_SideSet2); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet2,bit_level_interface,meshset_SideSet2,MBTRI,true,3); CHKERRQ(ierr);
  EntityHandle meshset_SideSet3; //Dirihlet BC is there
  ierr = mField.get_msId_meshset(3,SideSet,meshset_SideSet3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet3,bit_level_interface,meshset_SideSet3,MBTRI,true,3); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(1);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(meshset_level_interface,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  BitRefLevel bit_level1;
  bit_level1.set(2);
  ierr = mField.add_verices_in_the_middel_of_edges(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_TET(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_PRISM(meshset_level0,bit_level1); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet1,bit_level1,meshset_SideSet1,MBTRI,true,3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet2,bit_level1,meshset_SideSet2,MBTRI,true,3); CHKERRQ(ierr);
  ierr = mField.refine_get_childern(meshset_SideSet3,bit_level1,meshset_SideSet3,MBTRI,true,3); CHKERRQ(ierr);

  BitRefLevel problem_bit_level = bit_level1;

  /***/
  //Define problem

  //Fields
  ierr = mField.add_BitFieldId("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_MoFEMFE("ELASTIC"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = mField.modify_MoFEMFE_row_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_col_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_data_add_bit("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  //FE Interface
  ierr = mField.add_MoFEMFE("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_row_add_bit("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_col_add_bit("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_MoFEMFE_data_add_bit("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_BitProblemId("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_MoFEMFE_add_bit("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_MoFEMFE_add_bit("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_MoFEMFE_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_MoFEMFE_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",5); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  /****/
  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning 

  //partition
  ierr = mField.partition_problems("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2,SideSet3;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(3,SideSet,2,SideSet3,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 3 : %u\n",SideSet3.size());

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  const double alpha = 0.05;
  MyElasticFEMethod MyFE(moab,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),SideSet1,SideSet2,SideSet3);
  InterfaceFEMethod IntMyFE(moab,Aij,F,YoungModulus*alpha,SideSet1,SideSet2,SideSet3);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);


  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcDisplacementsEntMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PostProcCohesiveForces fe_post_proc_prisms(moab,YoungModulus*alpha);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
  }


  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);


  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

