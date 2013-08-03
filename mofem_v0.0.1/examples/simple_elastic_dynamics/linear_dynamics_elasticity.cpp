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
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
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

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("VELOCITIES",L2,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("STIFFNESS"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("MASS"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COPUPLING_VV"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("COPUPLING_VU"); CHKERRQ(ierr);


  //Define rows/cols and element data
  //STIFFNESS
  ierr = mField.modify_finite_element_add_field_row("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("STIFFNESS","DISPLACEMENT"); CHKERRQ(ierr);

  //MASS
  ierr = mField.modify_finite_element_add_field_row("MASS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MASS","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MASS","DISPLACEMENT"); CHKERRQ(ierr);

  //COPUPLING
  ierr = mField.modify_finite_element_add_field_row("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COPUPLING_VV","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("COPUPLING_VU","VELOCITIES"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("COPUPLING_VU","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","STIFFNESS"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","MASS"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VV"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","COPUPLING_VU"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"VELOCITIES"); CHKERRQ(ierr);


  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"STIFFNESS",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MASS",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VV",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"COPUPLING_VU",MBTET); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  const int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"VELOCITIES",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"VELOCITIES",1); CHKERRQ(ierr);

  /****/
  //build database

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

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
  Vec D,F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&D); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());

  /*
   * M*u'' + K*u' - F = 0
   *
   * F( t, [ dot_u, u], [ dot_u', u'] ) = [ 0 1 ][ dot_u' ] + [ -1 0 ][ dot_u ] + [ 0    ] = [ 0 ]
   *                                      [ M 0 ][ u'     ]   [  0 K ][ u     ]   [ F(t) ]   [ 0 ]
   * 
   * Fu  = [ -1 0 ][ dot_u ]
   *       [  0 K ][ u     ] 
   * 
   * Fu' = [ 0 1 ][ dot_u' ]
   *       [ M 0 ][ u'     ]
   *
   * G = Fu + a*Fu'
   *
   * dG/(du,ddot_u) = [ -1 0 ] + a*[ 0 1 ]
   *                  [  0 K ] +   [ M 0 ]
   *
   */

  struct DynamicElasticFEMethod: public ElasticFEMethod {


    double rho;
    DynamicElasticFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,
      double _lambda,double _mu,double _rho,Range &_SideSet1,Range &_SideSet2): 
      ElasticFEMethod(_moab,_Aij,_F,_lambda,_mu,
      _SideSet1,_SideSet2),rho(_rho) {};

    /// Set Neumann Boundary Conditions on SideSet2
    PetscErrorCode NeumannBC() {
      PetscFunctionBegin;
      ublas::vector<FieldData,ublas::bounded_array<double,3> > traction(3);
      //Set Direction of Traction On SideSet2
      traction[0] = 1; //X
      traction[1] = 0; //Y 
      traction[2] = 0; //Z
      //ElasticFEMethod::NeumannBC(...) function calulating external forces (see file ElasticFEMethod.hpp)
      ierr = ElasticFEMethod::NeumannBC(ts_F,traction,SideSet2); CHKERRQ(ierr);
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

    ublas::vector<FieldData> dot_velocities;
    vector<ublas::vector<FieldData> > velocities;
    PetscErrorCode GetVelocities_Form_TS_u_t() {
      PetscFunctionBegin;
      Vec u_t_local;
      ierr = VecGhostGetLocalForm(ts_u_t,&u_t_local); CHKERRQ(ierr);
      int local_size;
      ierr = VecGetLocalSize(u_t_local,&local_size); CHKERRQ(ierr);
      double *array;
      ierr = VecGetArray(u_t_local,&array); CHKERRQ(ierr);
      dot_velocities.resize(VelColLocal.size());
      int ii = 0;
      vector<DofIdx>::iterator it = VelColLocal.begin();
      for(;it!=VelColLocal.end();it++,ii++) {
	  if(*it < 0) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(*it >= local_size) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  dot_velocities[ii] = array[*it];
      }
      int cc = 0;
      velocities.resize(col_mat);
      for(;cc<col_mat;cc++) {
	velocities[cc].resize(ColLocal.size());
	int iii = 0;
	vector<DofIdx>::iterator iit = ColLocal[cc].begin();
	for(;iit!=ColLocal[cc].end();iit++) {
	  (velocities[cc])[iii] = array[*iit];
	}
      }
      ierr = VecRestoreArray(u_t_local,&array); CHKERRQ(ierr);
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
      Mass.resize(row_mat);
      int g_dim = g_NTET.size()/4;
      for(int cc = 0;cc<col_mat;cc++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  ublas::matrix<FieldData> &row_Mat = VelRowNMatrix[gg];
	  ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	  double w = rho*V*G_TET_W45[gg];
	  if(gg == 0) {
	    VU[cc] = w*prod( trans(row_Mat), col_Mat );
	  } else {
	    VU[cc] += w*prod( trans(row_Mat), col_Mat );
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

      if(fe_ent_ptr->get_name()=="STIFFNESS") {
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
	  case ctx_TSSetIFunction:
	    ierr = VecZeroEntries(ts_F); CHKERRQ(ierr);
	    ierr = VecGhostUpdateBegin(ts_F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(ts_F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
	    break;
	  case ctx_TSSetIJacobian:
	    ierr = MatZeroEntries(*ts_B); CHKERRQ(ierr);
	    ierr = VecDuplicate(ts_F,&Diagonal); CHKERRQ(ierr);
	    break;
	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="MASS") {
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="COPUPLING_VV") {
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="COPUPLING_VU") {
	PetscFunctionReturn(0);
      }

      SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      if(fe_ent_ptr->get_name()=="STIFFNESS") {
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="MASS") {
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="COPUPLING_VV") {
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="COPUPLING_VU") {

	switch (ts_ctx) {
	  case ctx_TSSetIFunction:
	    ierr = VecGhostUpdateBegin(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecGhostUpdateEnd(ts_F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
	    ierr = VecAssemblyBegin(ts_F); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(ts_F); CHKERRQ(ierr);
	    break;
	  case ctx_TSSetIJacobian:
	    ierr = VecAssemblyBegin(Diagonal); CHKERRQ(ierr);
	    ierr = VecAssemblyEnd(Diagonal); CHKERRQ(ierr);
	    ierr = MatDiagonalSet(*ts_B,Diagonal,ADD_VALUES); CHKERRQ(ierr);
	    ierr = VecDestroy(&Diagonal); CHKERRQ(ierr);
	    Diagonal = PETSC_NULL;
	    ierr = MatAssemblyBegin(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    ierr = MatAssemblyEnd(*ts_B,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	    break;
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

      if(fe_ent_ptr->get_name()=="STIFFNESS") {
	ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
	ierr = GetMatrices(); CHKERRQ(ierr);
	//ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	//Dirihlet Boundary Condition
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetIFunction: {
	    ierr = NeumannBC(); CHKERRQ(ierr);
	    ierr = Fint(ts_F); CHKERRQ(ierr);
	  } break;
	  case ctx_TSSetIJacobian: {
	    if(Diagonal==PETSC_NULL) SETERRQ(PETSC_COMM_SELF,1,"DrihletBC on diagonal imposible to set"); 
	    if(DirihletBC.size()>0) {
	      DirihletBCDiagVal.resize(DirihletBC.size());
	      fill(DirihletBCDiagVal.begin(),DirihletBCDiagVal.end(),1);
	      ierr = VecSetValues(Diagonal,DirihletBC.size(),&(DirihletBC[0]),&DirihletBCDiagVal[0],INSERT_VALUES); CHKERRQ(ierr);
	    }
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
      if(fe_ent_ptr->get_name()=="MASS") {
	ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
	ierr = GetMatrices(); CHKERRQ(ierr);
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	ierr = MassLhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetIFunction: {
	    //Inercia
	    ierr = MassLhs(); CHKERRQ(ierr);
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    for(int rr = 0;rr<row_mat;rr++) {
	      ublas::vector<FieldData> Mu_t = prod(Mass[rr],dot_velocities);
	      if(RowGlob.size()!=Mu_t.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(RowGlob[rr].size()==0) continue;
	      ierr = VecSetValues(ts_F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(Mu_t.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
	  case ctx_TSSetIJacobian: {
	    Mass *= ts_a;
	    for(int rr = 0;rr<row_mat;rr++) {
	      if(RowGlob[rr].size()==0) continue;
	      if(RowGlob[rr].size()!=Mass[rr].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(VelColGlob.size()!=Mass[rr].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = MatSetValues(*ts_B,RowGlob[rr].size(),&(RowGlob[rr])[0],VelColGlob.size(),&VelColGlob[0],&(Mass[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
  	  default:
	    SETERRQ(PETSC_COMM_SELF,1,"sorry... I don't know what to do");
	}
	ierr = OpStudentEnd(); CHKERRQ(ierr);
	PetscFunctionReturn(0);
      }
      if(fe_ent_ptr->get_name()=="COPUPLING_VV") {
	ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = VVLhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetIFunction: {
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    ublas::vector<FieldData> VVu_t = prod(VV,dot_velocities);
	    if(VelRowGlob.size()!=VVu_t.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = VecSetValues(ts_F,VelRowGlob.size(),&VelRowGlob[0],&(VVu_t.data()[0]),ADD_VALUES); CHKERRQ(ierr);
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
      if(fe_ent_ptr->get_name()=="COPUPLING_VU") {
	ierr = OpStudentStart_TET(g_NTET); CHKERRQ(ierr);
	ierr = GetMatrices(); CHKERRQ(ierr);
	ierr = GetMatricesVelocities(); CHKERRQ(ierr);
	ierr = ApplyDirihletBC(); CHKERRQ(ierr);
	ierr = VULhs(); CHKERRQ(ierr);
	switch (ts_ctx) {
	  case ctx_TSSetIFunction: {
	    ierr = GetVelocities_Form_TS_u_t(); CHKERRQ(ierr);
	    int cc = 0;
	    for(;cc<col_mat;cc++) {
	      ublas::vector<FieldData> VUu_t = prod(VU[cc],velocities[cc]);
	      if(VelRowGlob.size()!=VUu_t.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = VecSetValues(ts_F,VelRowGlob.size(),&VelRowGlob[0],&(VUu_t.data()[0]),ADD_VALUES); CHKERRQ(ierr);
	    }
	  } break;
	  case ctx_TSSetIJacobian: {
	    VU *= ts_a;
	    int cc = 0;	
	    for(;cc<<col_mat;cc++) {
	      if(VelRowGlob.size()!=VU[cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(ColGlob[cc].size()!=VU[cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
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

  //TS
  moabTsCtx TsCtx(mField,"ELASTIC_MECHANICS");

  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  ierr = TSSetIFunction(ts,PETSC_NULL,f_TSSetIFunction,&TsCtx); CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,Aij,Aij,f_TSSetIJacobian,&TsCtx); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,0.0,.001); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,D); CHKERRQ(ierr);

  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  //detroy matrices
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

