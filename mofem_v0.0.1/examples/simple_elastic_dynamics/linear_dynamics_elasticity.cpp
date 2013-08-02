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
  ierr = mField.add_field("VELOCITIES",H1,3); CHKERRQ(ierr);


  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","VELOCITIES"); CHKERRQ(ierr);


  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

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
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  //Get SideSet 1 and SideSet 2 defined in CUBIT
  Range SideSet1,SideSet2;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  ierr = mField.get_Cubit_msId_entities_by_dimension(2,SideSet,2,SideSet2,true); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 1 : %u\n",SideSet1.size());
  PetscPrintf(PETSC_COMM_WORLD,"Nb. faces in SideSet 2 : %u\n",SideSet2.size());

  struct MyElasticFEMethod: public ElasticFEMethod {
    double rho;
    MyElasticFEMethod(Interface& _moab,Mat &_Aij,Vec& _F,
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
      ierr = ElasticFEMethod::NeumannBC(F,traction,SideSet2); CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    PetscErrorCode MassLhs() {
      PetscFunctionBegin;
      ublas::matrix<FieldData> K[row_mat][col_mat];
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int cc = 0;cc<col_mat;cc++) {
	  for(int gg = 0;gg<g_dim;gg++) {
	    ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	    ublas::matrix<FieldData> &col_Mat = (colNMatrices[cc])[gg];
	    ///K matrices
	    double w = rho*V*G_TET_W45[gg];
	    if(gg == 0) {
	      K[rr][cc] = w*prod( trans(row_Mat), col_Mat );
	    } else {
	      K[rr][cc] += w*prod( trans(row_Mat), col_Mat );
	    }
	  }
	}
	for(int cc = 0;cc<col_mat;cc++) {
	  if(ColGlob[cc].size()==0) continue;
	  if(RowGlob[rr].size()!=K[rr][cc].size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  if(ColGlob[cc].size()!=K[rr][cc].size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	  ierr = MatSetValues(Aij,RowGlob[rr].size(),&(RowGlob[rr])[0],ColGlob[cc].size(),&(ColGlob[cc])[0],&(K[rr][cc].data())[0],ADD_VALUES); CHKERRQ(ierr);
	}
      }
      PetscFunctionReturn(0);
    }

    PetscErrorCode MassRhs() {
      PetscFunctionBegin;
      //ublas::vector<FieldData> f[row_mat];
      int g_dim = g_NTET.size()/4;
      for(int rr = 0;rr<row_mat;rr++) {
	for(int gg = 0;gg<g_dim;gg++) {
	  //ublas::matrix<FieldData> &row_Mat = (rowNMatrices[rr])[gg];
	  //ublas::matrix<FieldData> &col_Mat = (colNMatrices[rr])[gg];
	  //if(gg == 0) f[rr] = ublas::zero_vector<FieldData>(row_Mat.size2());
	  //double w = rho*V*G_TET_W45[gg];
	  //ublas::vector<FieldData> velocities(col_Mat.size2());


	}
	//if(RowGlob[rr].size()!=f[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	//if(RowGlob[rr].size()==0) continue;
	//ierr = VecSetValues(F,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f[rr].data())[0],ADD_VALUES); CHKERRQ(ierr);
      }
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

  };

  //TE
  moabTsCtx TsCtx(mField,"ELASTIC_MECHANICS");

  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

  //Assemble F and Aij
  const double YoungModulus = 1;
  const double PoissonRatio = 0.0;
  const double rho = 1;
  MyElasticFEMethod MyFE(moab,Aij,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio),rho,SideSet1,SideSet2);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::string wait;
  //std::cin >> wait;

  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

}

