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

#include "FieldInterface.hpp"
#include "FieldCore.hpp"
#include "FEMethod_UpLevelStudent.hpp"
#include "cholesky.hpp"
#include <petscksp.h>

#include "ElasticFEMethod.hpp"
#include "ElasticFEMethodInterface.hpp"

#include "SurfacePressure.hpp"
#include "BodyForce.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct ExampleDiriheltBC: public BaseDirihletBC {
  Range SideSet1_;

  string field_name;
  ExampleDiriheltBC(Interface &moab,Range& SideSet1): field_name("DISPLACEMENT") {
	ErrorCode rval;
	Range SideSet1Edges,SideSet1Nodes;
	rval = moab.get_adjacencies(SideSet1,1,false,SideSet1Edges,Interface::UNION); CHKERR_THROW(rval);
	rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_THROW(rval);
	SideSet1_.insert(SideSet1.begin(),SideSet1.end());
	SideSet1_.insert(SideSet1Edges.begin(),SideSet1Edges.end());
	SideSet1_.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
  }

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
      FieldInterface::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = SideSet1_.begin();
      for(;siit1!=SideSet1_.end();siit1++) {
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_riit = fe_method_ptr->row_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;riit!=hi_riit;riit++) {
	    if(riit->get_name()!=field_name) continue;
	    // all fixed
	    // if some ranks are selected then we could apply BC in particular direction
	    DirihletBC.push_back(riit->get_petsc_gloabl_dof_idx());
	    for(unsigned int rr = 0;rr<RowGlob.size();rr++) {
	      vector<DofIdx>::iterator it = find(RowGlob[rr].begin(),RowGlob[rr].end(),riit->get_petsc_gloabl_dof_idx());
	      if( it!=RowGlob[rr].end() ) *it = -1; // of idx is set -1 row is not assembled
	    }
	  }
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().lower_bound(*siit1);
	  FENumeredDofMoFEMEntity_multiIndex::index<MoABEnt_mi_tag>::type::iterator hi_ciit = fe_method_ptr->col_multiIndex->get<MoABEnt_mi_tag>().upper_bound(*siit1);
	  for(;ciit!=hi_ciit;ciit++) {
	    if(ciit->get_name()!=field_name) continue;
	    for(unsigned int cc = 0;cc<ColGlob.size();cc++) {
	      vector<DofIdx>::iterator it = find(ColGlob[cc].begin(),ColGlob[cc].end(),ciit->get_petsc_gloabl_dof_idx());
	      if( it!=ColGlob[cc].end() ) *it = -1; // of idx is set -1 column is not assembled
	    }
	  }
      }
      PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_ElementIndiciesFace(vector<DofIdx>& DirihletBC,
      vector<DofIdx> &FaceNodeIndices,
      vector<vector<DofIdx> > &FaceEdgeIndices,
      vector<DofIdx> &FaceIndices) {
      PetscFunctionBegin;
      vector<DofIdx>::iterator dit = DirihletBC.begin();
      for(;dit!=DirihletBC.end();dit++) {
	vector<DofIdx>::iterator it = find(FaceNodeIndices.begin(),FaceNodeIndices.end(),*dit);
	if(it!=FaceNodeIndices.end()) *it = -1; // of idx is set -1 row is not assembled
	for(int ee = 0;ee<3;ee++) {
	  it = find(FaceEdgeIndices[ee].begin(),FaceEdgeIndices[ee].end(),*dit);
	  if(it!=FaceEdgeIndices[ee].end()) *it = -1; // of idx is set -1 row is not assembled
	}
	it = find(FaceIndices.begin(),FaceIndices.end(),*dit);
	if(it!=FaceIndices.end()) *it = -1; // of idx is set -1 row is not assembled
      }
      PetscFunctionReturn(0);
  }

};

int main(int argc, char *argv[]) {

  try {

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
  ierr = PetscTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  //Create MoFEM (Joseph) database
  FieldCore core(moab);
  FieldInterface& mField = core;

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,BitRefLevel().set(0)); CHKERRQ(ierr);
  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(0));

  int ll = 1;
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SideSet|InterfaceSet,cit)) {
    //for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,SideSet,cit)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Insert Interface %d\n",cit->get_msId()); CHKERRQ(ierr);
    EntityHandle cubit_meshset = cit->get_meshset();
    {
      //get tet enties form back bit_level
      EntityHandle ref_level_meshset = 0;
      rval = moab.create_meshset(MESHSET_SET,ref_level_meshset); CHKERR_PETSC(rval);
      ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBTET,ref_level_meshset); CHKERRQ(ierr);
      ierr = mField.get_entities_by_type_and_ref_level(bit_levels.back(),BitRefLevel().set(),MBPRISM,ref_level_meshset); CHKERRQ(ierr);
      Range ref_level_tets;
      rval = moab.get_entities_by_handle(ref_level_meshset,ref_level_tets,true); CHKERR_PETSC(rval);
      //get faces and test to split
      ierr = mField.get_msId_3dENTS_sides(cubit_meshset,bit_levels.back(),true,0); CHKERRQ(ierr);
      //set new bit level
      bit_levels.push_back(BitRefLevel().set(ll++));
      //split faces and 
      ierr = mField.get_msId_3dENTS_split_sides(ref_level_meshset,bit_levels.back(),cubit_meshset,true,true,0); CHKERRQ(ierr);
      //clean meshsets
      rval = moab.delete_entities(&ref_level_meshset,1); CHKERR_PETSC(rval);
    }
    //update cubit meshsets
    for(_IT_CUBITMESHSETS_FOR_LOOP_(mField,ciit)) {
      EntityHandle cubit_meshset = ciit->meshset; 
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBVERTEX,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBEDGE,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTRI,true); CHKERRQ(ierr);
      ierr = mField.update_meshset_by_entities_children(cubit_meshset,bit_levels.back(),cubit_meshset,MBTET,true); CHKERRQ(ierr);
    }
  }
    
  BitRefLevel problem_bit_level = bit_levels.back();

  /***/
  //Define problem

  //Fields
  ierr = mField.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  //FE Interface
  ierr = mField.add_finite_element("INTERFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("INTERFACE","DISPLACEMENT"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","INTERFACE"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);

  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"ELASTIC",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(problem_bit_level,"INTERFACE",MBPRISM); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 4;
  }
  ierr = mField.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  PetscBool oscillationBool;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_interface_order_reduction",&oscillationBool,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    oscillationBool=PETSC_FALSE;}
  
  if(oscillationBool==PETSC_TRUE) {
    //reduce level of approximation for entities on inetrface
    Range prims;
    ierr = mField.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims); CHKERRQ(ierr);
    Range prims_faces;
    rval = mField.get_moab().get_adjacencies(prims,2,false,prims_faces,Interface::UNION); CHKERR_PETSC(rval);
    Range prims_faces_edges;
    rval = mField.get_moab().get_adjacencies(prims_faces,1,false,prims_faces_edges,Interface::UNION); CHKERR_PETSC(rval);
    ierr = mField.set_field_order(prims_faces,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);
    ierr = mField.set_field_order(prims_faces_edges,"DISPLACEMENT",order>1 ? order-1 : 0); CHKERRQ(ierr);
  }
  /****/

  ierr = MetaNeummanForces::addNeumannBCElements(mField,"ELASTIC_MECHANICS","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("BODY_FORCE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("BODY_FORCE","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","BODY_FORCE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Block_BodyForcesSet,it)) {
    Range tets;
    rval = mField.get_moab().get_entities_by_type(it->meshset,MBTET,tets,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TETs(tets,"BODY_FORCE"); CHKERRQ(ierr);
  }

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
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //print bcs
  ierr = mField.printCubitDisplacementSet(); CHKERRQ(ierr);
  ierr = mField.printCubitForceSet(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.printCubitMaterials(); CHKERRQ(ierr);

  //create matrices
  Vec D,F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Col,&D); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",Row,&F); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _mField,BaseDirihletBC *_dirihlet_ptr,
      Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu): 
      ElasticFEMethod(_mField,_dirihlet_ptr,_Aij,_D,_F,_lambda,_mu) {};

    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
	if(RowGlob[rr].size()==0) continue;
	if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
	ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  //Assemble F and Aij
  const double YoungModulus = 1000;
  const double PoissonRatio = 0.2;
  const double alpha = 0.1;
  CubitDisplacementDirihletBC myDirihletBC(mField,"ELASTIC_MECHANICS","DISPLACEMENT");
  ierr = myDirihletBC.Init(); CHKERRQ(ierr);
  MyElasticFEMethod MyFE(mField,&myDirihletBC,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  InterfaceFEMethod IntMyFE(mField,&myDirihletBC,Aij,D,F,YoungModulus*alpha);

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",IntMyFE);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
  //Neuman BCs and BodyForces
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setFiniteElementOperators(mField,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("ELASTIC_MECHANICS",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  BodyFroceConstantField body_forces_methods(mField);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,BlockSet|Block_BodyForcesSet,it)) {
    ierr = body_forces_methods.addBlock("DISPLACEMENT",F,it->get_msId()); CHKERRQ(ierr);
  }
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","BODY_FORCE",body_forces_methods.getLoopFe()); CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
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

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PostProcVertexMethod ent_method(moab);
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",Row,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }
  
  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(
    mField,"DISPLACEMENT",LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }
  

  PostProcCohesiveForces fe_post_proc_prisms(mField,YoungModulus*alpha);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  }

}

