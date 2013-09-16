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
#include <petscksp.h>

#include "moabSnes.hpp"
#include "PostProcVertexMethod.hpp"
#include "PostProcDisplacementAndStrainOnRefindedMesh.hpp"

#include "moabFEMethod_DriverComplexForLazy.hpp"

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct MyMeshSmoothing_ElasticFEMethod: public FEMethod_DriverComplexForLazy_MeshSmoothingProjected {

  MyMeshSmoothing_ElasticFEMethod(moabField& _mField,matPROJ_ctx &_proj_all_ctx,BaseDirihletBC *_dirihlet_bc_method_ptr,int _verbose = 0):
    FEMethod_DriverComplexForLazy_MeshSmoothingProjected(_mField,_proj_all_ctx,_dirihlet_bc_method_ptr,_verbose) {
    set_qual_ver(2);
  }

};

struct materialDirihletBC: public BaseDirihletBC {

  Interface& moab;
  Range &CornersNodes;
  string field_name;
  materialDirihletBC(Interface &_moab,Range& _CornerNodes): moab(_moab),CornersNodes(_CornerNodes),field_name("MESH_NODE_POSITIONS") {}

  PetscErrorCode SetDirihletBC_to_ElementIndicies(
      moabField::FEMethod *fe_method_ptr,vector<vector<DofIdx> > &RowGlob,vector<vector<DofIdx> > &ColGlob,vector<DofIdx>& DirihletBC) {
      PetscFunctionBegin;
      //Dirihlet form SideSet1
      DirihletBC.resize(0);
      Range::iterator siit1 = CornersNodes.begin();
      for(;siit1!=CornersNodes.end();siit1++) {
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

  PetscErrorCode SetDirihletBC_to_MatrixDiagonal(
    moabField::FEMethod *fe_method_ptr,Mat Aij) {
    PetscFunctionBegin;

    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);

    for(_IT_NUMEREDDOFMOFEMENTITY_ROW_BY_LOCIDX_FOR_LOOP_(fe_method_ptr->problem_ptr,dit)) {
      if(dit->get_part()!=pcomm->rank()) continue;
      if(find(CornersNodes.begin(),CornersNodes.end(),dit->get_ent()) == CornersNodes.end()) continue;
      ierr = MatSetValue(Aij,dit->get_petsc_gloabl_dof_idx(),dit->get_petsc_gloabl_dof_idx(),1.,INSERT_VALUES); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
  }

  PetscErrorCode SetDirihletBC_to_RHS(moabField::FEMethod *fe_method_ptr,Vec F) {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};



int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }
 
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  PetscLogDouble t1,t2;
  PetscLogDouble v1,v2;
  ierr = PetscGetTime(&v1); CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t1); CHKERRQ(ierr);

  moabField_Core core(moab);
  moabField& mField = core;

  Range CubitSideSets_meshsets;
  ierr = mField.get_CubitBCType_meshsets(SideSet,CubitSideSets_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.refine_get_ents(bit_level0,meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_SURFACE",H1,1); CHKERRQ(ierr);
  ierr = mField.add_field("LAMBDA_CORNER",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("MESH_SMOOTHER"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.add_finite_element("CTC_CORNER_ELEM"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("MESH_SMOOTHER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_SURFACE_ELEM","LAMBDA_SURFACE"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("C_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  ierr = mField.modify_finite_element_add_field_row("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("CTC_CORNER_ELEM","LAMBDA_CORNER"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.add_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.add_problem("C_ALL_MATRIX"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("MESH_SMOOTHING","MESH_SMOOTHER"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("CCT_ALL_MATRIX","CTC_SURFACE_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_CORNER_ELEM"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("C_ALL_MATRIX","C_SURFACE_ELEM"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("MESH_SMOOTHING",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("CCT_ALL_MATRIX",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("C_ALL_MATRIX",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"MESH_SMOOTHER",MBTET); CHKERRQ(ierr);

  //add tets on corners
  EntityHandle CornersNodesMeshset,SurfacesFacesMeshset;
  {
    Range CornersEdges,CornersNodes,SurfacesFaces;
    ierr = mField.get_Cubit_msId_entities_by_dimension(100,SideSet,1,CornersEdges,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
    ierr = mField.get_Cubit_msId_entities_by_dimension(102,SideSet,2,SurfacesFaces,true); CHKERRQ(ierr);
    Range SurfacesTets;
    rval = moab.get_adjacencies(SurfacesFaces,3,false,SurfacesTets,Interface::UNION); CHKERR_PETSC(rval);
    {
      Range CornersEdgesNodes;
      rval = moab.get_adjacencies(CornersEdges,0,false,CornersEdgesNodes,Interface::UNION); CHKERR_PETSC(rval);
      rval = moab.create_meshset(MESHSET_SET,CornersNodesMeshset); CHKERR_PETSC(rval);	
      CornersNodes.insert(CornersEdgesNodes.begin(),CornersEdgesNodes.end());
      Range SideSet1;
      ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
      Range SideSet1Nodes;
      rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_PETSC(rval);
      CornersNodes.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
      rval = moab.add_entities(CornersNodesMeshset,CornersNodes); CHKERR_PETSC(rval);
      //add surface elements
      Range CornersTets;
      rval = moab.get_adjacencies(CornersNodes,3,false,CornersTets,Interface::UNION); CHKERR_PETSC(rval);
      EntityHandle CornersTetsMeshset;
      rval = moab.create_meshset(MESHSET_SET,CornersTetsMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(CornersTetsMeshset,CornersTets); CHKERR_PETSC(rval);
      CornersTets = intersect(CornersTets,SurfacesTets);
      ierr = mField.add_ents_to_finite_element_by_TETs(CornersTetsMeshset,"C_CORNER_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(CornersTetsMeshset,"CTC_CORNER_ELEM"); CHKERRQ(ierr);
      rval = moab.delete_entities(&CornersTetsMeshset,1); CHKERR_PETSC(rval);
    }
    {
      Range SurfacesNodes;
      rval = moab.get_adjacencies(SurfacesFaces,0,false,SurfacesNodes,Interface::UNION); CHKERR_PETSC(rval);
      Range CornersEdgesNodes;
      rval = moab.get_connectivity(CornersEdges,CornersEdgesNodes,true); CHKERR_PETSC(rval);
      Range SideSet1;
      ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
      Range SideSet1Nodes;
      rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_PETSC(rval);
      SurfacesNodes = subtract(SurfacesNodes,CornersNodes);
      SurfacesNodes = subtract(SurfacesNodes,CornersEdgesNodes);
      SurfacesNodes = subtract(SurfacesNodes,SideSet1Nodes);
      rval = moab.create_meshset(MESHSET_SET,SurfacesFacesMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(SurfacesFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
      rval = moab.add_entities(SurfacesFacesMeshset,SurfacesNodes); CHKERR_PETSC(rval);
      //add surface elements
      EntityHandle SurfacesTetsMeshset;
      rval = moab.create_meshset(MESHSET_SET,SurfacesTetsMeshset); CHKERR_PETSC(rval);	
      rval = moab.add_entities(SurfacesTetsMeshset,SurfacesTets); CHKERR_PETSC(rval);
      ierr = mField.add_ents_to_finite_element_by_TETs(SurfacesTetsMeshset,"C_SURFACE_ELEM"); CHKERRQ(ierr);
      ierr = mField.add_ents_to_finite_element_by_TETs(SurfacesTetsMeshset,"CTC_SURFACE_ELEM"); CHKERRQ(ierr);
      rval = moab.delete_entities(&SurfacesTetsMeshset,1); CHKERR_PETSC(rval);
    }
  }

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(SurfacesFacesMeshset,"LAMBDA_SURFACE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_VERTICEs(CornersNodesMeshset,"LAMBDA_CORNER"); CHKERRQ(ierr);

  //NOTE: always order should be 1
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_SURFACE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"LAMBDA_CORNER",1); CHKERRQ(ierr);

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("MESH_SMOOTHING"); CHKERRQ(ierr);
  //partition
  ierr = mField.partition_problem("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("CCT_ALL_MATRIX"); CHKERRQ(ierr);
  //partition
  ierr = mField.compose_problem("C_ALL_MATRIX","CCT_ALL_MATRIX","MESH_SMOOTHING"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("C_ALL_MATRIX"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("C_ALL_MATRIX"); CHKERRQ(ierr);

  {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_MOABFIELD_BY_NAME_FOR_LOOP_(mField,"MESH_NODE_POSITIONS",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      int dof_rank = dof_ptr->get_dof_rank();
      double &fval = dof_ptr->get_FieldData();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      fval = coords[dof_rank];
    }
  }

  matPROJ_ctx proj_all_ctx(mField,"MESH_SMOOTHING","C_ALL_MATRIX");
  Mat precK;
  ierr = mField.MatCreateMPIAIJWithArrays("MESH_SMOOTHING",&precK); CHKERRQ(ierr);
  ierr = MatDuplicate(precK,MAT_DO_NOT_COPY_VALUES,&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = mField.MatCreateMPIAIJWithArrays("C_ALL_MATRIX",&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("C_ALL_MATRIX",Row,&proj_all_ctx.g); CHKERRQ(ierr);

  Mat CTC_QTKQ;
  int M,N,m,n;
  ierr = MatGetSize(proj_all_ctx.K,&M,&N); CHKERRQ(ierr);
  ierr = MatGetLocalSize(proj_all_ctx.K,&m,&n); CHKERRQ(ierr);
  ierr = MatCreateShell(PETSC_COMM_WORLD,m,n,M,N,&proj_all_ctx,&CTC_QTKQ); CHKERRQ(ierr);
  ierr = MatShellSetOperation(CTC_QTKQ,MATOP_MULT,(void(*)(void))matCTC_QTKQ_mult_shell); CHKERRQ(ierr);
  Vec F;
  ierr = mField.VecCreateGhost("MESH_SMOOTHING",Row,&F); CHKERRQ(ierr);

  Range CornersNodes;
  ierr = mField.get_Cubit_msId_entities_by_dimension(101,NodeSet,0,CornersNodes,true); CHKERRQ(ierr);
  Range SideSet1;
  ierr = mField.get_Cubit_msId_entities_by_dimension(1,SideSet,2,SideSet1,true); CHKERRQ(ierr);
  Range SideSet1Nodes;
  rval = moab.get_connectivity(SideSet1,SideSet1Nodes,true); CHKERR_PETSC(rval);
  CornersNodes.insert(SideSet1Nodes.begin(),SideSet1Nodes.end());
  materialDirihletBC myDirihletBC(moab,CornersNodes);

  MyMeshSmoothing_ElasticFEMethod MyFE(mField,proj_all_ctx,&myDirihletBC);

  moabSnesCtx SnesCtx(mField,"MESH_SMOOTHING");

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,CTC_QTKQ,precK,SnesMat,&SnesCtx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  moabSnesCtx::loops_to_do_type& loops_to_do_Rhs = SnesCtx.get_loops_to_do_Rhs();
  loops_to_do_Rhs.push_back(moabSnesCtx::loop_pair_type("MESH_SMOOTHER",&MyFE));
  moabSnesCtx::loops_to_do_type& loops_to_do_Mat = SnesCtx.get_loops_to_do_Mat();
  loops_to_do_Mat.push_back(moabSnesCtx::loop_pair_type("MESH_SMOOTHER",&MyFE));

  Vec D;
  ierr = mField.VecCreateGhost("MESH_SMOOTHING",Col,&D); CHKERRQ(ierr);
  ierr = mField.set_local_VecCreateGhost("MESH_SMOOTHING",Col,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("MESH_SMOOTHING",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  PostProcVertexMethod ent_method(moab,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_SMOOTHING","MESH_NODE_POSITIONS",Col,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("MESH_SMOOTHING","MESH_SMOOTHER",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  ierr = proj_all_ctx.DestroyQorP(); CHKERRQ(ierr);
  ierr = proj_all_ctx.DestroyQTKQ(); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.K); CHKERRQ(ierr);
  ierr = MatDestroy(&precK); CHKERRQ(ierr);
  ierr = MatDestroy(&proj_all_ctx.C); CHKERRQ(ierr);
  ierr = VecDestroy(&proj_all_ctx.g); CHKERRQ(ierr);
  ierr = MatDestroy(&CTC_QTKQ); CHKERRQ(ierr);

  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);

  PetscFinalize();

  return 0;
}



