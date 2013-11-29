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

#include <petscsys.h> 
#include <petsctime.h>

#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

using namespace MoFEM;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

//Rounding
#define RND_EPS 1e-6
double roundn(double n) {

    //break n into fractional part (fract) and integral part (intp)
    double fract, intp;
    fract = modf(n,&intp);
    
    // case where n approximates zero, set n to "positive" zero
    if (abs(intp)==0) {
      if(abs(fract)<=RND_EPS) {
	n=0.000;
      }
    }

    return n;
}


int main(int argc, char *argv[]) {

  try {

  PetscInitialize(&argc,&argv,(char *)0,help);

  Core mb_instance;
  Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  /*if(rank==0) {
    EntityHandle dummy_meshset;
    rval = moab.create_meshset(MESHSET_SET,dummy_meshset); CHKERR_PETSC(rval);
  }*/

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  FieldCore core(moab);
  FieldInterface& mField = core;

  //add filds
  ierr = mField.add_field("FIELD_HDIV",Hdiv,1); CHKERRQ(ierr);
  ierr = mField.add_field("FIELD_L2",L2,3); CHKERRQ(ierr);

  //add finite elements
  ierr = mField.add_finite_element("ELEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_HDIV","FIELD_HDIV"); CHKERRQ(ierr);

  ierr = mField.add_finite_element("ELEM_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_L2","FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_L2","FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_L2","FIELD_L2"); CHKERRQ(ierr);

  ierr = mField.add_finite_element("ELEM_L2HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("ELEM_L2HDIV","FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELEM_L2HDIV","FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELEM_L2HDIV","FIELD_HDIV"); CHKERRQ(ierr);

  //add problems 
  ierr = mField.add_problem("PROBLEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROBLEM_L2"); CHKERRQ(ierr);
  ierr = mField.add_problem("PROBLEM_L2HDIV"); CHKERRQ(ierr);

  //define problems and finite elements
  ierr = mField.modify_problem_add_finite_element("PROBLEM_HDIV","ELEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROBLEM_L2","ELEM_L2"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("PROBLEM_L2HDIV","ELEM_L2HDIV"); CHKERRQ(ierr);

  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //add ents to field and set app. order
  ierr = mField.add_ents_to_field_by_TETs(0,"FIELD_HDIV"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"FIELD_HDIV",4); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"FIELD_HDIV",4); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"FIELD_L2"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"FIELD_L2",5); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELEM_HDIV",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELEM_L2",MBTET); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELEM_L2HDIV",MBTET); CHKERRQ(ierr);

  //set problem level
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_HDIV",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_L2",bit_level0); CHKERRQ(ierr);
  ierr = mField.modify_problem_ref_level_add_bit("PROBLEM_L2HDIV",bit_level0); CHKERRQ(ierr);

  //build fields
  ierr = mField.build_fields(); CHKERRQ(ierr);
  //build finite elements
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition problems
  ierr = mField.simple_partition_problem("PROBLEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROBLEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROBLEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.simple_partition_problem("PROBLEM_L2"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROBLEM_L2"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROBLEM_L2"); CHKERRQ(ierr);
  ierr = mField.compose_problem("PROBLEM_L2HDIV","PROBLEM_L2","PROBLEM_HDIV"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("PROBLEM_L2HDIV"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("PROBLEM_L2HDIV"); CHKERRQ(ierr);

  struct ApproxAnaliticalFunction {

    FieldData scalar(ublas::vector<FieldData> coords) {
      return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2];
    }
  
    ublas::vector<FieldData> dx_field;
    ublas::vector<FieldData>& dx_scalar(ublas::vector<double,ublas::bounded_array<double, 3> > coords) {
      dx_field.resize(3);
      dx_field[0] = pow(coords[0],2);
      dx_field[1] = pow(coords[1],2);
      dx_field[2] = pow(coords[2],2);
      return dx_field;
    }

  };

  struct HdivApprox: public FEMethod_UpLevelStudent,ApproxAnaliticalFunction {

    Mat A;
    Vec F;
    HdivApprox(Interface& _moab,Mat _A,Vec _F): FEMethod_UpLevelStudent(_moab),A(_A),F(_F) {}; 

    const double *G_TET_W;
    vector<double> g_NTET;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      g_NTET.resize(45*4);
      ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      G_TET_W = G_TET_W45;
      PetscFunctionReturn(0);
    }

    vector<ublas::matrix<FieldData> > rowNMatrices_Volume; //rowDiffNMatrices_Volume;
    vector<vector<DofIdx> > RowGlob_Faces;
    vector<DofIdx> RowGlob_Volume;

    vector<vector<ublas::matrix<FieldData> > > colNMatrices_Faces;
    vector<ublas::matrix<FieldData> > colNMatrices_Volume;
    vector<vector<DofIdx> > ColGlob_Faces;
    vector<DofIdx> ColGlob_Volume;

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      try {
	fe_ent_ptr = fe_ptr->fe_ptr;
	ierr = InitDataStructures(); CHKERRQ(ierr);
	ierr = GlobIndices(); CHKERRQ(ierr);
	ierr = LocalIndices(); CHKERRQ(ierr);
	ierr = DataOp(); CHKERRQ(ierr);
	ierr = ShapeFunctions_TET(g_NTET); CHKERRQ(ierr);
	ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
	ierr = GetRowNMatrix_at_GaussPoint(); CHKERRQ(ierr); //ierr = GetRowDiffNMatrix_at_GaussPoint(); CHKERRQ(ierr);
	ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);
	EntityHandle fe_handle = fe_ptr->get_ent();

	V = Shape_intVolumeMBTET(diffNTET,&*coords.data().begin()); 
	if( V <= 0 ) SETERRQ1(PETSC_COMM_SELF,1,"V < 0 for EntityHandle = %lu\n",fe_handle);

	const int g_dim = get_dim_gNTET();
	coords_at_Gauss_nodes.resize(g_dim);
	for(int gg = 0;gg<g_dim;gg++) {
	  coords_at_Gauss_nodes[gg].resize(3);
	  for(int dd = 0;dd<3;dd++) {
	    (coords_at_Gauss_nodes[gg])[dd] = cblas_ddot(4,&coords[dd],3,&get_gNTET()[gg*4],1);
	  }
	}

	ierr = GetGaussRowNMatrix("FIELD_L2",MBTET,rowNMatrices_Volume); CHKERRQ(ierr);
	//ierr = GetGaussRowDiffNMatrix("FIELD_L2",MBTET,rowDiffNMatrices_Volume); CHKERRQ(ierr);
	colNMatrices_Faces.resize(4);
	for(int ff = 0;ff<4;ff++) { //faces matrices
	  ierr = GetGaussColNMatrix("FIELD_HDIV",MBTRI,colNMatrices_Faces[ff],ff); CHKERRQ(ierr);
	}
	ierr = GetGaussColNMatrix("FIELD_HDIV",MBTET,colNMatrices_Volume); CHKERRQ(ierr);

	ierr = GetRowGlobalIndices("FIELD_L2",MBTET,RowGlob_Volume); CHKERRQ(ierr);
	ColGlob_Faces.resize(4);
	for(int ff = 0;ff<4;ff++) { //faces matrices
	  ierr = GetColGlobalIndices("FIELD_HDIV",MBTRI,ColGlob_Faces[ff],ff); CHKERRQ(ierr);
	}
	ierr = GetColGlobalIndices("FIELD_HDIV",MBTET,ColGlob_Volume); CHKERRQ(ierr);

	for(int ff = 0;ff<4;ff++) {
	  if(RowGlob_Volume.size()>0) {
	    for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) { 
	      ublas::matrix<double> NTN;
	      NTN = V*G_TET_W[gg]*prod(trans(rowNMatrices_Volume[gg]),colNMatrices_Faces[ff][gg]);
	      //NTN = V*G_TET_W[gg]*prod(trans(rowDiffNMatrices_Volume[gg]),colNMatrices_Faces[ff][gg]);
	      if(RowGlob_Volume.size()!=NTN.size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      if(ColGlob_Faces[ff].size()!=NTN.size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	      ierr = MatSetValues(A,RowGlob_Volume.size(),&(RowGlob_Volume)[0],ColGlob_Faces[ff].size(),&(ColGlob_Faces[ff])[0],&(NTN.data())[0],ADD_VALUES); CHKERRQ(ierr);
	    }
	  }
	}

	for(unsigned int gg = 0;gg<g_NTET.size()/4;gg++) { 
	  if(ColGlob_Volume.size()>0) {
	    ublas::matrix<double> NTN;
	    NTN = V*G_TET_W[gg]*prod(trans(rowNMatrices_Volume[gg]),colNMatrices_Volume[gg]);
	    //NTN = V*G_TET_W[gg]*prod(trans(rowDiffNMatrices_Volume[gg]),colNMatrices_Volume[gg]);
	    if(RowGlob_Volume.size()!=NTN.size1()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    if(ColGlob_Volume.size()!=NTN.size2()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
	    ierr = MatSetValues(A,RowGlob_Volume.size(),&(RowGlob_Volume)[0],ColGlob_Volume.size(),&(ColGlob_Volume)[0],&(NTN.data())[0],ADD_VALUES); CHKERRQ(ierr);
	  }
	  ublas::vector<FieldData> t = dx_scalar(coords_at_Gauss_nodes[gg]);
	  //ublas::vector<double> f = V*G_TET_W[gg]*prod(trans(rowDiffNMatrices_Volume[gg]),t);
	  ublas::vector<double> f = V*G_TET_W[gg]*prod(trans(rowNMatrices_Volume[gg]),t);
	  ierr = VecSetValues(F,RowGlob_Volume.size(),&(RowGlob_Volume)[0],&(f.data())[0],ADD_VALUES); CHKERRQ(ierr);
	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

  };

  struct HdivApprox_Check: public FEMethod_UpLevelStudent,ApproxAnaliticalFunction {

    HdivApprox_Check(Interface& _moab): FEMethod_UpLevelStudent(_moab) {}; 

    ofstream myfile;
    vector<double> g_NTET;
    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      //g_NTET.resize(45*4);
      //ShapeMBTET(&g_NTET[0],G_TET_X45,G_TET_Y45,G_TET_Z45,45);
      g_NTET.resize(4*4);
      ShapeMBTET(&g_NTET[0],G_TET_X4,G_TET_Y4,G_TET_Z4,4);
      myfile.open("hdiv_approximation.txt");
      PetscFunctionReturn(0);
    }

    vector<vector<ublas::matrix<FieldData> > > colNMatrices_Faces;
    vector<ublas::matrix<FieldData> > colNMatrices_Volume;
    vector<ublas::vector<FieldData> > dofsFace;
    ublas::vector<FieldData> dofsVolume;
  

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      try {
	fe_ent_ptr = fe_ptr->fe_ptr;
	ierr = InitDataStructures(); CHKERRQ(ierr);
	ierr = GlobIndices(); CHKERRQ(ierr);
	ierr = LocalIndices(); CHKERRQ(ierr);
	ierr = DataOp(); CHKERRQ(ierr);
	ierr = ShapeFunctions_TET(g_NTET); CHKERRQ(ierr);
	ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
	ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);
      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      try {

	const int g_dim = get_dim_gNTET();
	coords_at_Gauss_nodes.resize(g_dim);
	for(int gg = 0;gg<g_dim;gg++) {
	  coords_at_Gauss_nodes[gg].resize(3);
	  for(int dd = 0;dd<3;dd++) {
	    (coords_at_Gauss_nodes[gg])[dd] = cblas_ddot(4,&coords[dd],3,&get_gNTET()[gg*4],1);
	  }
	}

	colNMatrices_Faces.resize(4);
	for(int ff = 0;ff<4;ff++) { //faces matrices
	  ierr = GetGaussColNMatrix("FIELD_HDIV",MBTRI,colNMatrices_Faces[ff],ff); CHKERRQ(ierr);
	}
	ierr = GetGaussColNMatrix("FIELD_HDIV",MBTET,colNMatrices_Volume); CHKERRQ(ierr);

	dofsFace.resize(4);
	for(int ff = 0;ff<4;ff++) {
	  ierr = GetDataVector("FIELD_HDIV",MBTRI,dofsFace[ff],ff); CHKERRQ(ierr);
	}
	ierr = GetDataVector("FIELD_HDIV",MBTET,dofsVolume); CHKERRQ(ierr);

	for(int gg = 0;gg<g_dim;gg++) {
	  vector< ublas::matrix<FieldData> > field_hdiv_data;
	  ierr = GetGaussDataVector_HcurlHdiv("FIELD_HDIV",field_hdiv_data); CHKERRQ(ierr);
	  ublas::vector<FieldData> t0 = dx_scalar(coords_at_Gauss_nodes[gg]);
	  ublas::matrix_row<ublas::matrix<FieldData> > t(field_hdiv_data[gg],0);

	  ublas::vector<FieldData> t_from_NMatrix;
	  t_from_NMatrix = ublas::zero_vector<FieldData>(3);

	  for(int ff = 0;ff<4;ff++) {
	    t_from_NMatrix += prod(colNMatrices_Faces[ff][gg],dofsFace[ff]);
	  }
	  t_from_NMatrix += prod(colNMatrices_Volume[gg],dofsVolume);

	  for(int dd = 0;dd<3;dd++) {
	    cout << boost::format("%.6lf") % roundn(t0[dd]-t[dd]) << " ";
	    myfile << boost::format("%.6lf") % roundn(t0[dd]-t[dd]) << " ";
	  }
	  for(int dd = 0;dd<3;dd++) {
	    cout << boost::format("%.6lf") % roundn(t0[dd]-t_from_NMatrix[dd]) << " ";
	    myfile << boost::format("%.6lf") % roundn(t0[dd]-t_from_NMatrix[dd]) << " ";
	  }
	  for(int dd = 0;dd<3;dd++) {
	    cout << boost::format("%.6lf") % roundn(t[dd]) << " ";
	    myfile << boost::format("%.6lf") % roundn(t[dd]) << " ";
	  }
	  cout << endl;
	  myfile << endl;

	  /*cerr << endl << endl;
	  cerr << t0 << endl;
	  cerr << t << endl;
	  cerr << t-t0 << endl;
	  cerr << t_from_NMatrix-t0 << endl;*/
	  //cerr << colNMatrices_Volume[gg] << endl;

	}

      } catch (const std::exception& ex) {
	ostringstream ss;
	ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << endl;
	SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }


      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;
      myfile.close();
      PetscFunctionReturn(0);
    }

  };

  struct HdivApproxPostProc: public PostProcDisplacementsOnRefMesh,ApproxAnaliticalFunction {

    //vector<vector<Tag> > th_tags_faces;
    //vector<Tag> th_tags_volume;
    Tag th_field,th_field_err;
    HdivApproxPostProc(Interface& _moab): PostProcDisplacementsOnRefMesh(_moab) {

      double def_VAL[3] = {0,0,0};

      /*th_tags_faces.resize(4);
      for(int ff = 0; ff<4;ff++) {
	th_tags_faces[ff].resize(NBFACE_Hdiv(5));
	for(int ii = 0;ii<NBFACE_Hdiv(5);ii++) {
	  ostringstream ss;
	  ss << "HDiv_App_Face_" << ff << "_" << ii;
	  rval = moab_post_proc.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,(th_tags_faces[ff])[ii],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
	}
      }

      th_tags_volume.resize(NBVOLUME_Hdiv(5));
      for(int ii = 0;ii<NBVOLUME_Hdiv(5);ii++) {
	ostringstream ss;
	ss << "HDiv_App_Volume_" << ii;
	rval = moab_post_proc.tag_get_handle(ss.str().c_str(),3,MB_TYPE_DOUBLE,th_tags_volume[ii],MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      }*/

      rval = moab_post_proc.tag_get_handle("HDIV_APP_FIELD",3,MB_TYPE_DOUBLE,th_field,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);
      rval = moab_post_proc.tag_get_handle("HDIV_APP_FIELD_ERROR",3,MB_TYPE_DOUBLE,th_field_err,MB_TAG_CREAT|MB_TAG_SPARSE,def_VAL); CHKERR_THROW(rval);

    }; 

    vector<vector<ublas::matrix<FieldData> > > colNMatrices_Faces;
    vector<ublas::matrix<FieldData> > colNMatrices_Volume;

    PetscErrorCode do_operator() {
      PetscFunctionBegin;

      fe_ent_ptr = fe_ptr->fe_ptr;
      ierr = InitDataStructures(); CHKERRQ(ierr);
      ierr = GlobIndices(); CHKERRQ(ierr);
      ierr = LocalIndices(); CHKERRQ(ierr);
      ierr = DataOp(); CHKERRQ(ierr);
      ierr = ShapeFunctions_TET(g_NTET); CHKERRQ(ierr);
      ierr = Data_at_GaussPoints(); CHKERRQ(ierr);
      ierr = GetColNMatrix_at_GaussPoint(); CHKERRQ(ierr);

      const int g_dim = get_dim_gNTET();
      coords_at_Gauss_nodes.resize(g_dim);
      for(int gg = 0;gg<g_dim;gg++) {
	coords_at_Gauss_nodes[gg].resize(3);
	for(int dd = 0;dd<3;dd++) {
	  (coords_at_Gauss_nodes[gg])[dd] = cblas_ddot(4,&coords[dd],3,&get_gNTET()[gg*4],1);
	}
      }

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

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      ierr = do_operator(); CHKERRQ(ierr);

      colNMatrices_Faces.resize(4);
      for(int ff = 0;ff<4;ff++) { 
	ierr = GetGaussColNMatrix("FIELD_HDIV",MBTRI,colNMatrices_Faces[ff],ff); CHKERRQ(ierr);
      }
      ierr = GetGaussColNMatrix("FIELD_HDIV",MBTET,colNMatrices_Volume); CHKERRQ(ierr);

      if(colNMatrices_Volume.size()>0) {
	if(colNMatrices_Volume.size()!=node_map.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }
      for(int ff = 0;ff<4;ff++) {
	if(colNMatrices_Faces[ff].size()!=node_map.size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
      }

      vector< ublas::matrix<FieldData> > field_hdiv_data;
      ierr = GetGaussDataVector_HcurlHdiv("FIELD_HDIV",field_hdiv_data); CHKERRQ(ierr);

      int gg =0;
      map<EntityHandle,EntityHandle>::iterator mit = node_map.begin();
      for(;mit!=node_map.end();mit++,gg++) {

	/*for(int ff = 0;ff<4;ff++) {
	  for(int jj = 0;jj<colNMatrices_Faces[ff][0].size2();jj++) {
	    double t[3] = {0,0,0};
	    for(int dd = 0;dd<3;dd++) {
	      t[dd] = (colNMatrices_Faces[ff][gg])(dd,jj);
	      assert(t[dd]==t[dd]);
	     //if(fabs(t[dd])<1e-12) t[dd] = 0;
	    }
	    rval = moab_post_proc.tag_set_data(th_tags_faces[ff][jj],&mit->second,1,t); CHKERR_PETSC(rval);
	  }
	}

	if(colNMatrices_Volume.size()>0) {
	  for(int jj = 0;jj<colNMatrices_Volume[0].size2();jj++) {
	    double t[3] = {0,0,0};
	    for(int dd = 0;dd<3;dd++) {
	      t[dd] = (colNMatrices_Volume[gg])(dd,jj);
	      assert(t[dd]==t[dd]);
	    }
	    rval = moab_post_proc.tag_set_data(th_tags_volume[jj],&mit->second,1,t); CHKERR_PETSC(rval);
	  }
	}*/

	ublas::vector<FieldData> t0 = dx_scalar(coords_at_Gauss_nodes[gg]);
	ublas::matrix_row<ublas::matrix<double> > t(field_hdiv_data[gg],0);
	double _t[3],_err_t[3];
	for(int dd = 0;dd<3;dd++) {
	  _t[dd] = t[dd];
	  _err_t[dd] = t[dd]-t0[dd];
	}
	rval = moab_post_proc.tag_set_data(th_field,&mit->second,1,_t); CHKERR_PETSC(rval);
	rval = moab_post_proc.tag_set_data(th_field_err,&mit->second,1,_err_t); CHKERR_PETSC(rval);

      }

      PetscFunctionReturn(0);
    }

  };

  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("PROBLEM_L2HDIV",&Aij); CHKERRQ(ierr);
  Vec F;
  ierr = mField.VecCreateGhost("PROBLEM_L2HDIV",Row,&F); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);
  ierr = VecZeroEntries(F); CHKERRQ(ierr);

  HdivApprox fe(moab,Aij,F);
  ierr = mField.loop_finite_elements("PROBLEM_L2HDIV","ELEM_L2HDIV",fe);  CHKERRQ(ierr);

  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //Matrix View
  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //std::string wait;
  //std::cin >> wait;

  Mat B;
  ierr = MatTransposeMatMult(Aij,Aij,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B); CHKERRQ(ierr); 
  Vec AijTF,D;
  ierr = mField.VecCreateGhost("PROBLEM_HDIV",Row,&AijTF); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("PROBLEM_HDIV",Col,&D); CHKERRQ(ierr);
  ierr = MatMultTranspose(Aij,F,AijTF); CHKERRQ(ierr);
  //MatView(B,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,B,B,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,AijTF,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //ierr = VecSet(D,1); CHKERRQ(ierr);
  //Save data on mesh
  ierr = mField.set_global_VecCreateGhost("PROBLEM_HDIV",Row,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = MatDestroy(&B); CHKERRQ(ierr);
  ierr = VecDestroy(&AijTF); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);

  HdivApprox_Check fe_check(moab);
  ierr = mField.loop_finite_elements("PROBLEM_HDIV","ELEM_HDIV",fe_check);  CHKERRQ(ierr);

  HdivApproxPostProc fe_post_proc(moab);
  ierr = mField.loop_finite_elements("PROBLEM_HDIV","ELEM_HDIV",fe_post_proc);  CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    rval = fe_post_proc.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }

  PetscFinalize();
  return 0;

  } catch (const char* msg) {
    SETERRQ(PETSC_COMM_SELF,1,msg);
  } catch (const std::exception& ex) {
    ostringstream ss;
    ss << "thorw in method: " << ex.what() << endl;
    SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
  }




}
