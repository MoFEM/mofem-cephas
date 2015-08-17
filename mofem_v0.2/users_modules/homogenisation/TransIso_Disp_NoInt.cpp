/* Copyright (C) 2014, Zahur Ullah (Zahur.Ullah AT glasgow.ac.uk)
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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include <ElasticFEMethodTransIso.hpp>

#include <ElasticFE_RVELagrange_Disp.hpp>
#include <ElasticFE_RVELagrange_Homogenized_Stress_Disp.hpp>
#include <RVEVolume.hpp>

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
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

  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 1;
  }

  char outName[PETSC_MAX_PATH_LEN]="out.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_out",outName,sizeof(outName),&flg); CHKERRQ(ierr);

  char outName2[PETSC_MAX_PATH_LEN]="out_post_proc.vtk";
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_post_out",outName2,sizeof(outName2),&flg); CHKERRQ(ierr);

  //Applied strain on the RVE (vector of length 6) strain=[xx, yy, zz, xy, xz, zy]^T
  double myapplied_strain[6];
  int nmax=6;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-myapplied_strain",myapplied_strain,&nmax,&flg); CHKERRQ(ierr);
  ublas::vector<FieldData> applied_strain;
  applied_strain.resize(6);
  cblas_dcopy(6, &myapplied_strain[0], 1, &applied_strain(0), 1);
  //    cout<<"applied_strain ="<<applied_strain<<endl;

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
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  Tag th_phi;
  //    double def_val  = 0;
  rval = moab.tag_get_handle("PHI",1,MB_TYPE_DOUBLE,th_phi); CHKERR_PETSC(rval);

  Tag th_meshset_info;
  int def_meshset_info[2] = {0,0};
  rval = moab.tag_get_handle("MESHSET_INFO",2,MB_TYPE_INTEGER,th_meshset_info,MB_TAG_CREAT|MB_TAG_SPARSE,&def_meshset_info);

  int meshset_data[2];
  EntityHandle root = moab.get_root_set();
  rval = moab.tag_get_data(th_meshset_info,&root,1,meshset_data); CHKERR_PETSC(rval);

  vector<BitRefLevel> bit_levels;
  bit_levels.push_back(BitRefLevel().set(meshset_data[0]-1));

  //    const clock_t begin_time = clock();
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field.build_adjacencies(bit_levels.back()); CHKERRQ(ierr);
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  //    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl;

	//Build FE

  EntityHandle out_meshset;
  rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
  //    ierr = m_field.get_problem_finite_elements_entities("POTENTIAL_PROBLEM","POTENTIAL_ELEM",out_meshset); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_levels.back(),BitRefLevel().set(),out_meshset); CHKERRQ(ierr);
  Range LatestRefinedTets;
  rval = moab.get_entities_by_type(out_meshset, MBTET,LatestRefinedTets,true); CHKERR_PETSC(rval);

  Range LatestRefinedPrisms;
  rval = moab.get_entities_by_type(out_meshset, MBPRISM,LatestRefinedPrisms,true); CHKERR_PETSC(rval);

  cout<<"No of Prisms/Interfaces = "<<LatestRefinedPrisms.size()<<endl;

  BitRefLevel problem_bit_level = bit_levels.back();

  EntityHandle meshset_Elastic, meshset_Trans_ISO;
  rval = moab.create_meshset(MESHSET_SET,meshset_Elastic); CHKERR_PETSC(rval);
  rval = moab.create_meshset(MESHSET_SET,meshset_Trans_ISO); CHKERR_PETSC(rval);

	///Getting No. of Fibres to be used for Potential Flow Problem
	int noOfFibres=0;
	for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,BLOCKSET|UNKNOWNCUBITNAME,it)) {

		std::size_t found=it->get_name().find("PotentialFlow");
		if (found==std::string::npos) continue;
		noOfFibres += 1;
	}
	cout<<"No. of Fibres for Potential Flow : "<<noOfFibres<<endl;

	vector<int> fibreList(noOfFibres,0);
	for (int aa=0; aa<noOfFibres; aa++) {
		fibreList[aa] = aa + 1;
	}

	vector<Range> RangeFibre(noOfFibres);
	vector<EntityHandle> fibre_meshset(noOfFibres);

	for (int ii=0; ii<noOfFibres; ii++) {
		ostringstream sss;
		sss << "POTENTIAL_ELEM" << ii+1;
		for(_IT_GET_FES_BY_NAME_FOR_LOOP_(m_field, sss.str().c_str() ,it)){
			RangeFibre[ii].insert(it->get_ent());
			rval = moab.create_meshset(MESHSET_SET,fibre_meshset[ii]); CHKERR_PETSC(rval);
			rval = moab.add_entities(fibre_meshset[ii],RangeFibre[ii]); CHKERR_PETSC(rval);
			rval = moab.unite_meshset(meshset_Trans_ISO,fibre_meshset[ii]); CHKERR_PETSC(rval);
		}
	}

  rval = moab.write_file("meshset_Trans_ISO.vtk","VTK","",&meshset_Trans_ISO,1); CHKERR_PETSC(rval);

	for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)){

		if(it->get_name() == "MAT_ELASTIC_1") {
			Range TetsInBlock;
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock,true); CHKERR_PETSC(rval);
			Range block_rope_bit_level = intersect(LatestRefinedTets,TetsInBlock);

      cout<<"=============  TetsInBlock  "<< TetsInBlock.size() <<endl;

			rval = moab.add_entities(meshset_Elastic,block_rope_bit_level);CHKERR_PETSC(rval);

		}
	}
  ierr = m_field.seed_finite_elements(meshset_Elastic); CHKERRQ(ierr);



	Range prims_on_problem_bit_level;
	ierr = m_field.get_entities_by_type_and_ref_level(problem_bit_level,BitRefLevel().set(),MBPRISM,prims_on_problem_bit_level); CHKERRQ(ierr);
  //to create meshset from range
  EntityHandle meshset_prims_on_problem_bit_level;
  rval = moab.create_meshset(MESHSET_SET,meshset_prims_on_problem_bit_level); CHKERR_PETSC(rval);
	rval = moab.add_entities(meshset_prims_on_problem_bit_level,prims_on_problem_bit_level); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_MESHSET(meshset_prims_on_problem_bit_level,BitRefLevel().set()); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
  int field_rank=3;
  ierr = m_field.add_field("DISPLACEMENT",H1,field_rank); CHKERRQ(ierr);
  ierr = m_field.add_field("Lagrange_mul_disp",H1,field_rank); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("Lagrange_elem"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  //FE Transverse Isotropic
  ierr = m_field.modify_finite_element_add_field_row("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","DISPLACEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TRAN_ISOTROPIC_ELASTIC","POTENTIAL_FIELD"); CHKERRQ(ierr);


  //C row as Lagrange_mul_disp and col as DISPLACEMENT
  ierr = m_field.modify_finite_element_add_field_row("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);

  //CT col as Lagrange_mul_disp and row as DISPLACEMENT
  ierr = m_field.modify_finite_element_add_field_col("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);

  //As for stress we need both displacement and temprature (Lukasz)
  ierr = m_field.modify_finite_element_add_field_data("Lagrange_elem","Lagrange_mul_disp"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("Lagrange_elem","DISPLACEMENT"); CHKERRQ(ierr);


  /***/
  //Declare problem

  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = m_field.add_ents_to_finite_element_by_TETs(meshset_Elastic,"ELASTIC",true); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TETs(meshset_Trans_ISO,"TRAN_ISOTROPIC_ELASTIC",true); CHKERRQ(ierr);

  Range SurfacesFaces;
  ierr = m_field.get_cubit_msId_entities_by_dimension(103,SIDESET,2,SurfacesFaces,true); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of SIDESET 103 = %d\n",SurfacesFaces.size()); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_TRIs(SurfacesFaces,"Lagrange_elem"); CHKERRQ(ierr);

  //to create meshset from range
  EntityHandle BoundFacesMeshset;
  rval = moab.create_meshset(MESHSET_SET,BoundFacesMeshset); CHKERR_PETSC(rval);
	rval = moab.add_entities(BoundFacesMeshset,SurfacesFaces); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_MESHSET(BoundFacesMeshset,BitRefLevel().set()); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_field_by_TRIs(BoundFacesMeshset,"Lagrange_mul_disp",2); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTRI,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"Lagrange_mul_disp",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"Lagrange_mul_disp",1); CHKERRQ(ierr);

  /****/
  //build database

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(problem_bit_level); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","Lagrange_elem"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",problem_bit_level); CHKERRQ(ierr);
  //build problem
  ierr = m_field.build_problems(); CHKERRQ(ierr);

  /****/
  //mesh partitioning

  //partition
  ierr = m_field.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F,D;
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",ROW,&F); CHKERRQ(ierr);
  ierr = m_field.VecCreateGhost("ELASTIC_MECHANICS",COL,&D); CHKERRQ(ierr);

  Mat Aij;
  ierr = m_field.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  struct MyElasticFEMethod: public ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _m_field,
                      Mat &_Aij,Vec &_D,Vec& _F,double _lambda,double _mu):
    ElasticFEMethod(_m_field,_Aij,_D,_F,_lambda,_mu) {};

    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
  };


  for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"POTENTIAL_FIELD",dof_ptr)) {
    if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
    EntityHandle ent = dof_ptr->get_ent();
    double &fval = dof_ptr->get_FieldData();
    double phi;
    rval = moab.tag_get_data(th_phi,&ent,1,&phi); CHKERR_PETSC(rval);
    fval = phi;
  }

  //Assemble F and Aij
//  double YoungModulusP;
//  double PoissonRatioP;
//  double YoungModulusZ;
//  double PoissonRatioPZ;
//  double ShearModulusZP;
  double YoungModulus;
  double PoissonRatio;
  double alpha;

  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it))
  {
    cout << endl << *it << endl;

    //Get block name
    string name = it->get_name();

//    if (name.compare(0,20,"MAT_ELASTIC_TRANSISO") == 0)
//    {
//      Mat_Elastic_TransIso mydata;
//      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
//      cout << mydata;
//      YoungModulusP=mydata.data.Youngp;
//      YoungModulusZ=mydata.data.Youngz;
//      PoissonRatioP=mydata.data.Poissonp;
//      PoissonRatioPZ=mydata.data.Poissonpz;
//      if (mydata.data.Shearzp!=0) {
//        ShearModulusZP=mydata.data.Shearzp;
//      }else{
//        ShearModulusZP=YoungModulusZ/(2*(1+PoissonRatioPZ));}
//    }
    if (name.compare(0,11,"MAT_ELASTIC") == 0)
    {
      Mat_Elastic mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      YoungModulus=mydata.data.Young;
      PoissonRatio=mydata.data.Poisson;
    }
    else if (name.compare(0,10,"MAT_INTERF") == 0)
    {
      Mat_Interf mydata;
      ierr = it->get_attribute_data_structure(mydata); CHKERRQ(ierr);
      cout << mydata;
      alpha = mydata.data.alpha;
    }
  }

  //alpha = 1e8;
  cout<<"alpha   = "<<alpha<<endl;
  MyElasticFEMethod MyFE(m_field,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
  TranIsotropicFibreDirRotElasticFEMethod MyTIsotFE(m_field,Aij,D,F);
  ElasticFE_RVELagrange_Disp MyFE_RVELagrange(m_field,Aij,D,F,applied_strain,"DISPLACEMENT","Lagrange_mul_disp",field_rank);


  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyFE);  CHKERRQ(ierr);
	ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyTIsotFE);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVELagrange);  CHKERRQ(ierr);

  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);


  //    //Matrix View
  //    MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
  //    std::string wait;
  //    std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);

  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = m_field.set_global_ghost_vector("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);


  //Calculation of Homogenized stress
  //====================================================================================================================================
  double RVE_volume;    RVE_volume=0.0;  //RVE volume for full RVE We need this for stress calculation
  Vec RVE_volume_Vec;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 1, pcomm->size(), &RVE_volume_Vec);  CHKERRQ(ierr);
  ierr = VecZeroEntries(RVE_volume_Vec); CHKERRQ(ierr);

  RVEVolume MyRVEVol(m_field,Aij,D,F,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio), RVE_volume_Vec);
  RVEVolumeTrans MyRVEVolTrans(m_field,Aij,D,F, RVE_volume_Vec);

  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",MyRVEVol);  CHKERRQ(ierr);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",MyRVEVolTrans);  CHKERRQ(ierr);

  ierr = VecView(RVE_volume_Vec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = VecSum(RVE_volume_Vec, &RVE_volume);  CHKERRQ(ierr);
  cout<<"Final RVE_volume = "<< RVE_volume <<endl;
  cout<<"Actual RVE_volume = "<< 3*0.3*0.78<<endl;  //Lx=3, Ly=0.3; Lz=0.78

  //create a vector for 6 components of homogenized stress
  Vec Stress_Homo;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, 6, 6*pcomm->size(), &Stress_Homo);  CHKERRQ(ierr);
  ierr = VecZeroEntries(Stress_Homo); CHKERRQ(ierr);

  //    ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ElasticFE_RVELagrange_Homogenized_Stress_Disp MyFE_RVEHomoStressDisp(m_field,Aij,D,F,&RVE_volume, applied_strain, Stress_Homo,"DISPLACEMENT","Lagrange_mul_disp",field_rank);
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","Lagrange_elem",MyFE_RVEHomoStressDisp);  CHKERRQ(ierr);

  if(pcomm->rank()) cout<< " Stress_Homo =  "<<endl;
  ierr = VecView(Stress_Homo,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //====================================================================================================================================


  PostProcVertexMethod ent_method(moab);
  ierr = m_field.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = m_field.get_problem_finite_elements_entities("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = m_field.get_problem_finite_elements_entities("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",out_meshset); CHKERRQ(ierr);
    ierr = m_field.get_problem_finite_elements_entities("ELASTIC_MECHANICS","INTERFACE",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file(outName,"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  //    TranIso_PostProc_FibreDirRot_OnRefMesh fe_post_proc_method( m_field, LAMBDA(YoungModulusP,PoissonRatioP),MU(YoungModulusP,PoissonRatioP), YoungModulusP,YoungModulusZ,PoissonRatioP,PoissonRatioPZ,ShearModulusZP);
  //
  //    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  //    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","TRAN_ISOTROPIC_ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  //
  //    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //    if(pcomm->rank()==0) {
  //        rval = fe_post_proc_method.moab_post_proc.write_file(outName2,"VTK",""); CHKERR_PETSC(rval);
  //    }
  //
  //    PostProcCohesiveForces fe_post_proc_prisms(m_field,YoungModulus*alpha);
  //    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","INTERFACE",fe_post_proc_prisms);  CHKERRQ(ierr);
  //    PetscSynchronizedFlush(PETSC_COMM_WORLD);
  //    if(pcomm->rank()==0) {
  //        rval = fe_post_proc_prisms.moab_post_proc.write_file("out_post_proc_prisms.vtk","VTK",""); CHKERR_PETSC(rval);
  //    }


  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  ierr = VecDestroy(&RVE_volume_Vec); CHKERRQ(ierr);
  ierr = VecDestroy(&Stress_Homo); CHKERRQ(ierr);


  ierr = PetscTime(&v2);CHKERRQ(ierr);
  ierr = PetscGetCPUTime(&t2);CHKERRQ(ierr);

  PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Total Rank %d Time = %f CPU Time = %f\n",pcomm->rank(),v2-v1,t2-t1);

  PetscFinalize();

}
