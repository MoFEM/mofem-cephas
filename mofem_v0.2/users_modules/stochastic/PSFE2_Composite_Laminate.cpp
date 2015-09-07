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

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>
#include <petsctime.h>

#include <SurfacePressure.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

//#include <PostProcOnRefMesh.hpp>
//#include <PostProcHookStresses.hpp>
//#include <PostProcHookStressesLaminates.hpp>

#include <ElasticFEMethod.hpp>
#include <FE2_ElasticFEMethod.hpp>

using namespace boost::numeric;
using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

const double young_modulus = 1;
const double poisson_ratio = 0.0;


PetscErrorCode Dmat_Transfermation(double theta, ublas::matrix<FieldData> &Dmat_in, ublas::matrix<FieldData> &Dmat_out) {
  PetscFunctionBegin;
  
  double l1, l2, l3, m1, m2, m3, n1, n2, n3;
  l1=cos(theta);
  m1=sin(theta);
  n1=0.0;
  l2=-sin(theta);
  m2=cos(theta);
  n2=0.0;
  l3=0.0;
  m3=0.0;
  n3=1.0;
  
  //  cout<<"l1 = "<<l1<<endl;
  //  cout<<"m1 = "<<m1<<endl;
  //  cout<<"l2 = "<<l2<<endl;
  //  cout<<"m2 = "<<m2<<endl;
  
  ublas::matrix<FieldData> T_strain;    T_strain.resize(6,6);
  T_strain(0,0)=l1*l1;    T_strain(0,1)=m1*m1;   T_strain(0,2)=n1*n1;    T_strain(0,3)=l1*m1;        T_strain(0,4)=l1*n1;        T_strain(0,5)=m1*n1;
  T_strain(1,0)=l2*l2;    T_strain(1,1)=m2*m2;   T_strain(1,2)=n2*n2;    T_strain(1,3)=l2*m2;        T_strain(1,4)=l2*n2;        T_strain(1,5)=m2*n2;
  T_strain(2,0)=l3*l3;    T_strain(2,1)=m3*m3;   T_strain(2,2)=n3*n3;    T_strain(2,3)=l3*m3;        T_strain(2,4)=l3*n3;        T_strain(2,5)=m3*n3;
  T_strain(3,0)=2*l1*l2;  T_strain(3,1)=2*m1*m2; T_strain(3,2)=2*n1*n2;  T_strain(3,3)=l1*m2+m1*l2;  T_strain(3,4)=l1*n2+n1*l2;  T_strain(3,5)=m1*n2+n1*m2;
  T_strain(4,0)=2*l1*l3;  T_strain(4,1)=2*m1*m3; T_strain(4,2)=2*n1*n3;  T_strain(4,3)=l1*m3+m1*l3;  T_strain(4,4)=l1*n3+n1*l3;  T_strain(4,5)=m1*n3+n1*m3;
  T_strain(5,0)=2*l2*l3;  T_strain(5,1)=2*m2*m3; T_strain(5,2)=2*n2*n3;  T_strain(5,3)=l2*m3+m2*l3;  T_strain(5,4)=l2*n3+n2*l3;  T_strain(5,5)=m2*n3+n2*m3;
//  cout<<"\n\nT_strain = "<<T_strain<<endl;
  
  ublas::matrix<FieldData> Mat1=prod(Dmat_in,T_strain);
  Dmat_out=prod(trans(T_strain), Mat1);
//  cout<<"\n\n Dmat_out = "<<Dmat_out<<endl;

  PetscFunctionReturn(0);
}

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

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  
  
  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& mField = core;
  
  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  
  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);
  
  //Read the saved Dmat mechancial (from the computational homgenisaiton of the 0deg RVE)
  ublas::matrix<FieldData> Dmat_0deg;
  Dmat_0deg.resize(6,6);  Dmat_0deg.clear();

  cout<<"Dmat_0deg Before =  "<<Dmat_0deg<<endl;
  int fd;
  PetscViewer view_in;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Mechanical_Dmat.out",FILE_MODE_READ,&view_in);
  PetscViewerBinaryGetDescriptor(view_in,&fd);
  PetscBinaryRead(fd,&Dmat_0deg(0,0),36,PETSC_DOUBLE);
  PetscViewerDestroy(&view_in);
  cout<<"Dmat_0deg After =   "<<Dmat_0deg<<endl;

  //Find Dmat_90deg by rotating Dmat_0deg at an angle of +90 deg about the z-axis
  double theta; theta=90*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_90deg;  Dmat_90deg.resize(6,6);   Dmat_90deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_90deg); CHKERRQ(ierr);
  cout<<"\n\nDmat_90deg = "<<Dmat_90deg<<endl;

  //Find Dmat_pos25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
  theta=25*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_pos25deg;  Dmat_pos25deg.resize(6,6);   Dmat_pos25deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_pos25deg); CHKERRQ(ierr);
  cout<<"\n\n Dmat_pos25deg = "<<Dmat_pos25deg<<endl;

  //Find Dmat_neg25deg by rotating Dmat_pos25deg at an angle of +25 deg about the z-axis
  theta=-25*(M_PI/180.0); //rotation angle about the Z-axis
  ublas::matrix<FieldData> Dmat_neg25deg;  Dmat_neg25deg.resize(6,6);   Dmat_neg25deg.clear();
  ierr = Dmat_Transfermation(theta, Dmat_0deg, Dmat_neg25deg); CHKERRQ(ierr);
  cout<<"\n\n Dmat_neg25deg = "<<Dmat_neg25deg<<endl;

  
  //Define problem
  
  //Fields
  ierr = mField.add_field("DISP",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  
  //FE
  ierr = mField.add_finite_element("ELASTIC_90deg",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("ELASTIC_pos25deg",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("ELASTIC_neg25deg",MF_ZERO); CHKERRQ(ierr);
  ierr = mField.add_finite_element("ELASTIC_0_90_post_process",MF_ZERO); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_90deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_90deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_90deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_90deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_pos25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_pos25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_pos25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_pos25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_neg25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_neg25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_neg25deg","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_neg25deg","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  
  
  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC_0_90_post_process","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC_0_90_post_process","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","DISP"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC_0_90_post_process","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  
  //set finite elements for problem
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC_neg25deg"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","ELASTIC_0_90_post_process"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_PROB",bit_level0); CHKERRQ(ierr);
  
  //Declare problem
  
  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"DISP"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);

  Range TetsInBlock_90deg, TetsInBlock_pos25deg, TetsInBlock_neg25deg;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(mField,BLOCKSET,it)){
    if(it->get_name() == "MAT_ELASTIC_90deg") {
      rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_90deg,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Pos25deg") {
      rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_pos25deg,true); CHKERR_PETSC(rval);
    }
    if(it->get_name() == "MAT_ELASTIC_Neg25deg") {
      rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_neg25deg,true); CHKERR_PETSC(rval);
    }
  }
  cout<<"===============================   TetsInBlock_90deg "<<TetsInBlock_90deg<<endl;
  cout<<"===============================   TetsInBlock_pos25deg "<<TetsInBlock_pos25deg<<endl;
  cout<<"===============================   TetsInBlock_neg25deg "<<TetsInBlock_neg25deg<<endl;
  
  
  ierr = mField.add_ents_to_finite_element_by_TETs(TetsInBlock_90deg, "ELASTIC_90deg"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(TetsInBlock_pos25deg,"ELASTIC_pos25deg"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(TetsInBlock_neg25deg,"ELASTIC_neg25deg"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_finite_element_by_TETs(0, "ELASTIC_0_90_post_process"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  //int order = 5;
  ierr = mField.set_field_order(0,MBTET,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"DISP",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"DISP",1); CHKERRQ(ierr);
  //
  ierr = mField.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);
  
  ierr = MetaNeummanForces::addNeumannBCElements(mField,"DISP"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_PROB","FORCE_FE"); CHKERRQ(ierr);
  
  //build database
  
  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);
  
  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);
  
  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);
  
  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);
  
  //mesh partitioning
  
  //partition
  ierr = mField.partition_problem("ELASTIC_PROB"); CHKERRQ(ierr);
  //PetscBarrier(PETSC_NULL);
  ierr = mField.partition_finite_elements("ELASTIC_PROB"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = mField.partition_ghost_dofs("ELASTIC_PROB"); CHKERRQ(ierr);
  
  //mField.list_dofs_by_field_name("DISP",true);
  
  //print bcs
  ierr = mField.print_cubit_displacement_set(); CHKERRQ(ierr);
  ierr = mField.print_cubit_force_set(); CHKERRQ(ierr);
  //print block sets with materials
  ierr = mField.print_cubit_materials_set(); CHKERRQ(ierr);
  
  //create matrices
  Vec F,D;
  ierr = mField.VecCreateGhost("ELASTIC_PROB",ROW,&F); CHKERRQ(ierr);
  ierr = mField.VecCreateGhost("ELASTIC_PROB",COL,&D); CHKERRQ(ierr);
  
  Mat A;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_PROB",&A); CHKERRQ(ierr);
  
  
  struct MyElasticFEMethod: public FE2_ElasticFEMethod {
    MyElasticFEMethod(FieldInterface& _mField,Mat _A,Vec _D,Vec& _F, ublas::matrix<FieldData> _Dmat,string _field_name):
    FE2_ElasticFEMethod(_mField,_A,_D,_F, _Dmat, _field_name) {};
    
    PetscErrorCode Fint(Vec F_int) {
      PetscFunctionBegin;
      ierr = FE2_ElasticFEMethod::Fint(); CHKERRQ(ierr);
      for(int rr = 0;rr<row_mat;rr++) {
        if(RowGlob[rr].size()!=f_int[rr].size()) SETERRQ(PETSC_COMM_SELF,1,"data inconsistency");
        if(RowGlob[rr].size()==0) continue;
        f_int[rr] *= -1; //This is not SNES we solve K*D = -RES
        ierr = VecSetValues(F_int,RowGlob[rr].size(),&(RowGlob[rr])[0],&(f_int[rr].data()[0]),ADD_VALUES); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }
    
  };
  
  
  Projection10NodeCoordsOnField ent_method_material(mField,"MESH_NODE_POSITIONS");
  ierr = mField.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
  
  //Assemble F and A
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(mField,"DISP",A,D,F);
  //preproc
  ierr = mField.problem_basic_method_preProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
  ierr = mField.set_global_ghost_vector("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
 
  
  MyElasticFEMethod my_fe_90deg   (mField,A,D,F,Dmat_90deg,   "DISP");
  MyElasticFEMethod my_fe_pos25deg(mField,A,D,F,Dmat_pos25deg,"DISP");
  MyElasticFEMethod my_fe_neg25deg(mField,A,D,F,Dmat_neg25deg,"DISP");

  
  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  
  //loop elems
  //PetscBarrier(PETSC_NULL);
  ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC_90deg",my_fe_90deg);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC_pos25deg",my_fe_pos25deg);  CHKERRQ(ierr);
  ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC_neg25deg",my_fe_neg25deg);  CHKERRQ(ierr);

  //forces and preassures on surface
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(mField,neumann_forces,F,"DISP"); CHKERRQ(ierr);
  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
  for(;mit!=neumann_forces.end();mit++) {
    ierr = mField.loop_finite_elements("ELASTIC_PROB",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
  }
  
  //postproc
  ierr = mField.problem_basic_method_postProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
  
  
  //set matrix possitives define and symetric for cholesky and icc preceonditionser
  ierr = MatSetOption(A,MAT_SPD,PETSC_TRUE); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

//  Matrix View
//  MatView(A,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  std::string wait;
//  std::cin >> wait;

  //Solver
  KSP solver;
  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
  ierr = KSPSetOperators(solver,A,A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
  ierr = KSPSetUp(solver); CHKERRQ(ierr);
  
  // elastic analys
  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = mField.problem_basic_method_preProcess("ELASTIC_PROB",my_dirichlet_bc); CHKERRQ(ierr);
//  ierr = VecView(D,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  //Save data on mesh
  ierr = mField.set_global_ghost_vector("ELASTIC_PROB",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

  /*PostPocOnRefinedMesh post_proc(mField);
  ierr = post_proc.generateReferenceElementMesh(); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("DISP"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = post_proc.addFieldValuesGradientPostProc("DISP"); CHKERRQ(ierr);
 
  //add postpocessing for sresses
  post_proc.getOpPtrVector().push_back(
    new PostPorcStressLaminates(
      mField,
      post_proc.postProcMesh,
      post_proc.mapGaussPts,
      "DISP",
      post_proc.commonData,Dmat_neg25deg)
  );

  
  ierr = mField.loop_finite_elements("ELASTIC_PROB","ELASTIC_neg25deg",post_proc); CHKERRQ(ierr);
  rval = post_proc.postProcMesh.write_file("out_macro.h5m","MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
*/
  
  //Destroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = KSPDestroy(&solver); CHKERRQ(ierr);
  PetscFinalize();
  
}

