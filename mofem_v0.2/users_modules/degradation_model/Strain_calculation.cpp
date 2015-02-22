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

#include <MoFEM.hpp>
using namespace MoFEM;

#include <DirichletBC.hpp>

#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <FluidPressure.hpp>
#include <BodyForce.hpp>
#include <ThermalStressElement.hpp>

#include <FEMethod_LowLevelStudent.hpp>
#include <FEMethod_UpLevelStudent.hpp>

#include <PostProcVertexMethod.hpp>
#include <PostProcDisplacementAndStrainOnRefindedMesh.hpp>

#include <ElasticFEMethod.hpp>
#include <ElasticFEMethod_strain_calculation.hpp>



using namespace boost::numeric;
using namespace ObosleteUsersModules;


ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";
int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,PETSC_NULL,help);

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
    order = 5;
  }
    

  //Read mesh to MOAB
  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);

  // stl::bitset see for more details
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  /***/
  //Define problem

  //Fields
//  ierr = m_field.add_field("DISPLACEMENT",H1,3); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("ELASTIC","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("ELASTIC","DISP_MACRO"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DISP_MACRO"); CHKERRQ(ierr);

  //define problems
  ierr = m_field.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  /***/
  //Declare problem


  //add entitities (by tets) to the field
//  ierr = m_field.add_ents_to_field_by_TETs(0,"DISPLACEMENT"); CHKERRQ(ierr);

  //add finite elements entities
  Range TetsInBlock_Ele_Strain_RVE;
  for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)){
		if(it->get_Cubit_name() == "Ele_Strain_RVE") {
			rval = moab.get_entities_by_type(it->meshset, MBTET,TetsInBlock_Ele_Strain_RVE,true); CHKERR_PETSC(rval);
//      cout<<"TetsInBlock in Ele_Strain_RVE   "<<TetsInBlock<<endl;
		}
	}

//  ierr = m_field.add_ents_to_finite_element_by_TETs(TetsInBlock_Ele_Strain_RVE,"ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

//  //set app. order
//  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
//  //int order = 5;
//  ierr = m_field.set_field_order(0,MBTET,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = m_field.set_field_order(0,MBTRI,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = m_field.set_field_order(0,MBEDGE,"DISPLACEMENT",order); CHKERRQ(ierr);
//  ierr = m_field.set_field_order(0,MBVERTEX,"DISPLACEMENT",1); CHKERRQ(ierr);
//
//  ierr = MetaNeummanForces::addNeumannBCElements(m_field,"DISPLACEMENT"); CHKERRQ(ierr);
//  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","FORCE_FE"); CHKERRQ(ierr);
//  ierr = m_field.modify_problem_add_finite_element("ELASTIC_MECHANICS","PRESSURE_FE"); CHKERRQ(ierr);

  /****/
  //build database

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

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


  //Assemble F and Aij
  const double young_modulus = 1;
  const double poisson_ratio = 0.0;
  DisplacementBCFEMethodPreAndPostProc my_dirichlet_bc(m_field,"DISP_MACRO",Aij,D,F);
  ElasticFEMethod_strain_calculation fe(m_field,Aij,D,F,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio),"DISP_MACRO");

  ierr = VecZeroEntries(F); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(F,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = MatZeroEntries(Aij); CHKERRQ(ierr);

  //preproc
  ierr = m_field.problem_basic_method_preProcess("ELASTIC_MECHANICS",my_dirichlet_bc); CHKERRQ(ierr);
  //loop elements
  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe);  CHKERRQ(ierr);
//  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
//  ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field,neumann_forces,F,"DISPLACEMENT"); CHKERRQ(ierr);
//  boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
//  for(;mit!=neumann_forces.end();mit++) {
//    ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS",mit->first,mit->second->getLoopFe()); CHKERRQ(ierr);
//  }
//  //postproc
//  ierr = m_field.problem_basic_method_postProcess("ELASTIC_MECHANICS",my_dirichlet_bc); CHKERRQ(ierr);
//
//  ierr = VecGhostUpdateBegin(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(F,ADD_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
//  ierr = VecAssemblyEnd(F); CHKERRQ(ierr);
//  ierr = MatAssemblyBegin(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//  ierr = MatAssemblyEnd(Aij,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//
//  //Matrix View
//  //MatView(Aij,PETSC_VIEWER_DRAW_WORLD);//PETSC_VIEWER_STDOUT_WORLD);
//  //std::string wait;
//  //std::cin >> wait;
//
//  //Solver
//  KSP solver;
//  ierr = KSPCreate(PETSC_COMM_WORLD,&solver); CHKERRQ(ierr);
//  ierr = KSPSetOperators(solver,Aij,Aij); CHKERRQ(ierr);
//  ierr = KSPSetFromOptions(solver); CHKERRQ(ierr);
//  ierr = KSPSetUp(solver); CHKERRQ(ierr);
//
//  ierr = KSPSolve(solver,F,D); CHKERRQ(ierr);
//  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
//
//  //Save data on mesh
//  ierr = m_field.set_global_VecCreateGhost("ELASTIC_MECHANICS",ROW,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
//  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
//
//  PostProcVertexMethod ent_method(moab);
//  ierr = m_field.loop_dofs("ELASTIC_MECHANICS","DISPLACEMENT",ROW,ent_method); CHKERRQ(ierr);

//  if(pcomm->rank()==0) {
//    EntityHandle out_meshset;
//    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
//    ierr = m_field.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
//    //rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//      rval = moab.write_file((string(mesh_file_name)+".vtk").c_str(),"VTK","",&out_meshset,1); CHKERR_PETSC(rval);
//    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
//  }

  //PostProcDisplacemenysAndStarinOnRefMesh fe_post_proc_method(moab);
//  PostProcDisplacemenysAndStarinAndElasticLinearStressOnRefMesh fe_post_proc_method(m_field,LAMBDA(YoungModulus,PoissonRatio),MU(YoungModulus,PoissonRatio));
//  ierr = m_field.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);

//  PetscSynchronizedFlush(PETSC_COMM_WORLD);
//  if(pcomm->rank()==0) {
//    //rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
//      rval = fe_post_proc_method.moab_post_proc.write_file((string(mesh_file_name)+".out.vtk").c_str(),"VTK",""); CHKERR_PETSC(rval);
//  }
//
//    //Open mesh_file_name.txt for writing
//    ofstream myfile;
//    myfile.open ((string(mesh_file_name)+".txt").c_str());
//    
//    //Output displacements
//    cout << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
//    myfile << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
//    
//    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field,"DISPLACEMENT",dof_ptr))
//    {
//        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
//        
//        if(dof_ptr->get_dof_rank()==0)
//        {
//            //Round and truncate to 3 decimal places
//            double fval = dof_ptr->get_FieldData();
//            cout << boost::format("%.3lf") % roundn(fval) << "  ";
//            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
//        }
//        if(dof_ptr->get_dof_rank()==1)
//        {
//            //Round and truncate to 3 decimal places
//            double fval = dof_ptr->get_FieldData();
//            cout << boost::format("%.3lf") % roundn(fval) << "  ";
//            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
//        }
//        if(dof_ptr->get_dof_rank()==2)
//        {
//            //Round and truncate to 3 decimal places
//            double fval = dof_ptr->get_FieldData();
//            cout << boost::format("%.3lf") % roundn(fval) << endl;
//            myfile << boost::format("%.3lf") % roundn(fval) << endl;
//        }
//        
//    }
//    
//    //Close mesh_file_name.txt
//    myfile.close();
//
//  //destroy matrices
//  ierr = VecDestroy(&F); CHKERRQ(ierr);
//  ierr = VecDestroy(&D); CHKERRQ(ierr);
//  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
//  ierr = KSPDestroy(&solver); CHKERRQ(ierr);

  PetscFinalize();
  return 0;

}

