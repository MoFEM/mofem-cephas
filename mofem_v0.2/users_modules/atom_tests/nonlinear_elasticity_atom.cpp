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

extern "C" {
  #include <complex_for_lazy.h>
}

#include <ArcLengthTools.hpp>
#include <FEMethod_ComplexForLazy.hpp>
#include <FEMethod_DriverComplexForLazy.hpp>

#include <SurfacePressureComplexForLazy.hpp>
#include <PostProcNonLinearElasticityStresseOnRefindedMesh.hpp>

using namespace ObosleteUsersModules;

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


struct NL_ElasticFEMethod: public NonLinearSpatialElasticFEMthod {

  NL_ElasticFEMethod(FieldInterface& _mField,double _lambda,double _mu,int _verbose = 0): 
      FEMethod_ComplexForLazy_Data(_mField,_verbose), 
      NonLinearSpatialElasticFEMthod(_mField,_lambda,_mu,_verbose)  {
    set_PhysicalEquationNumber(hooke);
  }

};

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
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

  MoFEM::Core core(moab);
  FieldInterface& mField = core;

  Range CubitSIDESETs_meshsets;
  ierr = mField.get_Cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = mField.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = mField.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = mField.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = mField.add_field("SPATIAL_POSITION",H1,3); CHKERRQ(ierr);

  //FE
  ierr = mField.add_finite_element("ELASTIC"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = mField.modify_finite_element_add_field_row("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);

  //define problems
  ierr = mField.add_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //set finite elements for problems
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","ELASTIC"); CHKERRQ(ierr);

  //set refinment level for problem
  ierr = mField.modify_problem_ref_level_add_bit("ELASTIC_MECHANICS",bit_level0); CHKERRQ(ierr);

  //add entitities (by tets) to the field
  ierr = mField.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

  //add finite elements entities
  ierr = mField.add_ents_to_finite_element_EntType_by_bit_ref(bit_level0,"ELASTIC",MBTET); CHKERRQ(ierr);

  //set app. order
  PetscInt order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    order = 4;
  }
  ierr = mField.set_field_order(0,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  ierr = mField.add_finite_element("NEUAMNN_FE"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_row("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_col("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("NEUAMNN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = mField.modify_problem_add_finite_element("ELASTIC_MECHANICS","NEUAMNN_FE"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = mField.add_ents_to_finite_element_by_TRIs(tris,"NEUAMNN_FE"); CHKERRQ(ierr);
  }

  #ifdef NONLINEAR_TEMPERATUTE

  ierr = mField.add_field("TEMPERATURE",H1,1); CHKERRQ(ierr);
  ierr = mField.modify_finite_element_add_field_data("ELASTIC","TEMPERATURE"); CHKERRQ(ierr);
  ierr = mField.add_ents_to_field_by_TETs(0,"TEMPERATURE"); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTET,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBTRI,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBEDGE,"TEMPERATURE",1); CHKERRQ(ierr);
  ierr = mField.set_field_order(0,MBVERTEX,"TEMPERATURE",1); CHKERRQ(ierr);

  #endif //NONLINEAR_TEMPERATUTE

  //build field
  ierr = mField.build_fields(); CHKERRQ(ierr);

  {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr)) {
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

  #ifdef NONLINEAR_TEMPERATUTE

  {
    EntityHandle node = 0;
    double coords[3];
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"TEMPERATURE",dof_ptr)) {
      if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
      EntityHandle ent = dof_ptr->get_ent();
      if(node!=ent) {
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
	node = ent;
      }
      double &fval = dof_ptr->get_FieldData();
      fval = coords[0];
    }
  }

  #endif //NONLINEAR_TEMPERATUTE

  //build finite elemnts
  ierr = mField.build_finite_elements(); CHKERRQ(ierr);

  //build adjacencies
  ierr = mField.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //build problem
  ierr = mField.build_problems(); CHKERRQ(ierr);

  //partition
  ierr = mField.partition_problem("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_finite_elements("ELASTIC_MECHANICS"); CHKERRQ(ierr);
  ierr = mField.partition_ghost_dofs("ELASTIC_MECHANICS"); CHKERRQ(ierr);

  //create matrices
  Vec F;
  ierr = mField.VecCreateGhost("ELASTIC_MECHANICS",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);
  Mat Aij;
  ierr = mField.MatCreateMPIAIJWithArrays("ELASTIC_MECHANICS",&Aij); CHKERRQ(ierr);

  const double young_modulus = 1.;
  const double poisson_ratio = 0.;
  NL_ElasticFEMethod my_fe(mField,LAMBDA(young_modulus,poisson_ratio),MU(young_modulus,poisson_ratio));
  #ifdef NONLINEAR_TEMPERATUTE
  my_fe.thermal_expansion = 0.1;
  #endif //NONLINEAR_TEMPERATUTE

  NeummanForcesSurfaceComplexForLazy neumann_forces(mField,Aij,F);
  NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_spatial = neumann_forces.getLoopSpatialFe();
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,NODESET|FORCESET,it)) {
    ierr = fe_spatial.addForce(it->get_msId()); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(mField,SIDESET|PRESSURESET,it)) {
    ierr = fe_spatial.addPreassure(it->get_msId()); CHKERRQ(ierr);
  }

  SpatialPositionsBCFEMethodPreAndPostProc MyDirichletBC(mField,"SPATIAL_POSITION",Aij,D,F);

  SnesCtx snes_ctx(mField,"ELASTIC_MECHANICS");
  
  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

  SnesCtx::loops_to_do_type& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
  snes_ctx.get_preProcess_to_do_Rhs().push_back(&MyDirichletBC);
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Rhs.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_spatial));
  snes_ctx.get_postProcess_to_do_Rhs().push_back(&MyDirichletBC);

  SnesCtx::loops_to_do_type& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
  snes_ctx.get_preProcess_to_do_Mat().push_back(&MyDirichletBC);
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("ELASTIC",&my_fe));
  loops_to_do_Mat.push_back(SnesCtx::loop_pair_type("NEUAMNN_FE",&fe_spatial));
  snes_ctx.get_postProcess_to_do_Mat().push_back(&MyDirichletBC);

  ierr = mField.set_local_VecCreateGhost("ELASTIC_MECHANICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  double step_size = -1.;
  ierr = neumann_forces.setForceScale(step_size); CHKERRQ(ierr);
  ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
  int its;
  ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  {

    //Open mesh_file_name.txt for writing
    ofstream myfile;
    #ifdef NONLINEAR_TEMPERATUTE
    myfile.open(("nonlinear_thermal_expansion_"+string(mesh_file_name)+".txt").c_str());
    #else
    myfile.open(("nonlinear_"+string(mesh_file_name)+".txt").c_str());
    #endif
    
    //Output displacements
    cout << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    myfile << "<<<< Displacements (X-Translation, Y-Translation, Z-Translation) >>>>>" << endl;
    
    for(_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(mField,"SPATIAL_POSITION",dof_ptr))
    {
        if(dof_ptr->get_ent_type()!=MBVERTEX) continue;
  
	double coords[3];
	EntityHandle ent = dof_ptr->get_ent();
	rval = moab.get_coords(&ent,1,coords); CHKERR_PETSC(rval);
        
        if(dof_ptr->get_dof_rank()==0)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==1)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << "  ";
            myfile << boost::format("%.3lf") % roundn(fval) << "  ";
        }
        if(dof_ptr->get_dof_rank()==2)
        {
            //Round and truncate to 3 decimal places
            double fval = dof_ptr->get_FieldData()-coords[dof_ptr->get_dof_rank()];
            cout << boost::format("%.3lf") % roundn(fval) << endl;
            myfile << boost::format("%.3lf") % roundn(fval) << endl;
        }
        
    }

  }

  /*//Save data on mesh
  ierr = mField.set_global_VecCreateGhost("ELASTIC_MECHANICS",Col,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  //ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  PostProcVertexMethod ent_method(moab,"SPATIAL_POSITION");
  ierr = mField.loop_dofs("ELASTIC_MECHANICS","SPATIAL_POSITION",Col,ent_method); CHKERRQ(ierr);

  if(pcomm->rank()==0) {
    EntityHandle out_meshset;
    rval = moab.create_meshset(MESHSET_SET,out_meshset); CHKERR_PETSC(rval);
    ierr = mField.problem_get_FE("ELASTIC_MECHANICS","ELASTIC",out_meshset); CHKERRQ(ierr);
    rval = moab.write_file("out.vtk","VTK","",&out_meshset,1); CHKERR_PETSC(rval);
    rval = moab.delete_entities(&out_meshset,1); CHKERR_PETSC(rval);
  }

  PostProcFieldsAndGradientOnRefMesh fe_post_proc_method(moab);
  ierr = mField.loop_finite_elements("ELASTIC_MECHANICS","ELASTIC",fe_post_proc_method);  CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  if(pcomm->rank()==0) {
    rval = fe_post_proc_method.moab_post_proc.write_file("out_post_proc.vtk","VTK",""); CHKERR_PETSC(rval);
  }*/

  //detroy matrices
  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  ierr = SNESDestroy(&snes); CHKERRQ(ierr);

  PetscFinalize();

  return 0;
}



