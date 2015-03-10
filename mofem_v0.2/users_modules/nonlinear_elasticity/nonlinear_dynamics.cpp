/* \file nonlinear_dynamics.cpp
 *
 * \brief Non-linear elastic dynamics.
 *
 * This file is part of MoFEM.
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
#include <SurfacePressureComplexForLazy.hpp>
#include <adolc/adolc.h> 
#include <ConvectiveMassElement.hpp>
#include <NonLienarElasticElement.hpp>

#include <PotsProcOnRefMesh.hpp>
#include <PostProcStresses.hpp>

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include <ElasticMaterials.hpp>

#define BLOCKED_PROBLEM

using namespace ObosleteUsersModules;

ErrorCode rval;
PetscErrorCode ierr;

static char help[] = "...\n\n";

struct MonitorPostProc: public FEMethod {

  FieldInterface &mField;
  PostPocOnRefinedMesh postProc;
  map<int,NonlinearElasticElement::BlockData> &setOfBlocks; 
  NonlinearElasticElement::MyVolumeFE &feElasticEnergy; ///< calculate elastic energy 
  ConvectiveMassElement::MyVolumeFE &feKineticEnergy; ///< calculate elastic energy 

  bool iNit;

  int pRT;
  int *step;

  MonitorPostProc(FieldInterface &m_field,
    map<int,NonlinearElasticElement::BlockData> &set_of_blocks,
    NonlinearElasticElement::MyVolumeFE &fe_elastic_energy,
    ConvectiveMassElement::MyVolumeFE &fe_kinetic_energy): 
    FEMethod(),mField(m_field),postProc(m_field),
    setOfBlocks(set_of_blocks),
    feElasticEnergy(fe_elastic_energy),
    feKineticEnergy(fe_kinetic_energy),
    iNit(false) { 
    
    ErrorCode rval;
    PetscErrorCode ierr;
    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    }

    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }

  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;
    PetscErrorCode ierr;
    ErrorCode rval;

    if(!iNit) {
      ierr = postProc.generateRefereneElemenMesh(); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_VELOCITY"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);

      map<int,NonlinearElasticElement::BlockData>::iterator sit = setOfBlocks.begin();
      for(;sit!=setOfBlocks.end();sit++) {
	postProc.get_op_to_do_Rhs().push_back(
	  new PostPorcStress(
	    postProc.postProcMesh,
	    postProc.mapGaussPts,
	    "SPATIAL_POSITION",
	    sit->second,
	    postProc.commonData));
      }

      iNit = true;
    }

    if((*step)%pRT==0) {
      ierr = mField.loop_finite_elements("DYNAMICS","MASS_ELEMENT",postProc); CHKERRQ(ierr);
      ostringstream sss;
      sss << "out_values_" << (*step) << ".h5m";
      rval = postProc.postProcMesh.write_file(sss.str().c_str(),"MOAB","PARALLEL=WRITE_PART"); CHKERR_PETSC(rval);
    }

    feElasticEnergy.snes_ctx = SnesMethod::CTX_SNESNONE;
    ierr = mField.loop_finite_elements("DYNAMICS","ELASTIC",feElasticEnergy); CHKERRQ(ierr);
    feKineticEnergy.ts_ctx = TSMethod::CTX_TSNONE;
    ierr = mField.loop_finite_elements("DYNAMICS","MASS_ELEMENT",feKineticEnergy); CHKERRQ(ierr);
    double E = feElasticEnergy.eNergy;
    double T = feKineticEnergy.eNergy;
    PetscPrintf(PETSC_COMM_WORLD,
      "%D Time %3.2e Elastic energy %3.2e Kinetic Energy %3.2e Total %3.2e\n",
      ts_step,ts_t,E,T,E+T);

    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

struct MonitorRestart: public FEMethod {

  double *time;
  int *step;
  FieldInterface &mField;
  int pRT;

  MonitorRestart(FieldInterface &m_field,TS ts): mField(m_field) {

    PetscErrorCode ierr;
    ErrorCode rval;
    double def_t_val = 0;

    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_time;
    rval = m_field.get_moab().tag_get_handle("_TsTime_",1,MB_TYPE_DOUBLE,th_time,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); CHKERR(rval);
      ierr = TSSetTime(ts,*time); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      rval = m_field.get_moab().tag_set_data(th_time,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); CHKERR(rval);
    }
    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val); 
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); CHKERR(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); CHKERR(rval);
    }

    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }


  }

  PetscErrorCode preProcess() {
    PetscFunctionBegin;

    //PetscErrorCode ierr;
    ErrorCode rval;
    (*time) = ts_t;
    if(pRT>0) {
      if((*step)%pRT==0) {
	ostringstream ss;
	ss << "restart_" << (*step) << ".h5m";
	rval = mField.get_moab().write_file(ss.str().c_str()/*,"MOAB","PARALLEL=WRITE_PART"*/); CHKERR_PETSC(rval);
      }
    }
    (*step)++;
    PetscFunctionReturn(0);
  }

  PetscErrorCode operator()() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

  PetscErrorCode postProcess() {
    PetscFunctionBegin;
    PetscFunctionReturn(0);
  }

};

//See file users_modules/elasticity/TimeForceScale.hpp
#include <TimeForceScale.hpp>

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);


  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // use this if your mesh is partotioned and you run code on parts, 
  // you can solve very big problems 
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-my_is_partitioned",&is_partitioned,&flg); CHKERRQ(ierr);
 
  if(is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option = "PARALLEL=BCAST_DELETE;PARALLEL_RESOLVE_SHARED_ENTS;PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
    rval = pcomm->resolve_shared_ents(0,3,0); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,1); CHKERR_PETSC(rval);
    rval = pcomm->resolve_shared_ents(0,3,2); CHKERR_PETSC(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval); 
  }

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  Range CubitSIDESETs_meshsets;
  ierr = m_field.get_Cubit_meshsets(SIDESET,CubitSIDESETs_meshsets); CHKERRQ(ierr);

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0,0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERR_PETSC(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);
  ierr = m_field.get_entities_by_ref_level(bit_level0,BitRefLevel().set(),meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

  bool check_if_spatial_field_exist = m_field.check_field("SPATIAL_POSITION");
  ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_POSITION"); CHKERRQ(ierr);

  //set app. order
  PetscInt disp_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_disp_order",&disp_order,&flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    disp_order = 1;	
  }
  PetscInt vel_order = disp_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-my_vel_order",&vel_order,&flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    vel_order = disp_order;	
  }
  
  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

  ierr = m_field.add_finite_element("NEUMANN_FE",MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("NEUMANN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("NEUMANN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUMANN_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("NEUMANN_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUMANN_FE"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERR_PETSC(rval);
    ierr = m_field.add_ents_to_finite_element_by_TRIs(tris,"NEUMANN_FE"); CHKERRQ(ierr);
  }
  //add nodal force element
  ierr = MetaNodalForces::addNodalForceElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);

  //Velocity
  ierr = m_field.add_field("SPATIAL_VELOCITY",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"SPATIAL_VELOCITY"); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  ierr = m_field.add_field("DOT_SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = m_field.add_field("DOT_SPATIAL_VELOCITY",H1,3,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(0,"DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  //set material model and mass element
  NonlinearElasticElement elastic(m_field,2);
  ElasticMaterials elastic_materials(m_field);
  ierr = elastic_materials.setBlocks(elastic.setOfBlocks); CHKERRQ(ierr);
  //NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<adouble> st_venant_kirchhoff_material_adouble;
  //NonlinearElasticElement::FunctionsToCalulatePiolaKirchhoffI<double> st_venant_kirchhoff_material_double;
  //ierr = elastic.setBlocks(&st_venant_kirchhoff_material_double,&st_venant_kirchhoff_material_adouble); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

  //set mass element
  ConvectiveMassElement inertia(m_field,1);
  //ierr = inertia.setBlocks(); CHKERRQ(ierr);
  ierr = elastic_materials.setBlocks(inertia.setOfBlocks); CHKERRQ(ierr);
  ierr = inertia.addConvectiveMassElement("MASS_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = inertia.addVelocityElement("VELOCITY_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);

  MonitorPostProc post_proc(m_field,elastic.setOfBlocks,elastic.getLoopFeEnergy(),inertia.getLoopFeEnergy());

  #ifdef BLOCKED_PROBLEM
    // elastic and mass element calculated in Kuu shell matrix problem. To
    // calculate Mass element, velocity field is needed.
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","SPATIAL_VELOCITY"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC","DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
  #endif

  //build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //10 node tets
  if(!check_if_spatial_field_exist) {
    Projection10NodeCoordsOnField ent_method_material(m_field,"MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS",ent_method_material); CHKERRQ(ierr);
    Projection10NodeCoordsOnField ent_method_spatial(m_field,"SPATIAL_POSITION");
    ierr = m_field.loop_dofs("SPATIAL_POSITION",ent_method_spatial); CHKERRQ(ierr);
  }

  //build finite elements
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //define problems
  #ifdef BLOCKED_PROBLEM
  {
    ierr = m_field.add_problem("Kuu",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","ELASTIC"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","NEUMANN_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","FORCE_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_add_bit("Kuu",bit_level0); CHKERRQ(ierr);
    if(is_partitioned) {
      ierr = m_field.build_partitioned_problem("Kuu"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("Kuu",true,0,pcomm->size(),1); CHKERRQ(ierr);
    } else {
      ierr = m_field.build_problem("Kuu"); CHKERRQ(ierr);
      ierr = m_field.partition_problem("Kuu"); CHKERRQ(ierr);
      ierr = m_field.partition_finite_elements("Kuu"); CHKERRQ(ierr);
    }
    ierr = m_field.partition_ghost_dofs("Kuu"); CHKERRQ(ierr);
  }
  #endif 

  ierr = m_field.add_problem("DYNAMICS",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","NEUMANN_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","MASS_ELEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","VELOCITY_ELEMENT"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("DYNAMICS",bit_level0); CHKERRQ(ierr);

  if(is_partitioned) {
    ierr = m_field.build_partitioned_problem("DYNAMICS"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("DYNAMICS",true,0,pcomm->size(),1); CHKERRQ(ierr);
  } else {
    ierr = m_field.build_problem("DYNAMICS"); CHKERRQ(ierr);
    ierr = m_field.partition_problem("DYNAMICS"); CHKERRQ(ierr);
    ierr = m_field.partition_finite_elements("DYNAMICS"); CHKERRQ(ierr);
  }
  ierr = m_field.partition_ghost_dofs("DYNAMICS"); CHKERRQ(ierr);

  //create tS
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.VecCreateGhost("DYNAMICS",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

  PetscBool linear;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-is_linear",&linear,&linear); CHKERRQ(ierr);

  #ifdef BLOCKED_PROBLEM
    //shell matrix
    ConvectiveMassElement::MatShellCtx *shellAij_ctx = new ConvectiveMassElement::MatShellCtx();
    ierr = m_field.MatCreateMPIAIJWithArrays("Kuu",&shellAij_ctx->K); CHKERRQ(ierr);
    ierr = MatDuplicate(shellAij_ctx->K,MAT_DO_NOT_COPY_VALUES,&shellAij_ctx->M); CHKERRQ(ierr);
    ierr = shellAij_ctx->iNit(); CHKERRQ(ierr);
    ierr = m_field.VecScatterCreate(D,"DYNAMICS",COL,shellAij_ctx->u,"Kuu",COL,&shellAij_ctx->scatterU); CHKERRQ(ierr);
    ierr = m_field.VecScatterCreate(D,"DYNAMICS","SPATIAL_VELOCITY",COL,shellAij_ctx->v,"Kuu","SPATIAL_POSITION",COL,&shellAij_ctx->scatterV); CHKERRQ(ierr);
    Mat shell_Aij;
    const MoFEMProblem *problem_ptr;
    ierr = m_field.get_problem("DYNAMICS",&problem_ptr); CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD,
      problem_ptr->get_nb_local_dofs_row(),problem_ptr->get_nb_local_dofs_col(),
      problem_ptr->get_nb_dofs_row(),problem_ptr->get_nb_dofs_row(),
      (void*)shellAij_ctx,&shell_Aij); CHKERRQ(ierr);
    ierr = MatShellSetOperation(shell_Aij,MATOP_MULT,(void(*)(void))ConvectiveMassElement::MultOpA); CHKERRQ(ierr);
    ierr = MatShellSetOperation(shell_Aij,MATOP_ZERO_ENTRIES,(void(*)(void))ConvectiveMassElement::ZeroEntriesOp); CHKERRQ(ierr);
    //blocked problem
    ConvectiveMassElement::ShellMatrixElement shell_matrix_element(m_field);
    SpatialPositionsBCFEMethodPreAndPostProc shell_dirihlet_bc(
      m_field,"SPATIAL_POSITION",shellAij_ctx->barK,PETSC_NULL,PETSC_NULL);
    shell_dirihlet_bc.fixFields.push_back("SPATIAL_VELOCITY");
    SpatialPositionsBCFEMethodPreAndPostProc my_dirihlet_bc(
      m_field,"SPATIAL_POSITION",PETSC_NULL,D,F);
    my_dirihlet_bc.fixFields.push_back("SPATIAL_VELOCITY");
    shell_matrix_element.problemName = "Kuu";
    shell_matrix_element.shellMatCtx = shellAij_ctx;
    shell_matrix_element.dirihletBcPtr = &shell_dirihlet_bc;
    shell_matrix_element.loopK.push_back(ConvectiveMassElement::ShellMatrixElement::LoopPairType("ELASTIC",&elastic.getLoopFeLhs()));
    ConvectiveMassElement::ShellResidualElement shell_matrix_residual(m_field);
    shell_matrix_residual.shellMatCtx = shellAij_ctx;
    //surface forces
    NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,shellAij_ctx->barK,F);
    NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_spatial = neumann_forces.getLoopSpatialFe();
    if(linear) {
      fe_spatial.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = fe_spatial.addForce(it->get_msId()); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      ierr = fe_spatial.addPreassure(it->get_msId()); CHKERRQ(ierr);
    }
    fe_spatial.methodsOp.push_back(new TimeForceScale());
    shell_matrix_element.loopK.push_back(ConvectiveMassElement::ShellMatrixElement::LoopPairType("NEUMANN_FE",&fe_spatial));

    ierr = inertia.setShellMatrixMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",linear); CHKERRQ(ierr);
    //element name "ELASTIC" is used, therefore M matrix is assembled as K
    //matrix. This is added to M is shell matrix. M matrix is a derivative of
    //inertia forces over spatial velocities
    shell_matrix_element.loopM.push_back(ConvectiveMassElement::ShellMatrixElement::LoopPairType("ELASTIC",&inertia.getLoopFeMassLhs()));
    //this calculate derivatives of inertia forces over spatial positions and add this to shell K matrix
    shell_matrix_element.loopAuxM.push_back(ConvectiveMassElement::ShellMatrixElement::LoopPairType("ELASTIC",&inertia.getLoopFeMassAuxLhs()));
  #else
    Mat Aij;
    ierr = m_field.MatCreateMPIAIJWithArrays("DYNAMICS",&Aij); CHKERRQ(ierr);
    SpatialPositionsBCFEMethodPreAndPostProc my_dirihlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);
    my_dirihlet_bc.fixFields.push_back("SPATIAL_VELOCITY");
  
    //surface forces
    NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,Aij,F);
    NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &fe_spatial = neumann_forces.getLoopSpatialFe();
    if(linear) {
      fe_spatial.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = fe_spatial.addForce(it->get_msId()); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      ierr = fe_spatial.addPreassure(it->get_msId()); CHKERRQ(ierr);
    }
    fe_spatial.methodsOp.push_back(new TimeForceScale());
  
    ierr = inertia.setConvectiveMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",false,linear); CHKERRQ(ierr);
    ierr = inertia.setVelocityOperators("SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  #endif

  //nodal forces
  boost::ptr_map<string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->get_msId(),true);  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }

  MonitorRestart monitor_restart(m_field,ts);
  ConvectiveMassElement::UpdateAndControl update_and_control(m_field,ts,"SPATIAL_VELOCITY","SPATIAL_POSITION");

  //TS
  TsCtx ts_ctx(m_field,"DYNAMICS");

  //right hand side
  //preprocess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_and_control);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirihlet_bc);
  //fe looops
  TsCtx::loops_to_do_type& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("ELASTIC",&elastic.getLoopFeRhs()));
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("NEUMANN_FE",&fe_spatial));
  boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type(fit->first,&fit->second->getLoopFe()));
  }
  loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&inertia.getLoopFeMassRhs()));
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_preProcess_to_do_IFunction().push_back(&shell_matrix_residual);
  #else
    loops_to_do_Rhs.push_back(TsCtx::loop_pair_type("VELOCITY_ELEMENT",&inertia.getLoopFeVelRhs()));
  #endif
  //postproc
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirihlet_bc);
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_postProcess_to_do_IFunction().push_back(&shell_matrix_residual);
  #endif

  //left hand side 
  //preprocess
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&update_and_control);
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirihlet_bc);
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_preProcess_to_do_IJacobian().push_back(&shell_matrix_element);
  #else 
     //fe loops
    TsCtx::loops_to_do_type& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("ELASTIC",&elastic.getLoopFeLhs()));
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("NEUMANN_FE",&fe_spatial));
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("VELOCITY_ELEMENT",&inertia.getLoopFeVelLhs()));
    loops_to_do_Mat.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&inertia.getLoopFeMassLhs()));
    //postrocess
    ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirihlet_bc);
  #endif
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&update_and_control);
  //monitor
  TsCtx::loops_to_do_type& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&post_proc));
  loops_to_do_Monitor.push_back(TsCtx::loop_pair_type("MASS_ELEMENT",&monitor_restart));

  ierr = TSSetIFunction(ts,F,f_TSSetIFunction,&ts_ctx); CHKERRQ(ierr);
  #ifdef BLOCKED_PROBLEM
    ierr = TSSetIJacobian(ts,shell_Aij,shell_Aij,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  #else 
    ierr = TSSetIJacobian(ts,Aij,Aij,f_TSSetIJacobian,&ts_ctx); CHKERRQ(ierr);
  #endif 
  ierr = TSMonitorSet(ts,f_TSMonitorSet,&ts_ctx,PETSC_NULL); CHKERRQ(ierr);

  double ftime = 1;
  ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
  ierr = TSSetSolution(ts,D); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
  #ifdef BLOCKED_PROBLEM
    //shell matrix pre-conditioner
    SNES snes;
    ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    KSP ksp;
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    PC pc;
    ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);
    ierr = PCSetType(pc,PCSHELL); CHKERRQ(ierr);
    ConvectiveMassElement::PCShellCtx pc_shell_ctx(shell_Aij);
    ierr = PCShellSetContext(pc,(void*)&pc_shell_ctx); CHKERRQ(ierr);
    ierr = PCShellSetApply(pc,ConvectiveMassElement::PCShellApplyOp); CHKERRQ(ierr);
    ierr = PCShellSetSetUp(pc,ConvectiveMassElement::PCShellSetUpOp); CHKERRQ(ierr);
    ierr = PCShellSetDestroy(pc,ConvectiveMassElement::PCShellDestroy);  CHKERRQ(ierr);
  #endif

  ierr = m_field.set_local_VecCreateGhost("DYNAMICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = TSSolve(ts,D); CHKERRQ(ierr);
  ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);

  PetscInt steps,snesfails,rejects,nonlinits,linits;
  ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
  ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
  ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
  ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
  ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
    "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
    steps,rejects,snesfails,ftime,nonlinits,linits);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);


  ierr = VecDestroy(&F); CHKERRQ(ierr);
  ierr = VecDestroy(&D); CHKERRQ(ierr);
  #ifdef BLOCKED_PROBLEM
    ierr = MatDestroy(&shellAij_ctx->K); CHKERRQ(ierr);
    ierr = MatDestroy(&shellAij_ctx->M); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&shellAij_ctx->scatterU); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&shellAij_ctx->scatterV); CHKERRQ(ierr);
    ierr = MatDestroy(&shell_Aij); CHKERRQ(ierr);
    delete shellAij_ctx;
  #else 
    ierr = MatDestroy(&Aij); CHKERRQ(ierr);
  #endif
 
  PetscFinalize();

  return 0;


}



