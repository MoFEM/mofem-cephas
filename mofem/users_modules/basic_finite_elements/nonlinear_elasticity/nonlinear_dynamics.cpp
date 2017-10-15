/** \file nonlinear_gs.cpp
 * \ingroup nonlinear_elastic_elem
 *
 * \brief Non-linear elastic dynamics.

 NOTE: For block solver is only for settings, some features are not implemented
 for this part.

 */

/*
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

#include <BasicFiniteElements.hpp>
using namespace MoFEM;

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include <ElasticMaterials.hpp>
#include <SurfacePressureComplexForLazy.hpp>

#define BLOCKED_PROBLEM




static char help[] = "...\n\n";

struct MonitorPostProc: public FEMethod {

  MoFEM::Interface &mField;
  PostProcVolumeOnRefinedMesh postProc;
  std::map<int, NonlinearElasticElement::BlockData> &setOfBlocks;
  NonlinearElasticElement::MyVolumeFE &feElasticEnergy;   ///< calculate elastic energy
  ConvectiveMassElement::MyVolumeFE &feKineticEnergy;     ///< calculate elastic energy

  bool iNit;

  int pRT;
  int *step;

  MonitorPostProc(
    MoFEM::Interface &m_field,
    std::map<int,NonlinearElasticElement::BlockData> &set_of_blocks,
    NonlinearElasticElement::MyVolumeFE &fe_elastic_energy,
    ConvectiveMassElement::MyVolumeFE &fe_kinetic_energy
  ):
  FEMethod(),mField(m_field),postProc(m_field),
  setOfBlocks(set_of_blocks),
  feElasticEnergy(fe_elastic_energy),
  feKineticEnergy(fe_kinetic_energy),
  iNit(false) {



    double def_t_val = 0;
    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val);
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); MOAB_THROW(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); MOAB_THROW(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); MOAB_THROW(rval);
    }

    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }

  }

  PetscErrorCode preProcess() {
    MoFEMFunctionBeginHot;



    if(!iNit) {
      ierr = postProc.generateReferenceElementMesh(); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("SPATIAL_VELOCITY"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesPostProc("MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = postProc.addFieldValuesGradientPostProc("SPATIAL_POSITION"); CHKERRQ(ierr);

      std::map<int,NonlinearElasticElement::BlockData>::iterator sit = setOfBlocks.begin();
      for (; sit != setOfBlocks.end(); sit++) {
        postProc.getOpPtrVector().push_back(
          new PostProcStress(
            postProc.postProcMesh,
            postProc.mapGaussPts,
            "SPATIAL_POSITION",
            sit->second,
            postProc.commonData)
          );
        }

        iNit = true;
      }

      if((*step)%pRT==0) {
        ierr = mField.loop_finite_elements("DYNAMICS","MASS_ELEMENT",postProc); CHKERRQ(ierr);
        std::ostringstream sss;
        sss << "out_values_" << (*step) << ".h5m";
        ierr = postProc.writeFile(sss.str().c_str()); CHKERRQ(ierr);
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

        MoFEMFunctionReturnHot(0);
      }

  PetscErrorCode operator()() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode postProcess() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

};

struct MonitorRestart: public FEMethod {

  double *time;
  int *step;
  MoFEM::Interface &mField;
  int pRT;

  MonitorRestart(MoFEM::Interface &m_field,TS ts): mField(m_field) {



    double def_t_val = 0;

    const EntityHandle root_meshset = mField.get_moab().get_root_set();

    Tag th_time;
    rval = m_field.get_moab().tag_get_handle("_TsTime_",1,MB_TYPE_DOUBLE,th_time,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val);
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); MOAB_THROW(rval);
      ierr = TSSetTime(ts,*time); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    } else {
      rval = m_field.get_moab().tag_set_data(th_time,&root_meshset,1,&def_t_val); MOAB_THROW(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_time,&root_meshset,1,(const void**)&time); MOAB_THROW(rval);
    }
    Tag th_step;
    rval = m_field.get_moab().tag_get_handle("_TsStep_",1,MB_TYPE_INTEGER,th_step,MB_TAG_CREAT|MB_TAG_EXCL|MB_TAG_MESH,&def_t_val);
    if(rval == MB_ALREADY_ALLOCATED) {
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); MOAB_THROW(rval);
    } else {
      rval = m_field.get_moab().tag_set_data(th_step,&root_meshset,1,&def_t_val); MOAB_THROW(rval);
      rval = m_field.get_moab().tag_get_by_ptr(th_step,&root_meshset,1,(const void**)&step); MOAB_THROW(rval);
    }

    PetscBool flg = PETSC_TRUE;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_output_prt",&pRT,&flg); CHKERRABORT(PETSC_COMM_WORLD,ierr);
    if(flg!=PETSC_TRUE) {
      pRT = 10;
    }


  }

  PetscErrorCode preProcess() {
    MoFEMFunctionBeginHot;

    //

    (*time) = ts_t;
    // if(pRT>0) {
    //   if((*step)%pRT==0) {
    //     std::ostringstream ss;
    //     ss << "restart_" << (*step) << ".h5m";
    //     rval = mField.get_moab().write_file(ss.str().c_str()/*,"MOAB","PARALLEL=WRITE_PART"*/); CHKERRQ_MOAB(rval);
    //   }
    // }
    (*step)++;
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode operator()() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

  PetscErrorCode postProcess() {
    MoFEMFunctionBeginHot;
    MoFEMFunctionReturnHot(0);
  }

};

//See file users_modules/elasticity/TimeForceScale.hpp
#include <TimeForceScale.hpp>

int main(int argc, char *argv[]) {

  PetscInitialize(&argc, &argv, (char *)0, help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if (pcomm == NULL) pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  // use this if your mesh is partitioned and you run code on parts,
  // you can solve very big problems
  PetscBool is_partitioned = PETSC_FALSE;
  ierr = PetscOptionsGetBool(
    PETSC_NULL,PETSC_NULL, "-my_is_partitioned", &is_partitioned, &flg
  ); CHKERRQ(ierr);

  PetscBool linear;
  ierr = PetscOptionsGetBool(PETSC_NULL,PETSC_NULL,"-is_linear",&linear,&linear); CHKERRQ(ierr);

  if (is_partitioned == PETSC_TRUE) {
    //Read mesh to MOAB
    const char *option;
    option =
      "PARALLEL=BCAST_DELETE;"
      "PARALLEL_RESOLVE_SHARED_ENTS;"
      "PARTITION=PARALLEL_PARTITION;";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  } else {
    const char *option;
    option = "";
    rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  }

  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //ref meshset ref level 0
  ierr = m_field.seed_ref_level_3D(0, 0); CHKERRQ(ierr);
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET, meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0, bit_level0); CHKERRQ(ierr);
  ierr = m_field.query_interface<Tools>()->getEntitiesByRefLevel(
    bit_level0, BitRefLevel().set(), meshset_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0, MBTET,"MESH_NODE_POSITIONS", 2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0, MBTRI,"MESH_NODE_POSITIONS", 2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0, MBEDGE,"MESH_NODE_POSITIONS", 2); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0, MBVERTEX,"MESH_NODE_POSITIONS", 1); CHKERRQ(ierr);

  bool check_if_spatial_field_exist = m_field.check_field("SPATIAL_POSITION");
  ierr = m_field.add_field("SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  //add entitities (by tets) to the field
  ierr = m_field.add_ents_to_field_by_type(0,MBTET, "SPATIAL_POSITION"); CHKERRQ(ierr);

  //set app. order
  PetscInt disp_order;
  ierr = PetscOptionsGetInt(
    PETSC_NULL,PETSC_NULL, "-my_disp_order", &disp_order, &flg); CHKERRQ(ierr);
  if(flg!=PETSC_TRUE) {
    disp_order = 1;
  }
  PetscInt vel_order = disp_order;
  ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-my_vel_order",&vel_order,&flg); CHKERRQ(ierr);
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
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET, it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI,"NEUMANN_FE"); CHKERRQ(ierr);
  }
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
    Range tris;
    rval = moab.get_entities_by_type(it->meshset,MBTRI,tris,true); CHKERRQ_MOAB(rval);
    ierr = m_field.add_ents_to_finite_element_by_type(tris,MBTRI, "NEUMANN_FE"); CHKERRQ(ierr);
  }
  // Add nodal force element
  ierr = MetaNodalForces::addElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);
  // Add fluid pressure finite elements
  FluidPressure fluid_pressure_fe(m_field);
  fluid_pressure_fe.addNeumannFluidPressureBCElements("SPATIAL_POSITION");
  fluid_pressure_fe.setNeumannFluidPressureFiniteElementOperators(
    "SPATIAL_POSITION",PETSC_NULL,false,true
  );

  // Velocity
  ierr = m_field.add_field("SPATIAL_VELOCITY",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"SPATIAL_VELOCITY"); CHKERRQ(ierr);

  ierr = m_field.set_field_order(0,MBTET,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  ierr = m_field.add_field("DOT_SPATIAL_POSITION",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_POSITION",disp_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_POSITION",1); CHKERRQ(ierr);
  ierr = m_field.add_field("DOT_SPATIAL_VELOCITY",H1,AINSWORTH_LEGENDRE_BASE,3,MB_TAG_SPARSE,MF_ZERO); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_type(0,MBTET,"DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTET,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBTRI,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBEDGE,"DOT_SPATIAL_VELOCITY",vel_order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(0,MBVERTEX,"DOT_SPATIAL_VELOCITY",1); CHKERRQ(ierr);

  // Set material model and mass element
  NonlinearElasticElement elastic(m_field, 2);
  ElasticMaterials elastic_materials(m_field);
  ierr = elastic_materials.setBlocks(elastic.setOfBlocks); CHKERRQ(ierr);
  //NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<adouble> st_venant_kirchhoff_material_adouble;
  //NonlinearElasticElement::FunctionsToCalculatePiolaKirchhoffI<double> st_venant_kirchhoff_material_double;
  //ierr = elastic.setBlocks(&st_venant_kirchhoff_material_double,&st_venant_kirchhoff_material_adouble); CHKERRQ(ierr);
  ierr = elastic.addElement("ELASTIC","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = elastic.setOperators("SPATIAL_POSITION"); CHKERRQ(ierr);

  //set mass element
  ConvectiveMassElement inertia(m_field,1);
  //ierr = inertia.setBlocks(); CHKERRQ(ierr);
  ierr = elastic_materials.setBlocks(inertia.setOfBlocks); CHKERRQ(ierr);
  ierr = inertia.addConvectiveMassElement("MASS_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  ierr = inertia.addVelocityElement("VELOCITY_ELEMENT","SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);

  //Add possibility to load accelerogram
  {
    string name = "-my_accelerogram";
    char time_file_name[255];
    PetscBool flg;
    ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,name.c_str(),time_file_name,255,&flg); CHKERRQ(ierr);
    if(flg == PETSC_TRUE) {
      inertia.methodsOp.push_back(new TimeAccelerogram(name));
    }
  }

  //damper element
  KelvinVoigtDamper damper(m_field);
  ierr = elastic_materials.setBlocks(damper.blockMaterialDataMap); CHKERRQ(ierr);
  {
    KelvinVoigtDamper::CommonData &common_data = damper.commonData;
    common_data.spatialPositionName = "SPATIAL_POSITION";
    common_data.spatialPositionNameDot = "DOT_SPATIAL_POSITION";
    ierr = m_field.add_finite_element("DAMPER",MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_row("DAMPER","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_col("DAMPER","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("DAMPER","SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("DAMPER","DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
    if(m_field.check_field("MESH_NODE_POSITIONS")) {
      ierr = m_field.modify_finite_element_add_field_data("DAMPER","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
    }
    std::map<int,KelvinVoigtDamper::BlockMaterialData>::iterator bit = damper.blockMaterialDataMap.begin();
    for(;bit!=damper.blockMaterialDataMap.end();bit++) {
      bit->second.lInear = linear;
      int id = bit->first;
      KelvinVoigtDamper::BlockMaterialData &material_data = bit->second;
      damper.constitutiveEquationMap.insert(
        id,new KelvinVoigtDamper::ConstitutiveEquation<adouble>(material_data)
      );
      ierr = m_field.add_ents_to_finite_element_by_type(bit->second.tEts,MBTET,"DAMPER"); CHKERRQ(ierr);
    }
    ierr = damper.setOperators(3); CHKERRQ(ierr);
  }

  MonitorPostProc post_proc(
    m_field,
    elastic.setOfBlocks,
    elastic.getLoopFeEnergy(),
    inertia.getLoopFeEnergy()
  );

  #ifdef BLOCKED_PROBLEM
    // elastic and mass element calculated in Kuu shell matrix problem. To
    // calculate Mass element, velocity field is needed.
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC", "SPATIAL_VELOCITY"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC", "DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.modify_finite_element_add_field_data("ELASTIC", "DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);
  #endif

  // build field
  ierr = m_field.build_fields(); CHKERRQ(ierr);
  //ierr = m_field.list_dofs_by_field_name("SPATIAL_POSITION"); CHKERRQ(ierr);

  // 10 node tets
  if (!check_if_spatial_field_exist) {
    Projection10NodeCoordsOnField ent_method_material(m_field, "MESH_NODE_POSITIONS");
    ierr = m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method_material); CHKERRQ(ierr);
    Projection10NodeCoordsOnField ent_method_spatial(m_field, "SPATIAL_POSITION");
    ierr = m_field.loop_dofs("SPATIAL_POSITION", ent_method_spatial); CHKERRQ(ierr);
  }

  //build finite elements
  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //define problems
  #ifdef BLOCKED_PROBLEM
  {
    ierr = m_field.add_problem("Kuu", MF_ZERO); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","ELASTIC"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","NEUMANN_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","FORCE_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_add_finite_element("Kuu","FLUID_PRESSURE_FE"); CHKERRQ(ierr);
    ierr = m_field.modify_problem_ref_level_add_bit("Kuu",bit_level0); CHKERRQ(ierr);

    ProblemsManager *prb_mng_ptr;
    ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
    if(is_partitioned) {
      ierr = prb_mng_ptr->buildProblemOnDistributedMesh("Kuu",true); CHKERRQ(ierr);
      ierr = prb_mng_ptr->partitionFiniteElements("Kuu",true,0,pcomm->size()); CHKERRQ(ierr);
    } else {
      ierr = prb_mng_ptr->buildProblem("Kuu",true); CHKERRQ(ierr);
      ierr = prb_mng_ptr->partitionProblem("Kuu"); CHKERRQ(ierr);
      ierr = prb_mng_ptr->partitionFiniteElements("Kuu"); CHKERRQ(ierr);
    }
    ierr = prb_mng_ptr->partitionGhostDofs("Kuu"); CHKERRQ(ierr);
    //ierr = m_field.partition_check_matrix_fill_in("Kuu",-1,-1,0); CHKERRQ(ierr);
  }
  #endif

  ierr = m_field.add_problem("DYNAMICS",MF_ZERO); CHKERRQ(ierr);
  //set finite elements for problems
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","ELASTIC"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","DAMPER"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","NEUMANN_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","FORCE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","FLUID_PRESSURE_FE"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","MASS_ELEMENT"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("DYNAMICS","VELOCITY_ELEMENT"); CHKERRQ(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("DYNAMICS",bit_level0); CHKERRQ(ierr);

  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  if(is_partitioned) {
    ierr = prb_mng_ptr->buildProblemOnDistributedMesh("DYNAMICS",true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("DYNAMICS",true,0,pcomm->size()); CHKERRQ(ierr);
  } else {
    ierr = prb_mng_ptr->buildProblem("DYNAMICS",true); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionProblem("DYNAMICS"); CHKERRQ(ierr);
    ierr = prb_mng_ptr->partitionFiniteElements("DYNAMICS"); CHKERRQ(ierr);
  }
  ierr = prb_mng_ptr->partitionGhostDofs("DYNAMICS"); CHKERRQ(ierr);

  Vec F;
  ierr = m_field.query_interface<VecManager>()->vecCreateGhost("DYNAMICS",COL,&F); CHKERRQ(ierr);
  Vec D;
  ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

  //create tS
  TS ts;
  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);

  #ifdef BLOCKED_PROBLEM
    //shell matrix
    ConvectiveMassElement::MatShellCtx *shellAij_ctx = new ConvectiveMassElement::MatShellCtx();
    ierr = m_field.MatCreateMPIAIJWithArrays("Kuu",&shellAij_ctx->K); CHKERRQ(ierr);
    ierr = MatDuplicate(shellAij_ctx->K,MAT_DO_NOT_COPY_VALUES,&shellAij_ctx->M); CHKERRQ(ierr);
    ierr = shellAij_ctx->iNit(); CHKERRQ(ierr);
    ierr = m_field.query_interface<VecManager>()->vecScatterCreate(D,"DYNAMICS",COL,shellAij_ctx->u,"Kuu",COL,&shellAij_ctx->scatterU); CHKERRQ(ierr);
    ierr = m_field.query_interface<VecManager>()->vecScatterCreate(
      D,"DYNAMICS","SPATIAL_VELOCITY",COL,shellAij_ctx->v,"Kuu","SPATIAL_POSITION",COL,&shellAij_ctx->scatterV
    ); CHKERRQ(ierr);
    Mat shell_Aij;
    const Problem *problem_ptr;
    ierr = m_field.get_problem("DYNAMICS",&problem_ptr); CHKERRQ(ierr);
    ierr = MatCreateShell(
      PETSC_COMM_WORLD,
      problem_ptr->getNbLocalDofsRow(),
      problem_ptr->getNbLocalDofsCol(),
      problem_ptr->getNbDofsRow(),
      problem_ptr->getNbDofsRow(),
      (void*)shellAij_ctx,&shell_Aij
    ); CHKERRQ(ierr);
    ierr = MatShellSetOperation(shell_Aij,MATOP_MULT,(void(*)(void))ConvectiveMassElement::MultOpA); CHKERRQ(ierr);
    ierr = MatShellSetOperation(shell_Aij,MATOP_ZERO_ENTRIES,(void(*)(void))ConvectiveMassElement::ZeroEntriesOp); CHKERRQ(ierr);
    //blocked problem
    ConvectiveMassElement::ShellMatrixElement shell_matrix_element(m_field);
    DirichletSpatialPositionsBc shell_dirichlet_bc(
      m_field,"SPATIAL_POSITION",shellAij_ctx->barK,PETSC_NULL,PETSC_NULL
    );
    DirichletSpatialPositionsBc my_dirichlet_bc(
      m_field,"SPATIAL_POSITION",PETSC_NULL,D,F
    );
    shell_matrix_element.problemName = "Kuu";
    shell_matrix_element.shellMatCtx = shellAij_ctx;
    shell_matrix_element.DirichletBcPtr = &shell_dirichlet_bc;
    shell_matrix_element.loopK.push_back(ConvectiveMassElement::ShellMatrixElement::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeLhs()));
    //damper
    shell_matrix_element.loopK.push_back(ConvectiveMassElement::ShellMatrixElement::PairNameFEMethodPtr("ELASTIC",&damper.feLhs));

    //surface forces
    NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,shellAij_ctx->barK,F);
    NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &surface_force = neumann_forces.getLoopSpatialFe();
    if(linear) {
      surface_force.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = surface_force.addForce(it->getMeshsetId()); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      ierr = surface_force.addPreassure(it->getMeshsetId()); CHKERRQ(ierr);
    }
    surface_force.methodsOp.push_back(new TimeForceScale());
    shell_matrix_element.loopK.push_back(ConvectiveMassElement::ShellMatrixElement::PairNameFEMethodPtr("NEUMANN_FE",&surface_force));

    ierr = inertia.setShellMatrixMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",linear); CHKERRQ(ierr);
    //element name "ELASTIC" is used, therefore M matrix is assembled as K
    //matrix. This is added to M is shell matrix. M matrix is a derivative of
    //inertia forces over spatial velocities
    shell_matrix_element.loopM.push_back(ConvectiveMassElement::ShellMatrixElement::PairNameFEMethodPtr("ELASTIC",&inertia.getLoopFeMassLhs()));
    //this calculate derivatives of inertia forces over spatial positions and add this to shell K matrix
    shell_matrix_element.loopAuxM.push_back(ConvectiveMassElement::ShellMatrixElement::PairNameFEMethodPtr("ELASTIC",&inertia.getLoopFeMassAuxLhs()));

    //Element to calualte shell matrix residual
    ConvectiveMassElement::ShellResidualElement shell_matrix_residual(m_field);
    shell_matrix_residual.shellMatCtx = shellAij_ctx;

  #else
    Mat Aij;
    ierr = m_field.MatCreateMPIAIJWithArrays("DYNAMICS",&Aij); CHKERRQ(ierr);
    DirichletSpatialPositionsBc my_dirichlet_bc(m_field,"SPATIAL_POSITION",Aij,D,F);
    //my_dirichlet_bc.fixFields.push_back("SPATIAL_VELOCITY");

    //surface forces
    NeummanForcesSurfaceComplexForLazy neumann_forces(m_field,Aij,F);
    NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE &surface_force = neumann_forces.getLoopSpatialFe();
    if(linear) {
      surface_force.typeOfForces = NeummanForcesSurfaceComplexForLazy::MyTriangleSpatialFE::NONCONSERVATIVE;
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
      ierr = surface_force.addForce(it->getMeshsetId()); CHKERRQ(ierr);
    }
    for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,SIDESET|PRESSURESET,it)) {
      ierr = surface_force.addPreassure(it->getMeshsetId()); CHKERRQ(ierr);
    }
    surface_force.methodsOp.push_back(new TimeForceScale());

    ierr = inertia.setConvectiveMassOperators("SPATIAL_VELOCITY","SPATIAL_POSITION","MESH_NODE_POSITIONS",false,linear); CHKERRQ(ierr);
    ierr = inertia.setVelocityOperators("SPATIAL_VELOCITY","SPATIAL_POSITION"); CHKERRQ(ierr);
  #endif

  //nodal forces
  boost::ptr_map<std::string,NodalForce> nodal_forces;
  string fe_name_str ="FORCE_FE";
  nodal_forces.insert(fe_name_str,new NodalForce(m_field));
  for(_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(m_field,NODESET|FORCESET,it)) {
    ierr = nodal_forces.at(fe_name_str).addForce("SPATIAL_POSITION",F,it->getMeshsetId(),true);  CHKERRQ(ierr);
    nodal_forces.at(fe_name_str).methodsOp.push_back(new TimeForceScale());
  }

  MonitorRestart monitor_restart(m_field,ts);
  ConvectiveMassElement::UpdateAndControl update_and_control(m_field,ts,"SPATIAL_VELOCITY","SPATIAL_POSITION");

  //TS
  TsCtx ts_ctx(m_field,"DYNAMICS");

  //right hand side
  //preprocess
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&update_and_control);
  ts_ctx.get_preProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  //fe looops
  TsCtx::FEMethodsSequence& loops_to_do_Rhs = ts_ctx.get_loops_to_do_IFunction();
  loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeRhs()));
  loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("DAMPER",&damper.feRhs));
  loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("NEUMANN_FE",&surface_force));
  boost::ptr_map<std::string,NodalForce>::iterator fit = nodal_forces.begin();
  for(;fit!=nodal_forces.end();fit++) {
    loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr(fit->first,&fit->second->getLoopFe()));
  }
  loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("FLUID_PRESSURE_FE",&fluid_pressure_fe.getLoopFe()));
  loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("MASS_ELEMENT",&inertia.getLoopFeMassRhs()));
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_preProcess_to_do_IFunction().push_back(&shell_matrix_residual);
  #else
    loops_to_do_Rhs.push_back(TsCtx::PairNameFEMethodPtr("VELOCITY_ELEMENT",&inertia.getLoopFeVelRhs()));
  #endif
  //postproc
  ts_ctx.get_postProcess_to_do_IFunction().push_back(&my_dirichlet_bc);
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_postProcess_to_do_IFunction().push_back(&shell_matrix_residual);
  #endif

  //left hand side
  //preprocess
  ts_ctx.get_preProcess_to_do_IJacobian().push_back(&update_and_control);
  #ifdef BLOCKED_PROBLEM
    ts_ctx.get_preProcess_to_do_IJacobian().push_back(&shell_matrix_element);
  #else
    //preprocess
    ts_ctx.get_preProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
    //fe loops
    TsCtx::FEMethodsSequence& loops_to_do_Mat = ts_ctx.get_loops_to_do_IJacobian();
    loops_to_do_Mat.push_back(TsCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeLhs()));
    loops_to_do_Mat.push_back(TsCtx::PairNameFEMethodPtr("DAMPER",&damper.feLhs));
    loops_to_do_Mat.push_back(TsCtx::PairNameFEMethodPtr("NEUMANN_FE",&surface_force));
    loops_to_do_Mat.push_back(TsCtx::PairNameFEMethodPtr("VELOCITY_ELEMENT",&inertia.getLoopFeVelLhs()));
    loops_to_do_Mat.push_back(TsCtx::PairNameFEMethodPtr("MASS_ELEMENT",&inertia.getLoopFeMassLhs()));
    //postrocess
    ts_ctx.get_postProcess_to_do_IJacobian().push_back(&my_dirichlet_bc);
  #endif
  ts_ctx.get_postProcess_to_do_IJacobian().push_back(&update_and_control);
  //monitor
  TsCtx::FEMethodsSequence& loops_to_do_Monitor = ts_ctx.get_loops_to_do_Monitor();
  loops_to_do_Monitor.push_back(TsCtx::PairNameFEMethodPtr("MASS_ELEMENT",&post_proc));
  loops_to_do_Monitor.push_back(TsCtx::PairNameFEMethodPtr("MASS_ELEMENT",&monitor_restart));

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
    // ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
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

  ierr = m_field.query_interface<VecManager>()->setLocalGhostVector("DYNAMICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateBegin(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  // Solve problem at time Zero
  PetscBool is_solve_at_time_zero = PETSC_FALSE;
  ierr = PetscOptionsGetBool(
    PETSC_NULL,PETSC_NULL, "-my_solve_at_time_zero",&is_solve_at_time_zero,&flg
  ); CHKERRQ(ierr);
  if(is_solve_at_time_zero) {

    #ifdef BLOCKED_PROBLEM

    Mat Aij = shellAij_ctx->K;
    Vec F;
    ierr = m_field.query_interface<VecManager>()->vecCreateGhost("Kuu",COL,&F); CHKERRQ(ierr);
    Vec D;
    ierr = VecDuplicate(F,&D); CHKERRQ(ierr);

    SnesCtx snes_ctx(m_field,"Kuu");

    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetApplicationContext(snes,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetFunction(snes,F,SnesRhs,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes,Aij,Aij,SnesMat,&snes_ctx); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    DirichletSpatialPositionsBc my_dirichlet_bc(
      m_field,"SPATIAL_POSITION",PETSC_NULL,D,F
    );

    SnesCtx::FEMethodsSequence& loops_to_do_Rhs = snes_ctx.get_loops_to_do_Rhs();
    snes_ctx.get_preProcess_to_do_Rhs().push_back(&my_dirichlet_bc);
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeRhs()));
    fluid_pressure_fe.getLoopFe().ts_t = 0;
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("FLUID_PRESSURE_FE",&fluid_pressure_fe.getLoopFe()));
    surface_force.ts_t = 0;
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("NEUMANN_FE",&surface_force));
    boost::ptr_map<std::string,NodalForce>::iterator fit = nodal_forces.begin();
    for(;fit!=nodal_forces.end();fit++) {
      fit->second->getLoopFe().ts_t = 0;
      loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr(fit->first,&fit->second->getLoopFe()));
    }
    inertia.getLoopFeMassRhs().ts_t = 0;
    loops_to_do_Rhs.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&inertia.getLoopFeMassRhs()));
    snes_ctx.get_postProcess_to_do_Rhs().push_back(&my_dirichlet_bc);

    SnesCtx::FEMethodsSequence& loops_to_do_Mat = snes_ctx.get_loops_to_do_Mat();
    snes_ctx.get_preProcess_to_do_Mat().push_back(&my_dirichlet_bc);
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&elastic.getLoopFeLhs()));
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("NEUMANN_FE",&surface_force));
    inertia.getLoopFeMassAuxLhs().ts_t = 0;
    inertia.getLoopFeMassAuxLhs().ts_a = 0;
    loops_to_do_Mat.push_back(SnesCtx::PairNameFEMethodPtr("ELASTIC",&inertia.getLoopFeMassAuxLhs()));
    snes_ctx.get_postProcess_to_do_Mat().push_back(&my_dirichlet_bc);

    ierr = m_field.field_scale(0,"SPATIAL_VELOCITY"); CHKERRQ(ierr);
    ierr = m_field.field_scale(0,"DOT_SPATIAL_POSITION"); CHKERRQ(ierr);
    ierr = m_field.field_scale(0,"DOT_SPATIAL_VELOCITY"); CHKERRQ(ierr);

    ierr = m_field.query_interface<VecManager>()->setLocalGhostVector("Kuu",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

    ierr = SNESSolve(snes,PETSC_NULL,D); CHKERRQ(ierr);
    int its;
    ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"number of Newton iterations = %D\n",its); CHKERRQ(ierr);

    ierr = m_field.query_interface<VecManager>()->setGlobalGhostVector("Kuu",COL,D,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);

    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = VecDestroy(&D); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    #endif //BLOCKED_PROBLEM
  }

  if(is_solve_at_time_zero) {
    ierr = m_field.query_interface<VecManager>()->setLocalGhostVector("DYNAMICS",COL,D,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  }

  #if PETSC_VERSION_GE(3,7,0)
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER); CHKERRQ(ierr);
  #endif
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


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  PetscFinalize();

  return 0;


}
