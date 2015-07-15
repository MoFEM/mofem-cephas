/** \file gel_analysis.cpp
  \brief Reads cubit file and solves problem with gel material
  \ingroup gel

  1) TODO: Current version is limited only to one material. If one like to have
  general problem for nonlinear elasticity should implement general time
  dependent problem. If inertia terms need to be considered, this material
  should be add to nonlinear dynamics problem.

  2) TODO: Internal history state variables need to be statically condensed. It
  can  be done  by implementing static condensation on finite element level or
  by implementing pre-conditioner.

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

#include <PostProcOnRefMesh.hpp>
#include <Projection10NodeCoordsOnField.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adolc/adolc.h>
#include <Gels.hpp>
#include <UserGelModel.hpp>

#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

#include <DirichletBC.hpp>
#include <SurfacePressure.hpp>
#include <NodalForce.hpp>
#include <EdgeForce.hpp>

// Elements for applying fluxes, convection or radiation
// are used to apply solvent concentration.
#include <ThermalElement.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>
namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

ErrorCode rval;
PetscErrorCode ierr;
static char help[] = "...\n\n";

struct BlockData {
  int oRder;
  BlockData():
  oRder(-1) {
  }
};

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;

  PetscBool flg_gel_config,flg_file;
  char mesh_file_name[255];
  char gel_config_file[255];
  PetscInt order = 2;
  PetscBool is_partitioned = PETSC_FALSE;

  {

    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","Elastic Config","none"); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-my_file",
      "mesh file name","",
      "mesh.h5m",mesh_file_name,255,&flg_file
    ); CHKERRQ(ierr);
    ierr = PetscOptionsInt(
      "-my_order",
      "default approximation order","",
      2,&order,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsBool(
      "-my_is_partitioned",
      "set if mesh is partitioned (this result that each process keeps only part of the mes","",
      PETSC_FALSE,&is_partitioned,PETSC_NULL
    ); CHKERRQ(ierr);
    ierr = PetscOptionsString(
      "-my_gel_config",
      "gel configuration file name","",
      "gel_config.in",gel_config_file,255,&flg_gel_config
    ); CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);

    //Reade parameters from line command
    if(flg_file != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
    }

    if(is_partitioned == PETSC_TRUE) {
      //Read mesh to MOAB
      const char *option;
      option =
      "PARALLEL=BCAST_DELETE;"
      "PARALLEL_RESOLVE_SHARED_ENTS;"
      "PARTITION=PARALLEL_PARTITION;";
      rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
    } else {
      const char *option;
      option = "";
      rval = moab.load_file(mesh_file_name, 0, option); CHKERR_PETSC(rval);
    }
    ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
    if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);
  }

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;

  // Seed all mesh entities to MoFEM database, those entities can be potentially used as finite elements
  // or as entities which carry some approximation field.
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  // Define fields and finite elements
  Gel gel(m_field);
  map<int,ThermalElement::FluxData> set_of_solvent_fluxes;
  {

    // Set approximation fields
    {

      // Add fields
      bool check_if_spatial_field_exist = m_field.check_field("SPATIAL_POSITION");
      ierr = m_field.add_field("SPATIAL_POSITION",H1,3,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_field("SOLVENT_CONCENTRATION",H1,1,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_field("HAT_EPS",L2,6,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_field("SPATIAL_POSITION_DOT",H1,3,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_field("SOLVENT_CONCENTRATION_DOT",H1,1,MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.add_field("HAT_EPS_DOT",L2,6,MF_ZERO); CHKERRQ(ierr);

      //Add field H1 space rank 3 to approximate geometry using hierarchical basis
      //For 10 node tetrahedral, before use, geometry is projected on that field (see below)
      ierr = m_field.add_field("MESH_NODE_POSITIONS",H1,3,MF_ZERO); CHKERRQ(ierr);

      //meshset consisting all entities in mesh
      EntityHandle root_set = moab.get_root_set();
      //Add entities to field (root_mesh, i.e. on all mesh entities fields are approx.)
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SPATIAL_POSITION_DOT"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SOLVENT_CONCENTRATION_DOT"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"HAT_EPS_DOT"); CHKERRQ(ierr);

      // Set approximation order. Solvent concentration has one order less than
      // order of of spatial position field. Tests need to be maid if that is
      // good enough to have stability for this type of problem. If not bubble
      // functions could be easily added by increasing approximate order for
      // volume.

      ierr = m_field.set_field_order(root_set,MBTET,"SPATIAL_POSITION",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTRI,"SPATIAL_POSITION",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBEDGE,"SPATIAL_POSITION",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBVERTEX,"SPATIAL_POSITION",1); CHKERRQ(ierr);

      ierr = m_field.set_field_order(root_set,MBTET,"SPATIAL_POSITION_DOT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTRI,"SPATIAL_POSITION_DOT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBEDGE,"SPATIAL_POSITION_DOT",order); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBVERTEX,"SPATIAL_POSITION_DOT",1); CHKERRQ(ierr);

      ierr = m_field.set_field_order(root_set,MBTET,"SOLVENT_CONCENTRATION",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTRI,"SOLVENT_CONCENTRATION",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBEDGE,"SOLVENT_CONCENTRATION",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBVERTEX,"SOLVENT_CONCENTRATION",1); CHKERRQ(ierr);

      ierr = m_field.set_field_order(root_set,MBTET,"SOLVENT_CONCENTRATION_DOT",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTRI,"SOLVENT_CONCENTRATION_DOT",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBEDGE,"SOLVENT_CONCENTRATION_DOT",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBVERTEX,"SOLVENT_CONCENTRATION_DOT",1); CHKERRQ(ierr);

      ierr = m_field.set_field_order(root_set,MBTET,"HAT_EPS",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTET,"HAT_EPS_DOT",order-1); CHKERRQ(ierr);

      //gemetry approximation is set to 2nd oreder
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

      try {

        map<int,BlockData> block_data;
        Gel::BlockMaterialData &material_data = gel.blockMaterialData;

        ifstream ini_file(gel_config_file);
        po::variables_map vm;
        po::options_description config_file_options;

        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
          ostringstream str_order;
          str_order << "block_" << it->get_msId() << ".oRder";
          config_file_options.add_options()
          (str_order.str().c_str(),po::value<int>(&block_data[it->get_msId()].oRder)->default_value(order));
        }

        config_file_options.add_options()
        ("gAlpha",po::value<double>(&material_data.gAlpha)->default_value(1));
        config_file_options.add_options()
        ("vAlpha",po::value<double>(&material_data.vAlpha)->default_value(0));
        config_file_options.add_options()
        ("gBeta",po::value<double>(&material_data.gBeta)->default_value(1));
        config_file_options.add_options()
        ("vBeta",po::value<double>(&material_data.vBeta)->default_value(0));
        config_file_options.add_options()
        ("gBetaHat",po::value<double>(&material_data.gBetaHat)->default_value(1));
        config_file_options.add_options()
        ("vBetaHat",po::value<double>(&material_data.vBetaHat)->default_value(0));
        config_file_options.add_options()
        ("oMega",po::value<double>(&material_data.oMega)->default_value(1));
        config_file_options.add_options()
        ("vIscosity",po::value<double>(&material_data.vIscosity)->default_value(1));
        config_file_options.add_options()
        ("pErmeability",po::value<double>(&material_data.pErmeability)->default_value(1));

        po::parsed_options parsed = parse_config_file(ini_file,config_file_options,true);
        store(parsed,vm);
        po::notify(vm);
        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
          if(block_data[it->get_msId()].oRder == -1) continue;
          if(block_data[it->get_msId()].oRder == order) continue;
          PetscPrintf(PETSC_COMM_WORLD,"Set block %d order to %d\n",it->get_msId(),block_data[it->get_msId()].oRder);
          Range block_ents;
          rval = moab.get_entities_by_handle(it->get_meshset(),block_ents,true); CHKERR_PETSC(rval);
          Range ents_to_set_order;
          rval = moab.get_adjacencies(block_ents,3,false,ents_to_set_order,Interface::UNION); CHKERR_PETSC(rval);
          ents_to_set_order = ents_to_set_order.subset_by_type(MBTET);
          rval = moab.get_adjacencies(block_ents,2,false,ents_to_set_order,Interface::UNION); CHKERR_PETSC(rval);
          rval = moab.get_adjacencies(block_ents,1,false,ents_to_set_order,Interface::UNION); CHKERR_PETSC(rval);
          ierr = m_field.synchronise_entities(ents_to_set_order); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"SPATIAL_POSITION",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"SOLVENT_CONCENTRATION",block_data[it->get_msId()].oRder-1); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"HAT_EPS",block_data[it->get_msId()].oRder-1); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"SPATIAL_POSITION_DOT",block_data[it->get_msId()].oRder); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"SOLVENT_CONCENTRATION_DOT",block_data[it->get_msId()].oRder-1); CHKERRQ(ierr);
          ierr = m_field.set_field_order(ents_to_set_order,"HAT_EPS_DOT",block_data[it->get_msId()].oRder-1); CHKERRQ(ierr);
        }
        vector<string> additional_parameters;
        additional_parameters = collect_unrecognized(parsed.options,po::include_positional);
        for(vector<string>::iterator vit = additional_parameters.begin();
        vit!=additional_parameters.end();vit++) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"** WARNING Unrecognised option %s\n",vit->c_str()); CHKERRQ(ierr);
        }
      } catch (const std::exception& ex) {
        ostringstream ss;
        ss << ex.what() << endl;
        SETERRQ(PETSC_COMM_SELF,MOFEM_STD_EXCEPTION_THROW,ss.str().c_str());
      }

      ierr = m_field.build_fields(); CHKERRQ(ierr);

      // Sett geometry approximation and initial spatial positions
      // 10 node tets
      if (!check_if_spatial_field_exist) {
        Projection10NodeCoordsOnField ent_method_material(m_field, "MESH_NODE_POSITIONS");
        ierr = m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method_material); CHKERRQ(ierr);
        Projection10NodeCoordsOnField ent_method_spatial(m_field, "SPATIAL_POSITION");
        ierr = m_field.loop_dofs("SPATIAL_POSITION", ent_method_spatial); CHKERRQ(ierr);
      }

    }

    //Set finite elements. The primary element is GEL_FE in addition elements
    //for applying tractions and fluxes are added.
    {
      ierr = m_field.add_finite_element("GEL_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SPATIAL_POSITION_DOT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SOLVENT_CONCENTRATION_DOT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","HAT_EPS_DOT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      EntityHandle root_set = moab.get_root_set();
      ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"GEL_FE"); CHKERRQ(ierr);

      // Add Neumann forces, i.e. on triangles, edges and nodes.
      ierr = MetaNeummanForces::addNeumannBCElements(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = MetaNodalForces::addElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = MetaEdgeForces::addElement(m_field,"SPATIAL_POSITION"); CHKERRQ(ierr);

      // Add solvent flux element
      {

        ierr = m_field.add_finite_element("SOLVENT_FLUX_FE",MF_ZERO); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_row("SOLVENT_FLUX_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_col("SOLVENT_FLUX_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data("SOLVENT_FLUX_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
        ierr = m_field.modify_finite_element_add_field_data("SOLVENT_FLUX_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);

        // Assume that boundary conditions are set in block containing surface
        // triangle elements and block name is "SOLVENT_FLUX"
        for(_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field,BLOCKSET,it)) {
          if(it->get_name().compare(0,12,"SOLVENT_FLUX") == 0) {
            vector<double> data;
            ierr = it->get_attributes(data); CHKERRQ(ierr);
            if(data.size()!=1) {
              SETERRQ(PETSC_COMM_SELF,1,"Data inconsistency");
            }
            // Here it set how block of for heat flux is set.  This is because
            // implementation from thermal element is used to enforce this
            // boundary condition.
            strcpy(set_of_solvent_fluxes[it->get_msId()].dAta.data.name,"HeatFlu");
            set_of_solvent_fluxes[it->get_msId()].dAta.data.flag1 = 1;
            set_of_solvent_fluxes[it->get_msId()].dAta.data.value1 = data[0];
            //cerr << set_of_solvent_fluxes[it->get_msId()].dAta << endl;
            rval = m_field.get_moab().get_entities_by_type(
              it->meshset,MBTRI,set_of_solvent_fluxes[it->get_msId()].tRis,true
            ); CHKERR_PETSC(rval);
            ierr = m_field.add_ents_to_finite_element_by_TRIs(
              set_of_solvent_fluxes[it->get_msId()].tRis,"SOLVENT_FLUX_FE"
            ); CHKERRQ(ierr);
          }
        }
      }

      //build finite elements
      ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
      //build adjacencies
      ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
    }

  }

  // Create gel instance and set operators.
  {

    gel.constitutiveEquationPtr = boost::shared_ptr<UserGelConstitutiveEquation<adouble> >(
      new UserGelConstitutiveEquation<adouble>(gel.blockMaterialData)
    );

    // Set name of fields which has been choose to approximate spatial
    // displacements, solvent concentration and internal state variables.
    Gel::CommonData &common_data = gel.commonData;
    common_data.spatialPositionName = "SPATIAL_POSITION";
    common_data.spatialPositionNameDot = "SPATIAL_POSITION_DOT";
    common_data.muName = "SOLVENT_CONCENTRATION";
    common_data.muNameDot = "SOLVENT_CONCENTRATION_DOT";
    common_data.strainHatName = "HAT_EPS";
    common_data.strainHatNameDot = "HAT_EPS_DOT";

    // Set operators to calculate field values at integration points, both for
    // left and right hand side elements.
    Gel::GelFE *fe_ptr[] = { &gel.feRhs, &gel.feLhs };
    for(int ss = 0;ss<2;ss++) {
      fe_ptr[ss]->getOpPtrVector().push_back(
        new Gel::OpGetDataAtGaussPts("SPATIAL_POSITION",common_data,false,true)
      );
      fe_ptr[ss]->getOpPtrVector().push_back(
        new Gel::OpGetDataAtGaussPts("SPATIAL_POSITION_DOT",common_data,false,true)
      );
      fe_ptr[ss]->getOpPtrVector().push_back(
        new Gel::OpGetDataAtGaussPts("SOLVENT_CONCENTRATION",common_data,true,true)
      );
      fe_ptr[ss]->getOpPtrVector().push_back(
        new Gel::OpGetDataAtGaussPts("HAT_EPS",common_data,true,false,MBTET)
      );
      fe_ptr[ss]->getOpPtrVector().push_back(
        new Gel::OpGetDataAtGaussPts("HAT_EPS_DOT",common_data,true,false,MBTET)
      );
    }

    // attach tags for each recorder
    vector<int> tags;
    tags.push_back(Gel::STRESSTOTAL); // ADOL-C tag used to calculate total stress
    tags.push_back(Gel::SOLVENTFLUX);
    tags.push_back(Gel::VOLUMERATE);
    tags.push_back(Gel::RESIDUALSTRAINHAT);

    // Right hand side operators
    gel.feRhs.getOpPtrVector().push_back(
      new Gel::OpJacobian(
        "SPATIAL_POSITION",tags,gel.constitutiveEquationPtr,gel.commonData,true,false
      )
    );
    gel.feRhs.getOpPtrVector().push_back(
      new Gel::OpRhsStressTotal(gel.commonData)
    );
    gel.feRhs.getOpPtrVector().push_back(
      new Gel::OpRhsSolventFlux(gel.commonData)
    );
    gel.feRhs.getOpPtrVector().push_back(
      new Gel::OpRhsVolumeDot(gel.commonData)
    );
    gel.feRhs.getOpPtrVector().push_back(
      new Gel::OpRhsStrainHat(gel.commonData)
    );

    // Left hand side operators
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpJacobian(
        "SPATIAL_POSITION",tags,gel.constitutiveEquationPtr,gel.commonData,false,true
      )
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdxdx(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdxdMu(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdxdStrainHat(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdStrainHatdStrainHat(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdStrainHatdx(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdMudMu(gel.commonData)
    );
    gel.feLhs.getOpPtrVector().push_back(
      new Gel::OpLhsdMudx(gel.commonData)
    );

  }

  // Create discrete manager instance
  DM dm;
  DMType dm_name = "DMGEL";
  {
    ierr = DMRegister_MoFEM(dm_name); CHKERRQ(ierr);
    ierr = DMCreate(PETSC_COMM_WORLD,&dm);CHKERRQ(ierr);
    ierr = DMSetType(dm,dm_name);CHKERRQ(ierr);
    ierr = DMMoFEMCreateMoFEM(dm,&m_field,dm_name,bit_level0); CHKERRQ(ierr);
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);
    //add elements to dm
    ierr = DMMoFEMAddElement(dm,"GEL_FE"); CHKERRQ(ierr);
    ierr = DMMoFEMAddElement(dm,"FORCE_FE"); CHKERRQ(ierr);
    ierr = DMMoFEMAddElement(dm,"PRESSURE_FE"); CHKERRQ(ierr);
    ierr = DMMoFEMAddElement(dm,"SOLVENT_FLUX_FE"); CHKERRQ(ierr);
    ierr = DMSetUp(dm); CHKERRQ(ierr);
  }

  // Create matrices and vectors used for analysis
  Vec T,F;
  Mat A;
  {
    ierr = DMCreateGlobalVector_MoFEM(dm,&T); CHKERRQ(ierr);
    ierr = VecDuplicate(T,&F); CHKERRQ(ierr);
    ierr = DMCreateMatrix_MoFEM(dm,&A); CHKERRQ(ierr);
  }

  // Setting finite element methods for Dirichelt boundary conditions
  SpatialPositionsBCFEMethodPreAndPostProc spatial_position_bc(
    m_field,"SPATIAL_POSITION",A,T,F
  );
  DirichletBCFromBlockSetFEMethodPreAndPostProc concentration_bc(
    m_field,"SOLVENT_CONCENTRATION","CONCENTRATION",A,T,F
  );

  // Setting finite element method for applying tractions
  boost::ptr_map<string,NeummanForcesSurface> neumann_forces;
  boost::ptr_map<string,NodalForce> nodal_forces;
  boost::ptr_map<string,EdgeForce> edge_forces;
  {
    //forces and pressures on surface
    ierr = MetaNeummanForces::setNeumannFiniteElementOperators(m_field,neumann_forces,F,"SPATIAL_POSITION"); CHKERRQ(ierr);
    //noadl forces
    ierr = MetaNodalForces::setOperators(m_field,nodal_forces,F,"SPATIAL_POSITION"); CHKERRQ(ierr);
    //edge forces
    ierr = MetaEdgeForces::setOperators(m_field,edge_forces,F,"SPATIAL_POSITION"); CHKERRQ(ierr);
  }

  // Solvent surface element, flux, convection radiation
  // TODO: add radiation and convection
  ThermalElement::MyTriFE solvent_surface_fe(m_field);
  {
    map<int,ThermalElement::FluxData>::iterator sit = set_of_solvent_fluxes.begin();
    for(;sit!=set_of_solvent_fluxes.end();sit++) {
      // add flux operator
      solvent_surface_fe.getOpPtrVector().push_back(
        new ThermalElement::OpHeatFlux("SOLVENT_CONCENTRATION",F,sit->second,true)
      );
    }
  }

  // Add finite elements to Time Stepping Solver, using Discrete Manager Interface
  {
    //Rhs
    ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&spatial_position_bc,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,&concentration_bc,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIFunction(dm,"GEL_FE",&gel.feRhs,NULL,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIFunction(dm,"SOLVENT_FLUX_FE",&solvent_surface_fe,NULL,NULL); CHKERRQ(ierr);
    {
      boost::ptr_map<string,NeummanForcesSurface>::iterator mit = neumann_forces.begin();
      for(;mit!=neumann_forces.end();mit++) {
        ierr = DMMoFEMTSSetIFunction(dm,mit->first.c_str(),&mit->second->getLoopFe(),NULL,NULL); CHKERRQ(ierr);
      }
    }
    {
      boost::ptr_map<string,NodalForce>::iterator fit = nodal_forces.begin();
      for(;fit!=nodal_forces.end();fit++) {
        ierr = DMMoFEMTSSetIFunction(dm,fit->first.c_str(),&fit->second->getLoopFe(),NULL,NULL); CHKERRQ(ierr);
      }
    }
    {
      boost::ptr_map<string,EdgeForce>::iterator fit = edge_forces.begin();
      for(;fit!=edge_forces.end();fit++) {
        ierr = DMMoFEMTSSetIFunction(dm,fit->first.c_str(),&fit->second->getLoopFe(),NULL,NULL); CHKERRQ(ierr);
      }
    }
    ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,NULL,&spatial_position_bc); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIFunction(dm,DM_NO_ELEMENT,NULL,NULL,&concentration_bc); CHKERRQ(ierr);

    //Lhs
    ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,&spatial_position_bc,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,&concentration_bc,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIJacobian(dm,"GEL_FE",&gel.feLhs,NULL,NULL); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,NULL,&spatial_position_bc); CHKERRQ(ierr);
    ierr = DMMoFEMTSSetIJacobian(dm,DM_NO_ELEMENT,NULL,NULL,&concentration_bc); CHKERRQ(ierr);

  }

  // Create Time Stepping solver
  TS ts;
  {
    ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
    ierr = TSSetType(ts,TSBEULER); CHKERRQ(ierr);
  }

  {
    ierr = DMoFEMMeshToLocalVector(dm,T,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = DMoFEMPreProcessFiniteElements(dm,&spatial_position_bc); CHKERRQ(ierr);
    ierr = DMoFEMPreProcessFiniteElements(dm,&concentration_bc); CHKERRQ(ierr);
    ierr = DMoFEMMeshToGlobalVector(dm,T,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  }

  // Solve problem
  {

    ierr = TSSetIFunction(ts,F,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);
    ierr = TSSetIJacobian(ts,A,A,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);

    //Monitor
    Gel::MonitorPostProc post_proc(m_field,"DMGEL","GEL_FE",gel.commonData);
    TsCtx *ts_ctx;
    DMMoFEMGetTsCtx(dm,&ts_ctx);
    {
      ts_ctx->get_postProcess_to_do_Monitor().push_back(&post_proc);
      ierr = TSMonitorSet(ts,f_TSMonitorSet,ts_ctx,PETSC_NULL); CHKERRQ(ierr);
    }

    double ftime = 1;
    ierr = TSSetDuration(ts,PETSC_DEFAULT,ftime); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);
    ierr = TSSetDM(ts,dm); CHKERRQ(ierr);
    ierr = TSSolve(ts,T); CHKERRQ(ierr);
    ierr = TSGetTime(ts,&ftime); CHKERRQ(ierr);
    PetscInt steps,snesfails,rejects,nonlinits,linits;
    ierr = TSGetTimeStepNumber(ts,&steps); CHKERRQ(ierr);
    ierr = TSGetSNESFailures(ts,&snesfails); CHKERRQ(ierr);
    ierr = TSGetStepRejections(ts,&rejects); CHKERRQ(ierr);
    ierr = TSGetSNESIterations(ts,&nonlinits); CHKERRQ(ierr);
    ierr = TSGetKSPIterations(ts,&linits); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,
      "steps %D (%D rejected, %D SNES fails), ftime %g, nonlinits %D, linits %D\n",
      steps,rejects,snesfails,ftime,nonlinits,linits
    );
  }

  // Clean and destroy
  {
    ierr = DMDestroy(&dm); CHKERRQ(ierr);
    ierr = VecDestroy(&T); CHKERRQ(ierr);
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
  }

  PetscFinalize();

  return 0;
}
