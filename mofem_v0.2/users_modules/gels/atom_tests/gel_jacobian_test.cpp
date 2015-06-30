/*s* \file gel_jacobian_test.cpp
  \brief Atom test testing calculation of element residual vectors and tangent matrices.
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

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <adolc/adolc.h>
#include <Gels.hpp>

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

int main(int argc, char *argv[]) {

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  Interface& moab = mb_instance;

  {
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
  }

  MoFEM::Core core(moab);
  FieldInterface& m_field = core;
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  // Define fields and finite elements
  {

    // Set approximation fields
    {
      // Seed all mesh entities to MoFEM database, those entities can be potentially used as finite elements
      // or as entities which carry some approximation field.

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
      //add entities to field (root_mesh, i.e. on all mesh entities fields are approx.)
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SPATIAL_POSITION_DOT"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"SOLVENT_CONCENTRATION_DOT"); CHKERRQ(ierr);
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"HAT_EPS_DOT"); CHKERRQ(ierr);

      PetscBool flg;
      PetscInt order;
      ierr = PetscOptionsGetInt(PETSC_NULL,"-my_order",&order,&flg); CHKERRQ(ierr);
      if(flg != PETSC_TRUE) {
        order = 2;
      }
      if(order < 2) {
        //SETERRQ()
      }

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

      ierr = m_field.set_field_order(root_set,MBTET,"HAT_EPS",order-1); CHKERRQ(ierr);
      ierr = m_field.set_field_order(root_set,MBTET,"HAT_EPS_DOT",order-1); CHKERRQ(ierr);

      //gemetry approximation is set to 2nd oreder
      ierr = m_field.add_ents_to_field_by_TETs(root_set,"MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTET,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBTRI,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBEDGE,"MESH_NODE_POSITIONS",2); CHKERRQ(ierr);
      ierr = m_field.set_field_order(0,MBVERTEX,"MESH_NODE_POSITIONS",1); CHKERRQ(ierr);

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

    //Set finite elements
    {
      ierr = m_field.add_finite_element("GEL_FE",MF_ZERO); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SPATIAL_POSITION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SPATIAL_POSITION_DOT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","SOLVENT_CONCENTRATION"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_row("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_col("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","HAT_EPS"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","HAT_EPS_DOT"); CHKERRQ(ierr);
      ierr = m_field.modify_finite_element_add_field_data("GEL_FE","MESH_NODE_POSITIONS"); CHKERRQ(ierr);
      EntityHandle root_set = moab.get_root_set();
      ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"GEL_FE"); CHKERRQ(ierr);

      //build finite elemnts
      ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
      //build adjacencies
      ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);
    }

  }

  // Create gel instance
  Gel gel(m_field);
  {

    Gel::BlockMaterialData &material_data = gel.blockMaterialData;

    // Set material parameters
    material_data.gAlpha = 1;
    material_data.vAlpha = 0.3;
    material_data.gBeta = 1;
    material_data.vBeta = 0.3;
    material_data.gBetaHat = 1;
    material_data.vBetaHat = 0.3;
    material_data.oMega = 1;
    material_data.vIscosity = 1;
    material_data.pErmeability = 2.;

    Gel::CommonData &common_data = gel.commonData;
    common_data.spatialPositionName = "SPATIAL_POSITION";
    common_data.spatialPositionNameDot = "SPATIAL_POSITION_DOT";
    common_data.muName = "SOLVENT_CONCENTRATION";
    common_data.strainHatName = "HAT_EPS";
    common_data.strainHatNameDot = "HAT_EPS_DOT";

    Gel::GelFE *fe_ptr[] = { &gel.feRhs, &gel.feLhs };
    for(int ss = 0;ss<2;ss++) {
      fe_ptr[ss]->getOpPtrVector().push_back(new Gel::OpGetDataAtGaussPts("SPATIAL_POSITION",common_data,false,true));
      fe_ptr[ss]->getOpPtrVector().push_back(new Gel::OpGetDataAtGaussPts("SPATIAL_POSITION_DOT",common_data,false,true));
      fe_ptr[ss]->getOpPtrVector().push_back(new Gel::OpGetDataAtGaussPts("SOLVENT_CONCENTRATION",common_data,true,true));
      fe_ptr[ss]->getOpPtrVector().push_back(new Gel::OpGetDataAtGaussPts("HAT_EPS",common_data,true,false,MBTET));
      fe_ptr[ss]->getOpPtrVector().push_back(new Gel::OpGetDataAtGaussPts("HAT_EPS_DOT",common_data,true,false,MBTET));

      // attach tags for each recorder
      vector<int> tags;
      tags.push_back(1);
      tags.push_back(2);
      tags.push_back(3);
      tags.push_back(4);

      gel.feRhs.getOpPtrVector().push_back(
        new Gel::OpJacobian("SPATIAL_POSITION", tags, gel.constitutiveEquation,gel.commonData,true,false)
      );
      gel.feRhs.getOpPtrVector().push_back(
        new Gel::OpRhsAssembleStrainTotal(gel.commonData)
      );
      gel.feRhs.getOpPtrVector().push_back(
        new Gel::OpRhsAssembleSolventFlux(gel.commonData)
      );
      gel.feLhs.getOpPtrVector().push_back(
        new Gel::OpJacobian("SPATIAL_POSITION", tags, gel.constitutiveEquation, gel.commonData,false,true)
      );
    }
  }

  // Create dm instance
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
    ierr = DMSetUp(dm); CHKERRQ(ierr);
  }

  // Make calculations
  Vec F;
  {
    ierr = DMCreateGlobalVector_MoFEM(dm,&F); CHKERRQ(ierr);
    gel.feRhs.snes_f = F; // Set right hand side vector manually
    ierr = DMoFEMLoopFiniteElements(dm,"GEL_FE",&gel.feRhs); CHKERRQ(ierr);
    ierr = DMoFEMLoopFiniteElements(dm,"GEL_FE",&gel.feLhs); CHKERRQ(ierr);
  }


  // Clean and destroy
  {
    ierr = VecDestroy(&F); CHKERRQ(ierr);
    ierr = DMDestroy(&dm); CHKERRQ(ierr);
  }

  PetscFinalize();

  return 0;
}
