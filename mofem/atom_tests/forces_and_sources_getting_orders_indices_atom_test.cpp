/** \file forces_and_sources_getting_orders_indices_atom_test.cpp
  \brief Atome test for getting orders on entities
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

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  ErrorCode rval;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  BARRIER_RANK_START(pcomm)
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);
  BARRIER_RANK_END(pcomm)

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRQ_MOAB(rval);
  ierr = m_field.seed_ref_level_3D(0,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,1); CHKERRQ(ierr);
  ierr = m_field.add_field("FIELD2",H1,3); CHKERRQ(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRQ(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","FIELD2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD2"); CHKERRQ(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRQ(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRQ(ierr);
  //set refinment level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRQ(ierr);


  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD1"); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_field_by_TETs(root_set,"FIELD2"); CHKERRQ(ierr);
  //add entities to finite element
  ierr = m_field.add_ents_to_finite_element_by_TETs(root_set,"TEST_FE"); CHKERRQ(ierr);


  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 5;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD2",1); CHKERRQ(ierr);

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
  ierr = m_field.partition_simple_problem("TEST_PROBLEM"); CHKERRQ(ierr);
  ierr = m_field.partition_finite_elements("TEST_PROBLEM"); CHKERRQ(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = m_field.partition_ghost_dofs("TEST_PROBLEM"); CHKERRQ(ierr);

  struct ForcesAndSurcesCore_TestFE: public ForcesAndSurcesCore {

    ErrorCode rval;
    PetscErrorCode ierr;

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs;
    TeeDevice my_tee;
    TeeStream my_split;

    DataForcesAndSurcesCore data;
    DerivedDataForcesAndSurcesCore derived_data;

    ForcesAndSurcesCore_TestFE(MoFEM::Interface &_m_field):
      ForcesAndSurcesCore(_m_field),
      ofs("forces_and_sources_getting_orders_indices_atom_test.txt"),
      my_tee(std::cout, ofs),my_split(my_tee),
      data(MBTET),
      derived_data(data) {};

    PetscErrorCode preProcess() {
      PetscFunctionBegin;
      PetscFunctionReturn(0);
    }

    PetscErrorCode operator()() {
      PetscFunctionBegin;

      my_split << "\n\nNEXT ELEM\n\n";

      ierr = getSpacesAndBaseOnEntities(data); CHKERRQ(ierr);

      ierr = getEdgesSense(data); CHKERRQ(ierr);
      ierr = getTrisSense(data); CHKERRQ(ierr);
      ierr = getEdgesDataOrder(data,H1); CHKERRQ(ierr);
      ierr = getTrisDataOrder(data,H1); CHKERRQ(ierr);
      ierr = getTetDataOrder(data,H1); CHKERRQ(ierr);
      ierr = getFaceTriNodes(data); CHKERRQ(ierr);

      MatrixDouble gauss_pts(4,4);
      for(int gg = 0;gg<4;gg++) {
        gauss_pts(0,gg) = G_TET_X4[gg];
        gauss_pts(1,gg) = G_TET_Y4[gg];
        gauss_pts(2,gg) = G_TET_Z4[gg];
        gauss_pts(3,gg) = G_TET_W4[gg];
      }
      ierr = TetPolynomialBase().getValue(
        gauss_pts,
        boost::shared_ptr<BaseFunctionCtx>(
          new EntPolynomialBaseCtx(data,H1,AINSWORTH_COLE_BASE)
        )
      ); CHKERRQ(ierr);

      ierr = getEdgesDataOrderSpaceAndBase(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTrisDataOrderSpaceAndBase(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTetDataOrderSpaceAndBase(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getRowNodesIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getEdgesRowIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTrisRowIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getTetsRowIndices(data,"FIELD1"); CHKERRQ(ierr);
      ierr = getNodesFieldData(data,"FIELD1"); CHKERRQ(ierr);
      data.dataOnEntities[MBVERTEX][0].getFieldData().resize(0);

      my_split << "FIELD1:\n";
      my_split << data << std::endl;

      derived_data.dataOnEntities[MBVERTEX][0].getBase() = AINSWORTH_COLE_BASE;
      ierr = getEdgesDataOrderSpaceAndBase(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getTrisDataOrderSpaceAndBase(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getTetDataOrderSpaceAndBase(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getColNodesIndices(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getEdgesColIndices(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getTrisColIndices(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getTetsColIndices(derived_data,"FIELD2"); CHKERRQ(ierr);
      ierr = getNodesFieldData(derived_data,"FIELD2"); CHKERRQ(ierr);
      derived_data.dataOnEntities[MBVERTEX][0].getFieldData().resize(0);

      my_split << "FIELD2:\n";
      my_split << derived_data << std::endl;

      PetscFunctionReturn(0);
    }

    PetscErrorCode postProcess() {
      PetscFunctionBegin;

      my_split.close();

      PetscFunctionReturn(0);
    }

  };

  ForcesAndSurcesCore_TestFE fe1(m_field);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",fe1);  CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
