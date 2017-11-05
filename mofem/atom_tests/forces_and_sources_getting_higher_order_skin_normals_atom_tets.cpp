/** \file  forces_and_sources_getting_higher_order_skin_normals_atom_tets.cpp

  \brief Atom test to calculate normals on faces approximated with HO geometry representation

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

  
  

  PetscInitialize(&argc,&argv,(char *)0,help);

  try {

  moab::Core mb_instance;
  moab::Interface& moab = mb_instance;
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscBool flg = PETSC_TRUE;
  char mesh_file_name[255];
  #if PETSC_VERSION_GE(3,6,4)
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRG(ierr);
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
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRG(rval);

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  EntityHandle meshset_level0;
  rval = moab.create_meshset(MESHSET_SET,meshset_level0); CHKERRG(rval);
  ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRG(ierr);

  //Fields
  ierr = m_field.add_field("FIELD1",H1,AINSWORTH_LEGENDRE_BASE,3); CHKERRG(ierr);

  //FE
  ierr = m_field.add_finite_element("TEST_FE"); CHKERRG(ierr);

  //Define rows/cols and element data
  ierr = m_field.modify_finite_element_add_field_row("TEST_FE","FIELD1"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_col("TEST_FE","FIELD1"); CHKERRG(ierr);
  ierr = m_field.modify_finite_element_add_field_data("TEST_FE","FIELD1"); CHKERRG(ierr);

  //Problem
  ierr = m_field.add_problem("TEST_PROBLEM"); CHKERRG(ierr);

  //set finite elements for problem
  ierr = m_field.modify_problem_add_finite_element("TEST_PROBLEM","TEST_FE"); CHKERRG(ierr);
  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM",bit_level0); CHKERRG(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"FIELD1"); CHKERRG(ierr);
  //add entities to finite element
  Range tets;
  rval = moab.get_entities_by_type(0,MBTET,tets,false); CHKERRG(rval);
  Skinner skin(&m_field.get_moab());
  Range tets_skin;
  rval = skin.find_skin(0,tets,false,tets_skin); CHKERRG(rval);
  ierr = m_field.add_ents_to_finite_element_by_type(tets_skin,MBTRI,"TEST_FE"); CHKERRG(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 3;
  ierr = m_field.set_field_order(root_set,MBTET,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBEDGE,"FIELD1",order); CHKERRG(ierr);
  ierr = m_field.set_field_order(root_set,MBVERTEX,"FIELD1",1); CHKERRG(ierr);

  /****/
  //build database
  //build field
  ierr = m_field.build_fields(); CHKERRG(ierr);
  //set FIELD1 from positions of 10 node tets
  Projection10NodeCoordsOnField ent_method(m_field,"FIELD1");
  ierr = m_field.loop_dofs("FIELD1",ent_method); CHKERRG(ierr);
  //build finite elemnts
  ierr = m_field.build_finite_elements(); CHKERRG(ierr);
  //build adjacencies
  ierr = m_field.build_adjacencies(bit_level0); CHKERRG(ierr);
  //build problem
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.getInterface(prb_mng_ptr); CHKERRG(ierr);
  ierr = prb_mng_ptr->buildProblem("TEST_PROBLEM",true); CHKERRG(ierr);

  /****/
  //mesh partitioning
  //partition
  ierr = prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM"); CHKERRG(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM"); CHKERRG(ierr);
  //what are ghost nodes, see Petsc Manual
  ierr = prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM"); CHKERRG(ierr);

  struct ForcesAndSourcesCore_TestFE: public ForcesAndSourcesCore {

    
    

    MatrixDouble hoCoords_at_GaussPt;
    MatrixDouble nOrmals_at_GaussPt;
    MatrixDouble tAngent1_at_GaussPt;
    MatrixDouble tAngent2_at_GaussPt;
    DataForcesAndSourcesCore data;
    OpGetCoordsAndNormalsOnFace op;

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;
    std::ofstream ofs;
    TeeDevice my_tee;
    TeeStream my_split;

    ForcesAndSourcesCore_TestFE(MoFEM::Interface &_m_field):
      ForcesAndSourcesCore(_m_field),data(MBTRI),
      op(
        hoCoords_at_GaussPt,nOrmals_at_GaussPt,tAngent1_at_GaussPt,tAngent2_at_GaussPt
      ),
      ofs("forces_and_sources_getting_higher_order_skin_normals_atom.txt"),
      my_tee(std::cout,ofs),my_split(my_tee) {};

    MoFEMErrorCode preProcess() {
      MoFEMFunctionBeginHot;
      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode operator()() {
      MoFEMFunctionBeginHot;

      ierr = getSpacesAndBaseOnEntities(data); CHKERRG(ierr);

      ierr = getEdgesSense(data); CHKERRG(ierr);
      ierr = getEdgesDataOrder(data,H1); CHKERRG(ierr);
      ierr = getTrisDataOrder(data,H1); CHKERRG(ierr);
      ierr = getEdgesDataOrderSpaceAndBase(data,"FIELD1"); CHKERRG(ierr);
      ierr = getTrisDataOrderSpaceAndBase(data,"FIELD1"); CHKERRG(ierr);
      ierr = getRowNodesIndices(data,"FIELD1"); CHKERRG(ierr);
      ierr = getEdgesRowIndices(data,"FIELD1"); CHKERRG(ierr);
      ierr = getTrisRowIndices(data,"FIELD1"); CHKERRG(ierr);
      ierr = getNodesFieldData(data,"FIELD1"); CHKERRG(ierr);
      ierr = getEdgesFieldData(data,"FIELD1"); CHKERRG(ierr);
      ierr = getTrisFieldData(data,"FIELD1"); CHKERRG(ierr);

      data.dataOnEntities[MBVERTEX][0].getN(NOBASE).resize(4,3,false);
      ierr = ShapeMBTRI(
        &*data.dataOnEntities[MBVERTEX][0].getN(NOBASE).data().begin(),
        G_TRI_X4,G_TRI_Y4,4
      ); CHKERRG(ierr);
      data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).resize(3,2,false);
      ierr = ShapeDiffMBTRI(
        &*data.dataOnEntities[MBVERTEX][0].getDiffN(NOBASE).data().begin()
      ); CHKERRG(ierr);

      MatrixDouble gauss_pts(2,4,false);
      for(int gg = 0;gg<4;gg++) {
        gauss_pts(0,gg) = G_TRI_X4[gg];
        gauss_pts(1,gg) = G_TRI_Y4[gg];
      }

      ierr = TriPolynomialBase().getValue(
        gauss_pts,
        boost::shared_ptr<BaseFunctionCtx>(
          new EntPolynomialBaseCtx(data,H1,AINSWORTH_LEGENDRE_BASE,NOBASE)
        )
      ); CHKERRG(ierr);

      nOrmals_at_GaussPt.resize(4,3);
      tAngent1_at_GaussPt.resize(4,3);
      tAngent2_at_GaussPt.resize(4,3);

      try {
        ierr = op.opRhs(data); CHKERRG(ierr);
        ierr = op.calculateNormals(); CHKERRG(ierr);
      } catch (std::exception& ex) {
        std::ostringstream ss;
        ss << "thorw in method: " << ex.what() << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
        SETERRQ(PETSC_COMM_SELF,1,ss.str().c_str());
      }

      my_split.precision(3);
      my_split << "coords: " << hoCoords_at_GaussPt << std::endl;
      my_split << "normals: " << nOrmals_at_GaussPt << std::endl;
      my_split << "tangent1: " << tAngent1_at_GaussPt << std::endl;
      my_split << "tangent2: " << tAngent2_at_GaussPt << std::endl;

      MoFEMFunctionReturnHot(0);
    }

    MoFEMErrorCode postProcess() {
      MoFEMFunctionBeginHot;
      MoFEMFunctionReturnHot(0);
    }


  };

  ForcesAndSourcesCore_TestFE fe1(m_field);
  ierr = m_field.loop_finite_elements("TEST_PROBLEM","TEST_FE",fe1);  CHKERRG(ierr);

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRG(ierr);

  return 0;

}
