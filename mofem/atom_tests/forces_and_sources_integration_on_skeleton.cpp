/** \file build_problems.cpp

  \brief Testing integration on skeleton

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
  ierr = PetscOptionsGetString(PETSC_NULL,"","-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #else
  ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-my_file",mesh_file_name,255,&flg); CHKERRQ(ierr);
  #endif
  if(flg != PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF,1,"*** ERROR -my_file (MESH FILE NEEDED)");
  }

  ParallelComm* pcomm = ParallelComm::get_pcomm(&moab,MYPCOMM_INDEX);
  if(pcomm == NULL) pcomm =  new ParallelComm(&moab,PETSC_COMM_WORLD);

  const char *option;
  option = "";//"PARALLEL=BCAST;";//;DEBUG_IO";
  rval = moab.load_file(mesh_file_name, 0, option); CHKERRQ_MOAB(rval);

  //Create MoFEM (Joseph) database
  MoFEM::Core core(moab);
  MoFEM::Interface& m_field = core;

  //set entitities bit level
  BitRefLevel bit_level0;
  bit_level0.set(0);
  ierr = m_field.query_interface<BitRefManager>()->setBitRefLevelByDim(0,3,bit_level0); CHKERRQ(ierr);

  //Fields
  ierr = m_field.add_field("F2",HDIV,AINSWORTH_LEGENDRE_BASE,1); CHKERRQ(ierr);

  //meshset consisting all entities in mesh
  EntityHandle root_set = moab.get_root_set();
  //add entities to field
  ierr = m_field.add_ents_to_field_by_type(root_set,MBTET,"F2"); CHKERRQ(ierr);

  //set app. order
  //see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes (Mark Ainsworth & Joe Coyle)
  int order = 1;
  ierr = m_field.set_field_order(root_set,MBTET,"F2",order); CHKERRQ(ierr);
  ierr = m_field.set_field_order(root_set,MBTRI,"F2",order); CHKERRQ(ierr);

  ierr = m_field.build_fields(); CHKERRQ(ierr);

  //add elements
  ierr = m_field.add_finite_element("V1"); CHKERRQ(ierr);
  ierr = m_field.add_finite_element("S2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("V1","F2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("V1","F2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("V1","F2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_row("S2","F2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_col("S2","F2"); CHKERRQ(ierr);
  ierr = m_field.modify_finite_element_add_field_data("S2","F2"); CHKERRQ(ierr);

  ierr = m_field.add_ents_to_finite_element_by_type(root_set,MBTET,"V1"); CHKERRQ(ierr);
  Range faces;
  ierr = m_field.get_entities_by_type_and_ref_level(
    bit_level0,BitRefLevel().set(),MBTRI,faces
  ); CHKERRQ(ierr);
  ierr = m_field.add_ents_to_finite_element_by_type(faces,MBTRI,"S2"); CHKERRQ(ierr);

  ierr = m_field.build_finite_elements(); CHKERRQ(ierr);
  ierr = m_field.build_adjacencies(bit_level0); CHKERRQ(ierr);

  //Problems
  ierr = m_field.add_problem("P1"); CHKERRQ(ierr);

  //set refinement level for problem
  ierr = m_field.modify_problem_ref_level_add_bit("P1",bit_level0); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("P1","V1"); CHKERRQ(ierr);
  ierr = m_field.modify_problem_add_finite_element("P1","S2"); CHKERRQ(ierr);

  //build problems
  ProblemsManager *prb_mng_ptr;
  ierr = m_field.query_interface(prb_mng_ptr); CHKERRQ(ierr);
  ierr = prb_mng_ptr->buildProblem("P1",true); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionProblem("P1"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionFiniteElements("P1"); CHKERRQ(ierr);
  ierr = prb_mng_ptr->partitionGhostDofs("P1"); CHKERRQ(ierr);

  struct SkeletonFE: public FaceElementForcesAndSourcesCore::UserDataOperator {

    VolumeElementForcesAndSourcesCoreOnSide volSideFe;
    struct OpVolSide: public VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator {
      OpVolSide():
      VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator("F2",UserDataOperator::OPROW) {
      }
      PetscErrorCode doWork(int side, EntityType type,DataForcesAndSurcesCore::EntData &data) {
        PetscFunctionBegin;
        std::cout << "\tVolume" << getFEMethod()->nInTheLoop << std::endl;
        std::cout << "\tGauss pts " << getGaussPts() << std::endl;
        std::cout << "\tCoords " << getCoordsAtGaussPts() << endl;
        MatrixDouble diff = getCoordsAtGaussPts() - getFaceCoordsAtGaussPts();
        std::cout << std::fixed << std::setprecision(3) << "\tDiff coords " << diff << endl;
        const double eps = 1e-12;
        if(norm_inf(diff)>eps) {
          SETERRQ(
            PETSC_COMM_WORLD,
            MOFEM_ATOM_TEST_INVALID,
            "coordinates at integration pts are different"
          );
        }
        PetscFunctionReturn(0);
      }
    };

    SkeletonFE(MoFEM::Interface &m_field):
    FaceElementForcesAndSourcesCore::UserDataOperator("F2",UserDataOperator::OPROW),
    volSideFe(m_field) {
      volSideFe.getOpPtrVector().push_back(new SkeletonFE::OpVolSide());
    }

    int getRule(int order) { return order; };

    PetscErrorCode doWork(int side,EntityType type,DataForcesAndSurcesCore::EntData &data) {
      
      PetscFunctionBegin;
      if(type == MBTRI && side == 0) {
        std::cout << "Face" << std::endl;
        std::cout << "Gauss pts " << getGaussPts() << std::endl;
        std::cout << "Coords " << getCoordsAtGaussPts() << endl;
        ierr = loopSideVolumes("V1",volSideFe); CHKERRQ(ierr);
      }
      PetscFunctionReturn(0);
    }

  };

  FaceElementForcesAndSourcesCore face_fe(m_field);
  face_fe.getOpPtrVector().push_back(new SkeletonFE(m_field));
  ierr = m_field.loop_finite_elements("P1","S2",face_fe);  CHKERRQ(ierr);


  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF,e.errorCode,e.errorMessage);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;

}
