/** \file forces_and_sources_integration_on_skeleton.cpp
 * \example forces_and_sources_integration_on_skeleton.cpp
 * \brief Testing integration on skeleton
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

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    const char *option;
    option = ""; 
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM database
    MoFEM::Core core(moab);
    // Access to database through interface
    MoFEM::Interface &m_field = core;

    // set entireties bit level
    BitRefLevel bit_level0 = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("F2", HDIV, AINSWORTH_LEGENDRE_BASE, 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F2");

    // set app. order
    int order = 1;
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order);

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("V1");
    CHKERR m_field.add_finite_element("S2");
    CHKERR m_field.modify_finite_element_add_field_row("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_data("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_row("S2", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("S2", "F2");
    CHKERR m_field.modify_finite_element_add_field_data("S2", "F2");

    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "V1");
    Range faces;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit_level0, BitRefLevel().set(), MBTRI, faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(faces, MBTRI, "S2");

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_level0);

    // Problems
    CHKERR m_field.add_problem("P1");

    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_level0);
    CHKERR m_field.modify_problem_add_finite_element("P1", "V1");
    CHKERR m_field.modify_problem_add_finite_element("P1", "S2");

    // build problems
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("P1", true);
    CHKERR prb_mng_ptr->partitionProblem("P1");
    CHKERR prb_mng_ptr->partitionFiniteElements("P1");
    CHKERR prb_mng_ptr->partitionGhostDofs("P1");

    struct SkeletonFE
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      VolumeElementForcesAndSourcesCoreOnSide volSideFe;
      struct OpVolSide
          : public VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator {
        OpVolSide()
            : VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator(
                  "F2", UserDataOperator::OPROW) {}
        MoFEMErrorCode doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
          MoFEMFunctionBeginHot;
          std::cout << "\tVolume" << getFEMethod()->nInTheLoop << std::endl;
          std::cout << "\tGauss pts " << getGaussPts() << std::endl;
          std::cout << "\tCoords " << getCoordsAtGaussPts() << endl;
          MatrixDouble diff = getCoordsAtGaussPts() - getFaceCoordsAtGaussPts();
          std::cout << std::fixed << std::setprecision(3) << "\tDiff coords "
                    << diff << endl;
          const double eps = 1e-12;
          if (norm_inf(diff) > eps) {
            SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                    "coordinates at integration pts are different");
          }
          MoFEMFunctionReturnHot(0);
        }
      };

      SkeletonFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "F2", UserDataOperator::OPROW),
            volSideFe(m_field) {
        volSideFe.getOpPtrVector().push_back(new SkeletonFE::OpVolSide());
      }

      int getRule(int order) { return order; };

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {

        MoFEMFunctionBeginHot;
        if (type == MBTRI && side == 0) {
          std::cout << "Face" << std::endl;
          std::cout << "Gauss pts " << getGaussPts() << std::endl;
          std::cout << "Coords " << getCoordsAtGaussPts() << endl;
          CHKERR loopSideVolumes("V1", volSideFe);
        }
        MoFEMFunctionReturnHot(0);
      }
    };

    FaceElementForcesAndSourcesCore face_fe(m_field);
    face_fe.getOpPtrVector().push_back(new SkeletonFE(m_field));
    CHKERR m_field.loop_finite_elements("P1", "S2", face_fe);

  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
