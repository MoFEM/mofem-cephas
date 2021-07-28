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
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // create one tet
    double tet_coords[] = {0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2};
    EntityHandle nodes[4];
    for (int nn = 0; nn < 4; nn++) {
      CHKERR moab.create_vertex(&tet_coords[3 * nn], nodes[nn]);
    }
    EntityHandle tet;
    CHKERR moab.create_element(MBTET, nodes, 4, tet);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("HDIV", HDIV, AINSWORTH_LEGENDRE_BASE, 1);

    // FE TET
    CHKERR m_field.add_finite_element("HDIV_TET_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("HDIV_TET_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_col("HDIV_TET_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_data("HDIV_TET_FE", "HDIV");

    // FE TRI
    CHKERR m_field.add_finite_element("HDIV_TRI_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("HDIV_TRI_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_col("HDIV_TRI_FE", "HDIV");
    CHKERR m_field.modify_finite_element_add_field_data("HDIV_TRI_FE", "HDIV");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "HDIV_TET_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "HDIV_TRI_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "HDIV");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "HDIV_TET_FE");

    Range tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTET, tets);
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    CHKERR skin.find_skin(0, tets, false, skin_faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(skin_faces, MBTRI,
                                                      "HDIV_TRI_FE");

    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "HDIV", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "HDIV", order);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    // build problem

    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);

    /****/
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_hdiv_approximation_functions.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct OpPrintingHdivApproximationFunctions
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpPrintingHdivApproximationFunctions(TeeStream &my_split)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "HDIV", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << std::endl
                << "type " << type << " side " << side << std::endl;
        mySplit.precision(5);

        const double eps = 1e-6;
        for (unsigned int dd = 0; dd < data.getN().data().size(); dd++) {
          if (fabs(data.getN().data()[dd]) < eps)
            data.getN().data()[dd] = 0;
        }
        for (unsigned int dd = 0; dd < data.getDiffN().data().size();
             dd++) {
          if (fabs(data.getDiffN().data()[dd]) < eps)
            data.getDiffN().data()[dd] = 0;
        }

        mySplit << std::fixed << data.getN() << std::endl;
        mySplit << std::fixed << data.getDiffN() << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyFE : public VolumeElementForcesAndSourcesCore {

      MyFE(MoFEM::Interface &m_field)
          : VolumeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 1; };
    };

    struct OpFacePrintingHdivApproximationFunctions
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpFacePrintingHdivApproximationFunctions(TeeStream &my_split)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "HDIV", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << std::endl
                << "type " << type << " side " << side << std::endl;
        mySplit.precision(5);

        const double eps = 1e-6;
        for (unsigned int dd = 0; dd < data.getN().data().size(); dd++) {
          if (fabs(data.getN().data()[dd]) < eps)
            data.getN().data()[dd] = 0;
        }
        for (unsigned int dd = 0; dd < data.getDiffN().data().size();
             dd++) {
          if (fabs(data.getDiffN().data()[dd]) < eps)
            data.getDiffN().data()[dd] = 0;
        }

        mySplit << std::fixed << data.getN() << std::endl;
        // mySplit << std::fixed << data.getDiffN() << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyFaceFE : public FaceElementForcesAndSourcesCore {

      MyFaceFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 1; };
    };

    MyFE tet_fe(m_field);
    MyFaceFE tri_fe(m_field);

    tet_fe.getOpPtrVector().push_back(
        new OpPrintingHdivApproximationFunctions(my_split));

    tri_fe.getOpPtrVector().push_back(
        new OpHOSetContravariantPiolaTransformOnFace(HDIV));
    tri_fe.getOpPtrVector().push_back(
        new OpFacePrintingHdivApproximationFunctions(my_split));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "HDIV_TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "HDIV_TRI_FE", tri_fe);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
