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
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("HCURL", HCURL, AINSWORTH_LEGENDRE_BASE, 1);

    // FE TET
    CHKERR m_field.add_finite_element("HCURL_TET_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("HCURL_TET_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_col("HCURL_TET_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_data("HCURL_TET_FE",
                                                        "HCURL");

    // FE TRI
    CHKERR m_field.add_finite_element("HCURL_TRI_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("HCURL_TRI_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_col("HCURL_TRI_FE", "HCURL");
    CHKERR m_field.modify_finite_element_add_field_data("HCURL_TRI_FE",
                                                        "HCURL");

    // FE EDGE
    CHKERR m_field.add_finite_element("HCURL_EDGE_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("HCURL_EDGE_FE",
                                                       "HCURL");
    CHKERR m_field.modify_finite_element_add_field_col("HCURL_EDGE_FE",
                                                       "HCURL");
    CHKERR m_field.modify_finite_element_add_field_data("HCURL_EDGE_FE",
                                                        "HCURL");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "HCURL_TET_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "HCURL_TRI_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM",
                                                     "HCURL_EDGE_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "HCURL");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "HCURL_TET_FE");

    Range tets;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        BitRefLevel().set(0), BitRefLevel().set(), MBTET, tets);
    Skinner skin(&moab);
    Range skin_faces; // skin faces from 3d ents
    CHKERR skin.find_skin(0, tets, false, skin_faces);
    CHKERR m_field.add_ents_to_finite_element_by_type(skin_faces, MBTRI,
                                                      "HCURL_TRI_FE");
    Range skin_edges;
    CHKERR moab.get_adjacencies(skin_faces, 1, false, skin_edges,
                                moab::Interface::UNION);
    CHKERR m_field.add_ents_to_finite_element_by_type(skin_edges, MBEDGE,
                                                      "HCURL_EDGE_FE");

    // set app. order
    int order = 4;
    CHKERR m_field.set_field_order(root_set, MBTET, "HCURL", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "HCURL", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "HCURL", order);

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

    // mesh partitioning

    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_hcurl_approximation_functions.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct OpPrintingHdivApproximationFunctions
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpPrintingHdivApproximationFunctions(TeeStream &my_split)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "HCURL", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << std::endl
                << "type " << type << " side " << side << std::endl;

        const double eps = 1e-4;
        for (unsigned int dd = 0; dd < data.getN().data().size(); dd++) {
          if (std::abs(data.getN().data()[dd]) < eps)
            data.getN().data()[dd] = 0;
        }
        for (unsigned int dd = 0; dd < data.getDiffN().data().size(); dd++) {
          if (std::abs(data.getDiffN().data()[dd]) < eps)
            data.getDiffN().data()[dd] = 0;
        }

        mySplit << std::fixed << std::setprecision(5) << data.getN()
                << std::endl;
        mySplit << std::fixed << std::setprecision(5) << data.getDiffN()
                << std::endl;

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
                "HCURL", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << std::endl
                << "type " << type << " side " << side << std::endl;
        mySplit.precision(5);

        const double eps = 1e-4;
        for (unsigned int dd = 0; dd < data.getN().data().size(); dd++) {
          if (std::abs(data.getN().data()[dd]) < eps)
            data.getN().data()[dd] = 0;
        }
        for (unsigned int dd = 0; dd < data.getDiffN().data().size(); dd++) {
          if (std::abs(data.getDiffN().data()[dd]) < eps)
            data.getDiffN().data()[dd] = 0;
        }

        mySplit << std::fixed << std::setprecision(5) << data.getN()
                << std::endl;
        mySplit << std::fixed << std::setprecision(5) << data.getDiffN()
                << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyFaceFE : public FaceElementForcesAndSourcesCore {

      MyFaceFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 1; };
    };

    struct OpEdgePrintingHdivApproximationFunctions
        : public EdgeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpEdgePrintingHdivApproximationFunctions(TeeStream &my_split)
          : EdgeElementForcesAndSourcesCore::UserDataOperator(
                "HCURL", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBeginHot;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << std::endl
                << "type " << type << " side " << side << std::endl;
        mySplit.precision(5);

        const double eps = 1e-4;
        for (unsigned int dd = 0; dd < data.getN().data().size(); dd++) {
          if (std::abs(data.getN().data()[dd]) < eps)
            data.getN().data()[dd] = 0;
        }

        mySplit << std::fixed << std::setprecision(5) << data.getN()
                << std::endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    struct MyEdgeFE : public EdgeElementForcesAndSourcesCore {

      MyEdgeFE(MoFEM::Interface &m_field)
          : EdgeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return 1; };
    };

    MyFE tet_fe(m_field);
    MyFaceFE tri_fe(m_field);
    MyEdgeFE edge_fe(m_field);

    tet_fe.getOpPtrVector().push_back(
        new OpPrintingHdivApproximationFunctions(my_split));
    tri_fe.getOpPtrVector().push_back(
        new OpHOSetCovariantPiolaTransformOnFace3D(HCURL));
    tri_fe.getOpPtrVector().push_back(
        new OpFacePrintingHdivApproximationFunctions(my_split));
    edge_fe.getOpPtrVector().push_back(
        new OpEdgePrintingHdivApproximationFunctions(my_split));

    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "HCURL_TET_FE", tet_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "HCURL_TRI_FE", tri_fe);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "HCURL_EDGE_FE",
                                        edge_fe);

    /*PostProcVolumeOnRefinedMesh post_proc(m_field);
    CHKERR post_proc.generateReferenceElementMesh();
    CHKERR post_proc.addHdivFunctionsPostProc("HCURL");
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM","HCURL_TET_FE",post_proc);
    CHKERR post_proc.postProcMesh.write_file("out.vtk","VTK",""); */
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
