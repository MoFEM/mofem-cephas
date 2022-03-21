/** \file forces_and_sources_calculate_jacobian.cpp

  \brief Atom test checking with blessed file how Jacobian's are calculated on
  elements

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

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const char *option;
    option = ""; 
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);

    // FE
    CHKERR m_field.add_finite_element("TEST_FE");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE", "FIELD1");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TEST_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD1");
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET,
                                                      "TEST_FE");

    // set app. order
    // see Hierarchic Finite Element Bases on Unstructured Tetrahedral Meshes
    // (Mark Ainsworth & Joe Coyle)
    int order = 5;
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD1", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD1", 1);

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

    struct ForcesAndSourcesCore_TestFE : public ForcesAndSourcesCore {

      typedef tee_device<std::ostream, std::ofstream> TeeDevice;
      typedef stream<TeeDevice> TeeStream;

      std::ofstream ofs;
      TeeDevice my_tee;
      TeeStream my_split;

      struct PrintJacobian : public DataOperator {

        TeeStream &my_split;
        PrintJacobian(TeeStream &_my_split) : my_split(_my_split) {}

        ~PrintJacobian() { my_split.close(); }

        MoFEMErrorCode doWork(int side, EntityType type,
                              EntitiesFieldData::EntData &data) {
          MoFEMFunctionBeginHot;
          const double eps = 1e-6;
          for (unsigned int ii = 0; ii != data.getDiffN().size1(); ii++) {
            for (unsigned int jj = 0; jj != data.getDiffN().size2(); jj++) {
              if (fabs(data.getDiffN()(ii, jj)) < eps) {
                data.getDiffN()(ii, jj) = 0;
              }
            }
          }
          my_split << "side: " << side << " type: " << type << std::fixed
                   << std::setprecision(4) << data.getDiffN() << std::endl;
          MoFEMFunctionReturnHot(0);
        }
      };

      MatrixDouble3by3 invJac;
      MatrixDouble dataFIELD1;
      MatrixDouble dataDiffFIELD1;
      VectorDouble coords;
      PrintJacobian opPrintJac;
      OpSetInvJacH1 opSetInvJac;
      OpGetDataAndGradient<1, 3> opGetData_FIELD1;

      ForcesAndSourcesCore_TestFE(MoFEM::Interface &_m_field)
          : ForcesAndSourcesCore(_m_field),
            ofs("forces_and_sources_calculate_jacobian.txt"),
            my_tee(std::cout, ofs), my_split(my_tee), invJac(3, 3),
            opPrintJac(my_split), opSetInvJac(invJac),
            opGetData_FIELD1(dataFIELD1, dataDiffFIELD1), data(MBTET) {}

      MoFEMErrorCode preProcess() {
        MoFEMFunctionBeginHot;
        MoFEMFunctionReturnHot(0);
      }

      EntitiesFieldData data;

      MoFEMErrorCode operator()() {
        MoFEMFunctionBegin;

        CHKERR getSpacesAndBaseOnEntities(data);

        CHKERR getEntitySense<MBEDGE>(data);
        CHKERR getEntitySense<MBTRI>(data);
        CHKERR getEntityDataOrder<MBEDGE>(data, H1);
        CHKERR getEntityDataOrder<MBTRI>(data, H1);
        CHKERR getEntityDataOrder<MBTET>(data, H1);
        CHKERR getFaceNodes(data);

        CHKERR getRowNodesIndices(data, "FIELD1");
        CHKERR getEntityRowIndices(data, "FIELD1", MBEDGE);

        CHKERR getNodesFieldData(data, "FIELD1");
        CHKERR getEntityFieldData(data, "FIELD1", MBEDGE);

        MatrixDouble gauss_pts(4, 4);
        for (int gg = 0; gg < 4; gg++) {
          gauss_pts(0, gg) = G_TET_X4[gg];
          gauss_pts(1, gg) = G_TET_Y4[gg];
          gauss_pts(2, gg) = G_TET_Z4[gg];
          gauss_pts(3, gg) = G_TET_W4[gg];
        }
        CHKERR TetPolynomialBase().getValue(
            gauss_pts,
            boost::shared_ptr<BaseFunctionCtx>(
                new EntPolynomialBaseCtx(data, H1, AINSWORTH_LEGENDRE_BASE)));

        EntityHandle ent = numeredEntFiniteElementPtr->getEnt();
        int num_nodes;
        const EntityHandle *conn;
        CHKERR mField.get_moab().get_connectivity(ent, conn, num_nodes, true);
        coords.resize(num_nodes * 3);
        CHKERR mField.get_moab().get_coords(conn, num_nodes,
                                            &*coords.data().begin());

        CHKERR ShapeJacMBTET(
            &*data.dataOnEntities[MBVERTEX][0].getDiffN().data().begin(),
            &*coords.begin(), &*invJac.data().begin());
        CHKERR ShapeInvJacVolume(&*invJac.data().begin());

        CHKERR opSetInvJac.opRhs(data);
        CHKERR opPrintJac.opRhs(data);
        CHKERR opGetData_FIELD1.opRhs(data);

        my_split << "data FIELD1:" << std::endl;
        my_split << dataFIELD1 << std::endl;
        my_split << "data diff FIELD1:" << std::endl;
        my_split << dataDiffFIELD1 << std::endl;

        MoFEMFunctionReturn(0);
      }

      MoFEMErrorCode postProcess() {
        MoFEMFunctionBeginHot;
        MoFEMFunctionReturnHot(0);
      }
    };

    ForcesAndSourcesCore_TestFE fe1(m_field);
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TEST_FE", fe1);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
