/** \file forces_and_sources_getting_mult_H1_H1_atom_test.cpp
  \brief Atom test verifying forces and sources operator on H1 approx. space

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

constexpr double eps = 1e-6;
template <typename T> void zero_entries(T &t) {
  for (auto &v : t)
    if (std::abs(v) < eps)
      v = 0;
}

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

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const char *option;
    option = ""; ;
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // set entitities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    CHKERR moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERRG(rval);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);

    // Fields
    CHKERR m_field.add_field("FIELD1", H1, AINSWORTH_LEGENDRE_BASE, 1);
    CHKERR m_field.add_field("FIELD2", H1, AINSWORTH_LEGENDRE_BASE, 3);

    // FE
    CHKERR m_field.add_finite_element("TEST_FE");

    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("TEST_FE", "FIELD2");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("TEST_FE", "FIELD2");

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
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD2");
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
    CHKERR m_field.set_field_order(root_set, MBTET, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD2", order);
    CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD2", 1);

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    // build problem
    //
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    // const Problem_multiIndex *problems_ptr;
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

      struct my_mult_H1_H1 : public DataOperator {

        std::ofstream ofs;
        TeeDevice my_tee;
        TeeStream my_split;

        my_mult_H1_H1()
            : ofs("forces_and_sources_getting_mult_H1_H1_atom_test.txt"),
              my_tee(std::cout, ofs), my_split(my_tee){};

        ~my_mult_H1_H1() { my_split.close(); }

        ublas::matrix<FieldData> NN;

        MoFEMErrorCode doWork(int row_side, int col_side, EntityType row_type,
                              EntityType col_type,
                              DataForcesAndSourcesCore::EntData &row_data,
                              DataForcesAndSourcesCore::EntData &col_data) {
          MoFEMFunctionBegin;

          row_data.getBase() = AINSWORTH_LEGENDRE_BASE;
          col_data.getBase() = AINSWORTH_LEGENDRE_BASE;
          int nb_row_dofs = row_data.getN().size2();
          int nb_col_dofs = col_data.getN().size2();

          my_split << row_side << " " << col_side << " " << row_type << " "
                   << col_type << std::endl;
          my_split << "nb_row_dofs " << nb_row_dofs << " nb_col_dofs "
                   << nb_col_dofs << std::endl;
          NN.resize(nb_row_dofs, nb_col_dofs);

          zero_entries(row_data.getN().data());
          zero_entries(row_data.getDiffN().data());

          my_split << std::setprecision(3);
          my_split << std::fixed;
          my_split << row_data.getN() << std::endl;
          my_split << col_data.getN() << std::endl;

          for (unsigned int gg = 0; gg < row_data.getN().size1(); gg++) {

            bzero(&*NN.data().begin(),
                  nb_row_dofs * nb_col_dofs * sizeof(FieldData));

            cblas_dger(CblasRowMajor, nb_row_dofs, nb_col_dofs, 1,
                       &row_data.getN()(gg, 0), 1, &col_data.getN()(gg, 0), 1,
                       &*NN.data().begin(), nb_col_dofs);

            my_split << "gg " << gg << " : ";
            my_split << std::setprecision(3);
            my_split << std::fixed;

            MatrixDouble difference =
                NN - outer_prod(row_data.getN(gg), col_data.getN(gg));
            zero_entries(difference.data());

            my_split << difference << std::endl;
            if (row_type != MBVERTEX) {
              my_split << row_data.getDiffN(gg) << std::endl;
            }

            if (row_type == MBVERTEX) {
              my_split << row_data.getDiffN() << std::endl;
            } else {
              typedef ublas::array_adaptor<FieldData> storage_t;
              storage_t st(nb_row_dofs * 3, &row_data.getDiffN()(gg, 0));
              ublas::matrix<FieldData, ublas::row_major, storage_t>
                  digNatGaussPt(nb_row_dofs, 3, st);
              my_split << std::endl << digNatGaussPt << std::endl;
            }
          }

          my_split << std::endl;

          MoFEMFunctionReturn(0);
        }
      };

      my_mult_H1_H1 op;

      ForcesAndSourcesCore_TestFE(MoFEM::Interface &_m_field)
          : ForcesAndSourcesCore(_m_field), data_row(MBTET), data_col(MBTET){};

      MoFEMErrorCode preProcess() {
        MoFEMFunctionBeginHot;
        MoFEMFunctionReturnHot(0);
      }

      DataForcesAndSourcesCore data_row, data_col;

      MoFEMErrorCode operator()() {
        MoFEMFunctionBeginHot;

        CHKERR getSpacesAndBaseOnEntities(data_row);
        CHKERR getSpacesAndBaseOnEntities(data_col);

        CHKERR getEntitySense<MBEDGE>(data_row);
        CHKERR getEntitySense<MBTRI>(data_row);
        CHKERR getEntitySense<MBEDGE>(data_col);
        CHKERR getEntitySense<MBTRI>(data_col);

        CHKERR getEntityDataOrder<MBEDGE>(data_row, H1);
        CHKERR getEntityDataOrder<MBEDGE>(data_col, H1);
        CHKERR getEntityDataOrder<MBTRI>(data_row, H1);
        CHKERR getEntityDataOrder<MBTRI>(data_col, H1);
        CHKERR getEntityDataOrder<MBTET>(data_row, H1);
        CHKERR getEntityDataOrder<MBTET>(data_col, H1);
        data_row.dataOnEntities[MBVERTEX][0].getBase() =
            AINSWORTH_LEGENDRE_BASE;
        CHKERR getEntityFieldData(data_row, "FIELD1", MBEDGE);
        data_col.dataOnEntities[MBVERTEX][0].getBase() =
            AINSWORTH_LEGENDRE_BASE;
        CHKERR getEntityFieldData(data_col, "FIELD2", MBEDGE);    
        CHKERR getRowNodesIndices(data_row, "FIELD1");
        CHKERR getColNodesIndices(data_col, "FIELD2");
        CHKERR getEntityRowIndices(data_row, "FIELD1", MBEDGE);
        CHKERR getEntityColIndices(data_col, "FIELD2", MBEDGE);
        CHKERR getFaceTriNodes(data_row);
        CHKERR getFaceTriNodes(data_col);

        MatrixDouble gauss_pts(4, 4);
        for (int gg = 0; gg < 4; gg++) {
          gauss_pts(0, gg) = G_TET_X4[gg];
          gauss_pts(1, gg) = G_TET_Y4[gg];
          gauss_pts(2, gg) = G_TET_Z4[gg];
          gauss_pts(3, gg) = G_TET_W4[gg];
        }
        CHKERR TetPolynomialBase().getValue(
            gauss_pts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                data_row, H1, AINSWORTH_LEGENDRE_BASE)));
        CHKERR TetPolynomialBase().getValue(
            gauss_pts,
            boost::shared_ptr<BaseFunctionCtx>(new EntPolynomialBaseCtx(
                data_col, H1, AINSWORTH_LEGENDRE_BASE)));

        CHKERR op.opLhs(data_row, data_col);

        MoFEMFunctionReturnHot(0);
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
