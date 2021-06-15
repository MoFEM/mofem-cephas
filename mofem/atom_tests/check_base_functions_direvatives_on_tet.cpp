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

static const double eps = 1e-6;
static const double eps_diff = 1e-8;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Select space
    enum spaces { H1TET, HDIVTET, HCURLTET, LASTSPACEOP };

    const char *list_spaces[] = {"h1tet", "hdivtet", "hcurltet"};

    PetscBool flg;
    PetscInt choise_space_value = H1TET;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                LASTSPACEOP, &choise_space_value, &flg);
    
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "space not set");
    }

    FieldSpace space = LASTSPACE;
    if (choise_space_value == H1TET) {
      space = H1;
    } else if (choise_space_value == HDIVTET) {
      space = HDIV;
    } else if (choise_space_value == HCURLTET) {
      space = HCURL;
    }

    // Select base
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };

    const char *list_bases[] = {"ainsworth", "demkowicz"};

    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);
    
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    FieldApproximationBase base = NOBASE;
    if (choice_base_value == AINSWORTH) {
      base = AINSWORTH_LEGENDRE_BASE;
    } else if (choice_base_value == DEMKOWICZ) {
      base = DEMKOWICZ_JACOBI_BASE;
    }

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // create one tet
    double tet_coords[] = {0, 0, 0, 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5};
    EntityHandle nodes[4];
    for (int nn = 0; nn < 4; nn++) {
      CHKERR moab.create_vertex(&tet_coords[3 * nn], nodes[nn]);
    }
    EntityHandle tet;
    CHKERR moab.create_element(MBTET, nodes, 4, tet);

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    

    // Fields
    CHKERR m_field.add_field("FIELD", space, base, 1);

    // FE TET
    CHKERR m_field.add_finite_element("TET_FE");
    
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TET_FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_col("TET_FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_data("TET_FE", "FIELD");
    

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");
    

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TET_FE");
    
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);
    

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "FIELD");
    
    // add entities to finite element
    CHKERR
        m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "TET_FE");
    

    // set app. order
    int order = 5;
    if (space == H1) {
      CHKERR m_field.set_field_order(root_set, MBTET, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD", 1);
      
    }
    if (space == HCURL) {
      CHKERR m_field.set_field_order(root_set, MBTET, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD", order);
      
    }
    if (space == HDIV) {
      CHKERR m_field.set_field_order(root_set, MBTET, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
    }

    /****/
    // build database
    // build field
    CHKERR m_field.build_fields();
    
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);
    

    /****/
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    
    // build problem
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    
    // mesh partitioning
    // partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    
    // what are ghost nodes, see Petsc Manual
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;

    std::ofstream ofs("forces_and_sources_checking_direvatives.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    struct OpCheckingDirevatives
        : public VolumeElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpCheckingDirevatives(TeeStream &my_split)
          : VolumeElementForcesAndSourcesCore::UserDataOperator(
                "FIELD", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBegin;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);

        mySplit << "Type  " << type << " side " << side << endl;

        if (data.getFieldDofs()[0]->getSpace() == H1) {

          // mySplit << std::fixed << data.getN() << std::endl;
          // mySplit << std::fixed << data.getDiffN() << std::endl;

          const int nb_dofs = data.getN().size2();
          for (int dd = 0; dd != nb_dofs; dd++) {
            const double dksi = (data.getN()(1, dd) - data.getN()(0, dd)) / eps;
            const double deta = (data.getN()(3, dd) - data.getN()(2, dd)) / eps;
            const double dzeta =
                (data.getN()(5, dd) - data.getN()(4, dd)) / eps;
            mySplit << "DKsi " << dksi - data.getDiffN()(6, 3 * dd + 0)
                    << std::endl;
            mySplit << "DEta " << deta - data.getDiffN()(6, 3 * dd + 1)
                    << std::endl;
            mySplit << "DZeta " << dzeta - data.getDiffN()(6, 3 * dd + 2)
                    << std::endl;
            if (fabs(dksi - data.getDiffN()(6, 3 * dd + 0)) > eps_diff) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "H1 inconsistent dKsi derivative");
            }
            if (fabs(deta - data.getDiffN()(6, 3 * dd + 1)) > eps_diff) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "H1 inconsistent dEta derivative");
            }
            if (fabs(dzeta - data.getDiffN()(6, 3 * dd + 2)) > eps_diff) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "H1 inconsistent dZeta derivative");
            }
          }
        }

        if (data.getFieldDofs()[0]->getSpace() == HDIV ||
            data.getFieldDofs()[0]->getSpace() == HCURL) {

          FTensor::Tensor1<double *, 3> base_ksi_m(
              &data.getN()(0, 0), &data.getN()(0, 1), &data.getN()(0, 2), 3);
          FTensor::Tensor1<double *, 3> base_ksi_p(
              &data.getN()(1, 0), &data.getN()(1, 1), &data.getN()(1, 2), 3);
          FTensor::Tensor1<double *, 3> base_eta_m(
              &data.getN()(2, 0), &data.getN()(2, 1), &data.getN()(2, 2), 3);
          FTensor::Tensor1<double *, 3> base_eta_p(
              &data.getN()(3, 0), &data.getN()(3, 1), &data.getN()(3, 2), 3);
          FTensor::Tensor1<double *, 3> base_zeta_m(
              &data.getN()(4, 0), &data.getN()(4, 1), &data.getN()(4, 2), 3);
          FTensor::Tensor1<double *, 3> base_zeta_p(
              &data.getN()(5, 0), &data.getN()(5, 1), &data.getN()(5, 2), 3);

          FTensor::Tensor2<double *, 3, 3> diff_base(
              &data.getDiffN()(6, 0), &data.getDiffN()(6, 3),
              &data.getDiffN()(6, 6), &data.getDiffN()(6, 1),
              &data.getDiffN()(6, 4), &data.getDiffN()(6, 7),
              &data.getDiffN()(6, 2), &data.getDiffN()(6, 5),
              &data.getDiffN()(6, 8), 9);

          FTensor::Index<'i', 3> i;
          FTensor::Number<0> N0;
          FTensor::Number<1> N1;
          FTensor::Number<2> N2;

          const int nb_dofs = data.getN().size2() / 3;
          for (int dd = 0; dd != nb_dofs; dd++) {
            FTensor::Tensor1<double, 3> dksi;
            dksi(i) = (base_ksi_p(i) - base_ksi_m(i)) / eps;
            FTensor::Tensor1<double, 3> deta;
            deta(i) = (base_eta_p(i) - base_eta_m(i)) / eps;
            FTensor::Tensor1<double, 3> dzeta;
            dzeta(i) = (base_zeta_p(i) - base_zeta_m(i)) / eps;

            dksi(i) -= diff_base(i, N0);
            deta(i) -= diff_base(i, N1);
            dzeta(i) -= diff_base(i, N2);

            mySplit << "dKsi " << dksi(0) << "  " << dksi(1) << " " << dksi(2)
                    << " " << sqrt(dksi(i) * dksi(i)) << endl;
            mySplit << "dEta " << deta(0) << "  " << deta(1) << " " << deta(2)
                    << " " << sqrt(deta(i) * deta(i)) << endl;
            mySplit << "dZeta " << dzeta(0) << "  " << dzeta(1) << " "
                    << dzeta(2) << " " << sqrt(dzeta(i) * dzeta(i)) << endl;

            if (sqrt(dksi(i) * dksi(i)) > eps_diff) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "%s inconsistent dKsi derivative  for type %d",
                       FieldSpaceNames[data.getFieldDofs()[0]->getSpace()],
                       type);
            }
            if (sqrt(deta(i) * deta(i)) > eps_diff) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "%s inconsistent dEta derivative for type %d",
                       FieldSpaceNames[data.getFieldDofs()[0]->getSpace()],
                       type);
            }
            if (sqrt(dzeta(i) * dzeta(i)) > eps_diff) {
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "%s inconsistent dZeta derivative for type %d",
                       FieldSpaceNames[data.getFieldDofs()[0]->getSpace()],
                       type);
            }

            ++base_ksi_m;
            ++base_ksi_p;
            ++base_eta_m;
            ++base_eta_p;
            ++base_zeta_m;
            ++base_zeta_p;
            ++diff_base;
          }
        }

        MoFEMFunctionReturn(0);
      }
    };

    struct MyFE : public VolumeElementForcesAndSourcesCore {

      MyFE(MoFEM::Interface &m_field)
          : VolumeElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return -1; };

      MoFEMErrorCode setGaussPts(int order) {
        MoFEMFunctionBeginHot;

        const double ksi = G_TET_X1[0];
        const double eta = G_TET_Y1[0];
        const double zeta = G_TET_Z1[0];

        gaussPts.resize(4, 7);
        gaussPts.clear();

        gaussPts(0, 0) = ksi - eps;
        gaussPts(1, 0) = eta;
        gaussPts(2, 0) = zeta;
        gaussPts(0, 1) = ksi + eps;
        gaussPts(1, 1) = eta;
        gaussPts(2, 1) = zeta;

        gaussPts(0, 2) = ksi;
        gaussPts(1, 2) = eta - eps;
        gaussPts(2, 2) = zeta;
        gaussPts(0, 3) = ksi;
        gaussPts(1, 3) = eta + eps;
        gaussPts(2, 3) = zeta;

        gaussPts(0, 4) = ksi;
        gaussPts(1, 4) = eta;
        gaussPts(2, 4) = zeta - eps;
        gaussPts(0, 5) = ksi;
        gaussPts(1, 5) = eta;
        gaussPts(2, 5) = zeta + eps;

        gaussPts(0, 6) = ksi;
        gaussPts(1, 6) = eta;
        gaussPts(2, 6) = zeta;

        for (unsigned int ii = 0; ii != gaussPts.size2(); ii++) {
          gaussPts(3, ii) = 1;
        }

        cerr << gaussPts << endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    MyFE tet_fe(m_field);

    tet_fe.getOpPtrVector().push_back(new OpCheckingDirevatives(my_split));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TET_FE", tet_fe);
    
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  
}
