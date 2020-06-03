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
static const double eps_diff = 1e-4;

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    enum spaces { H1TRI, HCURLTRI, LASTOP };

    const char *list_spaces[] = {"h1tri", "hcurltri"};

    PetscBool flg;
    PetscInt choice_value = H1TRI;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces, LASTOP,
                                &choice_value, &flg);
    
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    }

    FieldSpace space = LASTSPACE;
    if (choice_value == H1TRI) {
      space = H1;
    } else if (choice_value == HCURLTRI) {
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
    double tri_coords[] = {0,   0,   0,

                           0.5, 0,   0,

                           0,   2., 0};
    EntityHandle nodes[3];
    for (int nn = 0; nn < 3; nn++) {
      CHKERR moab.create_vertex(&tri_coords[3 * nn], nodes[nn]);
      
    }
    EntityHandle tri;
    CHKERR moab.create_element(MBTRI, nodes, 3, tri);
    
    // Create adjacencies entities
    Range adj;
    CHKERR moab.get_adjacencies(&tri, 1, 1, true, adj);
    

    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // Create MoFEM database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    BitRefLevel bit_level0;
    bit_level0.set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, bit_level0);
    

    // Fields
    CHKERR m_field.add_field("FIELD", space, base, 1);
    

    // FE TET
    CHKERR m_field.add_finite_element("TRI_FE");
    
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("TRI_FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_col("TRI_FE", "FIELD");
    CHKERR m_field.modify_finite_element_add_field_data("TRI_FE", "FIELD");
    

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");

    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "TRI_FE");
    
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);
    

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTRI, "FIELD");
    
    // add entities to finite element
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTRI,
                                                      "TRI_FE");

    // set app. order
    int order = 5;
    if (space == H1) {
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBVERTEX, "FIELD", 1);
      
    }
    if (space == HCURL) {
      CHKERR m_field.set_field_order(root_set, MBTRI, "FIELD", order);
      CHKERR m_field.set_field_order(root_set, MBEDGE, "FIELD", order);
    }

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
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      TeeStream &mySplit;
      OpCheckingDirevatives(TeeStream &my_split)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "FIELD", UserDataOperator::OPROW),
            mySplit(my_split) {}

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
        MoFEMFunctionBegin;

        if (data.getFieldData().size() == 0)
          MoFEMFunctionReturnHot(0);
        mySplit << "type " << type << " side " << side << endl;

        if (data.getFieldDofs()[0].lock()->getSpace() == H1) {

          // mySplit << std::fixed << data.getN() << std::endl;
          // mySplit << std::fixed << data.getDiffN() << std::endl;

          const int nb_dofs = data.getN().size2();
          for (int dd = 0; dd != nb_dofs; dd++) {
            const double dksi =
                (data.getN()(1, dd) - data.getN()(0, dd)) / (eps);
            const double deta =
                (data.getN()(3, dd) - data.getN()(2, dd)) / (4 * eps);
            mySplit << "DKsi " << dksi << std::endl;
            mySplit << "DEta " << deta << std::endl;
            mySplit << "diffN " << data.getDiffN()(4, 2 * dd + 0) << std::endl;
            mySplit << "diffN " << data.getDiffN()(4, 2 * dd + 1) << std::endl;
            if (fabs(dksi - data.getDiffN()(4, 2 * dd + 0)) > eps_diff) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "H1 inconsistent dKsi derivative");
            }
            if (fabs(deta - data.getDiffN()(4, 2 * dd + 1)) > eps_diff) {
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "H1 inconsistent dEta derivative");
            }
          }
        }

        if (data.getFieldDofs()[0].lock()->getSpace() == HCURL) {

          FTensor::Tensor1<double *, 3> base_ksi_m(&data.getN()(0, HVEC0),
                                                   &data.getN()(0, HVEC1),
                                                   &data.getN()(0, HVEC2), 3);
          FTensor::Tensor1<double *, 3> base_ksi_p(&data.getN()(1, HVEC0),
                                                   &data.getN()(1, HVEC1),
                                                   &data.getN()(1, HVEC2), 3);
          FTensor::Tensor1<double *, 3> base_eta_m(&data.getN()(2, HVEC0),
                                                   &data.getN()(2, HVEC1),
                                                   &data.getN()(2, HVEC2), 3);
          FTensor::Tensor1<double *, 3> base_eta_p(&data.getN()(3, HVEC0),
                                                   &data.getN()(3, HVEC1),
                                                   &data.getN()(3, HVEC2), 3);

          // cerr << data.getN() << endl;
          // cerr << data.getDiffN() << endl;

          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> diff_base(
              &data.getDiffN()(4, HVEC0_0), &data.getDiffN()(4, HVEC0_1),
              &data.getDiffN()(4, HVEC1_0), &data.getDiffN()(4, HVEC1_1),
              &data.getDiffN()(4, HVEC2_0), &data.getDiffN()(4, HVEC2_1));

          FTensor::Index<'i', 3> i;
          FTensor::Number<0> N0;
          FTensor::Number<1> N1;

          const int nb_dofs = data.getN().size2() / 3;
          for (int dd = 0; dd != nb_dofs; dd++) {

            mySplit << "MoFEM " << diff_base(0, 0) << " " << diff_base(1, 0)
                    << " " << diff_base(2, 0) << endl;
            mySplit << "MoFEM " << diff_base(0, 1) << " " << diff_base(1, 1)
                    << " " << diff_base(2, 1) << endl;

            FTensor::Tensor1<double, 3> dksi;
            dksi(i) = (base_ksi_p(i) - base_ksi_m(i)) / (eps);
            FTensor::Tensor1<double, 3> deta;
            deta(i) = (base_eta_p(i) - base_eta_m(i)) / (4 * eps);
            mySplit << "Finite difference dKsi " << dksi(0) << "  " << dksi(1)
                    << " " << dksi(2) << endl;
            mySplit << "Finite difference dEta " << deta(0) << "  " << deta(1)
                    << " " << deta(2) << endl;

            dksi(i) -= diff_base(i, N0);
            deta(i) -= diff_base(i, N1);
            if (sqrt(dksi(i) * dksi(i)) > eps_diff) {
              // mySplit << "KSI ERROR\n";
              SETERRQ2(
                  PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "%s inconsistent dKsi derivative  for type %d",
                  FieldSpaceNames[data.getFieldDofs()[0].lock()->getSpace()],
                  type);
            } else {
              mySplit << "OK" << std::endl;
            }
            if (sqrt(deta(i) * deta(i)) > eps_diff) {
              // mySplit << "ETA ERROR\n";
              SETERRQ2(
                  PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "%s inconsistent dEta derivative for type %d",
                  FieldSpaceNames[data.getFieldDofs()[0].lock()->getSpace()],
                  type);
            } else {
              mySplit << "OK" << std::endl;
            }

            ++base_ksi_m;
            ++base_ksi_p;
            ++base_eta_m;
            ++base_eta_p;
            ++diff_base;
          }
        }

        mySplit << endl;

        MoFEMFunctionReturn(0);
      }
    };

    struct MyFE : public FaceElementForcesAndSourcesCore {

      MyFE(MoFEM::Interface &m_field)
          : FaceElementForcesAndSourcesCore(m_field) {}
      int getRule(int order) { return -1; };

      MoFEMErrorCode setGaussPts(int order) {
        MoFEMFunctionBeginHot;

        const double ksi = G_TRI_X1[0];
        const double eta = G_TRI_Y1[0];

        gaussPts.resize(3, 5);
        gaussPts.clear();

        gaussPts(0, 0) = ksi - eps;
        gaussPts(1, 0) = eta;
        gaussPts(0, 1) = ksi + eps;
        gaussPts(1, 1) = eta;

        gaussPts(0, 2) = ksi;
        gaussPts(1, 2) = eta - eps;
        gaussPts(0, 3) = ksi;
        gaussPts(1, 3) = eta + eps;

        gaussPts(0, 4) = ksi;
        gaussPts(1, 4) = eta;

        for (unsigned int ii = 0; ii != gaussPts.size2(); ii++) {
          gaussPts(2, ii) = 1;
        }

        // cerr << gaussPts << endl;

        MoFEMFunctionReturnHot(0);
      }
    };

    MyFE tri_fe(m_field);

    MatrixDouble inv_jac;
    tri_fe.getOpPtrVector().push_back(new OpCalculateInvJacForFace(inv_jac));
    if (space == H1) {
      tri_fe.getOpPtrVector().push_back(new OpSetInvJacH1ForFace(inv_jac));
    }
    if (space == HCURL) {
      tri_fe.getOpPtrVector().push_back(new OpSetInvJacHcurlFace(inv_jac));
    }
    tri_fe.getOpPtrVector().push_back(new OpCheckingDirevatives(my_split));
    CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "TRI_FE", tri_fe);
    
    cerr << inv_jac << endl;
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
  
}
