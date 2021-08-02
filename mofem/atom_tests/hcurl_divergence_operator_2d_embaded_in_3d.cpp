/**
 * \file hcurl_divergence_operator_2d.cpp
 * \example hcurl_divergence_operator_2d.cpp
 *
 * Testing Hcurl base, transfromed to Hdiv base in 2d using Green theorem
 * 
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

using FaceEle = MoFEM::FaceElementForcesAndSourcesCore;

using EdgeEle = MoFEM::EdgeElementForcesAndSourcesCore;

using FaceEleOp = FaceEle::UserDataOperator;
using EdgeEleOp = EdgeEle::UserDataOperator;

struct OpDivergence : public FaceEleOp {

  double &dIv;
  OpDivergence(double &div) : FaceEleOp("FIELD1", OPROW), dIv(div) {}

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();
    auto t_diff_base_fun = data.getFTensor2DiffN<3, 3>();
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getArea() * getGaussPts()(2, gg);
      for (int bb = 0; bb != nb_dofs; bb++) {
        dIv += val * (t_diff_base_fun(0, 0) + t_diff_base_fun(1, 1) +
                      t_diff_base_fun(2, 2));
        ++t_diff_base_fun;
      }
    }
    MoFEMFunctionReturn(0);
  }
};

struct OpFlux : public EdgeEleOp {

  double &fLux;
  OpFlux(double &flux) : EdgeEleOp("FIELD1", OPROW), fLux(flux) {}

  FTensor::Index<'i', 3> i;

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {
    MoFEMFunctionBegin;
    const int nb_dofs = data.getIndices().size();
    if (nb_dofs == 0)
      MoFEMFunctionReturnHot(0);
    const int nb_gauss_pts = data.getN().size1();

    auto t_base_fun = data.getFTensor1N<3>();
    FTensor::Index<'i', 2> i;
    for (int gg = 0; gg != nb_gauss_pts; gg++) {
      const double val = getGaussPts()(1, gg);
      for (int bb = 0; bb != nb_dofs; bb++) {
        fLux += val * t_base_fun(0);
        ++t_base_fun;
      }
    }

    MoFEMFunctionReturn(0);
  }
};

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    PetscBool flg_file = PETSC_TRUE;
    char mesh_file_name[255];
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg_file);
    if (flg_file != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA,
              "*** ERROR -my_file (MESH FILE NEEDED)");

    // Read mesh to MOAB
    CHKERR moab.load_file(mesh_file_name, 0, "");

    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // set entities bit level
    BitRefLevel bit_level0 = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 2, bit_level0);

    // Declare elements
    enum bases { AINSWORTH, DEMKOWICZ, LASBASETOP };
    const char *list_bases[] = {"ainsworth", "demkowicz"};
    PetscBool flg;
    PetscInt choice_base_value = AINSWORTH;
    CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                LASBASETOP, &choice_base_value, &flg);
    if (flg != PETSC_TRUE)
      SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "base not set");
    FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
    if (choice_base_value == AINSWORTH)
      base = AINSWORTH_LEGENDRE_BASE;
    else if (choice_base_value == DEMKOWICZ)
      base = DEMKOWICZ_JACOBI_BASE;
    int order = 5;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-order", &order, PETSC_NULL);

    CHKERR m_field.add_field("FIELD1", HCURL, base, 1);
    CHKERR m_field.add_finite_element("FACE_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("FACE_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("FACE_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("FACE_FE", "FIELD1");

    CHKERR m_field.add_finite_element("EDGE_FE");
    // Define rows/cols and element data
    CHKERR m_field.modify_finite_element_add_field_row("EDGE_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_col("EDGE_FE", "FIELD1");
    CHKERR m_field.modify_finite_element_add_field_data("EDGE_FE", "FIELD1");

    // Problem
    CHKERR m_field.add_problem("TEST_PROBLEM");
    // set finite elements for problem
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "FACE_FE");
    CHKERR m_field.modify_problem_add_finite_element("TEST_PROBLEM", "EDGE_FE");
    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("TEST_PROBLEM", bit_level0);

    // Add entities

    CHKERR m_field.add_ents_to_field_by_type(0, MBTRI, "FIELD1");
    CHKERR m_field.add_ents_to_field_by_type(0, MBQUAD, "FIELD1");
    // Set order
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD1", order);
    CHKERR m_field.set_field_order(0, MBQUAD, "FIELD1", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD1", order);

    // Add entities to elements
    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBTRI, "FACE_FE");
    CHKERR m_field.add_ents_to_finite_element_by_type(0, MBQUAD, "FACE_FE");

    auto set_edge_elements_entities_on_mesh_skin = [&]() {
      MoFEMFunctionBegin;
      Range faces;
      CHKERR moab.get_entities_by_dimension(0, 2, faces, false);
      Skinner skin(&m_field.get_moab());
      Range faces_skin;
      CHKERR skin.find_skin(0, faces, false, faces_skin);
      CHKERR m_field.add_ents_to_finite_element_by_type(faces_skin, MBEDGE,
                                                        "EDGE_FE");
      MoFEMFunctionReturn(0);
    };
    CHKERR set_edge_elements_entities_on_mesh_skin();

    // Build database
    CHKERR m_field.build_fields();
    // build finite elemnts
    CHKERR m_field.build_finite_elements();
    // build adjacencies
    CHKERR m_field.build_adjacencies(bit_level0);

    // build problem
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("TEST_PROBLEM", true);
    // Partition
    CHKERR prb_mng_ptr->partitionSimpleProblem("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionFiniteElements("TEST_PROBLEM");
    CHKERR prb_mng_ptr->partitionGhostDofs("TEST_PROBLEM");

    // integration rule
    auto rule = [&](int, int, int p) { return 2 * p; };

    auto calculate_divergence = [&]() {
      double div = 0;
      FaceEle fe_face(m_field);
      fe_face.getRuleHook = rule;
      MatrixDouble inv_jac(2, 2), jac(2, 2);

      fe_face.getOpPtrVector().push_back(
          new OpCalculateJacForFaceEmbeddedIn3DSpace(jac));
      fe_face.getOpPtrVector().push_back(
          new OpCalculateInvJacForFaceEmbeddedIn3DSpace(inv_jac));
      fe_face.getOpPtrVector().push_back(new OpMakeHdivFromHcurl());
      fe_face.getOpPtrVector().push_back(
          new OpSetContravariantPiolaTransformOnFace2DEmbeddedIn3DSpace(jac));
      fe_face.getOpPtrVector().push_back(
          new OpSetInvJacHcurlFaceEmbeddedIn3DSpace(inv_jac));
      fe_face.getOpPtrVector().push_back(
          new OpSetHOWeigthsOnFace());
      fe_face.getOpPtrVector().push_back(new OpDivergence(div));
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "FACE_FE", fe_face);
      return div;
    };

    auto calculate_flux = [&]() {
      double flux = 0;
      EdgeEle fe_edge(m_field);
      fe_edge.getRuleHook = rule;
      // fe_edge.getOpPtrVector().push_back(
      //     new OpSetContravariantPiolaTransformOnEdge2D());
      fe_edge.getOpPtrVector().push_back(new OpFlux(flux));
      CHKERR m_field.loop_finite_elements("TEST_PROBLEM", "EDGE_FE", fe_edge);
      return flux;
    };

    const double div = calculate_divergence();
    const double flux = calculate_flux();

    PetscPrintf(PETSC_COMM_WORLD, "Div = %4.3e Flux = %3.4e Error = %4.3e\n",
                div, flux, div + flux);

    constexpr double tol = 1e-8;
    if (std::abs(div + flux) > tol)
      SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
               "Test failed (div != flux) %3.4e != %3.4e", div, flux);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
