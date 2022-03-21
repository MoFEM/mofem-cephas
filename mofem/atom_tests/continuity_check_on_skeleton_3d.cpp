/** \file continuity_check_on_skeleton_3d.cpp
 * \example continuity_check_on_skeleton_3d.cpp
 *
 * \brief Testing integration on skeleton for 3D
 *
 * Checking continuity of hdiv and hcurl spaces on faces, and testing methods
 * for integration on the skeleton.
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

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  auto core_log = logging::core::get();
  core_log->add_sink(
      LogManager::createSink(LogManager::getStrmSelf(), "ATOM_TEST"));
  LogManager::setLog("ATOM_TEST");

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
    auto get_base = []() -> FieldApproximationBase {
      enum bases { AINSWORTH, DEMKOWICZ, LASTBASEOP };
      const char *list_bases[] = {"ainsworth", "demkowicz"};
      PetscBool flg;
      PetscInt choice_base_value = AINSWORTH;
      CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-base", list_bases,
                                  LASTBASEOP, &choice_base_value, &flg);
      if (flg == PETSC_TRUE) {
        FieldApproximationBase base = AINSWORTH_LEGENDRE_BASE;
        if (choice_base_value == AINSWORTH)
          base = AINSWORTH_LEGENDRE_BASE;
        else if (choice_base_value == DEMKOWICZ)
          base = DEMKOWICZ_JACOBI_BASE;
        return base;
      }
      return LASTBASE;
    };

    auto get_space = []() -> FieldSpace {
      enum spaces { HDIV, HCURL, LASTSPACEOP };
      const char *list_spaces[] = {"hdiv", "hcurl"};
      PetscBool flg;
      PetscInt choice_space_value = HDIV;
      CHKERR PetscOptionsGetEList(PETSC_NULL, NULL, "-space", list_spaces,
                                  LASTSPACEOP, &choice_space_value, &flg);
      if (flg == PETSC_TRUE) {
        FieldSpace space = FieldSpace::HDIV;
        if (choice_space_value == HDIV)
          space = FieldSpace::HDIV;
        else if (choice_space_value == HCURL)
          space = FieldSpace::HCURL;
        return space;
      }
      return LASTSPACE;
    };

    CHKERR m_field.add_field("F2", get_space(), get_base(), 1);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_dim(root_set, 3, "F2");

    // set app. order
    int order = 3;
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBHEX, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBQUAD, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "F2", order);

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

    CHKERR m_field.add_ents_to_finite_element_by_dim(root_set, 3, "V1");
    Range faces;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByDimAndRefLevel(
        bit_level0, BitRefLevel().set(), 2, faces);
    CHKERR m_field.add_ents_to_finite_element_by_dim(faces, 2, "S2");

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

    struct CommonData {
      VectorDouble dotNormalFace;
      VectorDouble dotNormalEleLeft;
      VectorDouble dotNormalEleRight;
    };

    struct SkeletonFE
        : public FaceElementForcesAndSourcesCore::UserDataOperator {

      VolumeElementForcesAndSourcesCoreOnSide volSideFe;
      struct OpVolSide
          : public VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator {

        CommonData &elemData;
        OpVolSide(CommonData &elem_data)
            : VolumeElementForcesAndSourcesCoreOnSide::UserDataOperator(
                  "F2", UserDataOperator::OPROW),
              elemData(elem_data) {}
        MoFEMErrorCode doWork(int side, EntityType type,
                              EntitiesFieldData::EntData &data) {
          MoFEMFunctionBeginHot;

          if (CN::Dimension(type) == 3) {
            MatrixDouble diff =
                getCoordsAtGaussPts() - getFaceCoordsAtGaussPts();

            MOFEM_LOG("ATOM_TEST", Sev::noisy)
                << "getCoordsAtGaussPts() " << getCoordsAtGaussPts();
            MOFEM_LOG("ATOM_TEST", Sev::noisy)
                << "getFaceCoordsAtGaussPts() " << getFaceCoordsAtGaussPts();

            constexpr double eps = 1e-12;
            if (norm_inf(diff) > eps) {
              MOFEM_LOG("ATOM_TEST", Sev::error) << "diff " << diff;
              SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                      "Coordinates at integration pts are different");
            }
          }

          const size_t nb_dofs = data.getFieldData().size();
          if (nb_dofs) {
            const size_t nb_integration_pts = getGaussPts().size2();

            FTensor::Index<'i', 3> i;
            auto t_to_do_dot = getFTensor1Normal();
            if (data.getSpace() == HCURL) {
              auto s1 = getFTensor1Tangent1();
              auto s2 = getFTensor1Tangent2();
              t_to_do_dot(i) = s1(i) + s2(i);
            }

            VectorDouble *ptr_dot_elem_data = nullptr;
            if (getFEMethod()->nInTheLoop == 0) 
              ptr_dot_elem_data = &elemData.dotNormalEleLeft;
            else
              ptr_dot_elem_data = &elemData.dotNormalEleRight;
            auto &dot_elem_data = *ptr_dot_elem_data;
            if (dot_elem_data.size() == 0) {
              dot_elem_data.resize(nb_integration_pts, false);
              dot_elem_data.clear();
            }

            auto t_hdiv_base = data.getFTensor1N<3>();
            for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
              auto t_data = data.getFTensor0FieldData();
              for (size_t bb = 0; bb != nb_dofs; ++bb) {
                dot_elem_data(gg) += (t_to_do_dot(i) * t_hdiv_base(i)) * t_data;
                ++t_hdiv_base;
                ++t_data;
              }
            }
          }
          MoFEMFunctionReturnHot(0);
        }
      };

      CommonData &elemData;
      SkeletonFE(MoFEM::Interface &m_field, CommonData &elem_data)
          : FaceElementForcesAndSourcesCore::UserDataOperator(
                "F2", UserDataOperator::OPROW),
            volSideFe(m_field), elemData(elem_data) {

        auto jac_ptr = boost::make_shared<MatrixDouble>();
        auto inv_jac_ptr = boost::make_shared<MatrixDouble>();
        auto det_ptr = boost::make_shared<VectorDouble>();

        volSideFe.getOpPtrVector().push_back(new OpCalculateHOJacVolume(jac_ptr));
        volSideFe.getOpPtrVector().push_back(
            new OpInvertMatrix<3>(jac_ptr, det_ptr, inv_jac_ptr));
        volSideFe.getOpPtrVector().push_back(new OpSetHOWeights(det_ptr));
        volSideFe.getOpPtrVector().push_back(
            new OpSetHOCovariantPiolaTransform(HCURL, inv_jac_ptr));
        volSideFe.getOpPtrVector().push_back(
            new OpSetHOContravariantPiolaTransform(HDIV, det_ptr, jac_ptr));
        volSideFe.getOpPtrVector().push_back(
            new SkeletonFE::OpVolSide(elemData));
      }

      MoFEMErrorCode doWork(int side, EntityType type,
                            EntitiesFieldData::EntData &data) {

        MoFEMFunctionBegin;

        const size_t nb_integration_pts = getGaussPts().size2();

        if (side == 0 && type == MBEDGE) {
          elemData.dotNormalFace.resize(nb_integration_pts, false);
          elemData.dotNormalFace.clear();
        }

        const size_t nb_dofs = data.getFieldData().size();
        if (nb_dofs) {
          FTensor::Index<'i', 3> i;
          auto t_to_do_dot = getFTensor1Normal();
          if (data.getSpace() == HCURL) {
            auto s1 = getFTensor1Tangent1();
            auto s2 = getFTensor1Tangent2();
            t_to_do_dot(i) = s1(i) + s2(i);
          }
          auto t_hdiv_base = data.getFTensor1N<3>();
          for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
            auto t_data = data.getFTensor0FieldData();
            for (size_t bb = 0; bb != nb_dofs; ++bb) {
              elemData.dotNormalFace(gg) +=
                  (t_to_do_dot(i) * t_hdiv_base(i)) * t_data;
              ++t_hdiv_base;
              ++t_data;
            }
          }
        }

        if (CN::Dimension(type) == 2) {

          elemData.dotNormalEleLeft.resize(0, false);
          elemData.dotNormalEleRight.resize(0, false);
          CHKERR loopSideVolumes("V1", volSideFe);

          auto check_continuity_of_base = [&](auto &vol_dot_data) {
            MoFEMFunctionBegin;
            if (vol_dot_data.size() != elemData.dotNormalFace.size())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of integration points");
            const double eps = 1e-12;
            for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
              const double error =
                  std::abs(vol_dot_data(gg) - elemData.dotNormalFace(gg));
              MOFEM_LOG("ATOM_TEST", Sev::noisy) << "Error: " << error;
              if (error > eps)
                SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                         "Inconsistency %3.4e != %3.4e", vol_dot_data(gg),
                         elemData.dotNormalFace(gg));
            }
            MoFEMFunctionReturn(0);
          };
          if (elemData.dotNormalEleLeft.size() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleLeft);
          else if (elemData.dotNormalEleRight.size() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleRight);
        }

        MoFEMFunctionReturn(0);
      }
    };

    CommonData elem_data;
    FaceElementForcesAndSourcesCore face_fe(m_field);
    face_fe.getRuleHook = [](int, int, int) { return 1; };
    face_fe.getOpPtrVector().push_back(
        new OpHOSetContravariantPiolaTransformOnFace3D(HDIV));
    face_fe.getOpPtrVector().push_back(
        new OpHOSetCovariantPiolaTransformOnFace3D(HCURL));
    face_fe.getOpPtrVector().push_back(new SkeletonFE(m_field, elem_data));

    auto field_ents_ptr = m_field.get_field_ents();

    auto cache_ptr = boost::make_shared<CacheTuple>();
    CHKERR m_field.cache_problem_entities("P1", cache_ptr);

    for (auto &ent_ptr : (*field_ents_ptr)) {
      MOFEM_LOG("ATOM_TEST", Sev::verbose) << *ent_ptr;
      for(auto &v : ent_ptr->getEntFieldData())
        v = 1;
      CHKERR m_field.loop_finite_elements("P1", "S2", face_fe, 0, 1, nullptr,
                                          MF_ZERO, cache_ptr);
      for (auto &v : ent_ptr->getEntFieldData())
        v = 0;
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
