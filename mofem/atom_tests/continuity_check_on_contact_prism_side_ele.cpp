/** \file continuity_check_on_contact_prism_side_ele.cpp
 * \example continuity_check_on_contact_prism_side_ele.cpp
 *
 * \brief Testing integration on volume side on contact element
 *
 * Checking continuity of H1  ( later hdiv and hcurl spaces on faces), and
 * testing methods for integration on volume side on contact element.
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

    enum side_contact { MASTERSIDE, SLAVESIDE };

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
      enum spaces { H1, HDIV, HCURL, LASTSPACEOP };
      const char *list_spaces[] = {"h1", "hdiv", "hcurl"};
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
        else if (choice_space_value == H1)
          space = FieldSpace::H1;
        return space;
      }
      return LASTSPACE;
    };

    CHKERR m_field.add_field("F2", get_space(), get_base(), 1);

    auto add_prism_interface = [&](Range &tets, Range &prisms,
                                   Range &master_tris, Range &slave_tris,
                                   EntityHandle &meshset_tets,
                                   EntityHandle &meshset_prisms,
                                   std::vector<BitRefLevel> &bit_levels) {
      MoFEMFunctionBegin;
      PrismInterface *interface;
      CHKERR m_field.getInterface(interface);

      int ll = 1;
      for (_IT_CUBITMESHSETS_BY_SET_TYPE_FOR_LOOP_(m_field, SIDESET, cit)) {
        if (cit->getName().compare(0, 11, "INT_CONTACT") == 0) {
          CHKERR PetscPrintf(PETSC_COMM_WORLD, "Insert %s (id: %d)\n",
                             cit->getName().c_str(), cit->getMeshsetId());
          EntityHandle cubit_meshset = cit->getMeshset();
          Range tris;
          CHKERR moab.get_entities_by_type(cubit_meshset, MBTRI, tris, true);
          master_tris.merge(tris);

          {
            // get tet entities from back bit_level
            EntityHandle ref_level_meshset = 0;
            CHKERR moab.create_meshset(MESHSET_SET, ref_level_meshset);
            CHKERR m_field.getInterface<BitRefManager>()
                ->getEntitiesByTypeAndRefLevel(bit_levels.back(),
                                               BitRefLevel().set(), MBTET,
                                               ref_level_meshset);
            CHKERR m_field.getInterface<BitRefManager>()
                ->getEntitiesByTypeAndRefLevel(bit_levels.back(),
                                               BitRefLevel().set(), MBPRISM,
                                               ref_level_meshset);
            Range ref_level_tets;
            CHKERR moab.get_entities_by_handle(ref_level_meshset,
                                               ref_level_tets, true);
            // get faces and tets to split
            CHKERR interface->getSides(cubit_meshset, bit_levels.back(), true,
                                       0);
            // set new bit level
            bit_levels.push_back(BitRefLevel().set(ll++));
            // split faces and tets
            CHKERR interface->splitSides(ref_level_meshset, bit_levels.back(),
                                         cubit_meshset, true, true, 0);
            // clean meshsets
            CHKERR moab.delete_entities(&ref_level_meshset, 1);
          }
          // update cubit meshsets
          for (_IT_CUBITMESHSETS_FOR_LOOP_(m_field, ciit)) {
            EntityHandle cubit_meshset = ciit->meshset;
            CHKERR m_field.getInterface<BitRefManager>()
                ->updateMeshsetByEntitiesChildren(
                    cubit_meshset, bit_levels.back(), cubit_meshset, MBVERTEX,
                    true);
            CHKERR m_field.getInterface<BitRefManager>()
                ->updateMeshsetByEntitiesChildren(cubit_meshset,
                                                  bit_levels.back(),
                                                  cubit_meshset, MBEDGE, true);
            CHKERR m_field.getInterface<BitRefManager>()
                ->updateMeshsetByEntitiesChildren(cubit_meshset,
                                                  bit_levels.back(),
                                                  cubit_meshset, MBTRI, true);
            CHKERR m_field.getInterface<BitRefManager>()
                ->updateMeshsetByEntitiesChildren(cubit_meshset,
                                                  bit_levels.back(),
                                                  cubit_meshset, MBTET, true);
          }
        }
      }

      for (unsigned int ll = 0; ll != bit_levels.size() - 1; ll++) {
        CHKERR m_field.delete_ents_by_bit_ref(bit_levels[ll], bit_levels[ll],
                                              true);
      }
      CHKERR m_field.getInterface<BitRefManager>()->shiftRightBitRef(
          bit_levels.size() - 1);

      CHKERR moab.create_meshset(MESHSET_SET, meshset_tets);
      CHKERR moab.create_meshset(MESHSET_SET, meshset_prisms);

      CHKERR m_field.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(bit_levels[0], BitRefLevel().set(),
                                         MBTET, meshset_tets);
      CHKERR moab.get_entities_by_handle(meshset_tets, tets, true);

      CHKERR m_field.getInterface<BitRefManager>()
          ->getEntitiesByTypeAndRefLevel(bit_levels[0], BitRefLevel().set(),
                                         MBPRISM, meshset_prisms);
      CHKERR moab.get_entities_by_handle(meshset_prisms, prisms);

      Range tris;
      CHKERR moab.get_adjacencies(prisms, 2, false, tris,
                                  moab::Interface::UNION);
      tris = tris.subset_by_type(MBTRI);
      slave_tris = subtract(tris, master_tris);

      MoFEMFunctionReturn(0);
    };

    Range all_tets, contact_prisms, master_tris, slave_tris;
    EntityHandle meshset_tets, meshset_prisms;
    std::vector<BitRefLevel> bit_levels;

    bit_levels.push_back(BitRefLevel().set(0));
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_levels[0]);

    CHKERR add_prism_interface(all_tets, contact_prisms, master_tris,
                               slave_tris, meshset_tets, meshset_prisms,
                               bit_levels);

    // meshset consisting all entities in mesh
    EntityHandle root_set = moab.get_root_set();
    // add entities to field
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBTET, "F2");
    CHKERR m_field.add_ents_to_field_by_type(root_set, MBPRISM, "F2");

    // set app. order
    int order = 2;
    CHKERR m_field.set_field_order(root_set, MBTET, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBTRI, "F2", order);
    CHKERR m_field.set_field_order(root_set, MBEDGE, "F2", order);

    if (get_space() == FieldSpace::H1) {
      CHKERR m_field.set_field_order(root_set, MBVERTEX, "F2", 1);
    }

    CHKERR m_field.build_fields();

    // add elements
    CHKERR m_field.add_finite_element("V1");
    CHKERR m_field.add_finite_element("C2");
    CHKERR m_field.modify_finite_element_add_field_row("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_data("V1", "F2");
    CHKERR m_field.modify_finite_element_add_field_row("C2", "F2");
    CHKERR m_field.modify_finite_element_add_field_col("C2", "F2");
    CHKERR m_field.modify_finite_element_add_field_data("C2", "F2");

    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBTET, "V1");

    Range prism;
    CHKERR m_field.getInterface<BitRefManager>()->getEntitiesByTypeAndRefLevel(
        bit_levels[0], BitRefLevel().set(), MBPRISM, prism);
    CHKERR m_field.add_ents_to_finite_element_by_type(root_set, MBPRISM, "C2");

    CHKERR m_field.build_finite_elements();
    CHKERR m_field.build_adjacencies(bit_levels[0]);

    // Problems
    CHKERR m_field.add_problem("P1");

    // set refinement level for problem
    CHKERR m_field.modify_problem_ref_level_add_bit("P1", bit_levels[0]);
    CHKERR m_field.modify_problem_add_finite_element("P1", "V1");
    CHKERR m_field.modify_problem_add_finite_element("P1", "C2");

    // build problems
    ProblemsManager *prb_mng_ptr;
    CHKERR m_field.getInterface(prb_mng_ptr);
    CHKERR prb_mng_ptr->buildProblem("P1", true);

    struct CommonData {
      MatrixDouble dotNormalFace;
      MatrixDouble dotNormalEleLeft;
      MatrixDouble dotNormalEleRight;
      MatrixDouble shapeFunH1Values;
      MatrixDouble shapeFunH1VolSide;
    };

    struct OnContactSideMaster
        : public ContactPrismElementForcesAndSourcesCore::UserDataOperator {

      VolumeElementForcesAndSourcesCoreOnContactPrismSide volSideFe;
      struct OpVolSide : public VolumeElementForcesAndSourcesCoreOnContactPrismSide::
                             UserDataOperator {

        CommonData &elemData;
        OpVolSide(CommonData &elem_data)
            : VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator(
                  "F2", UserDataOperator::OPROW),
              elemData(elem_data) {}
        MoFEMErrorCode doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
          MoFEMFunctionBeginHot;

          if (type == MBTRI && side == getFaceSideNumber()) {

            MatrixDouble diff =
                getCoordsAtGaussPts() - getMasterCoordsAtGaussPts();
            constexpr double eps = 1e-12;
            if (norm_inf(diff) > eps)
              SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                      "coordinates at integration pts are different");

            const size_t nb_dofs = data.getN().size2() / 3;
            const size_t nb_integration_pts = data.getN().size1();

            if (data.getSpace() == H1) {
              MatrixDouble *ptr_elem_data = nullptr;
              ptr_elem_data = &elemData.shapeFunH1VolSide;

              MatrixDouble &elem_data = *ptr_elem_data;
              elem_data.resize(nb_integration_pts, nb_dofs, false);
              for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
                FTensor::Tensor0<double *> t_base(&data.getN()(gg, 0));
                for (size_t bb = 0; bb != nb_dofs; ++bb) {
                  elem_data(gg, bb) = t_base;
                  ++t_base;
                }
              }
            }
          }
          MoFEMFunctionReturnHot(0);
        }
      };

      CommonData &elemData;
      OnContactSideMaster(MoFEM::Interface &m_field, CommonData &elem_data)
          : ContactPrismElementForcesAndSourcesCore::UserDataOperator(
                "F2", UserDataOperator::OPROW,
                ContactPrismElementForcesAndSourcesCore::UserDataOperator::
                    FACEMASTER),
            volSideFe(m_field), elemData(elem_data) {
        volSideFe.getOpPtrVector().push_back(
            new OnContactSideMaster::OpVolSide(elemData));
      }

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {

        MoFEMFunctionBeginHot;
        if (type == MBTRI && side == 0) {
          const size_t nb_dofs = data.getN().size2() / 3;
          const size_t nb_integration_pts = data.getN().size1();
          if (data.getSpace() == H1) {
            elemData.shapeFunH1Values.resize(nb_integration_pts, nb_dofs,
                                             false);
            elemData.shapeFunH1VolSide.resize(nb_integration_pts, nb_dofs,
                                              false);

            for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
              FTensor::Tensor0<double *> t_base(&data.getN()(gg, 0));
              for (size_t bb = 0; bb != nb_dofs; ++bb) {
                elemData.shapeFunH1Values(gg, bb) = t_base;
                ++t_base;
              }
            }
          }
          std::string side_fe_name = "V1";
          const EntityHandle tri_master = getSideEntity(3, MBTRI);
          CHKERR loopSideVolumes(side_fe_name, volSideFe, 3, tri_master);

          auto check_continuity_of_base = [&](auto &vol_dot_data) {
            MoFEMFunctionBegin;

            if (vol_dot_data.size1() != elemData.dotNormalFace.size1())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of integration points");

            if (vol_dot_data.size2() != elemData.dotNormalFace.size2())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of base functions");
            constexpr double eps = 1e-12;
            for (size_t gg = 0; gg != vol_dot_data.size1(); ++gg)
              for (size_t bb = 0; bb != vol_dot_data.size2(); ++bb) {
                const double error = std::abs(vol_dot_data(gg, bb) -
                                              elemData.dotNormalFace(gg, bb));
                if (error > eps)
                  SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                           "Inconsistency %3.4e != %3.4e", vol_dot_data(gg, bb),
                           elemData.dotNormalFace(gg, bb));
              }
            MoFEMFunctionReturn(0);
          };

          auto check_continuity_of_h1_base = [&](auto &vol_data) {
            MoFEMFunctionBegin;

            if (vol_data.size1() != elemData.shapeFunH1Values.size1())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of integration points");

            if (vol_data.size2() != elemData.shapeFunH1Values.size2())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of base functions");
            constexpr double eps = 1e-12;
            for (size_t gg = 0; gg != vol_data.size1(); ++gg)
              for (size_t bb = 0; bb != vol_data.size2(); ++bb) {
                const double error = std::abs(
                    vol_data(gg, bb) - elemData.shapeFunH1Values(gg, bb));
                if (error > eps)
                  SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                           "Inconsistency %3.4e != %3.4e", vol_data(gg, bb),
                           elemData.shapeFunH1Values(gg, bb));
              }
            MoFEMFunctionReturn(0);
          };

          if (elemData.dotNormalEleLeft.size2() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleLeft);
          else if (elemData.dotNormalEleRight.size2() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleRight);
          else if (elemData.shapeFunH1VolSide.size2() != 0)
            CHKERR check_continuity_of_h1_base(elemData.shapeFunH1VolSide);
        }
        MoFEMFunctionReturnHot(0);
      }
    };

    struct OnContactSideSlave
        : public ContactPrismElementForcesAndSourcesCore::UserDataOperator {

      VolumeElementForcesAndSourcesCoreOnContactPrismSide volSideFe;
      struct OpVolSide : public VolumeElementForcesAndSourcesCoreOnContactPrismSide::
                             UserDataOperator {

        CommonData &elemData;
        OpVolSide(CommonData &elem_data)
            : VolumeElementForcesAndSourcesCoreOnContactPrismSide::UserDataOperator(
                  "F2", UserDataOperator::OPROW),
              elemData(elem_data) {}
        MoFEMErrorCode doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
          MoFEMFunctionBeginHot;

          if (type == MBTRI && side == getFaceSideNumber()) {

            MatrixDouble diff =
                getCoordsAtGaussPts() - getMasterCoordsAtGaussPts();
            constexpr double eps = 1e-12;
            if (norm_inf(diff) > eps)
              SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                      "coordinates at integration pts are different");

            const size_t nb_dofs = data.getN().size2() / 3;
            const size_t nb_integration_pts = data.getN().size1();

            if (data.getSpace() == H1) {
              MatrixDouble *ptr_elem_data = nullptr;
              ptr_elem_data = &elemData.shapeFunH1VolSide;

              MatrixDouble &elem_data = *ptr_elem_data;
              elem_data.resize(nb_integration_pts, nb_dofs, false);
              for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
                FTensor::Tensor0<double *> t_base(&data.getN()(gg, 0));
                for (size_t bb = 0; bb != nb_dofs; ++bb) {
                  elem_data(gg, bb) = t_base;
                  ++t_base;
                }
              }
            }
          }
          MoFEMFunctionReturnHot(0);
        }
      };

      CommonData &elemData;
      OnContactSideSlave(MoFEM::Interface &m_field, CommonData &elem_data)
          : ContactPrismElementForcesAndSourcesCore::UserDataOperator(
                "F2", UserDataOperator::OPROW,
                ContactPrismElementForcesAndSourcesCore::UserDataOperator::
                    FACESLAVE),
            volSideFe(m_field), elemData(elem_data) {
        volSideFe.getOpPtrVector().push_back(
            new OnContactSideSlave::OpVolSide(elemData));
      }

      MoFEMErrorCode doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {

        MoFEMFunctionBeginHot;
        if (type == MBTRI && side == 0) {
          const size_t nb_dofs = data.getN().size2() / 3;
          const size_t nb_integration_pts = data.getN().size1();
          if (data.getSpace() == H1) {
            elemData.shapeFunH1Values.resize(nb_integration_pts, nb_dofs,
                                             false);
            elemData.shapeFunH1VolSide.resize(nb_integration_pts, nb_dofs,
                                              false);

            for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
              FTensor::Tensor0<double *> t_base(&data.getN()(gg, 0));
              for (size_t bb = 0; bb != nb_dofs; ++bb) {
                elemData.shapeFunH1Values(gg, bb) = t_base;
                ++t_base;
              }
            }
          }

          std::string side_fe_name = "V1";
          const EntityHandle tri_slave = getSideEntity(4, MBTRI);
          CHKERR loopSideVolumes(side_fe_name, volSideFe, 3, tri_slave);

          auto check_continuity_of_base = [&](auto &vol_dot_data) {
            MoFEMFunctionBegin;

            if (vol_dot_data.size1() != elemData.dotNormalFace.size1())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of integration points");

            if (vol_dot_data.size2() != elemData.dotNormalFace.size2())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of base functions");
            constexpr double eps = 1e-12;
            for (size_t gg = 0; gg != vol_dot_data.size1(); ++gg)
              for (size_t bb = 0; bb != vol_dot_data.size2(); ++bb) {
                const double error = std::abs(vol_dot_data(gg, bb) -
                                              elemData.dotNormalFace(gg, bb));
                if (error > eps)
                  SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                           "Inconsistency %3.4e != %3.4e", vol_dot_data(gg, bb),
                           elemData.dotNormalFace(gg, bb));
              }
            MoFEMFunctionReturn(0);
          };

          auto check_continuity_of_h1_base = [&](auto &vol_data) {
            MoFEMFunctionBegin;

            if (vol_data.size1() != elemData.shapeFunH1Values.size1())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of integration points");

            if (vol_data.size2() != elemData.shapeFunH1Values.size2())
              SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                      "Inconsistent number of base functions");
            constexpr double eps = 1e-12;
            for (size_t gg = 0; gg != vol_data.size1(); ++gg)
              for (size_t bb = 0; bb != vol_data.size2(); ++bb) {
                const double error = std::abs(
                    vol_data(gg, bb) - elemData.shapeFunH1Values(gg, bb));
                if (error > eps)
                  SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                           "Inconsistency %3.4e != %3.4e", vol_data(gg, bb),
                           elemData.shapeFunH1Values(gg, bb));
              }
            MoFEMFunctionReturn(0);
          };

          if (elemData.dotNormalEleLeft.size2() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleLeft);
          else if (elemData.dotNormalEleRight.size2() != 0)
            CHKERR check_continuity_of_base(elemData.dotNormalEleRight);
          else if (elemData.shapeFunH1VolSide.size2() != 0)
            CHKERR check_continuity_of_h1_base(elemData.shapeFunH1VolSide);
        }
        MoFEMFunctionReturnHot(0);
      }
    };

    CommonData elem_data;

    ContactPrismElementForcesAndSourcesCore contact_prism_fe_master(m_field);
    ContactPrismElementForcesAndSourcesCore contact_prism_fe_slave(m_field);

    // OnContactSideMaster
    contact_prism_fe_master.getOpPtrVector().push_back(
        new OnContactSideMaster(m_field, elem_data));
    // OnContactSideSlave
    contact_prism_fe_slave.getOpPtrVector().push_back(
        new OnContactSideSlave(m_field, elem_data));

    CHKERR m_field.loop_finite_elements("P1", "C2", contact_prism_fe_master);
    CHKERR m_field.loop_finite_elements("P1", "C2", contact_prism_fe_slave);
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();

  return 0;
}
