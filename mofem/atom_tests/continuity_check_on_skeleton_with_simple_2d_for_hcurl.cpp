/**
 * \file continuity_check_on_skeleton_with_simple_2d_for_hcurl.cpp
 * \ingroup mofem_simple_interface
 * \example continuity_check_on_skeleton_with_simple_2d_for_hcurl.cpp
 *
 * \brief Integration on skeleton for 2d
 * 
 * Teting integration on skeleton and checking of continuity of hcurl space on
 * edges.
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

using FaceEleOnSide = MoFEM::FaceElementForcesAndSourcesCoreOnSideSwitch<0>;

using EdgeEle = MoFEM::EdgeElementForcesAndSourcesCoreSwitch<
    EdgeElementForcesAndSourcesCore::NO_HO_GEOMETRY |
    EdgeElementForcesAndSourcesCoreBase::NO_COVARIANT_TRANSFORM_HCURL>;

using FaceEleOnSideOp = FaceEleOnSide::UserDataOperator;
using EdgeEleOp = EdgeEle::UserDataOperator;

struct CommonData {
  MatrixDouble dotEdge;
  MatrixDouble dotEleLeft;
  MatrixDouble dotEleRight;
};

struct SkeletonFE : public EdgeEleOp {

  struct OpFaceSide : public FaceEleOnSideOp {

    CommonData &elemData;
    OpFaceSide(CommonData &elem_data)
        : FaceEleOnSideOp("FIELD", UserDataOperator::OPROW), elemData(elem_data) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          DataForcesAndSourcesCore::EntData &data) {
      MoFEMFunctionBeginHot;

      if (type == MBEDGE && side == getEdgeSideNumber()) {

        MatrixDouble diff = getCoordsAtGaussPts() - getEdgeCoordsAtGaussPts();

        const double eps = 1e-12;
        if (norm_inf(diff) > eps)
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "coordinates at integration pts are different");

        const size_t nb_dofs = data.getN().size2() / 3;
        const size_t nb_integration_pts = data.getN().size1();

        auto t_tangent = getFTensor1Direction();
        auto t_hcurl_base = data.getFTensor1N<3>();
        FTensor::Index<'i', 3> i;
        MatrixDouble *ptr_dot_elem_data = nullptr;
        if (getFEMethod()->nInTheLoop == 0)
          ptr_dot_elem_data = &elemData.dotEleLeft;
        else
          ptr_dot_elem_data = &elemData.dotEleRight;
        MatrixDouble &dot_elem_data = *ptr_dot_elem_data;
        dot_elem_data.resize(nb_integration_pts, nb_dofs, false);

        for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
          for (size_t bb = 0; bb != nb_dofs; ++bb) {
            dot_elem_data(gg, bb) = t_tangent(i) * t_hcurl_base(i);
            ++t_hcurl_base;
          }
        }
      }
      MoFEMFunctionReturnHot(0);
    }
  };

  CommonData &elemData;
  FaceEleOnSide faceSideFe;

  SkeletonFE(MoFEM::Interface &m_field, CommonData &elem_data)
      : EdgeEle::UserDataOperator("FIELD", UserDataOperator::OPROW),
        faceSideFe(m_field), elemData(elem_data) {
    faceSideFe.getOpPtrVector().push_back(
        new OpHOSetCovariantPiolaTransformOnFace3D(HCURL));
    faceSideFe.getOpPtrVector().push_back(new SkeletonFE::OpFaceSide(elemData));
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBeginHot;
    if (type == MBEDGE) {

      const size_t nb_dofs = data.getN().size2() / 3;
      const size_t nb_integration_pts = data.getN().size1();

      auto t_tangent = getFTensor1Direction();
      auto t_hcurl_base = data.getFTensor1N<3>();
      FTensor::Index<'i', 3> i;
      elemData.dotEdge.resize(nb_integration_pts, nb_dofs, false);
      elemData.dotEleLeft.resize(nb_integration_pts, 0, false);
      elemData.dotEleRight.resize(nb_integration_pts, 0, false);

      for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
        for (size_t bb = 0; bb != nb_dofs; ++bb) {
          elemData.dotEdge(gg, bb) = t_tangent(i) * t_hcurl_base(i);
          ++t_hcurl_base;
        }
      }

      CHKERR loopSideFaces("dFE", faceSideFe);

      auto check_continuity_of_base = [&](auto &vol_dot_data) {
        MoFEMFunctionBegin;

        if (vol_dot_data.size1() != elemData.dotEdge.size1())
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Inconsistent number of integration points");

        if (vol_dot_data.size2() != elemData.dotEdge.size2())
          SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                  "Inconsistent number of base functions");
        const double eps = 1e-12;
        for (size_t gg = 0; gg != vol_dot_data.size1(); ++gg)
          for (size_t bb = 0; bb != vol_dot_data.size2(); ++bb) {
            const double error =
                std::abs(vol_dot_data(gg, bb) - elemData.dotEdge(gg, bb));
            if (error > eps)
              SETERRQ4(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "Inconsistency (%d, %d) %3.4e != %3.4e", gg, bb,
                       vol_dot_data(gg, bb), elemData.dotEdge(gg, bb));
            else
              MOFEM_LOG("ATOM", Sev::noisy) << "Ok";

          }
        MoFEMFunctionReturn(0);
      };

      if (elemData.dotEleLeft.size2() != 0)
        CHKERR check_continuity_of_base(elemData.dotEleLeft);
      else if (elemData.dotEleRight.size2() != 0)
        CHKERR check_continuity_of_base(elemData.dotEleRight);

    }
    MoFEMFunctionReturnHot(0);
  }
};

int main(int argc, char *argv[]) {

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Create MoAB database
    moab::Core moab_core;
    moab::Interface &moab = moab_core;

    // Create MoFEM database and link it to MoAB
    MoFEM::Core mofem_core(moab);
    MoFEM::Interface &m_field = mofem_core;

    auto core_log = logging::core::get();
    core_log->add_sink(
        LogManager::createSink(LogManager::getStrmSelf(), "ATOM"));
    LogManager::setLog("ATOM");

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Simple interface
    Simple *simple_interface;
    CHKERR m_field.getInterface(simple_interface);
    {
      // get options from command line
      CHKERR simple_interface->getOptions();
      // load mesh file
      CHKERR simple_interface->loadFile("");

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

      // add fields
      auto base = get_base();
      CHKERR simple_interface->addDomainField("FIELD", HCURL, base, 1);
      CHKERR simple_interface->addSkeletonField("FIELD", HCURL, base, 1);
      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD", 3);
      // setup problem
      CHKERR simple_interface->setUp();
      // get dm
      auto dm = simple_interface->getDM();

      // create elements
      CommonData elem_data;
      boost::shared_ptr<EdgeEle> skeleton_fe =
          boost::shared_ptr<EdgeEle>(new EdgeEle(m_field));

      skeleton_fe->getOpPtrVector().push_back(
          new OpHOSetContravariantPiolaTransformOnEdge3D(HCURL));
      skeleton_fe->getOpPtrVector().push_back(
          new SkeletonFE(m_field, elem_data));

      // iterate skeleton finite elements
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getSkeletonFEName(),
                                      skeleton_fe);
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}