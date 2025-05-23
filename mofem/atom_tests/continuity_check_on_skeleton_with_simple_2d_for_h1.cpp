/**
 * \file continuity_check_on_skeleton_with_simple_2d_for_h1.cpp
 * \ingroup mofem_simple_interface
  \example continuity_check_on_skeleton_with_simple_2d_for_h1.cpp
 *
 * \brief Integration on skeleton for 2d
 *
 * Teting integration on skeleton and checking of continuity of hcurl space on
 * edges.
 */

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

using FaceEleOnSide = MoFEM::FaceElementForcesAndSourcesCoreOnSide;
using EdgeEle = MoFEM::EdgeElementForcesAndSourcesCore;
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
        : FaceEleOnSideOp("FIELD", UserDataOperator::OPROW),
          elemData(elem_data) {}

    MoFEMErrorCode doWork(int side, EntityType type,
                          EntitiesFieldData::EntData &data) {
      MoFEMFunctionBeginHot;

      if (type == MBEDGE && side == getEdgeSideNumber()) {

        MatrixDouble diff = getCoordsAtGaussPts() - getEdgeCoordsAtGaussPts();

        const double eps = 1e-12;
        if (norm_inf(diff) > eps) {
          MOFEM_LOG_ATTRIBUTES("ATOM",
                               LogManager::BitLineID | LogManager::BitScope);
          MOFEM_LOG_FUNCTION();
          MOFEM_LOG("ATOM", Sev::error)
              << "Quad coords: " << getCoordsAtGaussPts();
          MOFEM_LOG("ATOM", Sev::error)
              << "Edge coords: " << getEdgeCoordsAtGaussPts();
          MOFEM_LOG("ATOM", Sev::error) << "Difference: " << diff;
          SETERRQ(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                  "Coordinates at integration pts are different");
        }

        const size_t nb_dofs = data.getN().size2();
        const size_t nb_integration_pts = data.getN().size1();

        auto t_base = data.getFTensor0N();
        MatrixDouble *ptr_dot_elem_data = nullptr;
        if (getFEMethod()->nInTheLoop == 0)
          ptr_dot_elem_data = &elemData.dotEleLeft;
        else
          ptr_dot_elem_data = &elemData.dotEleRight;
        MatrixDouble &dot_elem_data = *ptr_dot_elem_data;
        dot_elem_data.resize(nb_integration_pts, nb_dofs, false);

        for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
          for (size_t bb = 0; bb != nb_dofs; ++bb) {
            dot_elem_data(gg, bb) = t_base;
            ++t_base;
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
    faceSideFe.getOpPtrVector().push_back(new SkeletonFE::OpFaceSide(elemData));
  }

  MoFEMErrorCode doWork(int side, EntityType type,
                        EntitiesFieldData::EntData &data) {

    MoFEMFunctionBeginHot;
    if (type == MBEDGE) {

      const size_t nb_dofs = data.getN().size2();
      const size_t nb_integration_pts = data.getN().size1();

      MOFEM_LOG("ATOM", Sev::noisy)
          << "Cords at integration points" << getCoordsAtGaussPts();

      auto t_base = data.getFTensor0N();
      elemData.dotEdge.resize(nb_integration_pts, nb_dofs, false);
      elemData.dotEleLeft.resize(nb_integration_pts, 0, false);
      elemData.dotEleRight.resize(nb_integration_pts, 0, false);

      for (size_t gg = 0; gg != nb_integration_pts; ++gg) {
        for (size_t bb = 0; bb != nb_dofs; ++bb) {
          elemData.dotEdge(gg, bb) = t_base;
          ++t_base;
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
      // Add ghost cells needed for evaluation on skeleton
      simple_interface->getAddSkeletonFE() = true;

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
      CHKERR simple_interface->addDomainField("FIELD", H1, base, 1);
      // CHKERR simple_interface->addDomainField("TEST_FIELD", L2,
      //                                         AINSWORTH_LEGENDRE_BASE, 1);
      CHKERR simple_interface->addSkeletonField("FIELD", H1, base, 1);
      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD", 2);
      // CHKERR simple_interface->setFieldOrder("TEST_FIELD", 1);
      // setup problem
      CHKERR simple_interface->setUp();
      // get dm
      auto dm = simple_interface->getDM();

      // create elements
      CommonData elem_data;
      boost::shared_ptr<EdgeEle> skeleton_fe =
          boost::shared_ptr<EdgeEle>(new EdgeEle(m_field));
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