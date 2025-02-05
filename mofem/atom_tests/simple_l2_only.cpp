/**
 * \file simple_l2_only.cpp
 * \ingroup mofem_simple_interface
 * \example simple_l2_only.cpp
 *
 * Test iterating over boundary and skeleton elements only when L2 field is
 * presents  on the domain.
 *
 */

#include <MoFEM.hpp>
using namespace MoFEM;

static char help[] = "...\n\n";

constexpr int SPACE_DIM = 2;

template <int DIM> struct ElementsAndOps {};

template <> struct ElementsAndOps<2> {
  using DomainEle = PipelineManager::FaceEle;
  using DomainEleOp = DomainEle::UserDataOperator;
  using BoundaryEle = PipelineManager::EdgeEle;
  using BoundaryEleOp = BoundaryEle::UserDataOperator;
};

using DomainEle = ElementsAndOps<SPACE_DIM>::DomainEle;
using DomainEleOp = DomainEle::UserDataOperator;
using BoundaryEle = ElementsAndOps<SPACE_DIM>::BoundaryEle;
using BoundaryEleOp = BoundaryEle::UserDataOperator;
using EntData = EntitiesFieldData::EntData;

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
      CHKERR simple_interface->loadFile();

      CHKERR simple_interface->addDomainField("FIELD", L2,
                                              AINSWORTH_LEGENDRE_BASE, 1);

      CHKERR simple_interface->addMeshsetField("GLOBAL", NOFIELD, NOBASE, 1);
      CHKERR simple_interface->addMeshsetField("FIELD", L2,
                                               AINSWORTH_LEGENDRE_BASE, 1);

      // I add vols to meshset, that is to integrate on side of mFE. That make
      // "FIELD" adjacent to "GLOBAL", and all other "FIELD" adjacent to
      // "FIELD". That would create dense matrix. In principle you will add only
      // small subset of "FIELD" entities to "GLOBAL" meshset.
      Range vols;
      CHKERR m_field.get_moab().get_entities_by_dimension(0, SPACE_DIM, vols);
      simple_interface->getMeshsetFiniteElementEntities().push_back(
          vols); // create one meshset element

      simple_interface->getAddBoundaryFE() = true;
      simple_interface->getAddSkeletonFE() = true;

      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD", 1);
      // setup problem
      CHKERR simple_interface->setUp();

      int count_skeleton_fe;
      int count_side_fe;
      int count_meshset_fe;
      int count_meshset_side_fe;

      PipelineManager *pipeline_mng = m_field.getInterface<PipelineManager>();

      // Create OP for side FE
      auto op_side_fe = new DomainEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_side_fe->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                      EntityType type,
                                      EntitiesFieldData::EntData &data) {
        auto domain_op = static_cast<DomainEleOp *>(op_ptr);
        MoFEMFunctionBegin;

        MOFEM_LOG("SELF", Sev::verbose)
            << "Side element name [ " << count_side_fe << " ] "
            << domain_op->getFEName();

        ++count_side_fe;

        MoFEMFunctionReturn(0);
      };

      // Create side FE
      auto side_fe = boost::make_shared<DomainEle>(m_field);
      side_fe->getOpPtrVector().push_back(op_side_fe);

      // Create boundary FE OP

      auto do_work_rhs = [&](DataOperator *op_ptr, int side, EntityType type,
                             EntitiesFieldData::EntData &data) {
        auto bdy_op = static_cast<BoundaryEleOp *>(op_ptr);
        MoFEMFunctionBegin;

        MOFEM_LOG("SELF", Sev::verbose)
            << "Element name  [ " << count_skeleton_fe << " ] "
            << bdy_op->getFEName();

        CHKERR bdy_op->loopSide(simple_interface->getDomainFEName(),
                                side_fe.get(), SPACE_DIM);

        ++count_skeleton_fe;

        MoFEMFunctionReturn(0);
      };

      auto op_bdy_fe = new BoundaryEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_bdy_fe->doWorkRhsHook = do_work_rhs;

      auto op_skeleton_fe = new BoundaryEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_skeleton_fe->doWorkRhsHook = do_work_rhs;

      // create meshset fe
      auto op_meshset_side_fe = new DomainEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_meshset_side_fe->doWorkRhsHook =
          [&](DataOperator *op_ptr, int side, EntityType type,
              EntitiesFieldData::EntData &data) {
            auto domain_op = static_cast<DomainEleOp *>(op_ptr);
            MoFEMFunctionBegin;

            MOFEM_LOG("SELF", Sev::verbose)
                << "Side element name [ " << count_side_fe << " ] "
                << domain_op->getFEName();

            ++count_meshset_side_fe;

            MoFEMFunctionReturn(0);
          };

      auto meshset_side_fe = boost::make_shared<DomainEle>(m_field);
      meshset_side_fe->getOpPtrVector().push_back(op_meshset_side_fe);

      auto op_meshset_fe = new ForcesAndSourcesCore::UserDataOperator(
          "GLOBAL", ForcesAndSourcesCore::UserDataOperator::OPROW);
      op_meshset_fe->doWorkRhsHook = [&](DataOperator *op_ptr, int side,
                                         EntityType type,
                                         EntitiesFieldData::EntData &data) {
        MoFEMFunctionBegin;
        MOFEM_LOG("SELF", Sev::inform)
            << "Meshset element name " << data.getIndices();

        CHKERR op_meshset_fe->loopSide(simple_interface->getDomainFEName(),
                                       meshset_side_fe.get(), SPACE_DIM, 0,
                                       boost::make_shared<Range>(vols));

        ++count_meshset_fe;
        MoFEMFunctionReturn(0);
      };

      // Count boundary
      count_skeleton_fe = 0;
      count_side_fe = 0;
      count_meshset_fe = 0;
      count_meshset_side_fe = 0;

      pipeline_mng->getOpBoundaryRhsPipeline().push_back(op_bdy_fe);
      pipeline_mng->getOpSkeletonRhsPipeline().push_back(op_skeleton_fe);
      pipeline_mng->getOpMeshsetRhsPipeline().push_back(op_meshset_fe);
      pipeline_mng->loopFiniteElements();

      MOFEM_LOG("SELF", Sev::inform)
          << "Number of elements " << count_skeleton_fe;
      MOFEM_LOG("SELF", Sev::inform)
          << "Number of side elements " << count_side_fe;
      MOFEM_LOG("SELF", Sev::inform)
          << "Number of meshset elements " << count_meshset_fe;
      MOFEM_LOG("SELF", Sev::inform)
          << "Number of meshset side elements " << count_meshset_side_fe;

      if (count_skeleton_fe != 16)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of FEs");
      if (count_side_fe != 24)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of side FEs");
      if (count_meshset_fe != 1)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of side FEs");
      if (count_meshset_side_fe != 8)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of side FEs");
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
