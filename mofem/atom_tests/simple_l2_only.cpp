/**
 * \file simple_l2_only.cpp
 * \ingroup mofem_simple_interface
 * \example simple_l2_only.cpp
 *
 * Test iterating over boundary and skeleton elements only when L2 field is
 * presents  on the domain.
 *
 */

/* MIT License
 *
 * Copyright (c) 2022
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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

      simple_interface->getAddBoundaryFE() = true;
      simple_interface->getAddSkeletonFE() = true;

      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD", 1);
      // setup problem
      CHKERR simple_interface->setUp();

      int count_fe;
      int count_side_fe;

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

        MOFEM_LOG("SELF", Sev::verbose) << "Element name  [ " << count_fe
                                        << " ] " << bdy_op->getFEName();

        CHKERR bdy_op->loopSide(simple_interface->getDomainFEName(),
                                side_fe.get(), SPACE_DIM);

        ++count_fe;

        MoFEMFunctionReturn(0);
      };

      auto op_bdy_fe = new BoundaryEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_bdy_fe->doWorkRhsHook = do_work_rhs;

      auto op_skeleton_fe = new BoundaryEleOp(NOSPACE, DomainEleOp::OPSPACE);
      op_skeleton_fe->doWorkRhsHook = do_work_rhs;

      // Count boundary
      count_fe = 0;
      count_side_fe = 0;

      pipeline_mng->getOpBoundaryRhsPipeline().push_back(op_bdy_fe);
      pipeline_mng->getOpSkeletonRhsPipeline().push_back(op_skeleton_fe);
      pipeline_mng->loopFiniteElements();

      MOFEM_LOG("SELF", Sev::inform) << "Number of elements " << count_fe;
      MOFEM_LOG("SELF", Sev::inform)
          << "Number of side elements " << count_side_fe;

      if (count_fe != 16)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of FEs");
      if (count_side_fe != 24)
        SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                "Wrong numbers of side FEs");
    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
