/**
 * \file test_cache_entities.cpp
 * \example test_cache_entities.cpp
 *
 * Tetsing entities cache. Entities acache is data abstraction enabling user to
 * pass data between operators, finite elements, and problems, in convinient
 * way.
 *
 * It can be used to store indices, or history variables, and any other data,
 * which I can not imagine at that moment of time, and how to make use of it
 * depends on user creativity and imagination.
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

using VolEle = VolumeElementForcesAndSourcesCore;
using VolOp = VolEle::UserDataOperator;
using FaceEle = FaceElementForcesAndSourcesCore;
using FaceOp = FaceEle::UserDataOperator;

struct OpVolume : public VolOp {
  OpVolume(const std::string &field_name) : VolOp(field_name, OPROW) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;

    MoFEMFunctionReturn(0);
  }

};

struct OpFace : public FaceOp {
  OpFace(const std::string &field_name) : FaceOp(field_name, OPROW) {}
  MoFEMErrorCode doWork(int side, EntityType type,
                        DataForcesAndSourcesCore::EntData &data) {

    MoFEMFunctionBegin;
    MoFEMFunctionReturn(0);
  }
};

struct VolRule {
  int operator()(int, int, int) const { return 1; }
};
struct FaceRule {
  int operator()(int, int, int) const { return 1; }
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
      // add fields
      CHKERR simple_interface->addDomainField("MESH_NODE_POSITIONS", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      CHKERR simple_interface->addBoundaryField("MESH_NODE_POSITIONS", H1,
                                                AINSWORTH_LEGENDRE_BASE, 3);
      // set fields order
      CHKERR simple_interface->setFieldOrder("MESH_NODE_POSITIONS", 1);
      // setup problem
      CHKERR simple_interface->setUp();
      // Project mesh coordinate on mesh
      Projection10NodeCoordsOnField ent_method(m_field, "MESH_NODE_POSITIONS");
      CHKERR m_field.loop_dofs("MESH_NODE_POSITIONS", ent_method);

      // get dm
      auto dm = simple_interface->getDM();
      // create elements
      auto domain_fe = boost::make_shared<VolEle>(m_field);
      auto boundary_fe = boost::make_shared<FaceEle>(m_field);
      // set integration rule
      domain_fe->getRuleHook = VolRule();
      boundary_fe->getRuleHook = FaceRule();

      // set operator to the volume element
      domain_fe->getOpPtrVector().push_back(
          new OpVolume("MESH_NODE_POSITIONS"));
      // set operator to the face element
      boundary_fe->getOpPtrVector().push_back(
          new OpFace("MESH_NODE_POSITIONS"));
      // make integration in volume (here real calculations starts)
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getDomainFEName(),
                                      domain_fe);
      // make integration on boundary
      CHKERR DMoFEMLoopFiniteElements(dm, simple_interface->getBoundaryFEName(),
                                      boundary_fe);


    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();

  return 0;
}
