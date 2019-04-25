/** \file segments_distance.cpp
 * \brief test segments distance
 *
 * \ingroup mesh_cut
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

static char help[] = "testing mesh cut test\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // Read mesh and create moab and mofem data structures
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;

    MoFEM::Core core(moab, PETSC_COMM_WORLD);
    MoFEM::Interface &m_field = core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    EntityHandle root_set = moab.get_root_set();

    // Simple interface
    Simple *simple_interface;
    CHKERR m_field.getInterface(simple_interface);
    {

      // get options from command line
      CHKERR simple_interface->getOptions();
      // load mesh file
      CHKERR simple_interface->loadFile();
      // add fields
      CHKERR simple_interface->addDomainField("FIELD1", H1,
                                              AINSWORTH_LEGENDRE_BASE, 3);
      // set fields order
      CHKERR simple_interface->setFieldOrder("FIELD1", 1);
      // setup problem
      CHKERR simple_interface->setUp();

      CHKERR m_field.getInterface<FieldEvaluatorInterface>()->buildTree3D(
          simple_interface->getDomainFEName());

      std::array<double, 3> point = {0, 0, 0};
      const double dist = 0.1;

      std::array<double, 12> eval_points = {0.0,  0.0, 0.0, 0.0, 0.01, 0.0,
                                            -0.1, 0.0, 0.0, 0.1, 0.0,  0.0};

      DM dm;
      CHKERR simple_interface->getDM(&dm);
      const MoFEM::Problem *prb_ptr;
      CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);

      VolumeElementForcesAndSourcesCore vol_ele(m_field);

      using VolOp = VolumeElementForcesAndSourcesCore::UserDataOperator;

      struct MyOp: public VolOp {

        MyOp(): VolOp("FIELD1",OPROW) {}

        MoFEMErrorCode doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {
          MoFEMFunctionBegin;
          if (type == MBVERTEX) {

            std::cout << "Integration pts" << std::endl;
            std::cout << getGaussPts() << endl;

            std::cout << "Global coordinates " << endl;
            std::cout << getCoordsAtGaussPts() << std::endl;

          }
          MoFEMFunctionReturn(0);
        }
      };

      vol_ele.getOpPtrVector().push_back(new MyOp());

      CHKERR m_field.getInterface<FieldEvaluatorInterface>()
          ->evalFEAtThePoint3D(
              &point[0], dist, &eval_points[0], eval_points.size() / 3, 1e-12,
              prb_ptr->getName(), simple_interface->getDomainFEName(), vol_ele,
              m_field.get_comm_rank(), m_field.get_comm_rank(), MF_EXIST,
              QUIET);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
