/** \file loop_entities.cpp
  \example loop_entities.cpp
  \brief Make a loop over entities

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

struct TestEntityMethod : EntityMethod {
  TestEntityMethod(const Range &ents)
      : EntityMethod(), allEnts(ents) {}

  MoFEMErrorCode preProcess() {
    MoFEMFunctionBeginHot;
    entsIt = allEnts.begin();
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode operator()() {
    MoFEMFunctionBeginHot;
    if (entPtr->getEnt() != *entsIt)
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "Entity and entity iterator should be the same");
    ++entsIt;
    MoFEMFunctionReturnHot(0);
  }

  MoFEMErrorCode postProcess() {
    MoFEMFunctionBeginHot;
    if(entsIt != allEnts.end())
      SETERRQ(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
              "All entities should be iterated");
    MoFEMFunctionReturnHot(0);
  }

private:
  const Range allEnts;
  Range::const_iterator entsIt;
};

int main(int argc, char *argv[]) {
  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    // read mesh and create moab and mofem data structures
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
      CHKERR simple_interface->setFieldOrder("FIELD1", 4);
      // setup problem
      CHKERR simple_interface->setUp();

      Range all_ents;
      CHKERR m_field.get_moab().get_entities_by_handle(root_set, all_ents);
      all_ents = subtract(all_ents, all_ents.subset_by_type(MBENTITYSET));

      auto testingEntitiesInDatabase = [&]() {
        MoFEMFunctionBegin;
        TestEntityMethod method(all_ents);
        CHKERR m_field.loop_entities("FIELD1", method);
        MoFEMFunctionReturn(0);
      };

      CHKERR testingEntitiesInDatabase();

      auto testingEntitiesInDatabaseSubRange = [&]() {
        MoFEMFunctionBegin;
        Range edges = all_ents.subset_by_type(MBEDGE);
        TestEntityMethod method(edges);
        CHKERR m_field.loop_entities("FIELD1", method, &edges);
        MoFEMFunctionReturn(0);
      };

      CHKERR testingEntitiesInDatabaseSubRange();

      auto dm = simple_interface->getDM();

      auto testingEntitiesInDM = [&]() {
        MoFEMFunctionBegin;
        TestEntityMethod method(all_ents);
        const MoFEM::Problem *prb_ptr;
        CHKERR DMMoFEMGetProblemPtr(dm, &prb_ptr);
        CHKERR m_field.loop_entities(prb_ptr, "FIELD1", ROW, method, 0,
                                     m_field.get_comm_size());
        MoFEMFunctionReturn(0);
      };

      CHKERR testingEntitiesInDM();

    }
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, ect.
  MoFEM::Core::Finalize();

  return 0;
}
