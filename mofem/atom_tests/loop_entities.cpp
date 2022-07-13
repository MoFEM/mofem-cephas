/** \file loop_entities.cpp
  \example loop_entities.cpp
  \brief Make a loop over entities

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
