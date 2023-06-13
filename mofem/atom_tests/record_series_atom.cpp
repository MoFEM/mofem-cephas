

#include <MoFEM.hpp>

namespace bio = boost::iostreams;
using bio::stream;
using bio::tee_device;

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, PETSC_NULL, help);

  try {

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // Reade parameters from line command
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
    PetscInt order = 3;
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-my_order", &order, &flg);
#else
    CHKERR PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-my_order", &order,
                              &flg);
#endif

    // Read mesh to MOAB
    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM 
    // Starting from one  MoFEM::CoreTmp<1> for testing purposes
    MoFEM::CoreTmp<0> core(moab);
    MoFEM::Interface &m_field = core;

    // stl::bitset see for more details
    auto bit_level0 = BitRefLevel().set(0);
    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));

    /***/
    // Define problem

    // Fields
    CHKERR m_field.add_field("FIELD_A", H1, AINSWORTH_LEGENDRE_BASE, 3);
    CHKERR m_field.add_field("FIELD_B", H1, AINSWORTH_LEGENDRE_BASE, 3);

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_A");
    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_B");

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_A", order);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_A", 1);

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_B", 1);

    // build field
    CHKERR m_field.build_fields();

    CHKERR m_field.getInterface<FieldBlas>()->setField(0, MBVERTEX, "FIELD_B");
    CHKERR m_field.getInterface<FieldBlas>()->setField(1, MBVERTEX, "FIELD_A");

    SeriesRecorder *recorder_ptr;
    CHKERR m_field.getInterface(recorder_ptr);

    CHKERR recorder_ptr->add_series_recorder("TEST_SERIES1");
    CHKERR recorder_ptr->add_series_recorder("TEST_SERIES2");

    // initialize
    CHKERR recorder_ptr->initialize_series_recorder("TEST_SERIES1");

    CHKERR recorder_ptr->record_begin("TEST_SERIES1");
    CHKERR recorder_ptr->record_field("TEST_SERIES1", "FIELD_B", bit_level0,
                                      bit_level0);
    CHKERR recorder_ptr->record_end("TEST_SERIES1", 1);

    CHKERR m_field.getInterface<FieldBlas>()->fieldAxpy(1., "FIELD_A",
                                                        "FIELD_B");
    CHKERR recorder_ptr->record_begin("TEST_SERIES1");
    CHKERR recorder_ptr->record_field("TEST_SERIES1", "FIELD_B", bit_level0,
                                      bit_level0);

    CHKERR recorder_ptr->initialize_series_recorder("TEST_SERIES2");
    CHKERR recorder_ptr->record_begin("TEST_SERIES2");
    CHKERR recorder_ptr->record_field("TEST_SERIES2", "FIELD_A", bit_level0,
                                      bit_level0);
    CHKERR recorder_ptr->record_field("TEST_SERIES2", "FIELD_B", bit_level0,
                                      bit_level0);
    CHKERR recorder_ptr->record_end("TEST_SERIES2", 1);
    CHKERR recorder_ptr->finalize_series_recorder("TEST_SERIES2");

    CHKERR recorder_ptr->record_end("TEST_SERIES1", 2);

    // finalize
    CHKERR recorder_ptr->finalize_series_recorder("TEST_SERIES1");
    CHKERR recorder_ptr->print_series_steps();

    CHKERR m_field.getInterface<FieldBlas>()->fieldScale(2, "FIELD_A");

    // Setting MoFEM::CoreTmp<-1> for testing purposes
    MoFEM::CoreTmp<-1> core2(moab);
    MoFEM::Interface &m_field2 = core2;

    // build field
    CHKERR m_field2.build_fields();

    typedef tee_device<std::ostream, std::ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;
    std::ofstream ofs("record_series_atom.txt");
    TeeDevice my_tee(std::cout, ofs);
    TeeStream my_split(my_tee);

    SeriesRecorder *recorder2_ptr;
    CHKERR m_field2.getInterface(recorder2_ptr);
    CHKERR recorder2_ptr->print_series_steps();

    my_split << "TEST_SERIES1" << std::endl;
    for (_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder2_ptr, "TEST_SERIES1",
                                            sit)) {

      CHKERR recorder2_ptr->load_series_data("TEST_SERIES1",
                                             sit->get_step_number());

      my_split << "next step:\n";
      my_split << *sit << std::endl;

      {
        DofEntity_multiIndex_uid_view dofs_view;
        for (_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2, "FIELD_B", dof)) {
          dofs_view.insert(*dof);
        }
        for (DofEntity_multiIndex_uid_view::iterator dit = dofs_view.begin();
             dit != dofs_view.end(); dit++) {
          my_split << **dit << endl;
        }
      }
    }

    my_split << "TEST_SERIES2" << std::endl;
    for (_IT_SERIES_STEPS_BY_NAME_FOR_LOOP_(recorder2_ptr, "TEST_SERIES2",
                                            sit)) {

      CHKERR recorder2_ptr->load_series_data("TEST_SERIES2",
                                             sit->get_step_number());

      my_split << "next step:\n";
      {
        DofEntity_multiIndex_uid_view dofs_view;
        for (_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2, "FIELD_A", dof)) {
          dofs_view.insert(*dof);
        }
        for (DofEntity_multiIndex_uid_view::iterator dit = dofs_view.begin();
             dit != dofs_view.end(); dit++) {
          my_split << **dit << endl;
        }
      }
      {
        DofEntity_multiIndex_uid_view dofs_view;
        for (_IT_GET_DOFS_FIELD_BY_NAME_FOR_LOOP_(m_field2, "FIELD_B", dof)) {
          dofs_view.insert(*dof);
        }
        for (DofEntity_multiIndex_uid_view::iterator dit = dofs_view.begin();
             dit != dofs_view.end(); dit++) {
          my_split << **dit << endl;
        }
      }
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
