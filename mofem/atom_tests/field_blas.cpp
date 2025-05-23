/** \file field_blas.cpp
 * \example field_blas.cpp
 * \brief test field blas interface
 *
 * \ingroup mofem_field_algebra
 */



#include <MoFEM.hpp>

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
    int order = 1;
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-my_order", &order, PETSC_NULL);
#else
    CHKERR PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-my_order", &order,
                              PETSC_NULL);
#endif

    // Read mesh to MOAB
    const char *option;
    option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    // Create MoFEM database
    // Note: Is MoFEM::CoreTmp<1> for testing purposes only
    MoFEM::CoreTmp<-1> core(moab);
    MoFEM::Interface &m_field = core;

    CHKERR m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, BitRefLevel().set(0));

    // Create fields, add entities to field and set approximation order
    CHKERR m_field.add_field("FIELD_A", H1, AINSWORTH_LEGENDRE_BASE, 3,
                             MB_TAG_DENSE);
    CHKERR m_field.add_field("FIELD_B", H1, AINSWORTH_LEGENDRE_BASE, 3,
                             MB_TAG_DENSE);
    CHKERR m_field.add_field("FIELD_C", H1, AINSWORTH_LEGENDRE_BASE, 1,
                             MB_TAG_DENSE);

    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_A");
    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_B");
    CHKERR m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_C");

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_A", order + 1);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_A", order + 1);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_A", order + 1);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_A", 1);

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_B", order);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_B", 1);

    CHKERR m_field.set_field_order(0, MBTET, "FIELD_C", order);
    CHKERR m_field.set_field_order(0, MBTRI, "FIELD_C", order);
    CHKERR m_field.set_field_order(0, MBEDGE, "FIELD_C", order);
    CHKERR m_field.set_field_order(0, MBVERTEX, "FIELD_C", 1);

    // build field
    CHKERR m_field.build_fields();

    // get access to field agebra interface
    auto fb = m_field.getInterface<FieldBlas>();

    auto check_axpy = [&]() {
      MoFEMFunctionBegin;
      // set value to field
      CHKERR fb->setField(+1, MBVERTEX, "FIELD_A");
      CHKERR fb->setField(-2, MBVERTEX, "FIELD_B");

      // FIELD_A = FIELD_A + 0.5 * FIELD_B
      CHKERR fb->fieldAxpy(+0.5, "FIELD_B", "FIELD_A");
      // FIELD_B *= -0.5
      CHKERR fb->fieldScale(-0.5, "FIELD_B");

      auto dofs_ptr = m_field.get_dofs();
      for (auto dit : *dofs_ptr) {

        auto check = [&](const std::string name, const double expected) {
          MoFEMFunctionBegin;
          if (dit->getName() == name) {
            cout << name << " " << dit->getFieldData() << " " << expected
                 << endl;
            if (dit->getFieldData() != expected)
              SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "Wrong DOF value 0 != %4.3e for %s", dit->getFieldData(),
                       boost::lexical_cast<std::string>(*dit).c_str());
          }

          MoFEMFunctionReturn(0);
        };

        CHKERR check("FIELD_A", 0);
        CHKERR check("FIELD_B", 1);
      }
      MoFEMFunctionReturn(0);
    };

    auto check_set_vertex = [&]() {
      MoFEMFunctionBegin;

      auto set_distance = [&](VectorAdaptor &&field_data, double *xcoord,
                              double *ycoord, double *zcoord) {
        MoFEMFunctionBegin;
        FTensor::Index<'i', 3> i;
        FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_coord(
            xcoord, ycoord, zcoord);
        field_data[0] = sqrt(t_coord(i) * t_coord(i));
        MoFEMFunctionReturn(0);
      };

      CHKERR fb->setVertexDofs(set_distance, "FIELD_A");

      // Test set veryex
      struct TestMethod : EntityMethod {
        TestMethod(moab::Interface &moab) : EntityMethod(), moab(moab) {}

        MoFEMErrorCode preProcess() { return 0; }
        MoFEMErrorCode operator()() {
          MoFEMFunctionBegin;
          if (entPtr->getEntType() == MBVERTEX) {
            EntityHandle v = entPtr->getEnt();
            VectorDouble3 coords(3);
            CHKERR moab.get_coords(&v, 1, &coords[0]);
            if (std::abs(entPtr->getEntFieldData()[0] - norm_2(coords)) >
                std::numeric_limits<double>::epsilon())
              SETERRQ2(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                       "Wrong vals %3.4e != %3.4e",
                       entPtr->getEntFieldData()[0], norm_2(coords));
          }
          MoFEMFunctionReturn(0);
        }
        MoFEMErrorCode postProcess() { return 0; }

      private:
        moab::Interface &moab;
      };

      TestMethod method(moab);
      // checking if all is ok
      CHKERR m_field.loop_entities("FIELD_A", method);

      MoFEMFunctionReturn(0);
    };

    auto check_lambda = [&]() {
      MoFEMFunctionBegin;

      Range meshset_ents;
      for (_IT_CUBITMESHSETS_BY_BCDATA_TYPE_FOR_LOOP_(
               m_field, SIDESET | PRESSURESET, it)) {
        if (it->getMeshsetId() != 1)
          continue;
        Range ents, nodes;
        CHKERR it->getMeshsetIdEntitiesByType(m_field.get_moab(), MBTRI, ents,
                                              true);
        CHKERR m_field.get_moab().get_connectivity(ents, nodes, true);
        meshset_ents.merge(nodes);
      }
      // set value to field
      CHKERR fb->setField(+1, MBVERTEX, "FIELD_A");
      CHKERR fb->setField(-1, MBVERTEX, "FIELD_B");
      CHKERR fb->setField(1.5, MBVERTEX, "FIELD_C");

      auto field_axpy = [&](const double val_y, const double val_x) {
        return 0.5 * val_x + val_y;
      };
      auto field_scale = [&](const double val) { return -0.5 * val; };

      // note that y_ent is first
      // FIXME: this can be confusing?
      auto vector_times_scalar_field =
          [&](boost::shared_ptr<FieldEntity> ent_ptr_y,
              boost::shared_ptr<FieldEntity> ent_ptr_x) {
            MoFEMFunctionBeginHot;
            auto x_data = ent_ptr_x->getEntFieldData();
            auto y_data = ent_ptr_y->getEntFieldData();
            // const auto size_x = x_data.size(); // scalar
            const auto size_y = y_data.size(); // vector

            for (size_t dd = 0; dd != size_y; ++dd)
              y_data[dd] *= x_data[0];

            // y_data *= x_data[0]; //would work as well

            MoFEMFunctionReturnHot(0);
          };

      Range ents; // TODO: create test for subentities
      // FIELD_A = FIELD_A + 0.5 * FIELD_B
      CHKERR fb->fieldLambdaOnValues(field_axpy, "FIELD_B", "FIELD_A", true);
      // CHKERR fb->fieldAxpy(+0.5, "FIELD_B", "FIELD_A");

      // FIELD_B *= -0.5
      // CHKERR fb->fieldScale(-0.5, "FIELD_B");
      CHKERR fb->fieldLambdaOnValues(field_scale, "FIELD_B", &meshset_ents);

      // FIELD_B(i) *= FIELD_C
      CHKERR fb->fieldLambdaOnEntities(vector_times_scalar_field, "FIELD_C",
                                       "FIELD_B", true, &meshset_ents);

      auto dofs_ptr = m_field.get_dofs();
      constexpr double eps = std::numeric_limits<double>::epsilon();
      for (auto dit : *dofs_ptr) {
        auto check = [&](const std::string name, const double expected) {
          MoFEMFunctionBegin;
          if (dit->getName() == name && dit->getEntType() == MBVERTEX) {
            cout << name << " " << dit->getFieldData() << " " << expected
                 << endl;
            if (abs(dit->getFieldData() - expected) > eps)
              SETERRQ3(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                       "Wrong DOF value %4.3e != %4.3e for %s",
                       dit->getFieldData(), expected,
                       boost::lexical_cast<std::string>(*dit).c_str());
          }

          MoFEMFunctionReturn(0);
        };

        CHKERR check("FIELD_A", 0.5);

        if (meshset_ents.find(dit->getEnt()) != meshset_ents.end()) {
          CHKERR check("FIELD_B", 0.75);
          CHKERR check("FIELD_C", 1.5);
        }
      }

      MoFEMFunctionReturn(0);
    };

    CHKERR check_axpy();
    CHKERR check_set_vertex();
    CHKERR check_lambda();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
