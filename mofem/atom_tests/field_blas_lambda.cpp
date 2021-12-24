/** \file field_blas_lambda.cpp
 * \example field_blas_lambda.cpp
 * \brief test field blas interface
 *
 * \ingroup mofem_field_algebra
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
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-file_name", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-file_name",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -file_name (MESH FILE NEEDED)");
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
          const auto size_x = x_data.size(); // scalar
          const auto size_y = y_data.size(); // vector

          for (size_t dd = 0; dd != size_y; ++dd)
            y_data[dd] *= x_data[0];

          MoFEMFunctionReturnHot(0);
        };

    Range ents; // TODO: create test for subentities
    // FIELD_A = FIELD_A + 0.5 * FIELD_B
    CHKERR fb->fieldLambdaOnValues(field_axpy, "FIELD_B", "FIELD_A");
    // CHKERR fb->fieldAxpy(+0.5, "FIELD_B", "FIELD_A");

    // FIELD_B *= -0.5
    // CHKERR fb->fieldScale(-0.5, "FIELD_B");
    CHKERR fb->fieldLambdaOnValues(field_scale, "FIELD_B");

    // FIELD_B(i) *= FIELD_C
    CHKERR fb->fieldLambdaOnEntities(vector_times_scalar_field, "FIELD_C",
                                     "FIELD_B");

    auto dofs_ptr = m_field.get_dofs();
    constexpr double eps = std::numeric_limits<double>::epsilon();
    for (auto dit : *dofs_ptr) {

      auto check = [&](const std::string name, const double expected) {
        MoFEMFunctionBegin;
        if (dit->getName() == name) {
          cout << name << " " << dit->getFieldData() << " " << expected << endl;
          if (dit->getFieldData() - expected > eps)
            SETERRQ2(PETSC_COMM_WORLD, MOFEM_ATOM_TEST_INVALID,
                     "Wrong DOF value 0 != %4.3e for %s", dit->getFieldData(),
                     boost::lexical_cast<std::string>(*dit).c_str());
        }

        MoFEMFunctionReturn(0);
      };

      CHKERR check("FIELD_A", 0.5);
      CHKERR check("FIELD_B", 0.75);
      CHKERR check("FIELD_C", 1.5);
    }
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
