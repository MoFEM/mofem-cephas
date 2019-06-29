/** \file field_blas_set_vertex_dofs.cpp
 * \example field_blas_set_vertex_dofs.cpp
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

  // initialize petsc
  MoFEM::Core::Initialize(&argc, &argv, nullptr, help);

  try {

    // Create moab and mofem instances
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // Register DM Manager
    DMType dm_name = "DMMOFEM";
    CHKERR DMRegister_MoFEM(dm_name);

    // Simple interface
    Simple *simple_interface;
    CHKERR m_field.getInterface(simple_interface);
    CHKERR simple_interface->getOptions();
    CHKERR simple_interface->loadFile("");

    // add fields
    CHKERR simple_interface->addDomainField("u", H1, AINSWORTH_LEGENDRE_BASE,
                                            1);
    // set fields order
    CHKERR simple_interface->setFieldOrder("u", 1);

    // build fields (finite elements and problem not needed)
    CHKERR simple_interface->buildFields();

    FieldBlas *field_blas;
    CHKERR m_field.getInterface(field_blas);

    auto set_distance = [&](VectorAdaptor &field_data, double *xcoord,
                            double *ycoord, double *zcoord) {
      MoFEMFunctionBegin;
      FTensor::Index<'i', 3> i;
      FTensor::Tensor1<FTensor::PackPtr<double *, 0>, 3> t_coord(xcoord, ycoord,
                                                                 zcoord);
      field_data[0] = sqrt(t_coord(i) * t_coord(i));
      MoFEMFunctionReturn(0);
    };
    CHKERR field_blas->setVertexDofs(set_distance, "u");
    

    // This is only used for testing
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
                     "Wrong vals %3.4e != %3.4e", entPtr->getEntFieldData()[0],
                     norm_2(coords));
        }
        MoFEMFunctionReturn(0);
      }
      MoFEMErrorCode postProcess() { return 0; }

    private:
      moab::Interface &moab;
    };

    TestMethod method(moab);
    // checking if all is ok
    CHKERR m_field.loop_entities("u", method);
  }
  CATCH_ERRORS;

  // finish work cleaning memory, getting statistics, etc.
  MoFEM::Core::Finalize();
}