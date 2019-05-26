/** \file petsc_smart_ptr_objects.cpp
 * \example petsc_smart_ptr_objects.cpp
 * \brief Create and destroy PETSc objects with smart pointer
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

#define DEBG_PETSC_OBJ_INTRUSIVE

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    SmartPetscObj<Mat> m;
    auto check = [&](const int expected) {
      MoFEMFunctionBegin;
      if (m.use_count() != expected)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "use count should be %d but is %d", expected, m.use_count());
      MoFEMFunctionReturn(0);
    };

    CHKERR MatCreate(PETSC_COMM_SELF, m.get());
    CHKERR check(1);
    CHKERR MatSetSizes(m, 2, 2, 2, 2);

    {
      boost::intrusive_ptr<Mat> n = m;
      CHKERR check(2);
    }
    {
      boost::intrusive_ptr<Mat> n = m;
      CHKERR check(2);
    }

    CHKERR check(1);

    m.reset();
    CHKERR check(0);

  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}