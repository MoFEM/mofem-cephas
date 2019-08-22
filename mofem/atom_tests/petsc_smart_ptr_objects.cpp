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

#include <MoFEM.hpp>

using namespace MoFEM;

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    SmartPetscObj<Mat> m_ptr;

    // check is count is correct
    auto check = [&](const int expected) {
      MoFEMFunctionBegin;
      if (m_ptr.use_count() != expected)
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "use count should be %d but is %d", expected,
                 m_ptr.use_count());
      MoFEMFunctionReturn(0);
    };

    // check copy constructor
    {
      Mat m;
      CHKERR MatCreate(PETSC_COMM_SELF, &m);
      m_ptr = SmartPetscObj<Mat>(m);
    }

    // Check casting on PetscObject
    { PetscObject obj = static_cast<PetscObject>(m_ptr); }

    { 
      SmartPetscObj<Mat> n_ptr(m_ptr); 
      CHKERR check(2);
    }

    CHKERR check(1);

    // check if casting works well
    CHKERR MatSetSizes(m_ptr, 2, 2, 2, 2);

    {
      // nesting once
      SmartPetscObj<Mat> n_ptr = m_ptr;
      CHKERR check(2);
      {
        // again
        SmartPetscObj<Mat> i_ptr = m_ptr;
        CHKERR check(3);
      }
       CHKERR check(2);
    }

    
    {
      // nesting again
      SmartPetscObj<Mat> n_ptr = m_ptr;
      CHKERR check(2);
    }

    // only one now
    CHKERR check(1);

    // reset, i.e. delete matrix
    m_ptr.reset();

    // counts should be zero now
    CHKERR check(0);

    {
      SmartPetscObj<Mat> m;
    }

  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}