/** \file petsc_smart_ptr_objects.cpp
 * \example petsc_smart_ptr_objects.cpp
 * \brief Create and destroy PETSc objects with smart pointer
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