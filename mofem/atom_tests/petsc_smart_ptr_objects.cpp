/** \file petsc_smart_ptr_objects.cpp
 * \example petsc_smart_ptr_objects.cpp
 * \brief Create and destroy PETSc objects with smart pointer
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