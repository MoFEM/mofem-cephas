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
#define RND_EPS 1e-6

// Rounding
double roundn(double n) {
  // break n into fractional part (fract) and integral part (intp)
  double fract, intp;
  fract = modf(n, &intp);
  // case where n approximates zero, set n to "positive" zero
  if (std::abs(intp) == 0) {
    if (std::abs(fract) <= RND_EPS) {
      n = 0.000;
    }
  }
  return n;
}

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
    ierr = PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
    CHKERRG(ierr);
#else
    ierr = PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
    CHKERRG(ierr);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }
    PetscInt order;
#if PETSC_VERSION_GE(3, 6, 4)
    ierr = PetscOptionsGetInt(PETSC_NULL, "", "-my_order", &order, &flg);
    CHKERRG(ierr);
#else
    ierr =
        PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-my_order", &order, &flg);
    CHKERRG(ierr);
#endif
    if (flg != PETSC_TRUE) {
      order = 1;
    }

    // Read mesh to MOAB
    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    rval = moab.load_file(mesh_file_name, 0, option);
    CHKERRG(rval);
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);

    // Create MoFEM (Joseph) database
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    // stl::bitset see for more details
    BitRefLevel bit_level0;
    bit_level0.set(0);
    EntityHandle meshset_level0;
    rval = moab.create_meshset(MESHSET_SET, meshset_level0);
    CHKERRG(rval);
    ierr = m_field.getInterface<BitRefManager>()->setBitRefLevelByDim(
        0, 3, bit_level0);
    CHKERRG(ierr);
    ierr = m_field.getInterface<BitRefManager>()->getEntitiesByRefLevel(
        bit_level0, BitRefLevel().set(), meshset_level0);
    CHKERRG(ierr);

    /***/
    // Define problem

    // Fields
    ierr = m_field.add_field("FIELD_A", H1, AINSWORTH_LEGENDRE_BASE, 3);
    CHKERRG(ierr);
    ierr = m_field.add_field("FIELD_B", H1, AINSWORTH_LEGENDRE_BASE, 3);
    CHKERRG(ierr);

    ierr = m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_A");
    CHKERRG(ierr);
    ierr = m_field.add_ents_to_field_by_type(0, MBTET, "FIELD_B");
    CHKERRG(ierr);

    ierr = m_field.set_field_order(0, MBTET, "FIELD_A", order + 1);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBTRI, "FIELD_A", order + 1);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBEDGE, "FIELD_A", order + 1);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBVERTEX, "FIELD_A", 1);
    CHKERRG(ierr);

    ierr = m_field.set_field_order(0, MBTET, "FIELD_B", order);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBTRI, "FIELD_B", order);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBEDGE, "FIELD_B", order);
    CHKERRG(ierr);
    ierr = m_field.set_field_order(0, MBVERTEX, "FIELD_B", 1);
    CHKERRG(ierr);

    // build field
    ierr = m_field.build_fields();
    CHKERRG(ierr);

    ierr = m_field.getInterface<FieldBlas>()->setField(+1, MBVERTEX, "FIELD_A");
    CHKERRG(ierr);
    ierr = m_field.getInterface<FieldBlas>()->setField(-2, MBVERTEX, "FIELD_B");
    CHKERRG(ierr);

    ierr = m_field.getInterface<FieldBlas>()->fieldAxpy(+0.5, "FIELD_B",
                                                        "FIELD_A");
    CHKERRG(ierr);
    ierr = m_field.getInterface<FieldBlas>()->fieldScale(-0.5, "FIELD_B");
    CHKERRG(ierr);

    // Open mesh_file_name.txt for writing
    std::ofstream myfile;
    myfile.open("field_axpy_test.txt");

    const DofEntity_multiIndex *dofs_ptr;
    ierr = m_field.get_dofs(&dofs_ptr);
    CHKERRG(ierr);
    for (DofEntity_multiIndex::iterator dit = dofs_ptr->begin();
         dit != dofs_ptr->end(); dit++) {

      if ((*dit)->getEntType() != MBVERTEX)
        continue;

      if ((*dit)->getDofCoeffIdx() == 0) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
        myfile << boost::format("%.3lf") % roundn(fval) << "  ";
      }
      if ((*dit)->getDofCoeffIdx() == 1) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << "  ";
        myfile << boost::format("%.3lf") % roundn(fval) << "  ";
      }
      if ((*dit)->getDofCoeffIdx() == 2) {
        // Round and truncate to 3 decimal places
        double fval = (*dit)->getFieldData();
        std::cout << boost::format("%.3lf") % roundn(fval) << std::endl;
        myfile << boost::format("%.3lf") % roundn(fval) << std::endl;
      }
    }

    myfile.close();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();
  return 0;
}
