/** \file find_local_coordinates
 * \example find_local_coordinates
 * 
 * \brief testing finding local coordinates on tetrahedron
 *
 * \ingroup mesh_cut
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

static char help[] = "testing mesh cut test\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {

    auto test_tet = [&]() {
      MoFEMFunctionBegin;
      MatrixDouble elem_coords(4, 3);

      elem_coords(0, 0) = -1;
      elem_coords(0, 1) = -1;
      elem_coords(0, 2) = -1;
      elem_coords(1, 0) = 2;
      elem_coords(1, 1) = 0;
      elem_coords(1, 2) = 0;
      elem_coords(2, 0) = 0;
      elem_coords(2, 1) = 1;
      elem_coords(2, 2) = 0;
      elem_coords(3, 0) = 0;
      elem_coords(3, 1) = 0;
      elem_coords(3, 2) = 1;

      MatrixDouble init_local_coords(5, 3);
      init_local_coords(0, 0) = 0;
      init_local_coords(0, 1) = 0;
      init_local_coords(0, 2) = 0;
      init_local_coords(1, 0) = 0.5;
      init_local_coords(1, 1) = 0;
      init_local_coords(1, 2) = 0;
      init_local_coords(2, 0) = 0;
      init_local_coords(2, 1) = 0.5;
      init_local_coords(2, 2) = 0;
      init_local_coords(3, 0) = 0;
      init_local_coords(3, 1) = 0;
      init_local_coords(3, 2) = 0.5;
      init_local_coords(4, 0) = 1. / 3.;
      init_local_coords(4, 1) = 1. / 3.;
      init_local_coords(4, 2) = 1. / 3.;

      MatrixDouble shape(init_local_coords.size1(), 4);
      CHKERR Tools::shapeFunMBTET<3>(&shape(0, 0), &init_local_coords(0, 0),
                                     &init_local_coords(0, 1),
                                     &init_local_coords(0, 2), 5);

      MatrixDouble global_coords = prod(shape, elem_coords);

      MatrixDouble local_coords(init_local_coords.size1(), 3);
      CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
          &elem_coords(0, 0), &global_coords(0, 0), init_local_coords.size1(),
          &local_coords(0, 0));

      MatrixDouble residual = local_coords - init_local_coords;
      std::cout << residual << std::endl;
      for (auto v : residual.data())
        if (std::abs(v) > 1e-12)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Should be zer, but is v = %3.4e", v);

      MoFEMFunctionReturn(0);
    };

    CHKERR test_tet();
    
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
