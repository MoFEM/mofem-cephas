/** \file segments_distance.cpp
 * \brief test segments distance
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

    MatrixDouble elem_coords(4, 3);

    elem_coords(0, 0) = -1;
    elem_coords(0, 1) = 0;
    elem_coords(0, 2) = 0;
    elem_coords(1, 0) = 2;
    elem_coords(1, 1) = 0;
    elem_coords(1, 2) = 0;
    elem_coords(2, 0) = 0;
    elem_coords(2, 1) = 1;
    elem_coords(2, 2) = 0;
    elem_coords(3, 0) = 0;
    elem_coords(3, 1) = 0;
    elem_coords(3, 2) = 1;

    cerr << elem_coords << endl;

    MatrixDouble init_local_coords(4, 3);
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

    std::cerr << "init" << endl;
    std::cerr << init_local_coords << endl;

    MatrixDouble shape(4, 4);
    CHKERR Tools::nMBTET<3>(&shape(0, 0), &init_local_coords(0, 0),
                            &init_local_coords(0, 1), &init_local_coords(0, 2),
                            4);

    std::cerr << "shape" << endl;
    std::cerr << shape << endl;

    MatrixDouble global_coords = prod(shape, elem_coords);

    std::cerr << "glob" << endl;
    std::cerr << global_coords << endl;

    MatrixDouble local_coords(4, 3);
    CHKERR Tools::getLocalCoordinatesOnReferenceFourNodeTet(
        &elem_coords(0, 0), &global_coords(0, 0), 4, &local_coords(0, 0));

    std::cerr << "res" << endl;
    std::cerr << local_coords << std::endl;
    std::cerr << local_coords - init_local_coords << std::endl;
    
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
