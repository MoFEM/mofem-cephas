/** \file find_local_coordinates
 * \example find_local_coordinates
 * 
 * \brief testing finding local coordinates on tetrahedron
 *
 * \ingroup mesh_cut
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
      MOFEM_LOG("SELF", Sev::inform) << residual;
      for (auto v : residual.data())
        if (std::abs(v) > 1e-12)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Should be zer, but is v = %3.4e", v);

      MoFEMFunctionReturn(0);
    };

    auto test_tri = [&]() {
      MoFEMFunctionBegin;
      MatrixDouble elem_coords(3, 3);

      elem_coords(0, 0) = -1;
      elem_coords(0, 1) = -1;
      elem_coords(0, 2) = -1;
      elem_coords(1, 0) = 2;
      elem_coords(1, 1) = 0;
      elem_coords(1, 2) = 0;
      elem_coords(2, 0) = 0;
      elem_coords(2, 1) = 1;
      elem_coords(2, 2) = 0;

      MatrixDouble init_local_coords(4, 2);
      init_local_coords(0, 0) = 0;
      init_local_coords(0, 1) = 0;
      init_local_coords(1, 0) = 1;
      init_local_coords(1, 1) = 0;
      init_local_coords(2, 0) = 0;
      init_local_coords(2, 1) = 1;
      init_local_coords(3, 0) = 1. / 3.;
      init_local_coords(3, 1) = 1. / 3.;

      MatrixDouble shape(init_local_coords.size1(), 3);
      CHKERR Tools::shapeFunMBTRI<2>(&shape(0, 0), &init_local_coords(0, 0),
                                     &init_local_coords(0, 1),
                                     init_local_coords.size1());
      MOFEM_LOG("SELF", Sev::verbose) << "tri shape " << shape;

      MatrixDouble global_coords = prod(shape, elem_coords);
      MOFEM_LOG("SELF", Sev::verbose) << "tri global_coords " << global_coords;

      MatrixDouble local_coords(init_local_coords.size1(), 2);
      CHKERR Tools::getLocalCoordinatesOnReferenceTriNodeTri(
          &elem_coords(0, 0), &global_coords(0, 0), init_local_coords.size1(),
          &local_coords(0, 0));

      MatrixDouble residual = local_coords - init_local_coords;
      MOFEM_LOG("SELF", Sev::inform) << residual;
      for (auto v : residual.data())
        if (std::abs(v) > 1e-12)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Should be zer, but is v = %3.4e", v);

      MoFEMFunctionReturn(0);
    };

    auto test_edge = [&]() {
      MoFEMFunctionBegin;
      MatrixDouble elem_coords(2, 3);

      elem_coords(0, 0) = -1;
      elem_coords(0, 1) = -1;
      elem_coords(0, 2) = -1;
      elem_coords(1, 0) = 2;
      elem_coords(1, 1) = 0;
      elem_coords(1, 2) = 0;

      VectorDouble init_local_coords(6);
      init_local_coords[0] = 0;
      init_local_coords[1] = 1;
      init_local_coords[2] = 0.5;
      init_local_coords[4] = 0.25;
      init_local_coords[5] = 0.75;

      MatrixDouble shape(init_local_coords.size(), 2);
      CHKERR Tools::shapeFunMBEDGE<1>(&shape(0, 0), &init_local_coords[0],
                                      init_local_coords.size());
      MOFEM_LOG("SELF", Sev::verbose) << "edge shape " << shape;

      MatrixDouble global_coords = prod(shape, elem_coords);
      MOFEM_LOG("SELF", Sev::verbose) << "edge global_coords " << global_coords;

      VectorDouble local_coords(init_local_coords.size());
      CHKERR Tools::getLocalCoordinatesOnReferenceEdgeNodeEdge(
          &elem_coords(0, 0), &global_coords(0, 0), init_local_coords.size(),
          &local_coords[0]);

      VectorDouble residual = local_coords - init_local_coords;
      MOFEM_LOG("SELF", Sev::inform) << residual;
      for (auto v : residual.data())
        if (std::abs(v) > 1e-12)
          SETERRQ1(PETSC_COMM_SELF, MOFEM_ATOM_TEST_INVALID,
                   "Should be zer, but is v = %3.4e", v);

      MoFEMFunctionReturn(0);
    };

    CHKERR test_tet();
    CHKERR test_tri();
    CHKERR test_edge();
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}
