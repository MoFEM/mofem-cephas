/**
 * \file bernstein_bezier_genrate_base.cpp
 * \example bernstein_bezier_genrate_base.cpp
 *
 * Genarte Bernstein-Bezier base
 *
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

    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    // Create MoFEM instance
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    const int N = 4;

    // Edge

    VectorInt edge_alpha(2 * (2 + NBEDGE_H1(N)));
    CHKERR BernsteinBezier::generateIndicesVertexEdge(N, &edge_alpha[0]);
    CHKERR BernsteinBezier::generateIndicesEdgeEdge(N, &edge_alpha[4]);
    cerr << edge_alpha << endl;

    std::array<double, 6> edge_x_k = {0, 0, 0, 0, 1, 0};
    MatrixDouble edge_x_alpha(edge_alpha.size() / 2, 3);
    CHKERR BernsteinBezier::domainPoints3d(
        N, 2, edge_alpha.size() / 2, &*edge_alpha.data().begin(),
        edge_x_k.data(), &*edge_x_alpha.data().begin());
    cerr << edge_x_alpha << endl;
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}