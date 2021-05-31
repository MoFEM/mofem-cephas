/** \file gauss_points_on_quad.cpp
  \example gauss_points_on_quad.cpp
  \brief Testing gauss points coordinates and weights for a quad face

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
    MoFEM::Core core(moab);

    MatrixDouble pts_quad;
    int rule_ksi = 6;
    int rule_eta = 8;
    CHKERR Tools::outerProductOfEdgeIntegrationPtsForQuad(pts_quad, rule_ksi,
                                                          rule_eta);
    int nb_gauss_pts = pts_quad.size2();
    double sum_coords = 0, sum_gauss_pts = 0;
    for (int i = 0; i < nb_gauss_pts; ++i) {
      sum_coords += pts_quad(0, i) + pts_quad(1, i);
      sum_gauss_pts += pts_quad(2, i);
    }
    double eps = 1e-8;
    if (fabs(20.0 - sum_coords) > eps) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }
    if (fabs(1.0 - sum_gauss_pts) > eps) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong result");
    }
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
