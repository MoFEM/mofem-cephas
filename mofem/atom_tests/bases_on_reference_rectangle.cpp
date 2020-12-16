/**
 * \file bases_on_reference_rectangle.cpp
 * \example bases_on_reference_rectangle.cpp
 *
 * Print bases on reference rectangle
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

    const char *mesh_file_name = "reference_rectangle.h5m";
    const char *option = "";
    CHKERR moab.load_file(mesh_file_name, 0, option);

    std::vector<EntityHandle> verts;
    CHKERR moab.get_entities_by_type(0, MBVERTEX, verts);
    std::vector<double> coords(3 * verts.size());
    CHKERR moab.get_coords(&*verts.begin(), verts.size(), &*coords.begin());
    MatrixDouble base(verts.size(), 4);
    MatrixDouble diff_base(verts.size(), 4 * 2);

    const int nb_gauss_pts = verts.size();

    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      const double ksi = coords[3 * gg + 0];
      const double zeta = coords[3 * gg + 1];
      base(gg, 0) = N_MBQUAD0(ksi, zeta);
      base(gg, 1) = N_MBQUAD1(ksi, zeta);
      base(gg, 2) = N_MBQUAD2(ksi, zeta);
      base(gg, 3) = N_MBQUAD3(ksi, zeta);
      diff_base(gg, 0) = diffN_MBQUAD0x(zeta);
      diff_base(gg, 1) = diffN_MBQUAD0y(ksi);
      diff_base(gg, 2) = diffN_MBQUAD1x(zeta);
      diff_base(gg, 3) = diffN_MBQUAD1y(ksi);
      diff_base(gg, 4) = diffN_MBQUAD2x(zeta);
      diff_base(gg, 5) = diffN_MBQUAD2y(ksi);
      diff_base(gg, 6) = diffN_MBQUAD3x(zeta);
      diff_base(gg, 7) = diffN_MBQUAD3y(ksi);
    }

    const int p = 5;
    MatrixDouble face_family(
        2, 3 * NBFACEQUAD_DEMKOWICZ_FAMILY_QUAD_HCURL(p, p) * nb_gauss_pts);
    MatrixDouble diff_face_family(
        2, 3 * 2 * NBFACEQUAD_DEMKOWICZ_FAMILY_QUAD_HCURL(p, p) * nb_gauss_pts);

    int order[2] = {p, p};
    double *face_family_ptr[] = {&face_family(0, 0), &face_family(1, 0)};
    double *diff_face_family_ptr[] = {&diff_face_family(0, 0),
                                      &diff_face_family(1, 0)};
    int face_nodes[] = {0, 1, 2, 3};
    CHKERR DemkowiczHexAndQuad::Hcurl_FaceShapeFunctions_ONQUAD(
        face_nodes, order, &*base.data().begin(), &*diff_base.data().begin(),
        face_family_ptr, diff_face_family_ptr, nb_gauss_pts);

    int nb_dofs = NBFACEQUAD_DEMKOWICZ_HCURL(p);
    MatrixDouble face_n;
    MatrixDouble diff_face_n;
    face_n.resize(nb_gauss_pts, 3 * nb_dofs, false);
    diff_face_n.resize(nb_gauss_pts, 3 * 2 * nb_dofs, false);

    double *ptr_f0 = &face_family(0, 0);
    double *ptr_f1 = &face_family(1, 0);
    double *ptr = &face_n(0, 0);
    for (int n = 0; n != face_family.size2() / 3; ++n) {
      for (int j = 0; j != 3; ++j) {
        *ptr = *ptr_f0;
        ++ptr;
        ++ptr_f0;
      }
      for (int j = 0; j != 3; ++j) {
        *ptr = *ptr_f1;
        ++ptr;
        ++ptr_f1;
      }
    }
    double *diff_ptr_f0 = &diff_face_family(0, 0);
    double *diff_ptr_f1 = &diff_face_family(1, 0);
    double *diff_ptr = &diff_face_n(0, 0);
    for (int n = 0; n != diff_face_family.size2() / 6; ++n) {
      for (int j = 0; j != 6; ++j) {
        *diff_ptr = *diff_ptr_f0;
        ++diff_ptr;
        ++diff_ptr_f0;
      }
      for (int j = 0; j != 6; ++j) {
        *diff_ptr = *diff_ptr_f1;
        ++diff_ptr;
        ++diff_ptr_f1;
      }
    }

    auto create_tag = [&](auto name) {
      Tag th;
      double def_val[] = {0, 0, 0};
      CHKERR moab.tag_get_handle(name, 3, MB_TYPE_DOUBLE, th,
                                 MB_TAG_CREAT | MB_TAG_SPARSE, &def_val);
      return th;
    };

    auto th = create_tag("BASE_BUBBLE");

    int nb_base_bubble = 0;
    CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-nb_base_bubble",
                              &nb_base_bubble, PETSC_NULL);
    for (int gg = 0; gg != nb_gauss_pts; ++gg) {
      CHKERR moab.tag_set_data(th, &verts[gg], 1,
                               &face_n(gg, 3 * nb_base_bubble));
    }

    int sense[4];
    double *hcurl_edge_n[4];
    double *diff_hcurl_edge_n[4];
    MatrixDouble face_edge_n[4];
    MatrixDouble diff_face_edge_n[4];

    for (int ee = 0; ee != 4; ++ee) {

      sense[ee] = 1;
      const int nb_dofs = NBEDGE_DEMKOWICZ_HCURL(p);
      face_edge_n[ee].resize(nb_gauss_pts, 3 * nb_dofs, false);
      diff_face_edge_n[ee].resize(nb_gauss_pts, 3 * 2 * nb_dofs, false);
      hcurl_edge_n[ee] = &*face_edge_n[ee].data().begin();
      diff_hcurl_edge_n[ee] = &*diff_face_edge_n[ee].data().begin();

    }
    int pp[4] = {p, p, p, p};
    CHKERR DemkowiczHexAndQuad::Hcurl_EdgeShapeFunctions_ONQUAD(
        sense, pp, &*base.data().begin(), &*diff_base.data().begin(),
        hcurl_edge_n, diff_hcurl_edge_n, nb_gauss_pts);

    for (int ee = 0; ee != 4; ++ee) {

      auto th = create_tag(
          ("BASE_BUBBLE" + boost::lexical_cast<std::string>(ee)).c_str());

      int nb_base_edge = 0;
      CHKERR PetscOptionsGetInt(PETSC_NULL, "", "-nb_base_edge",
                                &nb_base_edge, PETSC_NULL);
      for (int gg = 0; gg != nb_gauss_pts; ++gg) {
        CHKERR moab.tag_set_data(th, &verts[gg], 1,
                                 &face_edge_n[ee](gg, 3 * nb_base_edge));
      }
    }

    CHKERR moab.write_file("out.vtk", "VTK", "");
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}