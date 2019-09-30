/** \file l2_base_on_prism.cpp
  \example l2_base_on_prism.cpp

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

static double sum_matrix(MatrixDouble &m) {
  double s = 0;
  for (unsigned int ii = 0; ii < m.size1(); ii++) {
    for (unsigned int jj = 0; jj < m.size2(); jj++) {
      s += m(ii, jj);
    }
  }
  return s;
}

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);
  try {
    moab::Core mb_instance;
    moab::Interface &moab = mb_instance;
    ParallelComm *pcomm = ParallelComm::get_pcomm(&moab, MYPCOMM_INDEX);
    if (pcomm == NULL)
      pcomm = new ParallelComm(&moab, PETSC_COMM_WORLD);
    MoFEM::Core core(moab);
    MoFEM::Interface &m_field = core;

    PetscBool flg = PETSC_TRUE;
    char mesh_file_name[255];
#if PETSC_VERSION_GE(3, 6, 4)
    CHKERR PetscOptionsGetString(PETSC_NULL, "", "-my_file", mesh_file_name,
                                 255, &flg);
#else
    CHKERR PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-my_file",
                                 mesh_file_name, 255, &flg);
#endif
    if (flg != PETSC_TRUE) {
      SETERRQ(PETSC_COMM_SELF, 1, "*** ERROR -my_file (MESH FILE NEEDED)");
    }

    const char *option;
    option = ""; //"PARALLEL=BCAST;";//;DEBUG_IO";
    CHKERR moab.load_file(mesh_file_name, 0, option);
    Range verts;
    CHKERR moab.get_entities_by_type(0, MBVERTEX, verts, true);
    int nb_gauss_pts = verts.size();
    MatrixDouble coords(nb_gauss_pts, 3);
    CHKERR moab.get_coords(verts, &coords(0, 0));

    cout << coords << endl;

    cout << "verts: " << verts.size() << endl;
    verts.print();

    Range prisms;
    CHKERR moab.get_entities_by_type(0, MBPRISM, prisms, true);
    Range faces, tris;
    Skinner skinner(&moab);
    CHKERR skinner.find_skin(0, prisms, 2, faces);
    tris = faces.subset_by_type(MBTRI);

    cout << "tris: " << tris.size() << endl;
    tris.print();

    Range tris_verts;
    CHKERR moab.get_connectivity(tris, tris_verts, true);
    tris_verts.print();

    double x[3];
    Range::iterator it = tris_verts.begin();
    while (it != tris_verts.end()) {
      CHKERR moab.get_coords(&*it, 1, x);
      if (x[2] > 0.5) {
        it = tris_verts.erase(it);
      } else {
        it++;
      }
    }
    tris_verts.print();

    int nb_gauss_pts_on_faces = tris_verts.size();
    int nb_gauss_pts_through_thickness = nb_gauss_pts / nb_gauss_pts_on_faces;

    MatrixDouble gauss_pts_triangles_only(nb_gauss_pts_on_faces, 3);
    CHKERR moab.get_coords(tris_verts, &gauss_pts_triangles_only(0, 0));
    gauss_pts_triangles_only = trans(gauss_pts_triangles_only);
    MatrixDouble triangles_only_N(nb_gauss_pts_on_faces, 3);
    CHKERR ShapeMBTRI(&*triangles_only_N.data().begin(),
                      &gauss_pts_triangles_only(0, 0),
                      &gauss_pts_triangles_only(1, 0), nb_gauss_pts_on_faces);

    MatrixDouble gauss_pts_through_thickness(2, nb_gauss_pts_through_thickness);
    gauss_pts_through_thickness.clear();
    double step = 1.0 / (nb_gauss_pts_through_thickness - 1);
    for (int ss = 0; ss < nb_gauss_pts_through_thickness; ss++) {
      gauss_pts_through_thickness(0, ss) = step * ss;
    }

    MatrixDouble N(nb_gauss_pts, 6);
    MatrixDouble diffN(nb_gauss_pts, 18);

    MatrixDouble triangles_only_diffN(1, 6);
    CHKERR ShapeDiffMBTRI(&*triangles_only_diffN.data().begin());

    // MatrixDouble gauss_pts(3, nb_gauss_pts);
    Range gauss_pts;
    Range::iterator gp_it;
    // Calculate "nobase" base functions on prism, this is cartesian product
    // of base functions on triangles with base functions through thickness
    for (int dd = 0; dd != 6; dd++) {
      gauss_pts.clear();
      gp_it = gauss_pts.begin();
      cout << "Base function " << dd << endl;
      int gg = 0;
      for (int ggf = 0; ggf < nb_gauss_pts_on_faces; ggf++) {
        int ddd = dd > 2 ? dd - 3 : dd;
        double tri_n = triangles_only_N(ggf, ddd);
        double dksi_tri_n = triangles_only_diffN(0, 2 * ddd + 0);
        double deta_tri_n = triangles_only_diffN(0, 2 * ddd + 1);
        for (int ggt = 0; ggt < nb_gauss_pts_through_thickness; ggt++, gg++) {
          double zeta = gauss_pts_through_thickness(0, ggt);
          double dzeta, edge_shape;
          if (dd < 3) {
            dzeta = diffN_MBEDGE0;
            edge_shape = N_MBEDGE0(zeta);
          } else {
            dzeta = diffN_MBEDGE1;
            edge_shape = N_MBEDGE1(zeta);
          }
          N(gg, dd) = tri_n * edge_shape;
          diffN(gg, 3 * dd + 0) = dksi_tri_n * edge_shape;
          diffN(gg, 3 * dd + 1) = deta_tri_n * edge_shape;
          diffN(gg, 3 * dd + 2) = tri_n * dzeta;

          double x = gauss_pts_triangles_only(0, ggf);
          double y = gauss_pts_triangles_only(1, ggf);
          double z = gauss_pts_through_thickness(0, ggt);

          Range::iterator verts_it = verts.begin();
          int ii;
          for (ii = 0; ii < nb_gauss_pts; ii++) {
            double eps = 1e-12;
            if (fabs(coords(ii, 0) - x) < eps &&
                fabs(coords(ii, 1) - y) < eps &&
                fabs(coords(ii, 2) - z) < eps) {
              verts_it += ii;
              gauss_pts.insert(gp_it, *verts_it);
              break;
            }
          }

          if (ii == nb_gauss_pts) {
            SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                    "GP data inconsistency");
          } 

          // cout << x << " " << y << " " << z << endl;
        }
      }
      gauss_pts.print();
    }

    int order_triangles_only = 1;
    int order_thickness = 1;

    int P = NBVOLUMEPRISM_L2(order_triangles_only, order_thickness);

    MatrixDouble prism_base(verts.size(), P);
    MatrixDouble prism_base_diff(verts.size(), P * 3);

    CHKERR L2_Ainsworth_ShapeFunctions_MBPRISM(
        order_triangles_only, order_thickness, &*N.data().begin(),
        &*diffN.data().begin(), &*prism_base.data().begin(),
        &*prism_base_diff.data().begin(), nb_gauss_pts, Legendre_polynomials);

    MatrixDouble trans_base = trans(prism_base);
    MatrixDouble trans_base_diff = trans(prism_base_diff);
    MatrixDouble trans_n = trans(N);
    for (int rr = 0; rr < prism_base.size2(); rr++) {
      Tag th_base;
      // Tag th_base_diff;
      double def_val[] = {0};
      CHKERR moab.tag_get_handle(
          ("base_" + boost::lexical_cast<std::string>(rr)).c_str(), 1,
          MB_TYPE_DOUBLE, th_base, MB_TAG_CREAT | MB_TAG_DENSE, def_val);
      CHKERR moab.tag_set_data(th_base, gauss_pts, &trans_base(rr, 0));
      if (rr < 6) {
        Tag th_N;
        CHKERR moab.tag_get_handle(
            ("N_" + boost::lexical_cast<std::string>(rr)).c_str(), 1,
            MB_TYPE_DOUBLE, th_N, MB_TAG_CREAT | MB_TAG_DENSE, def_val);
        CHKERR moab.tag_set_data(th_N, gauss_pts, &trans_n(rr, 0));
      }
      // for (int dd = 0; dd < 3; dd++) {
      //   CHKERR moab.tag_get_handle(("base_diff_" +
      //                               boost::lexical_cast<std::string>(rr) +
      //                               "_" +
      //                               boost::lexical_cast<std::string>(dd))
      //                                  .c_str(),
      //                              1, MB_TYPE_DOUBLE, th_base_diff,
      //                              MB_TAG_CREAT | MB_TAG_DENSE, def_val);
      //   CHKERR moab.tag_set_data(th_base_diff, verts,
      //                            &trans_base_diff(3 * rr + dd, 0));
      // }
    }

    CHKERR moab.write_file("l2_base_on_prism.h5m");
  }
  CATCH_ERRORS;

  CHKERR MoFEM::Core::Finalize();
}
