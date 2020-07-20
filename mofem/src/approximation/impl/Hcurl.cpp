/** \file Hcurl.cpp

  \brief Implementation of H-curl base

  Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle and by Demkowicz
  Shape functions for MBTRI/MBTET and HCurl space

*/

// #undef GENERATE_VTK_WITH_CURL_BASE
#ifndef GENERATE_VTK_WITH_CURL_BASE

using namespace MoFEM;

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_EdgeBaseFunctions_MBTET(
    int *sense, int *p, double *N, double *diffN, double *edge_n[],
    double *diff_edge_n[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBeginHot;

  const int edges_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                 {0, 3}, {1, 3}, {2, 3}};
  int P[6];
  for (int ee = 0; ee != 6; ee++)
    P[ee] = NBEDGE_AINSWORTH_HCURL(p[ee]);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};
  double edge_diff_ksi[6][3];
  FTensor::Tensor1<double *, 3> t_edge_diff_ksi[6] = {
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[0][0], &edge_diff_ksi[0][1],
                                    &edge_diff_ksi[0][2]),
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[1][0], &edge_diff_ksi[1][1],
                                    &edge_diff_ksi[1][2]),
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[2][0], &edge_diff_ksi[2][1],
                                    &edge_diff_ksi[2][2]),
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[3][0], &edge_diff_ksi[3][1],
                                    &edge_diff_ksi[3][2]),
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[4][0], &edge_diff_ksi[4][1],
                                    &edge_diff_ksi[4][2]),
      FTensor::Tensor1<double *, 3>(&edge_diff_ksi[5][0], &edge_diff_ksi[5][1],
                                    &edge_diff_ksi[5][2])};
  for (int ee = 0; ee != 6; ee++) {
    t_edge_diff_ksi[ee](i) = (t_node_diff_ksi[edges_nodes[ee][1]](i) -
                              t_node_diff_ksi[edges_nodes[ee][0]](i)) *
                             sense[ee];
  }

  FTensor::Tensor1<double *, 3> t_edge_n[6] = {
      FTensor::Tensor1<double *, 3>(&edge_n[0][0], &edge_n[0][1], &edge_n[0][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[1][0], &edge_n[1][1], &edge_n[1][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[2][0], &edge_n[2][1], &edge_n[2][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[3][0], &edge_n[3][1], &edge_n[3][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[4][0], &edge_n[4][1], &edge_n[4][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[5][0], &edge_n[5][1], &edge_n[5][2],
                                    3)};
  FTensor::Tensor2<double *, 3, 3> t_diff_edge_n[6] = {
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[0][0], &diff_edge_n[0][3], &diff_edge_n[0][6],
          &diff_edge_n[0][1], &diff_edge_n[0][4], &diff_edge_n[0][7],
          &diff_edge_n[0][2], &diff_edge_n[0][5], &diff_edge_n[0][8], 9),
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[1][0], &diff_edge_n[1][3], &diff_edge_n[1][6],
          &diff_edge_n[1][1], &diff_edge_n[1][4], &diff_edge_n[1][7],
          &diff_edge_n[1][2], &diff_edge_n[1][5], &diff_edge_n[1][8], 9),
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[2][0], &diff_edge_n[2][3], &diff_edge_n[2][6],
          &diff_edge_n[2][1], &diff_edge_n[2][4], &diff_edge_n[2][7],
          &diff_edge_n[2][2], &diff_edge_n[2][5], &diff_edge_n[2][8], 9),
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[3][0], &diff_edge_n[3][3], &diff_edge_n[3][6],
          &diff_edge_n[3][1], &diff_edge_n[3][4], &diff_edge_n[3][7],
          &diff_edge_n[3][2], &diff_edge_n[3][5], &diff_edge_n[3][8], 9),
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[4][0], &diff_edge_n[4][3], &diff_edge_n[4][6],
          &diff_edge_n[4][1], &diff_edge_n[4][4], &diff_edge_n[4][7],
          &diff_edge_n[4][2], &diff_edge_n[4][5], &diff_edge_n[4][8], 9),
      FTensor::Tensor2<double *, 3, 3>(
          &diff_edge_n[5][0], &diff_edge_n[5][3], &diff_edge_n[5][6],
          &diff_edge_n[5][1], &diff_edge_n[5][4], &diff_edge_n[5][7],
          &diff_edge_n[5][2], &diff_edge_n[5][5], &diff_edge_n[5][8], 9)};
  FTensor::Tensor1<double, 3> t_psi_e_0, t_psi_e_1;
  FTensor::Tensor2<double, 3, 3> t_diff_psi_e_0, t_diff_psi_e_1;

  for (int ii = 0; ii != nb_integration_pts; ii++) {

    const int node_shift = ii * 4;
    for (int ee = 0; ee != 6; ee++) {

      if (P[ee] == 0)
        continue;

      t_psi_e_0(i) = (N[node_shift + edges_nodes[ee][1]] *
                          t_node_diff_ksi[edges_nodes[ee][0]](i) -
                      N[node_shift + edges_nodes[ee][0]] *
                          t_node_diff_ksi[edges_nodes[ee][1]](i)) *
                     sense[ee];
      t_diff_psi_e_0(i, j) = (t_node_diff_ksi[edges_nodes[ee][1]](j) *
                                  t_node_diff_ksi[edges_nodes[ee][0]](i) -
                              t_node_diff_ksi[edges_nodes[ee][0]](j) *
                                  t_node_diff_ksi[edges_nodes[ee][1]](i)) *
                             sense[ee];

      t_psi_e_1(i) = N[node_shift + edges_nodes[ee][1]] *
                         t_node_diff_ksi[edges_nodes[ee][0]](i) +
                     N[node_shift + edges_nodes[ee][0]] *
                         t_node_diff_ksi[edges_nodes[ee][1]](i);
      t_diff_psi_e_1(i, j) = t_node_diff_ksi[edges_nodes[ee][1]](j) *
                                 t_node_diff_ksi[edges_nodes[ee][0]](i) +
                             t_node_diff_ksi[edges_nodes[ee][0]](j) *
                                 t_node_diff_ksi[edges_nodes[ee][1]](i);

      (t_edge_n[ee])(i) = t_psi_e_0(i);
      ++(t_edge_n[ee]);
      (t_edge_n[ee])(i) = t_psi_e_1(i);
      ++(t_edge_n[ee]);

      (t_diff_edge_n[ee])(i, j) = t_diff_psi_e_0(i, j);
      ++(t_diff_edge_n[ee]);
      (t_diff_edge_n[ee])(i, j) = t_diff_psi_e_1(i, j);
      ++(t_diff_edge_n[ee]);

      if (p[ee] > 1) {

        const double ksi_0i = (N[node_shift + edges_nodes[ee][1]] -
                               N[node_shift + edges_nodes[ee][0]]) *
                              sense[ee];
        double psi_l[p[ee] + 1], diff_psi_l[3 * p[ee] + 3];
        CHKERR base_polynomials(p[ee], ksi_0i, &edge_diff_ksi[ee][0], psi_l,
                                diff_psi_l, 3);

        FTensor::Tensor1<double *, 3> t_diff_psi_l(
            &diff_psi_l[0], &diff_psi_l[p[ee] + 1], &diff_psi_l[2 * p[ee] + 2],
            1);

        for (int ll = 2; ll != P[ee]; ll++) {

          const double a = (double)(2 * ll + 1) / (double)(ll + 1);
          const double b = (double)(ll) / (double)(ll + 1);

          (t_edge_n[ee])(i) = a * psi_l[ll - 1] * t_psi_e_1(i) -
                              b * psi_l[ll - 2] * t_psi_e_0(i);
          ++(t_edge_n[ee]);

          (t_diff_edge_n[ee])(i, j) =
              -b * (t_diff_psi_l(j) * t_psi_e_0(i) +
                    psi_l[ll - 2] * t_diff_psi_e_0(i, j));
          ++t_diff_psi_l;
          (t_diff_edge_n[ee])(i, j) +=
              a * (t_diff_psi_l(j) * t_psi_e_1(i) +
                   psi_l[ll - 1] * t_diff_psi_e_1(i, j));
          ++(t_diff_edge_n[ee]);
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_EDGE(
    int sense, int p, double *N, double *diffN, double *edge_n,
    double *diff_edge_n, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  MoFEMFunctionBegin;

  if (NBEDGE_AINSWORTH_HCURL(p) == 0)
    MoFEMFunctionReturnHot(0);
  if (diff_edge_n != NULL)
    SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
            "Calculation of derivatives not implemented");

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_node_diff_ksi[2];
  t_node_diff_ksi[0](0) = diffN[0];
  t_node_diff_ksi[0](1) = 0;
  t_node_diff_ksi[0](2) = 0;
  t_node_diff_ksi[1](0) = diffN[1];
  t_node_diff_ksi[1](1) = 0;
  t_node_diff_ksi[1](2) = 0;

  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_edge_n(
      &edge_n[0], &edge_n[1], &edge_n[2]);
  FTensor::Tensor1<double, 3> t_psi_e_0, t_psi_e_1;

  for (int ii = 0; ii != nb_integration_pts; ii++) {

    const int node_shift = ii * 2;

    t_psi_e_0(i) = (N[node_shift + 1] * t_node_diff_ksi[0](i) -
                    N[node_shift + 0] * t_node_diff_ksi[1](i)) *
                   sense;
    t_psi_e_1(i) = N[node_shift + 1] * t_node_diff_ksi[0](i) +
                   N[node_shift + 0] * t_node_diff_ksi[1](i);

    t_edge_n(i) = t_psi_e_0(i);
    ++t_edge_n;

    t_edge_n(i) = t_psi_e_1(i);
    ++t_edge_n;

    if (p > 1) {

      const double ksi_0i = (N[node_shift + 1] - N[node_shift + 0]) * sense;
      double psi_l[p + 1];
      CHKERR base_polynomials(p, ksi_0i, NULL, psi_l, NULL, 3);

      for (int ll = 2; ll != NBEDGE_AINSWORTH_HCURL(p); ll++) {
        const double a = (double)(2 * ll + 1) / (double)(ll + 1);
        const double b = (double)(ll) / (double)(ll + 1);
        t_edge_n(i) =
            a * psi_l[ll - 1] * t_psi_e_1(i) - b * psi_l[ll - 2] * t_psi_e_0(i);
        ++t_edge_n;
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_EdgeBaseFunctions_MBTET_ON_FACE(
    int *sense, int *p, double *N, double *diffN, double *edge_n[],
    double *diff_edge_n[], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  // TODO This is not by atom tests properly

  const int edges_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  int P[3];
  for (int ee = 0; ee < 3; ee++)
    P[ee] = NBEDGE_AINSWORTH_HCURL(p[ee]);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  FTensor::Tensor1<double, 3> t_node_diff_ksi[3] = {
      FTensor::Tensor1<double, 3>(diffN[0], diffN[1], 0.),
      FTensor::Tensor1<double, 3>(diffN[2], diffN[3], 0.),
      FTensor::Tensor1<double, 3>(diffN[4], diffN[5], 0.),
  };
  FTensor::Tensor1<double, 2> t_2d_diff_ksi[3] = {
      FTensor::Tensor1<double, 2>(diffN[0], diffN[1]),
      FTensor::Tensor1<double, 2>(diffN[2], diffN[3]),
      FTensor::Tensor1<double, 2>(diffN[4], diffN[5])};

  FTensor::Tensor1<double *, 3> t_edge_n[3] = {
      FTensor::Tensor1<double *, 3>(&edge_n[0][0], &edge_n[0][1], &edge_n[0][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[1][0], &edge_n[1][1], &edge_n[1][2],
                                    3),
      FTensor::Tensor1<double *, 3>(&edge_n[2][0], &edge_n[2][1], &edge_n[2][2],
                                    3)};
  FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_diff_edge_n[3] = {
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
          &diff_edge_n[0][HVEC0_0], &diff_edge_n[0][HVEC0_1],
          &diff_edge_n[0][HVEC1_0], &diff_edge_n[0][HVEC1_1],
          &diff_edge_n[0][HVEC2_0], &diff_edge_n[0][HVEC2_1]),
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
          &diff_edge_n[1][HVEC0_0], &diff_edge_n[1][HVEC0_1],
          &diff_edge_n[1][HVEC1_0], &diff_edge_n[1][HVEC1_1],
          &diff_edge_n[1][HVEC2_0], &diff_edge_n[1][HVEC2_1]),
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
          &diff_edge_n[2][HVEC0_0], &diff_edge_n[2][HVEC0_1],
          &diff_edge_n[2][HVEC1_0], &diff_edge_n[2][HVEC1_1],
          &diff_edge_n[2][HVEC2_0], &diff_edge_n[2][HVEC2_1])};

  FTensor::Tensor1<double, 3> t_psi_e_0, t_psi_e_1;
  FTensor::Tensor2<double, 3, 2> t_diff_psi_e_0, t_diff_psi_e_1;

  for (int ee = 0; ee != 3; ee++) {

    if (P[ee] == 0)
      continue;
    const int node0 = edges_nodes[ee][0];
    const int node1 = edges_nodes[ee][1];
    const int sense_edge = sense[ee];

    t_diff_psi_e_0(i, j) =
        (t_node_diff_ksi[node0](i) * t_2d_diff_ksi[node1](j) -
         t_node_diff_ksi[node1](i) * t_2d_diff_ksi[node0](j)) *
        sense_edge;
    t_diff_psi_e_1(i, j) = t_node_diff_ksi[node0](i) * t_2d_diff_ksi[node1](j) +
                           t_node_diff_ksi[node1](i) * t_2d_diff_ksi[node0](j);

    for (int ii = 0; ii != nb_integration_pts; ii++) {

      const int node_shift = ii * 3;
      const double n0 = N[node_shift + node0];
      const double n1 = N[node_shift + node1];

      t_psi_e_0(i) =
          (n1 * t_node_diff_ksi[node0](i) - n0 * t_node_diff_ksi[node1](i)) *
          sense_edge;
      t_psi_e_1(i) =
          n1 * t_node_diff_ksi[node0](i) + n0 * t_node_diff_ksi[node1](i);

      (t_edge_n[ee])(i) = t_psi_e_0(i);
      (t_diff_edge_n[ee])(i, j) = t_diff_psi_e_0(i, j);
      ++(t_edge_n[ee]);
      ++(t_diff_edge_n[ee]);
      (t_edge_n[ee])(i) = t_psi_e_1(i);
      (t_diff_edge_n[ee])(i, j) = t_diff_psi_e_1(i, j);
      ++(t_edge_n[ee]);
      ++(t_diff_edge_n[ee]);

      if (p[ee] > 1) {
        const double ksi_0i = (n1 - n0) * sense_edge;
        double diff_ksi_0i[] = {
            ((t_2d_diff_ksi[node1])(0) - (t_2d_diff_ksi[node0])(0)) *
                sense_edge,
            ((t_2d_diff_ksi[node1])(1) - (t_2d_diff_ksi[node0])(1)) *
                sense_edge};

        double psi_l[p[ee] + 1], diff_psi_l[2 * p[ee] + 2];
        CHKERR
        base_polynomials(p[ee], ksi_0i, diff_ksi_0i, psi_l, diff_psi_l, 2);

        FTensor::Tensor1<double *, 2> t_diff_psi_ll_m1(
            &diff_psi_l[0 + 2 - 1], &diff_psi_l[p[ee] + 1 + 2 - 1], 1);
        FTensor::Tensor1<double *, 2> t_diff_psi_ll_m2(
            &diff_psi_l[0 + 2 - 2], &diff_psi_l[p[ee] + 1 + 2 - 2], 1);
        for (int ll = 2; ll != P[ee]; ll++) {
          const double a = (double)(2 * ll + 1) / (double)(ll + 1);
          const double b = (double)(ll) / (double)(ll + 1);
          (t_edge_n[ee])(i) = a * psi_l[ll - 1] * t_psi_e_1(i) -
                              b * psi_l[ll - 2] * t_psi_e_0(i);
          (t_diff_edge_n[ee])(i, j) = a * t_psi_e_1(i) * t_diff_psi_ll_m1(j) +
                                      a * psi_l[ll - 1] * t_diff_psi_e_1(i, j) -
                                      b * t_psi_e_0(i) * t_diff_psi_ll_m2(j) -
                                      b * psi_l[ll - 2] * t_diff_psi_e_0(i, j);
          ++(t_edge_n[ee]);
          ++(t_diff_edge_n[ee]);
          ++t_diff_psi_ll_m1;
          ++t_diff_psi_ll_m2;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f_e[4][3],
    double *diff_phi_f_e[4][3], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;
  const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};
  FTensor::Tensor1<double, 3> t_edge_diff_ksi;
  FTensor::Tensor1<double, 3> t_diff_beta_e;

  for (int ff = 0; ff != 4; ff++) {

    const int o_nodes[3] = {faces_nodes[3 * ff + 2], faces_nodes[3 * ff + 0],
                            faces_nodes[3 * ff + 1]};
    FTensor::Tensor1<double *, 3> t_opposite_node_diff[3] = {
        FTensor::Tensor1<double *, 3>(&diffN[3 * o_nodes[0] + 0],
                                      &diffN[3 * o_nodes[0] + 1],
                                      &diffN[3 * o_nodes[0] + 2]),
        FTensor::Tensor1<double *, 3>(&diffN[3 * o_nodes[1] + 0],
                                      &diffN[3 * o_nodes[1] + 1],
                                      &diffN[3 * o_nodes[1] + 2]),
        FTensor::Tensor1<double *, 3>(&diffN[3 * o_nodes[2] + 0],
                                      &diffN[3 * o_nodes[2] + 1],
                                      &diffN[3 * o_nodes[2] + 2])};
    double psi_l[p[ff] + 1], diff_psi_l[3 * p[ff] + 3];

    const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_EDGE_HCURL(p[ff]);
    // cerr << nb_base_fun_on_face << " " << p[ff] << endl;
    if (nb_base_fun_on_face == 0)
      continue;

    for (int ee = 0; ee != 3; ee++) {

      FTensor::Tensor1<double *, 3> t_face_edge_base(
          &phi_f_e[ff][ee][0], &phi_f_e[ff][ee][1], &phi_f_e[ff][ee][2], 3);
      FTensor::Tensor2<double *, 3, 3> t_diff_face_edge_base(
          &diff_phi_f_e[ff][ee][0], &diff_phi_f_e[ff][ee][3],
          &diff_phi_f_e[ff][ee][6], &diff_phi_f_e[ff][ee][1],
          &diff_phi_f_e[ff][ee][4], &diff_phi_f_e[ff][ee][7],
          &diff_phi_f_e[ff][ee][2], &diff_phi_f_e[ff][ee][5],
          &diff_phi_f_e[ff][ee][8], 9);
      const int en[] = {faces_nodes[3 * ff + edges[ee][0]],
                        faces_nodes[3 * ff + edges[ee][1]]};
      t_edge_diff_ksi(i) =
          t_node_diff_ksi[en[1]](i) - t_node_diff_ksi[en[0]](i);

      for (int ii = 0; ii != nb_integration_pts; ii++) {

        const int node_shift = ii * 4;
        const double n[] = {N[node_shift + faces_nodes[3 * ff + edges[ee][0]]],
                            N[node_shift + faces_nodes[3 * ff + edges[ee][1]]]};
        const double ksi_0i = n[1] - n[0];
        CHKERR base_polynomials(p[ff], ksi_0i, &t_edge_diff_ksi(0), psi_l,
                                diff_psi_l, 3);

        FTensor::Tensor1<double *, 3> t_diff_psi_l(
            &diff_psi_l[0], &diff_psi_l[p[ff] + 1], &diff_psi_l[2 * p[ff] + 2],
            1);

        const double beta_e = n[0] * n[1];
        t_diff_beta_e(j) =
            t_node_diff_ksi[en[0]](j) * n[1] + n[0] * t_node_diff_ksi[en[1]](j);

        for (int ll = 0; ll != nb_base_fun_on_face; ll++) {
          // if(ll == nb_base_fun_on_face-1) cerr << psi_l[ll] << endl;

          t_face_edge_base(i) =
              beta_e * psi_l[ll] * t_opposite_node_diff[ee](i);
          ++t_face_edge_base;

          t_diff_face_edge_base(i, j) =
              (beta_e * t_diff_psi_l(j) + t_diff_beta_e(j) * psi_l[ll]) *
              t_opposite_node_diff[ee](i);

          ++t_diff_face_edge_base;
          ++t_diff_psi_l;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f_e[3],
    double *diff_phi_f_e[3], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_EDGE_HCURL(p);
  if (nb_base_fun_on_face == 0)
    MoFEMFunctionReturnHot(0);

  const int edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  const int o_nodes[3] = {2, 0, 1};
  FTensor::Tensor2<double, 3, 3> t_opposite_node_diff(
      diffN[2 * o_nodes[0] + 0], diffN[2 * o_nodes[0] + 1], 0.,
      diffN[2 * o_nodes[1] + 0], diffN[2 * o_nodes[1] + 1], 0.,
      diffN[2 * o_nodes[2] + 0], diffN[2 * o_nodes[2] + 1], 0.);
  double psi_l[p + 1];
  double diff_psi_l[2 * p + 2];

  FTensor::Tensor1<double *, 3> t_face_edge_base[3] = {
      FTensor::Tensor1<double *, 3>(&phi_f_e[0][HVEC0], &phi_f_e[0][HVEC1],
                                    &phi_f_e[0][HVEC2], 3),
      FTensor::Tensor1<double *, 3>(&phi_f_e[1][HVEC0], &phi_f_e[1][HVEC1],
                                    &phi_f_e[1][HVEC2], 3),
      FTensor::Tensor1<double *, 3>(&phi_f_e[2][HVEC0], &phi_f_e[2][HVEC1],
                                    &phi_f_e[2][HVEC2], 3),
  };
  FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
      t_diff_face_edge_base[3] = {
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[0][HVEC0_0], &diff_phi_f_e[0][HVEC0_1],
              &diff_phi_f_e[0][HVEC1_0], &diff_phi_f_e[0][HVEC1_1],
              &diff_phi_f_e[0][HVEC2_0], &diff_phi_f_e[0][HVEC2_1]),
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[1][HVEC0_0], &diff_phi_f_e[1][HVEC0_1],
              &diff_phi_f_e[1][HVEC1_0], &diff_phi_f_e[1][HVEC1_1],
              &diff_phi_f_e[1][HVEC2_0], &diff_phi_f_e[1][HVEC2_1]),
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[2][HVEC0_0], &diff_phi_f_e[2][HVEC0_1],
              &diff_phi_f_e[2][HVEC1_0], &diff_phi_f_e[2][HVEC1_1],
              &diff_phi_f_e[2][HVEC2_0], &diff_phi_f_e[2][HVEC2_1])};

  for (int ee = 0; ee != 3; ee++) {

    const int node0 = faces_nodes[edges[ee][0]];
    const int node1 = faces_nodes[edges[ee][1]];
    double diff_ksi0i[] = {diffN[2 * node1 + 0] - diffN[2 * node0 + 0],
                           diffN[2 * node1 + 1] - diffN[2 * node0 + 1]};

    for (int ii = 0; ii != nb_integration_pts; ii++) {

      const int node_shift = ii * 3;
      const double n0 = N[node_shift + node0];
      const double n1 = N[node_shift + node1];
      const double ksi_0i = n1 - n0;
      CHKERR base_polynomials(p, ksi_0i, diff_ksi0i, psi_l, diff_psi_l, 2);

      const double beta_e = n0 * n1;
      FTensor::Tensor1<double, 2> t_diff_beta_e(
          diffN[2 * node0 + 0] * n1 + n0 * diffN[2 * node1 + 0],
          diffN[2 * node0 + 1] * n1 + n0 * diffN[2 * node1 + 1]);
      FTensor::Tensor1<double *, 2> t_diff_psi_l(&diff_psi_l[0],
                                                 &diff_psi_l[p + 1], 1);

      for (int ll = 0; ll != nb_base_fun_on_face; ll++) {
        t_face_edge_base[ee](i) =
            beta_e * psi_l[ll] * t_opposite_node_diff(ee, i);
        t_diff_face_edge_base[ee](i, j) =
            beta_e * t_opposite_node_diff(ee, i) * t_diff_psi_l(j) +
            psi_l[ll] * t_opposite_node_diff(ee, i) * t_diff_beta_e(j);
        ++t_face_edge_base[ee];
        ++t_diff_face_edge_base[ee];
        ++t_diff_psi_l;
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_BubbleFaceFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f[4],
    double *diff_phi_f[4], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  // double coords[] = { 0,0,0, 1,0,0, 0,1,0, 0,0,1 };
  // FTensor::Tensor1<double*,3> t_coords[4] = {
  //   FTensor::Tensor1<double*,3>(&coords[0],&coords[ 1],&coords[ 2]),
  //   FTensor::Tensor1<double*,3>(&coords[3],&coords[ 4],&coords[ 5]),
  //   FTensor::Tensor1<double*,3>(&coords[6],&coords[ 7],&coords[ 8]),
  //   FTensor::Tensor1<double*,3>(&coords[9],&coords[10],&coords[11])
  // };
  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};
  FTensor::Tensor1<double, 3> t_diff_ksi0i, t_diff_ksi0j;
  FTensor::Tensor1<double, 3> diff_beta_0ij;

  FTensor::Tensor1<double, 3> tou_0i;
  FTensor::Tensor1<double, 3> tou_0j;

  for (int ff = 0; ff != 4; ff++) {

    if (NBFACETRI_AINSWORTH_FACE_HCURL(p[ff]) == 0)
      continue;

    int n0 = faces_nodes[3 * ff + 0];
    int n1 = faces_nodes[3 * ff + 1];
    int n2 = faces_nodes[3 * ff + 2];

    // tou_0i(i) = t_coords[n1](i)-t_coords[n0](i);
    // tou_0j(i) = t_coords[n2](i)-t_coords[n0](i);
    tou_0i(i) = t_node_diff_ksi[n1](i) - t_node_diff_ksi[n0](i);
    tou_0j(i) = t_node_diff_ksi[n2](i) - t_node_diff_ksi[n0](i);

    t_diff_ksi0i(i) = t_node_diff_ksi[n1](i) - t_node_diff_ksi[n0](i);
    t_diff_ksi0j(i) = t_node_diff_ksi[n2](i) - t_node_diff_ksi[n0](i);

    double psi_l_0i[p[ff] + 1], diff_psi_l_0i[3 * p[ff] + 3];
    double psi_l_0j[p[ff] + 1], diff_psi_l_0j[3 * p[ff] + 3];

    FTensor::Tensor1<double *, 3> t_phi_f(&phi_f[ff][0], &phi_f[ff][1],
                                          &phi_f[ff][2], 3);
    FTensor::Tensor2<double *, 3, 3> t_diff_phi_f(
        &diff_phi_f[ff][0], &diff_phi_f[ff][3], &diff_phi_f[ff][6],
        &diff_phi_f[ff][1], &diff_phi_f[ff][4], &diff_phi_f[ff][7],
        &diff_phi_f[ff][2], &diff_phi_f[ff][5], &diff_phi_f[ff][8], 9);
    FTensor::Tensor1<double, 3> t_b;

    for (int ii = 0; ii != nb_integration_pts; ii++) {

      const int node_shift = ii * 4;
      const double beta_0ij =
          N[node_shift + n0] * N[node_shift + n1] * N[node_shift + n2];
      diff_beta_0ij(i) =
          t_node_diff_ksi[n0](i) * N[node_shift + n1] * N[node_shift + n2] +
          N[node_shift + n0] * t_node_diff_ksi[n1](i) * N[node_shift + n2] +
          N[node_shift + n0] * N[node_shift + n1] * t_node_diff_ksi[n2](i);

      const double ksi_0i = N[node_shift + n1] - N[node_shift + n0];
      CHKERR base_polynomials(p[ff], ksi_0i, &t_diff_ksi0i(0), psi_l_0i,
                              diff_psi_l_0i, 3);

      const double ksi_0j = N[node_shift + n2] - N[node_shift + n0];
      CHKERR base_polynomials(p[ff], ksi_0j, &t_diff_ksi0j(0), psi_l_0j,
                              diff_psi_l_0j, 3);

      int cc = 0;
      for (int oo = 0; oo <= (p[ff] - 3); oo++) {
        FTensor::Tensor1<double *, 3> t_diff_psi_l_0i(
            &diff_psi_l_0i[0], &diff_psi_l_0i[p[ff] + 1],
            &diff_psi_l_0i[2 * p[ff] + 2], 1);
        for (int pp0 = 0; pp0 <= oo; pp0++) {
          const int pp1 = oo - pp0;
          if (pp1 >= 0) {
            FTensor::Tensor1<double *, 3> t_diff_psi_l_0j(
                &diff_psi_l_0j[pp1], &diff_psi_l_0j[p[ff] + 1 + pp1],
                &diff_psi_l_0j[2 * p[ff] + 2 + pp1], 1);
            // base functions
            const double a = beta_0ij * psi_l_0i[pp0] * psi_l_0j[pp1];
            t_phi_f(i) = a * tou_0i(i);
            ++t_phi_f;
            ++cc;
            t_phi_f(i) = a * tou_0j(i);
            ++t_phi_f;
            ++cc;
            // diff base functions
            t_b(j) = diff_beta_0ij(j) * psi_l_0i[pp0] * psi_l_0j[pp1] +
                     beta_0ij * t_diff_psi_l_0i(j) * psi_l_0j[pp1] +
                     beta_0ij * psi_l_0i[pp0] * t_diff_psi_l_0j(j);
            t_diff_phi_f(i, j) = t_b(j) * tou_0i(i);
            ++t_diff_phi_f;
            t_diff_phi_f(i, j) = t_b(j) * tou_0j(i);
            ++t_diff_phi_f;
            ++t_diff_psi_l_0i;
          }
        }
      }
      const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_FACE_HCURL(p[ff]);
      if (cc != nb_base_fun_on_face) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong number of base functions %d != %d", cc,
                 nb_base_fun_on_face);
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_BubbleFaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  double zero = 0;
  FTensor::Tensor2<double *, 3, 3> t_node_diff_ksi(&diffN[0], &diffN[1], &zero,
                                                   &diffN[2], &diffN[3], &zero,
                                                   &diffN[4], &diffN[5], &zero);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  if (NBFACETRI_AINSWORTH_FACE_HCURL(p) == 0)
    MoFEMFunctionReturnHot(0);

  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  FTensor::Number<2> N2;

  const int node0 = faces_nodes[0];
  const int node1 = faces_nodes[1];
  const int node2 = faces_nodes[2];

  FTensor::Tensor1<double, 3> tou_0i;
  FTensor::Tensor1<double, 3> tou_0j;

  tou_0i(i) = t_node_diff_ksi(N1, i) - t_node_diff_ksi(N0, i);
  ;
  tou_0j(i) = t_node_diff_ksi(N2, i) - t_node_diff_ksi(N0, i);
  ;

  double psi_l_0i[p + 1], psi_l_0j[p + 1];
  double diff_psi_l_0i[2 * p + 2], diff_psi_l_0j[2 * p + 2];
  FTensor::Tensor1<double *, 3> t_phi_f(&phi_f[0], &phi_f[1], &phi_f[2], 3);
  FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_diff_phi_f(
      &diff_phi_f[HVEC0_0], &diff_phi_f[HVEC0_1], &diff_phi_f[HVEC1_0],
      &diff_phi_f[HVEC1_1], &diff_phi_f[HVEC2_0], &diff_phi_f[HVEC2_1]);

  double diff_ksi_0i[] = {t_node_diff_ksi(node1, 0) - t_node_diff_ksi(node0, 0),
                          t_node_diff_ksi(node1, 1) -
                              t_node_diff_ksi(node0, 1)};
  double diff_ksi_0j[] = {t_node_diff_ksi(node2, 0) - t_node_diff_ksi(node0, 0),
                          t_node_diff_ksi(node2, 1) -
                              t_node_diff_ksi(node0, 1)};

  for (int ii = 0; ii != nb_integration_pts; ii++) {

    const int node_shift = ii * 3;
    const double n0 = N[node_shift + node0];
    const double n1 = N[node_shift + node1];
    const double n2 = N[node_shift + node2];

    const double beta_0ij = n0 * n1 * n2;
    FTensor::Tensor1<double, 2> diff_beta_0ij(
        t_node_diff_ksi(node0, 0) * n1 * n2 +
            n0 * t_node_diff_ksi(node1, 0) * n2 +
            n0 * n1 * t_node_diff_ksi(node2, 0),
        t_node_diff_ksi(node0, 1) * n1 * n2 +
            n0 * t_node_diff_ksi(node1, 1) * n2 +
            n0 * n1 * t_node_diff_ksi(node2, 1));

    const double ksi_0i = N[node_shift + node1] - N[node_shift + node0];
    CHKERR base_polynomials(p, ksi_0i, diff_ksi_0i, psi_l_0i, diff_psi_l_0i, 2);

    const double ksi_0j = N[node_shift + node2] - N[node_shift + node0];
    CHKERR base_polynomials(p, ksi_0j, diff_ksi_0j, psi_l_0j, diff_psi_l_0j, 2);

    int cc = 0;
    FTensor::Tensor1<double, 2> t_diff_a;
    for (int oo = 0; oo <= (p - 3); oo++) {
      for (int pp0 = 0; pp0 <= oo; pp0++) {
        const int pp1 = oo - pp0;
        if (pp1 >= 0) {
          FTensor::Tensor1<double, 2> t_diff_psi_l_0i(
              diff_psi_l_0i[0 + pp0], diff_psi_l_0i[p + 1 + pp0]);
          FTensor::Tensor1<double, 2> t_diff_psi_l_0j(
              diff_psi_l_0j[0 + pp1], diff_psi_l_0j[p + 1 + pp1]);
          const double a = beta_0ij * psi_l_0i[pp0] * psi_l_0j[pp1];
          t_diff_a(j) = diff_beta_0ij(j) * psi_l_0i[pp0] * psi_l_0j[pp1] +
                        beta_0ij * psi_l_0i[pp0] * t_diff_psi_l_0j(j) +
                        beta_0ij * psi_l_0j[pp1] * t_diff_psi_l_0i(j);

          t_phi_f(i) = a * tou_0i(i);
          t_diff_phi_f(i, j) = tou_0i(i) * t_diff_a(j);
          ++t_phi_f;
          ++t_diff_phi_f;
          ++cc;
          t_phi_f(i) = a * tou_0j(i);
          t_diff_phi_f(i, j) = tou_0j(i) * t_diff_a(j);
          ++t_phi_f;
          ++t_diff_phi_f;
          ++cc;
        }
      }
    }

    const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_FACE_HCURL(p);
    if (cc != nb_base_fun_on_face) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of base functions %d != %d", cc,
               nb_base_fun_on_face);
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_FaceInteriorFunctions_MBTET(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_v,
    double *diff_phi_v, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  if (NBVOLUMETET_AINSWORTH_FACE_HCURL(p) == 0)
    MoFEMFunctionReturnHot(0);

  const int face_opposite_nodes[] = {2, 0, 1, 3};

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};
  FTensor::Tensor1<double, 3> t_diff_ksi0i, t_diff_ksi0j;

  MatrixDouble m_psi_l_0i(4, p + 1);
  MatrixDouble m_psi_l_0j(4, p + 1);
  MatrixDouble m_diff_psi_l_0i(4, 3 * p + 3);
  MatrixDouble m_diff_psi_l_0j(4, 3 * p + 3);

  double *psi_l_0i[] = {&m_psi_l_0i(0, 0), &m_psi_l_0i(1, 0), &m_psi_l_0i(2, 0),
                        &m_psi_l_0i(3, 0)};
  double *psi_l_0j[] = {&m_psi_l_0j(0, 0), &m_psi_l_0j(1, 0), &m_psi_l_0j(2, 0),
                        &m_psi_l_0j(3, 0)};
  double *diff_psi_l_0i[] = {&m_diff_psi_l_0i(0, 0), &m_diff_psi_l_0i(1, 0),
                             &m_diff_psi_l_0i(2, 0), &m_diff_psi_l_0i(3, 0)};
  double *diff_psi_l_0j[] = {&m_diff_psi_l_0j(0, 0), &m_diff_psi_l_0j(1, 0),
                             &m_diff_psi_l_0j(2, 0), &m_diff_psi_l_0j(3, 0)};
  double beta_f[4];

  FTensor::Tensor1<double, 3> t_diff_beta_f[4];

  FTensor::Tensor1<double *, 3> t_phi_v(&phi_v[0], &phi_v[1], &phi_v[2], 3);
  FTensor::Tensor2<double *, 3, 3> t_diff_phi_v(
      &diff_phi_v[0], &diff_phi_v[3], &diff_phi_v[6], &diff_phi_v[1],
      &diff_phi_v[4], &diff_phi_v[7], &diff_phi_v[2], &diff_phi_v[5],
      &diff_phi_v[8], 9);

  for (int ii = 0; ii != nb_integration_pts; ii++) {

    for (int ff = 0; ff != 4; ff++) {

      t_diff_ksi0i(i) = t_node_diff_ksi[faces_nodes[3 * ff + 1]](i) -
                        t_node_diff_ksi[faces_nodes[3 * ff + 0]](i);
      t_diff_ksi0j(i) = t_node_diff_ksi[faces_nodes[3 * ff + 2]](i) -
                        t_node_diff_ksi[faces_nodes[3 * ff + 0]](i);

      const int node_shift = ii * 4;

      beta_f[ff] = N[node_shift + faces_nodes[3 * ff + 0]] *
                   N[node_shift + faces_nodes[3 * ff + 1]] *
                   N[node_shift + faces_nodes[3 * ff + 2]];

      t_diff_beta_f[ff](j) = t_node_diff_ksi[faces_nodes[3 * ff + 0]](j) *
                                 N[node_shift + faces_nodes[3 * ff + 1]] *
                                 N[node_shift + faces_nodes[3 * ff + 2]] +
                             N[node_shift + faces_nodes[3 * ff + 0]] *
                                 t_node_diff_ksi[faces_nodes[3 * ff + 1]](j) *
                                 N[node_shift + faces_nodes[3 * ff + 2]] +
                             N[node_shift + faces_nodes[3 * ff + 0]] *
                                 N[node_shift + faces_nodes[3 * ff + 1]] *
                                 t_node_diff_ksi[faces_nodes[3 * ff + 2]](j);

      const double ksi_0i = N[node_shift + faces_nodes[3 * ff + 1]] -
                            N[node_shift + faces_nodes[3 * ff + 0]];
      CHKERR base_polynomials(p, ksi_0i, &t_diff_ksi0i(0), psi_l_0i[ff],
                              diff_psi_l_0i[ff], 3);

      const double ksi_0j = N[node_shift + faces_nodes[3 * ff + 2]] -
                            N[node_shift + faces_nodes[3 * ff + 0]];
      CHKERR base_polynomials(p, ksi_0j, &t_diff_ksi0j(0), psi_l_0j[ff],
                              diff_psi_l_0j[ff], 3);
    }

    int cc = 0;
    for (int oo = 0; oo <= (p - 3); oo++) {
      FTensor::Tensor1<double *, 3> t_diff_psi_l_0i[] = {
          FTensor::Tensor1<double *, 3>(&diff_psi_l_0i[0][0],
                                        &diff_psi_l_0i[0][p + 1],
                                        &diff_psi_l_0i[0][2 * p + 2], 1),
          FTensor::Tensor1<double *, 3>(&diff_psi_l_0i[1][0],
                                        &diff_psi_l_0i[1][p + 1],
                                        &diff_psi_l_0i[1][2 * p + 2], 1),
          FTensor::Tensor1<double *, 3>(&diff_psi_l_0i[2][0],
                                        &diff_psi_l_0i[2][p + 1],
                                        &diff_psi_l_0i[2][2 * p + 2], 1),
          FTensor::Tensor1<double *, 3>(&diff_psi_l_0i[3][0],
                                        &diff_psi_l_0i[3][p + 1],
                                        &diff_psi_l_0i[3][2 * p + 2], 1),
      };
      for (int pp0 = 0; pp0 <= oo; pp0++) {
        const int pp1 = oo - pp0;
        if (pp1 >= 0) {
          for (int ff = 0; ff != 4; ff++) {
            FTensor::Tensor1<double *, 3> t_diff_psi_l_0j(
                &m_diff_psi_l_0j(ff, pp1), &m_diff_psi_l_0j(ff, p + 1 + pp1),
                &m_diff_psi_l_0j(ff, 2 * p + 2 + pp1), 1);
            const double t = psi_l_0i[ff][pp0] * psi_l_0j[ff][pp1];
            const double a = beta_f[ff] * t;
            t_phi_v(i) = a * t_node_diff_ksi[face_opposite_nodes[ff]](i);
            ++t_phi_v;
            ++cc;
            t_diff_phi_v(i, j) =
                (t_diff_beta_f[ff](j) * t +
                 beta_f[ff] * t_diff_psi_l_0i[ff](j) * psi_l_0j[ff][pp1] +
                 beta_f[ff] * psi_l_0i[ff][pp0] * t_diff_psi_l_0j(j)) *
                t_node_diff_ksi[face_opposite_nodes[ff]](i);
            ++t_diff_phi_v;
            ++t_diff_psi_l_0i[ff];
          }
        }
      }
    }

    const int nb_base_fun_on_face = NBVOLUMETET_AINSWORTH_FACE_HCURL(p);
    if (cc != nb_base_fun_on_face) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of base functions %d != %d", cc,
               nb_base_fun_on_face);
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_VolumeInteriorFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  if (NBVOLUMETET_AINSWORTH_TET_HCURL(p) == 0)
    MoFEMFunctionReturnHot(0);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  FTensor::Number<2> N2;

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};

  double diff_ksi0i[3], diff_ksi0j[3], diff_ksi0k[3];
  FTensor::Tensor1<double *, 3> t_diff_ksi0i(diff_ksi0i, &diff_ksi0i[1],
                                             &diff_ksi0i[2]);
  FTensor::Tensor1<double *, 3> t_diff_ksi0j(diff_ksi0j, &diff_ksi0j[1],
                                             &diff_ksi0j[2]);
  FTensor::Tensor1<double *, 3> t_diff_ksi0k(diff_ksi0k, &diff_ksi0k[1],
                                             &diff_ksi0k[2]);
  t_diff_ksi0i(i) = t_node_diff_ksi[1](i) - t_node_diff_ksi[0](i);
  t_diff_ksi0j(i) = t_node_diff_ksi[2](i) - t_node_diff_ksi[0](i);
  t_diff_ksi0k(i) = t_node_diff_ksi[3](i) - t_node_diff_ksi[0](i);

  std::vector<double> v_psi_l_0i(p + 1), v_diff_psi_l_0i(3 * p + 3);
  std::vector<double> v_psi_l_0j(p + 1), v_diff_psi_l_0j(3 * p + 3);
  std::vector<double> v_psi_l_0k(p + 1), v_diff_psi_l_0k(3 * p + 3);
  double *psi_l_0i = &*v_psi_l_0i.begin();
  double *diff_psi_l_0i = &*v_diff_psi_l_0i.begin();
  double *psi_l_0j = &*v_psi_l_0j.begin();
  double *diff_psi_l_0j = &*v_diff_psi_l_0j.begin();
  double *psi_l_0k = &*v_psi_l_0k.begin();
  double *diff_psi_l_0k = &*v_diff_psi_l_0k.begin();

  FTensor::Tensor1<double *, 3> t_phi_v(&phi_v[0], &phi_v[1], &phi_v[2], 3);
  FTensor::Tensor2<double *, 3, 3> t_diff_phi_v(
      &diff_phi_v[0], &diff_phi_v[3], &diff_phi_v[6], &diff_phi_v[1],
      &diff_phi_v[4], &diff_phi_v[7], &diff_phi_v[2], &diff_phi_v[5],
      &diff_phi_v[8], 9);
  FTensor::Tensor1<double, 3> t_b;

  for (int ii = 0; ii != nb_integration_pts; ii++) {

    const int node_shift = ii * 4;
    const int n0 = node_shift + 0;
    const int n1 = node_shift + 1;
    const int n2 = node_shift + 2;
    const int n3 = node_shift + 3;

    const double beta_v = N[n0] * N[n1] * N[n2] * N[n3];

    const double ksi_0i = N[n1] - N[n0];
    CHKERR base_polynomials(p, ksi_0i, diff_ksi0i, psi_l_0i, diff_psi_l_0i, 3);

    const double ksi_0j = N[n2] - N[n0];
    CHKERR base_polynomials(p, ksi_0j, diff_ksi0j, psi_l_0j, diff_psi_l_0j, 3);

    const double ksi_0k = N[n3] - N[n0];
    CHKERR base_polynomials(p, ksi_0k, diff_ksi0k, psi_l_0k, diff_psi_l_0k, 3);

    FTensor::Tensor1<double, 3> t_diff_beta_v;
    t_diff_beta_v(j) = t_node_diff_ksi[0](j) * N[n1] * N[n2] * N[n3] +
                       N[n0] * t_node_diff_ksi[1](j) * N[n2] * N[n3] +
                       N[n0] * N[n1] * t_node_diff_ksi[2](j) * N[n3] +
                       N[n0] * N[n1] * N[n2] * t_node_diff_ksi[3](j);

    int cc = 0;
    for (int oo = 0; oo <= (p - 4); oo++) {
      FTensor::Tensor1<double *, 3> t_diff_psi_l_0i(
          &diff_psi_l_0i[0], &diff_psi_l_0i[p + 1], &diff_psi_l_0i[2 * p + 2],
          1);
      for (int pp0 = 0; pp0 <= oo; pp0++) {
        FTensor::Tensor1<double *, 3> t_diff_psi_l_0j(
            &diff_psi_l_0j[0], &diff_psi_l_0j[p + 1], &diff_psi_l_0j[2 * p + 2],
            1);
        for (int pp1 = 0; (pp0 + pp1) <= oo; pp1++) {
          const int pp2 = oo - pp0 - pp1;
          if (pp2 >= 0) {
            FTensor::Tensor1<double *, 3> t_diff_psi_l_0k(
                &diff_psi_l_0k[0 + pp2], &diff_psi_l_0k[p + 1 + pp2],
                &diff_psi_l_0k[2 * p + 2 + pp2], 1);
            const double t = psi_l_0i[pp0] * psi_l_0j[pp1] * psi_l_0k[pp2];
            const double a = beta_v * t;
            t_phi_v(0) = a;
            t_phi_v(1) = 0;
            t_phi_v(2) = 0;
            ++t_phi_v;
            ++cc;
            t_phi_v(0) = 0;
            t_phi_v(1) = a;
            t_phi_v(2) = 0;
            ++t_phi_v;
            ++cc;
            t_phi_v(0) = 0;
            t_phi_v(1) = 0;
            t_phi_v(2) = a;
            ++t_phi_v;
            ++cc;
            t_b(j) =
                t_diff_beta_v(j) * t +
                beta_v * (t_diff_psi_l_0i(j) * psi_l_0j[pp1] * psi_l_0k[pp2] +
                          psi_l_0i[pp0] * t_diff_psi_l_0j(j) * psi_l_0k[pp2] +
                          psi_l_0i[pp0] * psi_l_0j[pp1] * t_diff_psi_l_0k(j));
            t_diff_phi_v(N0, j) = t_b(j);
            t_diff_phi_v(N1, j) = 0;
            t_diff_phi_v(N2, j) = 0;
            ++t_diff_phi_v;
            t_diff_phi_v(N0, j) = 0;
            t_diff_phi_v(N1, j) = t_b(j);
            t_diff_phi_v(N2, j) = 0;
            ++t_diff_phi_v;
            t_diff_phi_v(N0, j) = 0;
            t_diff_phi_v(N1, j) = 0;
            t_diff_phi_v(N2, j) = t_b(j);
            ++t_diff_phi_v;
          }
          ++t_diff_psi_l_0j;
        }
        ++t_diff_psi_l_0i;
      }
    }

    const int nb_base_fun_on_face = NBVOLUMETET_AINSWORTH_TET_HCURL(p);
    if (cc != nb_base_fun_on_face) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong number of base functions %d != %d", cc,
               nb_base_fun_on_face);
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_FaceFunctions_MBTET(
    int *face_nodes, int *p, double *N, double *diffN, double *phi_f[4],
    double *diff_phi_f[4], int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  try {

    MatrixDouble base_face_edge_functions[4];
    MatrixDouble diff_base_face_edge_functions[4];
    double *phi_f_e[4][3];
    double *diff_phi_f_e[4][3];
    for (int ff = 0; ff != 4; ff++) {
      if (NBFACETRI_AINSWORTH_EDGE_HCURL(p[ff]) == 0) {
        for (int ee = 0; ee != 3; ee++) {
          phi_f_e[ff][ee] = NULL;
          diff_phi_f_e[ff][ee] = NULL;
        }
      } else {
        base_face_edge_functions[ff].resize(
            3, 3 * NBFACETRI_AINSWORTH_EDGE_HCURL(p[ff]) * nb_integration_pts);
        diff_base_face_edge_functions[ff].resize(
            3, 9 * NBFACETRI_AINSWORTH_EDGE_HCURL(p[ff]) * nb_integration_pts);
        // base_face_edge_functions[ff].clear();
        // diff_base_face_edge_functions[ff].clear();
        for (int ee = 0; ee != 3; ee++) {
          phi_f_e[ff][ee] = &base_face_edge_functions[ff](ee, 0);
          diff_phi_f_e[ff][ee] = &diff_base_face_edge_functions[ff](ee, 0);
        }
      }
    }
    CHKERR Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET(
        face_nodes, p, N, diffN, phi_f_e, diff_phi_f_e, nb_integration_pts,
        base_polynomials);

    VectorDouble base_face_bubble_functions[4];
    VectorDouble diff_base_face_bubble_functions[4];
    double *phi_f_f[4];
    double *diff_phi_f_f[4];
    for (int ff = 0; ff != 4; ff++) {
      int nb_dofs = NBFACETRI_AINSWORTH_FACE_HCURL(p[ff]);
      if (nb_dofs == 0) {
        phi_f_f[ff] = NULL;
        diff_phi_f_f[ff] = NULL;
      } else {
        base_face_bubble_functions[ff].resize(3 * nb_dofs * nb_integration_pts,
                                              false);
        diff_base_face_bubble_functions[ff].resize(
            9 * nb_dofs * nb_integration_pts, false);
        phi_f_f[ff] = &*base_face_bubble_functions[ff].data().begin();
        diff_phi_f_f[ff] = &*diff_base_face_bubble_functions[ff].data().begin();
      }
    }
    CHKERR Hcurl_Ainsworth_BubbleFaceFunctions_MBTET(
        face_nodes, p, N, diffN, phi_f_f, diff_phi_f_f, nb_integration_pts,
        base_polynomials);

    FTensor::Index<'i', 3> i;
    FTensor::Index<'j', 3> j;

    for (int ff = 0; ff != 4; ff++) {

      if (NBFACETRI_AINSWORTH_EDGE_HCURL(p[ff]) == 0)
        continue;

      FTensor::Tensor1<double *, 3> t_face_edge_base[] = {
          FTensor::Tensor1<double *, 3>(&phi_f_e[ff][0][0], &phi_f_e[ff][0][1],
                                        &phi_f_e[ff][0][2], 3),
          FTensor::Tensor1<double *, 3>(&phi_f_e[ff][1][0], &phi_f_e[ff][1][1],
                                        &phi_f_e[ff][1][2], 3),
          FTensor::Tensor1<double *, 3>(&phi_f_e[ff][2][0], &phi_f_e[ff][2][1],
                                        &phi_f_e[ff][2][2], 3)};
      FTensor::Tensor2<double *, 3, 3> t_diff_face_edge_base[] = {
          FTensor::Tensor2<double *, 3, 3>(
              &diff_phi_f_e[ff][0][0], &diff_phi_f_e[ff][0][3],
              &diff_phi_f_e[ff][0][6], &diff_phi_f_e[ff][0][1],
              &diff_phi_f_e[ff][0][4], &diff_phi_f_e[ff][0][7],
              &diff_phi_f_e[ff][0][2], &diff_phi_f_e[ff][0][5],
              &diff_phi_f_e[ff][0][8], 9),
          FTensor::Tensor2<double *, 3, 3>(
              &diff_phi_f_e[ff][1][0], &diff_phi_f_e[ff][1][3],
              &diff_phi_f_e[ff][1][6], &diff_phi_f_e[ff][1][1],
              &diff_phi_f_e[ff][1][4], &diff_phi_f_e[ff][1][7],
              &diff_phi_f_e[ff][1][2], &diff_phi_f_e[ff][1][5],
              &diff_phi_f_e[ff][1][8], 9),
          FTensor::Tensor2<double *, 3, 3>(
              &diff_phi_f_e[ff][2][0], &diff_phi_f_e[ff][2][3],
              &diff_phi_f_e[ff][2][6], &diff_phi_f_e[ff][2][1],
              &diff_phi_f_e[ff][2][4], &diff_phi_f_e[ff][2][7],
              &diff_phi_f_e[ff][2][2], &diff_phi_f_e[ff][2][5],
              &diff_phi_f_e[ff][2][8], 9)};

      FTensor::Tensor1<double *, 3> t_face_base(&phi_f[ff][0], &phi_f[ff][1],
                                                &phi_f[ff][2], 3);
      FTensor::Tensor2<double *, 3, 3> t_diff_face_base(
          &diff_phi_f[ff][0], &diff_phi_f[ff][3], &diff_phi_f[ff][6],
          &diff_phi_f[ff][1], &diff_phi_f[ff][4], &diff_phi_f[ff][7],
          &diff_phi_f[ff][2], &diff_phi_f[ff][5], &diff_phi_f[ff][8], 9);

      if (NBFACETRI_AINSWORTH_FACE_HCURL(p[ff]) > 0) {
        FTensor::Tensor1<double *, 3> t_face_face_base(
            &phi_f_f[ff][0], &phi_f_f[ff][1], &phi_f_f[ff][2], 3);
        FTensor::Tensor2<double *, 3, 3> t_diff_face_face_base(
            &diff_phi_f_f[ff][0], &diff_phi_f_f[ff][3], &diff_phi_f_f[ff][6],
            &diff_phi_f_f[ff][1], &diff_phi_f_f[ff][4], &diff_phi_f_f[ff][7],
            &diff_phi_f_f[ff][2], &diff_phi_f_f[ff][5], &diff_phi_f_f[ff][8],
            9);
        for (int ii = 0; ii != nb_integration_pts; ii++) {
          int cc = 0;
          for (int oo = 0; oo <= p[ff]; oo++) {
            // Face-edge base
            if (oo > 1) {
              for (int ee = 0; ee != 3; ee++) {
                for (int ll = NBFACETRI_AINSWORTH_EDGE_HCURL(oo - 1);
                     ll != NBFACETRI_AINSWORTH_EDGE_HCURL(oo); ll++) {
                  t_face_base(i) = t_face_edge_base[ee](i);
                  ++cc;
                  ++t_face_base;
                  ++t_face_edge_base[ee];
                  t_diff_face_base(i, j) = t_diff_face_edge_base[ee](i, j);
                  ++t_diff_face_base;
                  ++t_diff_face_edge_base[ee];
                  // cerr << oo << " " << ll << " " << cc << " " <<
                  // NBFACETRI_AINSWORTH_EDGE_HCURL(oo) << endl;
                }
              }
            }
            // Face-face base
            if (oo > 2) {
              for (int ll = NBFACETRI_AINSWORTH_FACE_HCURL(oo - 1);
                   ll != NBFACETRI_AINSWORTH_FACE_HCURL(oo); ll++) {
                t_face_base(i) = t_face_face_base(i);
                ++cc;
                ++t_face_base;
                ++t_face_face_base;
                t_diff_face_base(i, j) = t_diff_face_face_base(i, j);
                ++t_diff_face_base;
                ++t_diff_face_face_base;
              }
            }
          }
          // check consistency
          const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_HCURL(p[ff]);
          if (cc != nb_base_fun_on_face) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong number of base functions %d != %d", cc,
                     nb_base_fun_on_face);
          }
        }
      } else {
        for (int ii = 0; ii != nb_integration_pts; ii++) {
          int cc = 0;
          for (int oo = 0; oo <= p[ff]; oo++) {
            // Face-edge base
            if (oo > 1) {
              for (int ee = 0; ee != 3; ee++) {
                for (int ll = NBFACETRI_AINSWORTH_EDGE_HCURL(oo - 1);
                     ll != NBFACETRI_AINSWORTH_EDGE_HCURL(oo); ll++) {
                  t_face_base(i) = t_face_edge_base[ee](i);
                  ++cc;
                  ++t_face_base;
                  ++t_face_edge_base[ee];
                  t_diff_face_base(i, j) = t_diff_face_edge_base[ee](i, j);
                  ++t_diff_face_base;
                  ++t_diff_face_edge_base[ee];
                  // cerr << oo << " " << ll << " " << cc << " " <<
                  // NBFACETRI_AINSWORTH_EDGE_HCURL(oo) << endl;
                }
              }
            }
          }
          // check consistency
          const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_HCURL(p[ff]);
          if (cc != nb_base_fun_on_face) {
            SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                     "Wrong number of base functions %d != %d", cc,
                     nb_base_fun_on_face);
          }
        }
      }
    }

  } catch (MoFEMException const &e) {
    SETERRQ(PETSC_COMM_SELF, e.errorCode, e.errorMessage);
  } catch (std::exception &ex) {
    std::ostringstream ss;
    ss << "thorw in method: " << ex.what() << " at line " << __LINE__
       << " in file " << __FILE__;
    SETERRQ(PETSC_COMM_SELF, MOFEM_STD_EXCEPTION_THROW, ss.str().c_str());
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_FaceFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  if (NBFACETRI_AINSWORTH_EDGE_HCURL(p) == 0)
    MoFEMFunctionReturnHot(0);

  MatrixDouble base_face_edge_functions, diff_base_face_edge_functions;
  double *phi_f_e[3];
  double *diff_phi_f_e[3];
  base_face_edge_functions.resize(3, 3 * NBFACETRI_AINSWORTH_EDGE_HCURL(p) *
                                         nb_integration_pts);
  diff_base_face_edge_functions.resize(
      3, 2 * 3 * NBFACETRI_AINSWORTH_EDGE_HCURL(p) * nb_integration_pts);
  // base_face_edge_functions.clear();
  for (int ee = 0; ee != 3; ee++) {
    phi_f_e[ee] = &base_face_edge_functions(ee, 0);
    diff_phi_f_e[ee] = &diff_base_face_edge_functions(ee, 0);
  }
  CHKERR Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET_ON_FACE(
      faces_nodes, p, N, diffN, phi_f_e, diff_phi_f_e, nb_integration_pts,
      base_polynomials);

  VectorDouble base_face_bubble_functions;
  VectorDouble diff_base_face_bubble_functions;
  double *phi_f_f, *diff_phi_f_f;
  base_face_bubble_functions.resize(3 * NBFACETRI_AINSWORTH_FACE_HCURL(p) *
                                    nb_integration_pts);
  diff_base_face_bubble_functions.resize(
      2 * 3 * NBFACETRI_AINSWORTH_FACE_HCURL(p) * nb_integration_pts);
  phi_f_f = &*base_face_bubble_functions.data().begin();
  diff_phi_f_f = &*diff_base_face_bubble_functions.data().begin();
  CHKERR Hcurl_Ainsworth_BubbleFaceFunctions_MBTET_ON_FACE(
      faces_nodes, p, N, diffN, phi_f_f, diff_phi_f_f, nb_integration_pts,
      base_polynomials);

  // cerr << diff_base_face_bubble_functions << endl;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;

  FTensor::Tensor1<double *, 3> t_face_edge_base[] = {
      FTensor::Tensor1<double *, 3>(&phi_f_e[0][HVEC0], &phi_f_e[0][HVEC1],
                                    &phi_f_e[0][HVEC2], 3),
      FTensor::Tensor1<double *, 3>(&phi_f_e[1][HVEC0], &phi_f_e[1][HVEC1],
                                    &phi_f_e[1][HVEC2], 3),
      FTensor::Tensor1<double *, 3>(&phi_f_e[2][HVEC0], &phi_f_e[2][HVEC1],
                                    &phi_f_e[2][HVEC2], 3)};
  FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
      t_diff_face_edge_base[] = {
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[0][HVEC0_0], &diff_phi_f_e[0][HVEC0_1],
              &diff_phi_f_e[0][HVEC1_0], &diff_phi_f_e[0][HVEC1_1],
              &diff_phi_f_e[0][HVEC2_0], &diff_phi_f_e[0][HVEC2_1]),
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[1][HVEC0_0], &diff_phi_f_e[1][HVEC0_1],
              &diff_phi_f_e[1][HVEC1_0], &diff_phi_f_e[1][HVEC1_1],
              &diff_phi_f_e[1][HVEC2_0], &diff_phi_f_e[1][HVEC2_1]),
          FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
              &diff_phi_f_e[2][HVEC0_0], &diff_phi_f_e[2][HVEC0_1],
              &diff_phi_f_e[2][HVEC1_0], &diff_phi_f_e[2][HVEC1_1],
              &diff_phi_f_e[2][HVEC2_0], &diff_phi_f_e[2][HVEC2_1])};

  FTensor::Tensor1<double *, 3> t_face_base(&phi_f[0], &phi_f[1], &phi_f[2], 3);
  FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_diff_face_base(
      &diff_phi_f[HVEC0_0], &diff_phi_f[HVEC0_1], &diff_phi_f[HVEC1_0],
      &diff_phi_f[HVEC1_1], &diff_phi_f[HVEC2_0], &diff_phi_f[HVEC2_1]);

  if (NBFACETRI_AINSWORTH_FACE_HCURL(p) > 0) {
    FTensor::Tensor1<double *, 3> t_face_face_base(
        &phi_f_f[HVEC0], &phi_f_f[HVEC1], &phi_f_f[HVEC2], 3);
    FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_diff_face_face_base(
        &diff_phi_f_f[HVEC0_0], &diff_phi_f_f[HVEC0_1], &diff_phi_f_f[HVEC1_0],
        &diff_phi_f_f[HVEC1_1], &diff_phi_f_f[HVEC2_0], &diff_phi_f_f[HVEC2_1]);
    for (int ii = 0; ii != nb_integration_pts; ii++) {
      int cc = 0;
      for (int oo = 0; oo <= p; oo++) {
        // Face-edge base
        if (oo > 1) {
          for (int ee = 0; ee != 3; ee++) {
            for (int ll = NBFACETRI_AINSWORTH_EDGE_HCURL(oo - 1);
                 ll != NBFACETRI_AINSWORTH_EDGE_HCURL(oo); ll++) {
              t_face_base(i) = t_face_edge_base[ee](i);
              t_diff_face_base(i, j) = t_diff_face_edge_base[ee](i, j);
              ++cc;
              ++t_face_base;
              ++t_face_edge_base[ee];
              ++t_diff_face_base;
              ++t_diff_face_edge_base[ee];
            }
          }
        }
        // Face-face base
        if (oo > 2) {
          for (int ll = NBFACETRI_AINSWORTH_FACE_HCURL(oo - 1);
               ll != NBFACETRI_AINSWORTH_FACE_HCURL(oo); ll++) {
            t_face_base(i) = t_face_face_base(i);
            t_diff_face_base(i, j) = t_diff_face_face_base(i, j);
            ++cc;
            ++t_face_base;
            ++t_face_face_base;
            ++t_diff_face_base;
            ++t_diff_face_face_base;
          }
        }
      }
      // check consistency
      const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_HCURL(p);
      if (cc != nb_base_fun_on_face) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong number of base functions %d != %d", cc,
                 nb_base_fun_on_face);
      }
    }
  } else {
    for (int ii = 0; ii != nb_integration_pts; ii++) {
      int cc = 0;
      for (int oo = 0; oo <= p; oo++) {
        // Face-edge base
        if (oo > 1) {
          for (int ee = 0; ee != 3; ee++) {
            for (int ll = NBFACETRI_AINSWORTH_EDGE_HCURL(oo - 1);
                 ll != NBFACETRI_AINSWORTH_EDGE_HCURL(oo); ll++) {
              t_face_base(i) = t_face_edge_base[ee](i);
              t_diff_face_base(i, j) = t_diff_face_edge_base[ee](i, j);
              ++cc;
              ++t_face_base;
              ++t_face_edge_base[ee];
              ++t_diff_face_base;
              ++t_diff_face_edge_base[ee];
              // cerr << oo << " " << ll << " " << cc << " " <<
              // NBFACETRI_AINSWORTH_EDGE_HCURL(oo) << endl;
            }
          }
        }
      }
      // check consistency
      const int nb_base_fun_on_face = NBFACETRI_AINSWORTH_HCURL(p);
      if (cc != nb_base_fun_on_face) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong number of base functions %d != %d", cc,
                 nb_base_fun_on_face);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hcurl_Ainsworth_VolumeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int nb_integration_pts,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

  VectorDouble base_face_inetrior_functions(
      3 * NBVOLUMETET_AINSWORTH_FACE_HCURL(p) * nb_integration_pts);
  VectorDouble diff_base_face_inetrior_functions(
      9 * NBVOLUMETET_AINSWORTH_FACE_HCURL(p) * nb_integration_pts);
  // base_face_inetrior_functions.clear();
  // diff_base_face_inetrior_functions.clear();
  double *phi_v_f = &*base_face_inetrior_functions.data().begin();
  double *diff_phi_v_f = &*diff_base_face_inetrior_functions.data().begin();
  int faces_nodes[] = {0, 1, 3, 1, 2, 3, 0, 2, 3, 0, 1, 2};
  CHKERR Hcurl_Ainsworth_FaceInteriorFunctions_MBTET(
      faces_nodes, p, N, diffN, phi_v_f, diff_phi_v_f, nb_integration_pts,
      base_polynomials);

  VectorDouble base_interior_functions(3 * NBVOLUMETET_AINSWORTH_TET_HCURL(p) *
                                       nb_integration_pts);
  VectorDouble diff_base_interior_functions(
      9 * NBVOLUMETET_AINSWORTH_TET_HCURL(p) * nb_integration_pts);
  // base_interior_functions.clear();
  // diff_base_interior_functions.clear();
  double *phi_v_v = &*base_interior_functions.data().begin();
  double *diff_phi_v_v = &*diff_base_interior_functions.data().begin();
  CHKERR Hcurl_Ainsworth_VolumeInteriorFunctions_MBTET(
      p, N, diffN, phi_v_v, diff_phi_v_v, nb_integration_pts, base_polynomials);

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double *, 3> t_face_interior(&phi_v_f[0], &phi_v_f[1],
                                                &phi_v_f[2], 3);
  FTensor::Tensor2<double *, 3, 3> t_diff_face_interior(
      &diff_phi_v_f[0], &diff_phi_v_f[3], &diff_phi_v_f[6], &diff_phi_v_f[1],
      &diff_phi_v_f[4], &diff_phi_v_f[7], &diff_phi_v_f[2], &diff_phi_v_f[5],
      &diff_phi_v_f[8], 9);

  FTensor::Tensor1<double *, 3> t_phi_v(&phi_v[0], &phi_v[1], &phi_v[2], 3);
  FTensor::Tensor2<double *, 3, 3> t_diff_phi_v(
      &diff_phi_v[0], &diff_phi_v[3], &diff_phi_v[6], &diff_phi_v[1],
      &diff_phi_v[4], &diff_phi_v[7], &diff_phi_v[2], &diff_phi_v[5],
      &diff_phi_v[8], 9);

  if (NBVOLUMETET_AINSWORTH_TET_HCURL(p) > 0) {
    FTensor::Tensor1<double *, 3> t_volume_interior(&phi_v_v[0], &phi_v_v[1],
                                                    &phi_v_v[2], 3);
    FTensor::Tensor2<double *, 3, 3> t_diff_volume_interior(
        &diff_phi_v_v[0], &diff_phi_v_v[3], &diff_phi_v_v[6], &diff_phi_v_v[1],
        &diff_phi_v_v[4], &diff_phi_v_v[7], &diff_phi_v_v[2], &diff_phi_v_v[5],
        &diff_phi_v_v[8], 9);
    for (int ii = 0; ii != nb_integration_pts; ii++) {
      int cc = 0;
      for (int oo = 0; oo <= p; oo++) {
        for (int ll = NBVOLUMETET_AINSWORTH_FACE_HCURL(oo - 1);
             ll != NBVOLUMETET_AINSWORTH_FACE_HCURL(oo); ll++) {
          t_phi_v(i) = t_face_interior(i);
          ++t_phi_v;
          ++t_face_interior;
          ++cc;
          t_diff_phi_v(i, j) = t_diff_face_interior(i, j);
          ++t_diff_phi_v;
          ++t_diff_face_interior;
        }
        for (int ll = NBVOLUMETET_AINSWORTH_TET_HCURL(oo - 1);
             ll != NBVOLUMETET_AINSWORTH_TET_HCURL(oo); ll++) {
          t_phi_v(i) = t_volume_interior(i);
          ++t_phi_v;
          ++t_volume_interior;
          ++cc;
          t_diff_phi_v(i, j) = t_diff_volume_interior(i, j);
          ++t_diff_phi_v;
          ++t_diff_volume_interior;
        }
      }
      // check consistency
      const int nb_base_fun_on_face = NBVOLUMETET_AINSWORTH_HCURL(p);
      if (cc != nb_base_fun_on_face) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong number of base functions %d != %d", cc,
                 nb_base_fun_on_face);
      }
    }
  } else {
    for (int ii = 0; ii != nb_integration_pts; ii++) {
      int cc = 0;
      for (int oo = 0; oo <= p; oo++) {
        for (int ll = NBVOLUMETET_AINSWORTH_FACE_HCURL(oo - 1);
             ll != NBVOLUMETET_AINSWORTH_FACE_HCURL(oo); ll++) {
          t_phi_v(i) = t_face_interior(i);
          ++t_phi_v;
          ++t_face_interior;
          ++cc;
          t_diff_phi_v(i, j) = t_diff_face_interior(i, j);
          ++t_diff_phi_v;
          ++t_diff_face_interior;
        }
      }
      // check consistency
      const int nb_base_fun_on_face = NBVOLUMETET_AINSWORTH_HCURL(p);
      if (cc != nb_base_fun_on_face) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "Wrong number of base functions %d != %d", cc,
                 nb_base_fun_on_face);
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

#endif // Not GENERATE_VTK_WITH_CURL_BASE

#ifdef GENERATE_VTK_WITH_CURL_BASE

#include <MoFEM.hpp>
#include <Hcurl.hpp>
using namespace MoFEM;
using namespace boost::numeric;

MoFEMErrorCode VTK_Ainsworth_Hcurl_MBTET(const string file_name) {
  MoFEMFunctionBeginHot;

  double base_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  EntityHandle nodes[4];
  for (int nn = 0; nn < 4; nn++) {
    CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
  }
  EntityHandle tet;
  CHKERR moab_ref.create_element(MBTET, nodes, 4, tet);

  MoFEM::Core m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  CHKERR m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
      0, 3, BitRefLevel().set(0));

  const int max_level = 4;
  for (int ll = 0; ll != max_level; ll++) {
    Range edges;
    CHKERR m_field_ref.getInterface<BitRefManager>()
        ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                       BitRefLevel().set(), MBEDGE, edges);
    Range tets;
    CHKERR m_field_ref.getInterface<BitRefManager>()
        ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                       BitRefLevel(ll).set(), MBTET, tets);
    // refine mesh
    MeshRefinement *m_ref;
    CHKERR m_field_ref.getInterface(m_ref);
    CHKERR m_ref->add_vertices_in_the_middle_of_edges(
        edges, BitRefLevel().set(ll + 1));
    CHKERR m_ref->refine_TET(tets, BitRefLevel().set(ll + 1));
  }

  Range tets;
  CHKERR m_field_ref.getInterface<BitRefManager>()
      ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(max_level),
                                     BitRefLevel().set(max_level), MBTET, tets);

  // Use 10 node tets to print base
  if (1) {

    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
    CHKERR moab_ref.add_entities(meshset, tets);
    CHKERR moab_ref.convert_entities(meshset, true, false, false);
    CHKERR moab_ref.delete_entities(&meshset, 1);
  }

  Range elem_nodes;
  CHKERR moab_ref.get_connectivity(tets, elem_nodes, false);

  const int nb_gauss_pts = elem_nodes.size();
  MatrixDouble gauss_pts(nb_gauss_pts, 4);
  gauss_pts.clear();
  Range::iterator nit = elem_nodes.begin();
  for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
    CHKERR moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
  }
  gauss_pts = trans(gauss_pts);

  MatrixDouble shape_fun;
  shape_fun.resize(nb_gauss_pts, 4);
  CHKERR ShapeMBTET(&*shape_fun.data().begin(), &gauss_pts(0, 0),
                    &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);

  double diff_shape_fun[12];
  CHKERR ShapeDiffMBTET(diff_shape_fun);

  // int edge_sense[6] = { 1,1,1, 1,1,1 };
  const int order = 5;
  // int edge_order[6] = { order,order,order, order,order,order };
  double def_val[] = {0, 0, 0, 0, 0, 0};
  int faces_order[] = {order, order, order, order};
  int faces_nodes[] = {0, 1, 3, 1, 2, 3, 0, 2, 3, 0, 1, 2};

  // cout << "NBEDGE_AINSWORTH_HCURL " <<  NBEDGE_AINSWORTH_HCURL(order) <<
  // endl; MatrixDouble
  // base_edge_functions(6,3*nb_gauss_pts*NBEDGE_AINSWORTH_HCURL(order));
  // double* edge_n[] = {
  //   &base_edge_functions(0,0),
  //   &base_edge_functions(1,0),
  //   &base_edge_functions(2,0),
  //   &base_edge_functions(3,0),
  //   &base_edge_functions(4,0),
  //   &base_edge_functions(5,0)
  // };
  //
  // MatrixDouble
  // diff_base_edge_functions(6,9*nb_gauss_pts*NBEDGE_AINSWORTH_HCURL(order));
  // double* diff_edge_n[] = {
  //   &diff_base_edge_functions(0,0),
  //   &diff_base_edge_functions(1,0),
  //   &diff_base_edge_functions(2,0),
  //   &diff_base_edge_functions(3,0),
  //   &diff_base_edge_functions(4,0),
  //   &diff_base_edge_functions(5,0)
  // };
  //
  // CHKERR Hcurl_Ainsworth_EdgeBaseFunctions_MBTET(
  //   edge_sense,
  //   edge_order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   edge_n,
  //   diff_edge_n,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  //
  //
  // for(int  ee = 0;ee!=6;ee++) {
  //   for(int ll = 0;ll!=NBEDGE_AINSWORTH_HCURL(order);ll++) {
  //     std::ostringstream ss;
  //     ss << "curl_edge_" << ee << "_" << ll;
  //     Tag th;
  //     CHKERR moab_ref.tag_get_handle(
  //       ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //     );
  //     std::ostringstream ss_grad;
  //     ss_grad << "grad_curl_edge_" << ee << "_" << ll;
  //     Tag th_grad;
  //     CHKERR moab_ref.tag_get_handle(
  //       ss_grad.str().c_str(),9,MB_TYPE_DOUBLE,th_grad,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //     );
  //
  //     int gg = 0;
  //     for(Range::iterator nit =
  //     elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //       CHKERR moab_ref.tag_set_data(
  //         th,&*nit,1,&(edge_n[ee][gg*3*NBEDGE_AINSWORTH_HCURL(order)+ll*3])
  //       );
  //       int sh = gg*9*NBEDGE_AINSWORTH_HCURL(order)+ll*9;
  //       double grad[9] = {
  //         diff_edge_n[ee][sh+0],diff_edge_n[ee][sh+3],diff_edge_n[ee][sh+6],
  //         diff_edge_n[ee][sh+1],diff_edge_n[ee][sh+4],diff_edge_n[ee][sh+7],
  //         diff_edge_n[ee][sh+2],diff_edge_n[ee][sh+5],diff_edge_n[ee][sh+8]
  //       };
  //       CHKERR moab_ref.tag_set_data(th_grad,&*nit,1,grad);
  //
  //     }
  //   }
  // }

  // cout << "NBFACETRI_AINSWORTH_EDGE_HCURL(order) " <<
  // NBFACETRI_AINSWORTH_EDGE_HCURL(order) << endl; MatrixDouble
  // base_face_edge_functions(
  //   4,3*3*NBFACETRI_AINSWORTH_EDGE_HCURL(order)*nb_gauss_pts
  // );
  // MatrixDouble diff_base_face_edge_functions(
  //   4,3*9*NBFACETRI_AINSWORTH_EDGE_HCURL(order)*nb_gauss_pts
  // );
  // double *phi_f_e[4][3];
  // double *diff_phi_f_e[4][3];
  // for(int ff = 0;ff!=4;ff++) {
  //   for(int ee = 0;ee!=3;ee++) {
  //     phi_f_e[ff][ee] =
  //     &base_face_edge_functions(ff,ee*3*NBFACETRI_AINSWORTH_EDGE_HCURL(order)*nb_gauss_pts);
  //     diff_phi_f_e[ff][ee] =
  //     &diff_base_face_edge_functions(ff,ee*9*NBFACETRI_AINSWORTH_EDGE_HCURL(order)*nb_gauss_pts);
  //   }
  // }
  //
  // CHKERR Hcurl_Ainsworth_EdgeBasedFaceFunctions_MBTET(
  //   faces_nodes,
  //   faces_order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   phi_f_e,
  //   diff_phi_f_e,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  //
  // for(int ff = 0;ff!=4;ff++) {
  //   for(int  ee = 0;ee!=3;ee++) {
  //     for(int ll = 0;ll!=NBFACETRI_AINSWORTH_EDGE_HCURL(order);ll++) {
  //       std::ostringstream ss;
  //       ss << "curl_face_edge_" << ff << "_" << ee << "_" << ll;
  //       Tag th;
  //       CHKERR moab_ref.tag_get_handle(
  //         ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //       );
  //
  //       std::ostringstream ss_grad;
  //       ss_grad << "grad_curl_face_edge_" << ff << "_" << ee << "_" << ll;
  //       Tag th_grad;
  //       CHKERR moab_ref.tag_get_handle(
  //         ss_grad.str().c_str(),9,MB_TYPE_DOUBLE,th_grad,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //       );
  //
  //       int gg = 0;
  //       for(Range::iterator nit =
  //       elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //
  //         int idx =
  //         3*NBFACETRI_AINSWORTH_EDGE_HCURL(order)*gg+ll*3;
  //         if(idx >= base_face_edge_functions.size2()) {
  //           cerr << ff << " " << ee << " " << ll << " " << gg << endl;
  //         }
  //
  //         CHKERR moab_ref.tag_set_data(th,&*nit,1,&(phi_f_e[ff][ee][idx]));
  //
  //
  //         int sh = gg*9*NBFACETRI_AINSWORTH_EDGE_HCURL(order)+ll*9;
  //         double grad[9] = {
  //           diff_phi_f_e[ff][ee][sh+0],diff_phi_f_e[ff][ee][sh+3],diff_phi_f_e[ff][ee][sh+6],
  //           diff_phi_f_e[ff][ee][sh+1],diff_phi_f_e[ff][ee][sh+4],diff_phi_f_e[ff][ee][sh+7],
  //           diff_phi_f_e[ff][ee][sh+2],diff_phi_f_e[ff][ee][sh+5],diff_phi_f_e[ff][ee][sh+8]
  //         };
  //         CHKERR moab_ref.tag_set_data(th_grad,&*nit,1,grad);
  //
  //
  //       }
  //     }
  //   }
  // }
  //
  cout << "NBFACETRI_AINSWORTH_FACE_HCURL "
       << NBFACETRI_AINSWORTH_FACE_HCURL(order) << endl;
  MatrixDouble base_face_bubble_functions(
      4, 3 * NBFACETRI_AINSWORTH_FACE_HCURL(order) * nb_gauss_pts);
  MatrixDouble diff_base_face_bubble_functions(
      4, 9 * NBFACETRI_AINSWORTH_FACE_HCURL(order) * nb_gauss_pts);
  double *phi_f[4];
  double *diff_phi_f[4];
  for (int ff = 0; ff != 4; ff++) {
    phi_f[ff] = &base_face_bubble_functions(ff, 0);
    diff_phi_f[ff] = &diff_base_face_bubble_functions(ff, 0);
  }

  CHKERR Hcurl_Ainsworth_BubbleFaceFunctions_MBTET(
      faces_nodes, faces_order, &*shape_fun.data().begin(), diff_shape_fun,
      phi_f, diff_phi_f, nb_gauss_pts, Legendre_polynomials);

  for (int ff = 0; ff != 4; ff++) {
    for (int ll = 0; ll != NBFACETRI_AINSWORTH_FACE_HCURL(order); ll++) {
      std::ostringstream ss;
      ss << "curl_face_bubble_" << ff << "_" << ll;
      Tag th;
      CHKERR moab_ref.tag_get_handle(ss.str().c_str(), 3, MB_TYPE_DOUBLE, th,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, def_val);
      std::ostringstream grad_ss;
      grad_ss << "grad_curl_face_bubble_" << ff << "_" << ll;
      Tag th_grad;
      CHKERR moab_ref.tag_get_handle(grad_ss.str().c_str(), 9, MB_TYPE_DOUBLE,
                                     th_grad, MB_TAG_CREAT | MB_TAG_SPARSE,
                                     def_val);

      int gg = 0;
      for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
           nit++, gg++) {
        int idx = 3 * NBFACETRI_AINSWORTH_FACE_HCURL(order) * gg + ll * 3;
        CHKERR moab_ref.tag_set_data(th, &*nit, 1, &(phi_f[ff][idx]));
        int sh = gg * 9 * NBFACETRI_AINSWORTH_FACE_HCURL(order) + ll * 9;
        double grad[9] = {diff_phi_f[ff][sh + 0], diff_phi_f[ff][sh + 3],
                          diff_phi_f[ff][sh + 6], diff_phi_f[ff][sh + 1],
                          diff_phi_f[ff][sh + 4], diff_phi_f[ff][sh + 7],
                          diff_phi_f[ff][sh + 2], diff_phi_f[ff][sh + 5],
                          diff_phi_f[ff][sh + 8]};
        CHKERR moab_ref.tag_set_data(th_grad, &*nit, 1, grad);
      }
    }
  }

  // cout << "NBVOLUMETET_AINSWORTH_FACE_HCURL " <<
  // NBVOLUMETET_AINSWORTH_FACE_HCURL(order) << endl; VectorDouble
  // base_face_inetrior_functions(3*NBVOLUMETET_AINSWORTH_FACE_HCURL(order)*nb_gauss_pts);
  // VectorDouble
  // diff_base_face_inetrior_functions(9*NBVOLUMETET_AINSWORTH_FACE_HCURL(order)*nb_gauss_pts);
  // double *phi_v_f = &base_face_inetrior_functions[0];
  // double *diff_phi_v_f = &diff_base_face_inetrior_functions[0];
  // CHKERR Hcurl_Ainsworth_FaceInteriorFunctions_MBTET(
  //   faces_nodes,
  //   order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   phi_v_f,
  //   diff_phi_v_f,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  // for(int ll = 0;ll!=NBVOLUMETET_AINSWORTH_FACE_HCURL(order);ll++) {
  //
  //   std::ostringstream ss;
  //   ss << "curl_face_interior_" << ll;
  //   Tag th;
  //   CHKERR moab_ref.tag_get_handle(
  //     ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //   );
  //
  //   std::ostringstream ss_grad;
  //   ss_grad << "grad_curl_face_interior_" << ll;
  //   Tag th_grad;
  //   CHKERR moab_ref.tag_get_handle(
  //     ss_grad.str().c_str(),9,MB_TYPE_DOUBLE,th_grad,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //   );
  //
  //   int gg = 0;
  //   for(Range::iterator nit =
  //   elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //     int idx = 3*NBVOLUMETET_AINSWORTH_FACE_HCURL(order)*gg+ll*3;
  //     CHKERR moab_ref.tag_set_data(th,&*nit,1,&(phi_v_f[idx]));
  //      int sh =
  //     gg*9*NBVOLUMETET_AINSWORTH_FACE_HCURL(order)+ll*9; double grad[9] = {
  //       diff_phi_v_f[sh+0],diff_phi_v_f[sh+3],diff_phi_v_f[sh+6],
  //       diff_phi_v_f[sh+1],diff_phi_v_f[sh+4],diff_phi_v_f[sh+7],
  //       diff_phi_v_f[sh+2],diff_phi_v_f[sh+5],diff_phi_v_f[sh+8]
  //     };
  //     CHKERR moab_ref.tag_set_data(th_grad,&*nit,1,grad);
  //   }
  // }

  // cout << "NBVOLUMETET_AINSWORTH_TET_HCURL " <<
  // NBVOLUMETET_AINSWORTH_TET_HCURL(order) << endl; VectorDouble
  // base_interior_functions(3*NBVOLUMETET_AINSWORTH_TET_HCURL(order)*nb_gauss_pts);
  // VectorDouble
  // diff_base_interior_functions(9*NBVOLUMETET_AINSWORTH_TET_HCURL(order)*nb_gauss_pts);
  // double *phi_v = &base_interior_functions[0];
  // double *diff_phi_v = &diff_base_interior_functions[0];
  // CHKERR Hcurl_Ainsworth_VolumeInteriorFunctions_MBTET(
  //   order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   phi_v,
  //   diff_phi_v,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  // for(int ll = 0;ll!=NBVOLUMETET_AINSWORTH_TET_HCURL(order);ll++) {
  //
  //   std::ostringstream ss;
  //   ss << "curl_interior_" << ll;
  //   Tag th;
  //   CHKERR moab_ref.tag_get_handle(
  //     ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //   );
  //
  //   std::ostringstream ss_gard;
  //   ss_gard << "grad_curl_interior_" << ll;
  //   Tag th_grad;
  //   CHKERR moab_ref.tag_get_handle(
  //     ss_gard.str().c_str(),9,MB_TYPE_DOUBLE,th_grad,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //   );
  //
  //   int gg = 0;
  //   for(Range::iterator nit =
  //   elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //     int idx = 3*NBVOLUMETET_AINSWORTH_TET_HCURL(order)*gg+ll*3;
  //     CHKERR moab_ref.tag_set_data(th,&*nit,1,&(phi_v[idx]));
  //      int sh =
  //     gg*9*NBVOLUMETET_AINSWORTH_TET_HCURL(order)+ll*9; double grad[9] = {
  //       diff_phi_v[sh+0],diff_phi_v[sh+3],diff_phi_v[sh+6],
  //       diff_phi_v[sh+1],diff_phi_v[sh+4],diff_phi_v[sh+7],
  //       diff_phi_v[sh+2],diff_phi_v[sh+5],diff_phi_v[sh+8]
  //     };
  //     CHKERR moab_ref.tag_set_data(th_grad,&*nit,1,grad);
  //   }
  // }

  // cout << "NBFACETRI_AINSWORTH_HCURL(order) " <<
  // NBFACETRI_AINSWORTH_HCURL(order) << endl; MatrixDouble base_face_functions(
  //   4,3*NBFACETRI_AINSWORTH_HCURL(order)*nb_gauss_pts
  // );
  // MatrixDouble diff_base_face_functions(
  //   4,9*NBFACETRI_AINSWORTH_HCURL(order)*nb_gauss_pts
  // );
  // for(int ff=0;ff!=4;ff++) {
  //   phi_f[ff] = &base_face_functions(ff,0);
  //   diff_phi_f[ff] = &diff_base_face_functions(ff,0);
  // }
  // CHKERR  Hcurl_Ainsworth_FaceFunctions_MBTET(
  //   faces_nodes,
  //   faces_order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   phi_f,
  //   diff_phi_f,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  // for(int ff = 0;ff!=4;ff++) {
  //   for(int ll = 0;ll!=NBFACETRI_AINSWORTH_HCURL(order);ll++) {
  //     std::ostringstream ss;
  //     ss << "curl_face_" << ff << "_" << ll;
  //     Tag th;
  //     CHKERR moab_ref.tag_get_handle(
  //       ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //     );
  //
  //     int gg = 0;
  //     for(Range::iterator nit =
  //     elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //       int idx = 3*NBFACETRI_AINSWORTH_HCURL   (order)*gg+ll*3;
  //       CHKERR moab_ref.tag_set_data(th,&*nit,1,&(phi_f[ff][idx]));
  //
  //     }
  //   }
  // }
  //
  // cout << "NBVOLUMETET_AINSWORTH_TET_HCURL(order) " <<
  // NBVOLUMETET_AINSWORTH_HCURL(order) << endl; VectorDouble
  // base_volume_functions(3*NBVOLUMETET_AINSWORTH_HCURL(order)*nb_gauss_pts);
  // VectorDouble
  // diff_base_volume_functions(9*NBVOLUMETET_AINSWORTH_HCURL(order)*nb_gauss_pts);
  // phi_v = &base_volume_functions[0]; diff_phi_v =
  // &diff_base_volume_functions[0]; CHKERR
  // MoFEM::Hcurl_Ainsworth_VolumeFunctions_MBTET(
  //   order,
  //   &*shape_fun.data().begin(),
  //   diff_shape_fun,
  //   phi_v,
  //   diff_phi_v,
  //   nb_gauss_pts,
  //   Legendre_polynomials
  // );
  // for(int ll = 0;ll!=NBVOLUMETET_AINSWORTH_HCURL(order);ll++) {
  //   std::ostringstream ss;
  //   ss << "curl_volume_" << ll;
  //   Tag th;
  //   CHKERR moab_ref.tag_get_handle(
  //     ss.str().c_str(),3,MB_TYPE_DOUBLE,th,MB_TAG_CREAT|MB_TAG_SPARSE,def_val
  //   );
  //
  //   int gg = 0;
  //   for(Range::iterator nit =
  //   elem_nodes.begin();nit!=elem_nodes.end();nit++,gg++) {
  //     int idx = 3*NBVOLUMETET_AINSWORTH_HCURL(order)*gg+ll*3;
  //     CHKERR moab_ref.tag_set_data(th,&*nit,1,&(phi_v[idx]));
  //
  //   }
  // }

  EntityHandle meshset;
  CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
  CHKERR moab_ref.add_entities(meshset, tets);
  CHKERR moab_ref.write_file(file_name.c_str(), "VTK", "", &meshset, 1);

  MoFEMFunctionReturnHot(0);
}

#endif // GENERATE_VTK_WITH_CURL_BASE

#ifndef GENERATE_VTK_WITH_CURL_BASE

struct HcurlEdgeBase {

  FTensor::Index<'i', 3> i;

  FTensor::Tensor1<double, 3> tGradN0pN1;
  FTensor::Tensor1<double, 3> tPhi0;
  FTensor::Tensor1<double, 3> tDiffb;

  VectorDouble fI;
  MatrixDouble diffFi;

  template <int DIM, bool CALCULATE_DIRVATIVES>
  MoFEMErrorCode
  calculate(int p, int nb_integration_pts, int n0_idx, int n1_idx, double n[],
            FTensor::Tensor1<double, 3> t_grad_n[],
            FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> &t_phi,
            FTensor::Tensor2<FTensor::PackPtr<double *, 3 * DIM>, 3, DIM>
                *t_diff_phi_ptr) {

    FTensor::Index<'j', DIM> j;

    MoFEMFunctionBegin;

    FTensor::Tensor1<double, 3> &t_grad_n0 = t_grad_n[n0_idx];
    FTensor::Tensor1<double, 3> &t_grad_n1 = t_grad_n[n1_idx];
    tGradN0pN1(i) = t_grad_n0(i) + t_grad_n1(i);

    fI.resize(p + 1);
    diffFi.resize(3, p + 1);
    diffFi.clear();

    FTensor::Tensor2<double, 3, DIM> t_diff_phi0;

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      const int shift_n = (DIM + 1) * gg;
      const double n0 = n[shift_n + n0_idx];
      const double n1 = n[shift_n + n1_idx];

      tPhi0(i) = n0 * t_grad_n1(i) - n1 * t_grad_n0(i);
      t_phi(i) = tPhi0(i);

      ++t_phi;

      if (CALCULATE_DIRVATIVES) {

        t_diff_phi0(i, j) =
            t_grad_n0(j) * t_grad_n1(i) - t_grad_n1(j) * t_grad_n0(i);
        (*t_diff_phi_ptr)(i, j) = t_diff_phi0(i, j);
        ++(*t_diff_phi_ptr);
      }

      if (p > 1) {

        if (CALCULATE_DIRVATIVES)
          CHKERR Jacobi_polynomials(p, 0, n1, n0 + n1, &t_grad_n1(0),
                                    &tGradN0pN1(0), &*fI.data().begin(),
                                    &*diffFi.data().begin(), DIM);
        else
          CHKERR Jacobi_polynomials(p, 0, n1, n0 + n1, nullptr, nullptr,
                                    &*fI.data().begin(), nullptr, DIM);

        FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3> t_diff_fi(
            &diffFi(0, 1), &diffFi(1, 1), &diffFi(2, 1));

        for (int oo = 1; oo <= p - 1; ++oo) {

          const double b = pow(n0 + n1, oo);
          t_phi(i) = b * fI[oo] * tPhi0(i);

          if (CALCULATE_DIRVATIVES) {

            tDiffb(i) =
                oo * pow(n0 + n1, oo - 1) * (t_grad_n0(i) + t_grad_n1(i));
            (*t_diff_phi_ptr)(i, j) = (b * fI[oo]) * t_diff_phi0(i, j) +
                                      (b * t_diff_fi(j)) * tPhi0(i) +
                                      tDiffb(j) * fI[oo] * tPhi0(i);
            ++t_diff_fi;
            ++(*t_diff_phi_ptr);
          }

          ++t_phi;
        }
      }
    }

    MoFEMFunctionReturn(0);
  }
};

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_EdgeBaseFunctions_MBTET(
    int *sense, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts) {

  constexpr int e_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                 {0, 3}, {1, 3}, {2, 3}};

  MoFEMFunctionBegin;

  FTensor::Tensor1<double, 3> t_grad_n[4];
  for (int nn = 0; nn != 4; ++nn)
    t_grad_n[nn] = FTensor::Tensor1<double, 3>(
        diff_n[3 * nn + 0], diff_n[3 * nn + 1], diff_n[3 * nn + 2]);

  HcurlEdgeBase h_curl_base_on_edge;

  for (int ee = 0; ee != 6; ++ee) {

    auto t_phi = getFTensor1FromPtr<3>(phi[ee]);
    auto t_diff_phi = getFTensor2FromPtr<3, 3>(diff_phi[ee]);

    int n0_idx = e_nodes[ee][0];
    int n1_idx = e_nodes[ee][1];
    if (sense[ee] == -1) {
      int n_tmp = n0_idx;
      n0_idx = n1_idx;
      n1_idx = n_tmp;
    }

    CHKERR h_curl_base_on_edge.calculate<3, true>(p[ee], nb_integration_pts,
                                                  n0_idx, n1_idx, n, t_grad_n,
                                                  t_phi, &t_diff_phi);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_EdgeBaseFunctions_MBTRI(
    int *sense, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts) {

  constexpr int e_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};

  MoFEMFunctionBegin;

  FTensor::Tensor1<double, 3> t_grad_n[3];
  for (int nn = 0; nn != 3; ++nn)
    t_grad_n[nn] =
        FTensor::Tensor1<double, 3>(diff_n[2 * nn + 0], diff_n[2 * nn + 1], 0.);

  HcurlEdgeBase h_curl_base_on_edge;

  for (int ee = 0; ee != 3; ++ee) {

    if (p[ee]) {

      auto t_phi = getFTensor1FromPtr<3>(phi[ee]);
      auto t_diff_phi = getFTensor2FromPtr<3, 2>(diff_phi[ee]);

      int n0_idx = e_nodes[ee][0];
      int n1_idx = e_nodes[ee][1];
      if (sense[ee] == -1) {
        int n_tmp = n0_idx;
        n0_idx = n1_idx;
        n1_idx = n_tmp;
      }

      CHKERR h_curl_base_on_edge.calculate<2, true>(p[ee], nb_integration_pts,
                                                    n0_idx, n1_idx, n, t_grad_n,
                                                    t_phi, &t_diff_phi);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_EdgeBaseFunctions_MBEDGE(
    int sense, int p, double *n, double *diff_n, double *phi, double *diff_phi,
    int nb_integration_pts) {
  MoFEMFunctionBegin;

  FTensor::Tensor1<double, 3> t_grad_n[2];
  for (int nn = 0; nn != 2; ++nn)
    t_grad_n[nn] = FTensor::Tensor1<double, 3>(diff_n[nn], 0., 0.);

  HcurlEdgeBase h_curl_base_on_edge;

  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_phi(
      &phi[HVEC0], &phi[HVEC1], &phi[HVEC2]);

  if (diff_phi != NULL)
    SETERRQ(PETSC_COMM_WORLD, MOFEM_DATA_INCONSISTENCY,
            "Not implemented direvatives for edge for Hcurl Demkowicz base");

  int n0_idx = 0;
  int n1_idx = 1;
  if (sense == -1) {
    int n_tmp = n0_idx;
    n0_idx = n1_idx;
    n1_idx = n_tmp;
  }

  CHKERR h_curl_base_on_edge.calculate<1, false>(
      p, nb_integration_pts, n0_idx, n1_idx, n, t_grad_n, t_phi, nullptr);

  MoFEMFunctionReturn(0);
}

struct HcurlFaceBase {
  HcurlEdgeBase hCurlBaseOnEdge;
  VectorDouble f0PhiII, diffF0PhiII;
  VectorDouble f1PhiII, diffF1PhiII;
  VectorDouble iFiF0, diffIFiF0;
  VectorDouble iFiF1, diffIFiF1;

  FTensor::Index<'i', 3> i;

  template <int DIM>
  MoFEMErrorCode calculateOneFamily(
      int p, int nb_integration_pts, int n0f0_idx, int n1f0_idx, int n2f0_idx,
      double n[], FTensor::Tensor1<double, 3> t_grad_n[],
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> &t_phi,
      FTensor::Tensor2<FTensor::PackPtr<double *, DIM * 3>, 3, DIM>
          &t_diff_phi) {

    FTensor::Index<'j', DIM> j;

    MoFEMFunctionBegin;
    f0PhiII.resize(3 * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts, false);
    diffF0PhiII.resize(3 * DIM * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts,
                       false);

    // edge base for family I
    double *f0_phi_ii = &*f0PhiII.data().begin();
    double *diff_f0_phi_ii = &*diffF0PhiII.data().begin();
    auto t_f0_phi_ii = getFTensor1FromPtr<3>(f0_phi_ii);
    auto t_diff_f0_phi_ii = getFTensor2FromPtr<3, DIM>(diff_f0_phi_ii);

    CHKERR hCurlBaseOnEdge.calculate<DIM, true>(p - 1, nb_integration_pts,
                                                n0f0_idx, n1f0_idx, n, t_grad_n,
                                                t_f0_phi_ii, &t_diff_f0_phi_ii);

    FTensor::Tensor1<double, 3> &t_grad_n0f0 = t_grad_n[n0f0_idx];
    FTensor::Tensor1<double, 3> &t_grad_n1f0 = t_grad_n[n1f0_idx];
    FTensor::Tensor1<double, 3> &t_grad_n2f0 = t_grad_n[n2f0_idx];
    FTensor::Tensor1<double, 3> t_grad_n0f0_p_n1f0;
    t_grad_n0f0_p_n1f0(i) = t_grad_n0f0(i) + t_grad_n1f0(i) + t_grad_n2f0(i);

    iFiF0.resize(p + 1, false);
    diffIFiF0.resize(3 * p + 3, false);
    diffIFiF0.clear();

    double *ifif0 = &*iFiF0.data().begin();
    double *diff_ifif0 = &*diffIFiF0.data().begin();

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      const int shift_n = (DIM + 1) * gg;
      const double n0f0 = n[shift_n + n0f0_idx];
      const double n1f0 = n[shift_n + n1f0_idx];
      const double n2f0 = n[shift_n + n2f0_idx];

      int phi_shift = 3 * NBEDGE_DEMKOWICZ_HCURL(p - 1) * gg;
      int diff_phi_shift = (3 * DIM) * NBEDGE_DEMKOWICZ_HCURL(p - 1) * gg;

      for (int oo = 2; oo <= p; ++oo) {

        auto t_f0_phi_ii = getFTensor1FromPtr<3>(&f0_phi_ii[phi_shift]);
        auto t_diff_f0_phi_ii =
            getFTensor2FromPtr<3, DIM>(&diff_f0_phi_ii[diff_phi_shift]);

        for (int ii = 0; ii <= oo - 2; ii++) {

          int jj = oo - 2 - ii;

          // family I
          CHKERR IntegratedJacobi_polynomials(
              jj + 1, 2 * ii + 1, n2f0, n0f0 + n1f0 + n2f0, &t_grad_n2f0(0),
              &t_grad_n0f0_p_n1f0(0), ifif0, diff_ifif0, DIM);
          FTensor::Tensor1<double, 3> t_diff_ifif0(
              diff_ifif0[0 + jj], diff_ifif0[(jj + 1) + jj],
              diff_ifif0[2 * (jj + 1) + jj]);
          t_phi(i) = ifif0[jj] * t_f0_phi_ii(i);
          t_diff_phi(i, j) = ifif0[jj] * t_diff_f0_phi_ii(i, j) +
                             t_diff_ifif0(j) * t_f0_phi_ii(i);
          ++t_phi;
          ++t_diff_phi;
          ++t_f0_phi_ii;
          ++t_diff_f0_phi_ii;
        }
      }
    }

    MoFEMFunctionReturn(0);
  }

  template <int DIM>
  MoFEMErrorCode
  calculateTwoFamily(int p, int nb_integration_pts, int n0f0_idx, int n1f0_idx,
                     int n2f0_idx, int n0f1_idx, int n1f1_idx, int n2f1_idx,
                     double n[], FTensor::Tensor1<double, 3> t_grad_n[],
                     FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> &t_phi,
                     FTensor::Tensor2<FTensor::PackPtr<double *, 3 * DIM>, 3,
                                      DIM> &t_diff_phi) {

    FTensor::Index<'j', DIM> j;

    MoFEMFunctionBegin;

    f0PhiII.resize(3 * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts, false);
    diffF0PhiII.resize(3 * DIM * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts,
                       false);
    f1PhiII.resize(3 * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts, false);
    diffF1PhiII.resize(3 * DIM * NBEDGE_DEMKOWICZ_HCURL(p) * nb_integration_pts,
                       false);

    // edge base for family I
    double *f0_phi_ii = &*f0PhiII.data().begin();
    double *diff_f0_phi_ii = &*diffF0PhiII.data().begin();
    auto t_f0_phi_ii = getFTensor1FromPtr<3>(f0_phi_ii);
    auto t_diff_f0_phi_ii = getFTensor2FromPtr<3, DIM>(diff_f0_phi_ii);
    CHKERR hCurlBaseOnEdge.calculate<DIM, true>(p - 1, nb_integration_pts,
                                                n0f0_idx, n1f0_idx, n, t_grad_n,
                                                t_f0_phi_ii, &t_diff_f0_phi_ii);

    // edge base for family II
    double *f1_phi_ii = &*f1PhiII.data().begin();
    double *diff_f1_phi_ii = &*diffF1PhiII.data().begin();
    auto t_f1_phi_ii = getFTensor1FromPtr<3>(f1_phi_ii);
    auto t_diff_f1_phi_ii = getFTensor2FromPtr<3, DIM>(diff_f1_phi_ii);
    CHKERR hCurlBaseOnEdge.calculate<DIM, true>(p - 1, nb_integration_pts,
                                                n0f1_idx, n1f1_idx, n, t_grad_n,
                                                t_f1_phi_ii, &t_diff_f1_phi_ii);

    FTensor::Tensor1<double, 3> &t_grad_n0f0 = t_grad_n[n0f0_idx];
    FTensor::Tensor1<double, 3> &t_grad_n1f0 = t_grad_n[n1f0_idx];
    FTensor::Tensor1<double, 3> &t_grad_n2f0 = t_grad_n[n2f0_idx];
    FTensor::Tensor1<double, 3> t_grad_n0f0_p_n1f0;
    t_grad_n0f0_p_n1f0(i) = t_grad_n0f0(i) + t_grad_n1f0(i);

    FTensor::Tensor1<double, 3> &t_grad_n0f1 = t_grad_n[n0f1_idx];
    FTensor::Tensor1<double, 3> &t_grad_n1f1 = t_grad_n[n1f1_idx];
    FTensor::Tensor1<double, 3> &t_grad_n2f1 = t_grad_n[n2f1_idx];
    FTensor::Tensor1<double, 3> t_grad_n0f1_p_n1f1;
    t_grad_n0f1_p_n1f1(i) = t_grad_n0f1(i) + t_grad_n1f1(i);

    iFiF0.resize(p + 1, false);
    diffIFiF0.resize(3 * p + 3, false);
    diffIFiF0.clear();
    double *ifif0 = &*iFiF0.data().begin();
    double *diff_ifif0 = &*diffIFiF0.data().begin();
    iFiF1.resize(p + 1, false);
    diffIFiF1.resize(3 * p + 3, false);
    diffIFiF1.clear();
    double *ifif1 = &*iFiF1.data().begin();
    double *diff_ifif1 = &*diffIFiF1.data().begin();

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      const int shift_n = (DIM + 1) * gg;
      const double n0f0 = n[shift_n + n0f0_idx];
      const double n1f0 = n[shift_n + n1f0_idx];
      const double n2f0 = n[shift_n + n2f0_idx];
      const double n0f1 = n[shift_n + n0f1_idx];
      const double n1f1 = n[shift_n + n1f1_idx];
      const double n2f1 = n[shift_n + n2f1_idx];

      int phi_shift = 3 * NBEDGE_DEMKOWICZ_HCURL(p - 1) * gg;
      int diff_phi_shift = 3 * DIM * NBEDGE_DEMKOWICZ_HCURL(p - 1) * gg;

      int kk = 0;
      for (int oo = 2; oo <= p; ++oo) {

        auto t_f0_phi_ii = getFTensor1FromPtr<3>(&f0_phi_ii[phi_shift]);
        auto t_diff_f0_phi_ii =
            getFTensor2FromPtr<3, DIM>(&diff_f0_phi_ii[diff_phi_shift]);
        auto t_f1_phi_ii = getFTensor1FromPtr<3>(&f1_phi_ii[phi_shift]);
        auto t_diff_f1_phi_ii =
            getFTensor2FromPtr<3, DIM>(&diff_f1_phi_ii[diff_phi_shift]);

        for (int ii = 0; ii <= oo - 2; ii++) {

          int jj = oo - 2 - ii;

          // family I
          CHKERR IntegratedJacobi_polynomials(
              jj + 1, 2 * ii + 1, n2f0, n0f0 + n1f0, &t_grad_n2f0(0),
              &t_grad_n0f0_p_n1f0(0), ifif0, diff_ifif0, 3);
          FTensor::Tensor1<double, 3> t_diff_ifif0(
              diff_ifif0[0 + jj], diff_ifif0[(jj + 1) + jj],
              diff_ifif0[2 * (jj + 1) + jj]);

          t_phi(i) = ifif0[jj] * t_f0_phi_ii(i);
          t_diff_phi(i, j) = ifif0[jj] * t_diff_f0_phi_ii(i, j) +
                             t_diff_ifif0(j) * t_f0_phi_ii(i);

          ++t_phi;
          ++t_diff_phi;
          ++t_f0_phi_ii;
          ++t_diff_f0_phi_ii;
          ++kk;

          // family II
          CHKERR IntegratedJacobi_polynomials(
              jj + 1, 2 * ii + 1, n2f1, n0f1 + n1f1, &t_grad_n2f1(0),
              &t_grad_n0f1_p_n1f1(0), ifif1, diff_ifif1, 3);
          FTensor::Tensor1<double, 3> t_diff_ifif1(
              diff_ifif1[0 + jj], diff_ifif1[(jj + 1) + jj],
              diff_ifif1[2 * (jj + 1) + jj]);
          t_phi(i) = ifif1[jj] * t_f1_phi_ii(i);
          t_diff_phi(i, j) = ifif1[jj] * t_diff_f1_phi_ii(i, j) +
                             t_diff_ifif1(j) * t_f1_phi_ii(i);
          ++t_phi;
          ++t_diff_phi;
          ++t_f1_phi_ii;
          ++t_diff_f1_phi_ii;
          ++kk;
        }
      }
      if (kk != NBFACETRI_DEMKOWICZ_HCURL(p))
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of base functions");
    }
    MoFEMFunctionReturn(0);
  }
};

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_FaceBaseFunctions_MBTET(
    int *faces_nodes, int *p, double *n, double *diff_n, double *phi[],
    double *diff_phi[], int nb_integration_pts) {
  MoFEMFunctionBegin;

  FTensor::Tensor1<double, 3> t_grad_n[4];
  for (int nn = 0; nn != 4; ++nn) {
    t_grad_n[nn] = FTensor::Tensor1<double, 3>(
        diff_n[3 * nn + 0], diff_n[3 * nn + 1], diff_n[3 * nn + 2]);
  };

  HcurlFaceBase h_curl_face_base;

  for (int ff = 0; ff != 4; ++ff) {

    if (p[ff] > 1) {

      auto t_phi = getFTensor1FromPtr<3>(phi[ff]);
      auto t_diff_phi = getFTensor2FromPtr<3, 3>(diff_phi[ff]);

      // f0, f1 - family I and family II
      const int n0f0_idx = faces_nodes[3 * ff + 0];
      const int n1f0_idx = faces_nodes[3 * ff + 1];
      const int n2f0_idx = faces_nodes[3 * ff + 2];
      // family II
      const int n0f1_idx = faces_nodes[3 * ff + 1];
      const int n1f1_idx = faces_nodes[3 * ff + 2];
      const int n2f1_idx = faces_nodes[3 * ff + 0];

      CHKERR h_curl_face_base.calculateTwoFamily<3>(
          p[ff], nb_integration_pts, n0f0_idx, n1f0_idx, n2f0_idx, n0f1_idx,
          n1f1_idx, n2f1_idx, n, t_grad_n, t_phi, t_diff_phi);
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_FaceBaseFunctions_MBTRI(
    int *faces_nodes, int p, double *n, double *diff_n, double *phi,
    double *diff_phi, int nb_integration_pts) {
  MoFEMFunctionBegin;

  FTensor::Tensor1<double, 3> t_grad_n[3];
  for (int nn = 0; nn != 3; ++nn) {
    t_grad_n[nn] =
        FTensor::Tensor1<double, 3>(diff_n[2 * nn + 0], diff_n[2 * nn + 1], 0.);
  };

  HcurlFaceBase h_curl_face_base;

  if (p > 1) {

    auto t_phi = getFTensor1FromPtr<3>(phi);
    auto t_diff_phi = getFTensor2FromPtr<3, 2>(diff_phi);

    // f0, f1 - family I and family II
    const int n0f0_idx = faces_nodes[0];
    const int n1f0_idx = faces_nodes[1];
    const int n2f0_idx = faces_nodes[2];
    // family II
    const int n0f1_idx = faces_nodes[1];
    const int n1f1_idx = faces_nodes[2];
    const int n2f1_idx = faces_nodes[0];

    CHKERR h_curl_face_base.calculateTwoFamily<2>(
        p, nb_integration_pts, n0f0_idx, n1f0_idx, n2f0_idx, n0f1_idx, n1f1_idx,
        n2f1_idx, n, t_grad_n, t_phi, t_diff_phi);
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode MoFEM::Hcurl_Demkowicz_VolumeBaseFunctions_MBTET(
    int p, double *n, double *diff_n, double *phi, double *diff_phi,
    int nb_integration_pts) {

  constexpr int family[3][4] = {{0, 1, 2, 3}, {1, 2, 3, 0}, {2, 3, 0, 1}};
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  MoFEMFunctionBegin;

  if (p > 2) {

    auto t_phi = getFTensor1FromPtr<3>(phi);
    auto t_diff_phi = getFTensor2FromPtr<3, 3>(diff_phi);

    FTensor::Tensor1<double, 3> t_grad_n[4];
    for (int nn = 0; nn != 4; ++nn) {
      t_grad_n[nn] = FTensor::Tensor1<double, 3>(
          diff_n[3 * nn + 0], diff_n[3 * nn + 1], diff_n[3 * nn + 2]);
    };

    int nb_face_functions = (NBFACETRI_DEMKOWICZ_HCURL(p - 1)) / 2;
    MatrixDouble phi_ij(3, 3 * nb_face_functions * nb_integration_pts);
    MatrixDouble diff_phi_ij(3, 9 * nb_face_functions * nb_integration_pts);
    MatrixDouble fi_k(3, p + 1);
    MatrixDouble diff_fi_k(3, 3 * p + 3);
    HcurlFaceBase h_curl_face_base;

    // calate face base for each family
    for (int ff = 0; ff != 3; ++ff) {
      double *phi_ij_ptr = &phi_ij(ff, 0);
      double *diff_phi_ij_ptr = &diff_phi_ij(ff, 0);

      auto t_phi_ij = getFTensor1FromPtr<3>(phi_ij_ptr);
      auto t_diff_phi_ij = getFTensor2FromPtr<3, 3>(diff_phi_ij_ptr);

      const int n0_idx = family[ff][0];
      const int n1_idx = family[ff][1];
      const int n2_idx = family[ff][2];

      CHKERR h_curl_face_base.calculateOneFamily<3>(
          p - 1, nb_integration_pts, n0_idx, n1_idx, n2_idx, n, t_grad_n,
          t_phi_ij, t_diff_phi_ij);
    }

    FTensor::Tensor1<double, 3> &t_grad_n3f0 = t_grad_n[family[0][3]];
    FTensor::Tensor1<double, 3> &t_grad_n3f1 = t_grad_n[family[1][3]];
    FTensor::Tensor1<double, 3> &t_grad_n3f2 = t_grad_n[family[2][3]];

    FTensor::Tensor1<double, 3> t_sum_f0;
    t_sum_f0(i) = -t_grad_n3f0(i);
    FTensor::Tensor1<double, 3> t_sum_f1;
    t_sum_f1(i) = -t_grad_n3f1(i);
    FTensor::Tensor1<double, 3> t_sum_f2;
    t_sum_f2(i) = -t_grad_n3f2(i);

    for (int gg = 0; gg != nb_integration_pts; ++gg) {

      int shift_n = 4 * gg;

      double n3f0 = n[shift_n + family[0][3]];
      double n3f1 = n[shift_n + family[1][3]];
      double n3f2 = n[shift_n + family[2][3]];

      int kk = 0;
      for (int oo = 3; oo <= p; ++oo) {

        int phi_shift = 3 * nb_face_functions * gg;
        int diff_phi_shift = 9 * nb_face_functions * gg;

        auto t_phi_face_f0 = getFTensor1FromPtr<3>(&phi_ij(0, phi_shift));
        auto t_diff_phi_face_f0 =
            getFTensor2FromPtr<3, 3>(&diff_phi_ij(0, diff_phi_shift));
        auto t_phi_face_f1 = getFTensor1FromPtr<3>(&phi_ij(1, phi_shift));
        auto t_diff_phi_face_f1 =
            getFTensor2FromPtr<3, 3>(&diff_phi_ij(1, diff_phi_shift));
        auto t_phi_face_f2 = getFTensor1FromPtr<3>(&phi_ij(2, phi_shift));
        auto t_diff_phi_face_f2 =
            getFTensor2FromPtr<3, 3>(&diff_phi_ij(2, diff_phi_shift));

        int ij = 0;
        for (int oo_ij = 2; oo_ij != oo; ++oo_ij) {
          int k = oo - oo_ij;

          CHKERR IntegratedJacobi_polynomials(k, 2 * oo_ij, n3f0, 1 - n3f0,
                                              &t_grad_n3f0(0), &t_sum_f0(0),
                                              &fi_k(0, 0), &diff_fi_k(0, 0), 3);
          CHKERR IntegratedJacobi_polynomials(k, 2 * oo_ij, n3f1, 1 - n3f1,
                                              &t_grad_n3f1(0), &t_sum_f1(0),
                                              &fi_k(1, 0), &diff_fi_k(1, 0), 3);
          CHKERR IntegratedJacobi_polynomials(k, 2 * oo_ij, n3f2, 1 - n3f2,
                                              &t_grad_n3f2(0), &t_sum_f2(0),
                                              &fi_k(2, 0), &diff_fi_k(2, 0), 3);

          FTensor::Tensor1<double, 3> t_diff_fi_k_f0(
              diff_fi_k(0, 0 + k - 1), diff_fi_k(0, k + k - 1),
              diff_fi_k(0, 2 * k + k - 1));
          FTensor::Tensor1<double, 3> t_diff_fi_k_f1(
              diff_fi_k(1, 0 + k - 1), diff_fi_k(1, k + k - 1),
              diff_fi_k(1, 2 * k + k - 1));
          FTensor::Tensor1<double, 3> t_diff_fi_k_f2(
              diff_fi_k(2, 0 + k - 1), diff_fi_k(2, k + k - 1),
              diff_fi_k(2, 2 * k + k - 1));

          for (; ij != NBFACETRI_DEMKOWICZ_HCURL(oo_ij) / 2; ++ij) {
            t_phi(i) = fi_k(0, k - 1) * t_phi_face_f0(i);
            t_diff_phi(i, j) = t_diff_fi_k_f0(j) * t_phi_face_f0(i) +
                               fi_k(0, k - 1) * t_diff_phi_face_f0(i, j);
            ++t_phi;
            ++t_diff_phi;
            ++t_phi_face_f0;
            ++t_diff_phi_face_f0;
            ++kk;

            t_phi(i) = fi_k(1, k - 1) * t_phi_face_f1(i);
            t_diff_phi(i, j) = t_diff_fi_k_f1(j) * t_phi_face_f1(i) +
                               fi_k(1, k - 1) * t_diff_phi_face_f1(i, j);
            ++t_phi;
            ++t_diff_phi;
            ++t_phi_face_f1;
            ++t_diff_phi_face_f1;
            ++kk;

            t_phi(i) = fi_k(2, k - 1) * t_phi_face_f2(i);
            t_diff_phi(i, j) = t_diff_fi_k_f2(j) * t_phi_face_f2(i) +
                               fi_k(2, k - 1) * t_diff_phi_face_f2(i, j);
            ++t_phi;
            ++t_diff_phi;
            ++t_phi_face_f2;
            ++t_diff_phi_face_f2;
            ++kk;
          }
        }
      }
      if (kk != NBVOLUMETET_DEMKOWICZ_HCURL(p))
        SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                "Wrong number of base functions");
    }
  }

  MoFEMFunctionReturn(0);
}

#endif

#ifdef GENERATE_VTK_WITH_CURL_BASE

MoFEMErrorCode VTK_Demkowicz_Hcurl_MBTET(const string file_name) {
  MoFEMFunctionBegin;

  double base_coords[] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  moab::Core core_ref;
  moab::Interface &moab_ref = core_ref;

  EntityHandle nodes[4];
  for (int nn = 0; nn < 4; nn++) {
    CHKERR moab_ref.create_vertex(&base_coords[3 * nn], nodes[nn]);
  }
  EntityHandle tet;
  CHKERR moab_ref.create_element(MBTET, nodes, 4, tet);

  MoFEM::Core m_core_ref(moab_ref, PETSC_COMM_SELF, -2);
  MoFEM::Interface &m_field_ref = m_core_ref;

  CHKERR m_field_ref.getInterface<BitRefManager>()->setBitRefLevelByDim(
      0, 3, BitRefLevel().set(0));

  const int max_level = 3;
  for (int ll = 0; ll != max_level; ll++) {
    Range edges;
    CHKERR m_field_ref.getInterface<BitRefManager>()
        ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                       BitRefLevel().set(), MBEDGE, edges);
    Range tets;
    CHKERR m_field_ref.getInterface<BitRefManager>()
        ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(ll),
                                       BitRefLevel(ll).set(), MBTET, tets);
    // refine mesh
    MeshRefinement *m_ref;
    CHKERR m_field_ref.getInterface(m_ref);
    CHKERR m_ref->add_vertices_in_the_middle_of_edges(
        edges, BitRefLevel().set(ll + 1));
    CHKERR m_ref->refine_TET(tets, BitRefLevel().set(ll + 1));
  }

  Range tets;
  CHKERR m_field_ref.getInterface<BitRefManager>()
      ->getEntitiesByTypeAndRefLevel(BitRefLevel().set(max_level),
                                     BitRefLevel().set(max_level), MBTET, tets);

  // Use 10 node tets to print base
  if (1) {
    EntityHandle meshset;
    CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
    CHKERR moab_ref.add_entities(meshset, tets);
    CHKERR moab_ref.convert_entities(meshset, true, false, false);
    CHKERR moab_ref.delete_entities(&meshset, 1);
  }

  Range elem_nodes;
  CHKERR moab_ref.get_connectivity(tets, elem_nodes, false);

  const int nb_gauss_pts = elem_nodes.size();
  MatrixDouble gauss_pts(nb_gauss_pts, 4);
  gauss_pts.clear();
  Range::iterator nit = elem_nodes.begin();
  for (int gg = 0; nit != elem_nodes.end(); nit++, gg++) {
    CHKERR moab_ref.get_coords(&*nit, 1, &gauss_pts(gg, 0));
  }
  gauss_pts = trans(gauss_pts);

  MatrixDouble shape_fun;
  shape_fun.resize(nb_gauss_pts, 4);
  CHKERR ShapeMBTET(&*shape_fun.data().begin(), &gauss_pts(0, 0),
                    &gauss_pts(1, 0), &gauss_pts(2, 0), nb_gauss_pts);

  double diff_shape_fun[12];
  CHKERR ShapeDiffMBTET(diff_shape_fun);

  int edge_sense[6] = {1, 1, 1, 1, 1, 1};
  const int order = 4;
  int edge_order[6] = {order, order, order, order, order, order};

  MatrixDouble edge_phi(6, 3 * NBEDGE_DEMKOWICZ_HCURL(order) * nb_gauss_pts);
  MatrixDouble edge_diff_phi(6,
                             9 * NBEDGE_DEMKOWICZ_HCURL(order) * nb_gauss_pts);

  edge_phi.clear();
  edge_diff_phi.clear();

  double *edge_phi_ptr[] = {&edge_phi(0, 0), &edge_phi(1, 0), &edge_phi(2, 0),
                            &edge_phi(3, 0), &edge_phi(4, 0), &edge_phi(5, 0)};
  double *edge_diff_phi_ptr[] = {&edge_diff_phi(0, 0), &edge_diff_phi(1, 0),
                                 &edge_diff_phi(2, 0), &edge_diff_phi(3, 0),
                                 &edge_diff_phi(4, 0), &edge_diff_phi(5, 0)};

  CHKERR Hcurl_Demkowicz_EdgeBaseFunctions_MBTET(
      edge_sense, edge_order, &*shape_fun.data().begin(), diff_shape_fun,
      edge_phi_ptr, edge_diff_phi_ptr, nb_gauss_pts);

  // cerr << edge_phi << endl;

  for (int ee = 0; ee != 6; ++ee) {
    for (int ll = 0; ll != NBEDGE_DEMKOWICZ_HCURL(order); ++ll) {
      double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      std::string tag_name = "E" + boost::lexical_cast<std::string>(ee) + "_" +
                             boost::lexical_cast<std::string>(ll);
      Tag th;
      CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

      int gg = 0;
      for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
           nit++, gg++) {
        int idx = 3 * NBEDGE_DEMKOWICZ_HCURL(order) * gg + 3 * ll;
        CHKERR moab_ref.tag_set_data(th, &*nit, 1, &edge_phi(ee, idx));
      }
    }
  }

  int faces_order[] = {order, order, order, order};
  int faces_nodes[] = {0, 1, 3, 1, 2, 3, 0, 2, 3, 0, 1, 2};

  MatrixDouble face_phi(4, 3 * NBFACETRI_DEMKOWICZ_HCURL(order) * nb_gauss_pts);
  MatrixDouble face_diff_phi(4, 9 * NBFACETRI_DEMKOWICZ_HCURL(order) *
                                    nb_gauss_pts);
  face_phi.clear();
  face_diff_phi.clear();

  double *face_phi_ptr[] = {&face_phi(0, 0), &face_phi(1, 0), &face_phi(2, 0),
                            &face_phi(3, 0)};
  double *face_diff_phi_ptr[] = {&face_diff_phi(0, 0), &face_diff_phi(1, 0),
                                 &face_diff_phi(2, 0), &face_diff_phi(3, 0)};

  CHKERR Hcurl_Demkowicz_FaceBaseFunctions_MBTET(
      faces_nodes, faces_order, &*shape_fun.data().begin(), diff_shape_fun,
      face_phi_ptr, face_diff_phi_ptr, nb_gauss_pts);

  for (int ff = 0; ff != 4; ++ff) {
    for (int ll = 0; ll != NBFACETRI_DEMKOWICZ_HCURL(order); ++ll) {
      double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      std::string tag_name = "F" + boost::lexical_cast<std::string>(ff) + "_" +
                             boost::lexical_cast<std::string>(ll);
      Tag th;
      CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
                                     MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

      int gg = 0;
      for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
           nit++, gg++) {
        int idx = 3 * NBFACETRI_DEMKOWICZ_HCURL(order) * gg + 3 * ll;
        CHKERR moab_ref.tag_set_data(th, &*nit, 1, &face_phi(ff, idx));
      }
    }
  }

  MatrixDouble vol_phi(nb_gauss_pts, 3 * NBVOLUMETET_DEMKOWICZ_HCURL(order));
  MatrixDouble diff_vol_phi(nb_gauss_pts,
                            9 * NBVOLUMETET_DEMKOWICZ_HCURL(order));

  CHKERR Hcurl_Demkowicz_VolumeBaseFunctions_MBTET(
      order, &*shape_fun.data().begin(), diff_shape_fun, &vol_phi(0, 0),
      &diff_vol_phi(0, 0), nb_gauss_pts);

  for (int ll = 0; ll != NBVOLUMETET_DEMKOWICZ_HCURL(order); ++ll) {
    double def_val[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::string tag_name = "V_" + boost::lexical_cast<std::string>(ll);
    Tag th;
    CHKERR moab_ref.tag_get_handle(tag_name.c_str(), 3, MB_TYPE_DOUBLE, th,
                                   MB_TAG_CREAT | MB_TAG_SPARSE, def_val);

    int gg = 0;
    for (Range::iterator nit = elem_nodes.begin(); nit != elem_nodes.end();
         nit++, gg++) {
      int idx = 3 * ll;
      CHKERR moab_ref.tag_set_data(th, &*nit, 1, &vol_phi(gg, idx));
    }
  }

  EntityHandle meshset;
  CHKERR moab_ref.create_meshset(MESHSET_SET, meshset);
  CHKERR moab_ref.add_entities(meshset, tets);
  CHKERR moab_ref.write_file(file_name.c_str(), "VTK", "", &meshset, 1);

  MoFEMFunctionReturn(0);
}

static char help[] = "...\n\n";

int main(int argc, char *argv[]) {

  MoFEM::Core::Initialize(&argc, &argv, (char *)0, help);

  try {
    // CHKERR
    // VTK_Ainsworth_Hcurl_MBTET("out_curl_vtk_ainsworth_base_on_tet.vtk");
    CHKERR VTK_Demkowicz_Hcurl_MBTET("out_curl_vtk_demkowicz_base_on_tet.vtk");
  }
  CATCH_ERRORS;

  MoFEM::Core::Finalize();

  return 0;
}

#endif // GENERATE_VTK_WITH_CURL_BASE
