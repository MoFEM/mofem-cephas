/** \file Hdiv.cpp

  \brief Implementation of H-curl base

  Based on Hierarchic Finite Element Bases on Unstructured Tetrahedral
  Meshes, by Mark Ainsworth and Joe Coyle
  Shape functions for MBTRI/MBTET and HCurl space

*/

using namespace MoFEM;

boost::function<int(int)> AinsworthOrderHooks::broken_nbfacetri_edge_hdiv =
    [](int p) { return p; };
boost::function<int(int)> AinsworthOrderHooks::broken_nbfacetri_face_hdiv =
    [](int p) { return p; };
boost::function<int(int)> AinsworthOrderHooks::broken_nbvolumetet_edge_hdiv =
    [](int p) { return p; };
boost::function<int(int)> AinsworthOrderHooks::broken_nbvolumetet_face_hdiv =
    [](int p) { return p; };
boost::function<int(int)> AinsworthOrderHooks::broken_nbvolumetet_volume_hdiv =
    [](int p) { return p; };

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f_e[4][3],
    double *diff_phi_f_e[4][3], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;
#ifndef NDEBUG
  if (!diff_phi_f_e)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "expected to return derivatives");
#endif

  for (int ff = 0; ff < 4; ff++) {
    CHKERR Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET_ON_FACE(
        &faces_nodes[3 * ff], p[ff], N, diffN, phi_f_e[ff], diff_phi_f_e[ff],
        gdim, 4, base_polynomials);
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_EdgeFaceShapeFunctions_MBTET_ON_FACE(
    int *faces_nodes, int p, double *N, double *diffN, double *phi_f_e[3],
    double *diff_phi_f_e[3], int gdim, int nb,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  constexpr int face_edges_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  constexpr int face_opposite_edges_node[] = {2, 0, 1};
  FTENSOR_INDEX(3, i);
  FTENSOR_INDEX(3, j);
  FTENSOR_INDEX(3, k);

  MoFEMFunctionBeginHot;
  if (p < 1)
    MoFEMFunctionReturnHot(0);

  FTensor::Tensor1<double, 3> t_edge_cross[3];
  FTensor::Tensor1<double, 3> t_node_diff_ksi[4];
  FTensor::Tensor1<double, 3> t_diff_ksi[3];
  if (diffN) {
    t_node_diff_ksi[0] =
        FTensor::Tensor1<double, 3>(diffN[0], diffN[1], diffN[2]);
    t_node_diff_ksi[1] =
        FTensor::Tensor1<double, 3>(diffN[3], diffN[4], diffN[5]);
    t_node_diff_ksi[2] =
        FTensor::Tensor1<double, 3>(diffN[6], diffN[7], diffN[8]);
    t_node_diff_ksi[3] =
        FTensor::Tensor1<double, 3>(diffN[9], diffN[10], diffN[11]);
    for (int ee = 0; ee < 3; ee++) {
      const int n0 = faces_nodes[face_edges_nodes[ee][0]];
      const int n1 = faces_nodes[face_edges_nodes[ee][1]];
      t_diff_ksi[ee](i) = t_node_diff_ksi[n1](i) - t_node_diff_ksi[n0](i);
      t_edge_cross[ee](i) = levi_civita(i, j, k) * t_node_diff_ksi[n0](j) *
                            t_node_diff_ksi[n1](k);
    }
  } else {
    for (int ee = 0; ee < 3; ee++) {
      t_edge_cross[ee](0) = 1;
      t_edge_cross[ee](1) = 0;
      t_edge_cross[ee](2) = 0;
    }
  }
  double psi_l[p + 1];

  for (int ee = 0; ee != 3; ee++) {
    const int i0 = faces_nodes[face_edges_nodes[ee][0]];
    const int i1 = faces_nodes[face_edges_nodes[ee][1]];
    const int iO = faces_nodes[face_opposite_edges_node[ee]];
    auto t_psi_f_e = getFTensor1FromPtr<3>(phi_f_e[ee]);
    for (int ii = 0; ii != nb * gdim; ii += nb) {
      const double n0 = N[ii + i0];
      const double n1 = N[ii + i1];
      const double lambda = N[ii + iO];
      const double ksi = n1 - n0;
      base_polynomials(p, ksi, NULL, psi_l, NULL, 3);
      auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
      int jj = 0;
      for (int l = 0; l <= p - 1; l++) {
        t_psi_f_e(i) = (lambda * t_psi_l) * t_edge_cross[ee](i);
        ++t_psi_f_e;
        ++t_psi_l;
        ++jj;
      }
#ifndef NDEBUG
      if (jj != NBFACETRI_AINSWORTH_EDGE_HDIV(p)) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", jj, NBFACETRI_AINSWORTH_EDGE_HDIV(p));
      }
#endif
    }
  }

  if (diff_phi_f_e) {
    double diff_psi_l[3 * (p + 1)];
    for (int ee = 0; ee != 3; ee++) {
      const int i0 = faces_nodes[face_edges_nodes[ee][0]];
      const int i1 = faces_nodes[face_edges_nodes[ee][1]];
      const int iO = faces_nodes[face_opposite_edges_node[ee]];
      auto t_diff_phi_f_e = getFTensor2HVecFromPtr<3, 3>(diff_phi_f_e[ee]);
      for (int ii = 0; ii != nb * gdim; ii += nb) {
        const double n0 = N[ii + i0];
        const double n1 = N[ii + i1];
        const double lambda = N[ii + iO];
        const double ksi = n1 - n0;
        base_polynomials(p, ksi, &t_diff_ksi[ee](0), psi_l, diff_psi_l, 3);
        auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
        auto t_diff_psi_l = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>{
            &diff_psi_l[0], &diff_psi_l[p + 1], &diff_psi_l[2 * p + 2]};
        for (int l = 0; l <= p - 1; l++) {
          t_diff_phi_f_e(i, j) =
              (t_node_diff_ksi[iO](j) * t_psi_l + lambda * t_diff_psi_l(j)) *
              t_edge_cross[ee](i);
          ++t_diff_psi_l;
          ++t_diff_phi_f_e;
          ++t_psi_l;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_FaceBubbleShapeFunctions(
    int *faces_nodes, int *p, double *N, double *diffN, double *phi_f[],
    double *diff_phi_f[], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;

#ifndef NDEBUG
  if (!diff_phi_f)
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "expected to return derivatives");
#endif

  for (int ff = 0; ff < 4; ff++) {
    CHKERR Hdiv_Ainsworth_FaceBubbleShapeFunctions_ON_FACE(
        &faces_nodes[3 * ff], p[ff], N, diffN, phi_f[ff], diff_phi_f[ff], gdim,
        4, base_polynomials);
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_FaceBubbleShapeFunctions_ON_FACE(
    int *face_nodes, int p, double *N, double *diffN, double *phi_f,
    double *diff_phi_f, int gdim, int nb,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {
  FTENSOR_INDEX(3, i);
  FTENSOR_INDEX(3, j);
  FTENSOR_INDEX(3, k);

  MoFEMFunctionBeginHot;
  if (p < 3)
    MoFEMFunctionReturnHot(0);

  const int vert_i = face_nodes[1];
  const int vert_j = face_nodes[2];
  const int i0 = face_nodes[0];
  FTensor::Tensor1<double, 3> t_cross;
  FTensor::Tensor1<double, 3> t_node_diff_ksi[4];
  FTensor::Tensor1<double, 3> t_diff_ksi0i;
  FTensor::Tensor1<double, 3> t_diff_ksi0j;

  if (diffN) {
    t_node_diff_ksi[0] =
        FTensor::Tensor1<double, 3>(diffN[0], diffN[1], diffN[2]);
    t_node_diff_ksi[1] =
        FTensor::Tensor1<double, 3>(diffN[3], diffN[4], diffN[5]);
    t_node_diff_ksi[2] =
        FTensor::Tensor1<double, 3>(diffN[6], diffN[7], diffN[8]);
    t_node_diff_ksi[3] =
        FTensor::Tensor1<double, 3>(diffN[9], diffN[10], diffN[11]);
    t_diff_ksi0i(i) = t_node_diff_ksi[vert_i](i) - t_node_diff_ksi[i0](i);
    t_diff_ksi0j(i) = t_node_diff_ksi[vert_j](i) - t_node_diff_ksi[i0](i);
    t_cross(i) = levi_civita(i, j, k) * t_node_diff_ksi[vert_i](j) *
                 t_node_diff_ksi[vert_j](k);
  } else {
    t_cross(0) = 1;
    t_cross(1) = 0;
    t_cross(2) = 0;
  }

  double psi_l[p + 1], diff_psi_l[3 * (p + 1)];
  double psi_m[p + 1], diff_psi_m[3 * (p + 1)];
  auto t_psi_f = getFTensor1FromPtr<3>(phi_f);

  for (int ii = 0; ii < gdim; ii++) {

    const int node_shift = ii * nb;
    const double ni = N[node_shift + vert_i];
    const double nj = N[node_shift + vert_j];
    const double n0 = N[node_shift + i0];
    const double ksi0i = ni - n0;
    const double ksi0j = nj - n0;
    double beta_0ij = n0 * ni * nj;
    base_polynomials(p, ksi0i, NULL, psi_l, NULL, 3);
    base_polynomials(p, ksi0j, NULL, psi_m, NULL, 3);

    int jj = 0;
    int oo = 0;
    for (; oo <= p - 3; oo++) {
      auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
      for (int l = 0; l <= oo; l++) {
        int m = oo - l;
        if (m >= 0) {
          t_psi_f(i) = (beta_0ij * t_psi_l * psi_m[m]) * t_cross(i);
          ++t_psi_f;
        }
        ++t_psi_l;
        ++jj;
      }
    }
#ifndef NDEBUG
    if (jj != NBFACETRI_AINSWORTH_FACE_HDIV(p)) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong order %d != %d", jj, NBFACETRI_AINSWORTH_FACE_HDIV(p));
    }
#endif
  }

  if (diff_phi_f) {
    auto t_diff_phi_f = getFTensor2HVecFromPtr<3, 3>(diff_phi_f);

    for (int ii = 0; ii < gdim; ii++) {

      const int node_shift = ii * nb;
      const double ni = N[node_shift + vert_i];
      const double nj = N[node_shift + vert_j];
      const double n0 = N[node_shift + i0];
      const double ksi0i = ni - n0;
      const double ksi0j = nj - n0;
      double beta_0ij = n0 * ni * nj;
      FTensor::Tensor1<double, 3> t_diff_beta_0ij;
      t_diff_beta_0ij(i) = (ni * nj) * t_node_diff_ksi[i0](i) +
                           (n0 * nj) * t_node_diff_ksi[vert_i](i) +
                           (n0 * ni) * t_node_diff_ksi[vert_j](i);
      base_polynomials(p, ksi0i, &t_diff_ksi0i(0), psi_l, diff_psi_l, 3);
      base_polynomials(p, ksi0j, &t_diff_ksi0j(0), psi_m, diff_psi_m, 3);

      int jj = 0;
      int oo = 0;
      for (; oo <= p - 3; oo++) {
        auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
        auto t_diff_psi_l = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>{
            &diff_psi_l[0], &diff_psi_l[p + 1], &diff_psi_l[2 * p + 2]};
        for (int l = 0; l <= oo; l++) {
          int m = oo - l;
          if (m >= 0) {
            auto t_diff_psi_m =
                FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>{
                    &diff_psi_m[m], &diff_psi_m[p + 1 + m],
                    &diff_psi_m[2 * p + 2 + m]};
            t_diff_phi_f(i, j) = ((t_psi_l * psi_m[m]) * t_diff_beta_0ij(j) +
                                  (beta_0ij * psi_m[m]) * t_diff_psi_l(j) +
                                  (beta_0ij * t_psi_l) * t_diff_psi_m(j)) *
                                 t_cross(i);
            ++t_diff_phi_f;
          }
          ++t_psi_l;
          ++t_diff_psi_l;
          ++jj;
        }
      }
#ifndef NDEBUG
      if (jj != NBFACETRI_AINSWORTH_FACE_HDIV(p)) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", jj, NBFACETRI_AINSWORTH_FACE_HDIV(p));
      }
#endif
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_EdgeBasedVolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v_e[6],
    double *diff_phi_v_e[6], int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  constexpr int edges_nodes[6][2] = {{0, 1}, {1, 2}, {2, 0},
                                     {0, 3}, {1, 3}, {2, 3}};

  MoFEMFunctionBeginHot;
  if (p < 2)
    MoFEMFunctionReturnHot(0);
  if (diffN == NULL) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "data inconsistency");
  }

  FTensor::Tensor1<double, 3> t_coords[4] = {
      FTensor::Tensor1<double, 3>(0., 0., 0.),
      FTensor::Tensor1<double, 3>(1., 0., 0.),
      FTensor::Tensor1<double, 3>(0., 1., 0.),
      FTensor::Tensor1<double, 3>(0., 0., 1.)};
  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};

  FTENSOR_INDEX(3, i);
  FTENSOR_INDEX(3, j);

  FTensor::Tensor1<double, 3> t_tou_e;
  FTensor::Tensor1<double, 3> t_diff_ksi0i;
  FTensor::Tensor1<double, 3> t_diff_beta_e;

  double psi_l[p + 1];
  double diff_psi_l[3 * (p + 1)];

  for (int ee = 0; ee != 6; ee++) {
    t_tou_e(i) =
        t_coords[edges_nodes[ee][1]](i) - t_coords[edges_nodes[ee][0]](i);
    t_diff_ksi0i(i) = t_node_diff_ksi[edges_nodes[ee][1]](i) -
                      t_node_diff_ksi[edges_nodes[ee][0]](i);
    auto t_phi_v_e = getFTensor1FromPtr<3>(phi_v_e[ee]);
    for (int ii = 0; ii != gdim; ii++) {
      const int node_shift = ii * 4;
      const double ni = N[node_shift + edges_nodes[ee][1]];
      const double n0 = N[node_shift + edges_nodes[ee][0]];
      const double beta_e = ni * n0;
      const double ksi0i = ni - n0;
      base_polynomials(p, ksi0i, NULL, psi_l, NULL, 3);
      auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
      for (int l = 0; l <= p - 2; l++) {
        t_phi_v_e(i) = (beta_e * t_psi_l) * t_tou_e(i);
        ++t_phi_v_e;
        ++t_psi_l;
      }
    }
  }

  if (diff_phi_v_e) {
    for (int ee = 0; ee != 6; ee++) {
      t_tou_e(i) =
          t_coords[edges_nodes[ee][1]](i) - t_coords[edges_nodes[ee][0]](i);
      t_diff_ksi0i(i) = t_node_diff_ksi[edges_nodes[ee][1]](i) -
                        t_node_diff_ksi[edges_nodes[ee][0]](i);
      auto t_diff_phi_v_e = getFTensor2HVecFromPtr<3, 3>(diff_phi_v_e[ee]);
      for (int ii = 0; ii != gdim; ii++) {
        const int node_shift = ii * 4;
        const double ni = N[node_shift + edges_nodes[ee][1]];
        const double n0 = N[node_shift + edges_nodes[ee][0]];
        const double beta_e = ni * n0;
        const double ksi0i = ni - n0;
        t_diff_beta_e(i) = ni * t_node_diff_ksi[edges_nodes[ee][0]](i) +
                           t_node_diff_ksi[edges_nodes[ee][1]](i) * n0;
        base_polynomials(p, ksi0i, &t_diff_ksi0i(0), psi_l, diff_psi_l, 3);
        auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
        auto t_diff_psi_l = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
            diff_psi_l, &diff_psi_l[p + 1], &diff_psi_l[2 * p + 2]);
        for (int l = 0; l <= p - 2; l++) {
          t_diff_phi_v_e(i, j) =
              (t_diff_beta_e(j) * t_psi_l + beta_e * t_diff_psi_l(j)) *
              t_tou_e(i);
          ++t_diff_phi_v_e;
          ++t_diff_psi_l;
          ++t_psi_l;
        }
      }
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_FaceBasedVolumeShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v_f[], double *diff_phi_v_f[],
    int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  constexpr int faces_nodes[4][3] = {
      {0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1}};

  MoFEMFunctionBeginHot;
  if (p < 3)
    MoFEMFunctionReturnHot(0);

  FTensor::Tensor1<double, 3> t_coords[4] = {
      FTensor::Tensor1<double, 3>(0., 0., 0.),
      FTensor::Tensor1<double, 3>(1., 0., 0.),
      FTensor::Tensor1<double, 3>(0., 1., 0.),
      FTensor::Tensor1<double, 3>(0., 0., 1.)};

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};

  FTENSOR_INDEX(3, i);
  FTENSOR_INDEX(3, j);

  FTensor::Tensor1<double, 3> t_tau0i[4], t_tau0j[4];
  FTensor::Tensor1<double, 3> t_diff_ksi0i[4], t_diff_ksi0j[4];
  for (int ff = 0; ff != 4; ff++) {
    const int v0 = faces_nodes[ff][0];
    const int vi = faces_nodes[ff][1];
    const int vj = faces_nodes[ff][2];
    t_tau0i[ff](i) = t_coords[vi](i) - t_coords[v0](i);
    t_tau0j[ff](i) = t_coords[vj](i) - t_coords[v0](i);
    t_diff_ksi0i[ff](i) = t_node_diff_ksi[vi](i) - t_node_diff_ksi[v0](i);
    t_diff_ksi0j[ff](i) = t_node_diff_ksi[vj](i) - t_node_diff_ksi[v0](i);
  }

  double psi_l[p + 1], psi_m[p + 1];
  double diff_psi_l[3 * (p + 1)], diff_psi_m[3 * (p + 1)];
  for (int ff = 0; ff != 4; ff++) {
    const int v0 = faces_nodes[ff][0];
    const int vi = faces_nodes[ff][1];
    const int vj = faces_nodes[ff][2];
    auto t_phi_v_f = getFTensor1FromPtr<3>(phi_v_f[ff]);
    auto t_diff_phi_v_f = getFTensor2HVecFromPtr<3, 3>(diff_phi_v_f[ff]);
    for (int ii = 0; ii < gdim; ii++) {
      const int node_shift = 4 * ii;
      const double n0 = N[node_shift + v0];
      const double ni = N[node_shift + vi];
      const double nj = N[node_shift + vj];
      const double beta_f = n0 * ni * nj;
      FTensor::Tensor1<double, 3> t_diff_beta_f;
      t_diff_beta_f(i) = (ni * nj) * t_node_diff_ksi[v0](i) +
                         (n0 * nj) * t_node_diff_ksi[vi](i) +
                         (n0 * ni) * t_node_diff_ksi[vj](i);
      const double ksi0i = ni - n0;
      const double ksi0j = nj - n0;
      base_polynomials(p, ksi0i, &t_diff_ksi0i[ff](0), psi_l, diff_psi_l, 3);
      base_polynomials(p, ksi0j, &t_diff_ksi0j[ff](0), psi_m, diff_psi_m, 3);
      FTensor::Tensor1<double, 3> t_diff_a;
      int jj = 0;
      for (int oo = 0; oo <= p - 3; oo++) {
        auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
        auto t_diff_psi_l = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
            diff_psi_l, &diff_psi_l[p + 1], &diff_psi_l[2 * p + 2]);
        for (int l = 0; l <= oo; l++) {
          int m = oo - l;
          if (m >= 0) {
            auto t_diff_psi_m = FTensor::Tensor1<double, 3>{
                diff_psi_m[m], diff_psi_m[p + 1 + m],
                diff_psi_m[2 * p + 2 + m]};
            const double a = beta_f * t_psi_l * psi_m[m];
            t_phi_v_f(i) = a * t_tau0i[ff](i);
            ++t_phi_v_f;
            ++jj;
            t_phi_v_f(i) = a * t_tau0j[ff](i);
            ++t_phi_v_f;
            ++jj;

            t_diff_a(j) = (t_psi_l * psi_m[m]) * t_diff_beta_f(j) +
                          (beta_f * psi_m[m]) * t_diff_psi_l(j) +
                          (beta_f * t_psi_l) * t_diff_psi_m(j);
            t_diff_phi_v_f(i, j) = t_diff_a(j) * t_tau0i[ff](i);
            ++t_diff_phi_v_f;
            t_diff_phi_v_f(i, j) = t_diff_a(j) * t_tau0j[ff](i);
            ++t_diff_phi_v_f;

            ++t_psi_l;
            ++t_diff_psi_l;
          }
        }
      }
      if (jj != NBVOLUMETET_AINSWORTH_FACE_HDIV(p)) {
        SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
                 "wrong order %d != %d", jj,
                 NBVOLUMETET_AINSWORTH_FACE_HDIV(p));
      }
    }
  }
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Ainsworth_VolumeBubbleShapeFunctions_MBTET(
    int p, double *N, double *diffN, double *phi_v, double *diff_phi_v,
    int gdim,
    PetscErrorCode (*base_polynomials)(int p, double s, double *diff_s,
                                       double *L, double *diffL,
                                       const int dim)) {

  MoFEMFunctionBeginHot;
  if (p < 4)
    MoFEMFunctionReturnHot(0);

  FTensor::Tensor1<double *, 3> t_node_diff_ksi[4] = {
      FTensor::Tensor1<double *, 3>(&diffN[0], &diffN[1], &diffN[2]),
      FTensor::Tensor1<double *, 3>(&diffN[3], &diffN[4], &diffN[5]),
      FTensor::Tensor1<double *, 3>(&diffN[6], &diffN[7], &diffN[8]),
      FTensor::Tensor1<double *, 3>(&diffN[9], &diffN[10], &diffN[11])};

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  FTensor::Number<2> N2;

  FTensor::Tensor1<double, 3> t_diff_ksi0i;
  FTensor::Tensor1<double, 3> t_diff_ksi0j;
  FTensor::Tensor1<double, 3> t_diff_ksi0k;

  t_diff_ksi0i(i) = t_node_diff_ksi[1](i) - t_node_diff_ksi[0](i);
  t_diff_ksi0j(i) = t_node_diff_ksi[2](i) - t_node_diff_ksi[0](i);
  t_diff_ksi0k(i) = t_node_diff_ksi[3](i) - t_node_diff_ksi[0](i);

  double psi_l[p + 1];
  double diff_psi_l[3 * (p + 1)];
  double psi_m[p + 1];
  double diff_psi_m[3 * (p + 1)];
  double psi_n[p + 1];
  double diff_psi_n[3 * (p + 1)];

  auto t_phi_v = getFTensor1FromPtr<3>(phi_v);
  auto t_diff_phi_v = getFTensor2HVecFromPtr<3, 3>(diff_phi_v);

  FTensor::Tensor1<double, 3> t_diff_beta_v;
  for (int ii = 0; ii < gdim; ii++) {
    const int node_shift = ii * 4;
    const double n0 = N[node_shift + 0];
    const double ni = N[node_shift + 1];
    const double nj = N[node_shift + 2];
    const double nk = N[node_shift + 3];
    const double ksi0i = ni - n0;
    const double ksi0j = nj - n0;
    const double ksi0k = nk - n0;
    const double beta_v = n0 * ni * nj * nk;
    t_diff_beta_v(i) = (ni * nj * nk) * t_node_diff_ksi[0](i) +
                       (n0 * nj * nk) * t_node_diff_ksi[1](i) +
                       (n0 * ni * nk) * t_node_diff_ksi[2](i) +
                       (n0 * ni * nj) * t_node_diff_ksi[3](i);
    base_polynomials(p, ksi0i, &t_diff_ksi0i(0), psi_l, diff_psi_l, 3);
    base_polynomials(p, ksi0j, &t_diff_ksi0j(0), psi_m, diff_psi_m, 3);
    base_polynomials(p, ksi0k, &t_diff_ksi0k(0), psi_n, diff_psi_n, 3);

    FTensor::Tensor1<double, 3> t_diff_a;

    int jj = 0;
    for (int oo = 0; oo <= p - 4; oo++) {
      auto t_psi_l = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_l);
      auto t_diff_psi_l = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
          diff_psi_l, &diff_psi_l[p + 1], &diff_psi_l[2 * p + 2]);
      for (int l = 0; l <= oo; l++) {
        auto t_psi_m = FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(psi_m);
        auto t_diff_psi_m = FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
            diff_psi_m, &diff_psi_m[p + 1], &diff_psi_m[2 * p + 2]);
        for (int m = 0; (l + m) <= oo; m++) {
          int n = oo - l - m;
          if (n >= 0) {
            FTensor::Tensor1<double, 3> t_diff_psi_n(diff_psi_n[n],
                                                     diff_psi_n[p + 1 + n],
                                                     diff_psi_n[2 * p + 2 + n]);
            const double a = beta_v * t_psi_l * t_psi_m * psi_n[n];
            t_phi_v(0) = a;
            t_phi_v(1) = 0;
            t_phi_v(2) = 0;
            ++t_phi_v;
            t_phi_v(0) = 0;
            t_phi_v(1) = a;
            t_phi_v(2) = 0;
            ++t_phi_v;
            t_phi_v(0) = 0;
            t_phi_v(1) = 0;
            t_phi_v(2) = a;
            ++t_phi_v;
            t_diff_a(j) = (t_psi_l * t_psi_m * psi_n[n]) * t_diff_beta_v(j) +
                          (beta_v * t_psi_m * psi_n[n]) * t_diff_psi_l(j) +
                          (beta_v * t_psi_l * psi_n[n]) * t_diff_psi_m(j) +
                          (beta_v * t_psi_l * t_psi_m) * t_diff_psi_n(j);
            t_diff_phi_v(N0, j) = t_diff_a(j);
            t_diff_phi_v(N1, j) = 0;
            t_diff_phi_v(N2, j) = 0;
            ++t_diff_phi_v;
            t_diff_phi_v(N0, j) = 0;
            t_diff_phi_v(N1, j) = t_diff_a(j);
            t_diff_phi_v(N2, j) = 0;
            ++t_diff_phi_v;
            t_diff_phi_v(N0, j) = 0;
            t_diff_phi_v(N1, j) = 0;
            t_diff_phi_v(N2, j) = t_diff_a(j);
            ++t_diff_phi_v;
            ++jj;
          }
          ++t_psi_m;
          ++t_diff_psi_m;
        }
        ++t_psi_l;
        ++t_diff_psi_l;
      }
    }

    if (3 * jj != NBVOLUMETET_AINSWORTH_VOLUME_HDIV(p)) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong order %d != %d", jj,
               NBVOLUMETET_AINSWORTH_VOLUME_HDIV(p));
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode
MoFEM::Hdiv_Demkowicz_Face_MBTET_ON_FACE(int *faces_nodes, int p, double *N,
                                         double *diffN, double *phi_f,
                                         double *diff_phi_f, int gdim, int nb) {
  const int face_edges_nodes[3][2] = {{0, 1}, {1, 2}, {2, 0}};
  const int face_opposite_edges_node[] = {2, 0, 1};

  MoFEMFunctionBeginHot;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double, 3> t_cross[3];
  FTensor::Tensor2<double, 3, 3> t_diff_cross(0., 0., 0., 0., 0., 0., 0., 0.,
                                              0.);
  FTensor::Tensor1<double, 3> t_node_diff_ksi[4];
  FTensor::Tensor1<double, 3> t_node_diff_sum_n0_n1;

  const int i0 = faces_nodes[0];
  const int i1 = faces_nodes[1];
  const int i2 = faces_nodes[2];
  const int o[] = {faces_nodes[face_opposite_edges_node[0]],
                   faces_nodes[face_opposite_edges_node[1]],
                   faces_nodes[face_opposite_edges_node[2]]};

  FTensor::Tensor1<double, 3> t_diff_n0_p_n1;
  FTensor::Tensor1<double, 3> t_diff_n0_p_n1_p_n2;

  if (diff_phi_f) {
    t_node_diff_ksi[0] =
        FTensor::Tensor1<double, 3>(diffN[0], diffN[1], diffN[2]);
    t_node_diff_ksi[1] =
        FTensor::Tensor1<double, 3>(diffN[3], diffN[4], diffN[5]);
    t_node_diff_ksi[2] =
        FTensor::Tensor1<double, 3>(diffN[6], diffN[7], diffN[8]);
    t_node_diff_ksi[3] =
        FTensor::Tensor1<double, 3>(diffN[9], diffN[10], diffN[11]);
    t_diff_cross(i, j) = 0;
    for (int ee = 0; ee != 3; ee++) {
      int ei0 = faces_nodes[face_edges_nodes[ee][0]];
      int ei1 = faces_nodes[face_edges_nodes[ee][1]];
      t_cross[ee](0) = t_node_diff_ksi[ei0](1) * t_node_diff_ksi[ei1](2) -
                       t_node_diff_ksi[ei0](2) * t_node_diff_ksi[ei1](1);
      t_cross[ee](1) = t_node_diff_ksi[ei0](2) * t_node_diff_ksi[ei1](0) -
                       t_node_diff_ksi[ei0](0) * t_node_diff_ksi[ei1](2);
      t_cross[ee](2) = t_node_diff_ksi[ei0](0) * t_node_diff_ksi[ei1](1) -
                       t_node_diff_ksi[ei0](1) * t_node_diff_ksi[ei1](0);
      FTensor::Tensor1<double, 3> t_diff_o(
          diffN[3 * o[ee] + 0], diffN[3 * o[ee] + 1], diffN[3 * o[ee] + 2]);
      t_diff_cross(i, j) += t_cross[ee](i) * t_diff_o(j);
      // cerr << t_cross[ee](0) << " " << t_cross[ee](1) << " " <<
      // t_cross[ee](2) << endl;
    }
    // cerr << endl << endl;
    t_diff_n0_p_n1(i) = t_node_diff_ksi[i0](i) + t_node_diff_ksi[i1](i);
    t_diff_n0_p_n1_p_n2(i) = t_diff_n0_p_n1(i) + t_node_diff_ksi[i2](i);
  } else {
    for (int ee = 0; ee != 3; ee++) {
      t_cross[ee](0) = 1;
      t_cross[ee](1) = 0;
      t_cross[ee](2) = 0;
    }
  }

  FTensor::Tensor1<double *, 3> t_phi(&phi_f[HVEC0], &phi_f[HVEC1],
                                      &phi_f[HVEC2], 3);
  boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>> t_diff_phi_ptr;
  if (diff_phi_f) {
    t_diff_phi_ptr = boost::shared_ptr<FTensor::Tensor2<double *, 3, 3>>(
        new FTensor::Tensor2<double *, 3, 3>(
            &diff_phi_f[HVEC0_0], &diff_phi_f[HVEC0_1], &diff_phi_f[HVEC0_2],
            &diff_phi_f[HVEC1_0], &diff_phi_f[HVEC1_1], &diff_phi_f[HVEC1_2],
            &diff_phi_f[HVEC2_0], &diff_phi_f[HVEC2_1], &diff_phi_f[HVEC2_2],
            9));
  }

  double fi[p + 1], diff_fi[3 * p + 3];
  double fj[p + 1], diff_fj[3 * p + 3];
  double tmp_fj[p + 1], tmp_diff_fj[3 * p + 3];
  for (int ii = 0; ii != gdim; ii++) {
    const int shift = ii * nb;
    double n0 = N[shift + i0];
    double n1 = N[shift + i1];
    double n2 = N[shift + i2];
    double *diff_n1 = (diff_phi_f) ? &t_node_diff_ksi[i1](0) : NULL;
    double *diff_n0_p_n1 = (diff_phi_f) ? &t_diff_n0_p_n1(0) : NULL;
    ierr = Jacobi_polynomials(p, 0, n1, n0 + n1, diff_n1, diff_n0_p_n1, fi,
                              diff_phi_f ? diff_fi : NULL, 3);
    CHKERRG(ierr);
    for (int pp = 0; pp <= p; pp++) {
      double *diff_n2 = (diff_phi_f) ? &t_node_diff_ksi[i2](0) : NULL;
      double *diff_n0_p_n1_p_n2 = (diff_phi_f) ? &t_diff_n0_p_n1_p_n2(0) : NULL;
      ierr = Jacobi_polynomials(pp, 2 * pp + 1, n2, n0 + n1 + n2, diff_n2,
                                diff_n0_p_n1_p_n2, tmp_fj,
                                diff_phi_f ? tmp_diff_fj : NULL, 3);
      CHKERRG(ierr);
      fj[pp] = tmp_fj[pp];
      if (diff_phi_f) {
        diff_fj[0 * (p + 1) + pp] = tmp_diff_fj[0 * (pp + 1) + pp];
        diff_fj[1 * (p + 1) + pp] = tmp_diff_fj[1 * (pp + 1) + pp];
        diff_fj[2 * (p + 1) + pp] = tmp_diff_fj[2 * (pp + 1) + pp];
      }
    }
    double no0 = N[shift + o[0]];
    double no1 = N[shift + o[1]];
    double no2 = N[shift + o[2]];
    FTensor::Tensor1<double, 3> base0;
    base0(i) = no0 * t_cross[0](i) + no1 * t_cross[1](i) + no2 * t_cross[2](i);
    int jj = 0;
    for (int oo = 0; oo < p; oo++) {
      FTensor::Tensor0<double *> t_fi(fi);
      FTensor::Tensor1<double *, 3> t_diff_fi(&diff_fi[0], &diff_fi[p + 1],
                                              &diff_fi[2 * p + 2], 1);
      for (int ll = 0; ll <= oo; ll++) {
        const int mm = oo - ll;
        if (mm >= 0) {
          const double a = t_fi * fj[mm];
          // cerr << ll << " " << mm << " " << a << endl;
          t_phi(i) = a * base0(i);
          if (diff_phi_f) {
            FTensor::Tensor1<double *, 3> t_diff_fj(
                &diff_fj[0 + mm], &diff_fj[p + 1 + mm],
                &diff_fj[2 * p + 2 + mm], 1);
            (*t_diff_phi_ptr)(i, j) =
                a * t_diff_cross(i, j) +
                (t_diff_fi(j) * fj[mm] + t_fi * t_diff_fj(j)) * base0(i);
            ++(*t_diff_phi_ptr);
            ++t_diff_fi;
          }
          ++t_fi;
          ++t_phi;
          ++jj;
        }
      }
    }
    if (jj != NBFACETRI_DEMKOWICZ_HDIV(p)) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong number of base functions "
               "jj!=NBFACETRI_DEMKOWICZ_HDIV(p) "
               "%d!=%d",
               jj, NBFACETRI_DEMKOWICZ_HDIV(p));
    }
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode MoFEM::Hdiv_Demkowicz_Interior_MBTET(
    int p, double *N, double *diffN, int p_f[], double *phi_f[4],
    double *diff_phi_f[4], double *phi_v, double *diff_phi_v, int gdim) {

  const int opposite_face_node[4] = {2, 0, 1, 3};
  // list of zero node faces
  const int znf[] = {0, 2, 3};
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;

  FTensor::Tensor1<double, 3> t_node_diff_ksi[4];
  t_node_diff_ksi[0] =
      FTensor::Tensor1<double, 3>(diffN[0], diffN[1], diffN[2]);
  t_node_diff_ksi[1] =
      FTensor::Tensor1<double, 3>(diffN[3], diffN[4], diffN[5]);
  t_node_diff_ksi[2] =
      FTensor::Tensor1<double, 3>(diffN[6], diffN[7], diffN[8]);
  t_node_diff_ksi[3] =
      FTensor::Tensor1<double, 3>(diffN[9], diffN[10], diffN[11]);
  FTensor::Tensor1<double, 3> t_m_node_diff_ksi[4];
  for (int ff = 0; ff != 4; ++ff) {
    t_m_node_diff_ksi[ff](i) = -t_node_diff_ksi[ff](i);
  }

  FTensor::Tensor1<double *, 3> t_phi_v(&phi_v[HVEC0], &phi_v[HVEC1],
                                        &phi_v[HVEC2], 3);
  FTensor::Tensor2<double *, 3, 3> t_diff_phi_v(
      &diff_phi_v[HVEC0_0], &diff_phi_v[HVEC0_1], &diff_phi_v[HVEC0_2],
      &diff_phi_v[HVEC1_0], &diff_phi_v[HVEC1_1], &diff_phi_v[HVEC1_2],
      &diff_phi_v[HVEC2_0], &diff_phi_v[HVEC2_1], &diff_phi_v[HVEC2_2], 9);

  MatrixDouble fk(3, p + 1), diff_fk(3, 3 * p + 3);

  for (int ii = 0; ii != gdim; ii++) {
    const int shift = 4 * ii;

    for (int ff = 0; ff != 3; ff++) {
      const int fff = znf[ff];
      const int iO = opposite_face_node[fff];
      const double nO = N[shift + iO];
      for (int pp = 1; pp <= p; pp++) {
        CHKERR IntegratedJacobi_polynomials(
            pp, 2 * pp + 2, nO, 1 - nO, &t_node_diff_ksi[iO](0),
            &t_m_node_diff_ksi[iO](0), &fk(ff, 0), &diff_fk(ff, 0), 3);
      }
    }

    int jj = 0;
    for (int oo = 2; oo <= p; oo++) {
      for (int k = 1; k != oo; k++) {
        int OO = oo - k;
        if (OO >= 0) {
          int s = NBFACETRI_DEMKOWICZ_HDIV(OO - 1);
          // Note that we do faces 0,2,3, skipping 1. All the faces which have
          // zero node in it.
          int nb_dofs = NBFACETRI_DEMKOWICZ_HDIV(p_f[znf[0]]);
          int sp[] = {ii * 3 * nb_dofs + 3 * s, ii * 3 * nb_dofs + 3 * s,
                      ii * 3 * nb_dofs + 3 * s};
          FTensor::Tensor1<double *, 3> t_phi_f[] = {
              FTensor::Tensor1<double *, 3>(&phi_f[znf[0]][sp[0] + HVEC0],
                                            &phi_f[znf[0]][sp[0] + HVEC1],
                                            &phi_f[znf[0]][sp[0] + HVEC2], 3),
              FTensor::Tensor1<double *, 3>(&phi_f[znf[1]][sp[1] + HVEC0],
                                            &phi_f[znf[1]][sp[1] + HVEC1],
                                            &phi_f[znf[1]][sp[1] + HVEC2], 3),
              FTensor::Tensor1<double *, 3>(&phi_f[znf[2]][sp[2] + HVEC0],
                                            &phi_f[znf[2]][sp[2] + HVEC1],
                                            &phi_f[znf[2]][sp[2] + HVEC2], 3)};
          int sdp[] = {ii * 9 * nb_dofs + 9 * s, ii * 9 * nb_dofs + 9 * s,
                       ii * 9 * nb_dofs + 9 * s};
          FTensor::Tensor2<double *, 3, 3> t_diff_phi_f[] = {
              FTensor::Tensor2<double *, 3, 3>(
                  &diff_phi_f[znf[0]][sdp[0] + HVEC0_0],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC0_1],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC0_2],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC1_0],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC1_1],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC1_2],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC2_0],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC2_1],
                  &diff_phi_f[znf[0]][sdp[0] + HVEC2_2], 9),
              FTensor::Tensor2<double *, 3, 3>(
                  &diff_phi_f[znf[1]][sdp[1] + HVEC0_0],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC0_1],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC0_2],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC1_0],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC1_1],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC1_2],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC2_0],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC2_1],
                  &diff_phi_f[znf[1]][sdp[1] + HVEC2_2], 9),
              FTensor::Tensor2<double *, 3, 3>(
                  &diff_phi_f[znf[2]][sdp[2] + HVEC0_0],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC0_1],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC0_2],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC1_0],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC1_1],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC1_2],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC2_0],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC2_1],
                  &diff_phi_f[znf[2]][sdp[2] + HVEC2_2], 9)};
          for (int ij = s; ij != NBFACETRI_DEMKOWICZ_HDIV(OO); ij++) {
            for (int ff = 0; ff != 3; ff++) {
              FTensor::Tensor1<double, 3> t_diff_fk(diff_fk(ff, 0 * p + k - 1),
                                                    diff_fk(ff, 1 * p + k - 1),
                                                    diff_fk(ff, 2 * p + k - 1));
              t_phi_v(i) = fk(ff, k - 1) * t_phi_f[ff](i);
              t_diff_phi_v(i, j) = t_diff_fk(j) * t_phi_f[ff](i) +
                                   fk(ff, k - 1) * t_diff_phi_f[ff](i, j);
              ++t_phi_v;
              ++t_diff_phi_v;
              ++t_phi_f[ff];
              ++t_diff_phi_f[ff];
              ++jj;
            }
          }
        }
      }
    }
    if (jj != NBVOLUMETET_DEMKOWICZ_HDIV(p)) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong number of base functions "
               "jj!=NBVOLUMETET_DEMKOWICZ_HDIV(p) "
               "%d!=%d",
               jj, NBVOLUMETET_DEMKOWICZ_HDIV(p));
    }
  }

  MoFEMFunctionReturn(0);
}
