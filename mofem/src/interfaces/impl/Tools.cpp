/** \file Tools.cpp
 * \brief Auxiliary tools
 *
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
 * License along with MoFEM. If not, see <http://www.gnu.org/licenses/>
 */

#include <phg-quadrule/quad.h>

namespace MoFEM {

MoFEMErrorCode Tools::query_interface(const MOFEMuuid &uuid,
                                      UnknownInterface **iface) const {
  MoFEMFunctionBeginHot;
  *iface = NULL;
  if (uuid == IDD_MOFEMNodeMerger) {
    *iface = const_cast<Tools *>(this);
    MoFEMFunctionReturnHot(0);
  }
  SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "unknown interface");
  MoFEMFunctionReturnHot(0);
}

double Tools::volumeLengthQuality(const double *coords) {
  double lrms = 0;
  for (int dd = 0; dd != 3; dd++) {
    lrms += pow(coords[0 * 3 + dd] - coords[1 * 3 + dd], 2) +
            pow(coords[0 * 3 + dd] - coords[2 * 3 + dd], 2) +
            pow(coords[0 * 3 + dd] - coords[3 * 3 + dd], 2) +
            pow(coords[1 * 3 + dd] - coords[2 * 3 + dd], 2) +
            pow(coords[1 * 3 + dd] - coords[3 * 3 + dd], 2) +
            pow(coords[2 * 3 + dd] - coords[3 * 3 + dd], 2);
  }
  lrms = sqrt((1. / 6.) * lrms);
  double volume = tetVolume(coords);
  return 6. * sqrt(2.) * volume / pow(lrms, 3);
}

double Tools::tetVolume(const double *coords) {
  double diff_n[12];
  ShapeDiffMBTET(diff_n);
  FTensor::Tensor1<double *, 3> t_diff_n(&diff_n[0], &diff_n[1], &diff_n[2], 3);
  FTensor::Tensor1<const double *, 3> t_coords(&coords[0], &coords[1],
                                               &coords[2], 3);
  FTensor::Tensor2<double, 3, 3> jac;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  jac(i, j) = 0;
  for (int nn = 0; nn != 4; nn++) {
    jac(i, j) += t_coords(i) * t_diff_n(j);
    ++t_coords;
    ++t_diff_n;
  }
  return dEterminant(jac) / 6.;
}

MoFEMErrorCode
Tools::minTetsQuality(const Range &tets, double &min_quality, Tag th,
                      boost::function<double(double, double)> f) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  const EntityHandle *conn;
  int num_nodes;
  double coords[12];
  for (auto tet : tets) {
    CHKERR m_field.get_moab().get_connectivity(tet, conn, num_nodes, true);
    if (th) {
      CHKERR moab.tag_get_data(th, conn, num_nodes, coords);
    } else {
      CHKERR moab.get_coords(conn, num_nodes, coords);
    }
    double q = Tools::volumeLengthQuality(coords);
    if (!std::isnormal(q))
      q = -2;
    min_quality = f(q, min_quality);
  }
  MoFEMFunctionReturn(0);
}

constexpr std::array<double, 2> Tools::diffShapeFunMBEDGE;
constexpr std::array<double, 6> Tools::diffShapeFunMBTRI;
constexpr std::array<double, 12> Tools::diffShapeFunMBTET;
constexpr std::array<double, 4> Tools::shapeFunMBTETAt000;
constexpr std::array<double, 8> Tools::diffShapeFunMBQUADAtCenter;

MoFEMErrorCode Tools::getLocalCoordinatesOnReferenceFourNodeTet(
    const double *elem_coords, const double *global_coords, const int nb_nodes,
    double *local_coords) {
  FTensor::Index<'i', 4> i;
  MoFEMFunctionBeginHot;

  FTensor::Tensor1<FTensor::PackPtr<const double *, 1>, 4> t_elem_coords = {
      &elem_coords[0], &elem_coords[3], &elem_coords[6], &elem_coords[9]};

  FTensor::Tensor1<const double, 4> t_n = {
      shapeFunMBTETAt000[0], shapeFunMBTETAt000[1], shapeFunMBTETAt000[2],
      shapeFunMBTETAt000[3]};
  FTensor::Tensor1<double, 3> t_coords_at_0;

  // Build matrix and get coordinates of zero point
  // ii - global coordinates
  // jj - local direvatives
  MatrixDouble3by3 a(3, 3);
  for (auto ii : {0, 1, 2}) {
    FTensor::Tensor1<FTensor::PackPtr<const double *, 1>, 4> t_diff(
        &diffShapeFunMBTET[0], &diffShapeFunMBTET[3], &diffShapeFunMBTET[6],
        &diffShapeFunMBTET[9]);
    for (auto jj : {0, 1, 2}) {
      a(jj, ii) = t_diff(i) * t_elem_coords(i);
      ++t_diff;
    }
    t_coords_at_0(ii) = t_n(i) * t_elem_coords(i);
    ++t_elem_coords;
  }

  FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_global_coords = {
      &global_coords[0], &global_coords[1], &global_coords[2]};
  FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_local_coords = {
      &local_coords[0], &local_coords[1], &local_coords[2]};

  // Calculate right hand side
  FTensor::Index<'j', 3> j;
  for (int ii = 0; ii != nb_nodes; ++ii) {
    t_local_coords(j) = t_global_coords(j) - t_coords_at_0(j);
    ++t_local_coords;
    ++t_global_coords;
  }

  // Solve problem
  int IPIV[3];
  int info = lapack_dgesv(3, nb_nodes, &a(0, 0), 3, IPIV, local_coords, 3);
  if (info != 0)
    SETERRQ1(PETSC_COMM_SELF, 1, "info == %d", info);

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Tools::getTetsWithQuality(Range &out_tets, const Range &tets,
                                         Tag th,
                                         boost::function<bool(double)> f) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range to_write;
  const EntityHandle *conn;
  int num_nodes;
  double coords[12];
  for (auto tet : tets) {
    CHKERR m_field.get_moab().get_connectivity(tet, conn, num_nodes, true);
    if (th) {
      CHKERR moab.tag_get_data(th, conn, num_nodes, coords);
    } else {
      CHKERR moab.get_coords(conn, num_nodes, coords);
    }
    double q = Tools::volumeLengthQuality(coords);
    if (f(q)) {
      out_tets.insert(tet);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::writeTetsWithQuality(const char *file_name,
                                           const char *file_type,
                                           const char *options,
                                           const Range &tets, Tag th,
                                           boost::function<bool(double)> f) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  Range out_tets;
  CHKERR getTetsWithQuality(out_tets, tets, th, f);
  EntityHandle meshset;
  CHKERR moab.create_meshset(MESHSET_SET, meshset);
  CHKERR moab.add_entities(meshset, out_tets);
  CHKERR moab.write_file(file_name, file_type, options, &meshset, 1);
  CHKERR moab.delete_entities(&meshset, 1);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::checkIfPointIsInTet(const double tet_coords[],
                                          const double global_coord[],
                                          const double tol, bool &result) {
  double loc_coord[] = {0, 0, 0};
  double N[4], diffN[12];
  MoFEMFunctionBegin;
  CHKERR ShapeDiffMBTET(diffN);
  CHKERR ShapeMBTET(N, &loc_coord[0], &loc_coord[1], &loc_coord[2], 1);
  CHKERR ShapeMBTET_inverse(N, diffN, tet_coords, global_coord, loc_coord);
  CHKERR ShapeMBTET(N, &loc_coord[0], &loc_coord[1], &loc_coord[2], 1);
  result = true;
  for (int n = 0; n != 4; ++n) {
    if (N[n] < -tol || (N[n] - 1) > tol) {
      result = false;
      break;
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::checkVectorForNotANumber(const Problem *prb_ptr,
                                               const RowColData row_or_col,
                                               Vec v) {
  MoFEMFunctionBegin;
  int loc_size;
  CHKERR VecGetLocalSize(v, &loc_size);
  int prb_loc_size = 0;
  boost::shared_ptr<NumeredDofEntity_multiIndex> prb_dofs;
  switch (row_or_col) {
  case ROW:
    prb_loc_size = prb_ptr->getNbLocalDofsRow();
    prb_dofs = prb_ptr->getNumeredRowDofsPtr();
    break;
  case COL:
    prb_loc_size = prb_ptr->getNbLocalDofsCol();
    prb_dofs = prb_ptr->getNumeredColDofsPtr();
    break;
    break;
  default:
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Wrong argument, row_or_col should be row or column");
  }
  if (loc_size != prb_loc_size) {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Inconsistent size of vector and problem %d != %d", loc_size,
             prb_loc_size);
  }
  const double *a;
  CHKERR VecGetArrayRead(v, &a);
  MPI_Comm comm = PetscObjectComm((PetscObject)v);
  for (int ii = 0; ii != loc_size; ++ii) {
    if (!boost::math::isfinite(a[ii])) {
      NumeredDofEntityByLocalIdx::iterator dit =
          prb_dofs->get<PetscLocalIdx_mi_tag>().find(ii);
      std::ostringstream ss;
      ss << "Not a number " << a[ii] << " on dof: " << endl
         << **dit << endl
         << endl;
      PetscSynchronizedPrintf(comm, "%s", ss.str().c_str());
    }
  }
  CHKERR VecRestoreArrayRead(v, &a);
  PetscSynchronizedFlush(comm, PETSC_STDOUT);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::getTriNormal(const double *coords, double *normal) {
  MoFEMFunctionBegin;
  double diffN[6];
  CHKERR ShapeDiffMBTRI(diffN);
  CHKERR ShapeFaceNormalMBTRI(diffN, coords, normal);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::getTriNormal(const EntityHandle tri,
                                   double *normal) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  if (moab.type_from_handle(tri) != MBTRI) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "Works only for triangle");
  }
  const EntityHandle *conn;
  int num_nodes;
  double coords[9];
  CHKERR moab.get_connectivity(tri, conn, num_nodes, true);
  CHKERR moab.get_coords(conn, num_nodes, coords);
  CHKERR getTriNormal(coords, normal);
  MoFEMFunctionReturn(0);
}

double Tools::getTriArea(const EntityHandle tri) const {
  FTensor::Tensor1<double, 3> t_normal;
  ierr = getTriNormal(tri, &t_normal(0));
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  FTensor::Index<'i', 3> i;
  return sqrt(t_normal(i) * t_normal(i)) * 0.5;
}

double Tools::getEdgeLength(const double *edge_coords) {
  FTensor::Tensor1<double, 3> t_coords_n0(edge_coords[0], edge_coords[1],
                                          edge_coords[2]);
  FTensor::Tensor1<double, 3> t_coords_n1(edge_coords[3], edge_coords[4],
                                          edge_coords[5]);
  FTensor::Index<'i', 3> i;
  t_coords_n0(i) -= t_coords_n1(i);
  return sqrt(t_coords_n0(i) * t_coords_n0(i));
}

double Tools::getEdgeLength(const EntityHandle edge) {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  auto get_edge_coords = [edge, &moab](double *const coords) {
    MoFEMFunctionBegin;
    if (moab.type_from_handle(edge) != MBEDGE) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_INVALID_DATA, "Works only for edge");
    }
    const EntityHandle *conn;
    int num_nodes;
    CHKERR moab.get_connectivity(edge, conn, num_nodes, true);
    CHKERR moab.get_coords(conn, 2, coords);
    MoFEMFunctionReturn(0);
  };
  double coords[6];
  ierr = get_edge_coords(coords);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  return getEdgeLength(coords);
}

Tools::SEGMENT_MIN_DISTANCE
Tools::minDistancePointFromOnSegment(const double *w_ptr, const double *v_ptr,
                                     const double *p_ptr, double *const t_ptr) {
  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<const double *, 3> t_w(&w_ptr[0], &w_ptr[1], &w_ptr[2]);
  FTensor::Tensor1<const double *, 3> t_v(&v_ptr[0], &v_ptr[1], &v_ptr[2]);
  FTensor::Tensor1<const double *, 3> t_p(&p_ptr[0], &p_ptr[1], &p_ptr[2]);
  FTensor::Tensor1<double, 3> t_vw;
  t_vw(i) = t_v(i) - t_w(i);
  const double dot_vw = t_vw(i) * t_vw(i);
  if (std::abs(dot_vw) < std::numeric_limits<double>::epsilon()) {
    if (t_ptr)
      *t_ptr = 0;
    return SEGMENT_ONE_IS_POINT;
  }
  FTensor::Tensor1<double, 3> t_pw;
  t_pw(i) = t_p(i) - t_w(i);
  const double t = t_pw(i) * t_vw(i) / dot_vw;
  if (t_ptr)
    *t_ptr = t;
  return SOLUTION_EXIST;
}

Tools::SEGMENT_MIN_DISTANCE
Tools::minDistanceFromSegments(const double *w_ptr, const double *v_ptr,
                               const double *k_ptr, const double *l_ptr,
                               double *const tvw_ptr, double *const tlk_ptr) {

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<const double *, 3> t_w(&w_ptr[0], &w_ptr[1], &w_ptr[2]);
  FTensor::Tensor1<const double *, 3> t_v(&v_ptr[0], &v_ptr[1], &v_ptr[2]);
  FTensor::Tensor1<const double *, 3> t_k(&k_ptr[0], &k_ptr[1], &k_ptr[2]);
  FTensor::Tensor1<const double *, 3> t_l(&l_ptr[0], &l_ptr[1], &l_ptr[2]);

  // First segnent is a point
  FTensor::Tensor1<double, 3> t_vw;
  t_vw(i) = t_v(i) - t_w(i);
  double dot_vw = t_vw(i) * t_vw(i);
  if (std::abs(dot_vw) < std::numeric_limits<double>::epsilon()) {
    if (tvw_ptr)
      *tvw_ptr = 0;
    if (minDistancePointFromOnSegment(k_ptr, l_ptr, w_ptr, tlk_ptr) ==
        SEGMENT_ONE_IS_POINT)
      return SEGMENT_TWO_AND_TWO_ARE_POINT;
    else
      return SEGMENT_ONE_IS_POINT;
  }

  // Second segment is a point
  FTensor::Tensor1<double, 3> t_lk;
  t_lk(i) = t_l(i) - t_k(i);
  double dot_lk = t_lk(i) * t_lk(i);
  if (std::abs(dot_lk) < std::numeric_limits<double>::epsilon()) {
    if (tlk_ptr)
      *tlk_ptr = 0;
    if (minDistancePointFromOnSegment(w_ptr, v_ptr, k_ptr, tvw_ptr) ==
        SEGMENT_ONE_IS_POINT)
      return SEGMENT_TWO_AND_TWO_ARE_POINT;
    else
      return SEGMENT_TWO_IS_POINT;
  }

  const double a = t_vw(i) * t_vw(i);
  const double b = -t_vw(i) * t_lk(i);
  const double c = t_lk(i) * t_lk(i);

  const double det = a * c - b * b;
  if (std::abs(det) < std::numeric_limits<double>::epsilon()) {

    return NO_SOLUTION;

  } else {

    FTensor::Tensor1<double, 3> t_wk;
    t_wk(i) = t_w(i) - t_k(i);

    const double ft0 = t_vw(i) * t_wk(i);
    const double ft1 = -t_lk(i) * t_wk(i);
    const double t0 = (ft1 * b - ft0 * c) / det;
    const double t1 = (ft0 * b - ft1 * a) / det;

    if (tvw_ptr)
      *tvw_ptr = t0;
    if (tlk_ptr)
      *tlk_ptr = t1;

    return SOLUTION_EXIST;
  }
}

MoFEMErrorCode Tools::findMinDistanceFromTheEdges(
    const double *v_ptr, const int nb, Range edges, double *min_dist_ptr,
    double *o_ptr, EntityHandle *o_segments) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;

  FTensor::Index<'i', 3> i;

  auto get_point = [i](auto &t_w, auto &t_delta, auto t) {
    FTensor::Tensor1<double, 3> t_p;
    t = std::max(0., std::min(1., t));
    t_p(i) = t_w(i) + t * t_delta(i);
    return t_p;
  };

  auto get_distance = [i](auto &t_p, auto &t_n) {
    FTensor::Tensor1<double, 3> t_dist_vector;
    t_dist_vector(i) = t_p(i) - t_n(i);
    return sqrt(t_dist_vector(i) * t_dist_vector(i));
  };

  for (auto e : edges) {

    int num_nodes;
    const EntityHandle *conn_fixed;
    CHKERR moab.get_connectivity(e, conn_fixed, num_nodes, true);
    VectorDouble6 coords_fixed(6);
    CHKERR moab.get_coords(conn_fixed, num_nodes, &coords_fixed[0]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_f0(
        &coords_fixed[0], &coords_fixed[1], &coords_fixed[2]);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_f1(
        &coords_fixed[3], &coords_fixed[4], &coords_fixed[5]);

    FTensor::Tensor1<double, 3> t_edge_delta;
    t_edge_delta(i) = t_f1(i) - t_f0(i);

    FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_n(
        v_ptr, v_ptr + 1, v_ptr + 2);
    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_min_coords(
        o_ptr, o_ptr + 1, o_ptr + 2);
    FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_min_dist(min_dist_ptr);

    EntityHandle *colsest_segment_it = nullptr;
    if (o_segments)
      colsest_segment_it = o_segments;

    for (int n = 0; n != nb; ++n) {

      double t;
      if (Tools::minDistancePointFromOnSegment(&t_f0(0), &t_f1(0), &t_n(0),
                                               &t) == Tools::SOLUTION_EXIST) {
        auto t_p = get_point(t_f0, t_edge_delta, t);
        auto dist_n = get_distance(t_p, t_n);
        if (dist_n < t_min_dist || t_min_dist < 0) {
          t_min_dist = dist_n;
          if (o_ptr)
            t_min_coords(i) = t_p(i);
          if (o_segments)
            *colsest_segment_it = e;
        }
      }

      ++t_n;
      ++t_min_dist;
      if (o_ptr)
        ++t_min_coords;
      if (o_segments)
        ++colsest_segment_it;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::outerProductOfEdgeIntegrationPtsForQuad(
    MatrixDouble &gauss_pts, const int rule_ksi, const int rule_eta) {
  MoFEMFunctionBegin;

  auto check_rule_edge = [](int rule) {
    MoFEMFunctionBeginHot;
    if (rule < 0) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong integration rule: %d", rule);
    }
    if (rule > QUAD_1D_TABLE_SIZE) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_1D_TABLE_SIZE);
    }
    if (QUAD_1D_TABLE[rule]->dim != 1) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
    }
    if (QUAD_1D_TABLE[rule]->order < rule) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong order %d != %d", QUAD_1D_TABLE[rule]->order, rule);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR check_rule_edge(rule_ksi);
  CHKERR check_rule_edge(rule_eta);

  const int nb_gauss_pts_ksi = QUAD_1D_TABLE[rule_ksi]->npoints;
  const int nb_gauss_pts_eta = QUAD_1D_TABLE[rule_eta]->npoints;
  gauss_pts.resize(3, nb_gauss_pts_ksi * nb_gauss_pts_eta, false);

  int gg = 0;
  for (size_t i = 0; i != nb_gauss_pts_ksi; ++i) {
    const double wi = QUAD_1D_TABLE[rule_ksi]->weights[i];
    const double ksi = (QUAD_1D_TABLE[rule_ksi]->points[2 * i + 1]);
    for (size_t j = 0; j != nb_gauss_pts_eta; ++j, ++gg) {
      const double wk = wi * QUAD_1D_TABLE[rule_eta]->weights[j];
      const double eta = QUAD_1D_TABLE[rule_eta]->points[2 * j + 1];
      gauss_pts(0, gg) = ksi;
      gauss_pts(1, gg) = eta;
      gauss_pts(2, gg) = wk;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::outerProductOfEdgeIntegrationPtsForHex(
    MatrixDouble &gauss_pts, const int rule_ksi, const int rule_eta,
    const int rule_zeta) {
  MoFEMFunctionBegin;

  auto check_rule_edge = [](int rule) {
    MoFEMFunctionBeginHot;
    if (rule < 0) {
      SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "Wrong integration rule: %d", rule);
    }
    if (rule > QUAD_1D_TABLE_SIZE) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "rule > quadrature order %d < %d", rule, QUAD_1D_TABLE_SIZE);
    }
    if (QUAD_1D_TABLE[rule]->dim != 1) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "wrong dimension");
    }
    if (QUAD_1D_TABLE[rule]->order < rule) {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
               "wrong order %d != %d", QUAD_1D_TABLE[rule]->order, rule);
    }
    MoFEMFunctionReturnHot(0);
  };

  CHKERR check_rule_edge(rule_ksi);
  CHKERR check_rule_edge(rule_eta);
  CHKERR check_rule_edge(rule_zeta);

  const int nb_gauss_pts_ksi = QUAD_1D_TABLE[rule_ksi]->npoints;
  const int nb_gauss_pts_eta = QUAD_1D_TABLE[rule_eta]->npoints;
  const int nb_gauss_pts_zeta = QUAD_1D_TABLE[rule_zeta]->npoints;
  gauss_pts.resize(4, nb_gauss_pts_ksi * nb_gauss_pts_eta * nb_gauss_pts_zeta,
                   false);

  int gg = 0;
  for (size_t i = 0; i != nb_gauss_pts_ksi; ++i) {
    const double wi = QUAD_1D_TABLE[rule_ksi]->weights[i];
    const double ksi = QUAD_1D_TABLE[rule_ksi]->points[2 * i + 1];
    for (size_t j = 0; j != nb_gauss_pts_eta; ++j) {
      const double wj = wi * QUAD_1D_TABLE[rule_eta]->weights[j];
      const double eta = QUAD_1D_TABLE[rule_eta]->points[2 * j + 1];
      for (size_t k = 0; k != nb_gauss_pts_zeta; ++k, ++gg) {
        const double wk = wj * QUAD_1D_TABLE[rule_zeta]->weights[k];
        const double zeta = QUAD_1D_TABLE[rule_eta]->points[2 * k + 1];
        gauss_pts(0, gg) = ksi;
        gauss_pts(1, gg) = eta;
        gauss_pts(2, gg) = zeta;
        gauss_pts(3, gg) = wk;
      }
    }
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM
