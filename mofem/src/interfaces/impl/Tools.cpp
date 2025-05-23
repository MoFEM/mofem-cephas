/** \file Tools.cpp
 * \brief Auxiliary tools
 */

#include <phg-quadrule/quad.h>

namespace MoFEM {

MoFEMErrorCode Tools::query_interface(boost::typeindex::type_index type_index,
                                      UnknownInterface **iface) const {
  *iface = const_cast<Tools *>(this);
  return 0;
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
  return determinantTensor3by3(jac) / 6.;
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

constexpr double Tools::shapeFunMBEDGE0At00;
constexpr double Tools::shapeFunMBEDGE1At00;
constexpr std::array<double, 2> Tools::shapeFunMBEDGEAt00;

constexpr std::array<double, 2> Tools::diffShapeFunMBEDGE;
constexpr std::array<double, 6> Tools::diffShapeFunMBTRI;
constexpr std::array<double, 12> Tools::diffShapeFunMBTET;
constexpr double Tools::shapeFunMBTRI0At00;
constexpr double Tools::shapeFunMBTRI1At00;
constexpr double Tools::shapeFunMBTRI2At00;
constexpr std::array<double, 3> Tools::shapeFunMBTRIAt00;

constexpr std::array<double, 4> Tools::shapeFunMBTETAt000;
constexpr std::array<double, 8> Tools::diffShapeFunMBQUADAtCenter;
constexpr std::array<double, 24> Tools::diffShapeFunMBHEXAtCenter;

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
  // jj - local derivatives
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

template <typename T1, typename T2>
MoFEMErrorCode getLocalCoordinatesOnReferenceThreeNodeTriImpl(
    const T1 *elem_coords, const T2 *global_coords, const int nb_nodes,
    typename FTensor::promote<T1, T2>::V *local_coords) {

  FTensor::Index<'i', 3> i3;
  FTensor::Index<'j', 3> j3;
  FTensor::Index<'K', 2> K2;
  FTensor::Index<'L', 2> L2;
  MoFEMFunctionBeginHot;

  FTensor::Tensor1<const double, 3> t_n = {Tools::shapeFunMBTRIAt00[0],
                                           Tools::shapeFunMBTRIAt00[1],
                                           Tools::shapeFunMBTRIAt00[2]};
  auto t_diff = getFTensor2FromPtr<3, 2>(
      const_cast<double *>(Tools::diffShapeFunMBTRI.data()));
  auto t_elem_coords = getFTensor2FromPtr<3, 3>(const_cast<T1 *>(elem_coords));

  // Build matrix and get coordinates of zero point
  FTensor::Tensor2<T1, 2, 3> t_a;
  t_a(K2, i3) = t_diff(j3, K2) * t_elem_coords(j3, i3);
  FTensor::Tensor2_symmetric<T1, 2> t_b, t_inv_b;
  t_b(K2, L2) = t_a(K2, j3) ^ t_a(L2, j3);
  // Solve problem
  const auto inv_det = 1. / (t_b(0, 0) * t_b(1, 1) - t_b(0, 1) * t_b(1, 0));
  t_inv_b(0, 0) = t_b(1, 1) * inv_det;
  t_inv_b(0, 1) = -t_b(0, 1) * inv_det;
  t_inv_b(1, 1) = t_b(0, 0) * inv_det;

  // Coords at corner
  FTensor::Tensor1<T1, 3> t_coords_at_0;
  t_coords_at_0(i3) = t_n(j3) * t_elem_coords(j3, i3);

  auto t_global_coords = getFTensor1FromPtr<3>(const_cast<T2 *>(global_coords));
  auto t_local_coords = getFTensor1FromPtr<2>(local_coords);

  // Calculate right hand side
  for (int ii = 0; ii != nb_nodes; ++ii) {
    t_local_coords(L2) =
        t_inv_b(L2, K2) *
        (t_a(K2, j3) * (t_global_coords(j3) - t_coords_at_0(j3)));
    ++t_local_coords;
    ++t_global_coords;
  }

  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode Tools::getLocalCoordinatesOnReferenceThreeNodeTri(
    const double *elem_coords, const double *glob_coords, const int nb_nodes,
    double *local_coords) {
  return getLocalCoordinatesOnReferenceThreeNodeTriImpl<double, double>(
      elem_coords, glob_coords, nb_nodes, local_coords);
}

MoFEMErrorCode Tools::getLocalCoordinatesOnReferenceThreeNodeTri(
    const double *elem_coords, const std::complex<double> *glob_coords,
    const int nb_nodes, std::complex<double> *local_coords) {
  return getLocalCoordinatesOnReferenceThreeNodeTriImpl<double,
                                                        std::complex<double>>(
      elem_coords, glob_coords, nb_nodes, local_coords);
}

MoFEMErrorCode Tools::getLocalCoordinatesOnReferenceEdgeNodeEdge(
    const double *elem_coords, const double *global_coords, const int nb_nodes,
    double *local_coords) {

  FTensor::Index<'i', 3> i3;
  MoFEMFunctionBeginHot;

  FTensor::Tensor1<FTensor::PackPtr<const double *, 1>, 3> t_elem_coords = {
      &elem_coords[0], &elem_coords[3], &elem_coords[6]};

  FTensor::Tensor1<double, 3> t_coords_at_0;
  // Build matrix and get coordinates of zero point
  // ii - global coordinates
  FTensor::Tensor1<double, 3> t_a;
  for (auto ii : {0, 1, 2}) {
    t_a(ii) = diffShapeFunMBEDGE[0] * t_elem_coords(0) +
              diffShapeFunMBEDGE[1] * t_elem_coords(1);
    t_coords_at_0(ii) = shapeFunMBEDGEAt00[0] * t_elem_coords(0) +
                        shapeFunMBEDGEAt00[1] * t_elem_coords(1);
    ++t_elem_coords;
  }

  FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_global_coords = {
      &global_coords[0], &global_coords[1], &global_coords[2]};
  FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_local_coords = {
      &local_coords[0]};

  const double b = t_a(i3) * t_a(i3);
  const double inv_b = 1 / b;

  // Calculate right hand side
  for (int ii = 0; ii != nb_nodes; ++ii) {
    t_local_coords =
        inv_b * (t_a(i3) * (t_global_coords(i3) - t_coords_at_0(i3)));
    ++t_local_coords;
    ++t_global_coords;
  }

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

MoFEMErrorCode Tools::getTriNormal(const double *coords, double *normal,
                                   double *d_normal) {
  MoFEMFunctionBegin;
  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 3> j;
  FTensor::Index<'k', 3> k;
  FTensor::Index<'l', 3> l;
  FTensor::Index<'n', 3> n;
  FTensor::Index<'m', 3> m;
  FTensor::Index<'J', 2> J;
  FTensor::Number<0> N0;
  FTensor::Number<1> N1;
  auto diff_ptr = Tools::diffShapeFunMBTRI.data();
  auto t_diff_tensor = getFTensor2FromPtr<3, 2>(const_cast<double *>(diff_ptr));
  auto t_coords = getFTensor2FromPtr<3, 3>(const_cast<double *>(coords));
  FTensor::Tensor2<double, 3, 2> t_tangent;
  t_tangent(i, J) = t_coords(n, i) * t_diff_tensor(n, J);
  auto t_normal = getFTensor1FromPtr<3>(normal);
  t_normal(j) =
      (FTensor::levi_civita(i, j, k) * t_tangent(k, N0)) * t_tangent(i, N1);
  if (d_normal) {
    constexpr auto t_kd = FTensor::Kronecker_Delta<int>();
    FTensor::Tensor4<int, 3, 3, 3, 3> t_d_coords;
    t_d_coords(i, j, k, n) = t_kd(i, k) * t_kd(j, n);
    FTensor::Tensor4<double, 3, 3, 3, 2> t_d_tangent;
    t_d_tangent(i, k, l, J) = t_d_coords(n, i, k, l) * t_diff_tensor(n, J);
    auto t_d_normal = getFTensor3FromPtr<3, 3, 3>(d_normal);
    t_d_normal(j, m, n) = (FTensor::levi_civita(i, j, k) * t_tangent(i, N1)) *
                              t_d_tangent(k, m, n, N0)

                          +

                          (FTensor::levi_civita(i, j, k) * t_tangent(k, N0)) *
                              t_d_tangent(i, m, n, N1);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode Tools::getTriNormal(const EntityHandle tri,
                                   double *normal) const {
  MoFEM::Interface &m_field = cOre;
  moab::Interface &moab(m_field.get_moab());
  MoFEMFunctionBegin;
  if (type_from_handle(tri) != MBTRI) {
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
  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_normal;
  CHK_THROW_MESSAGE(getTriNormal(tri, &t_normal(0)), "calculate area");
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
    if (type_from_handle(edge) != MBEDGE) {
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
  if (gg != gauss_pts.size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong size of matrix");

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
        const double zeta = QUAD_1D_TABLE[rule_zeta]->points[2 * k + 1];
        gauss_pts(0, gg) = ksi;
        gauss_pts(1, gg) = eta;
        gauss_pts(2, gg) = zeta;
        gauss_pts(3, gg) = wk;
      }
    }
  }

  if (gg != gauss_pts.size2())
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY, "Wrong size of matrix");

  MoFEMFunctionReturn(0);
}

constexpr std::array<int, 12> Tools::uniformTriangleRefineTriangles;

std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
Tools::refineTriangle(int nb_levels) {

  std::vector<int> triangles{0, 1, 2, 3, 4, 5};
  std::vector<double> nodes{

      0.,  0.,  // 0
      1.,  0.,  // 1
      0.,  1.,  // 2
      0.5, 0.,  // 3
      0.5, 0.5, // 4
      0.,  0.5  // 5

  };
  std::map<std::pair<int, int>, int> edges{
      {{0, 1}, 3}, {{1, 2}, 4}, {{0, 2}, 5}};

  auto add_edge = [&](auto a, auto b) {
    if (a > b) {
      std::swap(a, b);
    }
    auto it = edges.find(std::make_pair(a, b));
    if (it == edges.end()) {
      int e = 3 + edges.size();
      edges[std::make_pair(a, b)] = e;
      for (auto n : {0, 1}) {
        nodes.push_back((nodes[2 * a + n] + nodes[2 * b + n]) / 2);
      }
#ifndef NDEBUG
      if (e != nodes.size() / 2 - 1)
        CHK_THROW_MESSAGE(MOFEM_DATA_INCONSISTENCY, "wrong node/edge index");
#endif
      return e;
    }
    return it->second;
  };

  auto add_triangle = [&](auto t) {
    for (auto tt : {0, 1, 2, 3}) {
      auto last = triangles.size() / 6;
      for (auto n : {0, 1, 2}) {
        // add triangle nodes
        triangles.push_back(
            triangles[6 * t + uniformTriangleRefineTriangles[3 * tt + n]]);
      }
      // add triangle edges
      auto cycle = std::array<int, 4>{0, 1, 2, 0};
      for (auto e : {0, 1, 2}) {
        triangles.push_back(add_edge(triangles[6 * last + cycle[e]],
                                     triangles[6 * last + cycle[e + 1]]));
      }
    }
  };

  std::vector<int> level_index{0, 1};
  auto l = 0;
  for (; l != nb_levels; ++l) {
    auto first_tet = level_index[l];
    auto nb_last_level_test = level_index.back() - level_index[l];
    for (auto t = first_tet; t != (first_tet + nb_last_level_test); ++t) {
      add_triangle(t);
    }
    level_index.push_back(triangles.size() / 6);
  }

  return std::make_tuple(nodes, triangles, level_index);
}

MatrixDouble Tools::refineTriangleIntegrationPts(
    MatrixDouble pts,
    std::tuple<std::vector<double>, std::vector<int>, std::vector<int>>
        refined) {

  auto [nodes, triangles, level_index] = refined;
  
  auto get_coords = [&](auto t) {
    auto [nodes, triangles, level_index] = refined;
    std::array<double, 9> ele_coords;
    for (auto n : {0, 1, 2}) {
      for (auto i : {0, 1}) {
        ele_coords[3 * n + i] = nodes[2 * triangles[6 * t + n] + i];
      }
      ele_coords[3 * n + 2] = 0;
    }
    return ele_coords;
  };

  auto get_normal = [](auto &ele_coords) {
    FTensor::Tensor1<double, 3> t_normal;
    Tools::getTriNormal(ele_coords.data(), &t_normal(0));
    return t_normal;
  };

  std::vector<double> ele_shape(3 * pts.size2());
  shapeFunMBTRI<1>(&*ele_shape.begin(), &pts(0, 0), &pts(1, 0), pts.size2());

  int nb_elems = level_index.back() - level_index[level_index.size() - 2];
  MatrixDouble new_pts(3, pts.size2() * nb_elems);
  FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2> t_gauss_pt{&new_pts(0, 0),
                                                                &new_pts(1, 0)};
  FTensor::Tensor0<FTensor::PackPtr<double *, 1>> t_gauss_weight{
      &new_pts(2, 0)};

  for (auto t = level_index[level_index.size() - 2]; t != level_index.back();
       ++t) {

    auto ele_coords = get_coords(t);
    auto t_normal = get_normal(ele_coords);
    auto area = t_normal.l2();

    FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_ele_shape{
        &ele_shape[0], &ele_shape[1], &ele_shape[2]};
    FTensor::Tensor2<double, 3, 2> t_ele_coords{ele_coords[0], ele_coords[1],
                                                ele_coords[3], ele_coords[4],
                                                ele_coords[6], ele_coords[7]};

    FTensor::Index<'i', 3> i;
    FTensor::Index<'J', 2> J;

    for (auto gg = 0; gg != pts.size2(); ++gg) {
      t_gauss_pt(J) = t_ele_shape(i) * t_ele_coords(i, J);
      t_gauss_weight = area * pts(2, gg);
      ++t_gauss_pt;
      ++t_gauss_weight;
      ++t_ele_shape;
    }

  }

  return new_pts;
}

MatrixDouble
Tools::refineTriangleIntegrationPts(int rule, RefineTrianglesReturn refined) {

  MatrixDouble gauss_pts;

  const size_t nb_gauss_pts = QUAD_2D_TABLE[rule]->npoints;
  gauss_pts.resize(3, nb_gauss_pts, false);
  cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[1], 3,
              &gauss_pts(0, 0), 1);
  cblas_dcopy(nb_gauss_pts, &QUAD_2D_TABLE[rule]->points[2], 3,
              &gauss_pts(1, 0), 1);
  cblas_dcopy(nb_gauss_pts, QUAD_2D_TABLE[rule]->weights, 1, &gauss_pts(2, 0),
              1);

  return Tools::refineTriangleIntegrationPts(gauss_pts, refined);
}

} // namespace MoFEM
