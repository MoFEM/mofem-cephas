/** \file Tools.cpp
 * \brief Auxilairy tools
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
    if(th) {
      CHKERR moab.tag_get_data(th, conn, num_nodes, coords);
    } else {
      CHKERR moab.get_coords(conn, num_nodes, coords);
    }
    double q = Tools::volumeLengthQuality(coords);
    min_quality = f(q, min_quality);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
Tools::getTetsWithQuality(Range &out_tets, const Range &tets,
                          Tag th, boost::function<bool(double)> f) {
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

MoFEMErrorCode
Tools::writeTetsWithQuality(const char *file_name, const char *file_type,
                            const char *options, const Range &tets,
                            Tag th, boost::function<bool(double)> f) {
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
                                          const double tol,bool &result) {
  double loc_coord[] = {0, 0, 0};
  double N[4],diffN[12];
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
  boost::shared_ptr< NumeredDofEntity_multiIndex > prb_dofs;
  switch(row_or_col) {
    case ROW:
      prb_loc_size = prb_ptr->getNbLocalDofsRow();
      prb_dofs = prb_ptr->getNumeredDofsRows();
      break;
    case COL:
      prb_loc_size = prb_ptr->getNbLocalDofsCol();
      prb_dofs = prb_ptr->getNumeredDofsCols();
      break;
    break;
   default:
     SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "Wrong argument, row_or_col should be row or column");
  }
  if(loc_size != prb_loc_size) {
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

}
