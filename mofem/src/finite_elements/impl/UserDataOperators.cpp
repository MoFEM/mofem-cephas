/** \file ForcesAndSourcesCore.cpp

\brief Implementation of Elements on Entities for Forces and Sources

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

namespace MoFEM {

MoFEMErrorCode
OpCalculateJacForFace::doWork(int side, EntityType type,
                              DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) 
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");
  
  if (type == MBVERTEX) {
    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();
    double j00 = 0, j01 = 0, j10 = 0, j11 = 0;

    // this is triangle, derivative of nodal shape functions is constant.
    // So only need to do one node.
    for (auto n : {0, 1, 2}) {
      j00 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j01 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 1];
      j10 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j11 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 1];
    }

    jac.resize(2, 2, false);
    jac(0, 0) = j00;
    jac(0, 1) = j01;
    jac(1, 0) = j10;
    jac(1, 1) = j11;

  }

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets = false;
  doPrisms = false;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpCalculateInvJacForFace::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) 
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");

  if (type == MBVERTEX) {
    VectorDouble &coords = getCoords();
    double *coords_ptr = &*coords.data().begin();
    double j00 = 0, j01 = 0, j10 = 0, j11 = 0;

    // this is triangle, derivative of nodal shape functions is constant.
    // So only need to do one node.
    for (auto n : {0, 1, 2}) {
      j00 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j01 += coords_ptr[3 * n + 0] * Tools::diffShapeFunMBTRI[2 * n + 1];
      j10 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 0];
      j11 += coords_ptr[3 * n + 1] * Tools::diffShapeFunMBTRI[2 * n + 1];
    }
    const double det = j00 * j11 - j01 * j10;
    
    invJac.resize(2, 2, false);
    invJac(0, 0) = j11 / det;
    invJac(0, 1) = -j01 / det;
    invJac(1, 0) = -j10 / det;
    invJac(1, 1) = j00 / det;
  }

  doVertices = true;
  doEdges = false;
  doQuads = false;
  doTris = false;
  doTets = false;
  doPrisms = false;

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacH1ForFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI &&
      getNumeredEntFiniteElementPtr()->getEntType() != MBQUAD) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");
  }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    unsigned int nb_functions = data.getN(base).size2();
    if (nb_functions) {
      unsigned int nb_gauss_pts = data.getN(base).size1();
      diffNinvJac.resize(nb_gauss_pts, 2 * nb_functions, false);

      if (type != MBVERTEX) {
        if (nb_functions != data.getDiffN(base).size2() / 2) {
          SETERRQ2(
              PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "data inconsistency nb_functions != data.diffN.size2()/2 ( %u != "
              "%u/2 )",
              nb_functions, data.getDiffN(base).size2());
        }
      }

      FTensor::Tensor2<double, 2, 2> t_inv_jac;
      t_inv_jac(0, 0) = invJac(0, 0);
      t_inv_jac(0, 1) = invJac(0, 1);
      t_inv_jac(1, 0) = invJac(1, 0);
      t_inv_jac(1, 1) = invJac(1, 1);

      switch (type) {
      case MBVERTEX:
      case MBEDGE:
      case MBTRI:
      case MBQUAD: {
        FTensor::Index<'i', 2> i;
        FTensor::Index<'j', 2> j;
        FTensor::Index<'k', 2> k;
        FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n(
            &diffNinvJac(0, 0), &diffNinvJac(0, 1));
        FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2> t_diff_n_ref(
            &data.getDiffN(base)(0, 0), &data.getDiffN(base)(0, 1));
        for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (unsigned int dd = 0; dd != nb_functions; ++dd) {
            t_diff_n(i) = t_inv_jac(k, i) * t_diff_n_ref(k);
            ++t_diff_n;
            ++t_diff_n_ref;
          }
        }
        data.getDiffN(base).data().swap(diffNinvJac.data());
      } break;
      default:
        SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpSetInvJacHcurlFace::doWork(int side, EntityType type,
                             DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI)
    MoFEMFunctionReturnHot(0);

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) 
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");

  FTensor::Tensor2<double *, 2, 2> t_inv_jac = FTensor::Tensor2<double *, 2, 2>(
      &invJac(0, 0), &invJac(0, 1), &invJac(1, 0), &invJac(1, 1));

  FTensor::Index<'i', 3> i;
  FTensor::Index<'j', 2> j;
  FTensor::Index<'k', 2> k;

  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getDiffN(base).size2() / 6;
    if (nb_base_functions) {
      const unsigned int nb_gauss_pts = data.getDiffN(base).size1();

      diffHcurlInvJac.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);

      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      double *t_inv_diff_n_ptr = &*diffHcurlInvJac.data().begin();
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2> t_inv_diff_n(
          t_inv_diff_n_ptr, &t_inv_diff_n_ptr[HVEC0_1],

          &t_inv_diff_n_ptr[HVEC1_0], &t_inv_diff_n_ptr[HVEC1_1],

          &t_inv_diff_n_ptr[HVEC2_0], &t_inv_diff_n_ptr[HVEC2_1]);

      for (unsigned int gg = 0; gg != nb_gauss_pts; gg++) {
        for (unsigned int bb = 0; bb != nb_base_functions; bb++) {
          t_inv_diff_n(i, j) = t_diff_n(i, k) * t_inv_jac(k, j);
          ++t_diff_n;
          ++t_inv_diff_n;
        }
      }

      data.getDiffN(base).data().swap(diffHcurlInvJac.data());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpMakeHdivFromHcurl::doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI)
    MoFEMFunctionReturnHot(0);

  if (getNumeredEntFiniteElementPtr()->getEntType() != MBTRI) 
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "This operator can be used only with element which is triangle");


  for (int b = AINSWORTH_LEGENDRE_BASE; b != USER_BASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getN(base).size2() / 3;
    if (nb_base_functions) {
      const unsigned int nb_gauss_pts = data.getDiffN(base).size1();

      auto t_n = data.getFTensor1N<3>(base);
      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {

          const double a = t_n(0);
          t_n(0) = -t_n(1);
          t_n(1) = a;

          for (auto n : {0, 1}) {
            const double b = t_diff_n(0, n);
            t_diff_n(0, n) = -t_diff_n(1, n);
            t_diff_n(1, n) = b;
          }

          ++t_n;
          ++t_diff_n;
        }
      }
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContravariantPiolaTransformFace ::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {

  MoFEMFunctionBegin;

  if (type != MBEDGE && type != MBTRI)
    MoFEMFunctionReturnHot(0);

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);

    const unsigned int nb_base_functions = data.getN(base).size2() / 3;
    if (!nb_base_functions)
      continue;

    const unsigned int nb_gauss_pts = data.getN(base).size1();

    double det;
    CHKERR determinantTensor2by2(tJac, det);

    piolaN.resize(nb_gauss_pts, data.getN(base).size2(), false);
    if (data.getN(base).size2() > 0) {
      auto t_n = data.getFTensor1N<3>(base);
      double *t_transformed_n_ptr = &*piolaN.data().begin();
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_transformed_n(
          t_transformed_n_ptr, // HVEC0
          &t_transformed_n_ptr[HVEC1], &t_transformed_n_ptr[HVEC2]);
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
          t_transformed_n(i) = tJac(i, k) * t_n(k) / det;
          ++t_n;
          ++t_transformed_n;
        }
      }
      data.getN(base).data().swap(piolaN.data());
    }

    piolaDiffN.resize(nb_gauss_pts, data.getDiffN(base).size2(), false);
    if (data.getDiffN(base).size2() > 0) {
      auto t_diff_n = data.getFTensor2DiffN<3, 2>(base);
      double *t_transformed_diff_n_ptr = &*piolaDiffN.data().begin();
      FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
          t_transformed_diff_n(
            t_transformed_diff_n_ptr, &t_transformed_diff_n_ptr[HVEC0_1],

            &t_transformed_diff_n_ptr[HVEC1_0], &t_transformed_diff_n_ptr[HVEC1_1],

            &t_transformed_diff_n_ptr[HVEC2_0], &t_transformed_diff_n_ptr[HVEC2_1]);
      for (unsigned int gg = 0; gg != nb_gauss_pts; ++gg) {
        for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
          t_transformed_diff_n(i, k) = tJac(i, j) * t_diff_n(j, k) / det;
          ++t_diff_n;
          ++t_transformed_diff_n;
        }
      }
      data.getDiffN(base).data().swap(piolaDiffN.data());
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpSetContrariantPiolaTransformOnEdge::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBeginHot;

  if (type != MBEDGE)
    MoFEMFunctionReturnHot(0);

  const auto &edge_direction = getDirection();

  FTensor::Index<'i', 3> i;
  FTensor::Tensor1<double, 3> t_m(-edge_direction[1], edge_direction[0], 0.);
  const double l0 = t_m(i) * t_m(i);

  std::vector<double> l1;
  {
    int nb_gauss_pts = getTangetAtGaussPts().size1();
    if (nb_gauss_pts) {
      l1.resize(nb_gauss_pts);
      const auto &edge_direction = getTangetAtGaussPts();
      FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_m_at_pts(
          &edge_direction(0, 0), &edge_direction(0, 1), &edge_direction(0, 2));
      for (int gg = 0; gg < nb_gauss_pts; ++gg) {
        l1[gg] = t_m_at_pts(i) * t_m_at_pts(i);
        ++t_m_at_pts;
      }
    }
  }

  for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {

    FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
    size_t nb_gauss_pts = data.getN(base).size1();
    size_t nb_dofs = data.getN(base).size2() / 3;
    if (nb_gauss_pts > 0 && nb_dofs > 0) {
      FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_h_div(
          &data.getN(base)(0, HVEC0), &data.getN(base)(0, HVEC1),
          &data.getN(base)(0, HVEC2));

      size_t cc = 0;
      const auto &edge_direction_at_gauss_pts = getTangetAtGaussPts();
      if (edge_direction_at_gauss_pts.size1() == nb_gauss_pts) {

        FTensor::Tensor1<FTensor::PackPtr<const double *, 3>, 3> t_m_at_pts(
            &edge_direction_at_gauss_pts(0, 1),
            &edge_direction_at_gauss_pts(0, 0),
            &edge_direction_at_gauss_pts(0, 2));

        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          const double l0 = l1[gg];
          for (int ll = 0; ll != nb_dofs; ++ll) {
            const double val = t_h_div(0);
            const double a = val / l0;
            t_h_div(i) = t_m_at_pts(i) * a;
            t_h_div(0) *= -1;
            ++t_h_div;
            ++cc;
          }
          ++t_m_at_pts;
        }

      } else {

        for (int gg = 0; gg != nb_gauss_pts; ++gg) {
          for (int ll = 0; ll != nb_dofs; ll++) {
            const double val = t_h_div(0);
            const double a = val / l0;
            t_h_div(i) = t_m(i) * a;
            ++t_h_div;
            ++cc;
          }
        }
        
      }

      if (cc != nb_gauss_pts * nb_dofs) 
        SETERRQ(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE, "Data inconsistency");
      
    }
  }

  MoFEMFunctionReturnHot(0);
}

}
