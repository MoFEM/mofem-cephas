/** \file UserDataOperators.cpp

\brief Generic user data operators for evaluate fields, and other common
purposes.

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
OpCalculateHoCoords::doWork(int side, EntityType type,
                            DataForcesAndSourcesCore::EntData &data) {
  FTensor::Index<'i', 3> i;
  MoFEMFunctionBegin;
  const auto nb_dofs = data.getFieldData().size() / 3;
  if (nb_dofs) {
    if (type == MBVERTEX)
      getCoordsAtGaussPts().clear();
    auto t_base = data.getFTensor0N();
    auto t_coords = getFTensor1CoordsAtGaussPts();
    const auto nb_integration_pts = data.getN().size1();
    const auto nb_base_functions = data.getN().size2();
    for (auto gg = 0; gg != nb_integration_pts; ++gg) {
      auto t_dof = data.getFTensor1FieldData<3>();
      size_t bb = 0;
      for (; bb != nb_dofs; ++bb) {
        t_coords(i) += t_base * t_dof(i);
        ++t_base;
      }
      for (; bb != nb_base_functions; ++bb)
        ++t_base;
      ++t_coords;
    }
  }
  MoFEMFunctionReturn(0);
};

// MoFEMErrorCode
// OpSetHoInvJacScalarBase::doWork(int side, EntityType type,
//                                 DataForcesAndSourcesCore::EntData &data) {
//   MoFEMFunctionBegin;

//   if (invHoJac.size2() != 9)
//     SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//              "It looks that ho inverse of Jacobian is not calculated %d != 9",
//              invHoJac.size2());

//   auto transform_base = [&](MatrixDouble &diff_n) {
//     MoFEMFunctionBeginHot;

//     unsigned int nb_gauss_pts = diff_n.size1();
//     if (nb_gauss_pts == 0)
//       MoFEMFunctionReturnHot(0);

//     if (invHoJac.size1() == nb_gauss_pts) {

//       unsigned int nb_base_functions = diff_n.size2() / 3;
//       if (nb_base_functions == 0)
//         MoFEMFunctionReturnHot(0);

//       double *t_inv_jac_ptr = &*invHoJac.data().begin();
//       FTensor::Tensor2<double *, 3, 3> t_inv_jac(
//           t_inv_jac_ptr, &t_inv_jac_ptr[1], &t_inv_jac_ptr[2],
//           &t_inv_jac_ptr[3], &t_inv_jac_ptr[4], &t_inv_jac_ptr[5],
//           &t_inv_jac_ptr[6], &t_inv_jac_ptr[7], &t_inv_jac_ptr[8], 9);

//       diffNinvJac.resize(nb_gauss_pts, 3 * nb_base_functions, false);

//       double *t_diff_n_ptr = &*diff_n.data().begin();
//       FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3> t_diff_n(
//           t_diff_n_ptr, &t_diff_n_ptr[1], &t_diff_n_ptr[2]);
//       double *t_inv_n_ptr = &*diffNinvJac.data().begin();
//       FTensor::Tensor1<double *, 3> t_inv_diff_n(t_inv_n_ptr, &t_inv_n_ptr[1],
//                                                  &t_inv_n_ptr[2], 3);

//       switch (type) {
//       case MBVERTEX:
//       case MBEDGE:
//       case MBTRI:
//       case MBTET: {
//         for (unsigned int gg = 0; gg < nb_gauss_pts; ++gg) {
//           for (unsigned int bb = 0; bb != nb_base_functions; ++bb) {
//             t_inv_diff_n(i) = t_diff_n(j) * t_inv_jac(j, i);
//             ++t_diff_n;
//             ++t_inv_diff_n;
//           }
//           ++t_inv_jac;
//         }
//       } break;
//       default:
//         SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED, "not implemented");
//       }

//       diff_n.data().swap(diffNinvJac.data());
//     }
//     MoFEMFunctionReturnHot(0);
//   };

//   for (int b = AINSWORTH_LEGENDRE_BASE; b != LASTBASE; b++) {
//     FieldApproximationBase base = static_cast<FieldApproximationBase>(b);
//     CHKERR transform_base(data.getDiffN(base));
//   }

//   switch (type) {
//   case MBVERTEX:
//     for (auto &m : data.getBBDiffNMap())
//       CHKERR transform_base(*(m.second));
//     break;
//   default:
//     for (auto &ptr : data.getBBDiffNByOrderArray())
//       if (ptr)
//         CHKERR transform_base(*ptr);
//   }

//   MoFEMFunctionReturn(0);
// }

MoFEMErrorCode OpMakeHighOrderGeometryWeightsOnFace::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_int_pts = getGaussPts().size2();
  if (getNormalsAtGaussPts().size1()) {
    if (getNormalsAtGaussPts().size1() == nb_int_pts) {
      const double a = getMeasure();
      auto t_w = getFTensor0IntegrationWeight();
      auto t_normal = getFTensor1NormalsAtGaussPts();
      FTensor::Index<'i', 3> i;
      for (size_t gg = 0; gg != nb_int_pts; ++gg) {
        t_w *= sqrt(t_normal(i) * t_normal(i)) / a;
        ++t_w;
        ++t_normal;
      }
    } else {
      SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
               "Number of rows in getNormalsAtGaussPts should be equal to "
               "number of integration points, but is not, i.e. %d != %d",
               getNormalsAtGaussPts().size1(), nb_int_pts);
    }
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode OpMakeHighOrderGeometryWeightsOnVolume::doWork(
    int side, EntityType type, DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;
  const size_t nb_int_pts = getGaussPts().size2();

  if (getHoGaussPtsDetJac().size() == nb_int_pts) {
    const double a = getMeasure();
    auto t_w = getFTensor0IntegrationWeight();
    auto t_w_ho = getFTensor0FromVec(getHoGaussPtsDetJac());

    for (size_t gg = 0; gg != nb_int_pts; ++gg) {
      t_w *= t_w_ho;
      ++t_w;
      ++t_w_ho;
    }
  } else {
    SETERRQ2(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
             "Number of rows in getHoGaussPtsDetJac should be equal to "
             "number of integration points, but is not, i.e. %d != %d",
             getHoGaussPtsDetJac().size(), nb_int_pts);
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM