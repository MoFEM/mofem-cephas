/** \file BaseDirevativesDataOperators.cpp


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

#include <cholesky.hpp>

namespace MoFEM {

MoFEMErrorCode
OpBaseDerivativesMass<1>::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

//   if (sPace != L2) {
//     SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//             "Space should be set to L2");
//   }

//   if (!data.getNSharedPtr(base)) {
//     auto fe_ptr = getPtrFE();
//     auto fe_type = getFEType();
//     // Set data structure to store bas
//     dataL2->dataOnEntities[fe_type].clear();
//     dataL2->dataOnEntities[MBVERTEX].push_back(
//         new DataForcesAndSourcesCore::EntData());
//     dataL2->dataOnEntities[fe_type].push_back(
//         new DataForcesAndSourcesCore::EntData());

//     auto &vertex_data = dataL2->dataOnEntities[MBVERTEX][0];
//     vertex_data.getNSharedPtr(NOBASE) =
//         fe_ptr->dataH1.dataOnEntities[MBVERTEX][0]->getNSharedPtr(NOBASE);
//     vertex_data.getDiffNSharedPtr(NOBASE) =
//         fe_ptr->dataH1.dataOnEntities[MBVERTEX][0]->getDiffNSharedPtr(NOBASE);

//     auto &ent_data = dataL2->dataOnEntities[fe_type][0];
//     ent_data.sEnse = 1;
//     ent_data.bAse = base;
//     ent_data.oRder = fe_ptr->getMaxDataOrder();

//     auto base_functions_generator = fe_ptr->getElementPolynomialBase();
//     CHKERR base_functions_generator->getValue(
//         getGaussPts(), boost::make_shared<EntPolynomialBaseCtx>(
//                            *dataL2, static_cast<FieldSpace>(L2),
//                            static_cast<FieldApproximationBase>(base), NOBASE));
//   }

//   auto &base = dataL2.getN(base);
//   const auto nb = base.size2();

//   if (&nb) {

//     const auto nb_integration_pts = base.size1();

// #ifndef NDEBUG
//     if (baseMassPtr)
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Mass matrix is null pointer");
// #endif

//     auto &nN = *baseMassPtr;
//     nN.resize(nb, nb, false);
//     nN.clear();

//     auto t_w = getFTensor0IntegrationWeight();
//     // get base function gradient on rows
//     auto t_row_base = dataL2.getFTensor0N();
//     // loop over integration points
//     for (int gg = 0; gg != nb_integration_pts; ++gg) {
//       // take into account Jacobian
//       const double alpha = t_w * vol;

//       for (int rr = 0; rr != nb; ++rr) {

//         // loop over rows base functions
//         auto a_mat_ptr = &*nN(rr, 0);
//         // get column base functions gradient at gauss point gg
//         auto t_col_base = dataL2.getFTensor0N(gg, 0);
//         // loop over columns
//         for (int cc = 0; cc <= rr; ++cc) {
//           // calculate element of local matrix
//           *a_mat_ptr += alpha * (t_row_base * t_col_base);
//           ++t_col_base;
//           ++a_mat_ptr;
//         }

//         ++t_row_base;
//       }
//       ++t_w; // move to another integration weight
//     }
//   }

//   cholesky_decompose(nN);

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpBaseDerivativesNext<1>::doWork(int side, EntityType type,
                                 DataForcesAndSourcesCore::EntData &data) {
  MoFEMFunctionBegin;

//   auto &base = data.getN(base);
//   const auto nb = base.size2();

//   if (&nb) {

//     auto fe_type = getFEType();
//     const auto nb_integration_pts = base.size1();

//     int base_direvative = 0;
//     for (auto base_direvative = 0; base_direvative != LastDerivative;
//          ++base_direvative) {
//       if (!data.getNSharedPtr(base, static_cast<BaseDirevatives>(d))) {
//         ++base_direvative;
//         break;
//       }
//     }

//     if (base_direvative == LastDerivative)
//       MoFEMFunctionReturnHot(0);

// #ifndef NDEBUG
//     if (baseMassPtr)
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Mass matrix is null pointer");
//     if (base_direvative < 2)
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Direvative should be larger than 2");
//     if (base_direvative < data.baseFunctionsAndBaseDirevatives.size())
//       SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//                "Direvative should be larger than %d",
//                data.baseFunctionsAndBaseDirevatives.size());
//     if (data.getNSharedPtr(base, base_direvative - 1))
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Derivative does not exist");
// #endif

//     const space_dim = data.getDiffN().size2() / nb;
//     auto &diff_base = *(data.getNSharedPtr(base, base_direvative - 1));

// #ifndef NDEBUG
//     if (diff_base.size1() != nb_integration_pts)
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Number of integration points is not consistent");
//     if (nb * space_dim == nb)
//       SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
//               "Number of base funtions derivatives is not consistent");
// #endif

//     auto &nN = *baseMassPtr;

//     int nb_direvatives = pow(space_dim, base_direvative - 1);
//     nF.resize(nb, nb_direvatives, false);
//     nF.clear();

//     auto diff_base_ptr = *diff_base.data().begin();

//     auto t_w = getFTensor0IntegrationWeight();
//     // get base function gradient on rows
//     auto &ent_data = dataL2->dataOnEntities[fe_type][0];
//     const int nb_base_functions = ent_data->getN().size2();
//     auto t_row_base = ent_data.getFTensor0N();
//     // loop over integration points
//     for (int gg = 0; gg != nb_integration_pts; ++gg) {
//       // take into account Jacobian
//       const double alpha = t_w * vol;

//       for (int rr = 0; rr != nb_base_functions; ++rr) {

//         auto f_vec_ptr = &*nF.data().begin();
//         for (int d = 0; d != nb * nb_direvatives; ++d) {
//           *f_vec_ptr += alpha * (t_row_base * (*diff_base_ptr));
//           ++diff_base_ptr;
//           ++f_vec_ptr;
//         }

//         ++t_row_base;
//       }
//       ++t_w; // move to another integration weight
//     }
//   }

//   for (auto rr = 0; rr != nb; ++rr) {
//     for (int d = 0; d != nb_direvatives; ++d) {
//       ublas::matrix_column<MatrixDouble> mc(mF, nb * nb_direvatives + d);
//       cholesky_solve(nN, nF, ublas::lower());
//     }
//   }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM