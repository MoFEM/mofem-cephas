/** \file BaseDerivativesDataOperators.cpp


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

OpBaseDerivativesBase::OpBaseDerivativesBase(
    boost::shared_ptr<MatrixDouble> base_mass_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s, int verb, Sev sev)
    : ForcesAndSourcesCore::UserDataOperator(s, OPSPACE), base(b),
      verbosity(verb), severityLevel(sev), baseMassPtr(base_mass_ptr),
      dataL2(data_l2) {}

MoFEMErrorCode
OpBaseDerivativesMass<1>::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  if (sPace != L2) {
    SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
            "Space should be set to L2");
  }

  auto fe_type = getFEType();

  auto calculate_base = [&]() {
    MoFEMFunctionBeginHot;
    auto fe_ptr = getPtrFE();
    // Set data structure to store base
    dataL2->dataOnEntities[fe_type].clear();
    dataL2->dataOnEntities[MBVERTEX].push_back(
        new EntitiesFieldData::EntData());
    dataL2->dataOnEntities[fe_type].push_back(new EntitiesFieldData::EntData());

    auto &vertex_data = dataL2->dataOnEntities[MBVERTEX][0];
    vertex_data.getNSharedPtr(NOBASE) =
        fe_ptr->getEntData(H1, MBVERTEX, 0).getNSharedPtr(NOBASE);
    vertex_data.getDiffNSharedPtr(NOBASE) =
        fe_ptr->getEntData(H1, MBVERTEX, 0).getDiffNSharedPtr(NOBASE);

    auto &ent_data = dataL2->dataOnEntities[fe_type][0];
    ent_data.getSense() = 1;
    ent_data.getBase() = base;
    ent_data.getOrder() = std::max(0, fe_ptr->getMaxDataOrder() - 1);

    CHKERR fe_ptr->getElementPolynomialBase()->getValue(
        getGaussPts(), boost::make_shared<EntPolynomialBaseCtx>(
                           *dataL2, static_cast<FieldSpace>(L2),
                           static_cast<FieldApproximationBase>(base), NOBASE));
    MoFEMFunctionReturnHot(0);
  };

  CHKERR calculate_base();

  auto &ent_data = dataL2->dataOnEntities[fe_type][0];
  auto &base_funcions = ent_data.getN(base);
  const auto nb = base_funcions.size2();

  if (nb) {

    const auto nb_integration_pts = getGaussPts().size2();

#ifndef NDEBUG
    if (!baseMassPtr)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Mass matrix is null pointer");
#endif

    auto &nN = *baseMassPtr;
    nN.resize(nb, nb, false);
    nN.clear();

    auto t_w = getFTensor0IntegrationWeight();
    // get base function gradient on rows
    auto t_row_base = ent_data.getFTensor0N();
    // loop over integration points
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      // take into account Jacobian
      const double alpha = t_w;

      for (int rr = 0; rr != nb; ++rr) {

        // loop over rows base functions
        auto a_mat_ptr = &nN(rr, 0);
        // get column base functions gradient at gauss point gg
        auto t_col_base = ent_data.getFTensor0N(gg, 0);
        // loop over columns
        for (int cc = 0; cc <= rr; ++cc) {
          // calculate element of local matrix
          *a_mat_ptr += alpha * (t_row_base * t_col_base);
          ++t_col_base;
          ++a_mat_ptr;
        }

        ++t_row_base;
      }
      ++t_w; // move to another integration weight
    }

    cholesky_decompose(nN);
  }

  MoFEMFunctionReturn(0);
}

OpBaseDerivativesNext<1>::OpBaseDerivativesNext(
    int derivative, boost::shared_ptr<MatrixDouble> base_mass_ptr,
    boost::shared_ptr<EntitiesFieldData> data_l2,
    const FieldApproximationBase b, const FieldSpace s, int verb, Sev sev)
    : OpBaseDerivativesBase(base_mass_ptr, data_l2, b, s, verb, sev),
      calcBaseDerivative(derivative) {}

template <int SPACE_DIM>
MoFEMErrorCode
OpBaseDerivativesNext<1>::setBase(EntitiesFieldData::EntData &data,
                                  EntitiesFieldData::EntData &ent_data) {
  MoFEMFunctionBegin;

  const int nb_gauss_pts = data.getN(base).size1();
  const int nb_approx_bases = data.getN(base).size2();
  const int nb_derivatives = std::pow(SPACE_DIM, calcBaseDerivative - 1);

  const int nb_prj_bases = ent_data.getN().size2();

  if (!data.getNSharedPtr(base,
                          static_cast<BaseDerivatives>(calcBaseDerivative))) {
    data.getNSharedPtr(base, static_cast<BaseDerivatives>(calcBaseDerivative)) =
        boost::make_shared<MatrixDouble>();
  }

  auto &nex_diff_base = *(data.getNSharedPtr(
      base, static_cast<BaseDerivatives>(calcBaseDerivative)));
  const int next_nb_derivatives = pow(SPACE_DIM, calcBaseDerivative);
  nex_diff_base.resize(nb_gauss_pts, nb_approx_bases * next_nb_derivatives);
  nex_diff_base.clear();

  FTensor::Index<'i', SPACE_DIM> i;
  auto next_diffs_ptr = &*nex_diff_base.data().begin();
  auto t_next_diff = getFTensor1FromPtr<SPACE_DIM>(next_diffs_ptr);

  for (int gg = 0; gg != nb_gauss_pts; ++gg) {

    auto ptr = &*nF.data().begin();
    for (auto r = 0; r != nb_approx_bases * nb_derivatives; ++r) {

        auto l2_diff_base = ent_data.getFTensor1DiffN<SPACE_DIM>(base, gg, 0);
        for (int rr = 0; rr != nb_prj_bases; ++rr) {
          t_next_diff(i) += l2_diff_base(i) * (*ptr);
          ++l2_diff_base;
          ++ptr;
        }

        ++t_next_diff;
    }
  }

  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
OpBaseDerivativesNext<1>::doWork(int side, EntityType type,
                                 EntitiesFieldData::EntData &data) {
  MoFEMFunctionBegin;

  auto &approx_base = data.getN(base);
  const auto nb_approx_bases = approx_base.size2();

  if (nb_approx_bases) {

    const auto fe_type = getFEType();
    const auto nb_integration_pts = approx_base.size1();

    const auto space_dim = data.getDiffN(base).size2() / nb_approx_bases;
    auto &diff_approx_base = *(data.getNSharedPtr(
        base, static_cast<BaseDerivatives>(calcBaseDerivative - 1)));
    int nb_derivatives = pow(space_dim, calcBaseDerivative - 1);

    auto &ent_data = dataL2->dataOnEntities[fe_type][0];
    const int nb_prj_bases = ent_data.getN().size2();

#ifndef NDEBUG
    if (!baseMassPtr)
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Mass matrix is null pointer");
#endif
    auto &nN = *baseMassPtr;

#ifndef NDEBUG
    if (diff_approx_base.size2() != nb_approx_bases * nb_derivatives) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Number of deriveraives and basses do not match");
    }
    if (ent_data.getN().size1() != nb_integration_pts) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Number of integration points is not consistent");
    }
    if (nN.size2() != nb_prj_bases) {
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Number of base functions and size of mass matrix does not math");
    }
#endif


    nF.resize(nb_approx_bases * nb_derivatives, nb_prj_bases, false);
    nF.clear();

    auto t_w = getFTensor0IntegrationWeight();
    // get base function gradient on rows

    auto diff_base_ptr = &*diff_approx_base.data().begin();
    // loop over integration points
    for (int gg = 0; gg != nb_integration_pts; ++gg) {
      // take into account Jacobian
      const double alpha = t_w;

      for (int r = 0; r != nb_approx_bases * nb_derivatives; ++r) {

        // Rows are base functions
        auto t_base = ent_data.getFTensor0N(base, gg, 0);
        for (int rr = 0; rr != nb_prj_bases; ++rr) {
          nF(r, rr) += alpha * (t_base * (*diff_base_ptr));
          ++t_base;
        }

        ++diff_base_ptr;
      }
      
      ++t_w; // move to another integration weight
    }

    for (auto r = 0; r != nb_approx_bases * nb_derivatives; ++r) {
      ublas::matrix_row<MatrixDouble> mc(nF, r);
      cholesky_solve(nN, mc, ublas::lower());
    }

    if (space_dim == 3)
      CHKERR setBase<3>(data, ent_data);
    else if (space_dim == 2)
      CHKERR setBase<2>(data, ent_data);
    // else if (space_dim == 1)
    //   CHKERR setBase<1>(data, ent_data);
    else
      SETERRQ1(PETSC_COMM_SELF, MOFEM_IMPOSIBLE_CASE,
               "Space dim can be only 1,2,3 but is %d", space_dim);
  }

  MoFEMFunctionReturn(0);
}

} // namespace MoFEM