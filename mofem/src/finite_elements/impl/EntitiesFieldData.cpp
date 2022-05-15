/** \file EntitiesFieldData.cpp
\brief Implementation for Data Structures in Forces and Sources

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

EntitiesFieldData::EntData::EntData(const bool allocate_base_matrices)
    : sEnse(0), oRder(0), bAse(NOBASE), entDataBitRefLevel(),
      N(baseFunctionsAndBaseDerivatives[ZeroDerivative]),
      diffN(baseFunctionsAndBaseDerivatives[FirstDerivative]) {
  if (allocate_base_matrices) {

    for (auto d = 0; d != LastDerivative; ++d) {
      for (int b = 0; b != LASTBASE; ++b) {
        baseFunctionsAndBaseDerivatives[d][b].reset(new MatrixDouble());
        N[b].reset(new MatrixDouble());
        diffN[b].reset(new MatrixDouble());
      }
    }
  }
}

int EntitiesFieldData::EntData::getSense() const { return sEnse; }

boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getNSharedPtr(const FieldApproximationBase base,
                                          const BaseDerivatives direvatie) {
  return baseFunctionsAndBaseDerivatives[direvatie][base];
}

boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getNSharedPtr(const FieldApproximationBase base) {
  return N[base];
}

boost::shared_ptr<MatrixDouble> &EntitiesFieldData::EntData::getDiffNSharedPtr(
    const FieldApproximationBase base) {
  return diffN[base];
}

static void constructor_data(EntitiesFieldData *data, const EntityType type) {

  using EntData = EntitiesFieldData::EntData;

  auto set_default = [&]() {
    std::array<size_t, MBMAXTYPE> count;
    std::fill(count.begin(), count.end(), 0);
    const int dim_type = moab::CN::Dimension(type);
    data->dataOnEntities[MBVERTEX].resize(1);
    if (type != MBVERTEX) {
      for (auto dd = dim_type; dd > 0; --dd) {
        int nb_ents = moab::CN::NumSubEntities(type, dd);
        for (int ii = 0; ii != nb_ents; ++ii) {
          auto sub_ent_type = moab::CN::SubEntityType(type, dd, ii);
          count[sub_ent_type] = nb_ents;
        }
        for (auto tt = moab::CN::TypeDimensionMap[dd].first;
             tt <= moab::CN::TypeDimensionMap[dd].second; ++tt) {
          data->dataOnEntities[tt].resize(count[tt]);
        }
      }
    }
  };

  switch (type) {
  case MBENTITYSET:
    break;

  default:
    set_default();
  }
}

EntitiesFieldData::EntitiesFieldData(EntityType type) {
  constructor_data(this, type);
}

MoFEMErrorCode EntitiesFieldData::setElementType(const EntityType type) {
  MoFEMFunctionBegin;
  constructor_data(this, type);
  MoFEMFunctionReturn(0);
}

static void
constructor_derived_data(DerivedEntitiesFieldData *derived_data,
                         const boost::shared_ptr<EntitiesFieldData> &data_ptr) {

  using EntData = EntitiesFieldData::EntData;
  using DerivedEntData = DerivedEntitiesFieldData::DerivedEntData;

  for (int tt = MBVERTEX; tt != MBMAXTYPE; ++tt) {
    auto &ent_data = data_ptr->dataOnEntities[tt];
    auto &derived_ent_data = derived_data->dataOnEntities[tt];
    for (auto c = derived_ent_data.size(); c < ent_data.size(); ++c) {
      boost::shared_ptr<EntData> ent_data_ptr(data_ptr, &ent_data[c]);
      derived_ent_data.push_back(new DerivedEntData(ent_data_ptr));
    }
    derived_ent_data.resize(ent_data.size());
  }
}

DerivedEntitiesFieldData::DerivedEntitiesFieldData(
    const boost::shared_ptr<EntitiesFieldData> &data_ptr)
    : EntitiesFieldData(), dataPtr(data_ptr) {
  constructor_derived_data(this, dataPtr);
}

MoFEMErrorCode DerivedEntitiesFieldData::setElementType(const EntityType type) {
  MoFEMFunctionBegin;
  constructor_derived_data(this, dataPtr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EntitiesFieldData::EntData::resetFieldDependentData() {
  MoFEMFunctionBeginHot;
  sPace = NOSPACE;
  bAse = NOBASE;
  iNdices.resize(0, false);
  localIndices.resize(0, false);
  dOfs.resize(0, false);
  fieldData.resize(0, false);
  MoFEMFunctionReturnHot(0);
}

MoFEMErrorCode EntitiesFieldData::resetFieldDependentData() {
  MoFEMFunctionBegin;
  for (EntityType t = MBVERTEX; t != MBMAXTYPE; t++)
    for (auto &e : dataOnEntities[t])
      CHKERR e.resetFieldDependentData();
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode
EntitiesFieldData::EntData::baseSwap(const std::string &field_name,
                                     const FieldApproximationBase base) {
  MoFEMFunctionBegin;
  auto make_swap = [](boost::shared_ptr<MatrixDouble> &ptr,
                      boost::shared_ptr<MatrixDouble> &ptrBB,
                      boost::shared_ptr<MatrixDouble> &swap_ptr) {
    if (swap_ptr) {
      ptr = swap_ptr;
      swap_ptr.reset();
    } else {
      swap_ptr = ptr;
      ptr = ptrBB;
    }
  };
  make_swap(getNSharedPtr(base), getBBNSharedPtr(field_name), swapBaseNPtr);
  make_swap(getDiffNSharedPtr(base), getBBDiffNSharedPtr(field_name),
            swapBaseDiffNPtr);
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode EntitiesFieldData::baseSwap(const std::string &field_name,
                                           const FieldApproximationBase base) {
  MoFEMFunctionBegin;
  // Note: Do not swap bases on entities sets
  for (int tt = MBVERTEX; tt != MBENTITYSET; ++tt) {
    auto &ent_data = dataOnEntities[tt];
    for (auto &side_data : ent_data)
      CHKERR side_data.baseSwap(field_name, base);
  }
  MoFEMFunctionReturn(0);
}

MoFEMErrorCode DerivedEntitiesFieldData::DerivedEntData::baseSwap(
    const std::string &field_name, const FieldApproximationBase base) {
  MoFEMFunctionBegin;
  auto make_swap = [](boost::shared_ptr<MatrixDouble> &ptr,
                      boost::shared_ptr<MatrixDouble> &ptrBB,
                      boost::shared_ptr<MatrixDouble> &swap_ptr) {
    if (swap_ptr) {
      ptr = swap_ptr;
      swap_ptr.reset();
    } else {
      swap_ptr = ptr;
      ptr = ptrBB;
    }
  };
  make_swap(getDerivedNSharedPtr(base), getBBNSharedPtr(field_name),
            swapBaseNPtr);
  make_swap(getDerivedDiffNSharedPtr(base), getBBDiffNSharedPtr(field_name),
            swapBaseDiffNPtr);
  MoFEMFunctionReturn(0);
}

DerivedEntitiesFieldData::DerivedEntData::DerivedEntData(
    const boost::shared_ptr<EntitiesFieldData::EntData> &ent_data_ptr)
    : EntitiesFieldData::EntData(false), entDataPtr(ent_data_ptr) {}

int DerivedEntitiesFieldData::DerivedEntData::getSense() const {
  return entDataPtr->getSense();
}

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getNSharedPtr(
    const FieldApproximationBase base, const BaseDerivatives derivative) {
  if (baseFunctionsAndBaseDerivatives[derivative][base])
    return baseFunctionsAndBaseDerivatives[derivative][base];
  else
    return entDataPtr->getNSharedPtr(base, derivative);
}

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getNSharedPtr(
    const FieldApproximationBase base) {
  if (N[base])
    return N[base];
  else
    return entDataPtr->getNSharedPtr(base);
}

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getDiffNSharedPtr(
    const FieldApproximationBase base) {
  if (diffN[base])
    return diffN[base];
  else
    return entDataPtr->getDiffNSharedPtr(base);
}
const boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getNSharedPtr(
    const FieldApproximationBase base) const {
  if (N[base])
    return N[base];
  else
    return entDataPtr->getNSharedPtr(base);
}
const boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getDiffNSharedPtr(
    const FieldApproximationBase base) const {
  if (diffN[base])
    return diffN[base];
  else
    return entDataPtr->getDiffNSharedPtr(base);
}

std::ostream &operator<<(std::ostream &os,
                         const EntitiesFieldData::EntData &e) {
  os << "sEnse: " << e.getSense() << std::endl
     << "oRder: " << e.getOrder() << std::endl
     << "global indices: " << e.getIndices() << std::endl
     << "local indices: " << e.getLocalIndices() << std::endl;
  // FIXME: precision should not be set here
  os << "fieldData: " << std::fixed << std::setprecision(2) << e.getFieldData()
     << std::endl;
  MatrixDouble base = const_cast<EntitiesFieldData::EntData &>(e).getN();
  MatrixDouble diff_base =
      const_cast<EntitiesFieldData::EntData &>(e).getDiffN();
  const double eps = 1e-6;
  for (unsigned int ii = 0; ii != base.size1(); ii++) {
    for (unsigned int jj = 0; jj != base.size2(); jj++) {
      if (fabs(base(ii, jj)) < eps)
        base(ii, jj) = 0;
    }
  }
  for (unsigned int ii = 0; ii != diff_base.size1(); ii++) {
    for (unsigned int jj = 0; jj != diff_base.size2(); jj++) {
      if (fabs(diff_base(ii, jj)) < eps)
        diff_base(ii, jj) = 0;
    }
  }
  os << "N: " << std::fixed << base << std::endl
     << "diffN: " << std::fixed << diff_base;
  return os;
}

std::ostream &operator<<(std::ostream &os, const EntitiesFieldData &e) {
  for (EntityType t = MBVERTEX; t != MBENTITYSET; ++t) {
    for (unsigned int nn = 0; nn < e.dataOnEntities[t].size(); nn++) {
      os << "dataOnEntities[" << moab::CN::EntityTypeName(t) << "][" << nn
         << "]" << std::endl
         << e.dataOnEntities[t][nn] << std::endl;
    }
  }
  return os;
}

/** \name Specializations for H1/L2 */

/**@{*/

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1FieldData<3>() {
  if (dOfs[0]->getNbOfCoeffs() != 3) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for tensor rank 1 dimension 3";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1FieldData<2>() {
  if (dOfs[0]->getNbOfCoeffs() != 2) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for tensor rank 1 dimension 2";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>(ptr, &ptr[1]);
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 1>
EntitiesFieldData::EntData::getFTensor1FieldData<1>() {
  if (dOfs[0]->getNbOfCoeffs() != 1) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for tensor rank 1 dimension 1";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 1>(ptr);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2FieldData<3, 3>() {
  if (dOfs[0]->getNbOfCoeffs() != 9) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for tensor rank 2 dimensions 3 by 3 so 9 coefficients "
         "is expected";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      ptr, &ptr[1], &ptr[2], &ptr[3], &ptr[4], &ptr[5], &ptr[6], &ptr[7],
      &ptr[8]);
}

template <>
FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>
EntitiesFieldData::EntData::getFTensor2SymmetricFieldData<3>() {
  if (dOfs[0]->getNbOfCoeffs() != 6) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for symmetric tensor rank 2 dimensions 3 by 3 so 6 "
         "coefficients "
         "is expected";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>(
      ptr, &ptr[1], &ptr[2], &ptr[3], &ptr[4], &ptr[5]);
}

template <>
FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>
EntitiesFieldData::EntData::getFTensor2SymmetricFieldData<2>() {
  if (dOfs[0]->getNbOfCoeffs() != 3) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for symmetric tensor rank 2 dimensions 2 by 2 so 3 "
         "coefficients "
         "is expected";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>(
      ptr, &ptr[1], &ptr[2]);
}

FTensor::Tensor0<FTensor::PackPtr<double *, 1>>
EntitiesFieldData::EntData::getFTensor0FieldData() {
  if (dOfs[0]->getNbOfCoeffs() != 1) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but expected scalar field, tensor of rank 0";
    THROW_MESSAGE(s.str());
  }
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1>>(
      &*fieldData.data().begin());
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1DiffN(
    const FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<double *, Tensor_Dim>();
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1DiffN() {
  return getFTensor1DiffN<Tensor_Dim>(bAse);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>() {
  return getFTensor1DiffN<3>(bAse);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>(ptr, &ptr[1]);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>() {
  return getFTensor1DiffN<2>(bAse);
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1DiffN(const FieldApproximationBase base,
                                             const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>();
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base, const int gg, const int bb) {
  double *ptr = &getDiffN(base)(gg, 3 * bb);
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(ptr, &ptr[1],
                                                            &ptr[2]);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1DiffN<3>(const int gg, const int bb) {
  return getFTensor1DiffN<3>(bAse, gg, bb);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base, const int gg, const int bb) {
  double *ptr = &getDiffN(base)(gg, 2 * bb);
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>(ptr, &ptr[1]);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
EntitiesFieldData::EntData::getFTensor1DiffN<2>(const int gg, const int bb) {
  return getFTensor1DiffN<2>(bAse, gg, bb);
}

/**@}*/

/** \name Specializations for HDiv/HCrul */

/**@{*/

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1N(FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>();
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>
EntitiesFieldData::EntData::getFTensor1N(FieldApproximationBase base,
                                         const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<FTensor::PackPtr<double *, Tensor_Dim>, Tensor_Dim>();
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
EntitiesFieldData::EntData::getFTensor2DiffN(FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << "x" << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<double *, Tensor_Dim0, Tensor_Dim1>();
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
EntitiesFieldData::EntData::getFTensor2DiffN(FieldApproximationBase base,
                                             const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << "x" << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<double *, Tensor_Dim0, Tensor_Dim1>();
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1N<3>(FieldApproximationBase base) {
  double *t_n_ptr = &*getN(base).data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(t_n_ptr, // HVEC0
                                                            &t_n_ptr[HVEC1],
                                                            &t_n_ptr[HVEC2]);
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
EntitiesFieldData::EntData::getFTensor1N<3>(FieldApproximationBase base,
                                            const int gg, const int bb) {
  double *t_n_ptr = &getN(base)(gg, 3 * bb);
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(t_n_ptr, // HVEC0
                                                            &t_n_ptr[HVEC1],
                                                            &t_n_ptr[HVEC2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 3>(
    FieldApproximationBase base) {
  double *t_diff_n_ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      t_diff_n_ptr, &t_diff_n_ptr[HVEC0_1], &t_diff_n_ptr[HVEC0_2],
      &t_diff_n_ptr[HVEC1_0], &t_diff_n_ptr[HVEC1_1], &t_diff_n_ptr[HVEC1_2],
      &t_diff_n_ptr[HVEC2_0], &t_diff_n_ptr[HVEC2_1], &t_diff_n_ptr[HVEC2_2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 3>(FieldApproximationBase base,
                                                   const int gg, const int bb) {
  double *t_diff_n_ptr = &getDiffN(base)(gg, 9 * bb);
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      t_diff_n_ptr, &t_diff_n_ptr[HVEC0_1], &t_diff_n_ptr[HVEC0_2],
      &t_diff_n_ptr[HVEC1_0], &t_diff_n_ptr[HVEC1_1], &t_diff_n_ptr[HVEC1_2],
      &t_diff_n_ptr[HVEC2_0], &t_diff_n_ptr[HVEC2_1], &t_diff_n_ptr[HVEC2_2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 2>(
    FieldApproximationBase base) {
  double *t_diff_n_ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      t_diff_n_ptr, &t_diff_n_ptr[HVEC0_1], &t_diff_n_ptr[HVEC1_0],
      &t_diff_n_ptr[HVEC1_1], &t_diff_n_ptr[HVEC2_0], &t_diff_n_ptr[HVEC2_1]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
EntitiesFieldData::EntData::getFTensor2DiffN<3, 2>(FieldApproximationBase base,
                                                   const int gg, const int bb) {
  double *t_diff_n_ptr = &getDiffN(base)(gg, 6 * bb);
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      t_diff_n_ptr, &t_diff_n_ptr[HVEC0_1], &t_diff_n_ptr[HVEC1_0],
      &t_diff_n_ptr[HVEC1_1], &t_diff_n_ptr[HVEC2_0], &t_diff_n_ptr[HVEC2_1]);
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
EntitiesFieldData::EntData::getFTensor2N(FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << ", " << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                          Tensor_Dim0, Tensor_Dim1>();
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2N<3, 3>(FieldApproximationBase base) {
  double *t_n_ptr = &*(getN(base).data().begin());
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(

      &t_n_ptr[HVEC0], &t_n_ptr[HVEC1], &t_n_ptr[HVEC2],

      &t_n_ptr[3 + HVEC0], &t_n_ptr[3 + HVEC1], &t_n_ptr[3 + HVEC2],

      &t_n_ptr[6 + HVEC0], &t_n_ptr[6 + HVEC1], &t_n_ptr[6 + HVEC2]

  );
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
EntitiesFieldData::EntData::getFTensor2N(FieldApproximationBase base,
                                         const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << ", " << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                          Tensor_Dim0, Tensor_Dim1>();
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
EntitiesFieldData::EntData::getFTensor2N<3, 3>(FieldApproximationBase base,
                                               const int gg, const int bb) {
  double *t_n_ptr = &getN(base)(gg, 9 * bb);
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(

      &t_n_ptr[HVEC0], &t_n_ptr[HVEC1], &t_n_ptr[HVEC2],

      &t_n_ptr[3 + HVEC0], &t_n_ptr[3 + HVEC1], &t_n_ptr[3 + HVEC2],

      &t_n_ptr[6 + HVEC0], &t_n_ptr[6 + HVEC1], &t_n_ptr[6 + HVEC2]

  );
}

/**@}*/

/** \name Bernstein-Bezier base only functions */

/**@{*/

boost::shared_ptr<MatrixInt> &
EntitiesFieldData::EntData::getBBAlphaIndicesSharedPtr(
    const std::string &field_name) {
  return bbAlphaIndices[field_name];
}

boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBNSharedPtr(const std::string &field_name) {
  return bbN[field_name];
}

/**
 * Get shared pointer to BB base base functions
 */
const boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBNSharedPtr(
    const std::string &field_name) const {
  return bbN.at(field_name);
}

/**
 * Get shared pointer to BB derivatives of base base functions
 */
boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBDiffNSharedPtr(const std::string &field_name) {
  return bbDiffN[field_name];
}

/**
 * Get shared pointer to derivatives of BB base base functions
 */
const boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBDiffNSharedPtr(
    const std::string &field_name) const {
  return bbDiffN.at(field_name);
}

std::map<std::string, boost::shared_ptr<MatrixInt>> &
EntitiesFieldData::EntData::getBBAlphaIndicesMap() {
  return bbAlphaIndices;
}

std::map<std::string, boost::shared_ptr<MatrixDouble>> &
EntitiesFieldData::EntData::getBBNMap() {
  return bbN;
}

std::map<std::string, boost::shared_ptr<MatrixDouble>> &
EntitiesFieldData::EntData::getBBDiffNMap() {
  return bbDiffN;
}

boost::shared_ptr<MatrixInt> &
EntitiesFieldData::EntData::getBBAlphaIndicesByOrderSharedPtr(const size_t o) {
  return bbAlphaIndicesByOrder[o];
}

boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBNByOrderSharedPtr(const size_t o) {
  return bbNByOrder[o];
}

boost::shared_ptr<MatrixDouble> &
EntitiesFieldData::EntData::getBBDiffNByOrderSharedPtr(const size_t o) {
  return bbDiffNByOrder[o];
}

std::array<boost::shared_ptr<MatrixInt>,
           EntitiesFieldData::EntData::MaxBernsteinBezierOrder> &
EntitiesFieldData::EntData::getBBAlphaIndicesByOrderArray() {
  return bbAlphaIndicesByOrder;
}

std::array<boost::shared_ptr<MatrixDouble>,
           EntitiesFieldData::EntData::MaxBernsteinBezierOrder> &
EntitiesFieldData::EntData::getBBNByOrderArray() {
  return bbNByOrder;
}

std::array<boost::shared_ptr<MatrixDouble>,
           EntitiesFieldData::EntData::MaxBernsteinBezierOrder> &
EntitiesFieldData::EntData::getBBDiffNByOrderArray() {
  return bbDiffNByOrder;
}

boost::shared_ptr<MatrixInt> &
DerivedEntitiesFieldData::DerivedEntData::getBBAlphaIndicesSharedPtr(
    const std::string &field_name) {
  return entDataPtr->getBBAlphaIndicesSharedPtr(field_name);
}

boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getBBNSharedPtr(
    const std::string &field_name) {
  return entDataPtr->getBBNSharedPtr(field_name);
}

const boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getBBNSharedPtr(
    const std::string &field_name) const {
  return entDataPtr->getBBNSharedPtr(field_name);
}

/**
 * Get shared pointer to BB derivatives of base base functions
 */
boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getBBDiffNSharedPtr(
    const std::string &field_name) {
  return entDataPtr->getBBDiffNSharedPtr(field_name);
}

/**
 * Get shared pointer to derivatives of BB base base functions
 */
const boost::shared_ptr<MatrixDouble> &
DerivedEntitiesFieldData::DerivedEntData::getBBDiffNSharedPtr(
    const std::string &field_name) const {
  return entDataPtr->getBBDiffNSharedPtr(field_name);
}

/**@}*/

BitRefLevel &EntitiesFieldData::EntData::getEntDataBitRefLevel() {
  return entDataBitRefLevel;
}

BitRefLevel &DerivedEntitiesFieldData::DerivedEntData::getEntDataBitRefLevel() {
  return entDataPtr->getEntDataBitRefLevel();
}


} // namespace MoFEM
