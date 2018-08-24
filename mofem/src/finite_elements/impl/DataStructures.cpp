/** \file DataStructures.cpp

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

#include <Includes.hpp>
#include <version.h>
#include <definitions.h>
#include <Common.hpp>

#include <h1_hdiv_hcurl_l2.h>

#include <UnknownInterface.hpp>

#include <MaterialBlocks.hpp>
#include <BCData.hpp>
#include <TagMultiIndices.hpp>
#include <CoordSysMultiIndices.hpp>
#include <FieldMultiIndices.hpp>
#include <EntsMultiIndices.hpp>
#include <DofsMultiIndices.hpp>
#include <FEMultiIndices.hpp>
#include <ProblemsMultiIndices.hpp>
#include <AdjacencyMultiIndices.hpp>
#include <BCMultiIndices.hpp>
#include <CoreDataStructures.hpp>
#include <SeriesMultiIndices.hpp>

#include <LoopMethods.hpp>
#include <Interface.hpp>
#include <Core.hpp>

#include <FTensor.hpp>
#include <DataStructures.hpp>

#include <boost/scoped_ptr.hpp>

#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

namespace MoFEM {

template <>
FTensor::Tensor0<FTensor::PackPtr<double *, 1> >
getFTensor0FromVec<double, DoubleAllocator>(
    ublas::vector<double, DoubleAllocator> &data) {
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1> >(
      &*data.data().begin());
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>
getFTensor1FromMat<3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 3) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0));
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>
getFTensor1FromMat<2, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 2) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>(&data(0, 0),
                                                            &data(1, 0));
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>
getFTensor2FromMat<3, 3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 9) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 3>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0));
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 2>
getFTensor2FromMat<3, 2, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &data) {
  if (data.size1() != 6) {
    THROW_MESSAGE("Wrong size of data matrix");
  }
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 3, 2>(
      &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
      &data(5, 0));
}

DataForcesAndSourcesCore::EntData::EntData()
    : sEnse(0), oRder(0), bAse(NOBASE) {
  N.resize(LASTBASE);
  diffN.resize(LASTBASE);
  for (ShapeFunctionBasesVector::iterator nit = N.begin(); nit != N.end();
       nit++) {
    nit->reset(new MatrixDouble());
  }
  for (ShapeFunctionBasesVector::iterator nit = diffN.begin();
       nit != diffN.end(); nit++) {
    nit->reset(new MatrixDouble());
  }
}

DataForcesAndSourcesCore::EntData::~EntData() {}

template <class T>
void cOnstructor(DataForcesAndSourcesCore *data, EntityType type, T) {

  data->dataOnEntities[MBENTITYSET].push_back(new T());

  switch (type) {
  case MBTET:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for (int ee = 0; ee < 6; ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    for (int ff = 0; ff < 4; ff++) {
      data->dataOnEntities[MBTRI].push_back(new T());
    }
    data->dataOnEntities[MBTET].push_back(new T());
    break;
  case MBTRI:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for (int ee = 0; ee < 3; ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    data->dataOnEntities[MBTRI].push_back(new T());
    break;
  case MBEDGE:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    data->dataOnEntities[MBEDGE].push_back(new T());
    break;
  case MBVERTEX:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    break;
  case MBPRISM:
    data->dataOnEntities[MBVERTEX].push_back(new T());
    for (int ee = 0; ee < 9; ee++) {
      data->dataOnEntities[MBEDGE].push_back(new T());
    }
    for (int ff = 0; ff < 5; ff++) {
      data->dataOnEntities[MBQUAD].push_back(new T());
    }
    for (int ff = 0; ff < 5; ff++) {
      data->dataOnEntities[MBTRI].push_back(new T());
    }
    data->dataOnEntities[MBPRISM].push_back(new T());
    break;
  default:
    throw MoFEMException(MOFEM_NOT_IMPLEMENTED);
  }
}

DataForcesAndSourcesCore::DataForcesAndSourcesCore(EntityType type) {
  cOnstructor(this, type, EntData());
}

DerivedDataForcesAndSourcesCore::DerivedDataForcesAndSourcesCore(
    DataForcesAndSourcesCore &data)
    : DataForcesAndSourcesCore() {

  boost::ptr_vector<EntData>::iterator iit;

  boost::ptr_vector<EntData>::iterator it;
  for (it = data.dataOnEntities[MBVERTEX].begin();
       it != data.dataOnEntities[MBVERTEX].end(); it++) {
    dataOnEntities[MBVERTEX].push_back(new DerivedEntData(*it));
  }
  for (it = data.dataOnEntities[MBEDGE].begin();
       it != data.dataOnEntities[MBEDGE].end(); it++) {
    dataOnEntities[MBEDGE].push_back(new DerivedEntData(*it));
  }
  for (it = data.dataOnEntities[MBTRI].begin();
       it != data.dataOnEntities[MBTRI].end(); it++) {
    dataOnEntities[MBTRI].push_back(new DerivedEntData(*it));
  }
  for (it = data.dataOnEntities[MBQUAD].begin();
       it != data.dataOnEntities[MBQUAD].end(); it++) {
    dataOnEntities[MBQUAD].push_back(new DerivedEntData(*it));
  }
  for (it = data.dataOnEntities[MBTET].begin();
       it != data.dataOnEntities[MBTET].end(); it++) {
    dataOnEntities[MBTET].push_back(new DerivedEntData(*it));
  }
  for (it = data.dataOnEntities[MBPRISM].begin();
       it != data.dataOnEntities[MBPRISM].end(); it++) {
    dataOnEntities[MBPRISM].push_back(new DerivedEntData(*it));
  }
}

std::ostream &operator<<(std::ostream &os,
                         const DataForcesAndSourcesCore::EntData &e) {
  os << "sEnse: " << e.getSense() << std::endl
     << "oRder: " << e.getOrder() << std::endl
     << "global indices: " << e.getIndices() << std::endl
     << "local indices: " << e.getLocalIndices() << std::endl;
  os.precision(2);
  os << "fieldData: " << std::fixed << e.getFieldData() << std::endl;
  MatrixDouble base = e.getN();
  MatrixDouble diff_base = e.getDiffN();
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

std::ostream &operator<<(std::ostream &os, const DataForcesAndSourcesCore &e) {
  for (unsigned int nn = 0; nn < e.dataOnEntities[MBVERTEX].size(); nn++) {
    os << "dataOnEntities[MBVERTEX][" << nn << "]" << std::endl
       << e.dataOnEntities[MBVERTEX][nn] << std::endl;
  }
  for (unsigned int ee = 0; ee < e.dataOnEntities[MBEDGE].size(); ee++) {
    os << "dataOnEntities[MBEDGE][" << ee << "]" << std::endl
       << e.dataOnEntities[MBEDGE][ee] << std::endl;
  }
  for (unsigned int ff = 0; ff < e.dataOnEntities[MBTRI].size(); ff++) {
    os << "dataOnEntities[MBTRI][" << ff << "] " << std::endl
       << e.dataOnEntities[MBTRI][ff] << std::endl;
  }
  for (unsigned int vv = 0; vv < e.dataOnEntities[MBTET].size(); vv++) {
    os << "dataOnEntities[MBTET][" << vv << "]" << std::endl
       << e.dataOnEntities[MBTET][vv] << std::endl;
  }
  return os;
}

/** \name Specializations for H1/L2 */

/**@{*/

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
DataForcesAndSourcesCore::EntData::getFTensor1FieldData<3>() {
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
DataForcesAndSourcesCore::EntData::getFTensor1FieldData<2>() {
  if (dOfs[0]->getNbOfCoeffs() != 2) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but you ask for tensor rank 1 dimension 3";
    THROW_MESSAGE(s.str());
  }
  double *ptr = &*fieldData.data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>(ptr, &ptr[1]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
DataForcesAndSourcesCore::EntData::getFTensor2FieldData<3, 3>() {
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

FTensor::Tensor0<FTensor::PackPtr<double *,1> >
DataForcesAndSourcesCore::EntData::getFTensor0FieldData() {
  if (dOfs[0]->getNbOfCoeffs() != 1) {
    std::stringstream s;
    s << "Wrong number of coefficients is " << dOfs[0]->getNbOfCoeffs();
    s << " but expected scalar field, tensor of rank 0";
    THROW_MESSAGE(s.str());
  }
  return FTensor::Tensor0<FTensor::PackPtr<double *, 1> >(
      &*fieldData.data().begin());
}

template <int Tensor_Dim>
FTensor::Tensor1<double *, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN(
    const FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<double *, Tensor_Dim>();
}

template <int Tensor_Dim>
FTensor::Tensor1<double *, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN() {
  return getFTensor1DiffN<Tensor_Dim>(bAse);
}

template <int Tensor_Dim>
FTensor::Tensor1<double *, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN(
    const FieldApproximationBase base, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<double *, Tensor_Dim>();
}

template <int Tensor_Dim>
FTensor::Tensor1<double *, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN(const int bb) {
  return getFTensor1DiffN<Tensor_Dim>(bAse, bb);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2], 3);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>() {
  return getFTensor1DiffN<3>(bAse);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base, const int bb) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<double *, 3>(
      &ptr[3 * bb], &ptr[3 * bb + 1], &ptr[3 * bb + 2], getDiffN(base).size2());
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>(const int bb) {
  return getFTensor1DiffN<3>(bAse, bb);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<double *, 2>(ptr, &ptr[1], 2);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>() {
  return getFTensor1DiffN<2>(bAse);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base, const int bb) {
  double *ptr = &*getDiffN(base).data().begin();
  return FTensor::Tensor1<double *, 2>(&ptr[2 * bb], &ptr[2 * bb + 1],
                                       getDiffN(base).size1());
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>(const int bb) {
  return getFTensor1DiffN<2>(bAse, bb);
}

template <int Tensor_Dim>
FTensor::Tensor1<double *, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN(
    const FieldApproximationBase base, const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<double *, Tensor_Dim>();
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>(
    const FieldApproximationBase base, const int gg, const int bb) {
  double *ptr = &getDiffN(base)(gg, 3 * bb);
  return FTensor::Tensor1<double *, 3>(ptr, &ptr[1], &ptr[2], 3);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 3d
 */
template <>
FTensor::Tensor1<double *, 3>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<3>(const int gg,
                                                       const int bb) {
  return getFTensor1DiffN<3>(bAse, gg, bb);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>(
    const FieldApproximationBase base, const int gg, const int bb) {
  double *ptr = &getDiffN(base)(gg, 2 * bb);
  return FTensor::Tensor1<double *, 2>(ptr, &ptr[1], 2);
}

/**
 * \brief Get spatial derivative of base function tensor for dimension 2d
 */
template <>
FTensor::Tensor1<double *, 2>
DataForcesAndSourcesCore::EntData::getFTensor1DiffN<2>(const int gg,
                                                       const int bb) {
  return getFTensor1DiffN<2>(bAse, gg, bb);
}

/**@}*/

/** \name Specializations for HDiv/HCrul */

/**@{*/

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1HdivN(
    FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, Tensor_Dim>();
}

template <int Tensor_Dim>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, Tensor_Dim>
DataForcesAndSourcesCore::EntData::getFTensor1HdivN(FieldApproximationBase base,
                                                    const int gg,
                                                    const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, Tensor_Dim>();
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN(
    FieldApproximationBase base) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << "x" << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<double *, Tensor_Dim0, Tensor_Dim1>();
}

template <int Tensor_Dim0, int Tensor_Dim1>
FTensor::Tensor2<FTensor::PackPtr<double *, Tensor_Dim0 * Tensor_Dim1>,
                 Tensor_Dim0, Tensor_Dim1>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN(
    FieldApproximationBase base, const int gg, const int bb) {
  std::stringstream s;
  s << "Template for tensor dimension " << Tensor_Dim0 << "x" << Tensor_Dim1
    << " not implemented";
  THROW_MESSAGE(s.str());
  return FTensor::Tensor2<double *, Tensor_Dim0, Tensor_Dim1>();
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
DataForcesAndSourcesCore::EntData::getFTensor1HdivN<3>(
    FieldApproximationBase base) {
  double *t_n_ptr = &*getHdivN(base).data().begin();
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(t_n_ptr, // HDIV0
                                                            &t_n_ptr[HDIV1],
                                                            &t_n_ptr[HDIV2]);
}

template <>
FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
DataForcesAndSourcesCore::EntData::getFTensor1HdivN<3>(
    FieldApproximationBase base, const int gg, const int bb) {
  double *t_n_ptr = &getHdivN(base)(gg, 3 * bb);
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>(t_n_ptr, // HDIV0
                                                            &t_n_ptr[HDIV1],
                                                            &t_n_ptr[HDIV2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN<3, 3>(
    FieldApproximationBase base) {
  double *t_diff_n_ptr = &*getDiffHdivN(base).data().begin();
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      t_diff_n_ptr, &t_diff_n_ptr[HDIV0_1], &t_diff_n_ptr[HDIV0_2],
      &t_diff_n_ptr[HDIV1_0], &t_diff_n_ptr[HDIV1_1], &t_diff_n_ptr[HDIV1_2],
      &t_diff_n_ptr[HDIV2_0], &t_diff_n_ptr[HDIV2_1], &t_diff_n_ptr[HDIV2_2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN<3, 3>(
    FieldApproximationBase base, const int gg, const int bb) {
  double *t_diff_n_ptr = &getDiffHdivN(base)(gg, 9 * bb);
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      t_diff_n_ptr, &t_diff_n_ptr[HDIV0_1], &t_diff_n_ptr[HDIV0_2],
      &t_diff_n_ptr[HDIV1_0], &t_diff_n_ptr[HDIV1_1], &t_diff_n_ptr[HDIV1_2],
      &t_diff_n_ptr[HDIV2_0], &t_diff_n_ptr[HDIV2_1], &t_diff_n_ptr[HDIV2_2]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN<3, 2>(
    FieldApproximationBase base) {
  double *t_diff_n_ptr = &*getDiffHdivN(base).data().begin();
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      t_diff_n_ptr, &t_diff_n_ptr[HCURL0_1], &t_diff_n_ptr[HCURL1_0],
      &t_diff_n_ptr[HCURL1_1], &t_diff_n_ptr[HCURL2_0],
      &t_diff_n_ptr[HCURL2_1]);
}

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>
DataForcesAndSourcesCore::EntData::getFTensor2DiffHdivN<3, 2>(
    FieldApproximationBase base, const int gg, const int bb) {
  double *t_diff_n_ptr = &getDiffHdivN(base)(gg, 6 * bb);
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      t_diff_n_ptr, &t_diff_n_ptr[HCURL0_1], &t_diff_n_ptr[HCURL1_0],
      &t_diff_n_ptr[HCURL1_1], &t_diff_n_ptr[HCURL2_0], &t_diff_n_ptr[HCURL2_1]);
}

/**@}*/

} // namespace MoFEM
