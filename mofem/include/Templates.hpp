/** \file Templates.hpp
 * \brief Templates declarations
 */

#ifndef __TEMPLATES_HPP__
#define __TEMPLATES_HPP__

namespace MoFEM {

template <typename T> using ShardVec = boost::shared_ptr<std::vector<T>>;

/**
 * @brief Get Vector adaptor
 *
 * \code
 *
 * double *a;
 * CHKERR VecGetArray(v,&a);
 *
 * for(int n = 0; n != nodes; ++n) {
 *
 *   auto a = getVectorAdaptor(&a[3*n], 3);
 *   double dot = inner_prod(a, a);
 *
 * }
 *
 * CHKERR VecRetsoreArray(v,&a);
 * \endcode
 *
 */
template <typename T1> inline auto getVectorAdaptor(T1 ptr, const size_t n) {
  typedef typename std::remove_pointer<T1>::type T;
  return VectorShallowArrayAdaptor<T>(n,
                                      ublas::shallow_array_adaptor<T>(n, ptr));
};

/**
 * @brief Get Matrix adaptor
 *
 * \code
 *
 * double *a;
 * CHKERR VecGetArray(v,&a);
 *
 * for(int n = 0; n != nodes; ++n) {
 *
 *   auto F = getMatrixAdaptor(&a[3*3*n], 3, 3);
 *   MatrixDouble C = prod(F, trans(F));
 *
 * }
 *
 * CHKERR VecRetsoreArray(v,&a);
 * \endcode
 *
 */
template <typename T1>
inline auto getMatrixAdaptor(T1 ptr, const size_t n, const size_t m) {
  typedef typename std::remove_pointer<T1>::type T;
  return MatrixShallowArrayAdaptor<T>(
      n, m, ublas::shallow_array_adaptor<T>(n * m, ptr));
};

/**
 * This small utility that cascades two key extractors will be
 * used throughout the boost example
 * <a
 * href=http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp>
 * http://www.boost.org/doc/libs/1_53_0/libs/multi_index/example/complex_structs.cpp
 * </a>
 */
template <class KeyExtractor1, class KeyExtractor2> struct KeyFromKey {
public:
  typedef typename KeyExtractor1::result_type result_type;

  KeyFromKey(const KeyExtractor1 &key1_ = KeyExtractor1(),
             const KeyExtractor2 &key2_ = KeyExtractor2())
      : key1(key1_), key2(key2_) {}

  template <typename Arg> result_type operator()(Arg &arg) const {
    return key1(key2(arg));
  }

private:
  KeyExtractor1 key1;
  KeyExtractor2 key2;
};

template <typename id_type> struct LtBit {
  inline bool operator()(const id_type &valueA, const id_type &valueB) const {
    return valueA.to_ulong() < valueB.to_ulong();
  }
};

template <typename id_type> struct EqBit {
  inline bool operator()(const id_type &valueA, const id_type &valueB) const {
    return valueA.to_ulong() == valueB.to_ulong();
  }
};

template <typename id_type> struct HashBit {
  inline unsigned int operator()(const id_type &value) const {
    return value.to_ulong();
  }
};

template <class X> inline std::string toString(X x) {
  std::ostringstream buffer;
  buffer << x;
  return buffer.str();
}

template <int S, class T, class A> struct GetFTensor0FromVecImpl {
  static inline auto get(ublas::vector<T, A> &data) {
    return FTensor::Tensor0<FTensor::PackPtr<T *, S>>(&*data.data().begin());
  }
};

/**
* \brief Get tensor rank 0 (scalar) form data vector

Example how to use it.
\code
VectorDouble vec;
vec.resize(nb_gauss_pts,false);
vec.clear();
auto t0 = getFTensor0FromData(data);
for(int gg = 0;gg!=nb_gauss_pts;gg++) {

  ++t0;
}
\endcode

*/
template <int S = 1, class T, class A>
static inline auto getFTensor0FromVec(ublas::vector<T, A> &data) {
  return GetFTensor0FromVecImpl<S, T, A>::get(data);
}

template <int Tensor_Dim, int S, class T, class L, class A>
struct GetFTensor1FromMatImpl {};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDDEBUG
    if (data.size1() != 3)
      THROW_MESSAGE(
          "getFTensor1FromMat<3>: wrong size of data matrix, number of "
          "rows should be 3 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 3>(
        &data(0, 0), &data(1, 0), &data(2, 0));
  }
};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<4, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDDEBUG
    if (data.size1() != 4)
      THROW_MESSAGE(
          "getFTensor1FromMat<4>: wrong size of data matrix, number of "
          "rows should be 4 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 4>(
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0));
  }
};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<6, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDDEBUG
    if (data.size1() != 6)
      THROW_MESSAGE(
          "getFTensor1FromMat<6>: wrong size of data matrix, number of "
          "rows should be 6 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 6>(
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
        &data(5, 0));
  }
};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<9, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDDEBUG
    if (data.size1() != 9)
      THROW_MESSAGE(
          "getFTensor1FromMat<9>: wrong size of data matrix, number of "
          "rows should be 9 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 9>(
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
        &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0));
  }
};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDDEBUG
    if (data.size1() != 2)
      THROW_MESSAGE(
          "getFTensor1FromMat<2>: wrong size of data matrix, number of "
          "rows should be 2 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 2>(&data(0, 0),
                                                              &data(1, 0));
  }
};

template <int S, class T, class A>
struct GetFTensor1FromMatImpl<1, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 1)
      THROW_MESSAGE(
          "getFTensor1FromMat<1>: wrong size of data matrix, number of "
          "rows should be 1 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor1<FTensor::PackPtr<double *, S>, 1>(&data(0, 0));
  }
};

/**
 * \brief Get tensor rank 1 (vector) form data matrix
 */
template <int Tensor_Dim, int S = 1, class T, class L, class A>
inline FTensor::Tensor1<FTensor::PackPtr<T *, S>, Tensor_Dim>
getFTensor1FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value, "not implemented");
}

/**
 * \brief Get tensor rank 1 (vector) form data matrix (specialization)
 */
template <int Tensor_Dim, int S = 1>
inline auto getFTensor1FromMat(MatrixDouble &data) {
  return GetFTensor1FromMatImpl<Tensor_Dim, S, double, ublas::row_major,
                                DoubleAllocator>::get(data);
}

template <int Tensor_Dim1, int Tensor_Dim2, int S, class T, class L, class A>
struct GetFTensor2FromMatImpl {
  static inline auto get(ublas::matrix<T, L, A> &data) {
#ifndef NDEBUG
    if (data.size1() != Tensor_Dim1 * Tensor_Dim2) {
      THROW_MESSAGE(
          "getFTensor2FromMat<" +
          boost::lexical_cast<std::string>(Tensor_Dim1) + "," +
          boost::lexical_cast<std::string>(Tensor_Dim2) +
          ">: wrong size of rows in data matrix, should be " +
          boost::lexical_cast<std::string>(Tensor_Dim1 * Tensor_Dim2) +
          " but is " + boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    std::array<double *, Tensor_Dim1 * Tensor_Dim2> ptrs;
    for (auto i = 0; i != Tensor_Dim1 * Tensor_Dim2; ++i)
      ptrs[i] = &data(i, 0);
    return FTensor::Tensor2<FTensor::PackPtr<double *, S>, Tensor_Dim1,
                            Tensor_Dim2>(ptrs);
  }
};

/**
 * \brief Get tensor rank 2 (matrix) form data matrix
 */
template <int Tensor_Dim1, int Tensor_Dim2>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, Tensor_Dim1, Tensor_Dim2>
getFTensor2FromMat(MatrixDouble &data) {
  return GetFTensor2FromMatImpl<Tensor_Dim1, Tensor_Dim2, 1, double,
                                ublas::row_major, DoubleAllocator>::get(data);
}

template <int Tensor_Dim1, int Tensor_Dim2>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, Tensor_Dim1, Tensor_Dim2>
getFTensor2FromVec(VectorDouble &data);

/**
 * Template specialization for getFTensor2FromMat
 */
template <>
inline FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 1, 1>
getFTensor2FromVec(VectorDouble &data) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 1>, 1, 1>(&*data.begin());
}

template <int Tensor_Dim, int S, class T, class L, class A>
struct GetFTensor2SymmetricFromMatImpl {};

template <int S, class T, class L, class A>
struct GetFTensor2SymmetricFromMatImpl<3, S, T, L, A> {
  static inline auto get(ublas::matrix<T, L, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 6)
      THROW_MESSAGE(
          "getFTensor2SymmetricFromMat<3>: wrong size of data matrix, numer "
          "of rows should be 6 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor2_symmetric<FTensor::PackPtr<T *, S>, 3>(
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
        &data(5, 0));
  }
};

template <int S, class T, class L, class A>
struct GetFTensor2SymmetricFromMatImpl<2, S, T, L, A> {
  static inline auto get(ublas::matrix<T, L, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 3)
      THROW_MESSAGE(
          "getFTensor2SymmetricFromMat<2>: wrong size of data matrix, number "
          "of rows should be 3 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor2_symmetric<FTensor::PackPtr<T *, S>, 2>(
        &data(0, 0), &data(1, 0), &data(2, 0));
  }
};

/**
 * \brief Get symmetric tensor rank 2 (matrix) form data matrix
 */
template <int Tensor_Dim, int S, class T, class L, class A>
static inline auto getFTensor2SymmetricFromMat(ublas::matrix<T, L, A> &data) {
  return GetFTensor2SymmetricFromMatImpl<Tensor_Dim, S, T, L, A>::get(data);
}

template <int Tensor_Dim, int S = 1>
static inline auto getFTensor2SymmetricFromMat(MatrixDouble &data) {
  return getFTensor2SymmetricFromMat<Tensor_Dim, S, double, ublas::row_major,
                                     DoubleAllocator>(data);
}

template <int Tensor_Dim01, int Tensor_Dim23, int S, class T, class L, class A>
struct GetFTensor4DdgFromMatImpl {};

template <int S, class T, class A>
struct GetFTensor4DdgFromMatImpl<1, 1, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 1)
      THROW_MESSAGE(
          "getFTensor4DdgFromMat<1, 1>: wrong size of data matrix, number "
          "of rows should be 1 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Ddg<FTensor::PackPtr<double *, S>, 1, 1>{&data(0, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor4DdgFromMatImpl<2, 2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 9) {
      THROW_MESSAGE(
          "getFTensor4DdgFromMat<2, 2>: wrong size of data matrix, number "
          "of rows should be 9 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Ddg<FTensor::PackPtr<double *, S>, 2, 2>{
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
        &data(5, 0), &data(6, 0), &data(7, 0), &data(8, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor4DdgFromMatImpl<3, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 36) {
      cerr << data.size1() << endl;
      THROW_MESSAGE(
          "getFTensor4DdgFromMat<3, 3>: wrong size of data matrix, number "
          "of rows should be 36 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Ddg<FTensor::PackPtr<double *, S>, 3, 3>{
        &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),  &data(4, 0),
        &data(5, 0),  &data(6, 0),  &data(7, 0),  &data(8, 0),  &data(9, 0),
        &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
        &data(15, 0), &data(16, 0), &data(17, 0), &data(18, 0), &data(19, 0),
        &data(20, 0), &data(21, 0), &data(22, 0), &data(23, 0), &data(24, 0),
        &data(25, 0), &data(26, 0), &data(27, 0), &data(28, 0), &data(29, 0),
        &data(30, 0), &data(31, 0), &data(32, 0), &data(33, 0), &data(34, 0),
        &data(35, 0)};
  }
};

/**
 * @brief Get symmetric tensor rank 4  on first two and last indices from
 * form data matrix
 *
 * @tparam Tensor_Dim01 dimension of first two indicies
 * @tparam Tensor_Dim23 dimension of second two indicies
 * @tparam Memory shift
 * @tparam T the type of object stored
 * @tparam L the storage organization
 * @tparam A 	the type of Storage array
 * @param data data container
 * @return FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
 */
template <int Tensor_Dim01, int Tensor_Dim23, int S = 1, class T, class L,
          class A>
static inline FTensor::Ddg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, Tensor_Dim23>
getFTensor4DdgFromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value,
                "Such getFTensor4DdgFromMat specialisation is not implemented");
}

template <int Tensor_Dim01, int Tensor_Dim23, int S = 1>
static inline auto getFTensor4DdgFromMat(MatrixDouble &data) {
  return GetFTensor4DdgFromMatImpl<Tensor_Dim01, Tensor_Dim23, S, double,
                                   ublas::row_major,
                                   DoubleAllocator>::get(data);
}

template <int Tensor_Dim01, int Tensor_Dim23, int S, class T>
struct GetFTensor4DdgFromPtrImpl;

template <int S, class T> struct GetFTensor4DdgFromPtrImpl<3, 3, S, T> {
  static inline auto get(T *ptr) {
    return FTensor::Ddg<FTensor::PackPtr<T *, S>, 3, 3>{

        ptr + 0,  ptr + 1,  ptr + 2,  ptr + 3,  ptr + 4,  ptr + 5,
        ptr + 6,  ptr + 7,  ptr + 8,  ptr + 9,  ptr + 10, ptr + 11,
        ptr + 12, ptr + 13, ptr + 14, ptr + 15, ptr + 16, ptr + 17,
        ptr + 18, ptr + 19, ptr + 20, ptr + 21, ptr + 22, ptr + 23,
        ptr + 24, ptr + 25, ptr + 26, ptr + 27, ptr + 28, ptr + 29,
        ptr + 30, ptr + 31, ptr + 32, ptr + 33, ptr + 34, ptr + 35};
  }
};

template <int Tensor_Dim01, int Tensor_Dim23, int S = 1, class T = double>
static inline auto getFTensor4DdgFromPtr(T *ptr) {
  return GetFTensor4DdgFromPtrImpl<Tensor_Dim01, Tensor_Dim23, S, T>::get(ptr);
}

template <int Tensor_Dim01, int Tensor_Dim2, int S, class T, class L, class A>
struct GetFTensor3DgFromMatImpl {};

template <int S, class T, class A>
struct GetFTensor3DgFromMatImpl<1, 1, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 1)
      THROW_MESSAGE(
          "getFTensor3DgFromMat<1, 1>: wrong size of data matrix, number "
          "of rows should be 1 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Dg<FTensor::PackPtr<double *, S>, 1, 1>{&data(0, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3DgFromMatImpl<2, 2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 6) {
      THROW_MESSAGE(
          "getFTensor4DdgFromMat<2, 2>: wrong size of data matrix, number "
          "of rows should be 6 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Dg<FTensor::PackPtr<double *, S>, 2, 2>{
        &data(0, 0), &data(1, 0), &data(2, 0),
        &data(3, 0), &data(4, 0), &data(5, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3DgFromMatImpl<3, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 18) {
      cerr << data.size1() << endl;
      THROW_MESSAGE(
          "getFTensor3DgFromMat<3, 3>: wrong size of data matrix, number "
          "of rows should be 18 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Dg<FTensor::PackPtr<double *, S>, 3, 3>{
        &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),  &data(4, 0),
        &data(5, 0),  &data(6, 0),  &data(7, 0),  &data(8, 0),  &data(9, 0),
        &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
        &data(15, 0), &data(16, 0), &data(17, 0)};
  }
};

/**
 * @brief Get symmetric tensor rank 3  on the first two indices from
 * form data matrix
 *
 * @tparam Tensor_Dim01 dimension of first two indicies
 * @tparam Tensor_Dim2 dimension of last index
 * @tparam T the type of object stored
 * @tparam L the storage organization
 * @tparam A 	the type of Storage array
 * @param data data container
 * @return FTensor::Dg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, TensorDim23>
 */
template <int Tensor_Dim01, int Tensor_Dim2, int S = 1, class T, class L,
          class A>
static inline FTensor::Dg<FTensor::PackPtr<T *, 1>, Tensor_Dim01, Tensor_Dim2>
getFTensor3DgFromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value,
                "Such getFTensor3DgFromMat specialisation is not implemented");
}

template <int Tensor_Dim01, int Tensor_Dim2, int S = 1>
static inline auto getFTensor3DgFromMat(MatrixDouble &data) {
  return GetFTensor3DgFromMatImpl<Tensor_Dim01, Tensor_Dim2, S, double,
                                  ublas::row_major, DoubleAllocator>::get(data);
}

template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int Tensor_Dim3,
          int S, class T, class L, class A>
struct GetFTensor4FromMatImpl {};

template <int S, class T, class A>
struct GetFTensor4FromMatImpl<1, 1, 1, 1, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 1)
      THROW_MESSAGE(
          "getFTensor4FromMat<1, 1, 1, 1>: wrong size of data matrix, number "
          "of rows should be 1 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor4<FTensor::PackPtr<double *, S>, 1, 1, 1, 1>{
        &data(0, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor4FromMatImpl<2, 2, 2, 2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 16) {
      THROW_MESSAGE(
          "getFTensor4FromMat<2, 2, 2, 2>: wrong size of data matrix, number "
          "of rows should be 16 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Tensor4<FTensor::PackPtr<double *, S>, 2, 2, 2, 2>{
        &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),
        &data(4, 0),  &data(5, 0),  &data(6, 0),  &data(7, 0),
        &data(8, 0),  &data(9, 0),  &data(10, 0), &data(11, 0),
        &data(12, 0), &data(13, 0), &data(14, 0), &data(15, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor4FromMatImpl<3, 3, 3, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 81) {
      cerr << data.size1() << endl;
      THROW_MESSAGE(
          "getFTensor4FromMat<3, 3, 3, 3>: wrong size of data matrix, number "
          "of rows should be 81 but is " +
          boost::lexical_cast<std::string>(data.size1()));
    }
#endif
    return FTensor::Tensor4<FTensor::PackPtr<double *, S>, 3, 3, 3, 3>{
        &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),  &data(4, 0),
        &data(5, 0),  &data(6, 0),  &data(7, 0),  &data(8, 0),  &data(9, 0),
        &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
        &data(15, 0), &data(16, 0), &data(17, 0), &data(18, 0), &data(19, 0),
        &data(20, 0), &data(21, 0), &data(22, 0), &data(23, 0), &data(24, 0),
        &data(25, 0), &data(26, 0), &data(27, 0), &data(28, 0), &data(29, 0),
        &data(30, 0), &data(31, 0), &data(32, 0), &data(33, 0), &data(34, 0),
        &data(35, 0), &data(36, 0), &data(37, 0), &data(38, 0), &data(39, 0),
        &data(40, 0), &data(41, 0), &data(42, 0), &data(43, 0), &data(44, 0),
        &data(45, 0), &data(46, 0), &data(47, 0), &data(48, 0), &data(49, 0),
        &data(50, 0), &data(51, 0), &data(52, 0), &data(53, 0), &data(54, 0),
        &data(55, 0), &data(56, 0), &data(57, 0), &data(58, 0), &data(59, 0),
        &data(60, 0), &data(61, 0), &data(62, 0), &data(63, 0), &data(64, 0),
        &data(65, 0), &data(66, 0), &data(67, 0), &data(68, 0), &data(69, 0),
        &data(70, 0), &data(71, 0), &data(72, 0), &data(73, 0), &data(74, 0),
        &data(75, 0), &data(76, 0), &data(77, 0), &data(78, 0), &data(79, 0),
        &data(80, 0)};
  }
};

/**
 * @brief Get tensor rank 4 (non symmetric) form data matrix
 *
 * @tparam Tensor_Dim0 dimension of frirst index
 * @tparam Tensor_Dim1 dimension of second index
 * @tparam Tensor_Dim2 dimension of third index
 * @tparam Tensor_Dim3 dimension of fourth index
 * @tparam T the type of object stored
 * @tparam L the storage organization
 * @tparam A 	the type of Storage array
 * @param data data container
 * @return FTensor::Tensor4<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>
 */
template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int Tensor_Dim3,
          int S = 1, class T, class L, class A>
static inline FTensor::Tensor4<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>
getFTensor4FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value,
                "Such getFTensor4FromMat specialisation is not implemented");
}

template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int Tensor_Dim3,
          int S = 1>
static inline auto getFTensor4FromMat(MatrixDouble &data) {
  return GetFTensor4FromMatImpl<Tensor_Dim0, Tensor_Dim1, Tensor_Dim2,
                                Tensor_Dim3, S, double, ublas::row_major,
                                DoubleAllocator>::get(data);
}

template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int S, class T,
          class L, class A>
struct GetFTensor3FromMatImpl {};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<1, 1, 1, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 1)
      THROW_MESSAGE(
          "getFTensor3FromMat<1, 1, 1>: wrong size of data matrix, number "
          "of rows should be 1 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 1, 1, 1>{
        &data(0, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<2, 2, 2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 8)
      THROW_MESSAGE(
          "getFTensor3FromMat<2, 2, 2>: wrong size of data matrix, number "
          "of rows should be 8 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 2, 2, 2>{
        &data(0, 0), &data(1, 0), &data(2, 0), &data(3, 0), &data(4, 0),
        &data(5, 0), &data(6, 0), &data(7, 0)

    };
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<3, 2, 2, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 12)
      THROW_MESSAGE(
          "getFTensor3FromMat<3, 2, 2>: wrong size of data matrix, number "
          "of rows should be 12 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 3, 2, 2>{
        &data(0, 0), &data(1, 0), &data(2, 0),  &data(3, 0),
        &data(4, 0), &data(5, 0), &data(6, 0),  &data(7, 0),
        &data(8, 0), &data(9, 0), &data(10, 0), &data(11, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<2, 2, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 12)
      THROW_MESSAGE(
          "getFTensor3FromMat<2, 2, 3>: wrong size of data matrix, number "
          "of rows should be 12 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 2, 2, 3>{
        &data(0, 0), &data(1, 0), &data(2, 0),  &data(3, 0),
        &data(4, 0), &data(5, 0), &data(6, 0),  &data(7, 0),
        &data(8, 0), &data(9, 0), &data(10, 0), &data(11, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<3, 3, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 27)
      THROW_MESSAGE(
          "getFTensor3FromMat<3, 3, 3>: wrong size of data matrix, number "
          "of rows should be 27 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 3, 3, 3>{
        &data(0, 0),  &data(1, 0),  &data(2, 0),  &data(3, 0),  &data(4, 0),
        &data(5, 0),  &data(6, 0),  &data(7, 0),  &data(8, 0),  &data(9, 0),
        &data(10, 0), &data(11, 0), &data(12, 0), &data(13, 0), &data(14, 0),
        &data(15, 0), &data(16, 0), &data(17, 0), &data(18, 0), &data(19, 0),
        &data(20, 0), &data(21, 0), &data(22, 0), &data(23, 0), &data(24, 0),
        &data(25, 0), &data(26, 0)};
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<6, 3, 3, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 54)
      THROW_MESSAGE(
          "getFTensor3FromMat<6, 3, 3>: wrong size of data matrix, number "
          "of rows should be 54 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    std::array<double *, 54> ptrs;
    for (auto i = 0; i != 54; ++i)
      ptrs[i] = &data(i, 0);
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 6, 3, 3>(ptrs);
  }
};

template <int S, class T, class A>
struct GetFTensor3FromMatImpl<3, 3, 6, S, T, ublas::row_major, A> {
  static inline auto get(ublas::matrix<T, ublas::row_major, A> &data) {
#ifndef NDEBUG
    if (data.size1() != 54)
      THROW_MESSAGE(
          "getFTensor3FromMat<3, 3, 6>: wrong size of data matrix, number "
          "of rows should be 54 but is " +
          boost::lexical_cast<std::string>(data.size1()));
#endif
    std::array<double *, 54> ptrs;
    for (auto i = 0; i != 54; ++i)
      ptrs[i] = &data(i, 0);
    return FTensor::Tensor3<FTensor::PackPtr<double *, S>, 3, 3, 6>(ptrs);
  }
};

/**
 * @brief Get tensor rank 3 (non symmetries) form data matrix
 *
 * @tparam Tensor_Dim0 dimension of first index
 * @tparam Tensor_Dim1 dimension of second index
 * @tparam Tensor_Dim2 dimension of third index
 * @tparam S shift size
 * @tparam T the type of object stored
 * @tparam L the storage organization
 * @tparam A 	the type of Storage array
 * @param data data container
 * @return FTensor::Tensor3<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1, Tensor_Dim2>
 */
template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int S = 1, class T,
          class L, class A>
static inline FTensor::Tensor3<FTensor::PackPtr<T *, 1>, Tensor_Dim0,
                               Tensor_Dim1, Tensor_Dim2>
getFTensor3FromMat(ublas::matrix<T, L, A> &data) {
  static_assert(!std::is_same<T, T>::value,
                "Such getFTensor3FromMat specialisation is not implemented");
}

template <int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2, int S = 1>
static inline auto getFTensor3FromMat(MatrixDouble &data) {
  return GetFTensor3FromMatImpl<Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, S,
                                double, ublas::row_major,
                                DoubleAllocator>::get(data);
}

template <int DIM, int S, typename T> struct GetFTensor1FromPtrImpl;

template <int S, typename T> struct GetFTensor1FromPtrImpl<1, S, T> {
  GetFTensor1FromPtrImpl() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 1>(ptr);
  }
};

template <int S, typename T> struct GetFTensor1FromPtrImpl<2, S, T> {
  GetFTensor1FromPtrImpl() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 2>(ptr + HVEC0,
                                                         ptr + HVEC1);
  }
};

template <int S, typename T> struct GetFTensor1FromPtrImpl<3, S, T> {
  GetFTensor1FromPtrImpl() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 3>(
        ptr + HVEC0, ptr + HVEC1, ptr + HVEC2);
  }
};

template <int S, typename T> struct GetFTensor1FromPtrImpl<4, S, T> {
  GetFTensor1FromPtrImpl() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 4>(ptr + 0, ptr + 1,
                                                         ptr + 2, ptr + 3);
  }
};

template <int S, typename T> struct GetFTensor1FromPtrImpl<6, S, T> {
  GetFTensor1FromPtrImpl() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 6>(
        ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5);
  }
};

/**
 * @brief Make Tensor1 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, 3 * DIM>, 3, DIM>
 */
template <int DIM, int S = DIM>
inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, DIM>
getFTensor1FromPtr(double *ptr) {
  return GetFTensor1FromPtrImpl<DIM, S, double>::get(ptr);
};

#ifdef WITH_ADOL_C
template <int DIM, int S = DIM>
inline FTensor::Tensor1<FTensor::PackPtr<adouble *, S>, DIM>
getFTensor1FromPtr(adouble *ptr) {
  return GetFTensor1FromPtrImpl<DIM, S, adouble>::get(ptr);
};
#endif

template <int DIM, int S = DIM>
inline FTensor::Tensor1<FTensor::PackPtr<std::complex<double> *, S>, DIM>
getFTensor1FromPtr(std::complex<double> *ptr) {
  return GetFTensor1FromPtrImpl<DIM, S, std::complex<double>>::get(ptr);
};

template <int DIM1, int DIM2, int S, typename T> struct GetFTensor2FromPtr;

template <int S, typename T> struct GetFTensor2FromPtr<3, 2, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 3, 2>(
        ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<6, 6, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 6, 6>(
        ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
        ptr + 8, ptr + 9, ptr + 10, ptr + 11, ptr + 12, ptr + 13, ptr + 14,
        ptr + 15, ptr + 16, ptr + 17, ptr + 18, ptr + 19, ptr + 20, ptr + 21,
        ptr + 22, ptr + 23, ptr + 24, ptr + 25, ptr + 26, ptr + 27, ptr + 28,
        ptr + 29, ptr + 30, ptr + 31, ptr + 32, ptr + 33, ptr + 34, ptr + 35);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<3, 3, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 3, 3>(
        ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
        ptr + 8);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<2, 2, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 2, 2>(&ptr[0], &ptr[1],
                                                            &ptr[2], &ptr[3]);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<1, 3, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 1, 3>(&ptr[0], &ptr[1],
                                                            &ptr[2]);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<1, 2, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 1, 2>(&ptr[0], &ptr[1]);
  }
};

template <int S, typename T> struct GetFTensor2FromPtr<1, 1, S, T> {
  GetFTensor2FromPtr() = delete;
  inline static auto get(T *ptr) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 1, 1>(&ptr[0]);
  }
};

/**
 * @brief Make Tensor2 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, S>, DIM1, DIM2>
 */
template <int DIM1, int DIM2, int S = DIM1 * DIM2>
inline auto getFTensor2FromPtr(double *ptr) {
  return GetFTensor2FromPtr<DIM1, DIM2, S, double>::get(ptr);
};

/**
 * @brief Make Tensor2 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, S>, DIM1, DIM2>
 */
template <int DIM1, int DIM2, int S = DIM1 * DIM2>
inline auto getFTensor2FromPtr(std::complex<double> *ptr) {
  return GetFTensor2FromPtr<DIM1, DIM2, S, std::complex<double>>::get(ptr);
};

/**
 * @brief Make Tensor2 for HVec base from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
 */
template <int DIM1, int DIM2>
inline FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
getFTensor2HVecFromPtr(double *ptr) {
  static_assert(DIM1 == DIM1 || DIM2 != DIM2,
                "Such getFTensor2HVecFromPtr is not implemented");
};

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3,
                 2> inline getFTensor2HVecFromPtr<3, 2>(double *ptr) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 6>, 3, 2>(
      ptr + HVEC0_0, ptr + HVEC0_1,

      ptr + HVEC1_0, ptr + HVEC1_1,

      ptr + HVEC2_0, ptr + HVEC2_1);
};

template <>
FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3,
                 3> inline getFTensor2HVecFromPtr<3, 3>(double *ptr) {
  return FTensor::Tensor2<FTensor::PackPtr<double *, 9>, 3, 3>(
      ptr + HVEC0_0, ptr + HVEC0_1, ptr + HVEC0_2,

      ptr + HVEC1_0, ptr + HVEC1_1, ptr + HVEC1_2,

      ptr + HVEC2_0, ptr + HVEC2_1, ptr + HVEC2_2);
};

/*
 * @brief Make Tensor3 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor3<FTensor::PackPtr<double *, DIM1 * DIM2* DIM3>, DIM1,
 * DIM2, DIM3>
 */
template <int DIM1, int DIM2, int DIM3>
inline FTensor::Tensor3<FTensor::PackPtr<double *, DIM1 * DIM2 * DIM3>, DIM1,
                        DIM2, DIM3>
getFTensor3FromPtr(double *ptr) {
  static_assert(DIM1 == DIM1 || DIM2 != DIM2 || DIM3 != DIM3,
                "Such getFTensor3FromPtr is not implemented");
};

template <>
inline FTensor::Tensor3<FTensor::PackPtr<double *, 12>, 3, 2, 2>
getFTensor3FromPtr<3, 2, 2>(double *ptr) {
  return FTensor::Tensor3<FTensor::PackPtr<double *, 12>, 3, 2, 2>(
      ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
      ptr + 8, ptr + 9, ptr + 10, ptr + 11);
};

template <>
inline FTensor::Tensor3<FTensor::PackPtr<double *, 27>, 3, 3, 3>
getFTensor3FromPtr<3, 3, 3>(double *ptr) {
  return FTensor::Tensor3<FTensor::PackPtr<double *, 27>, 3, 3, 3>(
      ptr + 0, ptr + 1, ptr + 2, ptr + 3, ptr + 4, ptr + 5, ptr + 6, ptr + 7,
      ptr + 8, ptr + 9, ptr + 10, ptr + 11, ptr + 12, ptr + 13, ptr + 14,
      ptr + 15, ptr + 16, ptr + 17, ptr + 18, ptr + 19, ptr + 20, ptr + 21,
      ptr + 22, ptr + 23, ptr + 24, ptr + 25, ptr + 26);
};

/**
 * @brief Make symmetric Tensor2 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
 */
template <int DIM>
inline FTensor::Tensor2_symmetric<
    FTensor::PackPtr<double *, (DIM * (DIM + 1)) / 2>, DIM>
getFTensor2SymmetricFromPtr(double *ptr) {
  static_assert(DIM, "Such getFTensor2SymmetricFromPtr is not implemented");
}

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>
getFTensor2SymmetricFromPtr<3>(double *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 6>, 3>(
      ptr + 0, ptr + 1, ptr + 2,

      ptr + 3, ptr + 4,

      ptr + 5);
};

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>
getFTensor2SymmetricFromPtr<2>(double *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 3>, 2>(
      &ptr[0], &ptr[1], &ptr[2]);
};

#ifdef WITH_ADOL_C

/**
 * @brief Make symmetric Tensor2 from pointer
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
 */
template <int DIM>
inline FTensor::Tensor2_symmetric<
    FTensor::PackPtr<adouble *, (DIM * (DIM + 1)) / 2>, DIM>
getFTensor2SymmetricFromPtr(adouble *ptr) {
  static_assert(DIM, "Such getFTensor2SymmetricFromPtr is not implemented");
}

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<adouble *, 6>, 3>
getFTensor2SymmetricFromPtr<3>(adouble *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<adouble *, 6>, 3>(
      ptr + 0, ptr + 1, ptr + 2,

      ptr + 3, ptr + 4,

      ptr + 5);
};

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<adouble *, 3>, 2>
getFTensor2SymmetricFromPtr<2>(adouble *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<adouble *, 3>, 2>(
      ptr + 0, ptr + 1, ptr + 2);
};

#endif

/**
 * @brief Make symmetric Tensor2 from pointer, taking lower triangle of matrix
 *
 * @tparam DIM
 * @param ptr
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, DIM1 * DIM2>, DIM1, DIM2>
 */
template <int DIM>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, DIM * DIM>, DIM>
getFTensor2SymmetricLowerFromPtr(double *ptr) {
  static_assert(DIM,
                "Such getFTensor2SymmetricLowerFromPtr is not implemented");
}

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 9>, 3>
getFTensor2SymmetricLowerFromPtr<3>(double *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 9>, 3>(
      ptr + HVEC0_0, ptr + HVEC0_1, ptr + HVEC0_2,

      ptr + HVEC1_0, ptr + HVEC1_1,

      ptr + HVEC2_2);
};

template <>
inline FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 4>, 2>
getFTensor2SymmetricLowerFromPtr<2>(double *ptr) {
  return FTensor::Tensor2_symmetric<FTensor::PackPtr<double *, 4>, 2>(
      ptr + 0, ptr + 1, ptr + 3);
};

template <int DIM, int S> struct GetFTensor1FromArray;

template <int S> struct GetFTensor1FromArray<1, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 1>{&data[0]};
  }
};

template <int S> struct GetFTensor1FromArray<2, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 2>{&data[0], &data[1]};
  }
};

template <int S> struct GetFTensor1FromArray<3, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 3>{&data[0], &data[1],
                                                         &data[2]};
  }
};

template <int S> struct GetFTensor1FromArray<4, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 4>{&data[0], &data[1],
                                                         &data[2], &data[3]};
  }
};

template <int S> struct GetFTensor1FromArray<6, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 6>{
        &data[0], &data[1], &data[2], &data[3], &data[4], &data[5]};
  }
};

template <int S> struct GetFTensor1FromArray<9, S> {
  GetFTensor1FromArray() = delete;
  template <typename V> static inline auto get(V &data) {
    using T = typename std::remove_reference<decltype(data[0])>::type;
    return FTensor::Tensor1<FTensor::PackPtr<T *, S>, 9>{
        &data[0], &data[1], &data[2], &data[3], &data[4],
        &data[5], &data[6], &data[7], &data[8]};
  }
};

/**
 * @brief Get FTensor1 from array
 *
 * \todo Generalise for different arrays and data types
 *
 * @tparam DIM
 * @param data
 * @return FTensor::Tensor1<FTensor::PackPtr<double *, S>, DIM>
 */
template <int DIM, int S> inline auto getFTensor1FromArray(VectorDouble &data) {
  return GetFTensor1FromArray<DIM, S>::get(data);
}

/** @copydoc getFTensor1FromArray */
template <int DIM, int S = 0>
inline auto getFTensor1FromArray(VectorDouble3 &data);

template <> inline auto getFTensor1FromArray<3, 0>(VectorDouble3 &data) {
  return GetFTensor1FromArray<3, 0>::get(data);
}

template <int DIM, int S>
inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, DIM>
getFTensor1FromMat(MatrixDouble &data, const size_t rr);

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>
getFTensor1FromMat(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 2>{&data(rr + 0, 0),
                                                            &data(rr + 1, 0)};
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>
getFTensor1FromMat(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 3>{
      &data(rr + 0, 0), &data(rr + 1, 0), &data(rr + 2, 0)};
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 4>
getFTensor1FromMat(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 4>{
      &data(rr + 0, 0), &data(rr + 1, 0), &data(rr + 2, 0), &data(rr + 3, 0)};
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 9>
getFTensor1FromMat(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 1>, 9>{
      &data(rr + 0, 0), &data(rr + 1, 0), &data(rr + 2, 0),
      &data(rr + 3, 0), &data(rr + 4, 0), &data(rr + 5, 0),
      &data(rr + 6, 0), &data(rr + 7, 0), &data(rr + 8, 0)};
}

/**
 * @brief Get FTensor1 from array
 *
 * \todo Generalise for diffrent arrays and data types
 *
 * @tparam DIM
 * @param data
 * @param rr
 * @return FTensor::Tensor1<FTensor::PackPtr<double *, DIM>, DIM>
 */
template <int DIM, int S>
inline FTensor::Tensor1<FTensor::PackPtr<double *, S>, DIM>
getFTensor1FromArrayDiag(MatrixDouble &data, const size_t rr) {
  static_assert(DIM != DIM, "not implemented");
  return FTensor::Tensor1<FTensor::PackPtr<double *, S>, DIM>();
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>
getFTensor1FromArrayDiag(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 2>, 2>{&data(rr + 0, 0),
                                                            &data(rr + 1, 1)};
}

template <>
inline FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>
getFTensor1FromArrayDiag(MatrixDouble &data, const size_t rr) {
  return FTensor::Tensor1<FTensor::PackPtr<double *, 3>, 3>{
      &data(rr + 0, 0), &data(rr + 1, 1), &data(rr + 2, 2)};
}

/**
 * @brief Get FTensor2 from array
 *
 * \note Generalise for other data types
 *
 * @tparam DIM1
 * @tparam DIM2
 * @tparam S
 * @param data
 * @return FTensor::Tensor2<FTensor::PackPtr<double *, S>, DIM1, DIM2>
 */
template <int DIM1, int DIM2, int S, class T, class L, class A>
struct GetFTensor2FromArrayImpl;

template <int DIM1, int DIM2, class T, class L, class A>
struct GetFTensor2FromArrayRawPtrImpl;

template <int S, class T, class L, class A>
struct GetFTensor2FromArrayImpl<2, 2, S, T, L, A> {
  GetFTensor2FromArrayImpl() = delete;
  inline static auto get(ublas::matrix<T, L, A> &data, const size_t rr,
                         const size_t cc) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 2, 2>{
        &data(rr + 0, cc + 0), &data(rr + 0, cc + 1),

        &data(rr + 1, cc + 0), &data(rr + 1, cc + 1)};
  }
};

template <int S, class T, class L, class A>
struct GetFTensor2FromArrayImpl<3, 3, S, T, L, A> {
  GetFTensor2FromArrayImpl() = delete;
  inline static auto get(ublas::matrix<T, L, A> &data, const size_t rr,
                         const size_t cc) {
    return FTensor::Tensor2<FTensor::PackPtr<T *, S>, 3, 3>{
        &data(rr + 0, cc + 0), &data(rr + 0, cc + 1), &data(rr + 0, cc + 2),
        &data(rr + 1, cc + 0), &data(rr + 1, cc + 1), &data(rr + 1, cc + 2),
        &data(rr + 2, cc + 0), &data(rr + 2, cc + 1), &data(rr + 2, cc + 2)};
  }
};

template <class T, class L, class A>
struct GetFTensor2FromArrayRawPtrImpl<2, 2, T, L, A> {
  GetFTensor2FromArrayRawPtrImpl() = delete;
  inline static auto get(ublas::matrix<T, L, A> &data, const size_t rr,
                         const size_t cc, const int ss = 0) {
    return FTensor::Tensor2<T *, 2, 2>(
        &data(rr + 0, cc + 0), &data(rr + 0, cc + 1),

        &data(rr + 1, cc + 0), &data(rr + 1, cc + 1), ss);
  }
};

template <class T, class L, class A>
struct GetFTensor2FromArrayRawPtrImpl<3, 3, T, L, A> {
  GetFTensor2FromArrayRawPtrImpl() = delete;
  inline static auto get(ublas::matrix<T, L, A> &data, const size_t rr,
                         const size_t cc, const int ss = 0) {
    return FTensor::Tensor2<T *, 3, 3>(
        &data(rr + 0, cc + 0), &data(rr + 0, cc + 1), &data(rr + 0, cc + 2),
        &data(rr + 1, cc + 0), &data(rr + 1, cc + 1), &data(rr + 1, cc + 2),
        &data(rr + 2, cc + 0), &data(rr + 2, cc + 1), &data(rr + 2, cc + 2),
        ss);
  }
};

template <int DIM1, int DIM2, int S>
inline FTensor::Tensor2<FTensor::PackPtr<double *, S>, DIM1, DIM2>
getFTensor2FromArray(MatrixDouble &data, const size_t rr, const size_t cc = 0) {
  return GetFTensor2FromArrayImpl<DIM1, DIM2, S, double, ublas::row_major,
                                  VecAllocator<double>>::get(data, rr, cc);
}

template <int DIM1, int DIM2>
inline FTensor::Tensor2<double *, DIM1, DIM2>
getFTensor2FromArray(MatrixDouble &data, const size_t rr, const size_t cc,
                     const int ss) {
  return GetFTensor2FromArrayRawPtrImpl<DIM1, DIM2, double, ublas::row_major,
                                        VecAllocator<double>>::get(data, rr, cc,
                                                                   ss);
}

template <int S, typename T, typename L, typename A>
inline auto getFTensor2FromArray2by2(ublas::matrix<T, L, A> &data,
                                     const FTensor::Number<S> &,
                                     const size_t rr, const size_t cc = 0) {
  return GetFTensor2FromArrayImpl<2, 2, S, T, L, A>::get(data, rr, cc);
}

template <int S, typename T, typename L, typename A>
inline auto getFTensor2FromArray3by3(ublas::matrix<T, L, A> &data,
                                     const FTensor::Number<S> &,
                                     const size_t rr, const size_t cc = 0) {
  return GetFTensor2FromArrayImpl<3, 3, S, T, L, A>::get(data, rr, cc);
}

#ifdef WITH_ADOL_C

template <int DIM1, int DIM2, int S>
inline auto getFTensor2FromArray(MatrixADouble &data, const size_t rr) {
  return GetFTensor2FromArrayImpl<DIM1, DIM2, S, adouble, ublas::row_major,
                                  VecAllocator<adouble>>::get(data, rr);
}

#endif

// list of lapack wrappers
/**
 * @brief compute matrix inverse with lapack dgetri
 *
 * @param mat input square matrix / output inverse matrix
 * @return MoFEMErrorCode
 */
inline MoFEMErrorCode computeMatrixInverse(MatrixDouble &mat) {
  MoFEMFunctionBegin;

  const size_t M = mat.size1();
  const size_t N = mat.size2();

  if (M != N)
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "The input matrix for inverse computation is not square %d != %d",
             M, N);

  int *ipv = new int[N];
  int lwork = N * N;
  double *work = new double[lwork];
  int info;
  info = lapack_dgetrf(N, N, &*mat.data().begin(), N, ipv);
  if (info != 0)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "lapack error info = %d", info);
  info = lapack_dgetri(N, &*mat.data().begin(), N, ipv, work, lwork);
  if (info != 0)
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "lapack error info = %d", info);

  delete[] ipv;
  delete[] work;

  MoFEMFunctionReturn(0);
}
/**
 * @brief solve linear system with lapack dgesv
 *
 * @param mat input lhs square matrix / output L and U from the factorization
 * @param f input rhs vector / output solution vector
 * @return MoFEMErrorCode
 */
inline MoFEMErrorCode solveLinearSystem(MatrixDouble &mat, VectorDouble &f) {
  MoFEMFunctionBegin;

  const size_t M = mat.size1();
  const size_t N = mat.size2();

  if (M == 0 || M != N)
    SETERRQ2(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
             "The input matrix for inverse computation is not square %d != %d",
             M, N);
  if (f.size() != M)
    f.resize(M, false);

  const int nrhs = 1;
  int info;
  int *ipiv = new int[M];
  info = lapack_dgesv(M, nrhs, &*mat.data().begin(), M, ipiv,
                      &*f.data().begin(), M);
  if (info != 0) {
    SETERRQ1(PETSC_COMM_SELF, MOFEM_OPERATION_UNSUCCESSFUL,
             "error lapack solve dgesv info = %d", info);
  }

  delete[] ipiv;
  MoFEMFunctionReturn(0);
}

/**
 * @brief Solve linear system of equations using Lapack
 *
 * @param mat
 * @param f
 * @return MoFEMErrorCode
 */
inline MoFEMErrorCode solveLinearSystem(const MatrixDouble &mat,
                                        VectorDouble &f) {
  MoFEMFunctionBegin;
  // copy matrix since on output lapack returns factorisation
  auto mat_copy = mat;
  CHKERR solveLinearSystem(mat_copy, f);
  MoFEMFunctionReturn(0);
}

/**
 * @brief compute eigenvalues of a symmetric matrix using lapack dsyev
 *
 * @param mat input symmetric matrix
 * @param eig output eigen values sorted
 * @param eigen_vec output matrix of row eigen vectors
 * @return MoFEMErrorCode
 */
inline MoFEMErrorCode computeEigenValuesSymmetric(const MatrixDouble &mat,
                                                  VectorDouble &eig,
                                                  MatrixDouble &eigen_vec) {
  MoFEMFunctionBegin;

  const size_t M = mat.size1();
  const size_t N = mat.size2();

  if (M == 0 || M != N)
    SETERRQ2(
        PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
        "The input matrix for eigen value computation is not square %d != %d",
        M, N);
  if (eig.size() != M)
    eig.resize(M, false);

  eigen_vec = mat;
  const int n = M;
  const int lda = M;
  const int size = (M + 2) * M;
  int lwork = size;
  double *work = new double[size];

  if (lapack_dsyev('V', 'U', n, &*eigen_vec.data().begin(), lda,
                   &*eig.data().begin(), work, lwork) > 0)
    PetscPrintf(PETSC_COMM_SELF,
            "The algorithm failed to compute eigenvalues.");

  delete[] work;
  MoFEMFunctionReturn(0);
}
/**
 * @brief compute eigenvalues of a symmetric matrix using lapack dsyev
 *
 * @tparam DIM
 * @param eigen_vec input / output DIM x DIM matrix of row eigen vectors
 * @param eig output eigen values sorted
 * @return MoFEMErrorCode
 */
template <int DIM, typename T1, typename T2>
inline MoFEMErrorCode
computeEigenValuesSymmetric(FTensor::Tensor2<T1, DIM, DIM> &eigen_vec,
                            FTensor::Tensor1<T2, DIM> &eig) {
  MoFEMFunctionBegin;

  const int n = DIM;
  const int lda = DIM;
  const int lwork = (DIM + 2) * DIM;
  std::array<double, (DIM + 2) * DIM> work;

  if (lapack_dsyev('V', 'U', n, &eigen_vec(0, 0), lda, &eig(0), work.data(),
                   lwork) > 0)
    PetscPrintf(PETSC_COMM_SELF,
            "The algorithm failed to compute eigenvalues.");
  MoFEMFunctionReturn(0);
}

/**
 * @brief compute eigenvalues of a symmetric tensor using lapack dsyev
 *
 * @tparam DIM
 * @param mat input tensor pointer of size DIM x DIM
 * @param eig output eigen values sorted
 * @param eigen_vec output matrix of row eigen vectors
 * @return MoFEMErrorCode
 */
template <int DIM, typename T1, typename T2, typename T3>
inline MoFEMErrorCode
computeEigenValuesSymmetric(const FTensor::Tensor2_symmetric<T1, DIM> &mat,
                            FTensor::Tensor1<T2, DIM> &eig,
                            FTensor::Tensor2<T3, DIM, DIM> &eigen_vec) {
  MoFEMFunctionBegin;
  for (int ii = 0; ii != DIM; ii++)
    for (int jj = 0; jj != DIM; jj++)
      eigen_vec(ii, jj) = mat(ii, jj);

  CHKERR computeEigenValuesSymmetric(eigen_vec, eig);

  MoFEMFunctionReturn(0);
}

/**
 * @brief Calculate the determinant of a 3x3 matrix or a tensor of rank 2
 *
 * @tparam T
 * @param t
 * @return double
 */
template <typename T> static inline auto determinantTensor3by3(T &t) {
  return t(0, 0) * t(1, 1) * t(2, 2) + t(1, 0) * t(2, 1) * t(0, 2) +
         t(2, 0) * t(0, 1) * t(1, 2) - t(0, 0) * t(2, 1) * t(1, 2) -
         t(2, 0) * t(1, 1) * t(0, 2) - t(1, 0) * t(0, 1) * t(2, 2);
}

/**
 * @brief Calculate the determinant of a 2x2 matrix or a tensor of rank 2
 *
 * @tparam T
 * @param t
 * @return double
 */
template <typename T> static inline auto determinantTensor2by2(T &t) {
  return t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0);
}

template <typename T, int DIM> struct DeterminantTensorImpl;

template <typename T> struct DeterminantTensorImpl<T, 3> {
  static inline auto get(T &t) { return determinantTensor3by3(t); }
};

template <typename T> struct DeterminantTensorImpl<T, 2> {
  static auto get(T &t) { return determinantTensor2by2(t); }
};

/**
 * @brief Calculate the determinant of a tensor of rank DIM
 */
template <typename T, int DIM>
static inline auto determinantTensor(FTensor::Tensor2<T, DIM, DIM> &t) {
  return DeterminantTensorImpl<FTensor::Tensor2<T, DIM, DIM>, DIM>::get(t);
}

/**
 * @brief Calculate the determinant of a tensor of rank DIM
 */
template <typename T, int DIM>
static inline auto determinantTensor(FTensor::Tensor2_symmetric<T, DIM> &t) {
  return DeterminantTensorImpl<FTensor::Tensor2_symmetric<T, DIM>, DIM>::get(t);
}

/**
 * \brief Calculate inverse of tensor rank 2 at integration points

 */
template <int Tensor_Dim, class T, class L, class A>
inline MoFEMErrorCode invertTensor3by3(ublas::matrix<T, L, A> &jac_data,
                                       ublas::vector<T, A> &det_data,
                                       ublas::matrix<T, L, A> &inv_jac_data) {
  MoFEMFunctionBeginHot;
  SETERRQ(PETSC_COMM_SELF, MOFEM_NOT_IMPLEMENTED,
          "Specialization for this template not yet implemented");
  MoFEMFunctionReturnHot(0);
}

template <>
inline MoFEMErrorCode
invertTensor3by3<3, double, ublas::row_major, DoubleAllocator>(
    MatrixDouble &jac_data, VectorDouble &det_data, MatrixDouble &inv_jac_data);

/**
 * \brief Calculate determinant 3 by 3

 */
template <class T1, class T2>
inline MoFEMErrorCode determinantTensor3by3(T1 &t, T2 &det) {
  MoFEMFunctionBeginHot;
  det = determinantTensor3by3(t);
  MoFEMFunctionReturnHot(0);
}

/**
 * \brief Calculate determinant 2 by 2

 */
template <class T1, class T2>
inline MoFEMErrorCode determinantTensor2by2(T1 &t, T2 &det) {
  MoFEMFunctionBeginHot;
  det = determinantTensor2by2(t);
  MoFEMFunctionReturnHot(0);
}

/**
 * \brief Calculate matrix inverse 3 by 3

 */
template <class T1, class T2, class T3>
inline MoFEMErrorCode invertTensor3by3(T1 &t, T2 &det, T3 &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 0) = (t(1, 2) * t(2, 0) - t(1, 0) * t(2, 2)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 0) = (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0)) * inv_det;
  inv_t(2, 1) = (t(0, 1) * t(2, 0) - t(0, 0) * t(2, 1)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

/**
 * \brief Calculate matrix inverse 2 by 2

 */
template <class T1, class T2, class T3>
inline MoFEMErrorCode invertTensor2by2(T1 &t, T2 &det, T3 &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = t(1, 1) * inv_det;
  inv_t(0, 1) = -t(0, 1) * inv_det;
  inv_t(1, 0) = -t(1, 0) * inv_det;
  inv_t(1, 1) = t(0, 0) * inv_det;
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_ADOL_C

/**
 * \brief Calculate matrix inverse, specialization for adouble tensor

 */
template <>
inline MoFEMErrorCode invertTensor3by3<FTensor::Tensor2<adouble, 3, 3>, adouble,
                                       FTensor::Tensor2<adouble, 3, 3>>(
    FTensor::Tensor2<adouble, 3, 3> &t, adouble &det,
    FTensor::Tensor2<adouble, 3, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) / det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) / det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) / det;
  inv_t(1, 0) = (t(1, 2) * t(2, 0) - t(1, 0) * t(2, 2)) / det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) / det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) / det;
  inv_t(2, 0) = (t(1, 0) * t(2, 1) - t(1, 1) * t(2, 0)) / det;
  inv_t(2, 1) = (t(0, 1) * t(2, 0) - t(0, 0) * t(2, 1)) / det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) / det;
  MoFEMFunctionReturnHot(0);
}

#endif

/**
 * \brief Calculate matrix inverse, specialization for symmetric tensor

 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<double, 3>, double,
                 FTensor::Tensor2_symmetric<double, 3>>(
    FTensor::Tensor2_symmetric<double, 3> &t, double &det,
    FTensor::Tensor2_symmetric<double, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

#ifdef WITH_ADOL_C

/**
 * \brief Calculate matrix inverse, specialization for adouble symmetric tensor

 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<adouble, 3>, adouble,
                 FTensor::Tensor2_symmetric<adouble, 3>>(
    FTensor::Tensor2_symmetric<adouble, 3> &t, adouble &det,
    FTensor::Tensor2_symmetric<adouble, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) / det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) / det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) / det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) / det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) / det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) / det;
  MoFEMFunctionReturnHot(0);
}

#endif

/**
 * \brief Calculate matrix inverse, specialization for symmetric (pointer)
 tensor

 */
template <>
inline MoFEMErrorCode
invertTensor3by3<FTensor::Tensor2_symmetric<double, 3>, double,
                 FTensor::Tensor2_symmetric<double *, 3>>(
    FTensor::Tensor2_symmetric<double, 3> &t, double &det,
    FTensor::Tensor2_symmetric<double *, 3> &inv_t) {
  MoFEMFunctionBeginHot;
  const auto inv_det = 1. / det;
  inv_t(0, 0) = (t(1, 1) * t(2, 2) - t(1, 2) * t(2, 1)) * inv_det;
  inv_t(0, 1) = (t(0, 2) * t(2, 1) - t(0, 1) * t(2, 2)) * inv_det;
  inv_t(0, 2) = (t(0, 1) * t(1, 2) - t(0, 2) * t(1, 1)) * inv_det;
  inv_t(1, 1) = (t(0, 0) * t(2, 2) - t(0, 2) * t(2, 0)) * inv_det;
  inv_t(1, 2) = (t(0, 2) * t(1, 0) - t(0, 0) * t(1, 2)) * inv_det;
  inv_t(2, 2) = (t(0, 0) * t(1, 1) - t(0, 1) * t(1, 0)) * inv_det;
  MoFEMFunctionReturnHot(0);
}

template <typename T1, typename T2, typename T3, int DIM>
struct InvertTensorImpl;

template <typename T1, typename T2, typename T3>
struct InvertTensorImpl<T1, T2, T3, 3> {
  inline static MoFEMErrorCode invert(T1 &t, T2 &det, T3 &inv_t) {
    return invertTensor3by3(t, det, inv_t);
  }
};

template <typename T1, typename T2, typename T3>
struct InvertTensorImpl<T1, T2, T3, 2> {
  inline static MoFEMErrorCode invert(T1 &t, T2 &det, T3 &inv_t) {
    return invertTensor2by2(t, det, inv_t);
  }
};

template <typename T1, typename T2, typename T3, int DIM>
static inline MoFEMErrorCode
invertTensor(FTensor::Tensor2<T1, DIM, DIM> &t, T2 &det,
             FTensor::Tensor2<T3, DIM, DIM> &inv_t) {
  return InvertTensorImpl<FTensor::Tensor2<T1, DIM, DIM>, T2,
                          FTensor::Tensor2<T3, DIM, DIM>, DIM>::invert(t, det,
                                                                       inv_t);
}

template <typename T1, typename T2, typename T3, int DIM>
static inline MoFEMErrorCode
invertTensor(FTensor::Tensor2_symmetric<T1, DIM> &t, T2 &det,
             FTensor::Tensor2_symmetric<T3, DIM> &inv_t) {
  return InvertTensorImpl<FTensor::Tensor2_symmetric<T1, DIM>, T2,
                          FTensor::Tensor2_symmetric<T3, DIM>,
                          DIM>::invert(t, det, inv_t);
}

/**
 * @brief Extract entity handle form multi-index container
 *
 */
struct RefEntExtractor {
  template <typename Iterator>
  static inline EntityHandle extract(const Iterator &it) {
    return (*it)->getEnt();
  }
};

/**
 * @brief Insert ordered mofem multi-index into range
 *
 * \note Inserted range has to be ordered.
 *
 * \code
 * auto hi_rit = refEntsPtr->upper_bound(start);
 * auto hi_rit = refEntsPtr->upper_bound(end);
 * Range to_erase;
 * insertOrdered(to_erase, RefEntExtractor(), rit, hi_rit);
 * \endcode
 *
 * @tparam Iterator
 * @param r
 * @param begin_iter
 * @param end_iter
 * @return moab::Range::iterator
 */
template <typename Extractor, typename Iterator>
moab::Range::iterator insertOrdered(Range &r, Extractor, Iterator begin_iter,
                                    Iterator end_iter) {
  moab::Range::iterator hint = r.begin();
  while (begin_iter != end_iter) {
    size_t j = 0;
    auto bi = Extractor::extract(begin_iter);
    Iterator pj = begin_iter;
    while (pj != end_iter && (bi + j) == Extractor::extract(pj)) {
      ++pj;
      ++j;
    }
    hint = r.insert(hint, bi, bi + (j - 1));
    begin_iter = pj;
  }
  return hint;
};

/**
 * @brief Do nothing, used to rebuild database
 *
 */
struct Modify_change_nothing {
  Modify_change_nothing() = default;
  template <typename T> inline void operator()(T &e) {}
};

/**
 * @brief Template used to reconstruct multi-index
 *
 * @tparam MI multi-index
 * @tparam Modifier
 * @param mi
 * @param mo
 * @return MoFEMErrorCode
 */
template <typename MI, typename MO = Modify_change_nothing>
inline MoFEMErrorCode reconstructMultiIndex(const MI &mi,
                                            MO &&mo = Modify_change_nothing()) {
  MoFEMFunctionBegin;
  for (auto it = mi.begin(); it != mi.end(); ++it) {
    if (!const_cast<MI &>(mi).modify(it, mo))
      SETERRQ(PETSC_COMM_SELF, MOFEM_DATA_INCONSISTENCY,
              "Houston we have a problem");
  }
  MoFEMFunctionReturn(0);
}

struct TempMeshset {
  TempMeshset(moab::Interface &moab) : moab(moab) {
    rval = moab.create_meshset(MESHSET_SET, meshset);
    MOAB_THROW(rval);
  }
  virtual ~TempMeshset() { delete_meshset(); }
  operator EntityHandle() const { return meshset; }
  auto get_ptr() { return &meshset; }

private:
  void delete_meshset() {
    rval = moab.delete_entities(&meshset, 1);
    MOAB_THROW(rval);
  }
  EntityHandle meshset;
  moab::Interface &moab;
};

/**
 * @brief  Create smart pointer to temporary meshset
 *
 */
inline auto get_temp_meshset_ptr(moab::Interface &moab) {
  return boost::make_shared<TempMeshset>(moab);
};

inline auto id_from_handle(const EntityHandle h) {
  return static_cast<EntityID>(h & MB_ID_MASK);
};

/**
 * @brief get type from entity handle
 *
 */
inline auto type_from_handle(const EntityHandle h) {
  return static_cast<EntityType>(h >> MB_ID_WIDTH);
};

/**
 * @brief get entity handle from type and id
 *
 */
inline auto ent_form_type_and_id(const EntityType type, const EntityID id) {
  return (static_cast<EntityHandle>(type) << MB_ID_WIDTH) | id;
};

/**
 * @brief get entity dimension form handle
 *
 */
inline auto dimension_from_handle(const EntityHandle h) {
  return moab::CN::Dimension(type_from_handle(h));
};

/**
 * @brief get entity type name from handle
 *
 */
inline auto type_name_from_handle(const EntityHandle h) {
  return moab::CN::EntityTypeName(type_from_handle(h));
};

/**
 * @brief get field bit id from bit number
 *
 */
inline auto field_bit_from_bit_number(const int bit_number) {
  return BitFieldId().set(bit_number - 1);
};

/**
 * @brief Insert ranges
 *
 * @tparam I
 * @param f
 * @param s
 * @param tester
 * @param inserter
 * @return auto
 */
template <typename I>
auto rangeInserter(const I f, const I s, boost::function<bool(I it)> tester,
                   boost::function<MoFEMErrorCode(I f, I s)> inserter) {
  MoFEMFunctionBegin;

  auto first = f;
  while (first != s)
    if (tester(first)) {

      auto second = first;
      ++second;

      while (second != s) {
        if (tester(second))
          ++second;
        else
          break;
      }

      CHKERR inserter(first, second);

      first = second;
      if (first != s)
        ++first;

    } else {
      ++first;
    }

  MoFEMFunctionReturn(0);
}

/**
 * @brief Create Array
 *
 * See:
 * <a
 * href="https://stackoverflow.com/questions/50942556/current-status-of-stdmake-array">See
 * stack overflow</a>
 *
 * @tparam Dest
 * @tparam Arg
 * @param arg
 * @return constexpr auto
 */
template <typename Dest = void, typename... Arg>
constexpr auto make_array(Arg &&...arg) {
  if constexpr (std::is_same<void, Dest>::value)
    return std::array<std::common_type_t<std::decay_t<Arg>...>, sizeof...(Arg)>{
        {std::forward<Arg>(arg)...}};
  else
    return std::array<Dest, sizeof...(Arg)>{{std::forward<Arg>(arg)...}};
}

template <int DIM> using i_FTIndex = FTensor::Index<'i', DIM>;
template <int DIM> using j_FTIndex = FTensor::Index<'j', DIM>;
template <int DIM> using k_FTIndex = FTensor::Index<'k', DIM>;
template <int DIM> using l_FTIndex = FTensor::Index<'l', DIM>;
template <int DIM> using m_FTIndex = FTensor::Index<'m', DIM>;
template <int DIM> using n_FTIndex = FTensor::Index<'n', DIM>;
template <int DIM> using o_FTIndex = FTensor::Index<'o', DIM>;
template <int DIM> using p_FTIndex = FTensor::Index<'p', DIM>;
template <int DIM> using I_FTIndex = FTensor::Index<'I', DIM>;
template <int DIM> using J_FTIndex = FTensor::Index<'J', DIM>;
template <int DIM> using K_FTIndex = FTensor::Index<'K', DIM>;
template <int DIM> using L_FTIndex = FTensor::Index<'L', DIM>;
template <int DIM> using M_FTIndex = FTensor::Index<'M', DIM>;
template <int DIM> using N_FTIndex = FTensor::Index<'N', DIM>;

#define FTENSOR_INDEX(DIM, I) I##_FTIndex<DIM> I;

} // namespace MoFEM

#endif //__TEMPLATES_HPP__