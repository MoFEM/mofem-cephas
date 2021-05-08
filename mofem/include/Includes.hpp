/** \file LibsIncludes
 * \brief Includes of header files form Blas/Lapack, Petsc, MOAB, Boost
 */

/* MoFEM is free software: you can redistribute it and/or modify it under
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

#ifndef __INCLUDES_HPP__
#define __INCLUDES_HPP__

//TETHEN
#ifdef WITH_TETGEN
  #ifdef PI
    #undef PI
  #endif //PI
  #include <tetgen.h>
  #undef REAL
#endif //WITH_TETGEN

//STD
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdarg.h>

//PETSc
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>
#include <petscdm.h>
#include <petscao.h>
#include <petscis.h>
#if PETSC_VERSION_GE(3, 14 , 0)
#include <petscsection.h>
#endif

//STL
#include <string>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <float.h>
#include <limits.h>
#include <bitset>
#include <exception>
#include <complex>
#include <cmath>
#include <regex>

//BOOST
#define BOOST_LOG_DYN_LINK
#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/log/sources/severity_channel_logger.hpp>
#include <boost/log/attributes.hpp>
#include <boost/log/attributes/scoped_attribute.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sources/severity_feature.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/core/null_deleter.hpp>
#include <boost/log/utility/setup/file.hpp>

#define BOOST_DISABLE_THREADS

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/global_fun.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/utility/string_ref.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>
  
#include <boost/move/unique_ptr.hpp>
#include <boost/move/make_unique.hpp>
#include <boost/move/move.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/format.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include <boost/typeof/typeof.hpp>
#include <boost/type_index.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/ref.hpp>

// Metaprogramming
// #include <boost/hana.hpp>
// #include <boost/hana/for_each.hpp>

//MOAB
#include <moab/MOABConfig.h>
#include <moab/ParallelComm.hpp>
#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/Range.hpp>
#include <MBTagConventions.hpp>
#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <moab/OrientedBoxTreeTool.hpp>
#include <moab/BVHTree.hpp>
#include <moab/MeshTopoUtil.hpp>
#include <moab/ReadUtilIface.hpp>
#include <moab/GeomUtil.hpp>
#include <moab/Util.hpp>

//LAPACK
#ifdef __cplusplus
extern "C" {
#endif
  #include <cblas.h>
  #include <lapack_wrap.h>
  #include <gm_rule.h>
#ifdef __cplusplus
}
#endif

//Name spaces
using namespace moab;
using namespace std;
using boost::multi_index_container;
using namespace boost::multi_index;
using namespace boost::multiprecision;
using namespace boost::numeric;

#endif //__INCLUDES_HPP__
