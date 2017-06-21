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
  #include <tetgen.h>
  #undef REAL
#endif //WITH_TETGEN

//STD
#include <fstream>
#include <iostream>

//PETSc
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscsnes.h>
#include <petscts.h>
#include <petscdm.h>
#include <petscao.h>
#include <petscis.h>

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

//BOOST
#include <boost/multi_index_container.hpp>
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

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/pointer_cast.hpp>

// #include <boost/move/unique_ptr.hpp>
// #include <boost/move/make_unique.hpp>
// #include <boost/move/move.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/format.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

// #include <boost/core/demangle.hpp>

//MOAB
#include <moab/ParallelComm.hpp>
#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/Range.hpp>
#include <MBTagConventions.hpp>
#include <moab/Skinner.hpp>
#include <moab/AdaptiveKDTree.hpp>
#include <moab/OrientedBoxTreeTool.hpp>
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
