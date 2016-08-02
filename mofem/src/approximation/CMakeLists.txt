
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# MoFEM is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MoFEM. If not, see <http://www.gnu.org/licenses/>

include_directories(${PROJECT_SOURCE_DIR}/src/ftensor/src)
include_directories(${PROJECT_SOURCE_DIR}/src/multi_indices)
include_directories(${PROJECT_SOURCE_DIR}/src/interfaces)
include_directories(${PROJECT_SOURCE_DIR}/src/finite_elements)
include_directories(${PROJECT_SOURCE_DIR}/src/petsc)

# Approximation lib
add_library(mofem_approx
  impl/EntityRefine.cpp
  impl/BaseFunction.cpp
  impl/LegendrePolynomial.cpp
  impl/LobattoPolynomial.cpp
  impl/EntPolynomialBaseCtx.cpp
  impl/TriPolynomialBase.cpp
  impl/TetPolynomialBase.cpp
  impl/EdgePolynomialBase.cpp
  impl/FlatPrismPolynomialBase.cpp
  impl/FatPrismPolynomialBase.cpp
  impl/Hcurl.cpp
  c_impl/base_functions.c
  c_impl/fem_tools.c
  c_impl/h1.c
  c_impl/hdiv.c
  c_impl/l2.c
)

install(TARGETS mofem_approx DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

add_executable(vtk_curl_base_on_tet impl/Hcurl.cpp)
set_target_properties(vtk_curl_base_on_tet PROPERTIES COMPILE_FLAGS "-DGENERATE_VTK_WITH_CURL_BASE")
target_link_libraries(vtk_curl_base_on_tet
  mofem_finite_elements
  mofem_interfaces
  mofem_multi_indices
  mofem_petsc
  mofem_approx
  mofem_cblas
  mofem_third_party
  ${PROJECT_LIBS}
)
