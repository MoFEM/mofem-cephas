
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

#approximation
add_library(mofem_approx
  impl/EntityRefine.cpp
  c_impl/fem_tools.c
  c_impl/h1.c
  c_impl/hdiv.c
  c_impl/l2.c
)

install(TARGETS mofem_approx DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

#install(FILES c/fem_tools.h c/h1_hdiv_hcurl_l2.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
