# Copyright (C) 2013, Lukasz Kaczmarczyk (likask AT wp.pl)
# The MoFEM package is copyrighted by Lukasz Kaczmarczyk.
# It can be freely used for educational and research purposes
# by other institutions. If you use this softwre pleas cite my work.
#
# MoFEM is free software: you can redistribute it and/or modify it under
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

include_directories("${PROJECT_SOURCE_DIR}/third_party/cblas")
link_directories("${PROJECT_BINARY_DIR}/third_party/cblas")

if(CBLAS_LIB AND BLAS_FOUND)
  add_library(mofem_cblas all.c)
elseif(BLAS_FOUND)
  add_library(mofem_cblas all.c)
else(CBLAS_LIB)
  add_library(mofem_cblas
    all.c
    cdotcsub.f
    cdotusub.f
    dasumsub.f
    ddotsub.f
    dnrm2sub.f
    dsdotsub.f
    dzasumsub.f
    dznrm2sub.f
    icamaxsub.f
    idamaxsub.f
    isamaxsub.f
    izamaxsub.f
    sasumsub.f
    scasumsub.f
    scnrm2sub.f
    sdotsub.f
    sdsdotsub.f
    snrm2sub.f
    zdotcsub.f
    zdotusub.f
  )
endif(CBLAS_LIB AND BLAS_FOUND)
# target_compile_options(mofem_cblas PUBLIC "-Wno-macro-redefined")

#install(FILES cblas.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
install(TARGETS mofem_cblas DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
