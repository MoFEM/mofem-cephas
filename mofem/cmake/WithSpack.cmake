option(WITHSPACK "Comopile with Spack" OFF)
if(WITH_SPACK)
  # enable @rpath in the install name for any shared library being built
  # note: it is planned that a future version of CMake will enable this by default
  set(CMAKE_MACOSX_RPATH 1)

  # Always use full RPATH
  # http://www.cmake.org/Wiki/CMake_RPATH_handling
  # http://www.kitware.com/blog/home/post/510

  # use, i.e. don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH  FALSE)

  # when building, don't use the install RPATH already
  # (but later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  endif("${isSystemDir}" STREQUAL "-1")
endif(WITH_SPACK)

