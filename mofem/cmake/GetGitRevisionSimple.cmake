find_package(Git)

function(get_git_hash GIT_DIR _hashvar)
  execute_process(COMMAND
    "${GIT_EXECUTABLE}"  rev-parse  HEAD
    WORKING_DIRECTORY ${GIT_DIR}
    OUTPUT_VARIABLE HEAD_HASH
    RESULT_VARIABLE res)
  if(NOT ${res})
    string(REGEX REPLACE "\n$" "" HEAD_HASH "${HEAD_HASH}")
    set(${_hashvar} "${HEAD_HASH}" PARENT_SCOPE)
  else(NOT ${res})
    set(${_hashvar} "SHA1-NOT FOUND" PARENT_SCOPE)
  endif(NOT ${res})
endfunction()

function(get_git_tag GIT_DIR FALLBACK _gittag) 
  execute_process(COMMAND
    "${GIT_EXECUTABLE}" describe --tags 
    WORKING_DIRECTORY ${GIT_DIR}
    OUTPUT_VARIABLE GIT_TAG
    RESULT_VARIABLE res)
  if(NOT ${res})
    string(REGEX REPLACE "\n$" "" GIT_TAG "${GIT_TAG}")
    set(${_gittag} "${GIT_TAG}" PARENT_SCOPE)
  else(NOT ${res})
    set(${_gittag} "${FALLBACK}-fallback" PARENT_SCOPE)
  endif(NOT ${res}) 
endfunction()

function(get_git_version
  GIT_TAG_VERSION _version_major _version_minor _version_build) 
  string(REGEX REPLACE 
    "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${GIT_TAG_VERSION}")
  string(REGEX REPLACE 
    "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${GIT_TAG_VERSION}")
  string(REGEX REPLACE 
    "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_BUILD "${GIT_TAG_VERSION}")
  set(${_version_major} "${VERSION_MAJOR}" PARENT_SCOPE)
  set(${_version_minor} "${VERSION_MINOR}" PARENT_SCOPE)
  set(${_version_build} "${VERSION_BUILD}" PARENT_SCOPE)
endfunction()