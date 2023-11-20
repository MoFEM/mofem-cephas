configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_build.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_build.sh"
  @ONLY
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_fast_check.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_fast_check.sh"
  @ONLY
  )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mofem_update.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mofem_update.sh"
  @ONLY
  )

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_build.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_fast_check.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/mofem_update.sh
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_EXECUTE
)
