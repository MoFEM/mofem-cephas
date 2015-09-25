# - Try to find LASSO

# If unset, try environment

if(NOT LASSO_DIR)
  set(LASSO_DIR $ENV{LASSO_DIR})
endif(NOT LASSO_DIR)

find_library(LASSO_LIBRARY NAMES iRel PATHS "${LASSO_DIR}/lib")
message(STATUS ${LASSO_LIBRARY})

if(LASSO_LIBRARY) 
  include_directories("${LASSO_DIR}/include")
endif(LASSO_LIBRARY)




