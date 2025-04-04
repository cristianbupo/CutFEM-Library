## Dependencies
find_package(MPI)


set(PARALLEL_SOURCES
  cfmpi.cpp
)

## ######################
## Build rules
## ######################
add_library(parallel SHARED ${PARALLEL_SOURCES})
add_library(parallel::parallel ALIAS parallel)


target_include_directories(parallel
                           PRIVATE 
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           ${MPI_CXX_INCLUDE_PATH})

# target_link_libraries(parallel /opt/homebrew/opt/libomp/lib/libomp.dylib)

