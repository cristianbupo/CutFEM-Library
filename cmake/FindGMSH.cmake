

find_path(GMSH_INCLUDE_DIR
    NAMES
        gmsh.h
    PATHS
    /usr/local/include
)

find_path(GMSH_LIBRARY_DIR
    NAMES 
        libgmsh.dylib
    PATHS
    /usr/local/lib
)

if(GMSH_LIBRARY_DIR)
    find_library(GMSH_LIBRARY
        NAMES 
            libgmsh.dylib
        PATHS 
            ${GMSH_LIBRARY_DIR}
            NO_DEFAULT_PATH
    )

message(STATUS "gmsh Found ")
message( STATUS "inlude directory: ${GMSH_INCLUDE_DIR}" )
message( STATUS "lib directory: ${GMSH_LIBRARY_DIR} ")
else()
message(STATUS " gmsh Not Found ")
endif()