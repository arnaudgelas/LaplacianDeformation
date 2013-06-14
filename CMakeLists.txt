project( LaplacianDeformation )

cmake_minimum_required( VERSION 2.8 )

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${LaplacianDeformation_BINARY_DIR}/bin
)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

include( CTest )

add_subdirectory( testing )

