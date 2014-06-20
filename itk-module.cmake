set(DOCUMENTATION "This module contains a class to perform Mesh Laplacian Deformation" )

itk_module( LaplacianDeformation
  DEPENDS
    ITKQuadEdgeMeshFiltering
 TEST_DEPENDS
    ITKTestKernel
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
)
