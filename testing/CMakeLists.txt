include_directories( ${LaplacianDeformation_SOURCE_DIR}/include )

add_executable(HeaderTest HeaderTest.cxx)
target_link_libraries(HeaderTest ${ITK_LIBRARIES})

add_executable(deformation deformation.cxx)
target_link_libraries(deformation ${ITK_LIBRARIES})

add_executable(deformation2 deformation2.cxx)
target_link_libraries(deformation2 ${ITK_LIBRARIES})
