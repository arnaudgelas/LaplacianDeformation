include_directories( ${LaplacianDeformation_SOURCE_DIR}/include )

add_executable(HeaderTest HeaderTest.cxx)
target_link_libraries(HeaderTest ${ITK_LIBRARIES})

add_test( NAME HeaderTest
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/HeaderTest
)

add_executable(softconstraints softconstraints.cxx)
target_link_libraries(softconstraints ${ITK_LIBRARIES})

add_test( NAME SoftConstraintsTest
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/softconstraints
          ${CMAKE_CURRENT_SOURCE_DIR}/input.vtk
          SoftConstraintsTestOutput.vtk
)

add_executable(hardconstraints hardconstraints.cxx)
target_link_libraries(hardconstraints ${ITK_LIBRARIES})

add_test( NAME HardConstraintsTest
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hardconstraints
          ${CMAKE_CURRENT_SOURCE_DIR}/input.vtk
          HardConstraintsTestOutput.vtk
)

