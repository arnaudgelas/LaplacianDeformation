itk_module_test()

CreateTestDriver( LaplacianDeformation "${LaplacianDeformation-Test_LIBRARIES}" "${LaplacianDeformationTests}" )

set( ITK_NO_IO_FACTORY_REGISTER_MANAGER "" )

add_executable(softconstraints softconstraints.cxx)
target_link_libraries(softconstraints ${LaplacianDeformation_LIBRARIES})

add_test( NAME SoftConstraintsTest
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/softconstraints
          ${CMAKE_CURRENT_SOURCE_DIR}/input.vtk
          SoftConstraintsTestOutput.vtk
)

add_executable(hardconstraints hardconstraints.cxx)
target_link_libraries(hardconstraints ${LaplacianDeformation_LIBRARIES})

add_test( NAME HardConstraintsTest
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hardconstraints
          ${CMAKE_CURRENT_SOURCE_DIR}/input.vtk
          HardConstraintsTestOutput.vtk
)

