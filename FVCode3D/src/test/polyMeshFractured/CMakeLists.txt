# Add all the source files to the add_executable command
# Replace exactSolution_test with the name of the executable to be generated
# (it will be both the name of the executable and the name of
# the target, so it must be unique within the whole library)
add_executable(polyMeshFractured_test main.cpp)

# This allows the test to be executed with ctest
add_test(NAME PolyMeshFracturedTest
         COMMAND polyMeshFractured_test) # this is the command that will be executed to run the test
target_link_libraries(polyMeshFractured_test ${LIBRARIES})

# If the test needs some files to be copied to the build directory
# we can specify them here. These files will be copied to this build
# directory regardless of their source directory
copy_files_to_binary_dir(polyMeshFractured_test_files # target name
    ./cubePolyF1.fvg
    ./cubePolyF2.fvg
    ./cubePolyF3.fvg
    ./cubePolyF4.fvg
    ./data.txt
)
