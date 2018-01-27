# Add all the source files to the add_executable command
# Replace exactSolution_test with the name of the executable to be generated
# (it will be both the name of the executable and the name of
# the target, so it must be unique within the whole library)
add_executable(tet_test main.cpp)

# This allows the test to be executed with ctest
add_test(NAME TetTest
         COMMAND tet_test) # this is the command that will be executed to run the test
target_link_libraries(tet_test ${LIBRARIES})

# If the test needs some files to be copied to the build directory
# we can specify them here. These files will be copied to this build
# directory regardless of their source directory
copy_files_to_binary_dir(tet_test_files # target name
#    ./cart100x100x1.fvg
    ./cubeFracturedTet.fvg
    ./cubeFracturedTet2.fvg
    ./cubeFracturedTet3.fvg
    ./data.txt
)
