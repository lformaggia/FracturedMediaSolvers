# Add all the source files to the add_executable command
# Replace exactSolution_test with the name of the executable to be generated
# (it will be both the name of the executable and the name of
# the target, so it must be unique within the whole library)
#add_executable(network1_test main.cpp)

# This allows the test to be executed with ctest
add_test(NAME Network1Test
  COMMAND ../../fvcode3d.exe) # this is the command that will be executed to run the test

#target_link_libraries(network1_test
#  ${LIBRARIES}
#  dl) #dl is for dynamic loading

#library for user defined functions
add_library(fun1 SHARED
  ./functions.cpp
  )

target_link_libraries(fun1 ${LIBRARIES})


# If the test needs some files to be copied to the build directory
# we can specify them here. These files will be copied to this build
# directory regardless of their source directory
copy_files_to_binary_dir(network1_test_files # target name
#    ./cart100x100x1.fvg
#    ./cart50x50x1.fvg
#    ./cart200x200x1.fvg
./compMFD.fvg
./data.txt
)
