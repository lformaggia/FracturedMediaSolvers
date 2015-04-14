# Set the minimum required version of cmake for a project
cmake_minimum_required(VERSION 2.8)

add_executable( polyGen.exe polyGen )
target_link_libraries( polyGen.exe fvcode3d tet ${LIBRARIES})

# Add data file
add_custom_command(
    TARGET polyGen.exe POST_BUILD

    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/tet2foam.txt
    ${CMAKE_CURRENT_BINARY_DIR}/tet2foam.txt

    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/foam2fvg.txt
    ${CMAKE_CURRENT_BINARY_DIR}/foam2fvg.txt

    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/tet2fvg.sh
    ${CMAKE_CURRENT_BINARY_DIR}/tet2fvg.sh

    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/data
    ${CMAKE_CURRENT_BINARY_DIR}/data


    COMMAND ${CMAKE_COMMAND} -E make_directory results

    )
