# Set the minimum required version of cmake for a project
cmake_minimum_required(VERSION 2.8)

# Add subdirectories with examples
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/converterToVTU)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/converterToSolver)
