# We create and install ${PROJECT_NAME}Config.cmake and
# ${PROJECT_NAME}ConfigVersion.cmake so that ${PROJECT_NAME} can be used
# by other projects using CMake
 
# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Generate a string containing this project libraries
STRING(REPLACE ";" " " CONF_LIBRARIES "${${PROJECT_NAME}_LIBS_LIBRARIES}")

# We use variables agnostic of the project name,
# so that we don't have to change these for every project
set(PROJECT_TPL_INCLUDE_DIRS ${${PROJECT_NAME}_TPL_INCLUDE_DIRS})
set(PROJECT_TPL_LIBRARY_DIRS ${${PROJECT_NAME}_TPL_LIBRARY_DIRS})
# Create the ${PROJECT_NAME}Config.cmake and ${PROJECT_NAME}ConfigVersion files
set(CMAKE_FILES_DIRECTORY "")
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH REL_LIBRARY_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_LIB_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/cmake/ProjectConfig.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${${PROJECT_NAME}_CMAKE_DIR}/${REL_INCLUDE_DIR}")
set(CONF_LIBRARY_DIRS "\${${PROJECT_NAME}_CMAKE_DIR}/${REL_LIBRARY_DIR}")
configure_file(${PROJECT_SOURCE_DIR}/cmake/ProjectConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake" @ONLY)
# ... for both
configure_file(${PROJECT_SOURCE_DIR}/cmake/ProjectConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake" @ONLY)
 
# Install the ${PROJECT_NAME}Config.cmake and ${PROJECT_NAME}ConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
