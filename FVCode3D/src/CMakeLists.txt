# Create static library
add_subdirectory(FVCode3D)
SET(LIBRARIES fvcode3d ${LIBRARIES} CACHE INTERNAL "")

add_executable(fvcode3d.exe main.cpp )
# Create executable and link libraries
target_link_libraries(fvcode3d.exe
    fvcode3d
    ${LIBRARIES}
    )

# Add subdirectories with tests
if(${PROJECT_NAME}_ENABLE_TESTS)
    add_subdirectory(test)
endif()

# Add tools
add_subdirectory (tools)

# Install
INSTALL(TARGETS fvcode3d.exe DESTINATION ${INSTALL_BIN_DIR})

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../lib/libtet.a DESTINATION ${INSTALL_LIB_DIR})

# Add data files
link_files_to_binary_dir(data_files data.txt ../data)

add_custom_command(
    TARGET fvcode.exe POST_BUILD                # Adds a post-build event to fvcode.exe
    COMMAND ${CMAKE_COMMAND} -E make_directory results
    )

