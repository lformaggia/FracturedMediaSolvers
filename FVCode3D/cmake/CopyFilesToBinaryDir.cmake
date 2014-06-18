# Copy a file into the current binary dir even if the file is elsewhere in the source tree
FUNCTION(COPY_FILES_TO_BINARY_DIR TARGET_NAME)
    set(FEST_FILES_LIST)
    foreach(f ${ARGN})
        # We want to get the file name (stripping the path)
        # because we want to copy it to the current binary directory
        # even if the file is elsewhere in the source tree
        get_filename_component(f_NAME ${f} NAME)
        add_custom_command(
        #    TARGET op PRE_BUILD
            OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${f_NAME}
            COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${f} ${CMAKE_CURRENT_BINARY_DIR}/${f_NAME}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${f}
            )
        SET(DEST_FILE_FULL "${CMAKE_CURRENT_BINARY_DIR}/${f_NAME}")
        LIST(APPEND DEST_FILES_LIST "${DEST_FILE_FULL}")
    endforeach()
    ADD_CUSTOM_TARGET(${TARGET_NAME} ALL DEPENDS ${DEST_FILES_LIST}) 

ENDFUNCTION()

