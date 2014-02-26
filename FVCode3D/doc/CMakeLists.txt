# find the Doxygen package add
find_package(Doxygen)

if(DOXYGEN_FOUND)
	# generate Doxyfile
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
               	   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

	# add doc target
	add_custom_target(doc
			${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			COMMENT "Generating API documentation with Doxygen" VERBATIM
			)       
endif(DOXYGEN_FOUND)
