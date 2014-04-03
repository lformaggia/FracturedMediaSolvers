# Set the minimum required version of cmake for a project
cmake_minimum_required(VERSION 2.8)

add_executable(cartGen.exe cartGen )
target_link_libraries( cartGen.exe fvcode tet )

# Add data file
add_custom_command(
	TARGET cartGen.exe POST_BUILD        	# Adds a post-build event to convToSolver.exe
	
	COMMAND ${CMAKE_COMMAND} -E copy	  	# which executes "cmake - E copy..."
	${CMAKE_CURRENT_SOURCE_DIR}/data.txt     	# <--this is in-file
	${CMAKE_CURRENT_BINARY_DIR}/data.txt		# <--this is out-file path
	
	COMMAND ${CMAKE_COMMAND} -E make_directory data

	COMMAND ${CMAKE_COMMAND} -E make_directory results

	)
