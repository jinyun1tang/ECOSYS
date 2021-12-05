function(add_ecosys_library lib)
  add_library(${lib} ${ARGN})
  if (BUILD_SHARED_LIBS)
	  target_link_libraries(${lib} ${ECOSYS_LIBRARIES})
  endif()
endfunction(add_ecosys_library)

