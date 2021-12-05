function(add_ecosys_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${ECOSYS_LIBRARIES})
endfunction(add_ecosys_executable)

