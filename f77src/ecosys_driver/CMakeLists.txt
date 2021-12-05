set(ECOSYS_DRIVER_SOURCES
	soil.f
	main.f
)

include_directories(${ECOSYS_INCLUDE_DIRS})
include_directories(${CMAKE_BINARY_DIR}/f77src/ecosys_core/)

include(add_ecosys_library)
add_ecosys_library(ecosys_driver ${ECOSYS_DRIVER_SOURCES})

set(ECOSYS_LIBRARIES ecosys_driver;${ECOSYS_LIBRARIES} PARENT_SCOPE)

link_libraries(${ECOSYS_LIBRARIES};${ECOSYS_TPLS})
include(add_ecosys_executable)
add_ecosys_executable(ecosys.x main.f soil.f)

install(TARGETS ecosys.x DESTINATION bin)
