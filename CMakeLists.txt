# Minimum CMake version, determined by hdf5.
cmake_minimum_required (VERSION 3.1)

# Version numbers.
set (ECOSYS_MAJOR_VERSION 0)
set (ECOSYS_MINOR_VERSION 1)
set (ECOSYS_PATCH_VERSION 0)
set (ECOSYS_VERSION "${ECOSYS_MAJOR_VERSION}.${ECOSYS_MINOR_VERSION}.${ECOSYS_PATCH_VERSION}")


# Adjust CMake's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Options for building Ecosys. These come from the xSDK compliance rules.
option(USE_XSDK_DEFAULTS "Set to use xSDK defaults for options [ON]." ON)
option(CMAKE_INSTALL_PREFIX "Sets installation prefix [/usr/local].")
option(XSDK_ENABLE_DEBUG "Enables Debug mode builds [OFF]." OFF)
option(BUILD_SHARED_LIBS "Builds shared libraries [ON]." ON)
option(XSDK_WITH_NETCDF "Enables support for netcdf [OFF]." ON)
option(TPL_NETCDF_LIBRARIES "List of absolute paths to netcdf link libraries [].")
option(TPL_NETCDF_INCLUDE_DIRS "List of absolute paths to netcdf include directories [].")



# For now, we disable shared libs on Macs.
if (APPLE)
  set(BUILD_SHARED_LIBS OFF)
endif()

if (BUILD_SHARED_LIBS)
	message("-- Ecosys will be built as a shared library.")
else()
	message("-- Ecosys will be built as a static library.")
endif()

if (NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()


# set the number of cores to something reasonable so we don't over
# subscribe a laptop or exceed hpc login node restrictions. could
# switch to something more sophisticated like polymec if necessary.
set(NUMBER_OF_CORES 4)
math(EXPR NUM_BUILD_THREADS "${NUMBER_OF_CORES} + 1")

include(set_up_platform)
include(set_up_compilers)
message("compilers set")
# Make sure compilers are set. This must be done before enabling languages.
if (NOT CMAKE_C_COMPILER)
  if (NOT $ENV{CC} STREQUAL "")
    set(CMAKE_C_COMPILER $ENV{CC})
  else()
    set(CMAKE_C_COMPILER cc)
  endif()
endif()
if (NOT CMAKE_C_FLAGS)
  set(CMAKE_C_FLAGS $ENV{CFLAGS})
endif()
if (NOT CMAKE_CXX_COMPILER)
  if (NOT $ENV{CXX} STREQUAL "")
    set(CMAKE_CXX_COMPILER $ENV{CXX})
  else()
    set(CMAKE_CXX_COMPILER c++)
  endif()
endif()
if (NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS $ENV{CXX_FLAGS})
endif()
if (NOT CMAKE_Fortran_COMPILER)
  if (NOT $ENV{FC} STREQUAL "")
    set(CMAKE_Fortran_COMPILER $ENV{FC})
  else()
    set(CMAKE_Fortran_COMPILER gfortran)
  endif()
endif()
if (NOT CMAKE_Fortran_FLAGS)
  set(CMAKE_Fortran_FLAGS $ENV{FCFLAGS})
endif()

enable_language(C)
enable_language(CXX)
message("set fortran")
enable_language(Fortran)
message("ok")

# We declare the project here.
project (ecosys)

message("-- C compiler is ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})")
message("-- CXX compiler is ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
message("-- Fortran compiler is ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})")

if (BGC)
  set (ECOSYS_BGC 1)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSYS_BGC")
else()
  set (ECOSYS_BGC 0)
endif()

if (ECOSYS)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DECOSYS")
endif()
# Figure out the system type.
set(ECOSYS_HAVE_BOOL 1) # All reasonable C99 compilers have this now.
if (APPLE EQUAL 1)
  set(SYS_FLAGS "-DAPPLE=1")
  set(DYLIB_SUFFIX "dylib")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Accelerate")
else ()
  if (LINUX EQUAL 1)
    set(SYS_FLAGS "-DLINUX=1")
    set(DYLIB_SUFFIX "so")
  else()
    if (WIN32 EQUAL 1)
      set(ECOSYS_HAVE_BOOL 0) # MS doesn't have reasonable C compilers.
      set(SYS_FLAGS "-DWINDOWS=1")
      set(DYLIB_SUFFIX "dll")
    endif()
  endif ()
endif ()

# Here we make sure CMake-installed binaries use the correct runpath, and
# that the path is not stripped during installation.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Basic libraries to be linked in.
set(ECOSYS_LIBRARIES m)

if (${NEED_LAPACK})
  # NEED_LAPACK is set by set_up_platform()
  include(FindBLAS)
  include(FindLAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if (${LAPACK_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${LAPACK_LIBRARY_DIR}")
  endif()
  if (${BLAS_LIBRARY_DIR})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${BLAS_LIBRARY_DIR}")
  endif()
  set(ECOSYS_LIBRARIES ${ECOSYS_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES})
endif()

set(ECOSYS_NEED_PETSC 0)

#if (XSDK_WITH_NETCDF AND MACHINE MATCHES "supported")

  if ((TPL_NETCDF_LIBRARIES MATCHES "lib") AND (TPL_NETCDF_INCLUDE_DIRS MATCHES "include"))
#    message(FATAL_ERROR "TPL_NETCDF_LIBRARIES option be set for netcdf support to be enabled.")
#    message("TPL_NETCDF_LIBRARIES option be set for netcdf support to be enabled.")
     message("ECOSYS has netcdf support from the system")
     message("TPL_NETCDF_LIBRARIES ${TPL_NETCDF_LIBRARIES}")
     message("TPL_NETCDF_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS}")
     list(APPEND ECOSYS_TPLS ${TPL_NETCDF_LIBRARIES})
     list(APPEND ECOSYS_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS})
     set(ECOSYS_HAVE_NETCDF 1)
  else()
     message("ECOSYS does not have netcdf support from the system")
     set(ECOSYS_HAVE_NETCDF 0)
  endif()
#X#  foreach(lib ${TPL_NETCDF_LIBRARIES})
#X#    if (NOT EXISTS ${lib})
#X#      message(FATAL_ERROR "netcdf library not found: ${lib}")
#X#    endif()
#X#  endforeach()

#  if (NOT TPL_NETCDF_INCLUDE_DIRS)
#    message(FATAL_ERROR "TPL_NETCDF_INCLUDE_DIRS option be set for netcdf support to be enabled.")
#  endif()
#X#  foreach(dir ${TPL_NETCDF_INCLUDE_DIRS})
#X#    if (NOT EXISTS ${dir})
#X#      message(FATAL_ERROR "netcdf include directory not found: ${dir}")
#X#    endif()
#X#  endforeach()

#  message("-- Enabled support for netcdf.")
#  list(APPEND ECOSYS_TPLS ${TPL_NETCDF_LIBRARIES})
#  list(APPEND ECOSYS_INCLUDE_DIRS ${TPL_NETCDF_INCLUDE_DIRS})
#  set(ECOSYS_HAVE_NETCDF 1)
#else()
#  set(ECOSYS_HAVE_NETCDF 0)
#endif()
# Other third-party libraries.
#add_subdirectory(3rd-party)

# Include the binary directory in the header file search path,
# since it's where we place the third-party libraries.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/include")
link_directories("${PROJECT_BINARY_DIR}/lib")
include_directories(${ECOSYS_INCDIRS})

# Unit testing.
enable_testing()

# Source code itself.
include_directories("${PROJECT_SOURCE_DIR}")
add_subdirectory(f77src)

# Drivers for benchmarks.
#X#add_subdirectory(drivers)

# Benchmarks.
#X#add_subdirectory(benchmarks)

# Now that we have gathered all our libraries, generate an ecosys.cmake
# file that contains all the vital information.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Templates/ecosys.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/ecosys.cmake"
  @ONLY
)

# Install miscellaneous build/test files.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ecosys.cmake DESTINATION share/ecosys)
