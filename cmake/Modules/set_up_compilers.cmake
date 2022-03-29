# This macro sets compiler specific flags
macro(set_up_compilers)
  #
  # General C compiler flags.
  #
  message("HOSTNAME ${HOSTNAME}")
  message("CMAKE_C_COMPILER_ID ${CMAKE_C_COMPILER_ID}")
  if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function -Wno-cast-qual")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-but-set-variable -Wno-overlength-strings")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-discarded-qualifiers")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-conversion -Wno-maybe-uninitialized")

    if (BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fPIC")
    endif()

    if (LINUX EQUAL 1)
      # Counter some of GCC's more recent stinginess on Linux.
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_POSIX_C_SOURCE=200809L")
    endif()

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic-errors -Wextra")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror-implicit-function-declaration ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-builtin")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-sign-compare -Wno-unused-parameter ")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-function")

  elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
    if (CMAKE_BUILD_TYPEi MATCHES "Debug")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -g")
    else()
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
    endif()
  endif()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SYS_FLAGS}")

  #
  # Fortran compiler flags.
  #
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -W -Wall -std=gnu -pedantic -finit-local-zero -Wall -cpp -O2")
#    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRGNU -fdefault-real-8")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-variable")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-parameter -Wno-unused-function")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wuninitialized -Wno-unused-dummy-argument")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=use-without-only -ffree-line-length-none")

  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCPRINTEL")

    if (HOSTNAME MATCHES "cori")
      message("Using set up for cori")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
    elseif (HOSTNAME MATCHES "[scs]")
      message("Using set up for lawrencium")
      if (CMAKE_BUILD_TYPEi MATCHES "Debug")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -debug -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
      endif()
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true")
    endif()
  endif()

endmacro()
