cmake_minimum_required(VERSION 3.20)

project(polikit
    VERSION 0.4
    DESCRIPTION "polikit -- Polyhedral analysis module"
    LANGUAGES Fortran
)

# Source files
file(GLOB SRC
    src/basic/*.f90
    src/func/*.f90
    src/io/*.f90
    src/main.f90
)

# Dependencies
add_subdirectory(ext/toml-f)

# Executable
add_executable(polikit ${SRC})

# Compiler options
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Building in Debug mode")

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(polikit PRIVATE
            -O0
            -g
            -Wall
            -Wextra
            -fcheck=all
            -fbacktrace
            -finit-real=snan
            -finit-integer=-999999
            -ffpe-trap=invalid,zero,overflow
        )
    endif()

else()
    message(STATUS "Building in Release mode")

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(polikit PRIVATE
            -O3
        )
    endif()

endif()

# Optional FFTW support
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

find_library(FFTW3_LIB NAMES fftw3)
find_path(
    FFTW3_INCLUDE
    NAMES fftw3.f03 fftw3.f
    PATH_SUFFIXES fftw3
)

if(FFTW3_LIB AND FFTW3_INCLUDE)
    message(STATUS "FFTW3 found -- enabling FFT support")

    target_link_libraries(polikit PRIVATE ${FFTW3_LIB})
    target_include_directories(polikit PRIVATE ${FFTW3_INCLUDE})
    target_compile_definitions(polikit PRIVATE HAVE_FFTW3)

else()
    message(STATUS "FFTW3 not found -- FFT support disabled")

endif()

# Link libraries
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

target_link_libraries(polikit PRIVATE toml-f)
