cmake_minimum_required(VERSION 3.20)

project(polikit
    VERSION 0.4
    DESCRIPTION "polikit -- Polyhedral analysis module"
    LANGUAGES Fortran)

enable_language(Fortran)

file(GLOB SRC
  src/basic/*.f90
  src/func/*.f90
  src/io/*.f90
  src/main.f90
)

add_subdirectory(ext/toml-f)

option(DEBUG "Debug compile options" OFF)

if(DEBUG)
    add_compile_options(-fbounds-check)
    message(STATUS "Building in Debug mode")
else()
    add_compile_options(-O3)
    message(STATUS "Building in Release mode")
endif()

add_executable(polikit ${SRC})

# Optional dependencies
find_library(FFTW3_LIB NAMES fftw3)
find_path(FFTW3_INCLUDE NAMES fftw3.f03 fftw3.f PATH_SUFFIXES fftw3)
find_package(OpenMP QUIET)

if(FFTW3_LIB AND FFTW3_INCLUDE)
    message(STATUS "FFTW3 found -- enabling FFT support")
    target_link_libraries(polikit PRIVATE ${FFTW3_LIB})
    target_include_directories(polikit PRIVATE ${FFTW3_INCLUDE})
    target_compile_definitions(polikit PRIVATE HAVE_FFTW3)
else()
    message(STATUS "FFTW3 not found -- FFT support disabled")
endif()

# Get the number of CPU cores
if(NOT DEFINED OMP_THREAD_COUNT)
    cmake_host_system_information(RESULT CORE_COUNT QUERY NUMBER_OF_LOGICAL_CORES)
    set(OMP_THREAD_COUNT ${CORE_COUNT})
    message(STATUS "Detected CPU cores: ${CORE_COUNT}")
endif()

if(OpenMP_Fortran_FOUND)
    message(STATUS "OpenMP found -- enabling OpenMP support")
    add_compile_definitions(OMP_THREAD_COUNT=${OMP_THREAD_COUNT})
    target_link_libraries(polikit PRIVATE OpenMP::OpenMP_Fortran)
    target_compile_definitions(polikit PRIVATE HAVE_OPENMP)
else()
    message(STATUS "OpenMP not found -- OpenMP support disabled")
endif()

target_link_libraries(polikit PRIVATE toml-f)


