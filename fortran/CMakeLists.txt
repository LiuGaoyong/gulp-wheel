cmake_minimum_required(VERSION 3.15)
project(GULP LANGUAGES Fortran C)
option(BUILD_DEBUG          "Build debug version"                       OFF)
option(USE_PLUMED           "Enable PLUMED support"                     OFF)
option(USE_FFTW3            "Enable FFTW3 support"                      OFF)
option(USE_KIM              "Enable OpenKIM support"                    OFF)
option(USE_ELPA             "Enable ELPA support"                       OFF)


# Set Compiler & Flags
set(CMAKE_Fortran_COMPILER "gfortran")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffpe-summary=invalid,zero,overflow")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wunused-parameter -fbounds-check -Wuninitialized -fbacktrace")
set(FORTRAN_FLAGS_EXTRA "-fallow-argument-mismatch")


# Global definations
add_compile_definitions(FLUSH)
if(NOT USE_KIM)
    add_compile_definitions(NOFOX)
endif()
if(USE_PLUMED)
    add_compile_definitions(PLUMED)
endif()
if(USE_FFTW3)
    add_compile_definitions(FFTW3)
endif()
if(USE_KIM)
    add_compile_definitions(KIM)
endif()
if(USE_ELPA)
    add_compile_definitions(ELPA)
endif()


# collect all source files
file(GLOB_RECURSE SOURCES "Src/*.F90" "Src/*.c")
file(GLOB_RECURSE LIB_SOURCES "Utils/pGFNFF/Src/*.F90")
list(FILTER SOURCES EXCLUDE REGEX ".*/FoX-.*")  # exclude FoX files


# build GFNFF library
add_library(gfnff STATIC ${LIB_SOURCES})
set_target_properties(gfnff PROPERTIES POSITION_INDEPENDENT_CODE ON)
# main executable files
add_executable(gulp ${SOURCES})
target_link_libraries(gulp PRIVATE gfnff ${MATH_LIBS})
if(USE_ELPA)
    target_link_libraries(gulp PRIVATE elpa)
endif()
if(USE_FFTW3)
    find_package(FFTW3 REQUIRED)
    target_link_libraries(gulp PRIVATE FFTW3::FFTW3)
endif()
if(USE_KIM)
    find_package(OpenKIM REQUIRED)
    target_link_libraries(gulp PRIVATE OpenKIM::OpenKIM)
endif()