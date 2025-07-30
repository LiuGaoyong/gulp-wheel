cmake_minimum_required(VERSION 3.15)
project(GULP LANGUAGES Fortran C)

option(BUILD_LIB            "Build as library"                          OFF)
option(BUILD_DEBUG          "Build debug version"                       OFF)
option(BUILD_LIB_NO_MATH    "Build library without math dependencies"   OFF)
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

# 包含目录
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    # 用户自定义路径示例:
    # /usr/local/include/elpa-2020.05.001/modules/
)

# 收集所有源文件
file(GLOB_RECURSE SOURCES "Src/*.F90" "Src/*.c")
file(GLOB_RECURSE LIB_SOURCES "Utils/pGFNFF/Src/*.F90")

# 排除特定文件（如FoX）
list(FILTER SOURCES EXCLUDE REGEX ".*/FoX-.*")

# 构建目标配置
if(BUILD_LIB OR BUILD_LIB_NO_MATH)
    if(BUILD_LIB_NO_MATH)
        add_library(gulp_nomaths STATIC ${SOURCES})
        target_link_libraries(gulp_nomaths PRIVATE gfnff)
        set_target_properties(gulp_nomaths PROPERTIES OUTPUT_NAME "gulp")
    else()
        add_library(gulp STATIC ${SOURCES})
        target_link_libraries(gulp PRIVATE ${MATH_LIBS})
    endif()
else()
    # 构建GFNFF库
    add_library(gfnff STATIC ${LIB_SOURCES})
    set_target_properties(gfnff PROPERTIES POSITION_INDEPENDENT_CODE ON)

    # 主可执行文件
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
endif()


# 安装规则
install(TARGETS gulp DESTINATION bin)
install(FILES ${CMAKE_BINARY_DIR}/gulp DESTINATION bin)