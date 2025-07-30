cmake_minimum_required(VERSION 3.15)
project(GULP LANGUAGES Fortran C)

add_subdirectory(fortran)
install(TARGETS gulp DESTINATION gulp)