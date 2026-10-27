cmake_minimum_required(VERSION 3.18)
project(quarcnonlinear LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python 3.9 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(nanobind CONFIG REQUIRED)

# src/core is a copy of the repository's cpp/ (tools/sync_core.sh).
file(GLOB CORE_SOURCES CONFIGURE_DEPENDS src/core/core_*.cpp)
nanobind_add_module(_core NB_STATIC src/bindings.cpp ${CORE_SOURCES})
target_include_directories(_core PRIVATE src/core)
# Keep a*b+c as two roundings, as MATLAB computes it; see cpp/README.md.
if(NOT MSVC)
  target_compile_options(_core PRIVATE -ffp-contract=off)
endif()
install(TARGETS _core LIBRARY DESTINATION quarcnonlinear)
