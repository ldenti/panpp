cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PANPP VERSION 2.0.0)
add_executable(PANPP main.cpp)
set_target_properties(PANPP PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(PANPP PRIVATE -Wall -Wextra)
add_definitions(-DMASSIVE_DATA_RLCSA)

# find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_library(SPDLOG libspdlog.so REQUIRED)

# rlcsa
message(STATUS "rlcsa will be built from source")
ExternalProject_Add(rlcsa
    GIT_REPOSITORY https://github.com/ldenti/rlcsa.git
    # GIT_TAG 14ae2b430a19a0f8602dc08c1d38906607e73507
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(rlcsa SOURCE_DIR)
SET(RLCSA_SOURCE_DIR ${SOURCE_DIR})
SET(RLCSA_INCLUDE_DIR ${RLCSA_SOURCE_DIR})
add_library(RLCSA STATIC IMPORTED)
set_target_properties(RLCSA PROPERTIES IMPORTED_LOCATION ${RLCSA_SOURCE_DIR}/librlcsa.a)
add_dependencies(RLCSA rlcsa)

add_dependencies(PANPP rlcsa)

target_include_directories(PANPP
  PRIVATE ${RLCSA_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(PANPP
  PUBLIC RLCSA
  PUBLIC spdlog
  PUBLIC z
  PUBLIC OpenMP::OpenMP_CXX
)
