cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PANPP VERSION 2.0.0)
add_executable(PANPP main_index.cpp main_search.cpp main_fmd.cpp main.cpp)
set_target_properties(PANPP PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(PANPP PRIVATE -Wall -Wextra)
add_definitions(-DMASSIVE_DATA_RLCSA)

find_package(OpenMP REQUIRED)
# find_library(SPDLOG libspdlog.so REQUIRED)

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

# spdlog
message(STATUS "spdlog will be built from source")
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG 7e635fca68d014934b4af8a1cf874f63989352b7) # v1.12.0
FetchContent_MakeAvailable(spdlog)

# ropebwt2
message(STATUS "ropebwt2 will be built from source")
ExternalProject_Add(ropebwt
    GIT_REPOSITORY https://github.com/lh3/ropebwt2.git
    GIT_TAG bd8dbd3db2e9e3cff74acc2907c0742c9ebbf033
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    # PATCH_COMMAND bash ${CMAKE_SOURCE_DIR}/patches/patch_ropebwt.sh <SOURCE_DIR>/mrope.c # I haven't been able to run sed here.. I could use my branch btw
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
      "CC=${CMAKE_C_COMPILER}"
      "CFLAGS=-g -Wall -O2 ${CONDA_INCLUDE}"
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(ropebwt SOURCE_DIR)
SET(ROPEBWT_SOURCE_DIR ${SOURCE_DIR})
SET(ROPEBWT_INCLUDE_DIR ${ROPEBWT_SOURCE_DIR})
set(ropebwt-src ${ROPEBWT_SOURCE_DIR}/mrope.c ${ROPEBWT_SOURCE_DIR}/rope.c ${ROPEBWT_SOURCE_DIR}/rld0.c ${ROPEBWT_SOURCE_DIR}/rle.c)
set_source_files_properties(${ropebwt-src} PROPERTIES GENERATED TRUE)
add_library(ROPEBWT OBJECT ${ropebwt-src})
add_dependencies(ROPEBWT ropebwt)

add_dependencies(PANPP rlcsa ropebwt spdlog)

target_include_directories(PANPP
  PRIVATE ${RLCSA_INCLUDE_DIR}
  PRIVATE ${ROPEBWT_INCLUDE_DIR}
  PRIVATE ${SPDLOG_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(PANPP
  PUBLIC RLCSA
  PUBLIC ROPEBWT
  PUBLIC spdlog::spdlog
  PUBLIC z
  PUBLIC OpenMP::OpenMP_CXX
)
