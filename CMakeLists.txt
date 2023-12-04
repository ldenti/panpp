cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PANPP VERSION 1.0.0)
add_executable(PANPP rlcsa2/alphabet.cpp rlcsa2/misc/utils.cpp rlcsa2/bits/bitbuffer.cpp rlcsa2/bits/bitvector.cpp rlcsa2/bits/deltavector.cpp rlcsa2/bits/rlevector.cpp rlcsa2/suffixarray.cpp rlcsa2/rlcsa.cpp rlcsa2/rlcsa_builder.cpp rlcsa2/fmd.cpp main_index.cpp main_search.cpp main.cpp) 
set_target_properties(PANPP PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(PANPP PRIVATE -Wall -Wextra)
add_definitions(-DMASSIVE_DATA_RLCSA)

find_package(OpenMP REQUIRED)

# spdlog
message(STATUS "spdlog will be built from source")
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG 7e635fca68d014934b4af8a1cf874f63989352b7) # v1.12.0
FetchContent_MakeAvailable(spdlog)

# sdsl
message(STATUS "sdsl-lite will be built from source")
ExternalProject_Add(sdsl
    GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
    # GIT_TAG 14ae2b430a19a0f8602dc08c1d38906607e73507
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND bash install.sh .
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(sdsl SOURCE_DIR)
SET(SDSL_SOURCE_DIR ${SOURCE_DIR})
SET(SDSL_INCLUDE_DIR ${SDSL_SOURCE_DIR}/include)
add_library(SDSL STATIC IMPORTED)
set_target_properties(SDSL PROPERTIES IMPORTED_LOCATION ${SDSL_SOURCE_DIR}/lib/libsdsl.a)

add_dependencies(PANPP spdlog sdsl)

target_include_directories(PANPP
  PRIVATE ${SDSL_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_libraries(PANPP
  PUBLIC spdlog::spdlog
  PUBLIC SDSL
  PUBLIC z
  PUBLIC OpenMP::OpenMP_CXX
)
