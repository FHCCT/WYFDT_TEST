cmake_minimum_required(VERSION 3.13)
project(dt LANGUAGES C CXX)

option(DT_BUILD_EXECUTABLE "Build the dt command-line program" ON)

# Source discovery requires the build-system regeneration check.
set(CMAKE_SUPPRESS_REGENERATION OFF)

find_package(OpenMP REQUIRED COMPONENTS CXX)

file(GLOB_RECURSE DT_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hxx"
)
list(REMOVE_ITEM DT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/dt_main.cpp")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/src" PREFIX src FILES ${DT_SOURCES})

add_library(dt STATIC ${DT_SOURCES}
    extern/geom/intersection_func.cpp
    extern/geom/predicates.cpp
    extern/geom/ttotr.c
)
add_library(dt::dt ALIAS dt)

target_compile_features(dt PUBLIC cxx_std_11)
target_compile_definitions(dt PUBLIC DT_LIBRARY)
target_include_directories(dt PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/geom
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/tree
)
target_include_directories(dt SYSTEM PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/spdlog/include
)
target_link_libraries(dt PUBLIC OpenMP::OpenMP_CXX)
set_target_properties(dt PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)

if(MSVC)
    target_compile_definitions(dt PUBLIC FMT_UNICODE=0)
endif()
if(WIN32)
    target_link_libraries(dt PRIVATE psapi)
endif()

if(DT_BUILD_EXECUTABLE)
    add_executable(dt_exec src/dt_main.cpp)
    target_link_libraries(dt_exec PRIVATE dt::dt)
    target_include_directories(dt_exec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extern/cli11)
    set_target_properties(dt_exec PROPERTIES OUTPUT_NAME dt ARCHIVE_OUTPUT_NAME dt_exec)
    if(MSVC)
        target_link_options(dt_exec PRIVATE /STACK:8388608)
    endif()
endif()
