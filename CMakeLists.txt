cmake_minimum_required(VERSION 3.15)
project(dt LANGUAGES C CXX)

option(BUILD_EXECUTABLE "Build dt executable" ON)
option(BUILD_LIBRARY "Build libdt library" ON)
option(BUILD_BW_TESTS "Build mesh regression tests" OFF)

# Production uses OpenMP, including the one-thread path.
find_package(OpenMP REQUIRED COMPONENTS CXX)
set(DT_SOURCES
    src/dt.cpp src/dt_adp.cpp src/dt_bw.cpp src/dt_bw_parallel.cpp
    src/dt_func.cpp src/dt_omp.cpp src/dt_opt.cpp src/dt_ort.cpp
    src/dt_parallel.cpp src/dt_sus.cpp src/mesh_io.cpp
    extern/geom/intersection_func.cpp extern/geom/predicates.cpp extern/geom/ttotr.c)

function(dt_configure_target target)
    target_include_directories(${target} PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/extern/geom
        ${PROJECT_SOURCE_DIR}/extern/cli11
        ${PROJECT_SOURCE_DIR}/extern/eigen
        ${PROJECT_SOURCE_DIR}/extern/spdlog/include
        ${PROJECT_SOURCE_DIR}/extern/tree)
    target_compile_features(${target} PRIVATE cxx_std_14)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W0)
    else()
        target_compile_options(${target} PRIVATE -w)
    endif()
    if(WIN32)
        target_link_libraries(${target} PRIVATE psapi)
    endif()
endfunction()

if(BUILD_LIBRARY OR BUILD_EXECUTABLE)
    # Compile algorithms once; only the API/CLI entry differs.
    add_library(dt_core OBJECT ${DT_SOURCES})
    dt_configure_target(dt_core)
    set_target_properties(dt_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
if(BUILD_LIBRARY)
    add_library(dt STATIC src/dt_main.cpp $<TARGET_OBJECTS:dt_core>)
    dt_configure_target(dt)
    target_compile_definitions(dt PUBLIC DT_LIBRARY)
    set_target_properties(dt PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
if(BUILD_EXECUTABLE)
    add_executable(dt_exec src/dt_main.cpp $<TARGET_OBJECTS:dt_core>)
    dt_configure_target(dt_exec)
    target_compile_definitions(dt_exec PRIVATE DT_EXEC)
    set_target_properties(dt_exec PROPERTIES OUTPUT_NAME dt)
    if(MSVC)
        target_link_options(dt_exec PRIVATE /STACK:8388608)
    endif()
endif()
if(BUILD_BW_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    if(NOT BUILD_LIBRARY)
        message(FATAL_ERROR "BUILD_BW_TESTS requires BUILD_LIBRARY")
    endif()
    enable_testing()
    add_subdirectory(tests)
endif()
