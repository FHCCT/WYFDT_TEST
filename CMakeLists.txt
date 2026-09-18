cmake_minimum_required(VERSION 3.13)
project(dt LANGUAGES C CXX)

# TiGER may also provide these functions before add_subdirectory(wyfdt).
set(_wyfdt_internal_project "${CMAKE_CURRENT_LIST_DIR}/../../cmake/functions/internal_project.cmake")
if(EXISTS "${_wyfdt_internal_project}")
    include("${_wyfdt_internal_project}")
endif()
if(TIGER_MANAGED_BUILD)
    foreach(_wyfdt_command tiger_configure_build_target tiger_link_required_target tiger_enable_openmp)
        if(NOT COMMAND ${_wyfdt_command})
            message(FATAL_ERROR "TiGER managed dt build requires ${_wyfdt_command}; provide internal_project.cmake first")
        endif()
    endforeach()
    option(WYFDT_BUILD_EXECUTABLE "Build dt executable" OFF)
    option(WYFDT_BUILD_LIBRARY "Build libdt library" ON)
    option(WYFDT_BUILD_BW_TESTS "Build mesh regression tests" OFF)
    option(WYFDT_USE_MMG "Enable legacy MMG integration (unavailable in this source tree)" OFF)
    set(_wyfdt_build_executable ${WYFDT_BUILD_EXECUTABLE})
    set(_wyfdt_build_library ${WYFDT_BUILD_LIBRARY})
    set(_wyfdt_build_tests ${WYFDT_BUILD_BW_TESTS})
    set(_wyfdt_use_mmg ${WYFDT_USE_MMG})
else()
    option(BUILD_EXECUTABLE "Build dt executable" ON)
    option(BUILD_LIBRARY "Build libdt library" ON)
    option(BUILD_BW_TESTS "Build mesh regression tests" OFF)
    option(USE_MMG "Enable legacy MMG integration (unavailable in this source tree)" OFF)
    set(_wyfdt_build_executable ${BUILD_EXECUTABLE})
    set(_wyfdt_build_library ${BUILD_LIBRARY})
    set(_wyfdt_build_tests ${BUILD_BW_TESTS})
    set(_wyfdt_use_mmg ${USE_MMG})
endif()
if(_wyfdt_use_mmg)
    message(FATAL_ERROR "This dt source tree has no MMG adapter or extern/mmg; disable USE_MMG/WYFDT_USE_MMG")
endif()

# Production uses OpenMP, including the one-thread path.
if(NOT TIGER_MANAGED_BUILD)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
endif()
set(DT_SOURCES
    src/dt.cpp src/dt_adp.cpp src/dt_bw.cpp src/dt_bw_parallel.cpp
    src/dt_func.cpp src/dt_omp.cpp src/dt_opt.cpp src/dt_ort.cpp
    src/dt_parallel.cpp src/dt_sus.cpp src/dt_surface.cpp src/mesh_io.cpp
    extern/geom/intersection_func.cpp extern/geom/predicates.cpp extern/geom/ttotr.c)

function(dt_configure_target target)
    if(TIGER_MANAGED_BUILD)
        tiger_configure_build_target(${target})
    endif()
    target_include_directories(${target} PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/extern/geom
        ${PROJECT_SOURCE_DIR}/extern/tree)
    if(TIGER_MANAGED_BUILD)
        tiger_link_required_target(${target} PUBLIC TiGER::cli11)
        tiger_link_required_target(${target} PUBLIC TiGER::eigen)
        tiger_link_required_target(${target} PUBLIC TiGER::spdlog)
        tiger_enable_openmp(${target})
    else()
        target_include_directories(${target} PUBLIC
            ${PROJECT_SOURCE_DIR}/extern/cli11
            ${PROJECT_SOURCE_DIR}/extern/eigen
            ${PROJECT_SOURCE_DIR}/extern/spdlog/include)
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
    endif()
    # Set after TiGER configuration so its default standard cannot replace 11.
    target_compile_features(${target} PUBLIC cxx_std_11)
    set_target_properties(${target} PROPERTIES
        CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
    if(MSVC)
        # Preserve the existing Windows code-page encoding with bundled fmt.
        target_compile_definitions(${target} PUBLIC FMT_UNICODE=0)
        target_compile_options(${target} PRIVATE /W0)
    else()
        target_compile_options(${target} PRIVATE -w)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            # Reject newer-language extensions even when other warnings are muted.
            target_compile_options(${target} PRIVATE
                "$<$<COMPILE_LANGUAGE:CXX>:-pedantic-errors>")
        endif()
    endif()
    if(WIN32)
        target_link_libraries(${target} PRIVATE psapi)
    endif()
endfunction()

if(_wyfdt_build_library)
    add_library(dt STATIC ${DT_SOURCES} src/dt_main.cpp)
    dt_configure_target(dt)
    target_compile_definitions(dt PUBLIC DT_LIBRARY)
    set_target_properties(dt PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
if(_wyfdt_build_executable)
    add_executable(dt_exec src/dt_main.cpp)
    dt_configure_target(dt_exec)
    target_compile_definitions(dt_exec PRIVATE DT_EXEC)
    if(_wyfdt_build_library)
        # Reuse the library's algorithms without a third intermediate target.
        target_link_libraries(dt_exec PRIVATE dt)
    else()
        target_sources(dt_exec PRIVATE ${DT_SOURCES})
    endif()
    # Keep the executable's import-library output distinct from the static library.
    # MSBuild otherwise treats dt.lib as its own output and may skip relinking.
    set_target_properties(dt_exec PROPERTIES OUTPUT_NAME dt ARCHIVE_OUTPUT_NAME dt_exec)
    if(MSVC)
        target_link_options(dt_exec PRIVATE /STACK:8388608)
    endif()
endif()
if(_wyfdt_build_tests AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    if(NOT _wyfdt_build_library)
        message(FATAL_ERROR "BUILD_BW_TESTS requires BUILD_LIBRARY")
    endif()
    enable_testing()
    add_subdirectory(tests)
endif()
