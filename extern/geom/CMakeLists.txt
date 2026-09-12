# 设置编译的源文件
set(LIBGEOM_SRC 
		intersection_func.cpp
		geom_func.cpp
		geom_func.h
		predicates.hpp
		ttotr.c  )

#AUX_SOURCE_DIRECTORY(. DIR_SRCS)

# 定义宏
add_definitions(
    #-D_LIB
    -D_OPT_INTERSEC_ALGO
    #-D_USE_STATIC_FILTER
    )

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffloat-store")

# 设置生成库
if(WIN32)  # 暂时windows生成静态库
    add_library(geom ${LIBGEOM_SRC})
else()
if(SHARED)
    add_library(geom SHARED ${LIBGEOM_SRC})
else()
    add_library(geom STATIC ${LIBGEOM_SRC})
endif()
endif()

# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
