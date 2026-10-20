# DT

三维四面体网格生成、优化与自适应库。需要 CMake 3.13+、C++11 和 OpenMP。

## 集成到其他项目

将 `CMakeLists.txt`、`src/`、`extern/` 一起放入宿主项目，例如 `third_party/dt/`：

```cmake
add_subdirectory(third_party/dt)
target_link_libraries(your_target PRIVATE dt::dt)
```

默认只构建静态库 `dt`。头文件路径、接口编译宏和 OpenMP 链接依赖会自动传递，无需宿主手动配置。不读取宿主特定的 CMake 脚本，也不修改宿主的全局编译选项或运行库设置。

对外接口见 `src/dt_API.h`，通过 `#include <dt_API.h>` 使用；需要内部网格对象时可以包含 `dt.h`。API 实现在 `src/dt_api.cpp`，命令行入口独立放在 `src/dt_main.cpp`。

## 单独构建静态库

```sh
cmake -S . -B build
cmake --build build --config Release --parallel 4
```

Windows / Visual Studio 输出 `build/Release/dt.lib`；其他生成器的输出位置遵循 CMake 配置。

## 可选命令行程序

唯一的项目开关为 `DT_BUILD_EXECUTABLE`，默认 `OFF`：

```sh
cmake -S . -B build -DDT_BUILD_EXECUTABLE=ON
cmake --build build --config Release --target dt_exec --parallel 4
```

命令行目标名为 `dt_exec`，输出程序名为 `dt`。Windows 下可运行：

```sh
build/Release/dt.exe --help
build/Release/dt.exe --input mesh.vtk --output result.vtk
```

## 随附依赖

- Eigen：`extern/eigen`，许可见其中的 `COPYING*`。
- spdlog / fmt：`extern/spdlog`，许可见 `LICENSE` 和 `include/spdlog/fmt/bundled/fmt.license.rst`。
- 几何谓词与空间树：`extern/geom`、`extern/tree`。
- CLI11：仅命令行程序使用，许可见 `extern/cli11/CLI11.hpp` 文件头。

保留第三方版权和许可文件。测试、基准、实验数据和历史报告不参与库构建，也无需复制到宿主项目。
