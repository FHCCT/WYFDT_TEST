# WYFDT_TEST

DT 三维四面体网格生成与优化库，提供 C++ 库和命令行程序。

## 构建

需要 CMake 3.15 或更高版本、支持 C++14 的编译器及 OpenMP。所需 Eigen、spdlog、CLI11 头文件随源码提供。

```powershell
cmake -S . -B build
cmake --build build --config Release --parallel 4
```

Windows / Visual Studio 构建生成 `build/Release/dt.exe` 和 `build/Release/dt.lib`。

仅构建库：

```powershell
cmake -S . -B build -DBUILD_EXECUTABLE=OFF -DBUILD_LIBRARY=ON
cmake --build build --config Release --parallel 4
```

## 运行

```powershell
./build/Release/dt.exe --help
./build/Release/dt.exe --input mesh.vtk --constrain 0 --nthread 4
```

对外 C++ 接口见 `src/dt_API.h`。CMake 库目标为 `dt`，命令行目标为 `dt_exec`。

## 依赖

第三方代码保留原版权说明。Eigen 的许可见 `extern/eigen/COPYING*`，spdlog 的许可见 `extern/spdlog/LICENSE`，CLI11 的许可位于 `extern/cli11/CLI11.hpp` 文件头。
