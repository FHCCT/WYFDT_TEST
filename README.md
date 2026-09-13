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

## 日志等级

`Args::infolevel` / `--infolevel` 默认为 `1`：

| 等级 | 输出 |
| --- | --- |
| 0 | 静默运行，仅输出 error；不输出 warning、步骤或统计信息。 |
| 1 | 主要阶段及耗时；每轮边界边翻转成功数/尝试数、拓扑成功数/候选数、光滑成功数/调用数、Quality 汇总，以及二面角信息。 |
| 2 | 包含 1 级内容，并输出子步骤耗时、迭代、边界恢复、线程配置和质量明细。 |

大于 2 的旧参数按 2 处理。Quality 汇总使用当前优化度量，不再输出步骤前后的 SUS 质量快照。infolevel 为 0 时跳过二面角打印函数及其统计；体积优化停止判据需要的二面角仍独立计算。日志等级不改变迭代停止条件或网格结果，也不修改宿主程序的全局日志等级。计时按阶段汇总，不输出逐点性能日志；默认不生成 CSV 或临时网格文件。

## 依赖

第三方代码保留原版权说明。Eigen 的许可见 `extern/eigen/COPYING*`，spdlog 的许可见 `extern/spdlog/LICENSE`，CLI11 的许可位于 `extern/cli11/CLI11.hpp` 文件头。
