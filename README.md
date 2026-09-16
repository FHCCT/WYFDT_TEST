# WYFDT_TEST

DT 三维四面体网格生成与优化库，提供 C++ 库和命令行程序。

## 构建

需要 CMake 3.15 或更高版本、支持 C++11 的编译器及 OpenMP。所需 Eigen、spdlog、CLI11 头文件随源码提供。

```powershell
cmake -S . -B build
cmake --build build --config Release --parallel 4
```

构建目标使用 C++11，并禁用编译器语言扩展；链接 `dt` 的 CMake 项目会继承 C++11 的最低语言要求。MSVC 没有严格的 `/std:c++11` 开关，其最低选项是 `/std:c++14`；严格兼容性使用 GCC 的 `-std=c++11 -pedantic-errors` 编译验证。

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

## 各向同性自适应质量优化

`--adptype 2` 使用统一默认参数：`--optloop 18 --optjacob 0.8 --optanglestrict 0 --adpangle 160`。显式传入的这四个选项优先，不再被入口覆盖。这里 `optjacob` 是当前阶段质量度量的候选阈值，不是最终平均雅可比的保证值；`optanglestrict=0` 仅关闭优化阶段的插点，按输入单元 ID 进行的自适应加密仍照常执行。

所有常规优化入口统一使用 `flipEdgPass → TopologicalPass → SmoothPass`，按各流程原有约束决定是否执行边界翻转。`TopologicalPass` 直接收集并稳定排序候选，再执行翻转及失败修复；不再设置中间转发 Pass。`QuicklyOptPass` 在默认 18 轮预算下，前 12 轮使用 SUS 质量、最后 6 轮使用最小二面角质量；内部点在 SUS 阶段使用 SUS 能量，在角度阶段使用最小二面角活跃集。SUS 和角度阶段各至少一轮：总预算最低为 2，角度轮数为 max(1, floor(总轮数 / 3))，余下分配给 SUS，不再限制角度阶段最多 6 轮。SUS 阶段无差单元时直接进入角度阶段，角度阶段无差单元时才整体退出。不再使用雅可比组合度量。退出使用轮数预算或无候选差单元，插点角度阈值不作为整体质量的收敛判据。

`SmoothPass` 每轮统一选择并染色边界点、内部点；染色覆盖所有相邻物理四面体（包括好单元），同色点并行、不同颜色依次执行。在同一个循环中，边界点调用切向光滑，内部点调用 `smoothInteriorPoint`，其中 `improve_Metric == 2` 分派到 `smooth_angle`，其余度量保留原来的 `smooth_sus` 路径，临时边界缓冲区按线程复用。`modifyBnd == false` 时染色不纳入边界点，边界光滑函数自身也直接拒绝移动；边界只在允许修改边界的 SUS/角度流程中参与；平面/直线特征移动必须位于全部相邻边界三角形的平面内，跳过角点、曲面点、周期点及锁定点/边/面。

每轮光滑开始时固定全局最低 SUS 质量作为共同下限，允许较好局部单元降低质量以改善整体；不承诺每个点邻域的最低质量不变，也不承诺整个拓扑阶段的 SUS 最低值不变。单独调用 `smoothInteriorPoint(node)`、`smooth_sus(node)` 或 `smooth_angle(node)` 仍默认保护局部最低 SUS 质量；可通过第二个参数显式传入质量下限。光滑完成后缓存仍使用当前阶段度量，角度阶段不会混入 SUS 质量缓存。

角度活跃集以各相邻四面体的六个二面角为独立约束，收集接近最小角的解析梯度，通过梯度凸包的最小范数点寻找共同上升方向，并回溯步长。每次接受必须提高邻域最小二面角，同时保持正体积和 SUS 下限，并保护本次点光滑开始时的邻域平均最小二面角、平均 SUS 质量，避免只改善最差角却损伤大量较好单元；这些局部保护不等同于整个拓扑流程的最终平均质量保证。优化插点后的光滑也经过 `smoothInteriorPoint`，随当前度量分派。边界恢复/Steiner 移除和体积均匀化所用的专门光滑算法保留原有层级。

库调用 `adaptation_by_Tetid` 时设置同一组 `Args` 即可复现该配置：`optloop=18, optTh=0.8, optanglestrict=0, adpangle=160`。库保留调用方的参数，不按模型名称选择策略。提高轮数会增加耗时；以上配置侧重最终质量，具体平均角度与雅可比仍取决于输入网格和边界约束。

## 日志等级

`Args::infolevel` / `--infolevel` 默认为 `1`：

| 等级 | 输出 |
| --- | --- |
| 0 | 静默运行，仅输出 error；不输出 warning、步骤或统计信息。 |
| 1 | 主要阶段及耗时；每轮边界边翻转成功数/尝试数、拓扑成功数/候选数、光滑成功数/调用数、Quality 汇总，以及二面角信息。 |
| 2 | 包含 1 级内容，并输出子步骤耗时、迭代、边界恢复、线程配置和质量明细。 |

大于 2 的旧参数按 2 处理。Quality 汇总使用当前优化度量，不再输出步骤前后的 SUS 质量快照。infolevel 为 0 时跳过二面角打印函数及其统计；体积优化停止判据需要的二面角仍独立计算。日志等级不改变迭代停止条件或网格结果，也不修改宿主程序的全局日志等级。计时按阶段汇总，不输出逐点性能日志；默认不生成 CSV 或临时网格文件。

`dt_init` 在 `infolevel > 0` 时汇总 `constrain`、`refine`、线程数及优化参数；`refine == 1` 时另列 `size`、`minEdge`、`maxEdge`、`growsize`。线程数区分配置上限、用户请求值和按输入规模选出的初始线程数，后续阶段仍会动态调整。

初始化同时输出输入点、线段、三角面和四面体数量。生成阶段在细化结束后汇总网格量与速度，计时涵盖边界点插入、边界恢复、区域分类及细化；优化阶段输出四面体数 `A->B` 和速度；`adaptation_by_YunBoSzieControl` 在导出网格后汇总最终四面体数、总耗时和速度。速度统一为 `W/s`（万个四面体每秒），按阶段结束四面体数除以耗时再除以 10000 计算。内部计数排除删除、虚拟和包围单元；自适应结束计数使用实际导出网格。文件写出只显示路径，不再计时。

## 依赖

第三方代码保留原版权说明。Eigen 的许可见 `extern/eigen/COPYING*`，spdlog 的许可见 `extern/spdlog/LICENSE`，CLI11 的许可位于 `extern/cli11/CLI11.hpp` 文件头。

日志依赖使用 [spdlog 1.17.0](https://github.com/gabime/spdlog/releases/tag/v1.17.0) 及该版本配套的 fmt 12.1.0，位于 `extern/spdlog`。仅保留官方完整头文件和许可证，不引入示例、测试、源码构建目录或额外库目标。fmt 许可证位于 `include/spdlog/fmt/bundled/fmt.license.rst`。MSVC 构建使用 `FMT_UNICODE=0` 兼容现有 Windows 源码编码；通过 CMake 链接 `dt` 时自动继承该配置。
