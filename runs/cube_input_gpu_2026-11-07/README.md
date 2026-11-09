# cube input：RTX 5090 工作流验证

2026-11-07 在独立服务器上实际安装提交 `3a50c65` 构建的 wheel，运行 `scripts/validate_cube_workflow.py --gpu`，**完整通过**。本次没有修改计算实现；验证的是新增 `$one-mol/$two-mol/$reorganization-energy` 调度、GPU 能量/梯度以及 vacuum 四点法。未上传完整 MD 轨迹，未运行真实 MD 的 100 个单体、550 个 dimer 全量量化。

## 环境与方法

- RTX 5090，32607 MiB，compute capability 12.0，驱动 615.71.09；CPU 为 Xeon Platinum 8380，SCF 使用 4 线程。
- Python 3.12.15、PySCF 2.15.0、GPU4PySCF-cuda13x 1.9.1、CuPy-cuda13x 13.6.0、NumPy 2.5.4、geomeTRIC 1.1.1。运行时 CUDA 13.0，驱动 API 13.4。
- B3LYP/6-31G(d,p)，def2-universal-jkfit 密度拟合，grid 3，SCF 收敛阈值 1e-10。CPU/GPU 使用相同设置和起始几何。
- 两个乙烯单体与一个 dimer 用于 DIPRO；F2 用于中性/阳离子/阴离子几何优化和空穴/电子四点法。F2 是小体系接口与数值验收案例，不用于判断此泛函、基组的实验精度。

版本、安装包和脚本 SHA256 见 [environment.json](environment.json)。测试安装包与此前本地 CPU 报告的 wheel 哈希一致。GPU4PySCF 提示使用 CuPy 张量收缩引擎；没有静默切换到 CPU 量化。

## 数值对照

| 量 / eV | CPU | GPU | 绝对差 / eV |
|---|---:|---:|---:|
| 乙烯空穴 DIPRO 耦合幅值 | 0.421764473455 | 0.421764326268 | 1.472e-7 |
| 乙烯电子 DIPRO 耦合幅值 | 0.474021357909 | 0.474021212677 | 1.452e-7 |
| F2 空穴四点重组能 | 0.405514348565 | 0.405514070297 | 2.783e-7 |
| F2 电子四点重组能 | 5.706538099647 | 5.706537728674 | 3.710e-7 |

这是数值容差内的一致性，不要求逐位相同。预设耦合容差为 atol=2e-5 eV、rtol=2e-4；重组能绝对差阈值为 5e-3 eV，实际差异远小于这些阈值。结果及原始四点 Hartree 能量见 [summary.json](summary.json)、[CPU 重组能](cpu/F2-reorganization.json)、[GPU 重组能](gpu/F2-reorganization.json)；耦合与后处理信息见 [CPU DIPRO](cpu/ethylene-transfer.json)、[GPU DIPRO](gpu/ethylene-transfer.json)。

重组能明确为 **intramolecular vacuum**：优化中性 R0 与带电 Rq，再独立计算 E0(R0)、Eq(R0)、Eq(Rq)、E0(Rq)。未加入 QM/MM、溶剂或外球贡献，也未做频率分析来证明全局极小值。开壳层、自旋及各步骤 SCF 信息保存在原始 JSON。

## GPU 执行、梯度与可靠性检查

- 中性和两个带电优化全部收敛。GPU 优化每一步通过运行时检查，确认 MO 是 CuPy 数组；trace 记录 GPU 后端和实际 GPU 设备。梯度扫描器为 `gpu4pyscf.grad.rhf.Gradients_Scanner`。优化器本身在 CPU，能量与梯度由 GPU4PySCF 计算。
- 在同一 F2 起始几何，对 q=0、+1、−1 分别用真实中心能量差分核对解析梯度。位移 1e-3 Bohr，GPU 最大绝对差 **1.345e-5 Hartree/Bohr**，CPU 最大差 **2.087e-5 Hartree/Bohr**，均通过预设 1e-4 阈值。保留有限网格/SCF 误差和日志中的简并轨道警告，不将门限解释为通用材料精度保证。
- 单体分别居中，经过已校验的纯平移 checkpoint 映射进入共同 dimer 坐标；DIPRO 后处理读取已完成的 SCF，不重新执行 SCF。
- CPU/GPU 各再次 `--resume`：checkpoint、优化 trace 和 transfer JSON 修改时间均不变，任务报告标记复用。
- 人为限制优化为一步，使实际优化失败；状态为 failed，没有成功的重组能，见 [失败状态](failure/reorganization.json)。
- 安装后的 `.venv/bin/cube .../input --plan` 与 `--resume` 也实际通过，日志见 [入口预检](installed-entry-plan.log)、[入口续算](installed-entry-resume.log)。

完整验证输出见 [validation.log](validation.log)。归档日志仅去除了 ANSI 颜色码和行尾空白；原始日志保留在服务器及本地 build 副本。这里只提交小分子输入、数值、优化 trace 和日志，不提交 wheel、checkpoint 或大型轨迹。归档 JSON 中的绝对路径是服务器原始运行位置，用于保留计算来源；服务器保留了全部运行产物。本地完整取回副本在 `build/cube-workflow-rtx5090/`。

## 耗时与复现

完整小体系调度单次观察：CPU **20.194 s**、GPU **20.528 s**；其中包括多个小 SCF、进程初始化、优化、checkpoint 与结果写入，不包括之后独立的梯度差分检查。续算分别为 0.489 s 与 0.486 s。该案例没有体现 GPU 加速收益，不能作为实际大分子或批量计算的加速比。未进行重复计时或大体系性能基准。

在兼容 GPU 环境中，从仓库根目录运行并使用新的输出目录：

```bash
uv sync --locked --extra workflow --extra quantum-cuda13
uv run --locked --extra workflow --extra quantum-cuda13 python scripts/validate_cube_workflow.py \
  --output build/cube-workflow-gpu-new --gpu
```

真实 MD 提取已在本地对参考脚本验证；本次 GPU 验证从小分子 manifest 开始，不重复声称覆盖完整 MD→全量量化流程。GTX 1650 上 CUDA 13 wheel 的不兼容记录仍见 [此前报告](../cube_input_2026-11-07/README.md)。
