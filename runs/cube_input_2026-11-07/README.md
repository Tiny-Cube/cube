# cube input 入口：实现与本轮验证

实现 `cube input` 的 `$get_dimer/$one-mol/$two-mol/$reorganization-energy ... $end` 分段控制。仓库根目录有可执行 `cube`、已填入实际数据路径的 `input`；安装包提供同名命令。量化和任务顺序按依赖调度，原有独立量化/转移积分接口保持可用。

## 真实 MD 提取

实际文件为 `anneal_a_fix.gro / anneal_a.tpr / anneal_a.xtc / core.txt`。轨迹 2001 帧，最后一帧编号 **2000**（零基）、时间 **20000 ps**，盒子约 45.76794 Å，100 个 86 原子中性分子。core 为用户给的 32 个局部原子编号。

以严格最短 core 原子距离 `< 6 Å`，得到 **100 个完整单体和 550 个 dimer**，550 个不同分子对，无重复镜像冲突。对参考脚本 `extract_dimers.py` 的 `contact_images()` 使用完全相同的已还原坐标和盒子交叉核对，分子对/镜像集合相同，最大距离差 **0 Å**。此检查比较距离与镜像算法，不声称输出文件格式/脚本的 Gaussian/Slurm 设置逐字一致。

`format=xyz,gjf` 写出 650 个 XYZ 和 650 个 GJF，另有距离 CSV、manifest 和 workflow 状态。GJF route 使用用户指定的 `b3lyp/6-31G(d,p)`。`format=gjf` 不写 XYZ 的情况通过独立测试。完整输出留在 `mol/frame2000/`；不提交原始轨迹、真实 MD 几何或大批结构文件。源数据和参考脚本哈希见 `extraction-summary.json`。

## CPU 调度、DIPRO 和 vacuum 四点法

`scripts/validate_cube_workflow.py` 实际从 input 控制调度：

1. 两个乙烯单体分别居中，运行单点 SCF。
2. 单体 checkpoint 经验证的纯平移映射到共同 dimer 坐标，二聚体 SCF 后读取已有结果进行 DIPRO，后处理不重新 SCF。
3. 对 F2 中性/阳离子/阴离子实际优化并计算四点重组能。F2 用作小体系接口/数值验收，不用来评价此泛函/基组的实验精度。
4. 使用 `--resume` 再次运行，checkpoint、优化 trace 和 transfer JSON 的修改时间均保持不变；另验证强制优化不收敛时没有成功的重组能。

统一 B3LYP/6-31G(d,p)、def2-universal-jkfit 密度拟合、grid 3、SCF conv_tol 1e-10。记录结果为：

| 项目 | 空穴 / eV | 电子 / eV |
|---|---:|---:|
| 乙烯 DIPRO 耦合幅值 | 0.421764473455 | 0.474021357909 |
| F2 vacuum 四点重组能 | 0.405515314669 | 5.706538025021 |

中性与带电优化均收敛，保存完整四点 Hartree 能量、两项松弛能、优化坐标、梯度和自旋信息。明确为 **intramolecular vacuum**，没有 QM/MM、溶剂或外球贡献，也未以频率分析证明全局极小值。实际能量来源见 `cpu-F2-reorganization.json`。

另对中性、阳离子和阴离子在同一 F2 几何核对解析梯度与真实中心能量差分，步长 1e-3 Bohr；最大绝对差约 **3.78e-5 Hartree/Bohr**，通过预设 1e-4 门限。该测试保留有限网格/SCF 数值误差，不把一致性门限称为任意材料的精度保证。原始数值见 `finite-difference-gradients.json`。

## GPU 状态与限制

本轮本机检测到 GTX 1650 Max-Q、4 GB、驱动 615.71.09。安装锁定的 GPU 可选依赖后，CuPy 设备初始化/数组计算成功；**GPU4PySCF 1.9.1 的 CUDA 13 wheel 实际运行报 `no kernel image is available for execution on the device`**，失败在 `block_diag` 内核。因此这台机器没有完成新调度器的 GPU SCF/梯度/四点法验证，程序也没有回退 CPU 或产生伪成功 checkpoint。日志与失败状态为 `local-gpu-failure.log/json`。

此前 RTX 5090 的直接计算与 checkpoint 控制接口已经通过验证，见 `../transfer_control_gpu_2026-11-07/`；这不替代本轮新增四点法 GPU 验证。远程测试因自动审批要求明确允许上传新代码制品而暂未执行。此处不把本机 CuPy 测试成功声称为 GPU4PySCF 全流程通过。

## 测试和使用

- **45 项 Python 测试全部通过，无跳过**；其中 GPU 测试是 CuPy Löwdin 代数测试。此项不证明 GPU4PySCF 内核兼容本机。
- 11 项新增专项测试覆盖严格解析/引擎继承/路径/依赖顺序、周期边界与三斜盒/多镜像、格式选择与哈希复用、checkpoint 纯平移与拒绝旋转、四点公式符号和优化失败。
- C++、Fortran、Rust 及四语言对照测试通过，uv wheel/sdist 构建通过。入口可执行权限和安装命令均已实际检查。

```bash
uv sync --locked --extra workflow
./cube input --plan
./cube input --only get_dimer --resume
uv run --locked --extra workflow python scripts/validate_cube_workflow.py --output build/cube-workflow-cpu-new
```

实际数据的 100 个单体、550 个 dimer 和 100 个重组能任务没有作为快速验证全部启动。用户可通过 `molecules`/`dimers` 选择任务，或执行完整 `./cube input --resume`；GPU 模式需要兼容的运行环境。详见 `docs/cube-input.md`。
