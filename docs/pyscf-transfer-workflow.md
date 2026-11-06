# XYZ → PySCF 输出 → 转移积分

量化计算与转移积分后处理是两个独立阶段。第一阶段读取任意单分子/dimer XYZ，生成各自的 PySCF 输出。第二阶段按 TOML 控制文件读取这些输出，计算转移积分，**不重新运行 SCF**。

## 1. 用独立脚本计算单体和 dimer

CPU 示例（在仓库根目录执行）：

```bash
uv sync --locked --extra quantum
uv run --locked --extra quantum python scripts/run_pyscf.py \
  examples/quantum/ethylene_a.xyz \
  examples/quantum/ethylene_b.xyz \
  examples/quantum/ethylene_ab.xyz \
  --output-dir examples/quantum/results \
  --backend cpu --basis def2-svp --xc b3lyp
```

换成自己的 `A.xyz B.xyz AB.xyz` 即可。也可只给一个 XYZ，分三次独立执行。每个文件按其 basename 生成：

- `.log`：便于阅读的 SCF 计算日志；
- `.chk`：分子、基组、完整 MO、占据、轨道能量及 cube 附加收敛/Fock 数据，供后处理读取；
- `.status.json`：本次成功/失败、设置、能量与实际计算类。

例如 `ethylene_a.xyz → ethylene_a.chk / ethylene_a.log / ethylene_a.status.json`。如果 SCF 不收敛，不保留同名旧成功 checkpoint。批量中某个任务失败时退出 1，已经完成的其他任务保留。脚本不优化、平移或旋转坐标。

`--xc hf` 使用 RHF；否则使用 RKS。辅助基、单位、线程、网格、收敛阈值等可显式设置，见 `python scripts/run_pyscf.py --help`。首版限中性闭壳层、球谐、全电子分子计算。

GPU 量化使用同一脚本：

```bash
uv sync --locked --extra quantum-cuda13
uv run --locked --extra quantum-cuda13 python scripts/run_pyscf.py \
  A.xyz B.xyz AB.xyz --output-dir qc_results \
  --backend gpu --basis def2-svp --xc b3lyp
```

请求 GPU 时调用 GPU4PySCF，不静默回退 CPU。保存的是可在 CPU 后处理的标准 PySCF checkpoint 数据。

## 2. 编写转移积分控制文件

已有完整示例 `examples/quantum/transfer.toml`。例如自己的 `transfer.toml`：

```toml
schema_version = 1

[input]
fragment_a = "qc_results/A.chk"
fragment_b = "qc_results/B.chk"
dimer = "qc_results/AB.chk"

[calculation]
methods = ["dipro"]
carriers = ["hole", "electron"]
backend = "cpu"

[output]
file = "couplings/AB.json"
```

所有相对路径以 **TOML 所在目录**为基准，和执行命令时的目录无关。后处理不需要重填 basis/functional，使用 checkpoint 中保存的基组与轨道；不允许把 SCF 参数写进 `[calculation]`。

默认 dimer 原子顺序为“全部 A 原子，然后全部 B 原子”。若 dimer 原子交错排列，明确提供完整分子的原子映射，顺序必须分别对应 A/B checkpoint 的原子顺序：

```toml
[fragments]
a_atoms = [1, 3, 5, 7, 9, 11]
b_atoms = [2, 4, 6, 8, 10, 12]
```

编号从 **1** 开始，也可写 `a_atoms = "1-6"`。必须覆盖整个 dimer，不能只列用于筛选平面距离的共轭原子。A/B 必须保留它们在 AB 中的实际坐标和基组；单体单独优化、重新居中或旋转后，不能直接用来投影。

默认空穴选择 HOMO 附近、电子选择 LUMO 附近 0.01 eV 内的完整近简并子空间。可用 `[calculation] degeneracy_tolerance_ev` 调整，也可在知道轨道含义时明确指定 **从 1 开始**的 MO 编号：

```toml
[orbitals.hole]
a = [8]
b = [8]

[orbitals.electron]
a = [9]
b = [9]
```

上述 8/9 是乙烯示例，不是所有分子的 HOMO/LUMO 编号。程序检查占据分区，多轨道输出耦合矩阵，不自动压成单个速率参数。JSON 中的 `orbital_indices_*` 使用 Python 从 0 开始的编号，控制文件元数据明确记录两者约定。

## 3. 独立运行后处理

```bash
# 只检查路径和设置，不读取矩阵，不运行量化或探测 GPU
uv run --locked python scripts/calculate_transfer.py examples/quantum/transfer.toml --check

# 读取 A/B/AB checkpoint，计算转移积分
uv run --locked --extra quantum python scripts/calculate_transfer.py examples/quantum/transfer.toml

# 等价的安装命令
uv run --locked --extra quantum cube-transfer examples/quantum/transfer.toml
```

示例输出 `examples/quantum/results/transfer.json`。结果包括空穴/电子耦合的 eV/meV、原始和正交化 Hamiltonian、片段重叠、源文件 SHA256、收敛状态来源，以及 `scf_executed=false`。

后处理设置 `backend="gpu"` 并用 `--extra quantum-cuda13` 运行时，用 CuPy 加速投影和正交化；读取 checkpoint、从完整 MO 恢复特征算符、AO 重叠/偶极积分在 CPU。量化 SCF 是否使用 GPU，由**第一阶段**的 `--backend` 决定，两阶段可以分别选择后端。

## 已有 PySCF 脚本如何保存结果

在你自己的脚本中，在计算结束后加入：

```python
from universal_kmc.quantum.checkpoints import save_transfer_checkpoint

mf.kernel()
save_transfer_checkpoint(mf, "A.chk")
```

这同样适用于 GPU4PySCF 的 `mf`。函数会保存完整轨道、收敛标记和独立 Fock 矩阵；不收敛则拒绝导出。当前模板不接受 QM/MM 或隐式溶剂对象，标准原生 checkpoint 也无法让后处理自动确认原始环境模型。

普通 `.out/.log` 一般缺少完整轨道矩阵，不能代替 `.chk`。原生 PySCF 可通过 `mf.chkfile="A.chk"` 保存标准检查点，但[标准字段](https://pyscf.org/_modules/pyscf/scf/chkfile.html)不含收敛标记。对这种已有文件，须先确认其 SCF 日志收敛，再在控制文件中明确写 `assume_converged=true`。程序会记录 `user_assumed_from_external_log`，不会把这个假设当成独立验证；显式未收敛标记不能被该选项绕过。标准 checkpoint 未保存 functional 等信息时，结果会标明模型一致性无法从文件元数据确认，用户应保证三份文件使用一致方法。

输出模板保存的 Fock 矩阵直接用于 DIPRO，同时用完整 MO 展开作一致性核对。对仅有标准字段的原生 checkpoint，从完整、S 正交的 dimer MO 恢复 `F = S C diag(epsilon) C.T S`，再执行 DIPRO。截断虚轨道、破坏的 MO、坐标/基组不对应等都会被拒绝。恢复恒等式本身不是独立 SCF 收敛测试，结果明确区分保存 Fock 与恢复特征算符两种来源。

此入口支持 `dipro`、受等价性/子空间检查约束的 `esid` 和单电子 `orbital_gmh`；ESID 还须明确 `assume_equivalent=true`。收敛 AB checkpoint 不能定义 frozen-density Hamiltonian，故此入口拒绝 `frozen_density`。各方法的物理定义见[转移积分理论说明](transfer-integrals.md)。目前不在此流程中隐含 QM/MM、带电开壳层或周期量化。

## 验证

[两阶段 CPU 实测报告](../runs/transfer_control_2026-11-07/README.md)包含实际输出和拒绝用例。本次新接口的 GPU 验证尚未执行，不能把之前直接计算接口的 GPU 报告视为新接口的完整验证。可复现命令：

```bash
uv run --locked --extra quantum python scripts/validate_transfer_control.py \
  --output build/transfer-control-validation
# 在有可用 GPU 的机器验证两阶段后端交叉组合：
uv run --locked --extra quantum-cuda13 python scripts/validate_transfer_control.py \
  --output build/transfer-control-gpu-validation --gpu
```
