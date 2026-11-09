# 模型 Surface hopping 使用说明

这是 `feat/surface-hopping` 分支的第一阶段实现，面向下一步有机材料输运接口。现在可以运行 Tully 1/2/3 和固定正交基的多态 LVC 模型；没有接入真实 MD 的非绝热量化、GPU、SOC、材料迁移率或激子扩散长度。理论与开源项目比较见 [调研](research/surface-hopping.md)。

## 安装与入口

```bash
uv sync --extra dynamics
uv run cube examples/surface-hopping/input --plan
uv run cube examples/surface-hopping/input
uv run cube examples/surface-hopping/input --resume
```

示例的三个顶层入口位于文件最前，只开启 `$surface-hopping`；下面的模块为：

```text
$surface-hopping-dynamics
input = tully1.toml
output = ../../build/surface-hopping/tully1
$end
```

`input` 是 TOML 控制文件，路径相对于当前 section 控制文件；模块 `output` 覆盖 TOML 的输出目录。若仅使用模型，无需添加结构提取段。完整入口缺少该模块时会在启动 MD/QC 前报错；`do=0` 跳过它。`--only` 仍受入口门控。若显式加入共享结构提取段，它会按既有规则提取，但模型不会自动使用该结构。

也可直接运行专用脚本：

```bash
uv run cube-surface-hopping examples/surface-hopping/tully1.toml --check
uv run cube-surface-hopping examples/surface-hopping/tully1.toml
uv run cube-surface-hopping examples/surface-hopping/lvc.toml
```

## 单位和初始条件

所有数值使用原子单位：坐标 Bohr、能量 Hartree、时间 `hbar/Eh`（约 0.02418884 fs）、质量为电子质量，动量为相应原子单位。没有自动从 Å、amu 或 fs 转换。原子笛卡尔 provider 应把每个原子的质量重复到三个坐标；质量加权／正常模坐标则必须提供与之相符的坐标质量，不能混用。

`initial.active-state` 是从 0 开始、按能量升序排列的 **adiabatic** 态。控制文件当前使用确定的初始坐标、动量和纯活跃态，轨迹间随机种子独立；温度／Wigner 分布没有自动生成。Python API 可以传入归一化的初始电子振幅。散射验证脚本另行显式使用匹配量子波包的 Wigner 采样。

TOML 示例包括模型、初始条件、核时间步、电子子步、轨迹数、步数、种子、输出间隔及可选散射边界。`decoherence="none"` 是基础 FSSH；`"edc"` 使用 Granucci–Persico Eq.17，`edc-constant-hartree=0.1`。必须比较两种模型选择及其收敛，不能默认认为打开修正即可得到真实材料结果。

## LVC 与 provider

LVC JSON 格式详见 [示例](../examples/surface-hopping/lvc.json)，其参数全部是**合成测试数据**：

```text
H(q) = H0 + sum_k linear[:,:,k] * q[k] + 0.5 * sum_k harmonic[k] * q[k]^2 * I
```

JSON 的 `schema="cube-lvc-1"`、`units="atomic"`，并包含 `h0[n,n]`、`linear[n,n,dof]`、`harmonic[dof]`。不直接从单帧转移积分拟合导数。一般 provider：

```python
from universal_kmc.surface_hopping.engine import FSSH
# provider.evaluate(q) -> H[n,n], dH[n,n,dof], real symmetric in a fixed orthonormal basis
engine = FSSH(provider, coordinate_masses, dt=0.5, electronic_substeps=4)
state = engine.initialize(coordinates, momenta, active=0)
state, event = engine.step(state, rng)
```

当前每个核时间步的末尾最多一次 hop，电子子步只改善传播和概率流积分。若概率超过 1、态重叠过小或出现近简并，则要求减小时间步或更换支持局域 diabatic 的方法；没有静默归一化概率。Frustrated hop 保留动量，不反向速度。

## 输出与验证

`result.json` 包含运行范围、控制与模型哈希、NumPy/RNG 信息、最终活跃态频率及标准误、平均电子占据、hop 计数、能量漂移、归一化误差和逐轨迹终态。`trajectories.jsonl` 记录坐标、动量、时间、活跃态、两种表象的电子占据与 diabatic IPR。活跃态分布与电子振幅占据是不同估计，不能相互替代。`max_steps` 是终止原因，散射未退出的轨迹使结果标为 `partial`（退出码 2）。无散射边界的定长模型运行完成可标为 `completed`，这不代表统计或时间步已收敛。

`--resume` 只复用同一控制／模型哈希且轨迹文件哈希一致的完整结果。失败／部分结果不继续拼接，需选择新输出目录；当前不恢复中断轨迹。现有成功结果不允许未经 `--resume` 覆盖。

```bash
uv run python -m unittest discover -s tests -p test_surface_hopping.py -v
uv run python scripts/validate_surface_hopping.py --trajectories 128
```

第二条命令用独立 FFT split-operator 解 Tully 1 量子散射，并比较具有同一 Wigner 初始分布的 FSSH 轨迹。量子参考做网格／时间步加密，随机结果提供二项统计误差。单个高动量工况的吻合不能证明所有 Tully 模型和凝聚相输运都准确；低动量、反射、多次穿越、去相干和真实材料仍需后续系统验证。
