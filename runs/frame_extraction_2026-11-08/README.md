# 多帧单体、dimer 与冻结 MM 簇验证

2026-11-08，使用本机实际 `anneal_a.tpr`、`anneal_a_fix.gro`、`anneal_a.xtc` 与 `core.txt` 验证新增 `$get_dimers_monomers`。仅提取结构，未启动真实材料的新量化任务。

设置：起始帧 1980（从 0 开始）、frame_step=20、dimers_cutoff=6 Å、monomers-qmmm-cutoff 留空（12 Å）、format=xyz,gjf。轨迹共 2001 帧，实际提取 1980、2000；每帧 100 个完整分子，每个 86 原子。

| 帧 | 单体 | dimer | 单体簇 | XYZ 文件 | GJF 文件 | MM 点电荷文件 |
|---|---:|---:|---:|---:|---:|---:|
| 1980 | 100 | 552 | 100 | 752 | 752 | 100 |
| 2000 | 100 | 550 | 100 | 750 | 750 | 100 |

- 帧集合索引能通过各帧 manifest 的 SHA256 校验。
- 第 2000 帧全部 550 个 dimer 的 ID、core 最短距离及完整坐标，与旧提取逐项一致（坐标绝对容差 1e-12 Å）。
- 第 2000 帧 1 号分子的簇包含 86 个 QM 原子、5160 个 MM 原子，MM 来自 60 个完整分子镜像；MM 坐标和电荷与之前 CPU/GPU 四点法基准使用的原环境逐项完全相同。
- 同一输入 `--resume` 成功，200 个簇 XYZ 的修改时间均保持不变。
- 新增测试还验证了多帧实际 CPU HF/STO-3G 单体、dimer SCF 与 DIPRO：各帧独立输出、改变分子间距导致耦合改变、续算不重写 checkpoint。该小体系验证调度与数值流程，不代表材料量化结果。
- 全部 59 项测试中 58 项通过，1 项需要显式 GPU 的 CuPy 检查跳过；包含四语言一致性与冻结 MM/LJ 的既有回归。新增工作流测试共 15 项通过。
- `uv --cache-dir /tmp/cube-uv-cache build` 成功构建 sdist 与 wheel。

完整原始结构、簇和参数均仅保存在忽略的 `build/frame-extraction-20261108/` 中，不提交材料坐标。可在同一数据机器上使用该目录中的 `input` 复查：

```bash
./cube build/frame-extraction-20261108/input --only get_dimers_monomers --resume
```

`*_cluster.xyz` 是完整 QM+MM 簇几何；`*_cluster.gjf` 是中央 QM 与外部冻结点电荷的静电单点模板，未运行 Gaussian，不含 LJ 或四点法优化。实际四点法仍通过 cube 的 PySCF/GPU4PySCF 流程计算，环境选择共用同一函数。
