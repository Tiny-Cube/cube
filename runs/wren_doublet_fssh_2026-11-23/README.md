# 有限乙烯二聚体阳离子的空穴轨迹

两个乙烯分子的平面间距为 3.6 Å，面内错位 0.4 Å 和 0.2 Å。总电荷 +1，自旋多重度 2，共 31 个电子。十四个空间轨道构成双占据核心，三个活性电子分布在两个前线轨道中。每一步做 CPU ROHF、两个双重态的等权 SA-CASSCF(3,2)/STO-3G，计算解析梯度、态间耦合和片段跃迁电荷矩阵。

```bash
uv run --extra dynamics --extra quantum wren examples/surface-hopping/ethene-hole.input
uv run --extra dynamics --extra quantum python runs/wren_doublet_fssh_2026-11-23/validate.py build/ethene-hole-validation
```

`validation.json` 中记录以下核验。

- 在二聚体的一个实际核构型，使用碳 2p$_z$ 子空间投影验证 ROHF 活性轨道。对于常规导数耦合（关闭 ETF），把 16 个 alpha、15 个 beta 电子的 Slater 行列式重叠直接与 CI 系数收缩，计算带符号的中心差分。它独立检验开壳层 NAC 中闭壳层核心的空间轨道导数处理；同一位移也用于检验两个势能面的能量梯度。
- 按两个原子组的对称 Löwdin 电荷构造局域初态。用具有相对复相位的 CI 叠加态独立重建全 AO 一阶密度，核对片段净电荷。轨迹检查全系统始终有 31 个电子，总净电荷始终为 +1。
- 保持起始核动量为零、初始空穴片段 A、种子和总时长相同，比较 6 × 2 与 12 × 1 au 的短轨迹。传播采用 ETF，报告每个时间点的相干电荷、活性面、能量漂移与态重叠。这是理想化孤立堆叠的短时演示，未采用周期边界、热采样或极化环境。

实测初始 A 净电荷为 +0.99990820 e。0.2903 fs 后，粗、细步长分别得到 +0.98664928 e、+0.98664930 e，两者相差 1.45 × 10⁻⁸ e。两条轨迹均未发生离散跃迁。最大总能量漂移分别为 7.05 × 10⁻⁹ Ha、1.70 × 10⁻⁹ Ha，电荷和误差在 9 × 10⁻¹⁴ e 内。

可选依赖另用系统 Python 3.14.7 验证。运行 `PYTHONPATH=src python3 -m unittest discover -s tests -p '*hopping*.py'` 时，安装了 NumPy、未安装 PySCF 的环境中通过 17 项，跳过 6 项量化测试。增加 `python3 -S` 关闭第三方包后，通过 2 项纯输入测试，跳过 21 项。项目的量化环境中，23 项均通过。
