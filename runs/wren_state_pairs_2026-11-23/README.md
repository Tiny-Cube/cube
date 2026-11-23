# Wren 多语言构建与电子态验证

以下命令从项目根目录运行，示例的配置、单位和适用物理模型见 [Wren 手册](../../docs/site/index.md)。

| 验证 | 输入与结果 |
|---|---|
| H₂ 的四点法与垂直梯度 | `validate_h2.py`、`h2.xyz`、`h2-vg.json`、`h2-validation.json`。同一参考几何下独立核验解析伸缩曲率，分别记录四点和 VG 结果 |
| 苯的 S₀/S₁ 重组能 | `uv run --extra workflow wren runs/wren_state_pairs_2026-11-23/benzene.input`。真空 B3LYP/6-31G + TDA，两个形变能及根诊断见 `benzene-validation.json` |
| GPU 解析激发态梯度实验 | `validate_gpu_tda.py` 和 `public-pyscf-probe.json`。直接使用 PySCF 和 GPU4PySCF 公共接口计算 H₂，CPU/GPU 各自优化再比较。Wren 的自动 S₀/S₁ 适配器当前只开放 CPU |
| 实际 MD 的垂直 ΔSCF 和两位点 PDOS | `gpu-dos.json`，同一 86 原子单体的三个电荷态计算，以及 CPU/CuPy 的投影和 IPR。记录了成功的 B3LYP 与首次 HF 阴离子未收敛的情况 |
| GROMACS 冻结簇往返 | `cluster-export.json`，第 2000 帧分子 1、6 Å。核对原始和回读的坐标、电荷、原子顺序、TPR 类型对 LJ 能量及梯度 |
| 分离输入的完整四点计算 | `examples/reorganization/separated.input`，以水为 QM、冻结 Na⁺ 为 MM；`separated-validation.json` 检查相同环境指纹、四个态总能量和优化结果。参数用于输入和求解流程验证 |

安装验证覆盖 Linux GNU C++/Fortran、Cargo、Ninja 的 Release 和 Debug，外部 CMake 工程可链接两套原生库。Python-only 配置也经过安装前缀重定位。完整构建命令在[安装文档](../../docs/site/installation.md)。实际测试摘要见 `cmake-validation.txt`。
