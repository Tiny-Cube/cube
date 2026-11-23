# Wren：分子材料的电子结构与动力学

Wren 从分子动力学结构出发，组织电子结构计算、重组能、载流子 KMC 和模型表面跳跃。
软件原名 cube。CMake 统筹多语言编译、测试和安装，Python 包由 uv 构建。主命令为 `wren`。

本手册区分**已验证实现**与**设计中的能力**。不能把孤立单体轨道直方图称为完整材料 DOS，
也不能把 Tully 模型 FSSH 验证称为真实材料的量化驱动非绝热动力学。

```{toctree}
:maxdepth: 2
:caption: 使用手册

installation
workflow
gpu-kmc
dos-guide
reorganization
fragment-hopping
modules
validation
```

```{toctree}
:maxdepth: 2
:caption: 理论与引用

dos-theory
references
```

## 按问题选择路线

| 目的 | 起点 | 结果与边界 |
|---|---|---|
| MD 的分子轨道变化与构象无序 | 提取 manifest 或 XYZ | PySCF 单体轨道 DOS，含时间/空间能级描述，不含分子间混合 |
| 有序晶体或无序超胞的电子谱 | 带晶格及原位置的完整帧 | PySCF PBC k 点 DOS；周期 KS/HF 谱，不自动等于准粒子谱 |
| 研究分子间耦合如何改变 DOS | 明确定义的有符号有效 Hamiltonian | Bloch 求和与对角化；必须收敛 k 网格/超胞 |
| 零场独立载流子输运 | 位点能、耦合、重组能、周期网络 | Marcus KMC、扩散/迁移率张量与诊断 |
| 分子两态的形变能或重组自由能 | XYZ、已有四点能量、梯度/振动模或两态能隙轨迹 | 按态和方法输出弛豫或系综估计 |
| 模型非绝热跃迁 | Tully、多态 LVC 或距离相关片段耦合 | FSSH、动量调整、可选退相干；真实量化接口待扩展 |

## 可复现实例

```bash
uv sync --locked --extra dos --extra dynamics --extra transport
uv run wren examples/dos/input
uv run wren examples/charge-carrier/input
uv run wren examples/surface-hopping/input
```

理论内容从通用 DOS 定义扩展至有序/无序有机半导体与 MD 采样，具体实现范围以
[DOS 使用说明](dos-guide.md)与[验收记录](validation.md)为准。
