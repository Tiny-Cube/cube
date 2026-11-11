# 通用重组能核心与当前工作流

`$reorganization-energy` 是可独立运行、供各主模块复用的子模块。目标是指定任意两个
电子态及重组算法，不能把重组能定义限定成 hole/electron。

当前自动量化工作流仍支持中性/阳离子、阴离子的 vacuum 或冻结 MM 点电荷+LJ 四点法。
通用 A/B 四点分析核心已实现，其他方法和任意态的自动量化尚未接通。
具体路线见 [设计](plans/2026-11-11-general-reorganization-design.md) 和
[方法/软件调研](research/reorganization-methods.md)。

## 分析来自任意程序的四个总能量

这个 Python API 不执行量化、不需要 NumPy/PySCF/GPU；数据来源可为任意已验证的量化结果。
示例值是两条不同曲率的合成谐势能面，不是材料的 S0/S1 计算结果。

```python
from universal_kmc.reorganization import (
    ElectronicState, FourPointEnergies, analyze_four_point_pair,
)

a = ElectronicState('S0', charge=0, multiplicity=1, solver='scf')
b = ElectronicState('S1', charge=0, multiplicity=1, root=1, solver='tddft')
energies = FourPointEnergies(
    e_a_ra=0.0, e_b_ra=15.5, e_b_rb=2.0, e_a_rb=4.5, unit='eV',
)
result = analyze_four_point_pair(a, b, energies)
print(result['relaxation'])
# {'on_a': 4.5, 'on_b': 13.5, 'sum': 18.0}
```

四个参数依次是 `E_A(R_A)`、`E_B(R_A)`、`E_B(R_B)`、`E_A(R_B)`。
显式单位仅接受 `hartree` 或 `eV`，输出保留该单位，不隐式转换或混用能量单位。
输入必须是总能量，TDDFT 的垂直激发能不能直接填作 `E_B`。
`root` 依赖求解器/自旋扇区的定义；声明 S1 或某个 solver 名称不会启动求解或验证态身份。

结果保留 `on_a`、`on_b`、`sum` 和垂直/绝热能隙，不自动把 sum 或 sum/2 命名为任何
反应的 Marcus 重组能。状态为 `analyzed`；极小点、态身份、环境一致性不能由四个标量
证明，因此对应验证标记为 false，调用方需要检查原量化数据。

## 现有输出兼容性

载流子流程中 A=中性、B=带电，映射为：

| 通用核心 | 已有载流子输出 |
|---|---|
| on_a | neutral_relaxation |
| on_b | charged_relaxation |
| sum | total |

既有 JSON 与 hole.csv / electron.csv 字段不变。自动计算仍检查负弛豫能及优化/SCF 收敛。
本次数学核心拆分会改变实现来源指纹；旧目录如果与当前源码不一致，`--resume` 将按现有
规则拒绝复用，应选择新输出目录。仅对已有四点结果做上述数学分析无需重新执行量化。
