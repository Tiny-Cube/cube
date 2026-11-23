# KMC、重组能与表面跳跃

## 通用 KMC 与载流子

通用连续时间 KMC 用事件率之和决定等待时间，再按相对速率选择事件。Wren 保留 Python、
C++、Rust、Fortran 四套实现与一致随机数验证；通用引擎不要求模型必须是电荷输运。

静态单载流子模型目前使用非绝热 Marcus 速率：

$$k_{ij}=\frac{2\pi}{\hbar}|J_{ij}|^2\frac{1}{\sqrt{4\pi\lambda_{ij}k_BT}}
\exp\left[-\frac{(\Delta G_{ij}+\lambda_{ij})^2}{4\lambda_{ij}k_BT}\right].$$

该模型需要弱耦合、适当的核自由度与重组能近似，并非任意有机晶体输运的统一极限。
零场平衡扩散可通过 Einstein 关系得到迁移率；独立轨迹的无包裹位移、协方差斜率和统计
不确定性要在扩散时间区间检查。原始定义见
[Marcus，1956](https://doi.org/10.1063/1.1742723)。

```bash
uv run --extra transport wren examples/charge-carrier/input
```

该示例使用合成周期链。在同一材料控制入口选择 `backend = gpu`，可以在 NVIDIA GPU 上并行推进独立单载流子轨迹。多载流子相互作用、有场漂移和动态 Hamiltonian 仍待扩展。[GPU KMC 用法](gpu-kmc.md)

## 通用两态重组能

对 A/B 态分别保留

$$\lambda_A=E_A(R_B)-E_A(R_A),\quad
\lambda_B=E_B(R_A)-E_B(R_B).$$

两项之和是否等于目标转移过程的 Marcus 参数取决于模型。完整反应的等曲率两势能面各自
的形变能通常均为 lambda，此时和为 2lambda；等价分子间跳跃常组合中性和带电分子的形变代价。
不能对任意态对都自动相加或除以二。

独立的重组能入口支持 S₀/S₁ 真空 TDA 四点优化、显式电荷/自旋扇区的 SCF 四点优化、外部总能量导入、垂直梯度和双系综线性响应。详见[重组能输入、理论与限制](reorganization.md)。

## Surface hopping

模型 FSSH 根据电子态振幅演化和活性势能面上的经典核轨迹进行随机跳跃；接受跳跃时沿
非绝热耦合方向调整动量，能量不足的跳跃受挫。方法源自
[Tully，1990](https://doi.org/10.1063/1.459170)。

Wren 已整合 Tully 1/2/3 和实固定正交基的多态 LVC provider，支持相位连续、电子子步、
能量守恒动量调整及可选 EDC。退相干近似参考
[Granucci、Persico，2007](https://doi.org/10.1063/1.2715585)。

```bash
uv run --extra dynamics wren examples/surface-hopping/input
uv run --extra dynamics wren-surface-hopping examples/surface-hopping/tully1.toml --check
```

对于 MD 抽取的有限分子簇，有两种描述可选。已有局域位点能和拟合耦合时，`fragment-sites` 推进片段中心和电子振幅。直接给分子几何时，`pyscf-sa-casscf` 在每个核步重新求相同自旋多重度下的态能、解析力和非绝热耦合，用相邻几何波函数重叠跟踪相位。单重态和双重态均有 CPU 算例。乙烯二聚体阳离子例子可指定空穴初始局域在分子 A，并用相干跃迁密度输出两个分子的净电荷。周期体相与长时间迁移率需独立研究。[控制文件、轨迹和量化验证](fragment-hopping.md)
