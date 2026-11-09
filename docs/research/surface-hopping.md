# Surface hopping：理论、开源实现与 cube 路线

检索日期：2026-11-09。目标是有机材料中的电荷／激子输运。下文的项目借鉴点是针对 cube 的工程判断，不代表已经安装、逐个复现或直接移植它们的源码。

## 推荐阅读顺序

| 文献 | 用途 | 对 cube 的具体作用 |
| --- | --- | --- |
| J. C. Tully, *Molecular dynamics with electronic transitions*, JCP **93**, 1061–1071 (1990), [DOI](https://doi.org/10.1063/1.459170)，[作者论文副本](https://xuv.scs.illinois.edu/chem540/GroupProjects/Benke_Hammer_SourcePaper1.pdf) | FSSH 原始算法；三个一维双态模型 | 首先核对概率符号、传播、动量调整和散射 |
| J. E. Subotnik, A. Jain, B. Landry, A. Petit, W. Ouyang, N. Bellonzi, *Understanding the Surface Hopping View of Electronic Transitions and Decoherence*, Annu. Rev. Phys. Chem. **67**, 387–417 (2016), [DOI](https://doi.org/10.1146/annurev-physchem-040215-112245) | 方法综述：波包分支、去相干、表象、详细平衡与局限 | 解释什么必须验证，避免把轨迹能运行当成物理可信 |
| G. Granucci, M. Persico, *Critical appraisal of the fewest switches algorithm for surface hopping*, JCP **126**, 134114 (2007), [DOI](https://doi.org/10.1063/1.2715585)，[作者全文](https://people.unipi.it/static/persico/art/surfhop.pdf) | 内部一致性与能量型去相干修正，Eq.17 | 可选 EDC；同时报告活跃态频率和电子振幅占据 |
| J. Spencer, F. Gajdos, J. Blumberger, *FOB-SH: Fragment orbital-based surface hopping for charge carrier transport in organic and biological molecules and materials*, JCP **145**, 064102 (2016), [DOI](https://doi.org/10.1063/1.4960144)，[UCL 原文记录](https://discovery.ucl.ac.uk/id/eprint/1516956/) | 面向凝聚相分子材料的片段轨道方法 | 电荷输运最贴近当前项目的材料路线 |
| A. V. Akimov, *Nonadiabatic Molecular Dynamics with Tight-Binding Fragment Molecular Orbitals*, JCTC **12**, 5719–5736 (2016), [DOI](https://doi.org/10.1021/acs.jctc.6b00955) | 基于片段轨道的电荷与能量转移 | 对照片段、正交化、去相干和界面输运的组织方法 |
| J. R. Mannouch, J. O. Richardson, *A mapping approach to surface hopping*, JCP **158**, 104111 (2023), [DOI](https://doi.org/10.1063/5.0139734)，[预印本](https://arxiv.org/abs/2212.11773) | MASH：以映射变量决定表面切换 | 后续独立算法对照；不是本次实现的 FSSH |

补充材料方向综述：[Charge Transport in Organic Semiconductors: The Perspective from Nonadiabatic Molecular Dynamics](https://doi.org/10.1021/acs.accounts.1c00675)。其讨论的瞬态离域提醒我们：局域速率 KMC 与带相干性的 surface hopping 有不同适用范围。激子需构造激子态的能量与耦合，不能把电子／空穴 HOMO/LUMO 转移积分直接当作一般激子耦合。

## 方法需要什么数据

KMC 给定状态及速率后抽取事件和等待时间。FSSH 则在核时间步上联合传播核和电子，所需信息包括多态能量、力、非绝热耦合或能够构造等价电子传播的波函数重叠，以及初始坐标、动量、电子态。一个静态 MD 帧与 `transfer_integrals.csv` 无法唯一确定这些时间演化数据。

在绝热表象、原子单位下：

\[
\dot c_i=-iE_i c_i-\sum_j(\mathbf v\cdot\mathbf d_{ij})c_j,\qquad
\mathbf d_{ij}=\langle\phi_i|\nabla_R\phi_j\rangle.
\]

核沿活跃态的力运动；跃迁概率由活跃态流向目标态的电子概率流确定。跃迁时需要能量守恒的动量调整。长时间电子相干可能带来偏差，因此去相干、frustrated hops、初始分布与时间步都属于模型定义，需要写入输出。

对于固定、正交、实 diabatic 基，当前接口由 `H(R)` 和 `dH/dR` 构造：

\[
\nabla E_i=U_i^T(\nabla H)U_i,\qquad
\mathbf d_{ij}=\frac{U_i^T(\nabla H)U_j}{E_j-E_i}\quad(i\ne j).
\]

该公式在近简并处需要更稳健的局域 diabatic 传播，且不能遗漏移动片段基自身的导数。当前实现对此明确限制，不对真实材料作未经验证的外推。

## 开源项目比较

| 项目与官方来源 | 已核实的内容 | 值得借鉴的部分 | 对 cube 的取舍 |
| --- | --- | --- | --- |
| [SHARC 4](https://github.com/sharc-md/sharc4)，[官方手册](https://sharc-md.org/?page_id=1454) | Fortran90 核心、Python 接口、wfoverlap、测试及量化实例；GPL-3.0-or-later | 量化程序与动力学分离，态追踪和多种耦合的统一接口 | 用作成熟参照；不把 GPL 源码复制到 cube 的 MIT 实现中 |
| [Libra](https://github.com/Quantum-Dynamics-Hub/libra-code)，[文档](https://libra-documentation.readthedocs.io/) | 模块化量子／经典动力学库；模型、态追踪、相位处理、NBRA 等；GPL-3.0 | 多哈密顿量、多传播算法、材料接口和基准组织 | 对通用 provider 和输运架构很有参考价值 |
| [Newton-X](https://newtonx.org/about/)，[下载及许可](https://newtonx.org/download/) | FSSH、局域 diabatic、Landau–Zener 等；NAC、overlap、曲率路线；GPL-3.0-or-later；当前官方下载需填表 | 初始采样→轨迹→统计分析；缺少 NAC 时的方法区分 | 不假定可匿名 git clone，也未填写下载表或发送邮件 |
| [PyRAI2MD](https://github.com/lopez-lab/PyRAI2MD) | Python/Cython FSSH，velocity Verlet，量化接口，机器学习模块；该仓库标为 MIT | 轨迹对象、propagator、量化接口分层 | 最适合参考 Python 原型组织；未移植其源代码 |
| [PySurf dynamics](https://github.com/DavidPicconi/pysurf_dynamics)，[对应方法论文](https://doi.org/10.1039/D5CP01194A) | 可扩展插件；FSSH、去相干、速度调整；论文扩展 mapping 方法；README 标为 Apache-2.0 | 模型 provider、算法插件、不同方法在同一模型上比较 | 适合后续 FSSH/MASH 对照，版本与文件许可仍应在实际复用时核对 |

这些项目的“开源”和 GPU 能力是不同问题。Python/Cython 或 GPU 神经网络并不自动意味着整套非绝热量化、NAC、QM/MM 以及轨迹传播都在 GPU 上运行。SHARC 的现有 PySCF 接口以 SA-CASSCF/PDFT 为主，官方手册列出了状态、NAC 和高级耦合支持限制；不能直接推断与我们现有 GPU4PySCF 基态 QM/MM 接口等价。

FOB-SH 论文很重要，但“论文使用 CP2K 内部实现”不能据此证明某段代码已经公开进入 CP2K 主线；本次没有把未核实的 FOB-SH 分支列为可直接安装的开源依赖。[相关材料论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11498209/)明确称其为 in-house 实现。

## cube 已实现与后续顺序

第一阶段已实现独立、原生的 Python/NumPy 实 FSSH 内核，Tully 1/2/3 与 JSON LVC 模型，支持任意态数／坐标数的 provider、可复现独立轨迹种子、幺正电子传播、相位符号连续、NAC 方向动量重标度，以及可选 EDC。完整运行方法与输出解释见 [使用说明](../surface-hopping.md)。这些能力是模型动力学，不代表已经完成材料输运。

下一阶段按以下顺序推进：

1. 片段电荷／激子哈密顿量接口，明确单位、基底、轨道相位、正交化和态追踪；电荷与激子的参数供应分开。
2. 对一小段真实 MD 序列检查位点能、带符号耦合、梯度、overlap 的连续性，分别定义含核反作用和 NBRA 模式。
3. GPU 量化检查每个所需电子态的能量、梯度、NAC／重叠支持，与 CPU 基准比较；不能以基态 SCF 成功替代该验证。
4. 轨迹集合的 GPU 批量传播，再做粒子／激子空间占据、IPR、MSD、扩散张量及迁移率。需要时间步、轨迹数、尺寸和时间窗的收敛报告。
5. 对模型与材料比较 FSSH 去相干方案及 MASH；多语言端移植在物理接口和基准稳定后进行。
