# 两个电子态的重组能

`$reorganization-energy` 可以单独执行。先确定过程中的两个电子态 A、B 和希望测量的物理量，再选择方法。

| 已有数据或任务 | 控制字段 | 输出 |
|---|---|---|
| 从 XYZ 或 MD 单体计算 S₀/S₁ | `re-method=4-points`，`state-a=S0`，`state-b=S1` | 两个优化结构、四个态总能量、两项形变能、根间隔与态特征诊断 |
| 比较不同电荷或自旋的两个 SCF 态 | 显式 `state-a`、`state-b`，或载流子简写 `carriers=hole,electron` | 真空或冻结 MM 中的四点结果。载流子简写另外生成 hole.csv 和 electron.csv |
| 已有其他量化程序算出的四个总能量 | `re-method=4-points`，`energy-source=external` | 按态对输出带单位的四个能量、两项形变能及绝热、垂直能隙 |
| 有参考态的频率、正则振动模和另一态的梯度 | `re-method=vertical-gradient` | 目标势能面的谐近似弛豫能、各模贡献和 Huang–Rhys 因子 |
| 有分别在 A、B 态平衡采样的两条轨迹 | `re-method=linear-response` | 能隙均值和方差、Stokes 和涨落估计、可选分块标准误 |

## 从一个 XYZ 计算 S₀/S₁

把下列内容保存为 `input`，并将 `molecule.xyz` 放在同一目录。启动命令为 `uv run --extra workflow wren input`。已提取 MD 帧时可将 `xyz` 行换成 `input = mol/frame2000/manifest.json` 和 `molecules = 1`。

```text
$reorganization-energy
do = 1
re-method = 4-points
state-a = S0
state-b = S1
xyz = molecule.xyz
qc-engine = pyscf
method = b3lyp
basis-set = 6-31G(d,p)
auxbasis = none
environment = vacuum
output = results/singlet
$end
```

这两个简写表示中性闭壳层的 SCF 基态与同电荷的最低单重 TDA 激发态。四个点采用同一基组、泛函与能量零点。每个几何上的激发态总能量为基态总能量加该几何的垂直激发能。当前适配器使用 CPU PySCF 的常规积分和解析梯度。优化时检查相邻评估几何间的跃迁振幅重叠，以及与相邻根的能量间隔。默认条件是重叠至少 0.7、最近能隙大于 0.05 eV；出现竞争根或态特征丢失时停止。优化收敛表明到达驻点，确认极小点还需要频率分析。[PySCF 激发态优化](https://pyscf.org/user/geomopt.html#excited-states) · [Hirata 与 Head-Gordon 的 TDA 论文](https://doi.org/10.1016/S0009-2614(99)01149-5)

`results/singlet/molecule/reorganization.json` 保存 `analysis.energies`（Hartree）、`relaxation_ev`（eV）及优化和根诊断。以 $R_A,R_B$ 表示对应态的优化几何，

$$
\lambda_A=E_A(R_B)-E_A(R_A),\qquad
\lambda_B=E_B(R_A)-E_B(R_B).
$$

`on_a` 是在 A 势能面上将几何从 $R_A$ 变为 $R_B$ 的能量代价；`on_b` 是在 B 面上将几何从 $R_B$ 变为 $R_A$ 的代价。`sum` 是两项之和，也等于四点模型中吸收与发射的垂直能隙之差。等价分子之间的一次局域激子交换可以组合供体回到基态与受体进入激发态的两项形变代价。异种分子须分别取相应供体和受体的项。对完整反应的两条等曲率 Marcus 自由能面，则每条面上的形变能已经是该反应的 $\lambda$。[Marcus 的原始模型](https://doi.org/10.1063/1.1742723)

仓库里的苯算例给出 $\lambda_{S_0}=0.140952$ eV 和 $\lambda_{S_1}=0.141546$ eV，采用 B3LYP/6-31G 和真空环境。复现用的分子坐标、控制文件及四个总能量均随验证记录保存。[苯的控制文件](https://github.com/Tiny-Cube/cube/blob/feat/wren/runs/wren_state_pairs_2026-11-23/benzene.input)

## 显式态定义与冻结 MM

不同的 SCF 电荷或自旋扇区可以直接声明。例如，中性单重态和阳离子双重态使用以下两行。

```text
state-a = {name="neutral", charge=0, multiplicity=1, solver="scf", root=0}
state-b = {name="cation",  charge=1, multiplicity=2, solver="scf", root=0}
```

`charge` 指 QM 区的总电荷。冻结环境的总电荷允许非零，四个点使用相同的 MM 原子、坐标、电荷与 LJ 参数。计算带净电荷的有限簇时，应检查不同簇大小和静电边界下结果的收敛。两份 SCF 定义只有名字不同、但电荷和多重度完全相同时，需要用可区分两个局域电子态的求解器产生能量，然后走外部输入路线。[PySCF 静电嵌入](https://pyscf.org/user/qmmm.html)

用显式簇时，`coordinates = coordinates.extxyz` 定义全体原子坐标、QM/MM 区域及冻结标记，`topology = topology.json` 定义原子身份、类型及部分电荷。可选的 `forcefield = forcefield.json` 给出成对 LJ 参数。这三个文件与 `environment = qmmm` 配合。已选定的簇会整体读入。另一种选择是 `input = manifest.json`、`mm-source = manifest` 和 `mm-vdw = tpr`。单独引用提取清单时，Wren 默认采用单体簇记录的全原子截断；可用 `mm-cutoff` 覆盖。只有缺少记录的历史清单才回到 12 Å。

从 GROMACS 提取清单导出这一组文件，可运行：

```bash
uv run --extra workflow python scripts/export_qmmm_cluster.py \
  mol/frame2000/manifest.json --molecule 1 --cutoff 6 --lj --output cluster-1
```

例子以 6 Å 全原子最短距离选择完整邻居分子镜像，并从 TPR 保留实际 QM–MM 类型对的 C6/C12。去掉 `--lj` 将只导出坐标与点电荷拓扑。做好切簇后，四点计算使用同一份冻结环境。

`carriers=hole,electron` 适合直接把结果交给载流子 KMC。显式态对留给相应过程的调用方组合。外部电子态可分别定义两态的求解器和根；自动 SCF 支持不同电荷或多重度，自动 TDA 目前支持闭壳层真空中的最低单重激发态。

## 导入四个总能量

```text
$reorganization-energy
do = 1
re-method = 4-points
energy-source = external
input = external-points.json
output = results/imported
$end
```

示例数据列出 F₂ 的中性/阳离子、中性/阴离子两个态对。输入中的 `energy_kind` 为 `total`，单位统一为 `hartree` 或 `eV`。每对提供 `E_A_RA`、`E_B_RA`、`E_B_RB`、`E_A_RB`，顺序由键名指定。结果写入 `results/imported/reorganization.json`。可以记录四点各自的源文件。Wren 计算文件哈希以便复核，原程序的收敛、优化和态身份需要结合原始输出判断。[仓库示例](https://github.com/Tiny-Cube/cube/blob/feat/wren/examples/reorganization/external-points.json)

## 垂直梯度法

`re-method=vertical-gradient` 配合 `input=vg-modes.json`。在参考态的平衡几何上取目标态梯度。若用参考态模的曲率 $k_i$ 近似目标态的曲率，质量加权的正则坐标中有

$$
\Delta Q_i=-g_i/k_i,\qquad
\lambda_i=g_i^2/(2k_i),\qquad S_i=\lambda_i/(\hbar\omega_i).
$$

输入指定参考态、正频率（cm⁻¹）及采用的模子空间。可以提供已投影的梯度，也可以提供按原子 ID 排列的质量（u）、笛卡尔能量梯度（eV/Å）和质量加权正交模向量。输出给出每个模的贡献与所选模上的总和。它估计**一个目标势能面**的弛豫。两态的曲率差异、明显的非谐性和未纳入的模都会影响结果。固定环境中应显式选择允许移动的自由度。[Avila Ferrer 与 Santoro 的谐振模型比较](https://doi.org/10.1039/C2CP41169E) · [ORCA ESD 的 VG 模型](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/esd.html#the-simplest-model-vertical-gradient-vg)

## 用两态 MD 平衡系综估计重组自由能

`re-method=linear-response` 接受两个独立采样文件，`ensemble-a` 来自态 A 的平衡轨迹，`ensemble-b` 来自态 B。每一行都要在**同一帧几何**上计算两个态的总能量。CSV 表头为 `frame,E_A_eV,E_B_eV`，帧号按时间等间隔递增。控制文件还需两个态定义、`temperature`（K）和 `environment`。只有一条中性 MD 轨迹时，仍需取得另一态的平衡采样才能使用这条双系综路线。

定义 $X=E_B-E_A$。在经典等曲率线性响应模型下，

$$
\lambda_{\rm St}=\tfrac12(\langle X\rangle_A-\langle X\rangle_B),\quad
\Delta G_{\rm LR}=\tfrac12(\langle X\rangle_A+\langle X\rangle_B),\quad
\lambda_{\rm var,s}=\frac{\mathrm{Var}_s(X)}{2k_BT}.
$$

输出分别报告两个系综的样本方差估计，并保留与 Stokes 值的差异。对轨迹先去除未平衡段。设置 `block-size`（样本数）可以根据两个独立轨迹的连续等长分块均值估计 Stokes 标准误；每条轨迹至少需要两个完整块，末尾不足一块的样本仍用于中心估计。用多个块长检查误差是否达到平台，并检查取样、能隙分布及线性响应假设。随环境一起采样得到的重组自由能与冻结 MM 中单体的四点形变能对应不同的自由度。[有机分子材料的能隙和涨落估计](https://pmc.ncbi.nlm.nih.gov/articles/PMC12805572/)

## 运行仓库示例

```bash
uv sync --locked --extra workflow
uv run wren examples/reorganization/singlet.input
uv run wren examples/reorganization/scf-pair.input
uv run wren examples/reorganization/external.input
uv run wren examples/reorganization/vertical-gradient.input
uv run wren examples/reorganization/linear-response.input
uv run wren examples/reorganization/separated.input
```

外部四点和两系综能隙分析只依赖 Python 标准库。笛卡尔振动模投影使用 NumPy。激发态 QM/MM 优化、激发态 GPU 解析梯度、可极化或 GDMA 多极嵌入以及双态 Hessian 的 Duschinsky 模式混合还需要独立实现和数值验证。
