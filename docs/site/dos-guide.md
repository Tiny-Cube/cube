# 从 MD 或 PySCF 计算 DOS

## 选择你要计算的 DOS

| `mode` | 输入 | 当前实现 |
|---|---|---|
| `monomer-orbitals` | 单体 XYZ、单帧 manifest 或 frames.json | 中性闭壳层单体的 HF/KS 轨道谱；CPU/GPU 分子 SCF |
| `vertical-dscf` | 单体 XYZ、单帧 manifest 或 frames.json | 固定输入几何，从三个电荷态的单点总能量分别构建第一电离、电子添加能级分布 |
| `periodic` | 带原始分子位置与晶格的完整 manifest | PySCF 或 GPU4PySCF 周期 KRHF/KRKS，完整晶胞与等权 k 网格 |
| `tight-binding` | 位点能、带符号矩阵元，可提供片段重叠 | 正交或非正交模型的 Bloch 谱、位点投影和参与数 |

单体轨道分布用于观察构象、化学环境模型下的能级变化；当前单体模式使用真空独立分子，
不包含静电/极化环境或分子间轨道混合。完整材料 DOS 应根据问题采用周期或经过验证的有效 Hamiltonian。

## 独立分子例子

```bash
uv run --extra dos wren examples/dos/input
uv run --extra dos wren-dos examples/dos/monomer.toml
```

第二种使用独立 `[dos]` TOML 控制，与主输入的 `$dos` 字段相同。

```toml
[dos]
mode = "monomer-orbitals"
xyz = "molecule.xyz"
method = "b3lyp"
basis-set = "6-31G(d,p)"
qc-engine = "pyscf"
threads = 4
sigma = 0.05
output = "results/dos"
```

`input=manifest.json/frames.json` 可替代 `xyz`，二者互斥。`molecules="1,2"`
选择部分单体；多帧要求同一有序分子集合。当前为有限基组内的全轨道谱，不是完整连续谱。
闭壳层空间轨道按自旋简并权重 2 计数，单帧再按分子数平均；各帧等权平均。
占据数另存在原始谱中，DOS 不乘费米占据；占据电子谱与未占据态密度应分别解释。

## MD 的实际使用流程

1. 平衡轨迹、确认分子连接及周期盒，丢弃非平衡初期。
2. 用 `$get_dimers_monomers` 提取需要的帧。先选少量相隔足够远的帧做功能与成本检查。
3. 将 `$dos` 的 `input` 指向生成的 `frames.json`，先限制 `molecules` 验证量化设置。
4. 检查各帧 SCF 收敛、原始 HOMO/LUMO、DOS 积分和参考能量后，扩大采样。
5. 对轨迹长度、帧间隔、独立初始条件、基组、泛函、环境模型和展宽做收敛检验。

对同一条 MD 轨迹，不要把每帧当作独立样本。当前报告给出有限样本的静态/动态方差分解，
尚未自动估计相关时间和置信区间，`sampling_uncertainty_estimated=false`。

### 从 MD 计算垂直电离和电子添加能级

每帧取相同的一组分子，在每个分子原始坐标上执行三次真空单点 SCF。中性态取闭壳层单重态，
阳离子和阴离子取双重态。可先执行小分子示例，再将 `xyz` 换为轨迹提取目录的 `input`。

```bash
uv run --extra dos wren-dos examples/dos/vertical-dscf.toml
```

```toml
[dos]
mode = "vertical-dscf"
input = "mol/frames.json"
molecules = ["1", "2"]
method = "b3lyp"
basis-set = "aug-cc-pvdz"
threads = 4
sigma = 0.05
output = "results/vertical-dscf"
```

`vertical-energies.csv` 保留每帧每个分子的三份总能量（Hartree）、垂直 IP 和 EA（eV），
以及电子能量轴上的移除能级 $-IP$ 和添加能级 $-EA$。`removal.csv` 和 `addition.csv`
分别给出这两个通道的逐帧及平均分布，每通道全能窗积分各为 1。每个分子在每个通道贡献一个跃迁，
帧内按选取分子数平均，随后按帧等权平均。`spectra.json` 记录原始能量、电子态和各 SCF 状态文件；
`result.json` 逐通道给出归一化与能级波动。

步长、轨迹帧数决定构象采样，基组中的弥散函数尤其影响阴离子。可检查各带电态的
$\langle S^2\rangle$（理想双重态为 0.75），并通过增加弥散函数和比较电子密度分布检查束缚性。
`nonpositive_vertical_ea_count` 统计计算得到 EA ≤ 0 的样本数。这一模式给出指定自旋态和有限
分子模型下的第一加减电子能量分布。材料的环境极化、分子间混合和其他多体激发需要单独建模。

## 周期 PySCF

```bash
uv run --extra dos wren-dos examples/dos/periodic.toml
```

完整帧中每个单体的居中坐标加回 `transport_position_angstrom`，重建同一晶胞。
旧清单缺此位置时明确拒绝，不把所有分子叠在原点。
使用 `kmesh=[n1,n2,n3]` 的完整均匀网格、等权求和；`basis-set`、`pseudo`、
`auxbasis`、`precision` 等必须针对周期计算收敛。`auxbasis="auto"` 使用 PySCF 自动辅助基。
默认 `max-atoms=500` 是防止意外启动巨大超胞量化的显式限额，不是对物理体系的裁剪规则。
在 TOML 中设置 `qc-engine = "gpu-pyscf"` 启用 GPU 周期 SCF，设置 `dos-backend = "gpu"` 在 GPU 上计算元素投影与展宽。两个选项可独立设置。周期 DFT 两端均采用 Gaussian density fitting 和 `grid-level` 指定的 Becke 网格。输出目录的 `spectra.json` 记录实际 SCF 后端、显卡、电子数、每个 k 点的轨道占据与带边位置。

随附的实验萘晶胞含两个分子，来自 298 K、0.4 GPa 的单晶衍射精修。用同一结构逐步加密三个倒空间方向的采样，可以分别观察总能、带隙和 DOS 形状的变化。[晶体结构与精修条件](https://doi.org/10.1107/S0108768106026814)

```bash
uv run --extra dos wren-dos examples/dos/naphthalene-periodic.toml --kmesh-study 1x1x1 2x2x2 3x3x3
```

也可在统一的 `input` 中直接指定网格序列，输入路径相对于控制文件。以下示例在仓库根目录运行 `wren input`。独立 TOML 的 `[dos]` 表接受同一字段，例如 `kmesh-study = ["1x1x1", "2x2x2", "3x3x3"]`。

```text
$dos
do = 1
mode = periodic
input = examples/dos/naphthalene-cell.json
output = results/naphthalene
qc-engine = gpu-pyscf
dos-backend = gpu
method = pbe
basis-set = gth-szv
pseudo = gth-pbe
grid-level = 2
precision = 1e-7
conv-tol = 1e-8
sigma = 0.1
element-projections = true
kmesh-study = 1x1x1,2x1x1,2x2x2,3x3x3,4x4x4
$end
```

`<output>-kmesh/kmesh-study.json` 给出每个网格的采样间隙、最小直接间隙，以及相对上一个网格的变化量。各网格目录保留总 DOS、元素 PDOS 和原始能级。采样间隙允许价带顶和导带底来自不同 k 点，直接间隙在同一 k 点取差。`dos_total_variation` 使用同一晶胞的原始能量零点、相同 Gaussian 宽度和覆盖所有轨道的公共能量网格，计算两条谱之差的绝对值积分除以两倍每晶胞态数。该量接近 0 表示两次展宽谱接近。带边、基组、结构采样仍需分别检查。使用 `--resume` 可复用输入和文件校验均通过的已完成网格。

### 从多个 MD 帧计算周期 DOS

用提取模块生成的 `frames.json` 作为 `input`，各帧保持相同的分子编号、原子种类及顺序。Wren 对完整晶胞逐帧求解，再在共同能量网格上平均。三个方向的 k 网格可先用代表帧收敛。对于中性三维有机晶胞，下面的设置可以直接从坐标和 GTH 赝势计算内部参考。

```text
$dos
do = 1
mode = periodic
input = mol/frames.json
output = results/md-periodic
qc-engine = gpu-pyscf
dos-backend = gpu
method = pbe
basis-set = gth-szv
pseudo = gth-pbe
kmesh = 2,2,2
periodic-reference = cell-mean-local
sigma = 0.10
element-projections = true
block-size-frames = 2
$end
```

每条轨道的绘图能量为 $E_\mathrm{raw}-E_\mathrm{ref}$。总 DOS 与元素 PDOS 使用相同的平移。平均对帧等权，对每帧的 k 点按其积分权重求和，单位为每 eV、每晶胞自旋态数。在 `spectra.json` 中，`orbital_energies_ev` 与 `band_edges` 保存原始 SCF 参考下的值，`energy_alignment` 保存减去的参考值及变换后的带边。

| `periodic-reference` | 参考值的来源及用途 |
|---|---|
| `cell-mean-local` | 使用 PySCF 三维周期约定中 Hartree 势与离子局域 GTH 赝势之和的晶胞平均值。固定组成、赝势与体积时该值恒定，结构变化造成的带边移动会保留。变胞时按实际体积逐帧更新。适合比较同一组成的体相内部谱。 |
| `external` | 从 `reference-file` 读取每帧在原 SCF 能量零点下测得的参考（eV）。例如从同一中性薄膜超胞的平坦真空区取平均静电势。可用 `reference-label` 记录取样位置和方法。 |
| `frame-vbm` | 每帧减去自己的价带顶。可以比较展宽、带形及带隙。共同的价带顶涨落随之移除。要求每一帧均存在分开的占据带与空带。 |

外部表格按真实 MD 帧编号匹配。可增加 `manifest_sha256` 列，将每行绑定到提取的几何文件。下面两行仅展示表格语法，实际计算需填入该帧测得的电子势能参考。

```csv
frame,reference_ev
1980,1.832
2000,1.794
```

### 轨迹的分块误差

在单体或周期 DOS 计算中，`block-size-frames = B` 将连续等间隔帧分为每块 $B$ 帧。至少需要两个完整块。`dos.csv` 保留全部帧的谱与均值，`dos-blocks.csv` 给出每一能点上各完整块的均值、完整块总体均值和该均值的标准误。末尾不足一块的帧数记录在 `result.json` 中。该误差描述所选轨迹的时间采样。

已有 DOS 结果时，可直接比较多个块长。将独立 TOML 的 `[dos]` 表中 `output` 指向已完成的计算目录，运行以下命令。

```bash
wren-dos dos.toml --block-study 1 2 4
```

`<output>-blocks/block-study.json` 记录每个块长、有效帧数和丢弃尾帧数，`block-2.csv` 等文件列出标准误随能量的变化。过短的块会保留块间相关，过长的块会导致独立统计量太少。应结合更长轨迹、多个独立初态和不同超胞检查有意义能量范围内的稳定性。[Flyvbjerg 与 Petersen，1989](https://doi.org/10.1063/1.457480)

## 有效 Hamiltonian

输入 JSON 的字段为 `schema="wren-dos-tb-1"`、
`coupling_convention="signed_real_consistent_orbital_phases"`、`site_energies_ev`、`pairs`。
每个 pair 是 `[i,j,[image_a,image_b,image_c],signed_coupling_ev]`；编号从 0 开始。
每条无向周期通道只列一次，程序补 Hermitian 共轭。这个版本的输入使用正交局域基。
`spin-degeneracy=1/2` 显式选择计数；k 点权重和为 1，输出每模型晶胞的态密度。
只有耦合绝对值的转移积分 CSV 不足以可靠重建带干涉的谱，不能直接当作此模型的输入。

对于有重叠的归一化片段轨道，使用 `schema="wren-dos-tb-2"` 和
`matrix_element_convention="raw_signed_H_and_S_in_same_orbital_gauge"`。分别在
`hamiltonian_pairs_ev` 与 `overlap_pairs` 中列出原始带符号的 $H_{ij}$（eV）与 $S_{ij}$（无量纲）。
两表可有不同的非零项，使用相同的轨道编号、相位和周期像约定。省略的非对角元为零。
`site_energies_ev` 给出原始 $H$ 的对角元，程序补单位现场重叠。
已经通过二聚体正交化得到的有效转移积分与原始 $H_{ij}$ 有不同的定义。

```bash
uv run --extra dos wren-dos examples/dos/nonorthogonal-chain.toml
```

### 查看每个位点对能谱的贡献

在 `mode="tight-binding"` 的控制文件中设置 `site-projections=true`，同时导出位点投影 DOS
和参与数。示例由两个位点构成模型晶胞，晶胞内和跨晶胞的耦合均为 −0.1 eV。

```bash
uv run --extra dos wren-dos examples/dos/projected-chain.toml
```

| 输出 | 如何读取 |
|---|---|
| `pdos.csv` | 第一列为 eV，`site_0`、`site_1` 等列对应输入中从 0 编号的位点。在每个能量点把各列相加即可得到 `dos.csv` 的总 DOS。 |
| `participation.csv` | 每一行对应一个 k 点的一组能级。`k_index` 对应 `spectra.json` 的 k 点次序；`first_band`、`last_band` 从 0 计数，均含端点。`rank=1` 时给出 IPR 和参与数 `1/IPR`。 |
| `result.json` 中的 `site_projections` | 记录每个位点的预期积分、实际能窗积分及位点和相对于总 DOS 的最大绝对残差。裁剪能窗会减少积分。 |

此例使用 32 个 k 点和自旋权重 2，两个位点的全能窗积分各为 2，总积分为 4。
在分数波矢 $(1/2,0,0)$ 处，两个耦合恰好抵消，`rank=2`。此行的单态 IPR 和参与数留空，
`subspace_profile_ipr=0.5` 表示整个二维子空间在模型晶胞内的两个位点均匀分布。

`degeneracy-tol-ev` 默认为 $10^{-8}$ eV。每个 k 点独立分组，相邻能级差小于或等于阈值
就归入同一组。该阈值只影响参与数的分组，PDOS 保留每个本征值自己的能量。
投影需要存储各 k 点的位点权重；位点数为 $N$、k 点数为 $N_k$ 时，这部分内存随 $N_kN^2$ 增长。

### 从 MD 帧查看各元素的贡献

在 `monomer-orbitals` 或 `periodic` 模式加入 `element-projections = true`。
例如对提取帧中的 1 号分子使用以下控制段，运行 `wren input`。

```text
$dos
do = 1
mode = monomer-orbitals
input = mol/frames.json
molecules = 1
element-projections = true
method = hf
basis-set = sto-3g
auxbasis = weigend
sigma = 0.1
output = dos-mol1
$end
```

`pdos.csv` 的第一列为 `energy_ev`，其余列以元素符号命名。`dos.csv` 保存总 DOS
及逐帧曲线。两者使用相同能量网格、自旋权重和展宽；逐行求和可核对元素分解。
分子模式在每帧内对选中分子等权，再对帧等权；周期模式按每晶胞的 k 点权重求和。
`result.json` 的 `element_projections` 保存各元素的离散谱权重、能窗内积分和逐点求和残差。
上述最小基组便于功能核验，材料计算需另做基组、构型和环境收敛。

## 大型无序网络的稀疏 KPM

已有正交有效哈密顿量时，可在同一 `$dos` 段选择 KPM。输入的位点能、带符号耦合、
分子编号与周期像沿用前文紧束缚格式。对真实材料，轨道相位约定和位点能的共同参考需来自参数化过程。

```text
$dos
do = 1
mode = tight-binding
solver = kpm
input = network.json
kmesh = 1,1,1
spin-degeneracy = 1
kpm-moments = 128
kpm-trace = stochastic
kpm-vectors = 16
kpm-seed = 314159
points = 512
dos-backend = cpu
output = dos-kpm
$end
```

运行 `wren input`。仓库另附可直接运行的双位点周期示例
`uv run --extra dos wren-dos examples/dos/kpm-chain.toml`，使用 `kpm-trace = exact` 对每个基向量求迹。
稀疏 KPM 输出总 DOS。需要非正交 $S$、单态能量、IPR 或位点投影时，使用默认的 `diagonalization`。

| 参数或文件 | 含义 |
|---|---|
| `kpm-moments` | 保留 $T_0$ 至 $T_{M-1}$。默认 128。逐步加倍检查关心的带边或峰形。 |
| `kpm-vectors`、`kpm-seed` | 随机迹的独立向量数和种子。默认 16、0。提高向量数并用独立种子检查随机波动。 |
| `kpm-bounds-ev` | 可选的缩放能量范围，须覆盖程序从矩阵元得到的全部 Gershgorin 界。程序再扩大半宽 1%。可用于让几个模型具有相同的谱分辨率；全平谱需指定有限范围。 |
| `emin`、`emax`、`points` | 输出窗口及等宽能量箱数。默认窗口覆盖缩放区间，默认箱数为 $\max(512,4M)$。裁剪窗口会减少积分，不改变缩放。 |
| `dos.csv` | 每行给出箱中心、箱平均 DOS、该箱的随机迹标准误以及左右边界。总态数等于 `sum(dos * (right_ev - left_ev))`。 |
| `moments.csv` | 按阶数记录含 k 点与自旋权重的原始迹矩和标准误。精确迹的随机误差为零，多项式截断仍影响分辨率。 |
| `result.json` | 记录核、缩放、最大稀疏矩阵非零元数、预期态数、窗口积分与输入哈希。 |

单个 MD 快照的随机迹标准误衡量数值估计误差。讨论材料无序还需增加快照、分子和超胞，并检查时间相关性。
CPU 使用 SciPy CSR。选择 `dos-backend = gpu` 时，CSR 转入 CuPy，Chebyshev 递推在 GPU 执行，
能量箱积分于主机完成。当前本模块的实测结果见[验证记录](validation.md)。

## Gaussian 展宽、单位与输出

- `sigma` 单位 eV，表示 Gaussian 标准差，不是 FWHM；也不是自动估计的物理无序宽度。
- `emin/emax` 可裁剪能窗；不指定时自动覆盖能级两侧 8 sigma。
- 不指定 `points` 时自动使用不大于 sigma/5 的间距；显式网格必须不大于 sigma/2。
- 不对截断能窗强制重新归一化。`finite_grid_integral` 和 `expected_integral` 供核验。
- 全部能量保留原始参考，不默认逐帧/逐分子平移。

轨道和紧束缚模式输出 `dos.csv`（eV 网格、平均 DOS、逐帧 DOS）、`spectra.json`（原始能级、占据、k 点或
分子 ID）、`result.json`（参数、来源/源码哈希、归一化、积分和采样诊断），以及 SCF 日志/检查点。
此处 DOS 是能级计数，不能把其积分当作直接计算了电子密度或光吸收强度。
