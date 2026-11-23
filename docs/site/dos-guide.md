# 从 MD 或 PySCF 计算 DOS

## 选择你要计算的 DOS

| `mode` | 输入 | 当前实现 |
|---|---|---|
| `monomer-orbitals` | 单体 XYZ、单帧 manifest 或 frames.json | 中性闭壳层单体的 HF/KS 轨道谱；CPU/GPU 分子 SCF |
| `vertical-dscf` | 单体 XYZ、单帧 manifest 或 frames.json | 固定输入几何，从三个电荷态的单点总能量分别构建第一电离、电子添加能级分布 |
| `periodic` | 带原始分子位置与晶格的完整 manifest | CPU PySCF PBC KRHF/KRKS、全 k 网格；不裁剪分子子集 |
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
当前周期 GPU SCF 不接受未验证请求，也不回退 CPU。

周期计算没有天然共同真空零点。不同帧/晶胞之间的谱需明确平均静电势、真空能级或合适
参考的对齐；未对齐原始谱不能解释为绝对输运位点无序，也不应逐分子将 HOMO 强行归零。
当前尚未接通跨帧势能零点对齐，因此明确拒绝直接平均多个周期帧；可先逐帧计算并分析所需的对齐方案。

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

## 展宽、单位与输出

- `sigma` 单位 eV，表示 Gaussian 标准差，不是 FWHM；也不是自动估计的物理无序宽度。
- `emin/emax` 可裁剪能窗；不指定时自动覆盖能级两侧 8 sigma。
- 不指定 `points` 时自动使用不大于 sigma/5 的间距；显式网格必须不大于 sigma/2。
- 不对截断能窗强制重新归一化。`finite_grid_integral` 和 `expected_integral` 供核验。
- 全部能量保留原始参考，不默认逐帧/逐分子平移。

轨道和紧束缚模式输出 `dos.csv`（eV 网格、平均 DOS、逐帧 DOS）、`spectra.json`（原始能级、占据、k 点或
分子 ID）、`result.json`（参数、来源/源码哈希、归一化、积分和采样诊断），以及 SCF 日志/检查点。
此处 DOS 是能级计数，不能把其积分当作直接计算了电子密度或光吸收强度。
