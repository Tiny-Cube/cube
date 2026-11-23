# 片段模型与量化驱动的表面跳跃

`fragment-sites` 把一个有限分子簇表示为若干可移动的分子中心，每个中心带一个局域电子态。给出耦合随间距变化的拟合关系，即可跟踪初始局域的电子布居如何传播。有机晶体中，热运动确实会调制近邻分子间的转移积分，实际幅度取决于堆积与振动模式。[Troisi 与 Orlandi，2006](https://doi.org/10.1021/jp055432g)

## 运行三片段算例

```bash
uv run --extra dynamics wren-surface-hopping examples/surface-hopping/fragment-sites.toml --check
uv run --extra dynamics wren-surface-hopping examples/surface-hopping/fragment-sites.toml
```

控制文件选择 `[model] type = "fragment-sites"`，`file = "fragment-sites.json"` 给出中心与接触参数。也可以在 Wren 的 `$surface-hopping-dynamics` 段中把 `input` 指向这个 TOML。示例轨迹写在 `build/surface-hopping/fragment-sites/`，使用同一组初始核位置和动量，指定在 A 处制备局域电子态。电子相干传播及随机势能面跳跃由 FSSH 引擎完成。[Tully，1990](https://doi.org/10.1063/1.459170)

## 从接触参数构造 Hamiltonian

记 $r_{ij}=|\mathbf R_i-\mathbf R_j|$，$r^0_{ij}$ 由模型文件中的参考中心计算，$\Delta r_{ij}=r_{ij}-r^0_{ij}$。每条无向接触只填写一次。位点能和转移积分服从下式。

$$
V_\mathrm{nuc}=\frac12\sum_{i<j}k_{ij}\Delta r_{ij}^2,\qquad
H_{ii}=\epsilon_i+V_\mathrm{nuc}+\sum_{j\sim i}g_{i,ij}\Delta r_{ij},
$$

$$
H_{ij}=H_{ji}=J^0_{ij}\exp(-\beta_{ij}\Delta r_{ij}) .
$$

未列出的接触具有零电子耦合与零弹簧常数。公共弹簧势进入每一条绝热势能面一次。电子耦合的径向导数是 $-\beta_{ij} H_{ij}$。位点能斜率和弹簧力也纳入解析梯度。每一对中心受到等大反向的径向贡献，整体平移与旋转保持 Hamiltonian 不变。引擎在固定正交基中对角化 $H$，从 $\langle m|\nabla H|n\rangle/(E_n-E_m)$ 得到非简并态之间的导数耦合。

| JSON 字段 | 物理量及单位 |
| --- | --- |
| `sites[].id`, `center_bohr` | 唯一片段名、参考中心的三维坐标（Bohr） |
| `mass_me`, `energy_hartree` | 片段总质量（电子质量）、参考位点能（Hartree） |
| `pairs[].sites` | 两个片段 ID，按此顺序对应两项位点能斜率 |
| `coupling_hartree`, `decay_per_bohr` | 参考有符号耦合 $J^0_{ij}$、衰减系数 $\beta_{ij}$（Bohr$^{-1}$） |
| `spring_hartree_per_bohr2` | 接触弹簧常数（Hartree/Bohr²） |
| `site_slopes_hartree_per_bohr` | 两个端点各自的 $g_{i,ij}$（Hartree/Bohr），省略时均为零 |

初始动量 `[initial].momenta` 按文件中片段顺序排列，每个中心占 $x,y,z$ 三项，单位为原子单位动量。可通过 `coordinates` 显式指定偏离参考构型的初始中心，弹簧平衡距离仍由 `center_bohr` 决定。选 `fragment-site = "A"` 时，初态振幅是局域态在绝热基中的投影，初始活性面按投影概率独立抽样。若要从能量排序的某条绝热面出发，则填写从 0 开始的 `active-state`。

## 如何使用自己的 MD 和量化结果

1. 从已重构完整分子且已展开周期边界的 MD 帧取得质量中心、总质量和中心速度，用 $p_i=M_iv_i$ 换算初始动量。截取有限簇后显式列出想保留的分子对。
2. 为每种局域态定义位点能和相同的能量零点。从邻近参考构型的多个有符号 dimer 量化结果拟合 $J^0_{ij}$ 与 $\beta_{ij}$。用逐帧轨道重叠等方法统一轨道相位。只有耦合绝对值的表格不能确定环路上的相对符号。若平移滑移、转动或节点使耦合过零，单一径向指数无法表达这些变化，应检查拟合残差和留出的构型。
3. 从力场势能或力在参考构型附近的响应拟合接触弹簧，按需要拟合片段位点能斜率。这里的参考几何是公共弹簧势的平衡几何。某条电子态上的最低能构型可以不同。
4. 对多个独立 MD 初态运行，分别收敛核步长、电子子步、簇大小和轨迹数。单一起点上的随机轨迹重复抽样仅覆盖跳跃随机性。示例中的参数用于验证程序，尚未对特定材料拟合。

`result.json` 的 `site_ids` 对应 `trajectories.jsonl` 中 `diabatic_populations` 的顺序。可随时间检查电荷在 A/B/C 上的布居以及活性面的变化。历史字段 `diabatic_ipr` 存的是 $1/\sum_i p_i^2$，即有效参与的位点数。波函数完全局域在一个中心时为 1。有限簇内的布居移动本身不足以提取体相迁移率。

该模型采用固定实正交电子基、标量径向接触、有限边界及 CPU 动力学，适用范围取决于拟合的结构区间。移动片段轨道的显式基导数、分子取向、分子内振动与周期边界需由更详细的模型描述。基于原子分辨率片段轨道的 FOB-SH 可参见其方法论文与晶体应用。[Spencer、Gajdos 与 Blumberger，2016](https://doi.org/10.1063/1.4960144)、[Giannini 等，2019](https://doi.org/10.1038/s41467-019-11775-9)

## 用 PySCF 在每个核构型上重新求电子态

`pyscf-sa-casscf` 接受一个孤立分子或有限分子簇的 XYZ。在每个核步，CPU PySCF 使用 RHF（单重态）或 ROHF（双重态）轨道启动等权态平均 CASSCF 计算，给出相同总自旋各态的总能量、各态的解析核梯度和两两导数耦合。Wren 用当前活性态的力移动原子，传播各绝热态的电子振幅，并在接受跳跃时沿导数耦合方向调整核动量。适用体系要求所选活性空间能描述关注的激发或电荷重排。[Tully 的 FSSH](https://doi.org/10.1063/1.459170)、[PySCF 的态平均 CASSCF](https://pyscf.org/user/mcscf.html)、[SA-CASSCF 解析态间耦合](https://doi.org/10.1021/acs.jctc.6b00384)

从仓库根目录运行下面的乙烯示例。第一条命令使用常规 `$surface-hopping` 主入口，第二条只检查 LiH 的独立动力学输入。

```bash
uv run --extra dynamics --extra quantum wren examples/surface-hopping/ethene.input
uv run --extra dynamics --extra quantum wren-surface-hopping examples/surface-hopping/lih-casscf.toml --check
```

`$surface-hopping-dynamics` 的 `input` 指向动力学 TOML。其中 `[model]` 设为 `type = "pyscf-sa-casscf"`，`file = "ethene-casscf.json"`。JSON 中的 `xyz` 相对于该 JSON，元素顺序确定核坐标和各数组排列。乙烯的输入如下。

```json
{
  "schema": "wren-sa-casscf-1",
  "xyz": "ethene.xyz",
  "masses_u": [12.0, 12.0, 1.00782503, 1.00782503, 1.00782503, 1.00782503],
  "basis": "sto-3g", "charge": 0,
  "active_electrons": 2, "active_orbitals": 2, "states": 2,
  "initial_active_mo_1based": [8, 9],
  "orbital_gradient_tolerance": 0.000001,
  "threads": 1, "use_etfs": true
}
```

- XYZ 坐标以 Å 为单位，`masses_u` 给出每个活动原子的质量。TOML 中的 `momenta` 为按原子排列的三维动量，以原子单位计。已有 MD 速度时，先换算 $p_i=M_i v_i$。`active-state` 从 0 编号。指定活性态时，波函数从该态纯态出发；也可按下文的片段电荷制备相干初态。轨迹之间各自抽取跳跃随机数。MD 的多个初始构型或速度应分别运行。
- `initial_active_mo_1based` 从 **1** 编号，对应初始 RHF 或 ROHF 轨道。省略时，按参考轨道顺序选取紧接双占据核心的连续活性轨道。随后优化活性轨道，并在下一核构型投影前一步的轨道作为初猜。不同分子或基组需要重新检查轨道性质和活性空间。可选 multiplicity=1 或 2；活性空间分别需要偶数或奇数电子。活性空间以外的占据轨道均为双占据，体系使用实波函数。
- 默认 `min_state_overlap` 为 0.7，`min_gap_hartree` 为 $10^{-6}$ Hartree，`max_cycles` 为 100，`max_memory_mb` 为 2000。输出检查各态的 $\langle S^2\rangle$、活性自然占据数和导数方程收敛情况。遇到小能隙、按能量排序的根发生特征交换或求解失败时，轨迹停止并写入错误原因。

### 运动基组中的相位与导数耦合

PySCF 在每个构型使用随原子移动的 AO 基组。Wren 用两套 AO 间的交叉重叠积分构造前后多电子态的重叠矩阵，包括 CI 组态及双占据核心。每个能量序号对应的对角重叠须超过阈值，并优于行列中与其他已跟踪态的重叠。随后将 CI 的符号对齐再计算解析耦合。未纳入活性空间或根集合的电子态仍会限制动力学的可解释范围。[PySCF 多参考态接口](https://pyscf.org/user/mcscf.html)

引擎使用 $d_{ij}=\langle\Psi_i|\nabla_R\Psi_j\rangle$。适配器请求能量加权的解析态间导数，再显式除以 $E_j-E_i$，保证电子传播和动量调整使用一致的下标和符号。`use_etfs = true` 为默认值，应用电子平移因子的修正以消除耦合中的整体平移分量。`false` 给出包含移动基组贡献的通常导数耦合，可与相邻构型的全波函数重叠差分直接比较。[PySCF NAC 接口](https://pyscf.org/pyscf_api_docs/pyscf.nac.html)、[Fatehi 与 Subotnik 的电子平移因子](https://doi.org/10.1021/jz3006173)

### 阅读输出和验证适用范围

`trajectories.jsonl` 的每个采样点给出按能量排序的 `energies_hartree`、`gradients_hartree_per_bohr`、`derivative_couplings_per_bohr`、绝热态布居和当前活性态。梯度取负号后即为力。每条轨迹的相位重叠和初始活性轨道诊断保存在 `quantum`；发生跳跃尝试时会额外保留该步。可在 `result.json` 检查整个轨迹的最小能隙、最小相邻态重叠、最大能量漂移和波函数归一化误差。此接口未定义固定局域片段基，因此输出的绝热态布居不能直接解释为某个分子的电荷或激子布居。

验证脚本 `runs/wren_casscf_fssh_2026-11-23/validate.py` 可复算 LiH 两态的解析梯度和**带符号**的导数耦合，用包括 Li 核心电子在内的独立行列式重叠差分作对照。还可复算轻微扭曲的乙烯、检查初始两个 RHF 轨道对碳原子局部 $p_\pi$ 子空间的投影，并在相同演化时间比较两种核步长。两个算例都使用小活性空间和 STO-3G，用于检验计算链。材料的非绝热输运还需要适当的态集合、环境、初态采样、局域布居或位移观测量及长时间统计。此处的分子接口在真空 CPU 上运行，单重态之间的轨迹不包含自旋轨道耦合、周期边界或 QM/MM。

| 验证项 | 实际结果 |
| --- | --- |
| LiH 键长 3.5 Bohr 的常规 $d_{01,H_z}$ | 解析 $-0.122485$ Bohr$^{-1}$，跨几何全波函数的中心差分 $-0.122466$ Bohr$^{-1}$。Li 和 H 上两个态的梯度与差分之差均小于 $3.3\times10^{-6}$ Hartree/Bohr。 |
| 乙烯起始轨道 | 两个 CH₂ 局部平面法向相差 24.88°。初始 RHF #8 和 #9 投向两个碳的局部 $2p_\pi$ 子空间，权重分别为 0.934 与 0.961；两个碳上的轨道系数分别同相和反相。 |
| 乙烯的轨道优化 | 优化活性空间对该局部 $\pi$ 子空间的两个重叠奇异值平方为 0.539 和 0.957，表明其中一个方向含显著其他轨道成分。两个态的活性自然占据分别为 (1.959, 0.041) 和 (1.003, 0.997)。 |
| 乙烯核步长减半 | 从激发态传播相同的 2 个原子时间单位。0.5 和 0.25 两种步长下，终点核坐标最大差 $9.63\times10^{-8}$ Bohr，绝热布居最大差 $3.61\times10^{-12}$；总能量最大漂移分别为 $1.91\times10^{-6}$、$1.98\times10^{-6}$ Hartree。 |

两种步长的这些短轨迹均未抽中随机跳跃。另一个数值测试把均匀随机数固定为 0，在真实 LiH 量化帧上选中具有非零概率的向下跳跃，检验调整核动量后电子与核的能量交换。材料迁移率、激发态寿命和跳跃频率需要单独的体系级验证。

## 从有机二聚体上的局域空穴出发

```bash
uv run --extra dynamics --extra quantum wren examples/surface-hopping/ethene-hole.input
```

该算例将两个乙烯分子相隔 3.6 Å 平行堆叠，并在平面内错开 0.4 Å、0.2 Å。体系总电荷为 $+1$，总电子数为 31。以 ROHF 轨道 15 和 16 启动三电子两轨道的 SA-CASSCF，对同一个自旋双重态空间内最低两个绝热根等权优化。计算检查两个根的 $\langle S^2\rangle=0.75$。每个构型均重算能量、梯度、导数耦合及电子密度。

量化 JSON 新增如下设置。用户可将 `xyz` 换成抽取的其他二聚体，同时更新电荷、质量、活性空间和原子分组。每个原子恰好属于一个片段，原子编号从 1 开始，与 XYZ 行序对应。

```json
{
  "charge": 1, "multiplicity": 2,
  "active_electrons": 3, "active_orbitals": 2,
  "initial_active_mo_1based": [15, 16],
  "fragments": {"A": [1, 2, 3, 4, 5, 6], "B": [7, 8, 9, 10, 11, 12]}
}
```

在动力学 TOML 的 `[initial]` 中设置 `hole-fragment = "A"`。Wren 在所跟踪绝热态张成的空间内寻找使 A 的净电荷最大的归一化振幅，再按 $|c_I|^2$ 抽取每条轨迹的初始活性面。阴离子的 `electron-fragment` 选择净电荷最小的态。若给定空间无法唯一选出极值，程序报错。体系电荷、投影定义和所选态共同决定实际能达到的局域程度。

### 用跃迁密度计算片段电荷

记 $S$ 为当前几何的 AO 重叠矩阵，$P_F$ 为属于片段 $F$ 全部原子 AO 的对角选择矩阵。在 $S^{1/2}C$ 表示中构造对称 Löwdin 正交基上的分片电子数算符。双占据核心贡献只进入态空间矩阵的对角元，活性部分用一阶态密度和态间跃迁密度计算，保持任意电子振幅的相干项。[Löwdin，1950](https://doi.org/10.1063/1.1747632)、[Q-Chem 的布居分析说明](https://manual.q-chem.com/6.0/subsec_pop_anal.html)

$$
Q^F_{IJ} = Z_F\delta_{IJ} - \langle I|\hat N_F|J\rangle,
\qquad Q_F(t)=\sum_{IJ} c_I^*(t)Q^F_{IJ}(R(t))c_J(t).
$$

完整分组满足 $\sum_F Q^F=Q_{\mathrm{mol}}I$。程序在每个构型检验此矩阵恒等式。每个采样行记录 `fragment_net_charges_e` 和 `fragment_electron_populations`，以及两态振幅与片段电荷矩阵，供独立重算。一次随机轨迹的活性势能面与相干波函数的片段电荷各有含义；只用活性面电荷或按 $|c_I|^2$ 混合各态电荷都会丢失干涉项。电子部分使用密度矩阵，片段电荷依赖所选布居划分及基组。

开壳层轨道来自 ROHF，完整 CASSCF 计算与响应保留指定的双重态 CI。在 PySCF 2.15 的态间导数实现中，双占据 inactive core 的扣除项使用闭壳层的空间密度收缩；Wren 为该次 NAC 求值提供相应的空间轨道梯度驱动。LiH⁺ 与乙烯阳离子二聚体均可用两相邻几何的多电子重叠独立检验这一处理。电子平移因子默认用于传播。符号差分检验使用未经 ETF 修正的导数耦合。[PySCF SA-CASSCF NAC 实现](https://pyscf.org/_modules/pyscf/nac/sacasscf.html)

在 PySCF 2.15.0 的实算中，片段 A、B 起始的净电荷分别为 +0.99990820 e 和 +0.00009180 e。传播 12 au（0.2903 fs）后的结果如下。

| 核步长（au） | A 净电荷（e） | B 净电荷（e） | 最大能量漂移（Ha） | 接受的跳跃次数 |
|---:|---:|---:|---:|---:|
| 2 | +0.98664928 | +0.01335072 | 7.05 × 10⁻⁹ | 0 |
| 1 | +0.98664930 | +0.01335070 | 1.70 × 10⁻⁹ | 0 |

两次计算的终点 A 电荷相差 1.45 × 10⁻⁸ e。两个根全程保持 $\langle S^2\rangle=0.75$，各时间点的总净电荷与 +1 e 的偏差小于 9 × 10⁻¹⁴ e。两条轨迹全程均在初始活性势能面上，片段之间的净电荷随电子振幅相干演化发生变化。完整数值与复算脚本位于仓库的 `runs/wren_doublet_fssh_2026-11-23/`。

该算例考察真空有限二聚体上的短时相干电荷重分布。单条轨迹、单一起始构型以及 Löwdin 净电荷随时间变化，均不足以给出体相迁移率、扩散系数或不可逆转移速率。晶体及无序薄膜还需要考虑环境、更多分子和构型采样，并为目标材料验证活性空间、态数、基组和动力学时间尺度。
