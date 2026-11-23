# 从分子间距到片段表面跳跃

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
