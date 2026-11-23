# 从 MD 或 PySCF 计算 DOS

## 选择你要计算的 DOS

| `mode` | 输入 | 当前实现 |
|---|---|---|
| `monomer-orbitals` | 单体 XYZ、单帧 manifest 或 frames.json | 中性闭壳层单体的 HF/KS 轨道谱；CPU/GPU 分子 SCF |
| `periodic` | 带原始分子位置与晶格的完整 manifest | CPU PySCF PBC KRHF/KRKS、全 k 网格；不裁剪分子子集 |
| `tight-binding` | 明确的位点能与有符号实耦合 | 正交局域基 Hamiltonian 的 Bloch 谱；CPU/CuPy 数值后端 |

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
每条无向周期通道只列一次，程序补 Hermitian 共轭。局域基必须正交。
`spin-degeneracy=1/2` 显式选择计数；k 点权重和为 1，输出每模型晶胞的态密度。
只有耦合绝对值的转移积分 CSV 不足以可靠重建带干涉的谱，不能直接当作此模型的输入。

## 展宽、单位与输出

- `sigma` 单位 eV，表示 Gaussian 标准差，不是 FWHM；也不是自动估计的物理无序宽度。
- `emin/emax` 可裁剪能窗；不指定时自动覆盖能级两侧 8 sigma。
- 不指定 `points` 时自动使用不大于 sigma/5 的间距；显式网格必须不大于 sigma/2。
- 不对截断能窗强制重新归一化。`finite_grid_integral` 和 `expected_integral` 供核验。
- 全部能量保留原始参考，不默认逐帧/逐分子平移。

输出：`dos.csv`（eV 网格、平均 DOS、逐帧 DOS）、`spectra.json`（原始能级、占据、k 点或
分子 ID）、`result.json`（参数、来源/源码哈希、归一化、积分和采样诊断），以及 SCF 日志/检查点。
此处 DOS 是能级计数，不能把其积分当作直接计算了电子密度或光吸收强度。
