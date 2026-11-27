# 控制文件与 MD 流程

`wren input` 解读 `$section ... $end` 控制文件。所有相对路径以控制文件所在目录为基准。
`--plan` 只检查配置、显示任务选择，不读取轨迹、执行 SCF 或探测 GPU。

```text
$get_dimers_monomers
do = 1
tpr = ../anneal_a.tpr
gro = ../anneal_a_fix.gro
xtc = ../anneal_a.xtc
frame = -1
frame_step = 20
dimers_cutoff = 6
fragment-file = ../core.txt
monomers-qmmm-cutoff = 6
output = mol
format = hdf5
$end
```

`frame=-1` 仅取最后一帧；否则从给定起始帧开始按 `frame_step` 提取到轨迹末尾。
`fragment-file` 用分子内部从 1 开始的编号定义 core。dimer 判据为两个 core 的最短原子距离严格小于阈值。
QM/MM 环境判据为完整中央分子与完整邻居分子的最短全原子距离，选中后保留整个邻居及周期镜像。

默认输出为 `mol/geometry.wren.h5`，含每个选中帧的单体坐标、dimer 和各自的 MM 索引。用 `wren-xyz frame2000.mol8 --root mol -o mol8.xyz` 导出指定单体，用 `xyz2gjf mol8.xyz mol8.gjf` 转为 Gaussian 输入。兼容输出可显式指定 `format = xyz`。
提取要求同种中性完整分子和一致的分子内原子顺序。

## 从 LAMMPS data 和 dump 提取

```text
$get_dimers_monomers
do = 1
lammps-data = hf.data
lammps-dump = hf.dump
lammps-units = real
element-map = elements.txt
fragment-file = core.txt
frame = 0
frame_step = 1
dimers_cutoff = 2
monomers-qmmm-cutoff = 2
format = hdf5
output = ../../build/lammps/mol
$end
```

`lammps-data` 提供 `atom_style full` 的 Atoms（atom-ID、molecule-ID、type、q、x、y、z）、Masses 与 Bonds。每个分子需要一个正的 molecule-ID，多原子分子需连成一个共价片段。使用固定拓扑及固定电荷。每行元素映射给出数字类型和化学元素，例如 `1 H`、`2 F`。同一元素可有多个力场类型。分子内以 atom-ID 递增排序，`core.txt` 写这个顺序下从 1 开始的局部序号。HDF5 的 `topology/molecule_id` 和 `atom_id` 保留原始编号，`global_atom_index` 是将全体系按 atom-ID 排序后的索引。[LAMMPS data 格式](https://docs.lammps.org/read_data.html)

输入轨迹为单个未压缩原生文本 `dump custom`，每帧含全体系的 `id type`，以及 `x y z`、`xs ys zs`、`xu yu zu` 或 `xsu ysu zsu` 中的一组坐标。如提供多组，按上列顺序取第一组。帧内原子行可乱序。支持三个方向均为周期的正交盒与 restricted triclinic 斜盒。解析器从斜盒的 bounding bounds 和 xy/xz/yz 还原晶格，按键连接恢复跨边界分子。每个分子以最小 atom-ID 的原子定位在参考晶胞，镜像平移沿该帧的晶格矢量。要求每条键在最短周期镜像上可唯一确定。`frames/origin` 保存盒原点，单体中心相对于该原点。[LAMMPS dump 坐标](https://docs.lammps.org/dump.html) · [三斜盒约定](https://docs.lammps.org/Howto_triclinic.html)

`lammps-units` 显式指定 `real` 或 `metal`。两者距离均为 Å，时间分别为 fs、ps。使用可变时间步时，可用以下 LAMMPS 命令保存实际经过时间。

```text
dump traj all custom 200 hf.dump id mol type xs ys zs
dump_modify traj time yes units yes
```

已有轨迹省略 `ITEM: TIME` 时须在控制文件提供 `lammps-timestep = 0.5` 等实际固定步长。此时采用“步号 × 步长”，适用于同一时间原点与恒定步长。输出分别保存零起始帧索引 `frame`、原始步号 `lammps_timestep` 和 `time_ps`。[时间和单位字段](https://docs.lammps.org/dump_modify.html) · [LAMMPS 单位](https://docs.lammps.org/units.html)

在仓库根目录执行以下命令即可从随附的小盒得到两帧提取文件和真空单体轨道 DOS。

```bash
uv run --locked --extra workflow wren examples/lammps/input
uv run --locked --extra workflow wren-dos examples/lammps/dos.toml
```

首帧的 dimer `7_19_img_0_-1_0` 由分子 7 和沿 −b 平移的分子 19 组成，两个 core 原子的距离为 1.382932 Å。两帧的物理时间为 0 和 0.1 ps。该示例的 HF/STO-3G 轨道 DOS 积分为每个 HF 分子 12 个自旋态。对首帧的独立 GROMACS TPR/GRO/XTC 输入，测试逐项比较了镜像、二聚体几何和冻结邻居电荷位置，坐标容差为 0.00002 Å。邻居编号写在 `clusters.csv`，原子电荷来自原始拓扑。QM/MM 的 `mm-source=manifest, mm-vdw=none` 可按同一原始文件哈希恢复 LAMMPS 静电环境。LAMMPS pair style、混合规则和 Pair/PairIJ Coeffs 目前没有转换，`mm-vdw=tpr` 需要真实 GROMACS TPR。带色散排斥的 LAMMPS 模型可用显式分离的坐标、拓扑和完整 LJ 参数输入。[冻结 QM/MM 输入](reorganization.md)

## 模块独立与主流程

主模块 `$charge-carrier`、`$surface-hopping`、`$density-of-states` 组织各自子模块。
子模块 `do=1` 独立请求运行；`do=0` 禁用；省略时由主流程组织，多个请求只运行一次。
主模块显式要求某子模块、而子模块显式关闭时会在科学计算前报冲突。

```text
$density-of-states
do = 1
$end
$dos
mode = monomer-orbitals
input = mol/geometry.wren.h5
qc-engine = pyscf
method = b3lyp
basis-set = 6-31G(d,p)
sigma = 0.05
output = results/dos
$end
```

DOS 当前要求先存在提取清单，可独立调用；不隐式替换已有材料的量化设置。
`--resume` 只复用输入、实现指纹及产物哈希一致的成功结果；失败任务请用新目录，
不能把存在旧输出当作已完成。

## 显式冻结簇输入

重组能模块可分别引用带区域、冻结标记和稳定原子编号的 Extended XYZ，含类型和电荷的拓扑 JSON，以及可选的 LJ 力场 JSON。目录示例和计算入口见[重组能](reorganization.md)。GDMA 多极环境仍待接入。
