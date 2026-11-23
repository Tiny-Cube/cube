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
format = xyz,gjf
$end
```

`frame=-1` 仅取最后一帧；否则从给定起始帧开始按 `frame_step` 提取到轨迹末尾。
`fragment-file` 用分子内部从 1 开始的编号定义 core。dimer 判据为两个 core 的最短原子距离严格小于阈值。
QM/MM 环境判据为完整中央分子与完整邻居分子的最短全原子距离，选中后保留整个邻居及周期镜像。

输出包括每帧的 `manifest.json`、单体、簇及 dimer；多帧还生成 `frames.json` 索引。
当前 GROMACS 提取要求同种中性完整分子、原子顺序一致；LAMMPS 自动导入尚未接通。

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
input = mol/frames.json
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
