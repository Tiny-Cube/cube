# 分离坐标、拓扑、力场与静电数据的 QM/MM 输入设计

## 用户已确定的方向

用户要求坐标、拓扑和力场拆开放，后续可能引入 GDMA 生成的数据。
用户进一步确定：坐标使用 Extended XYZ，同时在逐原子列中携带 QM/MM 区域和冻结状态。
此前提出的自包含 .qmmm.json 或 .qmmm 单文件不再作为唯一推荐输入。
本文件记录设计方向；下列新输入字段及 GDMA 后端尚未实现，不应作为现有命令运行示例。

## 文件职责

- 坐标：全部 QM/MM 原子的元素和坐标，使用同一坐标系；采用 Extended XYZ。
  必须同时记录逐原子的稳定 atom_id、region（QM/MM）和 frozen（T/F）。
- 拓扑：稳定的原子/分子 ID、原子顺序、力场类型、参考部分电荷、连接关系，
  以及回溯原始体系的原子编号和周期镜像。拓扑中的连接关系不等于已实现共价 QM/MM 边界。
- 力场：相互作用模型、参数、单位和参数来源。现有 QM/MM 模型仍使用从 TPR 解析出的
  成对 C6/C12，必须保留覆盖项，不通过新的混合规则重建原参数。
- 控制文件：引用坐标、拓扑、力场及可选多极数据，设置 QM 总电荷/自旋、量化方法、后端和任务。
  新入口不重复设置 qm-atoms、active-atoms 或 mm-frozen；逐原子区域及冻结状态以坐标文件为准。
  QM/MM 区域与冻结状态是两个独立属性。首版只允许完整中央 QM 分子全部 frozen=F、
  MM 原子全部 frozen=T；其他冻结组合需在优化器支持后开放，不能静默忽略标记。
- 可选静电数据：单独引用 GDMA 原始 punch 输出及位点映射。扩展名不是识别依据。
  拓扑默认电荷与 GDMA 静电模型由显式开关选择，不能同时重复计入。

拓扑和力场文本语法、原生 GROMACS/LAMMPS 输入适配器另行确定，不能仅凭 .top/.ff
扩展名承诺兼容。共享力场可以被多个簇引用；参数缓存独立于坐标和量化设置。
新入口应保留当前 manifest 及 xyz+mm-charges+mm-vdw-file 入口的兼容性。

## Extended XYZ 与控制文件示意

Extended XYZ 的 Properties 语法使用通用格式，region、frozen、atom_id 和 coordinate_unit
的名称及语义是 cube 对扩展属性的约定，不是所有程序都自动识别的 QM/MM 标准。
atom_id 采用从 1 开始的稳定 ID，与拓扑逐项核验；不得与从 0 开始的原始轨迹索引混淆。
不依赖原子行的位置推断区域，也不因扩展名为 .xyz 就丢弃附加列。

```text
6
Properties=species:S:1:pos:R:3:region:S:1:frozen:L:1:atom_id:I:1 coordinate_unit=Angstrom
O  0.000  0.000  0.000  QM F 1
H  0.960  0.000  0.000  QM F 2
H -0.240  0.930  0.000  QM F 3
O  3.500  0.000  0.000  MM T 4
H  4.460  0.000  0.000  MM T 5
H  3.260  0.930  0.000  MM T 6
```

`input` 新接口草案（并非当前可执行语法）：

```text
$reorganization-energy
do = 1
coordinates = coordinates.extxyz
topology = topology.top
forcefield = forcefield.ff
qm-charge = 0
qm-multiplicity = 1
environment = qmmm
electrostatics = point-charge
mm-vdw = forcefield
qc-engine = gpu-pyscf
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
carriers = hole,electron
output = results/reorganization
$end
```

坐标文件包含提取时已选定的完整簇，运行中不根据截断距离再次隐式删改原子。

## GDMA 扩展边界

GDMA 输出是分布多极矩，可包含单极、偶极、四极及更高阶分量，不仅是部分电荷。
读取时必须保留单位、球谐归一化和分量顺序、展开中心、坐标框架及原始文件来源。
扩展位点模型不能强制假定每个多极位点等于一个原子；允许额外位点。

若将单分子参考几何的多极参数复用于簇中的多个分子，需有明确的位点/原子映射和
参考几何，并依据每个分子的旋转平移变换位点与各阶张量；周期平移不能错误旋转张量。
柔性构象的参数转移规则另行定义，不能把刚体拟合默认当作任意构象的精确参数化。

GDMA 数据不提供完整 LJ 力场，也不自动提供可极化模型所需的极化率或诱导求解器。
永久多极矩嵌入与可极化嵌入需分别实现和命名。

当前 PySCF/GPU4PySCF 路径只接通固定点电荷和 QM-MM LJ。未来高阶多极嵌入必须同时验证
QM Hamiltonian 和解析梯度的 CPU/GPU 实现，不能声称读入文件即已支持物理计算。
在高阶后端完成前，显式请求高阶多极嵌入应在启动优化前报不支持；任何只保留 rank-0
的近似必须由用户显式选择并在输出中记录，不能静默截断。

## 数据验证

坐标/拓扑数量、元素和顺序一致；每个原子/静电位点具有明确映射；单位显式可检查；
QM 核心与 MM 外部静电位点不重复；QM 和 MM LJ 类型全部有实际使用的成对参数。
复用计算核验各输入文件哈希、区域选择、静电模型及坐标框架，环境在四个点保持一致。

## 文献与官方接口

- GDMA 作者仓库：https://gitlab.com/anthonyjs/gdma
- Psi4 GDMA 接口、分量顺序和单位：https://psi4.github.io/psi4docs/master/gdma.html
- PySCF 固定点电荷 QM/MM：https://pyscf.org/user/qmmm.html
- CP2K 分离输入示例：https://manual.cp2k.org/trunk/methods/qm_mm/builtin.html
