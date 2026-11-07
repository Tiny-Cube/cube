# cube 与 input 控制文件

`cube` 是统一调度入口，默认读取当前目录的 `input`。仓库根目录同时提供可执行文件 `cube` 和填写了本次实际数据路径的 `input`；安装后也提供同名命令。用户数据在仓库上一级，示例路径为 `../anneal_a.tpr` 等。其他机器上按实际路径修改，不要求重命名数据。

```bash
uv sync --locked --extra workflow
./cube input --plan                   # 解析/检查路径并显示任务顺序，不探测 GPU
./cube input --only get_dimer         # 只提取，本次真实数据已经实测
./cube input --resume                 # 执行其余阶段并复用一致且完整的成功结果
```

GPU 量化再安装可选 GPU 环境，原有 CUDA 13 验证条件继续适用：

```bash
uv sync --locked --extra workflow --extra quantum-cuda13
uv run --locked --extra workflow --extra quantum-cuda13 cube input --resume
```

仓库启动文件在没有激活其他虚拟环境时使用本仓库 `.venv/bin/python`；也可以在已激活环境中运行。只有提取时不会导入 PySCF/CuPy。`--plan` 只依赖标准库，不读取轨迹或启动量化。

## 控制语法

段落以 `$名字` 开始，以 `$end` 结束；字段为 `key = value`。支持 `#` 行尾注释；路径包含 `#` 时用引号。路径以 **input 所在目录**为基准。重复段落、重复字段、未知字段、缺失 `$end` 或不支持的计算引擎会报错。

```text
$get_dimer
tpr = ../anneal_a.tpr
gro = ../anneal_a_fix.gro
xtc = ../anneal_a.xtc
frame = -1
cutoff = 6  # 单位 Å
fragment-file = ../core.txt
output = mol
format = xyz,gjf
$end

$two-mol
method = b3lyp
basis-set = 6-31G(d,p)
transfer-integrals-method = DIPRO
$end

$one-mol
qc-engine = gpu-pyscf
method = b3lyp
basis-set = 6-31G(d,p)
$end

$reorganization-energy
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = vacuum
$end
```

执行顺序由依赖决定：提取→单体→dimer/转移积分→重组能；不会因为 `$two-mol` 写在 `$one-mol` 前面就先启动 dimer。只写所需段落即可。`--only` 可重复，依赖的提取/单体单点仍会执行；已有提取需要加 `--resume`。

## 提取定义和输出格式

- `frame` 从 **0** 开始；空值或 `-1` 表示最后一帧。GRO 用于核对 TPR 的原子/残基顺序，实际坐标来自所选 XTC 帧。
- `core.txt` 是以空白或逗号分隔的**分子局部、从 1 开始**的原子编号，支持 `#` 注释。当前一个 core 文件对应一致的分子类型/局部原子顺序。
- cutoff 单位为 **Å**，判据与提供的参考脚本一致：两组所选 core 原子的最短距离 **严格小于 cutoff**。这不是无限平面间距或质心距离；清单还记录所选 core 的质量中心距离。
- 用 TPR 键拓扑还原完整分子，以整体分子晶胞平移搜索接触，包括跨周期边界接触和多个有效镜像；支持非退化三斜盒。只枚举不同分子 ID 的组合，不包含同一分子的自镜像。
- 输出完整分子及侧链，不把 core 选择当成截断 QM 区域。当前要求电中性分子，量化任务采用闭壳层中性参考。
- `format = xyz` 只写 XYZ；`format = gjf` 只写 GJF；`format = xyz,gjf` 写两种格式。两种格式保持相同的坐标和原子顺序。另写 manifest/距离 CSV，不额外写未选的结构格式。
- GJF 的单体/二聚体 route 分别取 `$one-mol/$two-mol` 的泛函和基组；当前 GJF 为中性单重态冻结构型单点输入，不自动运行 Gaussian、Slurm 或 MOMAP。

例如真实输入最后一帧输出在：

```text
mol/frame2000/
  manifest.json                # 几何、原子编号、镜像、源文件/输出 SHA256
  dimers.csv                   # 每个有效镜像一行
  monomers/1.xyz, 1.gjf, ...
  dimers/1_8_img_0_-1_0.xyz     # 镜像保存在 ID 中，防止覆盖
  dimers/1_8_img_0_-1_0.gjf
  workflow.json                # 最近一次调用的任务状态
```

这份实际数据含 100 个 86 原子分子；最后一帧编号 2000，时间 20000 ps，6 Å 判据得到 550 个 dimer。与参考 `extract_dimers.py` 的分子对、镜像、距离逐项一致。完整原始轨迹不加入 Git。

## 单体与 dimer 量化

`qc-engine = pyscf` 使用 CPU；`gpu-pyscf` 使用 GPU4PySCF。`two-mol` 和重组能默认继承 `one-mol` 的引擎，也可各自明确设置；没有 one-mol 时默认 CPU。GPU 请求不静默回退 CPU。

`method` 为 RKS 泛函，或 `hf`；默认 `b3lyp`。`basis-set` 默认 `6-31G(d,p)`。可选公共字段：`auxbasis`（默认 def2-universal-jkfit，采用密度拟合）、`threads`（4）、`grid-level`（3）、`conv-tol`（1e-10）、`max-cycle`（100）、`max-memory-mb`（4000）、`output`。DIPRO 需要单体/dimer 的方法、基组、辅助基、网格、收敛阈值一致。

同一个真空单体只算一次。提取文件各自居中，pair 的两个片段与 dimer 在共同坐标框架内；对单体 checkpoint 只允许经逐原子检查的**整体平移**，保留 MO/Fock 并记录源 SHA256，不能通过该功能旋转、优化或改变单体内部构型。嵌入 Hamiltonian 不在这项平移复用的范围。

默认 `transfer-integrals-method = DIPRO`，也可逗号选择 `ESID`、`orbital_gmh`；ESID 还须明确 `assume-equivalent = true` 并通过适用性检查。`carriers = hole,electron` 可选择载流子。逐 pair 结果是 JSON；汇总 `qc/two-mol/transfer_integrals.csv` 的 VH/VL 单位为 meV。一维 DIPRO 值按固定显示轨道规范保留符号，多轨道子空间不压成标量，CSV 留空并保留 JSON 矩阵。

大体系可在 one-mol/重组能中用 `molecules = 1,2`，two-mol 中用 `dimers = 1_8_img_0_-1_0` 选择任务；未给选择时计算全部。各段还可用 `input = path/to/manifest.json` 读取此前提取，单体与重组能也可用 `xyz = A.xyz,B.xyz` 独立运行。数据源可独立于本次提取。

输出保存 `.chk`、`.log`、`.status.json`。SCF/优化失败会写失败状态，不把旧 checkpoint 当成成功。`--resume` 核对任务几何/设置/版本及成功输出哈希，优化 trace 或 checkpoint 改变时拒绝复用。更改任务物理设置应选择新的输出位置。

## Vacuum 四点法

`environment = vacuum` 使用孤立分子四点法，不含 QM/MM、溶剂或外球重组能。中性单重态优化得到 R0；阳离子（hole，q=+1）或阴离子（electron，q=−1）双重态优化得到 Rq。四点均使用相同泛函/基组和规定的电荷/自旋：

```text
lambda = [Eq(R0) - Eq(Rq)] + [E0(Rq) - E0(R0)]
```

MD 瞬时坐标只作为优化起点。geomeTRIC 在 CPU 控制优化；GPU 模式的每步能量/梯度来自 GPU4PySCF，并明确传回优化器。默认 `opt-maxsteps = 100`；SCF 或优化不收敛时不报告成功的重组能。保存优化 XYZ、逐步坐标/梯度、四个 Hartree 能量、两项松弛能及 eV 总和，另记录开壳层自旋平方。负松弛项超出 1e-5 eV 容差时拒绝给出成功结果，不取绝对值掩盖问题。未执行频率分析来证明所得驻点一定是全局极小值。

## 冻结 MM 的 QM/MM 四点法

设置 `environment = qmmm`：MM 原子的坐标、电荷在中性/带电优化和四个单点中完全相同，只优化整个 QM 分子。使用 PySCF/GPU4PySCF 的 **固定点电荷静电嵌入**，每步能量和 QM 核梯度均包含 QM 核–MM、QM 电子–MM 静电作用，GPU 引擎不回退 CPU。控制文件例子：

```text
$reorganization-energy
input = mol/frame2000/manifest.json
molecules = 1
output = reorganization-qmmm
qc-engine = gpu-pyscf
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = qmmm
mm-source = manifest
mm-cutoff = 12  # Angstrom
mm-frozen = true
carriers = hole,electron
$end
```

和 `$get_dimer` 一起运行时可以省略 `input`。仓库根目录 `input` 已改为此冻结 MM 模式；上面的 `molecules = 1` 用于先选一个分子，不填写则对所有选定单体逐一构建各自的环境。可执行：

```bash
./cube input --only reorganization-energy --resume
```

`mm-source = manifest` 从提取清单中找到原 TPR/GRO/XTC，核验文件哈希、所选帧、盒子、完整分子及原子编号后读取 TPR 电荷和 XTC 坐标。以**初始 QM 全分子与 MM 全分子的最短原子距离严格小于 mm-cutoff**选入完整 MM 分子，默认 12 Å。使用整体分子的三斜盒周期镜像，不截断分子；排除中央 QM 分子，非中央的周期复制分子可作为经典环境。MM 坐标随 QM 的初始居中平移一次，之后完全冻结，不跟随优化后的 QM 重新平移或重新筛选。应按材料检查 cutoff 收敛；这是有限簇静电环境，不是无限周期 Ewald 求和。

也可以使用一个独立 QM XYZ 和外部 MM 点电荷文件：

```text
$reorganization-energy
xyz = molecule.xyz
qc-engine = gpu-pyscf
environment = qmmm
mm-charges = environment.pc
mm-frozen = true
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
$end
```

点电荷文件每行四列 **x y z q**，坐标单位 Å，电荷单位元电荷；支持空行和 `#` 注释，无原子数行或表头。所有点与 QM XYZ 必须在**同一个坐标框架**中，且不能包含中央 QM 原子自身的力场电荷。文件模式只允许一个 QM XYZ，全部所给 MM 点均纳入，不应用 `mm-cutoff`。可将 LAMMPS 等程序的数据导出为该格式；当前未提供 LAMMPS 力场/轨迹自动读取。

```text
# x y z q; Angstrom/e
3.0 0.0 0.0  0.1
0.0 3.0 0.0 -0.1
```

四点使用相同环境 M：

```text
lambda(M) = [Eq(R0; M) - Eq(Rq; M)] + [E0(Rq; M) - E0(R0; M)]
```

结果保存 `mm-environment.json`（完整固定坐标、电荷、来源、镜像映射及指纹）、中性/带电优化轨迹、四个能量、带嵌入标记的 checkpoint。每步核验实际 SCF 中的 MM 坐标和电荷，优化向量只有 QM 原子。MM 内容或坐标改变时，`--resume` 拒绝复用已有结果；修改输出快照也会使续算失败。嵌入 checkpoint 不允许进入真空 DIPRO 或真空平移复用路径。

**模型范围**：这是固定电荷静电嵌入下的 QM 几何松弛重组能。未包含 QM–MM van der Waals、跨边界成键/链接原子、可极化 MM 或 MM 结构松弛，因此不是包含全部力场项的 additive QM/MM，也没有环境结构松弛的外球重组能。固定 MM 自身能量为常数，在差值中抵消；QM–MM van der Waals 一般不为常数，当前是明确省略，不能称为抵消。静电模型见 [PySCF 官方说明](https://pyscf.org/user/qmmm.html)。

小分子示例和实际验证入口：

```bash
./cube examples/quantum/qmmm/input --plan
uv run --locked --extra workflow python scripts/validate_qmmm_workflow.py \
  --output build/qmmm-validation-new
# 兼容 GPU 上再加 --extra quantum-cuda13，并在脚本末尾加 --gpu
```

CPU 与 RTX 5090 GPU 的实际四点法、解析梯度和固定 MM 检查结果见 [QM/MM 验证报告](../runs/qmmm_frozen_mm_2026-11-07/README.md)。

## 独立脚本和验证

`scripts/get_dimer.py`、`one_mol.py`、`two_mol.py`、`reorganization_energy.py` 可分别读取相同 input，只执行对应阶段及其依赖。它们与 cube 调度器共用任务实现；旧 `run_pyscf.py` / `calculate_transfer.py` 两阶段入口保留。

```bash
uv run --locked --extra workflow python scripts/get_dimer.py input
uv run --locked --extra workflow python scripts/validate_cube_workflow.py --output build/cube-workflow-cpu
# RTX 5090 等已配置 GPU 环境
uv run --locked --extra workflow --extra quantum-cuda13 python scripts/validate_cube_workflow.py \
  --output build/cube-workflow-gpu --gpu
```

本次真实 MD 提取与 CPU 验证记录见 [cube input 报告](../runs/cube_input_2026-11-07/README.md)；新增 GPU 梯度/四点法及安装后的 `cube` 命令已在 RTX 5090 实测通过，见 [GPU 工作流报告](../runs/cube_input_gpu_2026-11-07/README.md)。GPU 报告使用乙烯/F2 小体系验证接口与数值一致性，不代表已完成真实 MD 的全部量化任务，也不是大分子加速比测试。
