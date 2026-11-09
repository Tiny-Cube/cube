# Charge-carrier：控制文件到静态 KMC

本分支已接入主模块内部选择、参数网络导入、量化结果组网、Marcus 速率、独立单载流子 KMC、扩散／迁移率张量与收敛诊断。当前协议是**静态快照、零电场、经典非绝热 Marcus、独立单载流子**；没有材料 GPU KMC、多粒子相互作用、有场漂移或时间依赖形貌。

## 主模块

```text
$charge-carrier
do = 1
get_dimers_monomers = 1
one-mol = 1
two-mol = 1
reorganization-energy = 1
kmc = 1
$end
```

字段决定本主流程调用哪些子模块；1 请求，0 跳过，省略保留旧自动调度。相应模块仍需自己的设置段。关闭前处理后，下游必须有明确的已有输入；配置缺失在科学计算前报错。子模块自身 `do=1` 是独立请求，因此主模块的 `one-mol=0` 不会取消独立单体请求。子模块明确 do=0 而主流程显式请求 1 时报告冲突。根 input 的子模块 do 被注释，便于直接用主段控制流程；独立执行时再启用相应 do。

主模块只请求参数准备、没有 KMC 任务时，仍报告 prepared / kmc_not_requested，退出码 2。只有 KMC 与指定统计诊断完成才报告载流子 completed；欠采样报告 partial。不会把量化准备成功当成迁移率完成。

## 不运行量化的完整示例

```bash
uv sync --locked --extra transport
uv run --extra transport cube examples/charge-carrier/input --plan
uv run --extra transport cube examples/charge-carrier/input
uv run --extra transport cube examples/charge-carrier/input --resume
```

该例将四个前处理开关置 0、kmc 置 1，读取合成一维周期链。空穴／电子各 512 条轨迹，和解析 D=k a² 对照。参数没有材料测量含义。需要完整 MD/量化链时使用 workflow 和对应 GPU 量化 extra。

`$kmc` 可在三个主模块全关时通过自身 `do=1` 独立运行，沿用同一实现。

## KMC 设置

```text
$kmc
# do = 1
input = charge-network.json
output = results/charge-carrier
backend = python
carriers = hole,electron
temperature = 300
duration = 1e-9
trajectories = 128
samples = 51
seed = 42
max-events = 1000000
fit-start = 0.5
relative-tolerance = 0.25
trace-trajectories = 3
exact-max-sites = 512
$end
```

温度 K，时间 s；每条轨迹相同物理时长。samples 至少 5，fit-start 在 0 与 0.9 之间；小于四条轨迹可以产生轨迹，但不提供可靠统计区间，报告 partial。每条轨迹达到 max-events 而时间未到时记录截断，不把余下位移固定填满来伪造完成。relative-tolerance 是统计区间范数与时间窗变化的诊断目标，不能代替时间、轨迹数、尺寸和形貌收敛研究。

input 与 manifest/couplings/reorganization 两种来源互斥。材料 KMC 目前 backend 只接受 python；四语言通用 KMC1 核心仍可独立使用。

## 参数网络格式

见 [完整合成网络](../examples/charge-carrier/chain.json)。`schema=cube-charge-network-1`：

- `lattice_angstrom`：3×3，行向量为 a、b、c，有限右手非退化晶胞。
- sites：唯一 id、材料坐标中的 position_angstrom、各载流子的 free_energy_ev。
- pairs：唯一 id、source/target、三个整数 image、各载流子的 coupling_ev 与 reorganization_ev。

每个 pair 是一个无向周期通道，程序生成正反两条边；不要重复提供反向 pair。允许同一对位点不同镜像与自镜像。位移 `r_j-r_i+image@cell` 以 m 累积。缺失耦合不是 0；非连通正速率网络需先显式选择连通组分，避免未定义初态组分权重。当前不自动删除位点或补耦合。

位点能是该载流子形成的自由能**模型**，不是自动读取的 HOMO/LUMO 能。显式 uniform 必须在网络中写零值，或在量化组网中指定 uniform，不能依赖隐含填零。成对 lambda 正反共用；Marcus 速率满足所采用位点能模型的详细平衡。参见 [VOTCA 的输运理论](https://www.votca.org/xtp/transport_theory.html) 与本仓库 [符号契约](research/transport-theory.md)。

## 本次／已有量化结果到网络

完整链默认使用本次提取 manifest、two-mol 的 transfer_integrals.csv 和 reorganization-energy 的完整四点 JSON。跳过这些前处理时可显式指定：

```text
$kmc
manifest = mol/frame2000/manifest.json
couplings = mol/frame2000/qc/two-mol/transfer_integrals.csv
reorganization = reorganization-index.json
site-energy-model = four-point-adiabatic
pair-reorganization-model = four-point-equivalent
output = results/charge-carrier
$end
```

reorganization-index.json 内容为 `{"1":"path/1/reorganization.json","2":"path/2/reorganization.json"}`，其中路径相对该索引文件。couplings 使用 cube 输出列 dimer、mol_i、mol_j、VH_meV、VL_meV；dimer 与镜像必须逐项对应 manifest。所有所需通道／位点参数必须完整。

site-energy-model 可选 uniform（明确忽略位点无序）或 four-point-adiabatic。后者使用 `[Eq(Rq)-E0(R0)]` 的电子形成能近似，空穴为中性→阳离子，电子为中性→阴离子，单位转换为 eV；不包含熵和振动自由能，也不保证已有电子结构方法适合目标材料。

自动成对重组能采用 `charged_relaxation_i + neutral_relaxation_j`，并核对反向表达式在 rtol=1e-4、atol=1e-6 eV 内一致。库内允许的微小优化数值误差不取绝对值。一般不等价分子／环境未通过该检查时，要求在显式网络中给出物理定义的成对 lambda；不会直接平均两个总重组能。冻结 MM 结果只包含该模型中的 QM 松弛，不推断 MM 外球松弛。

新提取清单保存 transport_position_angstrom。旧单体已居中，不能用其 XYZ 均值重建材料位点；仅在原始 MD 的来源哈希一致时重读原坐标。为避免误用旧缓存，提取指纹版本更新为 3；可使用新输出目录重新提取，或跳过提取并指定已有清单及可核实 MD 来源。

多帧本次流程逐帧匹配量化结果，每帧独立 KMC，输出到 frameNNNN。不能把多个快照当作同一条动态轨迹拼接。外部单个耦合/重组索引配多帧不自动复用，须逐帧运行或由本次工作流提供逐帧结果。

## 输出和解释

result.json 保存输入／实现哈希、模型范围、轨迹事件数、正反速率诊断、周期输运秩、张量与收敛状态。另写：

- network.json：归档参数网络；`*_rates.json`：定向速率与 eV ΔG。
- `*.kmc`：SI 单位 KMC1 图，可供原四语言通用程序读取。
- `*_traces.json`：指定前几条轨迹，包含事件、位点、未折回位移及终止观测。
- `*_samples.npz`：等物理时间采样；截断后未知位置为 NaN，不用于迁移率拟合。
- `*_observations.csv`：均值和中心化协方差。

有限时间扩散由协方差斜率的一半得到，使用轨迹间 jackknife 给出标准误；迁移率张量 `mu[cm²/(V s)]=D[m²/s]*1e4/(kB[eV/K]*T)`，只用于这里的零场 Einstein 协议。没有把三维 trace/3 当作所有低维材料的标量迁移率。保留原始有限窗估计，不把负的噪声特征值静默截为零。

小图另给可逆 CTMC corrector 的渐近参考，支持非均匀 Boltzmann 平衡、自镜像及多个平行通道；超过 exact-max-sites 或参考矩阵病态时明确跳过。参考误差、KMC 统计误差和材料模型误差是不同问题。周期回路平移秩为 0 的封闭簇长期 D=0；其短时间非零 MSD／拟合斜率是有界运动，不能解释为材料迁移率。

`--resume` 校验输入和实现指纹以及所有输出文件哈希，复用相同的已结束 completed/partial 结果。它不延长已算轨迹、不恢复中断轨迹；改变 duration、轨迹数或实现后使用新输出目录。partial 表示需要检查或增加采样，退出码 2，不代表获得了收敛材料迁移率。
