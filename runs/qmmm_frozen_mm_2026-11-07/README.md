# 冻结 MM 的 QM/MM 四点法：实现与验证

已实现 `cube input` 的 `environment = qmmm`：MM 原子坐标、电荷在两次几何优化及所有单点中固定，只优化 QM 分子。CPU 和 RTX 5090 GPU 全流程实际通过，真空 DIPRO/四点法回归也通过。

## 物理模型

使用 PySCF/GPU4PySCF 的固定点电荷静电嵌入：QM 核–MM 点电荷、QM 电子–MM 点电荷作用同时进入能量和 QM 梯度。公式为：

```text
lambda(M) = [Eq(R0; M) - Eq(Rq; M)] + [E0(Rq; M) - E0(R0; M)]
```

M 在所有步骤完全相同；中性单重态、阳离子/阴离子双重态独立收敛。结果是固定环境下的 QM 几何松弛重组能，**不含 QM–MM van der Waals、跨边界成键、可极化 MM 或环境结构松弛的外球贡献**。固定 MM 自身能量抵消，未计算的 QM–MM van der Waals 能量一般不会抵消，当前明确省略。没有 Ewald 静电求和或频率分析。模型对应 [PySCF 官方静电嵌入说明](https://pyscf.org/user/qmmm.html)。

## 使用

仓库根目录 `input` 已将重组能段改成：

```text
$reorganization-energy
method = b3lyp
basis-set = 6-31G(d,p)
re-method = 4-points
environment = qmmm
mm-source = manifest
mm-cutoff = 12
mm-frozen = true
$end
```

GROMACS 模式从所选 MD 帧和 TPR 读取环境坐标、电荷，核验原始文件哈希与原子映射。以初始 QM–MM 最短原子距离 `< mm-cutoff` 选择**完整分子**及整体周期镜像，排除中央 QM 分子；将环境变换到对应 QM 坐标框架后冻结，不随优化重新筛选、居中或移动。cutoff 默认 12 Å，需要针对材料检查收敛。可在该段添加 `molecules = 1` 先运行一个目标分子。

也可指定 `xyz = molecule.xyz` 和 `mm-charges = environment.pc`，文件每行 `x y z q`（Å/元电荷），所有点与 QM 必须在同一坐标框架。该模式只接受一个 QM XYZ，包含文件中的全部点。示例见 [独立小分子 input](../../examples/quantum/qmmm/input)；详细字段与命令见 [控制说明](../../docs/cube-input.md#冻结-mm-的-qmmm-四点法)。

## 实际运行结果

测试 QM 区域为 F2，MM 为 4 个固定点电荷（对称、中性四极环境），输入见 [F2.xyz](F2.xyz)、[mm.pc](mm.pc)。这是接口和数值一致性测试，不是材料力场参数或实验精度的验证。B3LYP/6-31G(d,p)、def2-universal-jkfit 密度拟合、grid 3、SCF 阈值 1e-10。

| 重组能 / eV | 同机 CPU | RTX 5090 GPU | 绝对差 / eV |
|---|---:|---:|---:|
| 空穴 | 0.409628310571 | 0.409629871642 | 1.5611e-6 |
| 电子 | 5.758811614196 | 5.758802499930 | 9.1143e-6 |

中性和带电优化均收敛。差异低于预设 5e-3 eV 门限；不要求逐位一致。见 [summary.json](summary.json)、[CPU 四点能量及状态](cpu/reorganization.json)、[GPU 四点能量及状态](gpu/reorganization.json)。相同设置的真空结果为约 0.405515 / 5.706538 eV；加入非零 MM 环境后结果确实改变。

GPU 硬件为 RTX 5090、32607 MiB、compute capability 12.0、驱动 615.71.09。Python 3.12.15，PySCF 2.15.0，GPU4PySCF-cuda13x 1.9.1，CuPy-cuda13x 13.6.0，NumPy 2.5.4，geomeTRIC 1.1.1，SCF 使用 4 CPU 线程。版本、wheel、验证脚本与实际安装源文件哈希见 [environment.json](environment.json)，已核对服务器实际安装代码与工作区源文件一致。

CPU 流程单次观察 14.984 s，GPU 14.406 s；包括初始化、优化和多个单点，不包括之后单独的梯度差分。续算约 0.46 s。这不是重复计时或大体系性能基准。

## 梯度、冻结与失败处理

- 对 q=0、+1、−1 在同一 F2 几何做中心能量差分，步长 1e-3 Bohr。CPU/GPU 最大梯度绝对误差分别为 1.3453e-5 / 1.3448e-5 Hartree/Bohr，通过预设 1e-4 阈值。另在 HF 小体系用非对称点电荷验证三个电荷态梯度，以及共同平移 QM/MM 后能量不变。
- 优化向量只有 2 个 QM 原子，不含 4 个 MM 原子。所有优化评价和四个能量点均验证 MM 指纹一致；checkpoint 中另保存实际 MM 坐标、电荷。GPU MO 数组通过 CuPy 运行时检查，梯度继承链包含 `gpu4pyscf.qmmm.itrf.QMMMGrad`，确保包含 MM 力项。优化器由 CPU 控制。
- CPU/GPU `--resume` 后全部计算产物的哈希与修改时间保持不变。改变 MM 电荷再续算时必须失败，原结果不被覆盖，见 [GPU 拒绝复用日志](gpu/changed-mm-resume.log)。
- 限制最大优化步数为 1 的真实失败任务写出 failed 状态，不给出成功的重组能，见 [失败结果](failure/reorganization.json)。优化器初始化异常也会写失败状态。
- 嵌入 checkpoint 带独立的 `cube_reorganization` 标记，不能被当作真空 DIPRO 输入，即使设置 assume_converged 也拒绝。嵌入结果不能进行真空单体平移复用。
- 首轮发现带电 GPU 梯度继承 `GradientsBase`，没有继承限制性更强的 restricted `Gradients`；已修正检查。表中结果来自修正后的独立目录 v2，服务器保留首轮诊断结果。

52 项 Python 测试中 51 项通过、1 项显式 GPU opt-in 的 CuPy 代数测试跳过；本轮 GPU 量化全流程另外在 RTX 5090 实际运行。见 [单元测试日志](unit-tests.log)。旧真空单体/dimer、DIPRO、两种载流子的四点法、梯度和续算 CPU/GPU 回归全部通过，见 [回归汇总](vacuum-regression-summary.json)。安装后的 `cube` 命令 `--plan`、`--resume` 也通过，见 [入口日志](installed-entry-resume.log)。uv wheel/sdist 构建通过。

## 真实 MD 环境读取

对用户已有 `anneal_a` 最后一帧（frame 2000，20000 ps），构建第 1 个分子的 MM 环境：QM 为 86 原子，12 Å 判据选入 60 个完整邻近分子、5160 个 MM 原子，最近 QM–MM 距离约 1.85005 Å，中央 QM 分子已排除。TPR 电荷总和约 3.46e-6 e 的浮点残差原样保留，没有自动归一化。

这一项验证**环境构建与坐标映射**，没有启动该材料的完整 QM/MM SCF 或四点优化，也没有上传真实 MD 几何到服务器。这里只保存 [数量、源文件哈希与检查结果](md-environment-summary.json)，不提交真实环境坐标；本地快照保留在 `build/qmmm-anneal-frame2000/`。

## 复现与产物

```bash
uv sync --locked --extra workflow --extra quantum-cuda13
uv run --locked --extra workflow --extra quantum-cuda13 python scripts/validate_qmmm_workflow.py \
  --output build/qmmm-validation-new --gpu
```

使用新输出目录。原始运行位置为服务器独立测试目录；完整结果已取回 `build/qmmm-rtx5090/`。仓库保存小分子输入、环境快照、四点能量、优化坐标/轨迹和日志，不包含 wheel 或 checkpoint。归档日志仅去除 ANSI 颜色码和行尾空白，未改变计算数值；JSON 中的绝对路径保留了服务器计算来源。
