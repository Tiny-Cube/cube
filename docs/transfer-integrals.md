# 转移积分：理论、方法和 CPU/GPU 实现

转移积分是局域电荷态之间的电子耦合。它是能量量纲，不是跃迁速率。多电子定义为

\[
V_{DA}=\langle\Psi_D|\hat H|\Psi_A\rangle.
\]

实际计算必须同时说明局域态、Hamiltonian、电子数、自旋、基组和环境模型。当前实现提供四个明确区分的方法；中性闭壳层、冻结分子构型是已实现的计算范围。空穴默认取单体 HOMO 子空间，电子取 LUMO 子空间。这是局域单电子近似，不是任意材料中的精确多电子耦合。

## 1. DIPRO：收敛二聚体算符的片段投影

对单体 A、B 分别做 SCF，将单体轨道系数按原子/AO 映射嵌入完整二聚体 AO 基，组成矩阵 C。以收敛二聚体 Fock/KS 矩阵 F 和 AO 重叠 S_AO 构造

\[
H=C^TFC,\quad S=C^TS_{AO}C,\quad H_{\mathrm{eff}}=S^{-1/2}HS^{-1/2}.
\]

最后取正交化 Hamiltonian 的 A–B 非对角块。对于每个片段各一个实轨道、且 S_AA=S_BB=1 的情况，得到

\[
V_{AB}=\frac{J_{AB}-\tfrac12 S_{AB}(\epsilon_A+\epsilon_B)}{1-S_{AB}^2}.
\]

这里 J_AB 是未正交化的算符矩阵元，不能直接拿它替代有效 V。能量原点变换 H→H+cS 后，V 应保持不变。实现使用矩阵对称正交化，也支持多轨道子空间。重叠奇异或病态时拒绝计算，不静默删去基矢。

理论依据：Baumeier、Kirkpatrick、Andrienko，[PCCP 12, 11103–11113 (2010)](https://doi.org/10.1039/C002337J)，Sec. II，Eq. (7)；[作者机构记录](https://research.tue.nl/en/publications/density-functional-based-determination-of-intermolecular-charge-t/)。二态表达式也可在 [FHI-aims 官方手册 Eq. (3.259)](https://www.fhi-aims.org/uploads/manual/Ch3/S48.html)核对。该手册还讨论了不同 FO-DFT 电子数选择，不能把它们都称为同一个方法。

当前基组约定是**孤立单体自身基组**，然后零填充嵌入二聚体 AO 基；不使用另一片段的 ghost/counterpoise 基函数。完整二聚体 MO 展开另行重建投影矩阵，用残差检查 AO 映射和 SCF 一致性。所有 MO 都参与该诊断，不只取两个前线 MO。

## 2. Frozen density：中性单体密度叠加近似

保持相同单体参考，改用

\[
P=P_A\oplus P_B,\quad F_{\mathrm{frozen}}=F[P],
\]

只构造一次二聚体 Fock/KS 算符，再投影与正交化，不对这个二聚体密度进行 SCF。只请求此方法时，程序不运行二聚体 SCF。它与 DIPRO 的区别是 Hamiltonian 的密度来源，而非换一种矩阵乘法得到同一个定义。

这里明确采用 AO 密度矩阵叠加。叠加密度一般不是由正交化全部占据轨道构成的幂等行列式密度；不能混同二者。它也不是带电片段、occupation-reset、cDFT 或完整 frozen-density embedding 的实现。

Baumeier 等讨论了 noSCF 近似；Schober、Reuter、Oberhofer 的[原始论文](https://doi.org/10.1063/1.4940920)及[作者公开全文](https://arxiv.org/html/1512.00200)系统区分了电子数与片段态的选择，见预印本 Sec. II.3.1–II.3.3。当前代码不声称实现了该论文的所有 FO-DFT 变体。

## 3. ESID：有限适用域的半能级劈裂

正交两态 Hamiltonian 的本征值差满足

\[
\Delta E=\sqrt{(\epsilon_A-\epsilon_B)^2+4V^2}.
\]

仅在等价位点、适当的孤立两态近似下，才有 |V|=ΔE/2。程序要求显式 `--assume-equivalent`，还检查投影位点能之差和二聚体前线态对子空间的覆盖。这个选项是物理假设，不能绕过数值检查。非对称位点或多轨道前线子空间返回 `not_applicable`。

能级来自实际二聚体 SCF：空穴候选为 HOMO/HOMO−1，电子候选为 LUMO/LUMO+1；候选仍须通过子空间检查。两态假设不成立时不输出误导性的半能差。

理论与限制：Cave、Newton，[JCP 106, 9213–9226 (1997)](https://doi.org/10.1063/1.474023)，Sec. I；[作者公开记录](https://scholarship.claremont.edu/hmc_fac_pub/8/)。Valeev 等的 [JACS 128, 9882–9886 (2006)](https://doi.org/10.1021/ja061827h)强调极化导致的位点能差会使能级劈裂法失效；见[作者机构记录](https://experts.arizona.edu/en/publications/effect-of-electronic-polarization-on-charge-transport-parameters-/)。

## 4. 二态 GMH 与单电子轨道 GMH

对正交绝热态的能量和投影到指定方向的偶极矩阵，二态 GMH 为

\[
|V|=\frac{|\Delta E\,\mu_{12}|}{\sqrt{(\mu_{22}-\mu_{11})^2+4\mu_{12}^2}}.
\]

独立的 `gmh()` 代数函数接受两态能量和三维偶极矩阵。默认沿差偶极方向投影；差偶极因对称性消失时使用跃迁偶极方向。两者都为零时定义不足，程序拒绝给出数值。能量和偶极单位各自必须一致；返回与输入相同的能量单位。公共的对角偶极偏移抵消，因此结果应对坐标原点和刚体旋转不变。

GPU 工作流的 `orbital_gmh` 使用所选 KS/HF 轨道的能量与单电子偶极积分，**是轨道 GMH 近似**。它不生成多电子激发态或电荷转移态，也不把 KS 能差说成多电子激发能。要做多电子 GMH，须提供相应真实绝热态与跃迁密度/偶极，不能直接套用当前轨道工作流。

原始方法：Cave、Newton，[Chemical Physics Letters 249, 15–19 (1996)](https://doi.org/10.1016/0009-2614(95)01310-5)，[作者记录](https://scholarship.claremont.edu/hmc_fac_pub/816/)；1997 年上述论文作进一步比较。多态效应可能使二态结果不足，不能因两个方法输出接近就认为两态假设必然正确。

## 实现与安装

代码在 `src/universal_kmc/quantum/`，不改变通用 KMC 核心。方法旁的注释包含原论文 DOI、公式/章节和短原文摘句。发行物不包含论文全文或上游程序源码。

CPU 环境：

```bash
uv sync --locked --extra quantum
uv run --locked --extra quantum cube-transfer \
  --fragment-a examples/quantum/ethylene_a.xyz \
  --fragment-b examples/quantum/ethylene_b.xyz \
  --backend cpu --basis def2-svp --xc b3lyp \
  --method dipro --method frozen_density --output build/transfer-cpu.json
```

已实测的 CUDA 13 环境（需兼容的 NVIDIA 驱动，安装运行库时会下载较大的 CUDA wheel）：

```bash
uv sync --locked --extra quantum-cuda13
uv run --locked --extra quantum-cuda13 cube-transfer \
  --fragment-a examples/quantum/ethylene_a.xyz \
  --fragment-b examples/quantum/ethylene_b.xyz \
  --backend gpu --basis def2-svp --xc b3lyp \
  --method dipro --method frozen_density --method esid --method orbital_gmh \
  --assume-equivalent --output build/transfer-gpu.json
```

默认 XYZ 坐标为 Å；只有明确使用 Bohr 坐标时指定 `--unit Bohr`。`--xc hf` 选择 RHF，否则使用 RKS。CPU/GPU 都显式采用相同的密度拟合辅助基与收敛条件，完整设置写入 JSON。Python API 可用 `SCFConfig` 修改辅助基/精度，并通过 `orbital_indices={"hole": ([6,7], [6,7])}` 指定零基 MO 子空间。

GPU SCF、Fock 和矩阵投影使用 GPU4PySCF/CuPy；部分 AO 重叠/位置偶极积分由 CPU libcint 提供，小结果和 JSON 写出在 CPU 上。GPU 请求不静默转 CPU。返回的 `scf.*.class`、`gpu_runtime` 和实际设备名用于核查执行路径。

**本次 CUDA 13 CuPy wheel 缺少 cuTENSOR Python binding**；虽然单独安装 cuTENSOR 共享库不能补齐该接口。经实测采用 CuPy GPU 张量收缩，并在结果中记录 `tensor_contraction_engine="cupy"`；不声称使用了 cuTENSOR，也不把上游库告警吞掉。锁文件不安装这个环境下用不到的 cuTENSOR 包。未承诺其他 CUDA/显卡组合已通过同样测试。

## 输出、失败与验证

- JSON 原始 Hamiltonian 用 Hartree，有效耦合同时输出 eV/meV；一维单态返回标量 |V|，多轨道返回矩阵、奇异值和 Frobenius 范数。后者是诊断，不是自动用于任意多态速率的等效 V。
- 保存冻结坐标、参数、输入哈希、轨道索引、重叠、SCF 收敛、实际后端、依赖版本和同步墙钟时间。原始耦合的正负号依赖轨道相位；KMC 常用 |V|²，多态比较应考虑子空间规范。
- 任一 SCF 不收敛会失败；不输出伪有效耦合。CLI 成功退出 0，存在不适用方法退出 2，并保留逐方法原因；运行失败退出 1，输出明确的失败 JSON，替换同路径的旧成功结果。
- 本模块尚不实现 QM/MM、周期量化、带电开壳层或真正的多电子 GMH。这些不是当前结果的隐含物理条件；已有 QM/MM 规划仍是独立开发工作。

GPU 服务器完整复现命令：

```bash
uv run --locked --extra quantum-cuda13 python scripts/validate_transfer_gpu.py \
  --output build/transfer-validation
```

见[实际 GPU 验证报告](../runs/transfer_gpu_2026-11-06/README.md)与机器可读结果。对不同方法之间的差异，检查物理定义与近似；CPU/GPU 对照只在**同一方法和相同设置**之间判断数值一致性。
