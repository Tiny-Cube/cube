# 转移积分的 RTX 5090 实测

2026-11-06，在用户授权的独立 GPU 测试目录执行。**八类 CPU/GPU 量化案例全部通过**；最终验收使用 wheel 安装后的代码，而不是仅从源码目录导入。`summary.json` 的 `all_passed=true`、安装位置标志与源码哈希可核对该结论。

## 环境与数值结论

- RTX 5090，32607 MiB，驱动 615.71.09，计算能力 12.0；Python 3.12.15，uv 0.10.2。
- PySCF 2.15.0、GPU4PySCF 1.9.1、CuPy CUDA13 13.6.0、NumPy 2.5.4。依赖由根目录 `pyproject.toml` 和 `uv.lock` 配置，通过 uv 创建隔离环境；精确同步后移除了早期探测安装的无用 cuTENSOR 包。
- 实际 GPU 运行时报告 CUDA 13.0，驱动接口报告 13.4；完整已安装包版本见 [environment.json](environment.json)。CUDA 运行库由 uv 安装的 wheel 提供，本机不需要另行用 nvcc 编译 GPU4PySCF。
- CPU/GPU 最大 SCF 总能量差为 **5.730×10⁻¹¹ Hartree**。
- 同一方法 CPU/GPU 的最大耦合比较误差为 **6.231×10⁻⁷ eV**，约 **0.6231 μeV**。多轨道块比较奇异值，避免把简并子空间内的基变换误当成物理误差。
- 验收门限预先定义为总能量绝对差 10⁻⁶ Hartree，耦合差 `atol=2e-5 eV, rtol=2e-4`，实际误差更小。该误差是后端数值差，不是对实验或高等级多电子方法的物理误差估计。
- 最终完整驱动约 93.10 s。计时为同步墙钟时间，包含 SCF、参数后处理及验证任务；不是重复测量的性能基准，不据此承诺普遍 GPU 加速倍数。

当前 CUDA13 CuPy wheel 没有 GPU4PySCF 所用的 cuTENSOR Python binding，因此实际张量收缩引擎是 **CuPy GPU**。没有静默回退 CPU，也没有声称用到 cuTENSOR。每个 GPU SCF 记录了 `gpu4pyscf.*` 类，代码检查了设备端 MO 数组；投影矩阵使用 CuPy。AO 位置偶极积分使用 CPU libcint，这一边界也写入结果。

## 八类案例

| 案例 | 要验证的行为 | 结果 |
|---|---|---|
| cofacial | 3.5 Å 共面平行乙烯二聚体；四种方法、空穴/电子 | 通过 |
| slipped | 横向错位的乙烯；同一方法 CPU/GPU 对照 | 通过 |
| distorted | 一片段形变、位点不等价；拒绝不适用的 ESID | 通过 |
| separated | 8 Å 间距；DIPRO 耦合相对近距离下降超过一个数量级 | 通过 |
| hf | RHF/STO-3G，独立于 B3LYP/def2-SVP 路径 | 通过 |
| multiorbital | 显式两轨道/片段的 2×2 耦合块；拒绝两态方法 | 通过 |
| benzene | 自动识别苯的近简并 HOMO/LUMO 子空间，各返回 2×2 块 | 通过 |
| rigid_transform | 整个二聚体共同旋转和平移；耦合大小不变 | 通过 |

几何为公开、可复现的小分子验证构型，未使用原 e_006 的原始 MD 结构，不把此结果称为已经验证整个 MD/QM-MM/KMC 流水线。

## 同一几何上不同方法的示例

共面乙烯、3.5 Å、B3LYP/def2-SVP、相同密度拟合辅助基，GPU 结果约为：

| 方法 | 空穴 | 电子 |
|---|---:|---:|
| DIPRO | 422.000342 meV | 485.180386 meV |
| Frozen density | 478.127123 meV | 484.594097 meV |
| ESID | 421.957005 meV | 479.922663 meV |
| Orbital GMH | 421.957005 meV | 479.922663 meV |

不同方法采用的 Hamiltonian/局域态近似不同，**不要求彼此相等**。对称两态极限下 ESID 与 orbital GMH 接近，也不意味着对非对称或多态体系可互换。`orbital_gmh` 明确是单粒子近似，不是多电子绝热态 GMH。

## 其他已执行检查

- 8 项代数/输入/GPU 矩阵测试全部通过，0 项跳过，见 [algebra-tests.log](algebra-tests.log)。覆盖二态解析式、能量原点、相位/置换、广义本征谱、奇异度量、已知局域 Hamiltonian 的 GMH 恢复、偶极原点/空间旋转不变性和非法输入。
- 完整二聚体 MO 展开重建 DIPRO 投影，残差低于实现的 10⁻⁶ Hartree 门限。
- GPU frozen-density 单独运行只做两个单体 SCF，没有二聚体 SCF；其结果与四方法一起请求时一致到约 10⁻¹⁵ eV。
- CPU/GPU 的 SCF 不收敛案例均被拒绝。
- CLI 正常运行退出 0，未确认 ESID 等价假设时退出 2；SCF 失败退出 1，并以失败 JSON 替换同路径旧成功结果，见 [cli-failed.json](cli-failed.json)。
- 本地既有四语言构建、Rust/C++/Fortran 测试、19 项核心 Python 测试及 uv 打包通过。默认无 NumPy 的核心环境跳过可选量化测试；本目录服务器验收单独实际执行了全部 8 项量化测试。另用有 NumPy 的系统 Python 执行了两项既有分子分析测试。

## 复现与产物身份

从仓库根目录运行：

```bash
uv sync --locked --extra quantum-cuda13
uv run --locked --extra quantum-cuda13 python scripts/validate_transfer_gpu.py \
  --output build/transfer-validation
```

本次最终验收还检查了 wheel 安装：本地 `uv build --wheel`，服务器精确 `uv sync --locked --extra quantum-cuda13`，随后用 `uv pip install --python .venv/bin/python --no-deps --reinstall <wheel>` 安装到该环境，再由 `.venv/bin/python` 运行驱动，避免 `uv run` 自动切回 editable 安装。

[summary.json](summary.json)确认已安装模块与测试源码逐文件一致，`installed_from_source_tree=false`。uv 的 wheel `direct_url` 元数据没有提供归档哈希，故额外独立计算了 [wheel-artifact.json](wheel-artifact.json)；其校验和与本地交付 wheel 一致。全部量化模块、验证脚本、测试文件、依赖声明和锁文件的哈希也在 summary 中。

每个 `*-cpu.json` / `*-gpu.json` 保存实际输入几何、方法设置、重叠/耦合矩阵、SCF 记录、设备/版本和计时。理论与方法适用域见[转移积分文档](../../docs/transfer-integrals.md)，文献访问范围见[来源索引](../../docs/research/transfer-integral-sources.json)。发行源码只包含自行编写的实现和短引用，不附带论文全文、外部程序源码或登录凭据。
