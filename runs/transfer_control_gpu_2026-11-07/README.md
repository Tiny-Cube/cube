# RTX 5090：XYZ → PySCF 输出 → 控制文件 → 转移积分

2026-11-07 完成两阶段实际 GPU 验证。测试使用安装后的 wheel，已核对安装目录中的全部量化模块与上传源码 SHA256 一致，wheel SHA256 与本地构建一致。第一阶段实际分别运行 CPU/GPU SCF，第二阶段分别使用 CPU/GPU 读取计算结果，没有复用旧的量化结果。

硬件为 NVIDIA GeForce RTX 5090，32 GB，驱动 615.71.09，compute capability 12.0；Python 3.12.15、PySCF 2.15.0、GPU4PySCF 1.9.1、CuPy 13.6.0。CUDA runtime API 13.0，driver API 13.4，安装包详情见 `environment.json`。GPU4PySCF 明确提示使用 CuPy GPU 张量收缩；本测试没有把它称为 cuTENSOR。

## 计算与结果

仓库乙烯冻结二聚体，B 沿 z 平移 3.5 Å；中性闭壳层 RKS/B3LYP，def2-SVP、def2-universal-jkfit 密度拟合，grid level 3、conv_tol 1e-10、4 线程。测试 DIPRO、满足等价性/子空间条件的 ESID、单电子 orbital-GMH，两种载流子均计算。

GPU SCF 后使用 GPU 读取 checkpoint 得到：

| 载流子 | DIPRO / eV | ESID / eV | orbital-GMH / eV |
|---|---:|---:|---:|
| 空穴 | 0.422000342118 | 0.421957004673 | 0.421957004673 |
| 电子 | 0.485180385581 | 0.479922662991 | 0.479922662991 |

与同设置的独立 CPU 直接计算相比，所有载流子/方法的最大耦合差异约 **3.56e-7 eV = 0.356 µeV**，通过 `atol=2e-5 eV, rtol=2e-4` 的验证门限。A、B、AB 的 CPU/GPU SCF 总能最大差异约 **8.81e-13 Hartree**。这些是本示例的数值一致性检查，不是任意材料的物理精度或性能结论。

四种交叉组合全部成功，结果均 `status=completed`、`scf_executed=false`：

| 量化结果来源 | 后处理 | 与独立 CPU 直接参考的最大耦合差异 / eV |
|---|---|---:|
| CPU SCF | CPU | 约 5e-14 |
| CPU SCF | GPU | 约 5e-14 |
| GPU SCF | CPU | 约 3.56e-7 |
| GPU SCF | GPU | 约 3.56e-7 |

分离 SCF 差异之后，**同一份 checkpoint** 的 CPU/GPU 后处理差异更小：CPU 生成的 checkpoint 最大 **1.11e-15 eV**，GPU 生成的 checkpoint 最大 **2.22e-16 eV**。

## 执行路径和附加检查

- A/B/AB 的 GPU SCF 实际类均为 `gpu4pyscf.df.df_jk.DFRKS`；CPU 类为 `pyscf.df.df_jk.DFRKS`。每个状态 JSON 保存实际类、设置、能量和收敛状态。
- GPU 后处理记录 RTX 5090 和 `projection_array_module="cupy"`；额外拦截共享投影函数，确认输入的 AO 重叠及 Hamiltonian 均为真实 CuPy 数组。
- 同时禁止 CPU/GPU 两个 SCF kernel 后，GPU 后处理仍成功。第二阶段没有重新执行量化。
- dimer 原子交错、从 1 开始的显式 MO、原生 dimer checkpoint 的完整 MO 算符恢复均在 GPU 后处理上通过；与各自 CPU 后处理的最大差异分别约 1.11e-15、1.11e-15、1.67e-15 eV。
- 标准 checkpoint 缺少收敛确认、明确未收敛、破坏的 MO 正交性、截断轨道、错误坐标、不同方法和不一致 Fock 均被拒绝。ESID 未声明等价位点时退出 2，保留逐方法原因。
- CPU/GPU 第一阶段均验证了 `max_cycle=1` 不收敛路径：退出 1，生成失败状态，移除同名旧 checkpoint。
- 量化单元测试 **13 项全部通过，无跳过**，包括真实 GPU Löwdin 对照。日志为 `unit-tests.log`。

checkpoint 读取、AO 积分和完整 MO 特征算符恢复在 CPU；SCF 使用 GPU4PySCF，投影及正交化使用 CuPy。测试范围仍是中性闭壳层、冻结构型，尚不包含 QM/MM 或 GPU KMC，也没有进行通用 GPU 加速比基准测试。

## 复现与证据

```bash
uv sync --locked --extra quantum-cuda13
uv run --locked --extra quantum-cuda13 python scripts/validate_transfer_control.py \
  --output build/transfer-control-gpu-validation --gpu
CUBE_GPU_TESTS=1 uv run --locked --extra quantum-cuda13 python -m unittest discover \
  -s tests -p 'test_quantum*.py' -v
```

`summary.json`、`environment.json`、各后处理 JSON、`cpu/` 与 `gpu/` 中的 SCF 日志/状态构成本报告证据。日志仅去除行末空白。原始 HDF5 checkpoint 已下载到本机 `build/transfer-control-gpu-validation/`，不提交二进制输入。`capture_environment.py` 是本次环境核对脚本。源文件 SHA256 可与仓库核对；本次只增强验证脚本，生产计算代码无需修复。
