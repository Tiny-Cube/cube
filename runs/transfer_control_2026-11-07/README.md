# 两阶段 PySCF → 转移积分 CPU 验证

验证对象：`scripts/run_pyscf.py` 接受单体 A/B 和 dimer AB 的 XYZ，生成 `.chk/.log/.status.json`；`scripts/calculate_transfer.py transfer.toml` 按控制文件读取结果，计算转移积分，不执行 SCF。

使用仓库乙烯冻结二聚体，B 相对于 A 沿 z 方向平移 3.5 Å。RKS/B3LYP、def2-SVP、def2-universal-jkfit 密度拟合、grid level 3、conv_tol 1e-10、4 线程。依赖为锁文件中的 PySCF 2.15.0 / NumPy 2.5.4。本报告是 **CPU** 实测，生成本报告时尚未执行新接口的 GPU 交叉验证。本机 NVIDIA 驱动不可用；后续已完成的验证见 [RTX 5090 两阶段报告](../transfer_control_gpu_2026-11-07/README.md)。之前直接计算接口的 GPU 验证独立记录于 `../transfer_gpu_2026-11-06/`。

## 结果

| 载流子 | DIPRO / eV | ESID / eV | orbital-GMH / eV |
|---|---:|---:|---:|
| 空穴 | 0.422000611513 | 0.421956884429 | 0.421956884429 |
| 电子 | 0.485180696540 | 0.479922306969 | 0.479922306969 |

ESID 仅在本例显式确认等价位点且数值/子空间检查通过时成立；orbital-GMH 是单电子近似。这些结果不是对任意材料方法有效性的证明。

保存 Fock 的 checkpoint 后处理与独立直接 SCF→投影参考的最大耦合差异为 **3.624878175401136e-14 eV**。仅保留 PySCF 标准字段的 dimer checkpoint，通过完整 MO 恢复特征算符后，与保存 Fock 结果最大差异 **2.911948873141057e-7 eV**；该小差异反映最终 SCF 轨道与重新构造 Fock 的有限收敛残差。两者在结果中明确标注不同算符来源，不把恢复恒等式当成收敛证明。

通过项目：

- 三个 XYZ 批量生成独立结果，控制文件 CLI 从不同工作目录读取相对路径；仓库附带示例控制也实际执行成功。
- 在禁止 `pyscf.scf.hf.SCF.kernel` 调用的条件下，后处理正常完成，`scf_executed=false`。
- dimer 原子交错排序及对应 AO 映射，结果与 A+B 排序一致；显式从 1 开始的 MO 选择与自动前线态选择一致。
- 原生标准 checkpoint 在明确确认 SCF 收敛后读取，结果标记 `user_assumed_from_external_log`；无确认时拒绝。缺失模型元数据不冒充已验证的一致方法。
- 拒绝未收敛标记、坏的 S 正交性、截断 MO、坐标不对应、不同电子方法和不一致 Fock；失败结果替换旧成功 JSON。
- ESID 缺少等价性假设时部分失败并退出 2；第一阶段 SCF 不收敛时退出 1 并移除旧 checkpoint。
- 旧的 XYZ 命令行 frozen-density 模式仍正常，空穴 0.478127106949 eV、电子 0.484594059282 eV，没有执行 dimer SCF。
- 34 项 Python 测试运行完成，其中唯一 GPU opt-in 测试跳过；C++、Fortran、Rust（4 项）测试及四语言对照通过；uv wheel/sdist 构建通过。

## 证据和复现

`summary.json` 为机器可读报告；`cpu/` 保存真实计算日志与状态；其他 JSON 为实际后处理结果。HDF5 checkpoint 和完整构建产物保留在工作目录 `build/transfer-control-validation/`，不提交二进制输入。原始 XYZ 在 `examples/quantum/`。`environment.json` 记录本机运行环境、源码与 wheel 的 SHA256。

```bash
uv run --locked --extra quantum python scripts/validate_transfer_control.py \
  --output build/transfer-control-validation
uv run --locked python scripts/check.py
```

最初完整执行了第一阶段和独立参考；随后使用同一目录的 `--reuse` 复用这些成功输出，验证最终后处理代码。报告中的 `reused_existing_converged_outputs` 表示这一点，未把读取已有结果误报为重新 SCF。`--reuse` 只适用于已生成该目录的情况。
