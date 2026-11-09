# 冻结 MM 的静电＋LJ12-6：实现验证

根据用户选择，为 QM/MM 四点法补入原力场的 QM–MM van der Waals 排斥和色散。此目录记录**小体系实现验证**；86-QM/5160-MM 真实体系的完整 CPU/GPU 计时独立运行，不能用这里的秒数代替材料计时。

## 实现和参数

`mm-vdw = tpr` 从 `gmx dump` 读取完整、已解析的 C6/C12 原子类型对矩阵，通过 `gmx dump -sys` 对应全局原子类型。没有推断混合规则或使用通用替代参数。对用户 TPR，实际读到 10 种类型、100 组 LJ12-6 参数；QM 为 86 原子、MM 为 5160 原子。类型映射沿用已验证的完整分子及周期镜像，中央 QM 分子不重复作为 MM。

`V = C12/r^12-C6/r^6` 的单位从 kJ/mol、nm 转为 Hartree 和 Hartree/Bohr；只计算 QM–MM 项。MM 不动，LJ 参数不随 QM 电荷态改变。所有优化能量/梯度和四个单点都包含同一个经典 LJ 模型，并分别记录电子结构能量、LJ 能量和总能量。标准 SCF checkpoint 的 `e_tot` 仍是电子结构部分，总能量记录在嵌入标记/结果中。

LJ 在固定 R0/Rq 的四点总和中代数抵消，但会改变 R0/Rq，因此必须进入优化。模型为完整固定分子有限簇的裸 LJ12-6，不含原 MD 的 switch/shift、PME 或尾项修正，也不含可极化 MM、跨边界成键或 MM 结构松弛。公式见 [GROMACS 非键相互作用](https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html)。

控制文件和导出方式见 [使用说明](../../docs/cube-input.md#加入原力场的-qmmm-范德华作用)。根目录 `input` 已启用 `mm-vdw = tpr`。`scripts/export_qmmm_cluster.py` 可以从原清单导出单个 QM XYZ、冻结点电荷和绑定的 LJ JSON，无需上传完整轨迹。

## CPU/GPU 实际验证

小体系使用 F2、4 个固定 MM 点电荷，LJ 参数为明确标注的合成数值测试参数，不是材料力场。设置 B3LYP/6-31G(d,p)、def2-universal-jkfit、grid 3、SCF 1e-10；RTX 5090 与同机 4 线程 CPU。

| 重组能 / eV | CPU | GPU |
|---|---:|---:|
| 空穴 | 0.409058976245 | 0.409060319149 |
| 电子 | 5.627757877685 | 5.627750811105 |

最大 CPU/GPU 差约 7.07e-6 eV，低于预设 5e-3 eV 门限。所有中性/带电优化、四点能量、冻结环境指纹、checkpoint 和无重算续算检查通过。对三个电荷态的**总能量**作中心差分，最大梯度差约 1.346e-5 Hartree/Bohr，通过预设 1e-4 门限。

独立单元测试核对：LJ 能量的 nm/Å 单位、解析导数、平移不变性、类型矩阵完整性/对称性、拒绝非 12 次幂和 A/B 变参；实际 HF 嵌入 SCF 能量差分同时核对优化器接收的完整梯度。LJ 参数文件绑定初始 QM 几何、原子顺序和 MM 坐标/电荷，不能混用其他输入。

完整 Python 测试 55 项：54 项通过、1 项显式 GPU opt-in 的 CuPy 代数测试跳过。GPU 量化全流程在上述 RTX 5090 另行实际执行。见 [单元测试](unit-tests.log)、[运行汇总](summary.json)、[完整预检日志](validation.log)、[CPU 原始四点能量](cpu/reorganization.json)、[GPU 原始四点能量](gpu/reorganization.json)。日志仅去除颜色码和行尾空白。

版本、测试 wheel 和源码哈希见 [environment.json](environment.json)；已核对测试 wheel 内全部 workflow 源文件与工作区一致。归档保留小分子输入与快照，真实材料坐标、完整 MM 环境和 checkpoint 不加入 Git。

```bash
uv run --locked --extra workflow --extra quantum-cuda13 python scripts/validate_qmmm_workflow.py \
  --output build/qmmm-lj-validation-new --gpu --lj
```

## 真实体系最初的失败诊断

同一个 86-QM/5160-MM 体系在纯静电模型下，中性优化约 10 次几何评价后触发 0.5 Å QM–MM 距离保护，进程用时 **689.386 s**、退出码 1，没有成功的重组能。该失败时间不能作为完整重组能计算时间，也不能用于报告 CPU/GPU 加速比。

用户明确选择补齐原力场 van der Waals 后重测。新基准从相同原始 MD 构型重新开始，保留原失败目录；GPU 使用 4 主机线程，CPU 使用 40 线程，两者按顺序独立运行，避免相互争用资源。完整材料结果另行归档。
