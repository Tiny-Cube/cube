# G16 / GPU4PySCF：6 Å 冻结 MM 对照准备与数值验证

2026-11-08 提交此记录时，真实材料正式计时处于队列中，尚无 G16/GPU 速度倍数。已有 12 Å CPU 四点法仍在运行；独立后台队列会等待其退出及资源空闲，再依次执行两个引擎的初始能量/梯度与完整四点法。不停止、不修改已有计算。

材料选择 anneal_a 第 1980 帧第 85 个分子：86 QM 原子、2580 冻结 MM 原子（6 Å，30 个完整分子镜像）。三个必要输入为导出的 QM XYZ、MM 点电荷与原 TPR 的 LJ 文件，分别记录 SHA256；真实几何只保存在忽略的 build 和已授权计算服务器上，不提交。

## 模型与比较范围

用户选择固定点电荷＋原 TPR 的 QM–MM LJ12-6，两端使用同一个 LJ 实现、同一冻结环境和 geomeTRIC，只优化 QM。G16 通过 `Force Charge NoSymm` 及 formatted checkpoint 提供电子能量与 QM 梯度，坐标及元素顺序逐次检查。这是 **G16 外部电荷 QM＋cube LJ**，不是 G16 原生 ONIOM/Berny 的计时。

正式设置：B3LYP（Gaussian/VWN-RPA 版本）、6-31G(d,p)、球谐 5D/7F、SCF 目标 1e-10、名义 99×590 DFT 网格。G16 使用 UltraFine 与常规积分；GPU4PySCF 使用相同名义网格、原生剪枝和 DF/def2-universal-jkfit；网格实现和积分算法不同，需要一并报告数值误差。G16 用 40 个 CPU 线程，GPU 用 RTX 5090＋4 个主机线程；这代表整套生产配置的时间差，不单独归因于 GPU 硬件。

G16 输出含 MM 点电荷自能常数，必须扣除才能和 PySCF 约定直接比较。H2/STO-3G＋两个固定电荷的实际探针验证：扣除后能量差约 1.32e-10 Ha，梯度最大差约 1.31e-9 Ha/Bohr。梯度直接读取 Cartesian Gradient，不与 Gaussian 打印的 Forces 混淆，也不人为移除冻结外场对 QM 的净力。

## 已完成的预检

使用公开 F2 小体系、四个固定点电荷与相同合成 LJ，实际运行 G16 C.01 和 GPU4PySCF 1.9.1：

- 三个电荷态 q=0、+1、−1 的最大梯度差不超过 2.02e-5 Ha/Bohr。
- 完整四点法均收敛；空穴重组能 G16−GPU 约 +0.0001801 eV，电子约 −0.0003318 eV，均低于预设 0.005 eV 容差。
- G16 三个电荷态的原子顺序和坐标框架检查均通过。
- 小体系运行时旧 CPU 任务仍在运行，所有预检耗时均排除出性能结论，不能用它们推算材料加速比。
- 新增两项接口回归通过：点电荷自能扣除/梯度符号及冻结外场净力保留；失败计算和不完整梯度不被接受。

精确数值见 [summary.json](summary.json)。源实现见 [benchmark_qmmm_engines.py](../../scripts/benchmark_qmmm_engines.py)。真实任务后台目录名 `validation/g16-gpu-qmmm-6A-20261108`，状态文件 `comparison.json`；进程独立于 SSH 会话。

## 独立使用

该脚本是基准工具，尚未把 `g16` 添加为 cube input 的 `qc-engine`。先以已导出的同一个簇测三种电荷态初始能量＋梯度，输出目录必须是新的：

```bash
uv run --locked --extra workflow python scripts/benchmark_qmmm_engines.py \
  --engine g16 --task force --threads 40 \
  --xyz cluster/molecule.xyz --mm-charges cluster/environment.pc --lj cluster/lj.json \
  --g16 /path/to/g16 --output results/g16-force
```

GPU 改为 `--engine gpu --threads 4`，运行环境需 GPU4PySCF/CUDA。`--task four-point` 对同一起点独立完成中性/阳离子/阴离子优化及七个唯一单点。三个初始梯度通过数值对照后，后台队列才启动材料完整四点法；不以单步耗时推算优化总时间，不复用另一引擎的优化坐标。
