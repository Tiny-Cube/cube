# GPU 单载流子 KMC

在 `$kmc` 中设置 `backend = gpu`，Wren 会在 NVIDIA GPU 上同时推进多条独立的单载流子轨迹。
每条轨迹使用同一静态分子网络中的 Marcus 跃迁速率，分别统计空穴或电子在零场下的位移。
对多帧 MD 输入，每个结构快照独立计算。位点占据排斥、载流子间作用和随轨迹时间变化的速率
需要相应的多体或时变模型。

## 运行示例

```bash
uv sync --locked --extra transport-cuda13
uv run --no-sync wren examples/charge-carrier/input-gpu --plan
uv run --no-sync wren examples/charge-carrier/input-gpu
```

示例使用每个晶胞含一个位点的周期链，两个方向的跃迁都保留了晶格位移。`$kmc` 也可读取已有
的分子网络，或接收同一控制文件中结构提取、转移积分和重组能模块准备的数据。CPU 运行使用
`backend = python`，速率、初始分布、观测时刻和分析算法与 GPU 配置相同。

```text
$kmc
do = 1
input = network.json
output = results/gpu-kmc
backend = gpu
carriers = hole,electron
temperature = 300
duration = 1e-9
trajectories = 512
samples = 51
seed = 20261109
trace-trajectories = 3
max-events = 1000000
$end
```

`duration` 单位为秒，`samples` 包含起点和终点。任一轨迹在终止时刻前达到
`max-events` 时，其后观测值填充 NaN，载流子估计状态为 `partial`。报告保留
受影响的轨迹编号；提高步数上限后可在新的输出目录运行。

## 跃迁时间、位移和随机数

对当前位点 $i$，令 $K_i=\sum_j k_{ij}$。从开区间 $(0,1)$ 取均匀样本 $u$，停留时间为
$-\ln u/K_i$。若下一次跃迁位于观测窗口内，以第二个独立样本从累积速率中选择出边。
等待期间的坐标保持上一次跃迁后的位置；跨周期边直接累加有符号的位移。吸收位点保持原位直到
终止时刻。[Gillespie 直接法原文](https://doi.org/10.1021/j100540a008)

CPU 负责建立速率图、根据 Boltzmann 权重抽取起始位点、计算统计量并写出结果。设备端以一个
CUDA 线程推进一条轨迹，将有共同起点的通道连续存储。轨迹按批发送；同一轨迹被分配到哪个批次
不影响其随机流。GPU 直接导入 Python `random.Random` 的 MT19937 初始状态，在设备端生成相同的
53 位均匀数。轨迹种子按用户种子、载流子类型、轨迹编号和随机流用途派生。跳跃时间使用双精度，
平台的对数函数在浮点边界附近可能产生舍入差异。
[Matsumoto 与 Nishimura 的 MT19937 论文](https://doi.org/10.1145/272991.272995)

## 结果与核验

| 文件 | 内容 |
|---|---|
| `result.json` | CUDA 设备、驱动、CuPy 版本、源码哈希、逐轨迹事件数、精确网络参考以及张量估计和统计误差。 |
| `hole_samples.npz`、`electron_samples.npz` | 共用物理时间网格 `time_s` 与轨迹 × 时间 × 三方向的无包裹位移 `displacement_m`。 |
| `*_traces.json` | 前 `trace-trajectories` 条轨迹的起始位点、事件种子、逐事件时间、通道、位点及累计位移。 |
| `*_observations.csv` | 观测时刻的平均位移和协方差；扩散张量取协方差对时间的斜率的一半。 |

在有 NVIDIA GPU 的主机上运行下列测试。它会对相同种子的 CPU/GPU 轨迹比较初态、事件顺序、
事件计数、时间、位移和扩散张量，再通过正常 `$kmc` 入口检查持久化结果与续用。

```bash
WREN_TEST_GPU_KMC=1 uv run --no-sync python -m unittest discover -s tests -p 'test_transport_gpu.py' -v
```

示例周期链的解析结果为 $D_{xx}=ka^2$，其中 $k$ 为单方向跃迁率，$a$ 为周期长度。
对实际材料应另外检查轨迹数、模拟时长、拟合窗口及快照采样。计时比较时固定输入、轨迹数、
时长、事件上限和追踪数量，记录源码版本、硬件及首次内核编译开销。短轨迹的小批量计算中，
CUDA 初始化和主机设备数据传输可能占主要时间。[CuPy RawKernel 文档](https://docs.cupy.dev/en/stable/reference/generated/cupy.RawKernel.html)
