# 上游源码审读与可复现证据

2026-11-05。已获取并审读 7 个公开仓库的选定实现；这里的“审读”指实际读取下列核心函数和上下游接口，不是全代码形式化审计。源码保存在本机被 Git 忽略的 `build/research_sources/`，没有复制进 cube 源码或发行物。

精确提交、每个审读文件的 SHA-256、文件链接、用途和执行范围在 [sources.json](sources.json)。下文的固定版本索引由该清单生成。不同版本的目录或算法可能变化，不以 master 链接作为永久证据。

## 1. VOTCA-XTP

实际阅读 `rate_engine.cc`、`kmcwaitingtime.cc`、`kmcengine.cc`、`bulksimulation.cc` 的相关段落。

- `Rate_Engine` 把场、位点能和正反重组能接到速率；本次版本的 `Marcusrate` 使用 `(deltaG−reorg)^2` 形式，因此 cube 的 final-minus-initial 约定不能按变量名直接代入。
- 静态等待时间预计算位点出边累计率；引擎用候选事件队列与占据状态驱动运行。另有时变率 thinning/积分求解接口，说明时间依赖应是显式算法选择。
- `BulkSimulation::Exact` 检查闭合类，解稳态与 corrector、报告稀疏求解残差；驻留时间与观测窗口分开处理。

cube 的设计动作：统一自由能符号、CSR 出边与模型/时钟接口、以一般确定性求解交叉检验 KMC。没有编译或运行 VOTCA 全套，不宣称已完成与其数值回归。

## 2. MorphCT

阅读 `mobility_kmc.py`、`helper_functions.py`、`kmc_analyze.py` 与工作流 README。

其单载流子路径对各邻居生成指数等待时间并取最早者；在静态速率条件下与总率 direct SSA 在分布上等价，但随机数消耗不同，不要求逐轨迹同 seed 相同。代码显式携带邻居镜像，并区分链内/链间平均率选项；分析以三维 MSD 斜率/6 与 Einstein 关系生成标量迁移率。

cube 的设计动作：保留聚合物片段/链映射和不同通道，分析升级为完整张量与用户指定输运维数。模型中的各种近似应由配置明示。只做源码研究，未运行 MorphCT 流水线。

## 3. mcal

阅读 `hopping_mobility_model.py` 和 `rcal.py` 的入口、四点任务/能量组合。解析扩散函数构造对称跳跃生成元及伪逆修正项，描述具有周期连通路径的相关跳跃；与不考虑路径关联的独立跳跃近似不同。其标准化流程自动加入对称反向边，不是任意非对称/有场速率的通用求解器。

实际执行了该版本的 `diffusion_coefficient_tensor`，使用 cube 已发布的 `e_006` 几何、周期 image 和同一组速率。mcal 使用晶胞平移作为跳矢量，cube 使用含胞内位移的跳矢量；长期扩散张量应对这种有界胞内坐标 gauge 不变，本次对照也检验了该性质。

结果见 [mcal-crosscheck.json](mcal-crosscheck.json)：空穴/电子张量的相对 Frobenius 误差均小于 10⁻⁹（实际约 10⁻¹³/10⁻¹⁴ 量级）。这只验证**静态、对称、同一网络/同一速率**下的解析 D，未验证量化参数计算、完整 CLI、有限场或能量无序。

复现：

```bash
# 需要系统 Python 已装 NumPy，或用 uv run --with numpy 执行最后一行。
git clone https://github.com/matsui-lab-yamagata/mcal.git build/research_sources/mcal
# 按 sources.json 的 mcal commit 切换版本，再执行：
python3 scripts/research_crosscheck.py \
  --mcal-source build/research_sources/mcal \
  --output build/research_sources/mcal-crosscheck-reproduced.json
```

若目录已存在直接复用；脚本检查提交与目标函数哈希，不会自己下载或改动上游文件。NumPy 浮点环境不同可造成末位变化。输入哈希也已记录；Git 的 LF 规范化可能改变 CSV 字节哈希，而不改变解析数值。

## 4. Excimontec

阅读 `Polaron.h` 的 Hop 速率及 `OSC_Sim.cpp` 的事件应用与 `calculateMobilityData`。MA 与 Marcus 速率均通过应用层事件设置；ToF 迁移率使用器件长度、势差和 transit time，不能等同于零场扩散估计。README 和 `.gitmodules` 确认底层依赖独立的 KMC_Lattice 仓库；本次没有取得该子模块并运行完整引擎测试。

cube 的设计动作：边界条件与观测协议独立于事件内核；bulk、ToF、器件电流拥有不同的配置与结果标签。读取 MIT 文件许可，不自动推定其所有依赖都同许可。

## 5. CharTED-KMC

阅读主程序、`modules/event.f` 的 `myrate` 与 `modules/printings.f` 的通量/窗口分析。主程序使用 `−log(u)/rate_tot` 推进时钟，事件层含 Marcus/MA 分支，输出层计算通量与迁移率并检查后期窗口。

cube 的设计动作：四语言使用统一物理时间/通量定义，生产运行要有平衡和窗口分析。统计置信区间按 cube 独立验证的轨迹/分块方法实现，不因上游存在窗口分析就视为误差估计已正确复现。未编译/运行此软件。

## 6. ToFeT

阅读研究分支 `jfcoker/tofet` 的 `kmc.cc`、`graph.cc`、`vertex.cc` 和接口头。图层将分子对转换为正反边，静态时可一次计算速率，库仑变化时更新能量相关部分。顶点层处理空邻居选择与驻留时间。`FRM`/`FRM_FET` 表达不同运行边界和漂移/电流观测。

cube 的设计动作：分开不可变几何/耦合、依场温生成的静态率及随占据变化的率；新多粒子算法先与全量重算比较。该 fork 的主文档描述 Marcus，但实际源码也有 MA 分支；能力表据源码注明。源文件头为 LGPL-3.0-or-later，不能仅依 GitHub 根目录 GPL 自动标签判断。

## 7. AMSET

阅读 `core/transport.py`。实现从分机制散射率求寿命、通过能带速度/态密度积分得到电导，再按电子或空穴浓度换算迁移率；不是实空间跳跃的等待时间算法。

cube 的设计动作：同一个“mobility”字段必须附理论、载流子、温度、浓度与协议；未来 BTE 适配器使用不同的电子结构/散射参数规范。未执行 AMSET 的完整第一性原理流程。许可证是带额外贡献条款的 BSD 风格文本，以固定提交 LICENSE 为证据。

## 8. 文档研究的边界与转化规则

MOMAP、LightForge、Bumblebee 的核心源码未取得；Perturbo 核心源码需申请访问，本次没有提交申请；EPW 只阅读官方文档。没有声称对这些程序做源码审计。

公开理论可用来建立独立模型，但函数名、参数名、时钟算法、默认单位和输出“迁移率”的定义须分别核对。cube 本次只增加自己的研究文档与调用已审读 mcal 函数的对照脚本，没有引入上游函数实现或商业二进制。

## 固定版本索引

### votca `e2852159efc2c616ce3c99497caa05ec1e9f34ae`

- [xtp/src/libxtp/rate_engine.cc](https://github.com/votca/votca/blob/e2852159efc2c616ce3c99497caa05ec1e9f34ae/xtp/src/libxtp/rate_engine.cc) — Rate_Engine::Rate/Marcusrate; pair energy and field conventions
- [xtp/src/libxtp/kmcwaitingtime.cc](https://github.com/votca/votca/blob/e2852159efc2c616ce3c99497caa05ec1e9f34ae/xtp/src/libxtp/kmcwaitingtime.cc) — ExponentialWaitingTime and static-rate contract; thinning interface
- [xtp/src/libxtp/kmcengine.cc](https://github.com/votca/votca/blob/e2852159efc2c616ce3c99497caa05ec1e9f34ae/xtp/src/libxtp/kmcengine.cc) — carrier pending-event queue, occupancy and stopping limits
- [xtp/src/libxtp/bulksimulation.cc](https://github.com/votca/votca/blob/e2852159efc2c616ce3c99497caa05ec1e9f34ae/xtp/src/libxtp/bulksimulation.cc) — occupation-time sampling and Exact stationary/corrector solve
- [LICENSE](https://github.com/votca/votca/blob/e2852159efc2c616ce3c99497caa05ec1e9f34ae/LICENSE) — Apache-2.0

### morphct `cae0166d73994506f90645d33ede4cbe42a7be2d`

- [morphct/mobility_kmc.py](https://github.com/cmelab/morphct/blob/cae0166d73994506f90645d33ede4cbe42a7be2d/morphct/mobility_kmc.py) — Carrier hop: per-neighbor exponential competition, image, intra/inter choices
- [morphct/helper_functions.py](https://github.com/cmelab/morphct/blob/cae0166d73994506f90645d33ede4cbe42a7be2d/morphct/helper_functions.py) — get_hop_rate/get_event_tau
- [morphct/kmc_analyze.py](https://github.com/cmelab/morphct/blob/cae0166d73994506f90645d33ede4cbe42a7be2d/morphct/kmc_analyze.py) — calc_mobility/get_times_msds
- [README.md](https://github.com/cmelab/morphct/blob/cae0166d73994506f90645d33ede4cbe42a7be2d/README.md) — chromophore workflow, scope and license
- [LICENSE](https://github.com/cmelab/morphct/blob/cae0166d73994506f90645d33ede4cbe42a7be2d/LICENSE) — GPL-3.0

### excimontec `19a568f3a5871efad0db54f5bb950bb9e2bc6b65`

- [src/Polaron.h](https://github.com/MikeHeiber/Excimontec/blob/19a568f3a5871efad0db54f5bb950bb9e2bc6b65/src/Polaron.h) — Hop::calculateRateConstant MA/Marcus parameter conventions
- [src/OSC_Sim.cpp](https://github.com/MikeHeiber/Excimontec/blob/19a568f3a5871efad0db54f5bb950bb9e2bc6b65/src/OSC_Sim.cpp) — calculateMobilityData, event rate setup and ToF/steady-state application
- [README.md](https://github.com/MikeHeiber/Excimontec/blob/19a568f3a5871efad0db54f5bb950bb9e2bc6b65/README.md) — capabilities and KMC_Lattice dependency
- [.gitmodules](https://github.com/MikeHeiber/Excimontec/blob/19a568f3a5871efad0db54f5bb950bb9e2bc6b65/.gitmodules) — upstream engine submodule identity
- [LICENSE](https://github.com/MikeHeiber/Excimontec/blob/19a568f3a5871efad0db54f5bb950bb9e2bc6b65/LICENSE) — MIT

### charted-kmc `8624fbde8cd56b85c37f37c9894f0a8a6267594a`

- [code/CharTED-KMC.f90](https://github.com/thomasjamespope/CharTED-KMC/blob/8624fbde8cd56b85c37f37c9894f0a8a6267594a/code/CharTED-KMC.f90) — main event loop and total-rate clock
- [code/modules/event.f](https://github.com/thomasjamespope/CharTED-KMC/blob/8624fbde8cd56b85c37f37c9894f0a8a6267594a/code/modules/event.f) — myrate Marcus/MA branching
- [code/modules/printings.f](https://github.com/thomasjamespope/CharTED-KMC/blob/8624fbde8cd56b85c37f37c9894f0a8a6267594a/code/modules/printings.f) — flux mobility and recent-window statistics
- [README.md](https://github.com/thomasjamespope/CharTED-KMC/blob/8624fbde8cd56b85c37f37c9894f0a8a6267594a/README.md) — transport capabilities
- [LICENSE.md](https://github.com/thomasjamespope/CharTED-KMC/blob/8624fbde8cd56b85c37f37c9894f0a8a6267594a/LICENSE.md) — GPL-3.0

### mcal `282d0145c0509a967c4a8d8b6505d08cfcb9050b`

- [src/mcal/calculations/hopping_mobility_model.py](https://github.com/matsui-lab-yamagata/mcal/blob/282d0145c0509a967c4a8d8b6505d08cfcb9050b/src/mcal/calculations/hopping_mobility_model.py) — marcus_rate, mobility_tensor, analytic diffusion, standardization and MC reference
- [src/mcal/calculations/rcal.py](https://github.com/matsui-lab-yamagata/mcal/blob/282d0145c0509a967c4a8d8b6505d08cfcb9050b/src/mcal/calculations/rcal.py) — four-point calculation driver and energy combination
- [README.md](https://github.com/matsui-lab-yamagata/mcal/blob/282d0145c0509a967c4a8d8b6505d08cfcb9050b/README.md) — input and quantum backends
- [LICENSE](https://github.com/matsui-lab-yamagata/mcal/blob/282d0145c0509a967c4a8d8b6505d08cfcb9050b/LICENSE) — MIT

### amset `562904077533fba1de5d41a81930a84602590e67`

- [amset/core/transport.py](https://github.com/hackingmaterials/amset/blob/562904077533fba1de5d41a81930a84602590e67/amset/core/transport.py) — solve_boltzman_transport_equation, _calculate_mobility and scattering lifetimes
- [README.md](https://github.com/hackingmaterials/amset/blob/562904077533fba1de5d41a81930a84602590e67/README.md) — software identity and theory paper
- [LICENSE](https://github.com/hackingmaterials/amset/blob/562904077533fba1de5d41a81930a84602590e67/LICENSE) — BSD-style redistribution plus enhancement clause

### tofet `332a7d818e1a8722c41d1e06b87c573628d506d7`

- [source/kmc.cc](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/source/kmc.cc) — FRM/FRM_FET, drift mobility and time/output protocol
- [source/graph.cc](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/source/graph.cc) — paired edges, field and static vs interacting rate updates
- [source/vertex.cc](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/source/vertex.cc) — Marcus/MA rates, occupied destinations and residence-time weighting
- [source/kmc.h](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/source/kmc.h) — transport state/observables
- [README.md](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/README.md) — ToF/FET scope
- [COPYING](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/COPYING) — GPL-3.0 text shipped with LGPL
- [COPYING.LESSER](https://github.com/jfcoker/tofet/blob/332a7d818e1a8722c41d1e06b87c573628d506d7/COPYING.LESSER) — LGPL-3.0 terms; source headers say version 3 or later
