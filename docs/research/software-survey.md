# 载流子迁移率软件调研

查阅：2026-11-05。cube 首批范围已由用户确定为**有机分子晶体、非晶、聚合物**。下表记录本次查证的能力，不是对所有功能的穷举；未安装完整软件不等于已完成数值复现。源码证据见 [审读记录](source-review.md)，版本和文件校验和见 [清单](sources.json)。

## 1. 与 cube 最直接相关的软件

| 软件 | 迁移率与主要方法 | 输入和适用问题 | 本次源码/许可核对 | 对 cube 的启发 |
|---|---|---|---|---|
| [MOMAP](https://www.hzw.ai/momap/index) | 有机分子电荷转移、随机游走、扩散与 Einstein 关系；厂商公开说明包含考虑核隧穿的全量子电荷转移速率 | 量化参数、振动信息、分子排列；还覆盖发光性质 | 查阅公开理论/产品说明；未取得核心源码与可核验源码许可证 | 电子参数→速率→输运的工作流；量子振动效应作为独立模型阶段 |
| [VOTCA-XTP](https://github.com/votca/votca) | Marcus 网络；bulk KMC、ToF、单粒子确定性网络求解等 | 分子形貌、位点能、耦合、内外重组能、极化环境 | C++；已阅读固定提交；Apache-2.0 | 分层数据模型、场与环境能、KMC/主方程相互验证 |
| [MorphCT](https://github.com/cmelab/morphct) | 发色团网络、Marcus 类速率、KMC、MSD 迁移率 | 分子/聚合物形貌、发色团原子映射、量化参数 | Python；已阅读；GPL-3.0 | 链段划分、链内/链间通道、形貌到电子网络的自动化 |
| [mcal](https://github.com/matsui-lab-yamagata/mcal) | 有机晶体的相关跳跃网络扩散/迁移率张量；解析路径连续性修正 | CIF、重组能/耦合计算；Gaussian、PySCF 或 ORCA 后端 | Python；已阅读并执行其解析张量函数；MIT | 张量而非单个标量；确定性解可检验 KMC，不能只对独立跳跃求平均 |
| [Excimontec](https://github.com/MikeHeiber/Excimontec) | 晶格 KMC；MA/Marcus 类跃迁；ToF、稳态输运，亦含激子与复合 | 有机薄膜/混合形貌、DOS、场和模型参数 | C++；已阅读应用层；MIT；KMC_Lattice 为另一个子模块 | 将材料模型、事件、边界、测量协议分开；器件功能分阶段引入 |
| [CharTED-KMC](https://github.com/thomasjamespope/CharTED-KMC) | Fortran KMC；Marcus/MA、电子/空穴/激子；通量与漂移迁移率 | 位点、相对取向、无序及相关模型参数 | Fortran；已阅读；GPL-3.0 | 四语言后端可行；显式平衡期、漂移通量、统计窗口 |
| [ToFeT](https://github.com/jfcoker/tofet) | 分子图 KMC，ToF/FET 与漂移/电流；源码含 Marcus 和 MA 分支 | 分子位置、耦合/能量、场、电极及可选相互作用 | C++；已阅读此研究分支；核心源文件头为 LGPL-3.0-or-later，仓库同时有 COPYING/COPYING.LESSER | 聚合物与器件研究接口；占据、场、库仑引起的速率更新需独立设计 |
| [LightForge](https://people.bath.ac.uk/pysabw/extmosmanual/extmos/LightForge/LightForgeDocu.html) | 有机材料/器件多粒子 KMC，迁移率、I–V、激子等 | 多层器件、显式或模型生成的材料/跃迁参数 | 查阅公开功能文档；[作者论文](https://doi.org/10.1038/s41467-023-36748-x)说明可申请研究试用；未取得源码 | 电极、层结构和相互作用属于器件级工作流 |
| [Bumblebee](https://www.scm.com/amsterdam-modeling-suite/bumblebee/) | 3D KMC，有[体迁移率](https://www.scm.com/doc.2026/Bumblebee/Examples/Bulk.html)、器件电流、激子与光电过程 | 材料模型/数据库、分层器件、场、接触与统计控制 | 查阅厂商功能和理论文档；商业授权产品，本次未取得源码 | bulk 和器件实验分开；重复形貌、稳态与误差随运行显示 |

上表有 6 个软件的实际源码审读，另外 3 个是公开文档研究；加上下一节 AMSET，共审读 7 个软件的源码。不是声称这些软件均已安装运行。ToFeT 按具体 fork 和源文件许可描述，不能仅凭 GitHub 自动标签判断许可。

## 2. 也能算迁移率，但物理路线不同

| 软件 | 已查证方法与用途 | 可访问性/本次研究 | 对 cube 的定位 |
|---|---|---|---|
| [AMSET](https://hackingmaterials.lbl.gov/amset/) | 从能带和材料参数估计声学形变势、离化杂质、极性光学声子、压电散射，经 BTE 得迁移率 | Python 公开源码；已阅读 transport.py；BSD 风格许可证，含附加贡献条款，以仓库 LICENSE 为准 | 学习散射与温度/浓度数据分层；其能带迁移率不能替代分子 hopping 迁移率 |
| [EPW](https://epw-code.org/) | QE/Wannier 电子–声子耦合、BTE；[官方参数](https://docs.epw-code.org/Inputs/Inputs.html)明确包含输运与电子/空穴迁移率 | 官方提供公开源码；本次查阅官方理论/参数，未审读或运行完整内核 | 未来能带输运接口的参照，单独子项目评估 |
| [Perturbo](https://perturbo-code.github.io/) | 第一性原理电子–声子相互作用、BTE 电导/迁移率与非平衡动力学 | 官方称开源；[下载流程](https://perturbo-code.github.io/mydoc_installation.html)要求申请访问；本次未申请或取得核心源码 | 学习电子结构与输运求解分离、HDF5、并行；不把公开后处理仓库误当核心源码 |

AMSET 的具体 Born 近似散射模型见其[官方说明](https://hackingmaterials.lbl.gov/amset/)。EPW/Perturbo 的主要路线不是分子位点 Marcus-KMC。此分类按输运机制，不能简化为“有机必然 hopping、无机必然能带”。高迁移率有机晶体同样可能需要离域、动态无序与非绝热理论。

仅有能带、电导除以未知弛豫时间、迁移积分或跃迁速率的软件，也不能自动称为能预测绝对迁移率的软件。例如量化化学程序和 tcal 可以供应参数，但还需要输运模型、边界条件与物理时间统计。

## 3. 对 cube 理论范围的判断

推荐首个可验收产品：**分子/链段跳跃网络的材料级载流子输运软件**，同时提供 KMC 和单载流子确定性求解。先计算可解释、可复现的张量与场温依赖，再扩展多载流子和器件。

这一选择是结合当前数据、代码与用户材料范围作出的设计建议，不是认为所有有机材料均可用 Marcus 描述。对于强链内耦合和部分离域聚合物，要允许更大的电子活性片段，并检查分片敏感性；不能将每个单体强制当独立跳跃点。

对于高纯晶体，需保留紧束缚 Hamiltonian、时间序列耦合和未来非绝热/离域求解器接口。已有 [FOB-SH 研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC8928466/)与[晶体研究](https://www.nature.com/articles/s41467-019-11775-9)说明这类边界；本次只研究其公开物理依据，没有声称审读或实现对应动力学程序。

## 4. 读到的说明与待核实项

- 用户最初给的 [pyMINDs theory](https://pyminds.readthedocs.io/en/latest/theory.html) 页面自称 MOMAP，但本次不能验证它与当前商业 MOMAP 版本的对应关系。其“选中通道时间取 1/k”和“同一种分子 ΔG=0”的说明不能直接作为 cube 的算法定义：正确静态 CTMC 等待时间由总逃逸率决定；同种分子也可有构象、静电与极化导致的位点能无序。
- MOMAP 官方页面明确提到全量子速率。本次只查阅公开资料，不将经典 Marcus 或单模 MLJ 宣称为完整复现其全量子模块。完整版本差异、模式约定与官方手册需在有合法可访问资料后逐项验证。
- 商业/需申请产品的源代码不可访问是审读范围限制。cube 的独立实现不依赖取得这些源码；将来可用用户依法生成的标准化中间数据比较数值。
- 没有取得相同结构、能量、耦合、速率、边界与采样协议之前，不用不同软件输出的一个迁移率数字判断程序正确性。
