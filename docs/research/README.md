# cube 载流子输运调研与设计

日期：2026-11-05。用户确认首批材料为**有机分子晶体、非晶、聚合物**，产品命名为 **cube**。本次交付软件/源码调研、理论契约与实现规划，不代表已实现完整材料工作流。

| 用户需求 | 交付与核验证据 |
|---|---|
| 找到 MOMAP、VOTCA 之外的软件 | [软件比较](software-survey.md)：共 12 个，区分有机 hopping、器件与能带输运 |
| 学习理论 | [理论契约](transport-theory.md)：能量/场符号、Marcus/MA/MLJ、KMC、迁移率张量与漂移/ToF、确定性解、统计和适用域 |
| 学习公开源代码 | [审读记录](source-review.md)：7 个仓库的实际函数/接口审读；[固定提交和文件哈希](sources.json) |
| 给出独立数值证据 | [mcal 对照结果](mcal-crosscheck.json)，复现脚本为 `scripts/research_crosscheck.py`；只比较共同适用域的解析张量 |
| 规划 cube 实现 | [架构与里程碑](../plans/2026-11-05-cube-carrier-transport-design.md)：当前差距、三条路线、数据规范、四语言分工、M0–M6 与验收 |
| 现有功能与计划区分 | 项目 README 标为 cube；已有 universal-kmc/ukmc 可用，新 cube 工作流命令明确标为拟议 |

研究范围：7 个公开软件读取选定源码；mcal 的解析张量函数在已有真实体系上实际执行。其他程序的完整计算流程未安装或复现；MOMAP、LightForge、Bumblebee、Perturbo 核心源码本次未取得，EPW 仅研究官方文档。没有将文档功能说明当作完整源码验证或数值一致性结论。

源码检出目录为被 Git 忽略的 `build/research_sources/`。发行物保留研究证据、固定上游链接与 cube 自写脚本，不附带这些上游仓库。


## 后续构建讨论

用户已补充要求从 GROMACS/LAMMPS 单帧出发，按所选共轭区域提取 dimer、计算转移积分与 QM/MM 重组能，并同时支持量化和 KMC 的 GPU 路线。见 [构建讨论稿](../plans/2026-11-05-cube-construction-discussion.md)、[更新后的 MD/量化/GPU 工作流](../plans/2026-11-05-cube-md-quantum-gpu-workflow.md)。当前已接通 GROMACS 提取、CPU/GPU DIPRO 和固定 MM 坐标/电荷的四点法，后者只优化 QM 分子；已增加原 TPR 的 QM–MM LJ12-6 排斥/色散。原静电模型 GPU 验证见[报告](../../runs/qmmm_frozen_mm_2026-11-07/README.md)，新参数接口见[控制说明](../cube-input.md#加入原力场的-qmmm-范德华作用)。跨边界成键、可极化环境、LAMMPS 自动导入和 GPU KMC 仍属于后续目标。

[共轭区域距离示意与数值](dimer-distance-examples/README.md)为讨论用合成例子；[GPU4PySCF 补充源码证据](gpu-qmmm-sources.json)是早期静态调研，保留当时状态。后续实际 GPU 运行以对应测试报告为准。
