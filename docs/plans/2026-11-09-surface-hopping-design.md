# Surface hopping 第一阶段设计

用户已选择首先面向有机材料的电荷／激子输运，并授权开始实现。实现位于 `feat/surface-hopping` 工作树；其他计算分支不改动。

## 方法选择

首先实现可检查的标准 FSSH，而后再接真实材料。NBRA 适合复用预计算 MD，但会忽略电子对核的反作用，不能用同一名字掩盖模型差别。MASH 是后续可比较的算法，需要独立实现及验证，不作为 FSSH 的局部开关。

现阶段采用固定、正交、实 diabatic 基，provider 提供 `H[n,n]` 和 `dH[n,n,dof]`，全部使用原子单位。由本征分解得到能量、Hellmann–Feynman 梯度与 NAC；通过相邻时间步同一 adiabatic 态的重叠符号保持相位连续。近简并、复杂 SOC、快速态翻转不采用静默近似，而是明确报错。运动中的非正交片段轨道不能直接满足此接口，需要额外处理基底导数与正交化。

核采用 velocity Verlet，电子使用分步指数中点传播。对每个核时间步积分电子态间净流量，并在步末进行至多一次随机 hop；这是需要时间步收敛的离散算法。成功 hop 沿 NAC 调整动量，选动量变化最小的能量守恒根；能量不足时保留原动量并记录 frustrated hop。去相干可选 `none` 或明确指定的 Granucci–Persico EDC。

## 入口与输出

保留三个顶层入口。`$surface-hopping` 的 `do` 门控新模块 `$surface-hopping-dynamics`，该模块引用 TOML；不复用载流子基态量化或带电重组能作为激发态动力学。完整 surface-hopping 请求缺少 dynamics 模块时，在 MD/QC 前拒绝。模型示例不含结构提取段，因而不依赖 GRO、TPR、PySCF 或 GPU。

输出独立轨迹、电子/活跃态占据、能量及 hop 诊断；通过控制/模型哈希和轨迹文件哈希复用完整结果。散射在 max-steps 仍未出边界时报告 `partial`，不会冒充收敛。第一阶段没有中断轨迹恢复、真实材料迁移率或 GPU 加速。

## 验证与后续

使用三个 Tully 模型检查解析导数，独立 diabatic TDSE 对照检查电子传播，非等质量多坐标检查动量守恒，三态 LVC 检查接口，并完成控制门控与复用测试。另用独立 FFT split-operator 量子波包和相同 Wigner 初始分布的 FSSH 轨迹对照一个散射工况，明确统计误差与 FSSH 近似误差。

后续以 FOB-SH/片段激子哈密顿量为材料方向，先解决随构型变化的位点能、带符号耦合、梯度及基底重叠，再接现有 MD/QM/MM 流程。GPU 量化接口需单独确认各态能量、梯度、NAC/overlap 的实际支持；CPU/GPU 轨迹批量传播在数值对照通过后再开展。
