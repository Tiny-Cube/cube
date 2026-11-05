# 参考来源与独立设计

当前 cube 材料级输运方案见 [2026-11-05 软件调研与源码证据](research/README.md)。下面保留 2026-11-02 通用 KMC 内核的历史设计依据。

查阅日期：2026-11-02。上游页面与 master 分支会变化。本项目独立编写，没有复制上游实现，也不宣称与其文件格式/完整功能兼容。

- [VOTCA 源码仓库](https://github.com/votca/votca)：参考速率计算与事件图/模拟层的分离。具体阅读了 [KMCCalculator](https://github.com/votca/votca/blob/master/xtp/src/libxtp/kmccalculator.cc) 与 [Rate_Engine](https://github.com/votca/votca/blob/master/xtp/src/libxtp/rate_engine.cc)。没有引入其拓扑、量子化学、Eigen/Boost 或数据库依赖。
- [VOTCA Charge and Energy Transport Model](https://www.votca.org/xtp/transport_theory.html)：参考连续时间指数等待与相对速率选择的 KMC 表述。其材料模型属于本项目可选适配器的背景，不是内核的类型约束。
- [pyMINDs theory 页面](https://pyminds.readthedocs.io/en/latest/theory.html)（页面标题显示 MOMAP）：参考 Marcus 跃迁、主方程、随机游走和扩散统计的联系。该页提及按选中通道取 1/k 的时间，我们没有把它用作通用 CTMC 时钟；这里使用总离开率的指数分布。

本项目 Marcus 定义明确取 `delta_g = G_final - G_initial`，公式为

`k = (2*pi/hbar) * |V|^2 / sqrt(4*pi*lambda*kB*T) * exp(-(delta_g+lambda)^2/(4*lambda*kB*T))`。

不要直接套用上游变量名或符号；不同上游代码采用的 deltaG 约定可能相反。电场、位移与载流子电荷的符号应在模型层自行构造 delta_g，并通过 detailed balance 或无场解析结果验证。

通用化后的设计决策：核心只要求 rates+action；状态类型完全由调用方定义；事件可动态生成；空间和粒子数观测可选。四语言独立实现，以统一随机流、文件样例、算法契约与解析统计进行交叉验证。VOTCA 与 pyMINDs 作为领域参考，不决定本库的应用范围。
