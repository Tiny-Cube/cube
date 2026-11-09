# 顶层计算入口与分支开发

用户明确 charge-carrier、exciton、surface-hopping 是控制后续模块的入口，不是排在量化后的独立开关。三个段放 input 最前；do 严格 0/1，缺失为关闭，存在任意入口即启用路由。全关不执行模块，不验证不活动模块文件/物理设置；全无入口保留旧脚本兼容。

实现依赖标准库的 CalculationType 注册表，声明每类计算拥有的公共模块、solver 回调和未实现原因；路由/CLI/解析从注册表派生。当前 charge-carrier 运行现有四模块的参数准备，未接通迁移率则报告 prepared/退出2。exciton 和 surface-hopping 仅共享提取，完整求解预检明确失败，不能误用载流子量化。--only 限定入口/模块子集，不能绕过 do 或模块归属。共享模块只调度一次，隐式提取依赖同样检查归属。

公共基础分支 feat/calculation-entries，三类计算分支 feat/charge-carrier、feat/exciton、feat/surface-hopping 都从基础提交开始。使用独立 worktree 允许同时编辑；各计算分支仅开启对应入口，专用说明记录仍待实现的求解。共同基础不重写已有历史或破坏旧分支。

验证包括全关且数据文件不存在、入口/模块选择交集、别名与继承依赖、激子/SH 不走基态模块、新注册类型自动被解析与CLI识别、共享模块去重、回调求解与准备状态区分，以及现有真实小分子两帧SCF/DIPRO续算完整回归。
