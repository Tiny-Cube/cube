# 表面跳跃动力学开发分支

分支：`feat/surface-hopping`。公共基础提交：`813959bb2cb171d140d3308e8db72d23c6234ba3`。

本分支的 input 只开启 `$surface-hopping`，其他计算入口关闭。已实现 Python/NumPy 模型 FSSH：Tully 1/2/3、固定正交基多态 LVC、控制文件门控、轨迹输出、相位连续、NAC 动量调整与可选 EDC。运行示例使用 `examples/surface-hopping/input`，详见 [使用说明](surface-hopping.md) 与 [理论／开源调研](research/surface-hopping.md)。根 input 保留材料准备配置，但缺少 dynamics 模块时完整入口会明确报错；不会启动材料量化。

当前不是完整材料输运产品：真实位点能／耦合导数、移动片段基处理、量化 provider、GPU 与迁移率尚未接通。本次的 registry 扩展暂留本分支，待接口稳定后抽取公共部分到 `feat/calculation-entries`；详见 [分支协作说明](development-branches.md)。
