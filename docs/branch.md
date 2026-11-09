# 载流子迁移率开发分支

分支 `feat/charge-carrier`。当前重点是主模块编排与可验证的独立单载流子 KMC。

已接通主段内逐项 0/1 选择、已有网络导入、提取＋耦合＋四点结果组网、Marcus 速率、静态零场 KMC、扩散／迁移率张量、轨迹统计与参考解；见 [使用说明](charge-carrier.md)。验证示例为 `examples/charge-carrier/input`，无需 MD 或量化。

完整材料 input 仍需要实际电子结构计算，并满足所选位点能／成对重组模型；不会把示例结果说成实际材料结果。GPU 材料 KMC、多载流子、有场漂移、动态无序和一般非等价势能面模型尚未接通。公共框架先进入 feat/calculation-entries，科学代码留在本分支。
