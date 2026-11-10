# 重组能 CSV 汇总

按用户给定列顺序，在每帧的 reorganization 输出根目录生成 hole.csv 和 electron.csv。
每个计算完成的分子占一行，mol 使用 manifest 中的分子 ID（XYZ 输入使用文件名标识）。
只为 carriers 请求的通道生成表格；多帧分别汇总，避免重复分子 ID 混淆。

表头固定为：

```text
mol,E0_R0,Eq_R0,Eq_Rq,E0_Rq,charged_relaxation,neutral_relaxation,total
```

四点能量保留 JSON 的 Hartree 单位，三个弛豫能字段使用 relaxation_ev，单位 eV。
CSV 不插入额外注释行或单位行；文档及 workflow.json 的 csv_units 记录单位。

汇总直接读取计算返回的完整结果，不改变量化实现或其缓存指纹。每完成一个分子，
以临时文件原子替换本次选择的汇总表；复用已完成计算时也重新生成表格。
失败或尚未完成的分子不写零值、不写成功数据。运行中可能存在只有表头或部分分子的表，
是否全部完成仍以 workflow.json 和各分子的 reorganization.json 状态为准。

验证覆盖精确列顺序、JSON 到 CSV 的单位映射、分子选择、空表、单通道、无效结果拒绝，
以及工作流的多帧、自定义路径和断点复用。
