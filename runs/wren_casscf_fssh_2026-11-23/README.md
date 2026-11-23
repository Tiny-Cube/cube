# 基于 PySCF SA-CASSCF 的短轨迹验证

`validation.json` 保存 CPU 计算结果，`validate.py` 从 `examples/surface-hopping/` 读取受版本控制的 XYZ 和活性空间，重新运行 LiH 和扭曲乙烯 C₂H₄。命令在仓库根目录执行，输出到新的目录。

```bash
uv run --extra dynamics --extra quantum python runs/wren_casscf_fssh_2026-11-23/validate.py build/casscf-fssh-validation
uv run --extra dynamics --extra quantum wren examples/surface-hopping/ethene.input
```

- LiH 的两态解析梯度和不带 ETF 的带符号导数耦合与中心差分比较。检验脚本直接枚举含一个双占据 Li 核心轨道的行列式重叠，独立于传播器的 CI 基变换实现。
- 两个体系的核步长 0.5 和 0.25 au 在相同终点比较。传播使用 ETF，输出能量漂移、归一化误差、最小能隙与跨步多电子态重叠。短轨迹均无随机跳跃；另有受控非零概率向下跳跃的动量和能量回归测试。
- 乙烯的两个初始 RHF 活性轨道按局部 CH₂ 平面法向的碳 $2p$ 子空间投影核验。SA-CASSCF 的轨道优化引入混合。第一次更紧的乙烯轨道梯度阈值 $2\times10^{-7}$ 在 $3.1\times10^{-7}$ 左右停滞，PySCF 判为未收敛；记录中的乙烯输入显式使用 $10^{-6}$，两个步长的所有 SCF、CASSCF 和响应方程均通过收敛检查。
- 实际验证使用真空、小基组和单一起点。这些计算仅检验导数、态跟踪和积分器，不能据此估计晶体输运、光化学产率或寿命。
