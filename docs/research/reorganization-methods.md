# 通用重组能模块：方法与软件借鉴

调研日期：2026-11-11。范围为官方文档、公开项目接口及原始研究；未安装这些软件进行
数值对照，也没有复制其源码。此文描述可借鉴的科学方法和数据契约，不代表全部已实现。

## 方法及适用范围

| 路线 | 必要输入 | 主要输出与边界 | 借鉴来源 |
|---|---|---|---|
| 优化势能面四点法 | 两个明确定义的电子态、各自优化结构、四个一致模型下的总能量 | 分别输出两个势能面上的形变能；必须明确单体形变与完整转移反应的区别 | 现有 cube；[原始研究的四点/振动模比较](https://pubs.rsc.org/en/content/articlehtml/2024/tc/d4tc01744g) |
| 正常振动模/绝热 Hessian | 两态结构、质量、Hessian、允许移动的原子 | 模分辨形变能、Huang–Rhys 因子、Duschinsky 关系；谐近似不保证等于大幅非谐四点能差 | [ORCA ESD](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/esd.html)、[DUSHIN 项目](https://github.com/lenhanpham/dushin)、[MOMAP](https://www.hzw.ai/momap/index) |
| 垂直梯度近似 | 一个参考几何上的另一电子态梯度与指定参考 Hessian | 不做完整目标态优化的谐势能面估计；必须标记 VG 等近似，不冒充优化四点结果 | [ORCA ESD 的 VG/VH/AH/AHAS](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/esd.html) |
| 能隙采样/线性响应 | 明确 A/B 态和温度的平衡轨迹、同一构型上的垂直总能隙 | 自由能重组量及统计误差；需要检查线性响应、能隙分布与采样，不能拿单帧代替 | [能隙方差与四点定义比较的原始研究](https://doi.org/10.1021/acs.jctc.5c01753) |

ORCA ESD 将 PES 建模、态计算、梯度和 Hessian 组织为不同层次，说明计算近似应显式选择。
DUSHIN 的公开版本以外部 Gaussian log/fchk 为输入进行振动模重组能分解，说明后处理可以
独立于产生量化数据的程序。MOMAP 是相关振动耦合/光物理路线的参考，不能仅据产品页面
宣称其所有数据接口已被验证或可直接复用。

## 统一定义，避免因子二混淆

令 R_A、R_B 是两个态的参考平衡几何，定义：

```text
lambda_A = E_A(R_B) - E_A(R_A)
lambda_B = E_B(R_A) - E_B(R_B)
sum_relaxation = lambda_A + lambda_B
```

这两项不必相等。若 A/B 是完整反应的两条等曲率 Marcus 势能面，则通常 lambda_A=lambda_B=lambda，
两项之和是 2lambda；不能自动把 sum 当作任意反应的 Marcus lambda。
若 A/B 是单个分子的中性/带电态，等价分子间载流子跳跃的分子内重组参数常使用两项之和。
非等价供受体则需要组合具体供体和受体的形变代价，不能不加核对地平均单体 total。
S0/S1 单分子的 lambda_S0、lambda_S1 也应分别保留，光谱、局域激子跳跃的消费方式需显式说明。

TDDFT/TDA 的激发能不等于该态的绝对总能量。四点输入应使用同一约定下的态总能量。
参见 [PySCF 激发态接口](https://pyscf.org/user/tddft.html) 与
[激发态优化及态交换问题](https://pyscf.org/user/geomopt.html#excited-states)。

## 环境与硬件独立

固定 MM、环境弛豫、永久多极和可极化环境是不同的物理协议，而非可任意交换的文件格式。
固定 MM 模式四个点使用完全相同的坐标、电荷和静电位点。允许响应的环境必须另行规定
各能量点的平衡/非平衡响应和能量记账，不能据 QM 优化四点宣称得到了环境外球重组能。

PySCF 的标准 qmmm 接口提供外部点电荷嵌入；LJ、高阶多极或极化模型需有对应的能量和
梯度实现。[PySCF 官方说明](https://pyscf.org/user/qmmm.html)
CPU/GPU 能力要按具体态、能量、梯度、Hessian 和环境组合检查，不能由 GPU SCF 可用
推断 GPU 激发态 QM/MM 梯度也已可用。
