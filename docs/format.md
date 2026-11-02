# KMC1 内置模型交换格式

KMC1 是图/反应两个示例适配器的便携交换格式。任意用户模型通过库 API 接入，不必转换为图或反应文件。

空白分隔，支持 `#` 到行末的注释。整数只接受十进制非负数字，范围 0..9007199254740991；浮点接受普通十进制与 e/E 指数，不接受 NaN、inf、Fortran D 指数、十六进制、逗号、重复计数写法。截断、额外 token、非有限值和非法维度会失败。内存资源上限由主机/语言实现决定；此格式不是处理不可信大文件的沙箱。

## 有向图：单个自由跳跃实体

```text
KMC1 graph
N M INITIAL
SOURCE TARGET RATE DX DY DZ
... 共 M 条边
```

N>=1，M>=0，节点编号为 0..N-1。`RATE>=0`，零速率边禁用；当前节点的出边参与竞争。平行边、自环合法。无可用出边为吸收态。`DX DY DZ` 是每次事件累计的真实位移，而不是目标绝对坐标。周期系统应写跨界的最短物理跃迁位移以累计展开轨迹。若不需要空间信息，写 0 0 0。

## 质量作用反应

```text
KMC1 reaction
N_SPECIES M_REACTIONS
N0 N1 ...
C R0 R1 ... P0 P1 ...
... 共 M_REACTIONS 行
```

反应 j 的 propensity：`c_j * product_i choose(n_i, r_ji)`。R 是消耗数、P 是生成数（不是净变化），事件后 n'=n-R+P。人口不足时禁用。允许零阶生成、降解、催化剂、无变化通道。c 是**随机反应速率常数**，宏观速率、体积与阶乘因子由调用方转换。比如 `2A -> B` 对应 `a=c*A*(A-1)/2`。

## 命令行与 CSV

`ukmc MODEL [STEPS [SEED [UNTIL]]]`；原生命令名分别 `ukmc-cpp` / `ukmc-rust` / `ukmc-fortran`。

默认 STEPS=1000、SEED=1、UNTIL 无限制。若给出 UNTIL，应为有限非负十进制数。stdout 是 CSV；stderr 是 status/events/time 或错误说明。失败退出非零，可能已经写出一部分轨迹，调用方须检查退出码。

```text
step,time,event,x,y,z,state0,state1,...
```

初态 step=0、time=0、event=-1。事件记录包含事件执行后的状态与累计位移。事件 ID 对应输入中的第几条边/反应，从 0 起。截止时间和吸收态的末端观测使用 event=-1，不增加 step。图仅有 state0（当前节点）；反应有 N 个 state 列，xyz 恒为 0。Fortran 浮点列可能含前导空白，标准 CSV 数值读取时可正常解析。

输入数量、速率或种子非法时四版本均拒绝。参数错误 Fortran 使用 `error stop 2`；其他 CLI 返回 2。
