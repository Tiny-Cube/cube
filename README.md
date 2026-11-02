# Universal KMC — 通用动力学蒙特卡洛方法库

四套**独立实现**：Python、C++17、Rust、Fortran 2008。Python 用 **uv** 构建。

核心面向任意连续时间 Markov 跳跃过程：你定义状态、当前事件的速率和事件执行后的变化，引擎负责选择事件与推进物理时间。状态可以是对象、字符串队列、自旋数组、图、设备状态、包含实数属性的结构体等。**不要求粒子、坐标、化学计量或输运模型。**

## 从自己的问题开始

只需要两个方法：

```python
from universal_kmc import Engine

class Machine:
    def __init__(self):
        self.online = True
        self.failures = 0

    def propensities(self):
        # 故障 / 修复：当前不能发生的事件速率为 0。
        return [0.02, 0.0] if self.online else [0.0, 0.5]

    def fire(self, event):
        self.online = (event == 1)
        self.failures += (event == 0)

model = Machine()
engine = Engine(model, seed=42)
for record in engine.run(max_steps=10000, until=1000.0):
    print(record.time, model.online, model.failures)
```

也可以使用 `EventModel`，组合事件的 `rate`、`action` 和可选的 `enabled` 回调；事件工厂可以在每步生成不同数量的事件。示例见 [Python](examples/generic.py)、[C++](examples/generic.cpp)、[Rust](rust/examples/generic.rs)、[Fortran](examples/generic.f90)。四种语言的排队示例已做逐结果一致性检查；Python 还提供 Ising 自旋翻转示例。

## 安装、构建、验证

需要 Python ≥ 3.11、uv；原生版本需要 C++17 编译器、Cargo/Rust、gfortran 与 ar。没有第三方运行时依赖。

```bash
uv sync --locked
uv run python examples/generic.py
uv build                            # wheel + sdist

uv run python scripts/build.py      # 三种原生版本、静态库和示例
uv run python scripts/check.py      # 完整构建、全部测试、跨语言检查、打包
```

脚本默认使用 `c++` / `gfortran`，可通过 `CXX` / `FC` 指定兼容编译器。脚本的 Fortran 检查参数面向 gfortran；其他编译器可使用 CMake 或自行配置参数。`scripts/check.py` 默认把 uv 缓存放在系统临时目录；若直接运行 uv 遇到缓存目录无写权限，可设置 `UV_CACHE_DIR=/tmp/ukmc-uv-cache`。

可选 CMake 构建 C++/Fortran（本机未实测；已验证的构建入口为上面的 scripts/build.py）：

```bash
cmake -S . -B build/cmake
cmake --build build/cmake
ctest --test-dir build/cmake --output-on-failure
```

运行通用回调示例：

```bash
./build/generic-cpp
cargo run --manifest-path rust/Cargo.toml --release --example generic --offline
./build/generic-fortran
```

## 核心契约

| 能力 | Python | C++ | Rust | Fortran |
|---|---|---|---|---|
| 自定义模型 | `propensities/fire` 两方法 | 派生 `Model` | 实现 `Model` trait | 扩展 `model_t` |
| 回调适配器 | `EventModel[S]` | `EventModel<State>` | `EventModel<S>` | `callback_model_t` |
| 状态类型 | 任意 Python 对象 | 任意 C++ 类型 | 任意 Rust 类型 | `class(*)` 自定义派生类型 |
| 动态事件目录 | 支持 | 支持 | 支持 | 支持 |
| 替换随机源 | `uniform()` 对象 | `std::function<double()>` | `RandomSource` trait | 过程指针 `random` |
| 逐步执行 / 观察 | `step` / 生成器 | `step` / 回调 | `step` / 回调 | `step` / 回调 |

每一步速率可以依赖整个当前状态；事件间保持恒定。事件索引在一次计算与执行之间必须一致。自定义模型无需实现 `state()`、`position()`；它们只是内置模型的可选数值观测。任意字段可直接在观察函数中从模型读取，输出格式由你决定。具体扩展方式见 [API 文档](docs/api.md)。

算法为 Gillespie direct SSA / rejection-free KMC：

\[
A=\sum_i a_i(x),\quad \Delta t=-\log(u_1)/A,\quad P(i\mid x)=a_i(x)/A.
\]

零总速率返回吸收态；超过截止时间的事件不执行；负数、NaN、无限速率和不能推进的浮点时间被拒绝。内核每步 O(M) 扫描，M 为当前事件数，适合作为清晰、可扩展、经过核对的通用基础实现。

## 可选领域适配器与命令行

仓库另外附带有向图跳跃、质量作用反应网络，以及 Marcus / Arrhenius / Miller–Abrahams 速率函数。这些是应用层，不限制核心方法库。`KMC1` 文件目前只序列化这两个内置模型，**任意自定义模型通过库 API 接入**。

四个命令的参数相同：`MODEL [STEPS [SEED [UNTIL]]]`。默认 1000 事件、seed=1、无时间限制；`UNTIL` 是绝对模拟时间。CSV 写 stdout，结束状态写 stderr。

```bash
uv run ukmc examples/two_state.kmc 10000 42 10 > build/python.csv
./build/ukmc-cpp examples/two_state.kmc 10000 42 10 > build/cpp.csv
cargo run --quiet --manifest-path rust/Cargo.toml --release --bin ukmc-rust -- examples/two_state.kmc 10000 42 10 > build/rust.csv
./build/ukmc-fortran examples/two_state.kmc 10000 42 10 > build/fortran.csv
```

[文件格式](docs/format.md) · [算法与统计注意事项](docs/algorithm.md) · [参考来源与设计选择](docs/references.md) · [验证记录](docs/validation.md)

## 验证与当前边界

检查涵盖：任意状态与动态事件、禁用事件、指数等待时间与分支概率、截止时间、吸收态、非法输入、跨语言轨迹、通用排队例子、两态时间占据、反应守恒、固定时间扩散和出生死亡解析结果。

当前支持有限事件目录、分段恒定速率的 CTMC。持续随时间变化的 hazard、非指数等待时间、无限速率/Zeno 过程需要不同算法，不能直接套此内核。未实现 dependency graph/Fenwick tree 加速、MPI/GPU、通用检查点或 VOTCA 数据库导入。

默认 MINSTD 随机流主要用于四语言回归，可替换随机源；其周期和精度有限，长程生产模拟或并行独立流请接入合适的随机数生成器。Fortran 参数错误目前 `error stop 2`，其他版本分别使用异常或 `Result`。详见算法文档。

参考 VOTCA 的速率/事件分层与 pyMINDs 的输运理论，独立实现通用内核，没有复制上游源码。项目采用 MIT 许可。

## 真实体系验证

已使用 `e_006_ether_12-17_0` 的 100 个分子、550 对实际电子耦合，计算 Marcus 跃迁速率并完成四语言一致性和固定时间轨迹验证。输运是通用内核的一项应用验证；自定义问题仍通过上述两方法接口接入。

[体系报告与复现方法](runs/e_006_ether_12-17_0/README.md) · [机器可读结果](runs/e_006_ether_12-17_0/summary.json)

![真实体系 MSD 与精确 CTMC 对照](runs/e_006_ether_12-17_0/kmc_validation.png)
