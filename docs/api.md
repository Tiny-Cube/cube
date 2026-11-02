# 定义任意 KMC 模型

## 数学模型与程序接口

令状态为任意对象 `x`。每个事件 i 定义一个非负速率 `a_i(x)` 和变换 `T_i(x)`。引擎仅关心当前速率数组与按索引执行 `T_i` 的能力。状态可以包含连续数值，但它们只能在事件发生或外部规定的分段边界改变。若它们在事件间连续演化，需要积分 hazard 的算法。

`propensities()` 必须不改变物理状态，可以建立事件索引缓存；`fire(i)` 执行恰好一个事件。数组可为空，也可在下一步改变长度。禁用事件可以删除，或保留并返回零速率。概率不是速率：将所有速率乘以 c 会保持事件选择概率，但时间缩短为 1/c。

四种实现的事件 ID 都从 0 开始，包括 Fortran；数组在 Fortran 内部仍遵循通常的 1-based 下标。动态目录的 ID 只在当前步骤内有意义，如需稳定名称，请由模型记录。

## Python

最小模型实现 `propensities(self)` 和 `fire(self, event)`，不必继承任何类。引擎调用 `state()` / `position()` 只为创建可选数值观测；缺省返回空状态、零坐标。真实状态始终属于你的模型。

```python
from universal_kmc import Event, EventModel, Engine

state = {"available": ["a", "b", "c"], "completed": []}

def catalog(state):
    return [Event(
        rate=lambda s: 1.0,
        action=lambda s, name=name: (
            s["available"].remove(name), s["completed"].append(name)),
        name=f"complete-{name}",
    ) for name in state["available"]]

model = EventModel(state, catalog)
engine = Engine(model, seed=42)
for row in engine.run(max_steps=100, until=10.0):
    print(row.time, model.data)
```

事件的 action 原地修改 `data`；返回值被忽略。`enabled=False` 时不调用该事件的 rate。`EventModel` 会缓存本次工厂结果，执行后清空，避免对过期事件重复执行。`encode` 是可选的整型诊断向量，用于与内置 CSV 形式保持一致，不是对数据类型的限制。

`Engine.run` 是惰性生成器：必须迭代才能推进。初始记录和末端观测的 `event=-1`。`step(until=...)` 返回 `Record` 或 `None`；`engine.status` 为 `ready/running/absorbed/time_limit/step_limit`。`max_steps` 表示本次调用额外允许的事件数；`engine.steps` 是累计数。

若希望停止于任意状态条件，可以反复 `step()` 并检查自己的模型。若只关心最终结果，逐条丢弃记录即可，不必 `list()` 保存全部轨迹。

```python
while model.data['available']:
    record = engine.step(until=100.0)
    if record is None:
        break
```

替换随机源：

```python
import random
class RNG:
    def __init__(self, seed):
        self.r = random.Random(seed)
    def uniform(self):
        # 引擎要求严格开区间 (0,1)。
        while True:
            u = self.r.random()
            if u > 0.0:
                return u
engine = Engine(model, rng=RNG(42))
```

## C++17

包含 `ukmc.hpp`，派生 `ukmc::Model`，实现 `std::vector<double> propensities() const` 与 `void fire(std::size_t)`。模型生命周期必须覆盖引用它的 `Engine`。`run(max_steps, until, observer)` 流式调用观察函数；`snapshot` 是数值观测，任意状态直接从捕获的模型读取。

或包含 `ukmc_events.hpp` 使用 `EventModel<State>`：构造参数为初始状态、事件工厂、可选整数观察编码器。事件工厂返回 `std::vector<Event<State>>`，每个事件包含 rate/action/enabled/name。公开的 `model.data` 保留原始用户类型。

`Engine(model, seed, random_function)` 可接入自己的 RNG。模板库是 header-only，链接时无需 Python。`cpp/tests/test.cpp` 展示最小模型与动态任务目录；`examples/generic.cpp` 是独立排队系统。

## Rust

`Model` trait 只要求实现 `propensities(&self) -> Result<Vec<f64>>` 与 `fire(&mut self, usize) -> Result<()>`。引擎与模型分开存储，`engine.step(&mut model, until)` 返回 `Result<Option<Record>>`。`engine.run(&mut model, max_steps, until, observer)` 可观察记录；需要在每步读取任意用户字段时直接循环调用 `step`，避免 run 借用模型期间的 Rust 借用冲突。

`events::EventModel<S>` 用任意类型 S 保存数据，工厂返回闭包事件数组。`Event::new(rate, action)` 的 enabled 默认 true，也可替换。回调通过内部 RefCell 缓存目录；不是线程共享模型。每个线程使用自己的模型与引擎。

`Engine::with_rng(rng)` 接受实现 `RandomSource` 的类型。参见 `rust/examples/generic.rs` 与 `rust/tests/core.rs`。

## Fortran 2008

扩展抽象类型 `model_t`，实现 `propensities(self)`（allocatable real(dp) 数组）和 `fire(self,event)`（integer(i8) 的 0-based ID）。`state` 和 `position` 有默认实现，非必需。

也可以使用 `ukmc_events` 中的 `callback_model_t`：`allocate(model%data, source=your_state)`；分别设置 `model%rates => rates_callback` 与 `model%action => action_callback`。回调接受 `class(*)`，使用 `select type` 恢复你的类型。事件数和事件含义完全由回调决定。`encode` 是可选诊断回调。

`call engine%seed(42_i8)`；`advanced = engine%step(model, until, event)`；或 `call engine%run(model, max_steps, until, observer)`。所有公共事件 ID 从 0 开始。`engine%random => your_uniform` 可替换 RNG。回调错误目前使用 fail-fast 的 `error stop 2`；如果宿主应用要求错误恢复，需要在调用前自行验证，或扩展此错误策略。

参见 `examples/generic.f90`，它的状态包含动态分配的任务 ID 队列。

## 模型责任与观察

- rate/factory 不应调用随机数改变物理状态；否则会改变要模拟的随机过程。
- engine 不回滚自定义回调。action 抛错或失败后模型可能已部分改变，不应盲目继续。
- 回调可以在事件发生时采样额外随机标记，但应自行管理随机源与跨语言一致性。
- 不要在引擎执行期间从其他线程修改模型。每条轨迹独立拥有模型和 RNG。
- 事件后状态用于下一段驻留时间；计算时间平均时用前一状态乘以当前/前一记录的时间差。不能直接对事件记录做算术平均。
- 可选 `state()` 的整型向量只是简单观测通道；真实数值、字符串、图或对象可自行序列化，核心不会强制转换。
