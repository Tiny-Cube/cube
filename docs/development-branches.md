# 三类计算的并行开发

公共入口路由维护在 `feat/calculation-entries`。三类计算从同一个公共基础提交创建独立分支：

| 分支 | 用途 | 独立工作目录 |
|---|---|---|
| `feat/calculation-entries` | 入口解析、模块归属与复用、统一报告、注册协议 | `cube` |
| `feat/charge-carrier` | 载流子速率、KMC、迁移率与收敛报告 | `cube-charge-carrier` |
| `feat/exciton` | 激发态参数、能量转移与激子扩散 | `cube-exciton` |
| `feat/surface-hopping` | 多态能量/梯度、非绝热耦合与动力学 | `cube-surface-hopping` |

分支对应的 input 默认只开启本分支入口。建立独立 Git worktree 后，各目录拥有自己的工作文件和分支，可同时编辑与提交；每个目录执行 `uv sync --locked --extra workflow` 建立自己的环境，不共用指向其他工作树的可编辑安装。三个分支目前是公共框架及入口配置，不代表三种最终求解均已完成。

```bash
git worktree list
# 其他克隆可以自行创建独立工作目录，例如：
git worktree add ../cube-exciton feat/exciton
```

公共 API 的修复先提交到 `feat/calculation-entries`，各计算分支再合并公共分支；计算专用提交留在自己的分支，验证后再集成。数据、虚拟环境和量化输出不提交。并行编辑不同工作树不会覆盖文件，但涉及公共接口时仍应协调和回归测试。

## 注册新计算类型

`src/universal_kmc/workflow/calculations.py` 是入口注册表。新类型声明自己允许调用的模块，可提供延迟加载的 solver：

```python
from .calculations import CalculationType, register_calculation

# 示例仅演示接口，不是现有物理求解实现。
def solve_new_calculation(control, preparation_report, *, resume):
    from .new_calculation import solve
    return solve(control, preparation_report, resume=resume)

register_calculation('new-calculation', CalculationType(
    modules=('get_dimers_monomers', 'one-mol'),
    solver=solve_new_calculation,
))
```

注册在 CLI 构造参数之前完成，`$new-calculation / do=1 / $end` 和 `--only new-calculation` 自动可用，不需要逐处修改入口名字列表。模块按公共固定顺序去重；solver 在模块成功后执行，返回包含 `status=completed/partial` 的 JSON 可序列化字典。科学依赖在 solver 内导入，以保持解析与关闭入口无计算副作用。

只提供模块而未提供 solver 时，完整入口结束为 `prepared`；不支持整条计算时设置 `unavailable_reason` 并在预检说明。新增公共模块本身仍需实现其解析和调度，不能只在注册表中添加未知模块名称。
