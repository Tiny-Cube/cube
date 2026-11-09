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

`src/universal_kmc/workflow/calculations.py` 是主模块注册表。新类型声明自己自动组织的子模块，可提供延迟加载的 solver。这不是独立子模块的权限表：任何已有子模块都可通过自身 `do=1` 独立运行，不要求主模块开启。

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

子模块开关使用统一契约：`do=1` 添加独立执行请求；`do=0` 禁用并阻止隐式调用；省略 `do` 保留自动调度。主流程和独立请求共用科学实现、报告和续算。公共改动见 [主模块／子模块设计](plans/2026-11-09-main-submodule-design.md)。

主模块段内可使用其模块名作为 0/1 选择，例如 `one-mol=0` 使该主流程不调用单体；独立子模块 `do=1` 不受这个选择限制。显式主请求 1 与子模块 do=0 冲突时，在执行前报错。详见 [主模块内部选择](plans/2026-11-09-main-module-selectors-design.md)。
