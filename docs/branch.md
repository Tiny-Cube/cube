# 激子传输开发分支

分支：`feat/exciton`。公共基础提交：`813959bb2cb171d140d3308e8db72d23c6234ba3`。

本分支的 input 只开启 `$exciton`，其他计算入口关闭。结构准备可用；激发态参数、激子耦合与扩散求解待实现。

本次初始化只建立独立开发空间及默认入口，不代表最终物理求解已经完成。公共接口修改先进入 `feat/calculation-entries`，本分支合并公共更新后继续专用开发；详见 [分支协作说明](development-branches.md)。
