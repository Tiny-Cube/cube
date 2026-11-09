# 载流子迁移率开发分支

分支：`feat/charge-carrier`。公共基础提交：`813959bb2cb171d140d3308e8db72d23c6234ba3`。

本分支的 input 只开启 `$charge-carrier`，其他计算入口关闭。现有参数准备链可用；速率、迁移率及收敛求解待接通。

本次初始化只建立独立开发空间及默认入口，不代表最终物理求解已经完成。公共接口修改先进入 `feat/calculation-entries`，本分支合并公共更新后继续专用开发；详见 [分支协作说明](development-branches.md)。
