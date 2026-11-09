# 已有 PySCF 结果的控制文件接口

用户要求通过控制文件指定单体与 dimer 的 PySCF 输出，再计算转移积分，并要求另写通用 PySCF 脚本。这里把量化计算与后处理分开：`scripts/run_pyscf.py` 接受任意单体/dimer XYZ，生成各自 `.chk/.log/.status.json`；`scripts/calculate_transfer.py transfer.toml`（等价于 `cube-transfer transfer.toml`）只读取结果并执行投影，不调用 SCF kernel。

采用标准 PySCF HDF5 `.chk`，保持原生读取兼容性。普通 `.out/.log` 通常不包含完整 MO 系数，不能据两个前线轨道能量可靠恢复 DIPRO。单独设计 NPZ 交换格式会增加约定和维护成本，故首版不采用。输出模板在标准 checkpoint 中附加收敛标记、模型名称、独立 Fock 矩阵，避免仅凭文件存在判断 SCF 成功。标准原生 checkpoint 缺少收敛状态时，必须显式 `assume_converged=true`，并在结果中标明这是用户确认而非程序验证。

TOML 路径相对于控制文件目录解析；默认 A+B 原子顺序，也支持按单体原子顺序给出 dimer 的完整原子编号。编号从 1 开始；可选 MO 编号同样从 1 开始。检查实际坐标、元素、基组壳层、AO 顺序、闭壳层占据、MO 完整性和 S 正交性。重新优化/居中后的单体、ghost 基、ECP、开壳层、截断 MO 等不在首版范围内。

模板输出保存的 Fock 直接用于投影，完整 MO 展开独立核对其一致性。标准 checkpoint 没有保存 Fock 时，利用 `F = S C diag(epsilon) C.T S` 恢复二聚体算符，共用现有 DIPRO/ESID/orbital-GMH 投影代码。恢复恒等式不是独立的 SCF 收敛验证，结果明确说明。标准 checkpoint 未保存 functional 等元数据时，不用默认 B3LYP 设置冒充原始方法。冻结密度方法需要另一份 Hamiltonian，不能由收敛 AB checkpoint 代替，所以此入口拒绝 frozen_density。

`--check` 只做 TOML/路径/设置预检，不导入 NumPy/PySCF、不探测 GPU、不读取矩阵。真正计算再核对全部矩阵和原子映射。原命令行 XYZ→SCF 入口保持兼容。控制文件不允许混用命令行覆盖，防止参数来源不明；输出不得覆盖输入或控制文件，写出继续使用原子替换。

验证包括依赖无关的控制解析、错误输入与路径保护、真实 CPU/GPU 生成 checkpoint 后的读取、同一计算直接投影对照、原子交错顺序、显式 MO、多方法、不收敛和破坏矩阵等拒绝情形。GPU 用于后处理的投影/正交化，读取 HDF5、重建特征算符和 libcint AO 积分在 CPU；SCF 的 GPU 加速由输出模板控制。
