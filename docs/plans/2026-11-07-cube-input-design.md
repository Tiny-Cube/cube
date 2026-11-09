# cube 分段控制文件与任务调度

用户确定入口脚本名为 `cube`，默认读取名为 `input` 的文本控制文件，使用 `$get_dimer/$two-mol/$one-mol/$reorganization-energy` 与 `$end`。安装后提供同名 console script，仓库提供可执行启动文件。输入顺序不决定执行顺序；按提取→单体→dimer/转移积分→重组能的依赖调度。相对路径以 input 所在目录为准，拒绝未知字段、重复段落/键和不支持的引擎。

按用户给的 extract_dimers.py 采用分子局部 1-based core 编号，严格最短 core 原子距离小于 cutoff（Å）；全分子经 TPR 键拓扑还原 PBC，再枚举整体分子镜像。支持一般非退化晶胞，完整保存平移与镜像编号；多个镜像不覆盖文件。输出只写 format 选择的 xyz/gjf，另保存机器可读 manifest、原子映射和距离，GJF 使用相应量化段的真实泛函/基组。无需复制或调用参考脚本中的 Gaussian/Slurm/MOMAP 配置。

独立居中的单体与 dimer 不能直接传给现有 DIPRO 校验。单体真空 SCF 只计算一次；在每个 dimer 中通过检查过的整体平移重定位 checkpoint 的 Mole 坐标，保持完整轨道和 Fock 不变，记录来源与平移。这是孤立单体平移不变性的应用，不能用于未说明的 QM/MM 环境。每个 pair 的 A/B/AB 在同一坐标框架下投影。

量化通过可独立调用的脚本/同源 Python API 调度，允许 CPU PySCF 与 gpu-pyscf；two-mol 默认继承 one-mol 的引擎，单体/dimer DIPRO 的方法和基组必须一致。现有 cube-transfer TOML 入口保留。逐任务记录输入哈希、设置、真实后端、成功/失败；--resume 只复用哈希与文件一致的成功结果。--plan 不启动量化，--only 可单独提取或运行指定阶段。

用户已明确这次先实现 **vacuum 孤立分子四点法**。中性单重态与阳/阴离子双重态分别优化并检查 SCF/几何优化收敛，计算 E0(R0), Eq(R0), Eq(Rq), E0(Rq)，输出两项松弛能和总和，保留几何、梯度、单位与优化轨迹。禁止静默把 MD 构型当作优化平衡构型，不称 QM/MM 或含外球重组能。GPU 请求必须实际使用 GPU 能量/梯度，优化器控制仍在 CPU。

验证包括严格解析与依赖规划、周期整体镜像/临界距离、格式与原子坐标一致性、单体 checkpoint 平移与耦合、四点能量/收敛失败、小分子 CPU/GPU 工作流。真实 anneal_a 输入用于完整提取与参考算法交叉验证；100 分子体系的全部 DFT/优化任务不作为快速测试启动。原始 66 MB 轨迹不加入仓库。
