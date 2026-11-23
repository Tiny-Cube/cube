# 安装与命令

需要 Python 3.11 或更新版本，使用 uv 管理 Python 包和可选依赖。
新集成工作区位于 `wren`，命令均在该仓库根目录运行。

```bash
uv sync --locked --extra workflow --extra transport --extra dynamics --extra dos
uv run wren --help
uv run wren-dos --help
uv run wren-transfer --help
uv run wren-surface-hopping --help
uv run wren-kmc --help
```

GPU 分子量化使用锁定的 GPU4PySCF/CuPy/CUDA 依赖，需要实际兼容的 NVIDIA GPU 与驱动：

```bash
uv sync --locked --extra workflow --extra quantum-cuda13
```

请求 GPU 后不静默回退 CPU。当前周期 DOS 的 GPU SCF 尚未验证，输入要求显式使用 CPU；
单体 DOS 可指定 `qc-engine = gpu-pyscf`，DOS 数值展宽/有效 Hamiltonian 对角化可以独立选择
`dos-backend = gpu`，但其实际设备一致性验收仍须记录。

## 从源码构建与测试

```bash
uv build
uv run python scripts/build.py
uv run python -m unittest discover -s tests -v
```

四语言一致性测试需要 C++ 编译器、Rust/cargo 和 gfortran。Python 科学功能不要求构建这些程序。
历史输出中的 `cube-*` schema 和 KMC1 格式保留，不能简单按产品名称改写已有数据。

## 构建网页文档

```bash
uv sync --locked --extra docs
uv run sphinx-build -W --keep-going -b html docs/site build/html
uv run python -m http.server 38501 --directory build/html
```

网页含侧栏、搜索、代码块与数学公式。公式由 MathJax 显示，默认需要能访问其脚本来源。
源文档随软件版本管理。部署流程可使用相同静态构建产物，不需要将用户科学数据发布到网站。
