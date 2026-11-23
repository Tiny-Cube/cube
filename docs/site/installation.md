# 安装与命令

需要 Python 3.11 或更新版本，使用 uv 管理 Python 包和可选依赖。
新集成工作区位于 `wren`，命令均在该仓库根目录运行。

当前实现发布在现有 GitHub 仓库的 `feat/wren` 分支，原 `main` 和开发分支保留：

```bash
git clone --branch feat/wren https://github.com/Tiny-Cube/cube.git wren
cd wren
```

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

单体与周期 DOS 均可指定 `qc-engine = gpu-pyscf`。数值展宽、元素投影与有效 Hamiltonian 对角化可独立选择 `dos-backend = gpu`。周期量化使用 GPU4PySCF 的 PBC 接口，运行记录包含实际 SCF 类和显卡型号。

## 用 CMake 构建和安装 Wren

准备 CMake 3.25 或更新版本、C++17 编译器、Fortran 编译器、Cargo/Rust、Python 3.11+ 和 uv。默认构建编译 C++、Fortran 和 Rust，调用 uv 生成 Python wheel。已有 Ninja 时可在首条 CMake 命令添加 `-G Ninja`。在 Linux 仓库根目录运行：

```bash
uv sync --locked --extra workflow --extra dos --extra dynamics --extra transport
cmake -S . -B build/cmake -DCMAKE_BUILD_TYPE=Release -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build build/cmake --parallel
ctest --test-dir build/cmake --output-on-failure
cmake --install build/cmake --prefix "$HOME/.local"
```

安装前缀的 `bin/` 包含主命令 `wren`、Python 入口、`ukmc-cpp`、`ukmc-rust` 和 `ukmc-fortran`。Python 启动器使用 PATH 中的 `python3`；计算 DOS 或 QM/MM 时先激活装有相应科学依赖的环境。也可以直接在选定环境安装 `build/cmake/dist/` 下的 wheel。C++ 头文件与当前 Fortran 编译器生成的模块文件位于前缀的 `include/wren`。CMake 消费工程可通过 `find_package(Wren CONFIG REQUIRED)` 引用 `Wren::ukmc_cpp` 和 `Wren::ukmc_fortran`。

| 需要的产物 | 命令或配置 |
|---|---|
| wheel 和源码包 | `cmake --build build/cmake --target wren_package` |
| 只重新构建 Python wheel | `cmake --build build/cmake --target wren_python` |
| 四语言计时所需的原生程序 | `cmake --build build/cmake --target wren_benchmarks` |
| 同时构建并安装网页使用手册 | 配置时加 `-DWREN_BUILD_DOCS=ON` 并在所选 Python 环境安装 `docs` extra |
| 只构建 C++/Fortran | 配置时加 `-DWREN_BUILD_PYTHON=OFF -DWREN_BUILD_RUST=OFF -DBUILD_TESTING=OFF` |
| 只使用 Python | `uv sync --locked`，按计算类型选择 extras |

单配置生成器的原生程序在 `build/cmake/bin/`。可用 WREN_BUILD_CPP、WREN_BUILD_FORTRAN、WREN_BUILD_RUST 和 WREN_BUILD_PYTHON 四个 CMake 选项选择语言。不同编译器或 Debug/Release 安装请各用一个前缀，以保持 Fortran 模块与库匹配。多配置生成器构建时使用 `--config Release`、测试时使用 `-C Release`、安装时使用 `--config Release`。CTest 会选择本次构建的程序做 Python、C++、Rust、Fortran 的逐事件比较。原 `scripts/build.py` 和 `scripts/check.py` 调用同一套 CMake 构建。GPU 量化依赖仍由 uv 的 `quantum-cuda13` extra 提供。

## 构建网页文档

```bash
uv sync --locked --extra docs
cmake -S . -B build/manual -DWREN_BUILD_DOCS=ON -DBUILD_TESTING=OFF \
  -DWREN_BUILD_CPP=OFF -DWREN_BUILD_FORTRAN=OFF -DWREN_BUILD_RUST=OFF \
  -DWREN_BUILD_PYTHON=OFF -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build build/manual --target wren_docs
uv run --no-sync python -m http.server 38501 --directory build/manual/html
```

生成的入口为 `build/manual/html/index.html`。可选执行
`cmake --install build/manual --prefix "$HOME/.local" --component Documentation` 安装手册。
网页含侧栏、搜索、代码块与数学公式。公式由 MathJax 显示，默认需要能访问其脚本来源。
源文档随软件版本管理。部署流程可使用相同静态构建产物，不需要将用户科学数据发布到网站。
