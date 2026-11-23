# Wren 2026-11-23 集成与 DOS 验收

集成载流子和模型 FSSH、更名 Python 包/命令为 wren，并加入三种 DOS 路线与 Sphinx/Furo 网页。
详见 validation.json；图中 H2 为功能算例，不是实际有机材料的 DOS。

136 项测试：135 通过、1 项显式 GPU opt-in 跳过。Sphinx 严格构建通过，浏览器验证导航、公式、
图标和中文搜索。已构建并在隔离环境安装 wheel，检查无科学依赖时的主命令/配置预检。
网页尚未验证公网部署。实际 MD 两帧的单分子 CPU 功能算例在运行，不计为已完成材料 DOS。

原始本地输出保留在 build/dos，日志位于 /tmp/wren-full-tests.log 和 /tmp/wren-docs-build.log。
这些临时日志不是发行物，长期可追踪数值以 validation.json 与源测试为准。
