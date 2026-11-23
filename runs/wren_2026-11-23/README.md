# Wren 2026-11-23 集成与 DOS 验收

集成载流子和模型 FSSH、更名 Python 包/命令为 wren，并加入三种 DOS 路线与 Sphinx/Furo 网页。
详见 validation.json；图中 H2 为功能算例，不是实际有机材料的 DOS。

136 项测试：135 通过、1 项显式 GPU opt-in 跳过。Sphinx 严格构建通过，浏览器验证导航、公式、
图标和中文搜索。已构建并在隔离环境安装 wheel，检查无科学依赖时的主命令/配置预检。
网页尚未验证公网部署。实际 MD 的 frame1980/frame2000、分子 1（86 原子）已完成 CPU 功能算例，用时约 452 s，
积分为 604.0000000000001（期望 604）。这是两帧、一个分子的 HF/STO-3G + Weigend 辅助基
真空轨道谱，不是已收敛的材料 DOS；两帧不足以估计统计不确定性。

原始本地输出保留在 build/dos，日志位于 /tmp/wren-full-tests.log 和 /tmp/wren-docs-build.log。
这些临时日志不是发行物，长期可追踪数值以 validation.json 与源测试为准。

## GPU 对照

GPU 报告为 gpu/validation.json，同机 86 原子 CPU/GPU 对照为 gpu/md-comparison.json。
CPU 48.62 s、RTX5090 3.92 s（8 线程、HF/STO-3G+Weigend、单次、GPU 预热）；
最大轨道差 0.835 meV。不要将不同主机的 CPU 两帧 452 s 与服务器 GPU 单帧直接相除。
GitHub 文档构建也已通过：https://github.com/Tiny-Cube/cube/actions/runs/45168889194 。

## 在线文档

Sites owner-private 发布已成功，文档源码提交 446923fe4c8b1ebc8d7c9dcf68fde39f2175bf2e，
部署标识 appgdep_6b042e8513248191abee4d43d72a680f。该版本仅所有者可访问，未改变为公开站点。
