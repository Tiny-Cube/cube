# 簇查看器已独立迁移

簇查看器的实现、测试、使用说明与性能记录已迁至独立项目 [Tiny-Cube/cube-view](https://github.com/Tiny-Cube/cube-view)。cube 不再包含查看器代码。

在 cube-view 目录中执行：

```bash
uv sync --locked
uv run cube-view /path/to/structure.gro --cutoff 6 --open
```

[完整使用说明](https://github.com/Tiny-Cube/cube-view/blob/main/docs/usage.md) · [迁移验证](https://github.com/Tiny-Cube/cube-view/blob/main/docs/migration.md)
