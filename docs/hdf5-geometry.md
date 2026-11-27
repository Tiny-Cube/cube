# Wren 的 HDF5 提取格式

`$get_dimers_monomers` 默认使用 `format = hdf5`。给定 `output = mol`，所有选中帧写入 `mol/geometry.wren.h5`。也可以把 `output` 设为完整的 `.h5` 或 `.hdf5` 文件名。原始帧号以 0 为起点，`frame = -1` 仅保存最后一帧。即使某帧的 dimer 或 MM 环境为空，该帧仍有对应的索引范围。

根属性 `schema` 为 `wren-geometry/1`。位置和距离的单位是 Å，时间是 ps。一个归档使用同一套分子拓扑。每帧只保存一份完整分子坐标，已经按键还原跨周期边界的分子。

| 路径 | 内容 |
| --- | --- |
| `topology/molecule_id, molecule_offset` | 分子编号列表及各分子的原子切片边界。GROMACS 按拓扑 fragment 顺序从 1 编号，LAMMPS 保留 data 中的分子 ID。 |
| `topology/atom_id, global_atom_index, element, atom_name, charge_e, bonds` | 来源原子编号、从 0 开始的来源原子索引、元素、名称、电荷及键。键的端点按 HDF5 原子数组从 0 编号。 |
| `frames/id, time_ps, cell, origin, pbc` | 原始帧号、时间、按行排列的三个晶格向量、盒原点及三个周期轴标记。LAMMPS 另存 `lammps_timestep`。 |
| `geometry/positions` | 形状为（保存帧数，总原子数，3）。坐标保持输入读取器的浮点精度，按单帧分块，用 shuffle 和 gzip 无损压缩。 |
| `neighbors/frame_id` | 与 `frames/id` 指向同一数据集。 |
| `neighbors/monomer` | 每帧每个完整单体一行。含 `frame_ptr, qm_mol_id, mm_ptr, mm_mol_id, mm_image, mm_min_angstrom`。 |
| `neighbors/dimer` | 每帧每个符合 core 判据的分子对镜像一行。字段同 monomer，其中 `qm_mol_id` 每行含两个分子编号，另存 `qm_b_image, core_min_angstrom, core_com_angstrom`。 |

对任一张邻居表，第 `k` 个已保存帧的行范围为 `frame_ptr[k]:frame_ptr[k+1]`。第 `r` 行的 MM 范围为 `mm_ptr[r]:mm_ptr[r+1]`。两级指针均为整个文件内从 0 开始的绝对偏移，右端不包含在内。

每个 MM 邻居由分子编号和整数镜像向量共同确定。镜像坐标为 `positions + image @ cell`。对 dimer，第一个 QM 分子在原位，第二个应用 `qm_b_image`。MM 索引相对于第一个 QM 的晶胞，合并两个 QM 分子的环境并排除两个 QM 镜像本身。记录的 MM 距离是该完整邻居与任一完整 QM 分子间最近的原子距离。core 只用于判定 dimer。两个阈值均按严格小于比较。

以下示例读取第一个选中帧中第一个 dimer 的 QM 和 MM。如果该帧没有 dimer，行范围为空。

```python
import h5py

with h5py.File("mol/geometry.wren.h5") as h:
    d = h["neighbors/dimer"]
    start, stop = d["frame_ptr"][:2]
    if start < stop:
        qm = d["qm_mol_id"][start], d["qm_b_image"][start]
        a, b = d["mm_ptr"][start:start + 2]
        mm = list(zip(d["mm_mol_id"][a:b], d["mm_image"][a:b]))
```

`wren-xyz frame2000.mol8 --root mol -o mol8.xyz` 按帧号和分子编号直接读取单体。cube-view 可打开 HDF5，在界面中选择帧、单体或 dimer。量化工作流使用 `input = mol/geometry.wren.h5` 遍历已保存帧。原始 MD 和力场仍用于依赖它们的 QM/MM 计算。`--resume` 核对完整归档的来源、提取参数和内容校验和后复用结果。
