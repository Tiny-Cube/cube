# 标准文献与软件接口

理论章节在相关公式附近引用来源。以下均为论文原始出版信息或软件官方文档，链接有效性
与实现能力验证是两回事；引用某软件不表示 Wren 已实现其全部功能。

| 主题 | 文献 |
|---|---|
| 非正交片段轨道与元素 PDOS | P.-O. Löwdin, *On the Non-Orthogonality Problem Connected with the Use of Atomic Wave Functions in the Theory of Molecules and Crystals*, JCP **18**, 365–375 (1950), [DOI](https://doi.org/10.1063/1.1747632) |
| KS 与电子添加/移除能量 | Perdew, Parr, Levy, Balduz, *Density-Functional Theory for Fractional Particle Number: Derivative Discontinuities of the Energy*, PRL **49**, 1691 (1982), [DOI](https://doi.org/10.1103/PhysRevLett.49.1691) |
| 单粒子/多体激发与谱 | Onida, Reining, Rubio, *Electronic excitations: density-functional versus many-body Green's-function approaches*, RMP **74**, 601 (2002), [DOI](https://doi.org/10.1103/RevModPhys.74.601) |
| k 点积分 | Blöchl, Jepsen, Andersen, *Improved tetrahedron method for Brillouin-zone integrations*, PRB **49**, 16223 (1994), [DOI](https://doi.org/10.1103/PhysRevB.49.16223) |
| 无序局域化 | Anderson, *Absence of Diffusion in Certain Random Lattices*, PR **109**, 1492 (1958), [DOI](https://doi.org/10.1103/PhysRev.109.1492) |
| 有机无序输运 | Bässler, *Charge Transport in Disordered Organic Photoconductors a Monte Carlo Simulation Study*, pss(b) **175**, 15–56 (1993), [DOI](https://doi.org/10.1002/pssb.2221750102) |
| 动态耦合无序 | Troisi, Orlandi, *Electronic disorder in crystalline organic semiconductors: the role of thermal motions modulating the transfer integrals in polyacenes*, [作者预印本](https://arxiv.org/abs/cond-mat/0511314) |
| 大型稀疏 DOS | Weiße, Wellein, Alvermann, Fehske, *The kernel polynomial method*, RMP **78**, 275 (2006), [DOI](https://doi.org/10.1103/RevModPhys.78.275) |
| Jackson 核谱重建 | Silver, Roeder, Voter, Kress, *Kernel Polynomial Approximations for Densities of States and Spectral Functions*, J. Comput. Phys. **124**, 115–130 (1996), [DOI](https://doi.org/10.1006/jcph.1996.0048) |
| Marcus | Marcus, *On the Theory of Oxidation–Reduction Reactions Involving Electron Transfer. I*, JCP **24**, 966 (1956), [DOI](https://doi.org/10.1063/1.1742723) |
| FSSH | Tully, *Molecular dynamics with electronic transitions*, JCP **93**, 1061 (1990), [DOI](https://doi.org/10.1063/1.459170) |
| 退相干修正 | Granucci, Persico, *Critical appraisal of the fewest switches algorithm for surface hopping*, JCP **126**, 134114 (2007), [DOI](https://doi.org/10.1063/1.2715585) |
| 激发态四点与梯度 | Hirata, Head-Gordon, *Time-dependent density functional theory within the Tamm–Dancoff approximation*, [DOI](https://doi.org/10.1016/S0009-2614(99)01149-5) |
| 垂直梯度与谐模型 | Avila Ferrer, Santoro, *Comparison of vertical and adiabatic harmonic approaches for the calculation of the vibrational structure of electronic spectra*, [DOI](https://doi.org/10.1039/C2CP41169E) |
| 凝聚相能隙采样 | *Efficient Calculation of Electrostatic Energies for Large-Scale Nonadiabatic Molecular Dynamics in a Site Basis*, [DOI](https://doi.org/10.1021/acs.jctc.5c01753) |
| PySCF | Sun et al., *Recent developments in the PySCF program package*, JCP **153**, 024109 (2020), [DOI](https://doi.org/10.1063/5.0006074) |
| GPU4PySCF | *Introducing GPU Acceleration into the Python-Based Simulations of Chemistry Framework*, [DOI](https://doi.org/10.1021/acs.jpca.4c05876) |
| 周期高斯密度拟合 | Sun et al., *Gaussian and plane-wave mixed density fitting for periodic systems*, JCP **147**, 164119 (2017), [DOI](https://doi.org/10.1063/1.4998644) |
| 周期快照能量对齐 | Van de Walle 与 Martin，[平均势与能带对齐](https://doi.org/10.1103/PhysRevB.35.8154)；Hartwigsen 等，[局域与非局域赝势](https://doi.org/10.1103/PhysRevB.58.3641) |
| 时间采样与统计误差 | Flyvbjerg 与 Petersen，[Error estimates on averages of correlated data](https://doi.org/10.1063/1.457480) |
| 光学与电子–空穴谱 | Rohlfing 与 Louie，[Electron-hole excitations and optical spectra from first principles](https://doi.org/10.1103/PhysRevB.62.4927) |
| 有机晶体的电子–声子效应 | Brown-Altvater 等，[晶体萘的带隙重整化与电子–声子自能](https://doi.org/10.1103/PhysRevB.101.165102) |
| 萘晶胞的实验坐标 | Fabbiani et al., *Exploration of the high-pressure behaviour of polycyclic aromatic hydrocarbons*, Acta Cryst. B **62** (2006), [结构、压力与坐标表](https://doi.org/10.1107/S0108768106026814) |

官方接口：[PySCF PBC](https://pyscf.org/user/pbc/scf.html)、
[PySCF QM/MM](https://pyscf.org/user/qmmm.html)、
[ORCA ESD](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/esd.html)、
[GDMA 作者仓库](https://gitlab.com/anthonyjs/gdma)、
[VOTCA 文档](https://www.votca.org/)。
