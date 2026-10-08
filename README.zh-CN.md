# 边缘非线性 Hall 输运

[English](README.md)

本仓库整理了 **Nonlinear Hall Effect from Inequivalent Edge Bands in Bulk Insulators** 四张主图的 Python 代码与已保存的数值数据，包括十三个独立子图流程、Rice–Mele（RM）与量子自旋 Hall（QSH）共享输运后端，以及四张参考 PDF。运行所需输入均位于仓库内部。

## 用现成数据绘图

在仓库根目录运行：

```sh
python -m pip install -r requirements-plot.txt
python code/run_panel_workflow.py --stage plot
```

该命令绘制十三个子图并组装四张整图，读取 `data/panels/`，将 PDF、PNG、SVG 和元数据写入 `outputs/`。绘图只需要 NumPy 和 Matplotlib，不加载输运求解器。`figures/reference/` 中保留原参考图。

只画一个子图，或只组装一张整图：

```sh
python code/run_panel_workflow.py --figure 3 --panel c --stage plot
python code/run_panel_workflow.py --figure 4 --stage assemble
```

选择单个子图时不组装整图。数据缺失或与 JSON 定义不符时会报错；绘图入口不会自动启动计算。

## 执行阶段

| 阶段 | 操作 |
|---|---|
| `plot` | 读取现成数据绘图；未指定 `--panel` 时还会组装所选整图 |
| `assemble` | 直接从该图全部子图数据组装整图 |
| `data` | 从已保存的输入缓存准备子图数据，不绘图 |
| `all`（默认） | 准备数据、逐子图绘图、组装整图 |

```sh
python code/run_panel_workflow.py
```

默认执行缓存提取，会写入所选的 `data/panels/` NPZ/JSON 文件，不重跑输运扫描。只有显式指定 `--recompute` 才启动新计算。最近一次运行记录为 `outputs/last_run.json`。

## 将重算结果另存

先按[环境说明](docs/installation.md)安装计算依赖。建议通过子图入口将结果保存到独立位置：

```sh
python code/panels/fig2/calc_a.py --recompute --output outputs/recomputed/fig2/a.npz
python code/panels/fig2/plot_a.py --data outputs/recomputed/fig2/a.npz --output outputs/recomputed/fig2/a
```

计算同时生成 `a.npz` 和 `a.json`；绘图读取这对文件。总入口的 `--recompute` 使用默认 `data/panels/` 位置，会覆盖所选的随仓库提供的子图数据，因此比较新旧结果时请使用上述另存方式。输运、温度和无序扫描可能耗时较长；每个子图只执行自身需要的扫描。Fig. 1 重新生成定性示意图，不执行输运求解。

## 目录与图件

| 位置 | 内容 |
|---|---|
| `code/run_panel_workflow.py` | 选择图号、子图和执行阶段 |
| `code/panels/` | 各子图独立计算与绘图入口 |
| `code/*.py` | 共享模型与测量后端 |
| `code/checks/` | 补充数值诊断 |
| `data/panels/` | 可直接绘图的数组及定义 |
| `data/raw/` | 用于准备子图的原始缓存 |
| `figures/reference/` | 四张未改动的参考 PDF |
| `provenance/source_map.json` | 来源、路径迁移和 SHA-256 记录 |
| `scripts/` | 发布包验证工具 |
| `docs/` | 环境、流程、约定和数据来源说明 |
| `outputs/` | 新生成结果，Git 忽略此目录 |

Fig. 1 为装置示意；Fig. 2 展示 RM 能带、边缘谱密度与 LDOS、左入射 LPDOS 和二阶 Hall 响应；Fig. 3 展示 QSH 能带及温度、非磁无序、保自旋弹性虚探针扫描；Fig. 4 比较 RM/QSH 的探针耦合及弱探针无序响应。文稿源码和投稿信不包含在本包中。详细命令见[流程说明](docs/workflow.md)和[子图目录说明](code/panels/README.md)。

## 物理约定与来源

采用 `V1=+V/2`、`V2=-V/2`，其中 `V` 是完整源漏电压，横向电压为 `V3-V4`。无量纲二阶系数在 RM 中为 `t2*kappa2/e`，在 QSH 中为 `t*kappa2/e`。历史原始 `kH` 在计算侧仅转换一次：`-kH/4`；可绘图数据不再转换。Fig. 3(c) 误差棒为八个无序构型的 SEM，Fig. 4(c,d) 阴影为十六个构型的样本 SD。详见[约定说明](docs/conventions.md)。

参考 PDF 和保存的数值数组来自原工作区；独立包调整了路径及入口，原始与发布版本的哈希记录在 `provenance/source_map.json`。部分历史生成脚本与当前源码存在版本差异，JSON 保留了这些记录。缓存校验与提取复现已保存的输入，并不等于当前求解器的新数值验证。比较重算结果前请阅读[来源说明](docs/provenance.md)。

实际测试运行环境为 **Python 3.13.5（Anaconda）**，NumPy 2.1.3、SciPy 1.15.3、Matplotlib 3.10.0、Kwant 1.5.0、tinyarray 1.2.5、threadpoolctl 3.5.0。依赖文件记录这些版本。字体与绘图库版本可能影响新导出图件的外观和文件字节。

## 验证发布包

```sh
python scripts/verify_release.py
python scripts/verify_release.py --relocate --plot
```

第一条检查随包数据和发布记录。第二条还将仓库复制到独立临时目录，在禁止求解器导入的条件下准备十三份数据、绘制十三个子图并组装四张整图。报告写入 `outputs/validation/`。可选的 `--compute` 仅计算每个模型一个清洁弱探针代表点，需要完整计算依赖。详见[验证范围](docs/verification.md)；这些检查不重跑完整科学扫描。

## 许可证与引用

原始代码采用 [MIT 许可证](LICENSE)，版权为 2026 H. Geng and contributors。依赖另行安装并遵循各自许可证，本仓库不内置依赖库源码。

[CITATION.cff](CITATION.cff) 使用用户提供的四位作者姓名缩写，并记录[公开仓库地址](https://github.com/genghaophy/edge-nonlinear-hall)。软件 DOI 和文章发表信息尚未确定。
