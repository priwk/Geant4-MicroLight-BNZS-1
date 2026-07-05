Stage D 原点对齐轨迹图独立脚本说明

用途：
本文件夹用于单独生成 stageD_origin_aligned_trajectory_showcase_3d.png。
服务器端先从 Stage B 原始轨迹 CSV 中提取前 10000 行原点对齐轨迹，
生成同目录下的 origin_aligned_particle_trajectories_first10000.csv。
之后把整个文件夹交给 Windows 10 用户，绘图脚本默认只读取同目录 CSV，
不需要 Stage B 原始数据，也不需要项目目录结构。

第三方 Python 库：
- numpy
- matplotlib

Python 标准库：
- argparse
- csv
- math
- os
- re
- collections
- dataclasses
- pathlib

安装命令：

cd plots/stageD_zns_trajectories/origin_aligned_showcase_standalone
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

Windows 10 绘图命令：

python generate_origin_aligned_showcase.py

默认绘图输入与输出：
- 输入 CSV：origin_aligned_particle_trajectories_first10000.csv
- 输出 PNG：stageD_origin_aligned_trajectory_showcase_3d.png

上述两个文件默认都相对脚本所在文件夹，不使用绝对路径。

服务器端重新提取 CSV 的命令：

python3 plots/stageD_zns_trajectories/origin_aligned_showcase_standalone/generate_origin_aligned_showcase.py --extract-from-stageb

服务器端提取默认值：
- --project-root 默认由脚本所在位置向上推断到项目根目录
- --stageb-root 默认是 Output/stageB
- --ratios 默认是 1-2
- --showcase-ratio 默认是 1-2
- --thicknesses 默认是 1000
- --showcase-thickness 默认是 1000
- --row-limit 默认是 10000

如果服务器目录结构不是：
<project-root>/plots/stageD_zns_trajectories/origin_aligned_showcase_standalone
和
<project-root>/Output/stageB
这种形式，可以显式指定：

python3 generate_origin_aligned_showcase.py --extract-from-stageb --project-root /path/to/project-root

或者直接指定 Stage B 根目录：

python3 generate_origin_aligned_showcase.py --extract-from-stageb --stageb-root /path/to/Output/stageB

服务器端输入文件读取规则与 generate_stageD_zns_trajectories.py 一致：
- 在 Output/stageB 下扫描比例文件夹，例如 1-2。
- 对每个选中厚度，优先读取 *_alpha_li_steps.csv。
- 如果完整 alpha/Li 文件不存在，则回退读取 *_zns_track_steps.csv。
- --showcase-ratio 与 --showcase-thickness 决定用于生成独立 CSV 的具体文件。

字体路径：
1. 本地字体路径
   如果要把字体一起交给 Windows 10 用户，可以放在：
   fonts/Arial.ttf
   fonts/Arial-Bold.ttf
   fonts/Arial-Italic.ttf
   fonts/Arial-BoldItalic.ttf

   在完整项目目录中运行时，脚本也会尝试读取：
   ../fonts/Arial.ttf
   ../fonts/Arial-Bold.ttf
   ../fonts/Arial-Italic.ttf
   ../fonts/Arial-BoldItalic.ttf

2. 系统字体回退
   如果上述本地字体不存在，matplotlib 会按以下顺序回退：
   Arial, Helvetica, Liberation Sans, DejaVu Sans
