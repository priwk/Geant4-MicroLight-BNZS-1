# Geant4-MicroLight-BNZS

BN/ZnS 周期微结构 Geant4 模型。

## 运行模式

```text
/cfg/setRunMode StageA_NeutronPatch
/cfg/setRunMode StageB_ReplayAlphaLi
/cfg/setRunMode StageD_OpticalHomogenization
```

Stage C 的 Geant4 运行模式、动作实现、宏和批处理入口已经移除。Stage D 的名称保持为
`StageD_OpticalHomogenization`。

## 源码结构

```text
src/core      通用配置、周期 RVE 几何、物理列表、动作初始化
src/modes     A/B/D 运行模式分发
src/stageA    热中子局部区域输运与 Sigma_eff 汇总
src/stageB    从中子俘获 CSV 回放 alpha/Li7
src/stageD    周期 RVE 光学均匀化
```

## 周期 RVE

程序从 `Input/output_pbc` 选择 `format_version=3` 主 placement，并排除
`*_pbc_images.csv` 和 `*_radius_stats.csv`。主文件中的 `box_x_um`、`box_y_um`、
`box_z_um` 是 v3 几何尺寸来源；主文件用于唯一颗粒统计和源采样，PBC images 文件用于
Geant4 physical placement。

批处理脚本按目录中实际主 placement 数量运行，不设 50 套上限；`--max-placements 0`
表示全部。因此同一配比增加到 100 套时会自动处理 100 套，辅助 CSV 不计入套数。

默认模式：

```text
/cfg/setPlacementGeometryMode pbc_clipped
/cfg/setPeriodicImagesCsv <optional-explicit-path>
```

也可使用环境变量：

```text
BNZS_PLACEMENT_GEOMETRY_MODE=pbc_clipped
BNZS_PBC_IMAGES_CSV=<path>
```

`primary_only` 仅用于调试。默认关闭 Geant4 overlap check；需要诊断时可设置
`BNZS_CHECK_OVERLAPS=1`。启动时不再执行全颗粒周期重叠扫描、周期边界探针和重复 copy
映射复核，以缩短大 RVE 初始化时间。CSV 必填字段、ID、相位、半径类别、尺寸和统计量的
轻量校验仍保留。

## Stage A

对一个或多个配比目录运行全部排布：

```bash
STAGEA_EVENTS=100000 bash stageA_batch_placements.sh 1-2 1-3
```

输出位于：

```text
Output/stageA/<ratio>/neutron_transport_summary.csv
```

## Stage B

常用批处理入口：

```bash
bash batch_run.sh 1-2
```

推荐生产输出模式：

```bash
BNZS_STAGEB_OUTPUT_MODE=slim bash batch_run.sh 1-2
```

主要输入和输出：

```text
Input/stageA/<ratio>/neutron_capture_positions/*_neutron_capture_positions.csv
Input/output_pbc/<ratio>/placement_*.csv
Output/stageB/<ratio>/<thickness>_capture_anchors.csv
Output/stageB/<ratio>/<thickness>_zns_track_steps.csv
Output/stageB/<ratio>/<thickness>_boundary_stop_summary.csv
```

## Stage D

示例：

```bash
cd build
./Geant4-MicroLight-BNZS StageD_OpticalHomogenization.mac
```

运行模式名称与宏文件名称保持不变：

```text
StageD_OpticalHomogenization
StageD_OpticalHomogenization.mac
```

Stage D 默认使用 `periodic_wrap`：光子越过人工 RVE 面后，在对面同一周期坐标
继续，方向、偏振、能量、时间和权重保持不变。旧的 `same_phase_reentry`
统计重采样模式仍保留用于结果兼容。输出位于：

```text
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/
```

后处理：

```bash
python3 scripts/merge_stageD_optical_params.py --ratio 1-2
python3 scripts/calibrate_optical_params_with_Leff.py \
  --ratio 1-2 \
  --experimental-leff experimental_Leff.csv
```

## 构建

```bash
cmake -S . -B build -DWITH_GEANT4_UIVIS=OFF
cmake --build build -j$(nproc)
```
