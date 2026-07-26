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

`*_zns_track_steps.csv` records ZnS energy deposition from all transported
particles. The original `x_*_um` columns are folded RVE coordinates;
`unwrapped_*` columns are continuous across periodic cells, and `screen_*`
columns place those steps back in the finite-screen coordinate system.

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

逐次重入诊断默认关闭。调试时可通过宏命令开启并限制输出规模：

```text
/cfg/stageD/setWriteReentryDiagnostics true
/cfg/stageD/setReentryDiagnosticsSamplingRate 0.01
/cfg/stageD/setMaxDiagnosticRows 100000
```

对应环境变量为 `BNZS_STAGED_WRITE_REENTRY_DIAGNOSTICS`、
`BNZS_STAGED_REENTRY_DIAGNOSTICS_SAMPLING_RATE` 和
`BNZS_STAGED_MAX_DIAGNOSTIC_ROWS`。

逐光子 `stageD_events.csv` 在生产运行中默认关闭，summary 与相函数不受影响。
调试时可开启确定性抽样并限制行数：

```text
/cfg/stageD/setWriteEventCsv true
/cfg/stageD/setEventSamplingRate 0.01
/cfg/stageD/setMaxEventRows 100000
```

批处理命令对应 `--write-event-csv --event-sampling-rate 0.01 --max-event-rows 100000`。

```text
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/
  <source>__<boundary>__<config_hash>/run_<run_id>_seeds_<seed0>_<seed1>/
```

`server_parallel_run.py staged` 和 `stageD_run_batch.py` 会为每个 placement 写入不同的
Geant4 随机种子，避免重复批次误覆盖。Stage D CSV 使用 15 位有效数字。
正式均匀化默认使用 `uniform_all_phase + periodic_wrap + particle_encounter_no_threshold`。
单 run 同时输出 raw 与 thresholded 相函数；多 placement 后处理使用 ratio-of-sums，
并分别导出原生 RVE 参数和 `mu_s = mu_s_prime_raw, g = 0` 的各向同性等效参数。
`stageD_g0_equivalent_*` 保持 full-path 兼容语义，同时新增
`stageD_g0_full_path_*` 与 `stageD_g0_post_first_*`。后者描述完成一次完整颗粒
encounter 后的 moving-photon 子群，不能替代 ZnS 初始困光/提取模型。

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
