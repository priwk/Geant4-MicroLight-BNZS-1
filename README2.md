# StageD Notes

## 1. 当前结论

`StageD_OpticalHomogenization` 现在可以使用，但应理解为：

- 可以用于提取 `raw / provisional optical parameters`
- 目前输出的是：
  - `mu_a_raw`
  - `mu_s_raw`
  - `g_raw`
  - `mu_s_prime_raw`
- `n_eff_initial` 还没有放进 C++ 主程序输出，建议先在后处理里计算
- 目前还不应把这组参数直接当作“实验校正后的最终参数”

StageD 当前更适合做：

- 多个 placement 的 ensemble raw parameter extraction
- 不同 ratio 的相对比较
- 后续 Macro OpticalMC 的前置输入

不适合直接宣称：

- 已完成实验标定
- 已得到最终真实材料参数

## 2. StageD 现在到底做了什么

StageD 的目标不是宏观厚度出光模拟，而是对单个随机 RVE 做光学输运均匀化统计。

当前实现：

- `runMode = StageD_OpticalHomogenization`
- `source_mode = uniform_all_phase`
- `boundary_mode = periodic_wrap`
- `reentry_mode = state_matched`
- `particle_reentry_mode = sphere_q_mu`
- `matrix_reentry_mode = clearance_binned_portal`
- 输出逐 photon 事件表
- 输出单 run summary 表

兼容模式：

- `same_phase_reentry` 保留旧统计重入算法，仅用于复现旧结果

## 3. 边界模型说明

当前默认是严格周期坐标续传。

当前实现的是：

- 光子在人工 RVE 面离开后，从对面同一周期坐标继续
- 方向、偏振、能量、时间和权重保持不变

含义：

- 一个 Geant4 run 只使用一个 placement
- 不执行随机 portal 重采样
- 不跨 placement
- 要求 `format_version=3` 且 `pbc_clipped` 几何

因此：

- `one run = one placement`
- `many runs = many placements`
- 最终参数应来自多 placement 的 ensemble average

## 4. 当前使用建议

推荐按下面顺序使用：

1. 固定一个 ratio
2. 一次 run 只用一个 placement
3. 每个 placement 运行足够多 photon，例如 `beamOn 10000`
4. 对多个 placement 汇总 raw 参数
5. 后续如有实验，再做 `L_eff` 校正

## 5. 最常用命令

通常你的工作目录已经是项目根目录，例如：

```bash
cd B2
```

如果你通常是“跑一个配比下的所有 placement”，最常用的命令就是：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --beam-on 10000 \
  --run
```

跑完后，直接汇总这个 ratio 下所有已经完成的 placement：

```bash
python3 scripts/merge_stageD_optical_params.py \
  --ratio 1-2 \
  --project-root .
```

然后看合并结果：

```bash
cat Output/optical_params/1-2/rve_raw_optical_params_by_ratio.csv
```

如果想一条命令连续执行“全部 placement + 自动汇总”：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --beam-on 10000 \
  --run && \
python3 scripts/merge_stageD_optical_params.py \
  --ratio 1-2 \
  --project-root . && \
cat Output/optical_params/1-2/rve_raw_optical_params_by_ratio.csv
```

## 6. 单个 placement 的运行方法

在 `build/` 目录下运行。

示例：

```bash
cd build

cat >/tmp/stageD_run.mac <<'EOF'
/run/verbose 0
/event/verbose 0
/tracking/verbose 0

/cfg/setRunMode StageD_OpticalHomogenization
/cfg/setWeightRatio 1 2
/cfg/setPlacementFile ../Input/output_pbc/1-2/placement_phi0.63998_seed1245082203.csv

/cfg/stageD/setWavelengthNm 450
/cfg/stageD/setSourceMode uniform_all_phase
/cfg/stageD/setBoundaryMode periodic_wrap
/cfg/stageD/setReentryMode state_matched
/cfg/stageD/setParticleReentryMode sphere_q_mu
/cfg/stageD/setMatrixReentryMode clearance_binned_portal
/cfg/stageD/setTargetPrimaryScatter 0
/cfg/stageD/setThetaThresholdDeg 0.1
/cfg/stageD/setMaxReentry 10000
/cfg/stageD/setMaxSteps 100000
/cfg/stageD/setMaxPathLengthUm 1000000

/run/initialize
/run/beamOn 10000
EOF

./Geant4-MicroLight-BNZS /tmp/stageD_run.mac
```

如果要换 placement，只改：

```text
/cfg/setPlacementFile ../Input/output_pbc/<ratio>/<placement>.csv
```

如果要换 ratio，也同步改：

```text
/cfg/setWeightRatio <bnWt> <znsWt>
```

## 7. 批处理脚本用法

如果你不想手写宏，直接用：

```bash
python3 stageD_run_batch.py --ratio-tag 1-2
```

默认行为：

- 从 `Input/output_pbc/1-2/` 读取全部主 placement，并排除辅助 CSV
- 只生成宏
- 不实际执行

最常用参数：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --beam-on 10000 \
  --run
```

只跑前 5 个 placement：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --max-placements 5 \
  --beam-on 10000 \
  --run
```

跳过前 10 个，再跑后面的 5 个：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --start-index 10 \
  --max-placements 5 \
  --beam-on 10000 \
  --run
```

只跑指定 placement：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --placements placement_f_0.64_0004.csv placement_f_0.64_0008.csv \
  --beam-on 10000 \
  --run
```

打乱顺序后再选：

```bash
python3 stageD_run_batch.py \
  --ratio-tag 1-2 \
  --shuffle \
  --seed 20260507 \
  --max-placements 5 \
  --beam-on 10000 \
  --run
```

## 8. 输出位置

单个 run 输出：

```text
Output/stageD_optical_homogenization/<ratio>/<placement_stem>/<wavelength>/
  <source>__<boundary>__<config_hash>/run_<run_id>_seeds_<seed0>_<seed1>/
```

其中主要文件：

```text
optical_homogenization_events.csv
optical_homogenization_summary.csv
```

批处理脚本还会生成：

```text
Output/stageD_macros/<ratio>/
logs/stageD/<ratio>/
```

## 9. summary 里看什么

最关键列：

- `n_photons`
- `n_absorbed`
- `absorbed_fraction`
- `n_lost`
- `lost_fraction`
- `total_path_length_um`
- `total_real_scatter`
- `total_reentry`
- `mean_num_reentry`
- `mu_a_raw_per_um`
- `mu_s_raw_per_um`
- `g_raw`
- `mu_s_prime_raw_per_um`

当前推荐先检查：

- `n_lost` 是否接近 0
- `mu_s_raw_per_um` 是否非 0
- `g_raw` 是否非 0
- `mean_num_reentry` 是否处在合理范围

## 10. events 里看什么

如果要排查异常，重点列：

- `final_status`
- `num_reentry`
- `num_material_boundary`
- `num_real_scatter`
- `total_path_length_um`

如果大量出现这些状态，要小心：

- `lost`
- `reentry_failed`
- `max_reentry`
- `max_steps`
- `max_path_length`

## 11. summary 表和 merge 表的区别

单个 placement 的：

```text
Output/stageD_optical_homogenization/<ratio>/<placement>/<wavelength>/
  <config>/<run>/optical_homogenization_summary.csv
```

它给的是这个 placement 的 raw 参数：

- `mu_a_raw_per_um`
- `mu_s_raw_per_um`
- `g_raw`
- `mu_s_prime_raw_per_um`

它本身不直接给多 placement 误差条。

多个 placement merge 后的：

```text
Output/optical_params/<ratio>/rve_raw_optical_params_by_ratio.csv
```

这个表会给：

- `mu_a_mean_per_um`
- `mu_a_std_per_um`
- `mu_s_mean_per_um`
- `mu_s_std_per_um`
- `g_mean`
- `g_std`
- `mu_s_prime_mean_per_um`
- `mu_s_prime_std_per_um`

这里的 `std` 是 placement-to-placement 的离散度。

## 12. 发射角度和 re-entry 对 g 的影响

当前发射模型：

- 发射位置：在随机选中的 `ZnS` 球内均匀采样
- 发射方向：各向同性随机
- 偏振：随机，且垂直于动量方向

因此当前 source 没有额外方向偏置。

关于 `g`：

- re-entry 本身不直接改 `g`
- re-entry 不改方向
- re-entry 不记为 scatter
- `g` 只由后续真实方向变化统计得到

但 re-entry 会改变 photon 接下来遇到的局部微结构，因此会间接影响后续传播路径采样。

## 13. 当前验证情况

当前已经验证：

- `beamOn 10000`
- ratio `1-2`
- placement:
  - `placement_f_0.64_0001.csv`
  - `placement_f_0.64_0002.csv`
  - `placement_f_0.64_0004.csv`

这些 case 目前都能跑完，并且：

- `n_lost = 0`
- `mu_s_raw > 0`
- `g_raw > 0`

示例结果量级大致是：

- `mu_a_raw_per_um ~ 0.0181`
- `mu_s_raw_per_um ~ 0.67`
- `g_raw ~ 0.52`
- `mu_s_prime_raw_per_um ~ 0.32`
- `mean_num_reentry ~ 0.226`

这些值目前应理解为：

- simulation-derived raw parameters
- not yet experimentally calibrated

## 14. 后处理脚本

合并同一 ratio 下多个 placement：

```bash
python3 scripts/merge_stageD_optical_params.py \
  --ratio 1-2 \
  --project-root .
```

输出：

```text
Output/optical_params/1-2/rve_raw_optical_params_by_ratio.csv
```

如果后续有实验 `L_eff` 数据，再做校正：

```bash
python3 scripts/calibrate_optical_params_with_Leff.py \
  --ratio 1-2 \
  --project-root . \
  --experimental-leff /path/to/experimental_Leff.csv
```

## 15. 当前已知限制

- `n_eff_initial` 还没有在主程序里自动输出
- 没有实验时，只能得到 `raw/provisional` 参数
- `mu_s/g` 仍依赖当前几何、材料光学常数和统计定义
- 目前散射统计主要来自内部相界面方向改变
- 不应把当前结果直接当作最终真实样品参数

## 16. 推荐下一步

当前最推荐的工作流：

1. 对同一 ratio 跑多个 placement
2. 合并得到 `raw optical parameters`
3. 比较 placement 间离散度
4. 如有需要，扫描不同波长
5. 后续若拿到实验 `L_eff`，再做 calibrated parameters

如果只是当前阶段使用：

- 可以用
- 但请把输出理解为 `raw/provisional homogenized optical parameters`
