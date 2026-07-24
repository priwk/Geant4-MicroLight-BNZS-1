#ifndef AnalysisConfig_h
#define AnalysisConfig_h 1

#include <filesystem>
#include <string>

enum class RunMode
{
  StageA_NeutronPatch,  // 固定 50x50x30 um^3 微结构 patch 做热中子等效化
  StageB_ReplayAlphaLi, // 读取 capture CSV，重放 alpha / Li7
  StageD_OpticalHomogenization // 随机 RVE 统计均匀化光学输运参数
};

class AnalysisConfig
{
public:
  AnalysisConfig();
  ~AnalysisConfig();

  static const char *RunModeName(RunMode mode);
  static std::filesystem::path ProjectRootPath();
  static std::string PathForRecord(const std::filesystem::path &path);
  static std::string PathForRecord(const std::string &path);

public:
  // ---- 全局运行模式 ----
  RunMode runMode;

  // ---- 固定微结构基准参数（当前项目中默认不应被随意改动）----
  double patchXY_um;        // 固定 patch 横向尺寸，默认 50 um
  double microThickness_um; // 固定局部厚度，默认 30 um

  // ---- 配比信息 ----
  double bnWt;
  double znsWt;

  // ---- placement 输入控制 ----
  bool useRandomPlacement;
  std::string placementFilePath;
  std::string placementGeometryMode;
  std::string periodicImagesFilePath;

  // ---- Stage B 输入 ----
  std::string captureCsvPath;
  std::string captureInputDir;

  // ---- 共用光学材料参数 ----
  bool opticalParamsProvided;
  double opticalMatrixRIndex;
  double opticalMatrixAbsLengthUm;
  double opticalBnRIndex;
  double opticalBnAbsLengthUm;
  double opticalZnsRIndex;
  double opticalZnsAbsLengthUm;

  // ---- Stage D optical homogenization ----
  double stageD_wavelength_nm;
  std::string stageD_source_mode;
  std::string stageD_boundary_mode;
  std::string stageD_reentry_mode;
  std::string stageD_particle_reentry_mode;
  std::string stageD_matrix_reentry_mode;
  std::string stageD_scatter_metric;
  int stageD_target_primary_scatter;
  double stageD_theta_threshold_deg;
  int stageD_max_reentry;
  int stageD_max_steps;
  double stageD_max_path_length_um;
  std::string stageD_output_dir;
  int stageD_portal_nu;
  int stageD_portal_nv;
  double stageD_portal_margin_um;
  double stageD_clearance_bin0_um;
  double stageD_clearance_bin1_um;
  double stageD_clearance_bin2_um;
  int stageD_max_particle_reentry_trials;
  int stageD_max_portal_fallback_level;

  // ---- Stage B 深度映射兼容开关 ----
  // true: 允许 thickness_um == local patch thickness
  // 后面会用于替换当前 PrimaryGeneratorAction 里过严的限制
  bool allowThicknessEqualLocalPatch;
};

#endif
