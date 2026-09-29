# SI 线性不稳定性 MATLAB 求解器

## 1. 文件和运行

实际求解程序只有两个文件：

- `main_SI_instability.m`：主程序。可修改的无量纲参数、波数范围、数值精度及输出选项集中在前 60 行左右。
- `solve_SI_instability.m`：函数文件。负责离散矩阵、边界条件、特征值求解、波数扫描和局地峰值优化。其余辅助函数均包含在同一个文件内。

将两个 `.m` 文件放在同一文件夹，在 MATLAB 中运行：

```matlab
main_SI_instability
```

不需要任何预先生成的 `.mat`、`.csv`，也不读取经验公式来生成数值解。压缩包中的参考数据和检查记录只供核对，不是计算输入。

程序使用基本 MATLAB 函数（如 qz、eig、chol、fminbnd、table），不调用第三方工具箱。建议采用支持隐式数组扩展的 MATLAB 版本（R2016b 或更新）。这份 `.m` 程序尚未在 MATLAB 本体中执行；本次采用 Python/SciPy 对相同矩阵的移植实现进行了验证。请先运行主程序自带的四个固定波数基准检查，再进行大规模参数扫描。

## 2. 模型及变量约定

模型对应正文的固定热成风背景及其线性扰动方程，不是随混合系数变化而重新调整的 TTW 平衡背景。

- z ∈ [-1,0]；x 为沿锋方向，y 为跨锋方向。
- U(z)=z+0.5，V=W=0，背景 B_z=1、B_y=-1/Ri。
- 正常模为 exp[σt+i(kx+ℓy)]；Q=σ+ikU，D=d/dz。
- Ri=N_b²/S_0²；δ=H/L=f/S_0；E_v=ν_v/(fH²)；Pr=ν_v/κ_v；E_b=E_v/Pr。
- σ 及 Re(σ) 已经无量纲化。实际增长率为 f Re(σ)，不能再把输出除以 f。
- k、ℓ、K 都是以 L 为水平尺度的无量纲波数，K²=k²+ℓ²。
- 波向定义 k=K sinφ、ℓ=K cosφ。strict SI 对应 k=0、φ=0。
- 水平波长 λ/L=2π/K；strict SI 下 K=ℓ>0。

函数头部列出了所求解的五条原始扰动方程。

默认扩散仅含垂向黏性和浮力扩散。有限混合的扰动边界条件是 w=0、Du=Dv=Db=0。边界条件施加于扰动，而不是要求维持的背景剪切也为零。

## 3. 最常修改的主程序设置

```matlab
Ri_values = 0.70;
Ev_values = [1e-4, 1e-3, 3e-3];
p.delta = 0.10;
p.Pr = 1.00;
p.buoyancy_bc = 'neumann';
p.horizontal_ratio = 0;

opt.task = 'si';
opt.method = 'chebyshev';
opt.N = 96;
opt.scan_N = 64;
opt.K_grid = logspace(-1,log10(350),41);
opt.growth_tol = 1e-6;
```

`opt.scan_N=64` 仅用于寻找峰值所在的粗网格区间；局地优化和最终特征值/本征函数使用 `opt.N=96`。设置 `opt.scan_N=[]` 可在所有步骤使用同一个 N。

当前 Figure 3 的完整参数组合已在主程序中注释列出。取消下面两行的注释、替代默认向量，即可重新求解全部 88 个组合：

```matlab
Ri_values = [0.20,0.25,0.30,0.40,0.50,0.60,0.70,0.80];
Ev_values = [1e-5,3e-5,1e-4,3e-4,1e-3,1.5e-3,2e-3,3e-3,4e-3,6e-3,1e-2];
```

每次调用函数只处理一个 Ri、E_v 组合；主程序自动循环。delta 和 Pr 是当前整个扫描的标量设置，可直接修改后重跑。改变这些参数不会重新计算背景的非地转响应。

### 可选敏感性设置

- `p.buoyancy_bc='dirichlet'`：改为 b=0 的浮力边界条件；其他边界条件不变。
- `p.horizontal_ratio=1`：添加与垂向系数相等的**有量纲**水平黏性和扩散系数。对应无量纲动量的附加阻尼为 E_v δ²K²，浮力为 E_b δ²K²。默认 0 不包含这些项。
- `p.U_offset`：只改变 U(z)=z+U_offset 的常数偏移，不改变 U_z=1。k=0 增长率不受影响；k≠0 应只产生相位频率的 Doppler 偏移。

## 4. 三种任务

### (a) `opt.task='si'`：strict-SI 扫描

严格固定 k=0，在 `opt.K_grid` 中扫描 ℓ=K>0。每个波数均求全部通过有限性判别的特征值，保留最大实部。随后在正增长局地峰值两侧的数值区间内，用 fminbnd 在 log K 上优化。经验式或理论波长均不用于设置峰值位置。

### (b) `opt.task='point'`：任意指定波数

```matlab
opt.task = 'point';
opt.k = 0;
opt.ell = 25.5513726457298;
```

函数返回该点的特征值谱、正增长特征值，以及领先模态的速度/浮力/压力结构。此任务不寻找峰值，因此 `R.peak.lambda_over_L` 不被当作“所选波长”赋值；指定波数的几何波长保存在 `R.peak.candidate_lambda_over_L`。

### (c) `opt.task='spectrum'`：包含偏轴模态的扫描

```matlab
opt.task = 'spectrum';
opt.phi_deg = -90:5:90;
% 或专门检查近轴情况：opt.phi_deg = -10:0.5:10;
```

在每个输入波向上扫描/优化 K，比较所有被采样的角度。**这是对离散角度的比较，不是连续二维全局优化。** 要判断峰值是否略偏离 k=0，需要加密 phi_deg 并检查 N。完整波谱结果存放在 `R.scan`。正负增长都保留，不把所有非正增长改写成零。

## 5. 两种离散方法

### 默认：Chebyshev

`opt.method='chebyshev'; opt.N=96;`

移植现有 Figure 3 计算档案中 (ζ,w,β) 的广义本征矩阵。N 为区间数，节点数为 N+1。上下边界通过行替换实现，垂向速度质量矩阵不是单位阵。

为了保留原来计算的数值缩放，函数内部使用原档案的等价旋转坐标：

- k_c=δℓ，ℓ_c=-δk；M=-1/δ，S_2=Ri/δ²。
- β 为能量归一化浮力，b_man=β/√S_2。
- u_man=-δv_c，v_man=δu_c，w_man=w_c。

这些变换均由函数完成。用户只输入正文坐标下的 k、ℓ，不应自行再转换一次。

质量矩阵因边界约束而奇异，求解使用 QZ，而不是 inv(B)*A 或 B\A。保留满足 |b_qz|>qz_tol(|a_qz|+|b_qz|) 的齐次特征值对。qz_tol 默认 1e-9；排除数量也保存在结果中。该筛选不是对增长率或拟合指数的筛选。

### 兼容原结果：交错网格压力投影

```matlab
opt.method = 'staggered';
opt.N = 80;
opt.scan_N = [];
```

移植原来 Figure 1、2、4 使用的兼容原始方程求解器：u、v、b、p 在中心，w 在界面；压力梯度和散度互为负伴随；通过离散无散子空间投影消压，再用最小二乘恢复压力。求解的是普通特征值问题。

N=80 并不保证所有参数下的峰值波长都已收敛。尤其低 Ri、极弱混合的短尺度峰值，对网格精度较敏感。这也是当前 Figure 3 用 Chebyshev96 峰值的原因。

### E_v=0

若请求 Chebyshev 而 E_v=0，程序明确改用交错原始方程，并在 `R.method_used`、`R.notes` 中记录。这样避免把有限黏性/扩散的附加边界条件错误施加到降阶无黏问题。对于 inviscid strict SI，保留有限扫描的增长率，但不报告一个虚假的有限最快增长波长。

## 6. 收敛与接近零增长的处理

```matlab
opt.check_N = 128;        % 在同一个波数处复算
opt.check_peak = true;   % 在同一数值峰值区间重新优化峰位
```

固定波数增长率误差与重新优化后的波长误差分别输出，不能相互替代。此处峰值复查限定在已识别的局地区间，不等同于以更高 N 重扫整个谱。

如果谱方法给出的最大实部接近零，默认用 N=160 的交错求解器复查输入扫描点。当两种方法对是否超过 growth_tol 的判断不一致时，结果被标记为 unresolved，不自动把其中一种认定为正确。

峰值落在 K 扫描边界、局地优化未收敛或到达其边界、没有检测到正增长，以及两种方法分类不一致时，都不赋予 selected wavelength。候选波数、原始特征值和状态保留在结果中。

这些检查只针对输入范围；没有检测到正增长不代表数学上证明任意波数都稳定。

## 7. 输出

默认新建 `SI_solver_output`：

- `SI_results.mat`：所有参数、完整 `solutions{Ev_index,Ri_index}`、峰值矩阵、增长率与状态。
- `SI_peak_summary.csv`：每个参数组合的增长率、候选波数、所选波长、状态和残差。
- `SI_numerical_data.csv`：strict-SI 情形的紧凑数值表，可用于后续 Figure 3 绘图。
- `fixed_wavenumber_benchmarks.csv`：主程序在用户 MATLAB 中实际执行的基准检查。

主要结果字段：

| 字段 | 含义 |
|---|---|
| `R.peak.sigma` | 候选峰值处的复特征值 |
| `R.peak.growth` | Re(σ)，无量纲 |
| `R.peak.k / ell / K / phi_deg` | 候选峰值的波数和角度 |
| `R.peak.lambda_over_L` | 有正增长且峰位可接受时的所选波长；否则 NaN |
| `R.peak.selected_wavelength_defined` | 所选波长是否被赋值 |
| `R.eigenvalues` | 候选点通过有限性筛选的特征值，按实部降序排列 |
| `R.unstable_eigenvalues` | 其中实部超过 growth_tol 的特征值 |
| `R.mode` | 领先本征模的 u、v、w、b、p、ζ 和垂向坐标 |
| `R.scan` | K–φ 扫描的领先特征值与增长率；N 是扫描精度 |
| `R.check` | 固定波数和局地重优化的分辨率检查 |
| `R.confirmation` | 接近中性时的另一种离散方法检查 |

R.mode 只保存领先本征函数，不保存全部本征向量。交错网格的 z 按从底到顶排列，Chebyshev 节点从顶到底；w 的坐标使用 `z_w`，不要假定它与中心变量点数相同。

本征函数整体归一化为最大原始速度分量的幅度为 1，所有物理变量使用同一个复数缩放。`discrete_eigenvector` 保存相应矩阵坐标下的本征向量，`discrete_coordinates` 说明其未知量顺序。打开 `opt.return_matrices=true` 后可直接核对 A*x−σB*x。

能量诊断中的 shear_rate、frontal_rate、viscous_loss_rate、buoyancy_loss_rate 都按两倍总扰动能量归一化；后两项取正损失。其组合应接近 Re(σ)。m_eff2 是按扰动动能归一化的垂向梯度指标。它与实际垂向相位波数不是同一定义。

谱方法的 primitive_residual 是恢复压力后，在原始配点方程上的诊断残差（包括受行替换影响的边界附近点），没有用它来事后筛选符合期望的模态。交错投影的代数残差和谱方法的能量闭合也不是相同的误差指标；小残差不证明峰值波长已空间收敛。

## 8. 本次验证与限制

`validation/` 内均为 **Python/SciPy 矩阵移植检查，不是 MATLAB 执行记录**：

- `matrix_benchmark_checks.csv`：五个固定波数案例与现有 Figure 3 数据比较，包含本征值、原始方程、连续方程和能量诊断。
- `reoptimized_peak_checks.csv`：六个增长案例和一个无检测到增长案例的波数重扫描/优化比较。所测试增长案例的波长与既有存档最大差异约 0.2243%，增长率最大差异约 0.0170%；更强混合的平坦峰顶对求解器数值噪声和优化细节更敏感。
- `structure_checks.csv`：交错矩阵与原档案的一致性、无散基、共轭对称、均匀流 Doppler 偏移和 Pr=1 时水平扩散的统一衰减检查。
- `Figure3_reference_data.csv`：现有 Figure 3 数值表，原样保留，只作核对。

没有在本次中重新执行全部 88 个组合，也没有在 MATLAB 本体中执行这两个文件。主程序包含运行时基准检查。改变 Ri、Ev、delta、Pr 或边界条件后，应重新检查波数范围和分辨率，不能把本次有限案例的检查外推成对所有参数的保证。

## 9. 来源

物理与数值约定依据本次对话中的 manuscript、Supporting Information Text S1，以及原始 Figure3_expanded_Ri_theory_reproduction_package.zip 中的 source/mixing_core.py。改写的是程序组织与 MATLAB 实现，不是基本态、方程或物理假设。

MATLAB QZ 接口核对：MathWorks 官方 qz 文档，https://www.mathworks.com/help/matlab/ref/qz.html 。广义本征问题不能用奇异质量矩阵的逆替代，参见 https://www.mathworks.com/help/matlab/ref/eig.html 。
