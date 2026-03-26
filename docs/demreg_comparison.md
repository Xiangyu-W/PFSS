# DEM Code Examples 对比分析

对比两个 DEM（Differential Emission Measure）代码示例：

1. **demregpy 示例** — [`alasdairwilson/demregpy`](https://github.com/alasdairwilson/demregpy/blob/main/examples/running_demregpy_on_aia_data.py)
2. **demreg 示例** — [`ianan/demreg`](https://github.com/ianan/demreg/blob/master/python/example_demregpy_aiapxl.ipynb)

---

## 1. 基本信息

| 项目 | demregpy 示例 | demreg 示例 |
|------|--------------|-------------|
| **文件类型** | Python 脚本 (`.py`) | Jupyter Notebook (`.ipynb`) |
| **所属仓库** | `alasdairwilson/demregpy` | `ianan/demreg` |
| **核心函数** | `dn2dem`（from `demregpy`） | `dn2dem_pos`（from `dn2dem_pos`） |
| **分析对象** | 整个子区域（200×200 pixel region） | 单个像素 |
| **数据来源** | 在线下载 AIA 数据（Fido） | 预存本地子图（submaps） |

---

## 2. 使用的包与导入

### demregpy 示例
```python
from demregpy import dn2dem
from demregpy.tresp import aia_tresp  # 内置温度响应函数路径
from aiapy.calibrate import correct_degradation
from aiapy.calibrate.utils import get_correction_table
```

### demreg 示例
```python
from sys import path as sys_path
sys_path.append('/Users/iain/github/demreg/python')  # 手动添加路径
from dn2dem_pos import dn2dem_pos  # 原始 demreg 函数

from aiapy.calibrate import degradation
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.util import get_pointing_table, get_correction_table, get_error_table
```

**区别**：
- demregpy 使用 `pip install` 后直接 `import`，温度响应函数内置在包中
- demreg 需要手动设置 `sys.path`，函数名为 `dn2dem_pos`（而非 `dn2dem`）
- demreg 使用 `aiapy.calibrate.util`（旧版 API），demregpy 使用 `aiapy.calibrate.utils`（新版）

---

## 3. 温度响应函数（Temperature Response）

### demregpy 示例
```python
trin = io.readsav(aia_tresp)  # aia_tresp 是 demregpy.tresp 中的内置路径
```
- 使用 demregpy 包内置的 SSWIDL 生成的响应函数文件
- 文件由 `make_aiaresp_forpy.pro` 生成

### demreg 示例
```python
trin = io.readsav('aia_tresp_en.dat')
```
- 使用本地的 `aia_tresp_en.dat` 文件
- 同样由 SSWIDL 生成

**区别**：两者都依赖 SSWIDL 生成的响应函数，但 demregpy 将其打包到了库内部，更方便使用。

---

## 4. 数据获取与预处理

### demregpy 示例
```python
# 在线下载数据
q = Fido.search(
    attrs.Time(time_test, time_test + td),
    attrs.Instrument('AIA'),
    attrs.Wavelength(channels[0]) | ... | attrs.Wavelength(channels[5]),
)
files = Fido.fetch(q)

# 加载并排序
maps = [Map(f) for f in files]
maps = sorted(maps, key=lambda x: x.wavelength)

# 退化校正 + 曝光时间归一化（一步完成）
maps = [correct_degradation(m, correction_table=get_correction_table("JSOC")) / m.exposure_time 
        for m in maps]
```

### demreg 示例
```python
# 从本地加载预存的 submap FITS 文件
ffin = sorted(glob.glob(fdir + 'aia_smd_*.fits'))
aprep = sunpy.map.Map(ffin)

# 退化校正因子（独立计算）
cor_tab = get_correction_table("ssw")
degs[i] = degradation(channels[i], time, correction_table=cor_tab)

# DN 数据除以退化因子和曝光时间
dn_in = cor_data / durs  # cor_data = data / degs
```

**关键区别**：
| 方面 | demregpy 示例 | demreg 示例 |
|------|--------------|-------------|
| **数据下载** | 实时从 VSO/JSOC 下载 | 使用预存本地文件 |
| **AIA prep** | 使用 `correct_degradation` 一步完成 | 分步：`update_pointing` → `register` → 手动计算退化因子 |
| **退化校正** | `correct_degradation()` 返回校正后的 Map | 手动用 `degradation()` 获取退化因子，手动 `data/degs`（甚至硬编码数值避免重复请求） |
| **曝光时间** | `/ m.exposure_time` 归一化 | `/ durs` 手动归一化（从 `m.meta['exptime']` 取） |
| **Pointing 校正** | 未显式做（注释掉了 `aiaprep`） | 显式使用 `update_pointing` + `register`（代码被注释，使用预存 submap） |
| **PSF 反卷积** | 未做 | 未做 |
| **数据搜索方式** | VSO/Instrument 搜索：`attrs.Instrument('AIA')` | VSO 搜索（注释掉），主要使用预存本地 submap |

---

## 5. 误差估计

### demregpy 示例
```python
# 标准 AIA 误差模型
serr_per = 10.0  # 系统误差百分比
npix = 4096.**2 / (nx * ny)
gains = np.array([18.3, 17.6, 17.7, 18.3, 18.3, 17.6])
dn2ph = gains * [94, 131, 171, 193, 211, 335] / 3397.0
rdnse = 1.15 * np.sqrt(npix) / npix
drknse = 0.17
qntnse = 0.288819 * np.sqrt(npix) / npix

# 组合误差
etemp = np.sqrt(rdnse**2 + drknse**2 + qntnse**2 + (dn2ph[j]*abs(data[:,:,j]))/(npix*dn2ph[j]**2))
esys = serr_per * data[:,:,j] / 100.
edata[:,:,j] = np.sqrt(etemp**2 + esys**2)
```

### demreg 示例
```python
# 类似的误差模型，但只有光子噪声 + 读出噪声
gains = np.array([18.3, 17.6, 17.7, 18.3, 18.3, 17.6])
dn2ph = gains * np.array([94, 131, 171, 193, 211, 335]) / 3397.
rdnse = np.array([1.14, 1.18, 1.15, 1.20, 1.20, 1.18])  # 每个波段独立的读出噪声
num_pix = 1

# 光子噪声（注意：使用未退化校正的原始 DN 数据计算光子数）
shotnoise = (dn2ph * data * num_pix)**0.5 / dn2ph / num_pix / degs

# 组合误差（DN/px/s）
edn_in = (rdnse**2 + shotnoise**2)**0.5 / durs
```

**关键区别**：
- demregpy 示例包含 **暗电流噪声 (`drknse = 0.17`)**、**量化噪声 (`qntnse = 0.289`)**、**光子噪声** 和 **10% 系统误差 (`serr_per`)**
- demreg 示例只包含 **读出噪声 (`rdnse`)** 和 **光子噪声 (`shotnoise`)**，没有暗电流、量化噪声和系统误差
- demreg 示例中每个波段有独立的 `rdnse` 值（`[1.14, 1.18, 1.15, 1.20, 1.20, 1.18]`），而 demregpy 示例使用统一的 `1.15`
- demregpy 示例中 `npix = 4096² / (nx * ny)`，考虑了像素 rebinning 的情况
- demreg 示例提到可以额外加入 ~20% 的系统误差，但留给读者自行决定（原文："left as an exercise for the reader"）
- demreg 示例还讨论了使用 `aiapy.calibrate.estimate_error`（需要 aiapy >0.6）作为替代方案，并验证了手动计算与 aiapy 的结果一致

---

## 6. DEM 反演方法与参数

### demregpy 示例
```python
# 温度 bin 设置
nt = 16
temperatures = 10**np.linspace(5.7, 7.1, num=nt+1)

# DEM 权重（基于高斯 DEM 模型）
demwght0 = 10**np.interp(mlogt, tresp_logt, np.log10(dem_mod))
demwght0 /= max(demwght0)
dem_norm = np.array([[demwght0 for i in range(200)] for j in range(200)])

# 反演（单次调用，整个区域）
dem, edem, elogt, chisq, dn_reg = dn2dem(
    data[x1:x2, y1:y2, :], 
    edata[x1:x2, y1:y2, :],
    trmatrix, tresp_logt, temperatures, 
    dem_norm0=dem_norm, emd_int=True
)
```

### demreg 示例
```python
# 温度 bin 设置（更宽的范围）
temps = np.logspace(5.7, 7.6, num=42)

# 反演（三种模式，单个像素）
# 1. 默认模式（自动权重）
dem0, edem0, elogt0, chisq0, dn_reg0 = dn2dem(dn_in, edn_in, trmatrix, tresp_logt, temps)

# 2. EM Loci 权重模式
dem1, edem1, elogt1, chisq1, dn_reg1 = dn2dem(dn_in, edn_in, trmatrix, tresp_logt, temps, gloci=1)

# 3. EM Loci + EMD 内部计算模式
dem2, edem2, elogt2, chisq2, dn_reg2 = dn2dem(dn_in, edn_in, trmatrix, tresp_logt, temps, 
                                                 gloci=1, emd_int=True)
```

**关键区别**：
| 方面 | demregpy 示例 | demreg 示例 |
|------|--------------|-------------|
| **温度范围** | log T = 5.7 – 7.1 | log T = 5.7 – 7.6（更宽） |
| **温度 bin 数** | 16+1 = 17 个边界 | 42 个边界 |
| **空间维度** | 200×200 pixel 区域（3D） | 单个像素（1D） |
| **权重 (`dem_norm0`)** | 基于高斯 DEM 模型预计算 | 默认模式自动计算（两次正则化） |
| **`gloci` 参数** | 未使用（默认 0） | 测试了 `gloci=0` 和 `gloci=1` |
| **`emd_int` 参数** | `True` | 分别测试了 `False` 和 `True` |
| **反演次数** | 1 次 | 3 次（比较不同方法） |

---

## 7. 输出与可视化

### demregpy 示例
- 16 个 DEM 温度 bin 的 2D 图（200×200 区域，`imshow` + `inferno` colormap）
- 两个选定像素的 DEM 曲线（带误差棒）
- 171Å submap 图（`sqrt` 拉伸）
- 全盘 AIA 图

### demreg 示例
- 三种反演方法的 DEM 曲线分别绘制（带误差棒 + χ² 标注）
- DN_in vs DN_reg 散点图（log-log 坐标，含 6 个波段标注，三种方法叠加）
- DN_reg / DN_in 比值图（检查每个波段的重建质量）
- `(DN_in - DN_reg) / σ_DN_in` 残差图（归一化残差，判断拟合偏差）
- 综合 summary 图：DEM 曲线 + 残差并排

**区别**：demreg 示例更侧重于 **方法比较** 和 **反演质量评估**（χ²、残差分析），而 demregpy 示例更侧重于 **空间分布** 的可视化。

---

## 8. 核心函数接口差异

### `dn2dem`（demregpy）
```python
dem, edem, elogt, chisq, dn_reg = dn2dem(
    data,           # DN/px/s 数据，可以是多维数组
    edata,          # DN 误差
    trmatrix,       # 温度响应矩阵
    tresp_logt,     # 响应函数的 log T 网格
    temperatures,   # DEM 温度 bin 边界
    dem_norm0=...,  # 可选：初始权重
    emd_int=True    # 可选：在 EMD 空间内部计算
)
```

### `dn2dem_pos`（demreg 原版）
```python
dem, edem, elogt, chisq, dn_reg = dn2dem_pos(
    dn_in,          # DN/px/s 数据
    edn_in,         # DN 误差
    trmatrix,       # 温度响应矩阵
    tresp_logt,     # 响应函数的 log T 网格
    temps,          # DEM 温度 bin 边界
    gloci=0,        # 权重方案：0=自动（默认），1=EM Loci
    emd_int=False   # 是否在 EMD 空间计算
)
```

> ⚠️ **注意**：demreg 示例的 notebook 标题中写的是 `demregpy`，但实际使用的是 `dn2dem_pos`（原始 demreg 包的函数）。从示例代码看，`dn2dem`（demregpy版）和 `dn2dem_pos`（原版）的接口非常相似，返回值也相同。

---

## 9. 总结

| 维度 | demregpy 示例 | demreg 示例 |
|------|--------------|-------------|
| **适合场景** | 快速上手、区域 DEM 分析 | 深入理解方法、单像素分析 |
| **安装便捷性** | ✅ `pip install demregpy` | ❌ 需手动添加路径 |
| **数据预处理** | 较简单（`correct_degradation` 一步） | 更详细（分步 prep + 手动退化校正） |
| **误差模型** | 更完整（含暗电流、量化、系统误差） | 较简化（仅读出+光子噪声） |
| **反演方法** | 单一方法，但支持自定义权重 | 三种方法对比（默认/EM Loci/EMD） |
| **质量评估** | 仅 χ² | χ²、残差、DN 重建比较 |
| **文档/注释** | 较少 | 非常详细，附历史版本记录 |
| **代码风格** | 精简的脚本 | 交互式探索性分析 |
| **推荐用途** | 生产级 DEM 分析 | 学习 DEM 方法、对比不同参数 |

---

## 10. 实践建议

1. **如果你是初学者**：从 demreg 示例（notebook）开始，因为它提供了三种方法的对比，以及更详细的注释和质量评估。
2. **如果你要做批量分析**：使用 demregpy 包，因为它支持多维数组输入（整个区域一次反演）并且安装更方便。
3. **温度范围选择**：demregpy 用 5.7–7.1，demreg 用 5.7–7.6。选择取决于观测目标——如果关注热等离子体（如耀斑），使用更宽的范围。
4. **误差估计**：建议参考 demregpy 示例的更完整误差模型，或使用 `aiapy.calibrate.estimate_error`（需 aiapy ≥ 0.6）。
5. **权重选择**：
   - 默认自动权重（`gloci=0`）：最安全的选择，适合 AIA 数据
   - EM Loci 权重（`gloci=1`）：适合有尖锐温度响应的滤光片（如 X 射线观测）

---

## To do
1. 使用 aiapy 的误差估计
2. 使用 demreg_pos 的自动权重方法计算 2D DEM
3. DEM --> Temperature