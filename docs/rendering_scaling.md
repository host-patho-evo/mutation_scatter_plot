# Graphical Scaling and Rendering Architecture

This document describes the mathematical relationship between mutation frequencies
and their visual representation in both the Matplotlib and Bokeh backends.

## Three Scaling Modes

`mutation_scatter_plot` supports three circle-size scaling modes selected via
command-line flags.  The active mode is embedded in every output filename between
the matrix name and the colormap name:

```
{prefix}.{MATRIX}.{area_scaling|linear_scaling}.{colormap}.{ext}
```

This ensures outputs from different modes never silently overwrite each other.

---

### 1. Default — Area-proportional (`area_scaling`)

*No flag required.*  This is the perceptually recommended mode for bubble charts
(see Cleveland & McGill, 1984).

| Backend | Parameter | Formula | Proportionality |
|:---|:---|:---|:---|
| **Matplotlib** | `s` (area, pt²) | `s = freq × 5000` | area ∝ freq |
| **Bokeh** | `size` (diameter, px) | `size = √freq × 100` | area ∝ freq |

**Result in both backends:** area ∝ frequency, perceived radius ∝ √frequency.
A mutation twice as common occupies twice as much "ink" on the page.

> [!NOTE]
> Bokeh's raw `scatter(size=...)` parameter is a **diameter** in screen pixels,
> so `size² ∝ area`.  Without the square-root transform, a naive `size = freq × 100`
> would make area ∝ freq², causing high-frequency mutations to appear
> disproportionately large.  The default sqrt transform corrects this.

---

### 2. Linear-radius mode (`linear_scaling`)

Flag: **`--linear-circle-size`**

Renders circles with **radius proportional to frequency** in both backends.
A mutation twice as common appears exactly twice as wide.

| Backend | Parameter | Formula | Proportionality |
|:---|:---|:---|:---|
| **Matplotlib** | `s` (area, pt²) | `s = freq² × 5000` | radius ∝ freq |
| **Bokeh** | `size` (diameter, px) | `size = freq × 100` | radius ∝ freq |

**Derivation (Matplotlib):**
To achieve radius ∝ freq, we need area = π r² = π (freq × C)².
Choosing C = √(5000/π) ≈ 39.9 pt gives area = freq² × 5000 at the same maximum
circle size as the default mode (freq = 1 → s = 5000 in both formulas).

This mode **implies `--disable-bokeh-sqrt-size`** internally — `main()` sets
`bokeh_sqrt_size = False` whenever `linear_circle_size = True`, keeping both
figure types visually in sync without requiring an extra flag.

This matches the v0.2 Bokeh rendering (which used `size = freq × 5000`) while
correcting the Matplotlib side, which was always area-proportional in that version.

> [!TIP]
> Use `--linear-circle-size` when you want low-frequency mutations to appear
> *relatively* larger than in the default mode, or when you need to match the
> visual appearance of v0.2 figures.

---

### 3. Legacy Bokeh-only linear mode

Flag: **`--disable-bokeh-sqrt-size`**

Reverts Bokeh to diameter-proportional scaling without changing Matplotlib.

| Backend | Parameter | Formula | Proportionality |
|:---|:---|:---|:---|
| **Matplotlib** | `s` (area, pt²) | `s = freq × 5000` | area ∝ freq (unchanged) |
| **Bokeh** | `size` (diameter, px) | `size = freq × 100` | radius ∝ freq |

> [!WARNING]
> This leaves the two backends perceptually **inconsistent** with each other.
> Prefer `--linear-circle-size` if you want linear radius scaling in both outputs.

---

## Comparison Table

| Mode | Flag | `_mpl_s` formula | `_bokeh_size` formula | Matplotlib | Bokeh |
|:---|:---|:---|:---|:---|:---|
| **area_scaling** (default) | *(none)* | `freq × 5000` | `√freq × 100` | area ∝ freq | area ∝ freq |
| **linear_scaling** | `--linear-circle-size` | `freq² × 5000` | `freq × 100` | radius ∝ freq | radius ∝ freq |
| Legacy | `--disable-bokeh-sqrt-size` | `freq × 5000` | `freq × 100` | area ∝ freq | radius ∝ freq |

---

## Verification Data (Default Area-proportional Mode)

The following table shows internal scaling values for the standard legend
frequencies (0.1 %, 1 %, 10 %, 30 %).

**Input:** `tests/outputs/test2_full.x_after_count.frequencies.tsv`

| Frequency | Matplotlib `s` (pt²) | Area ratio | Bokeh `size` (px) | Calculated area (πr²) | Area ratio |
|:---|:---|:---|:---|:---|:---|
| 0.001 (0.1 %) | 5.0 | **1.0 ×** | 3.16 | 7.85 | **1.0 ×** |
| 0.01  (1.0 %) | 50.0 | **10.0 ×** | 10.00 | 78.54 | **10.0 ×** |
| 0.1  (10.0 %) | 500.0 | **100.0 ×** | 31.62 | 785.40 | **100.0 ×** |
| 0.3  (30.0 %) | 1500.0 | **300.0 ×** | 54.77 | 2356.19 | **300.0 ×** |

## Rasterized Verification

A quantitative pixel-area audit at 600 DPI (9 600 × 5 400 px PNG, OpenCV
contour analysis) confirmed the theoretical 1 : 10 : 100 : 300 area ratios
with > 95 % precision.  Small deviations are due to sub-pixel aliasing and
circular rasterization discretization.

## Coordinate Space Alignment

The frequency legend (`ax2`) and scatter plot (`ax1`) share a consistent
coordinate space.  Legend points are projected to Y = −400 (outside the
visible data range `[0, 1]`) so they can be captured by Matplotlib's legend
engine without interfering with the data visualization.
