# Symmetrized Colormap Mapping Architecture

This document describes the architectural implementation for zero-centered symmetric colormap bounds scaling natively applied in `mutation_scatter_plot`.

Per design requirements, zero *must* remain the perfect equilibrium (neutral grey `#dcdddd` or target center hue) across divergent colormaps regardless of empirical boundaries. This prevents red/blue tones from sliding over zero when asymmetric score minimum/maximums exist in the evaluated dataset.

## Implementation Mechanics

### 1. Symmetrization of Empirical Bounds
The pipeline natively calculates internal minimum and maximum boundaries iteratively by scoring exclusively genuine substitution ranges found inside the analyzed dataset (falling back symmetrically `[-bound_abs, +bound_abs]`).
```python
_bound_abs = max(abs(min(_emp_scores)), abs(max(_emp_scores)))
myoptions.cmap_actual_vmin = -_bound_abs
myoptions.cmap_actual_vmax = _bound_abs
```

### 2. Symmetrization of Theoretical Virtual Bounds (`load_matrix`)
For globally unified dataset tones natively enforced via `--spread-colormap-over-full-but-virtual-matrix-range`, the theoretical substitution matrix (e.g. BLOSUM80) limits are symmetrically evaluated identically:
```python
_bound_abs = max(abs(min(_theoretical)), abs(max(_theoretical)))
myoptions.matrix_min_theoretical = -_bound_abs
myoptions.matrix_max_theoretical = _bound_abs
```

### 3. Colormap Bound Unification
By forcing boundaries to be strictly symmetrized, the original geometrical `len(colors) // 2` midpoint derivations and continuous `[-cb_vmin, cb_vmax]` normalizations logically perfectly target index center `0`, naturally keeping Bokeh and Matplotlib rendering cleanly in sync with `coolwarm_r` divergent properties identically!

---

## Side-by-Side Symmetrized RGB Table Comparison

Assuming a hypothetical dataset empirically limited to `[-4, +4]` bound range dynamically, running under BLOSUM80 (`[-11, +11]` symmetric theoretic virtual):

| Score | Continuous HEAD | Cont. Virtual ($\pm11$) | Cont. Dynamic ($\pm4$) | Discrete HEAD | Discrete Sliced Vir. ($\pm11$) | Discrete Sliced Dyn. ($\pm4$) |
|-------|-----------------|-------------------------|------------------------|---------------|--------------------------------|-------------------------------|
| **-6**| `#ff0000` (Red Sent.) | `#ff0000` (Red Sent.) | `#ff0000` (Red Sent.) | `#cc0000` | `#cc0000` | `#cc0000` |
| **-5**| `#f6a283` | `#f6a283` | `#b40426` (Clamp)| `#ff0000` | `#ff0000` | `#ff0000` |
| **-4**| `#f7b396` | `#f7b396` | `#b40426` | `#ff4f00` | `#ff4f00` | `#ff4f00` |
| **-3**| `#f5c1a9` | `#f5c1a9` | `#de614d` | `#ff7c7c` | `#ff7c7c` | `#ff7c7c` |
| **-2**| `#f1cdba` | `#f1cdba` | `#f49a7b` | `#ff9999` | `#ff9999` | `#ff9999` |
| **-1**| `#e8d6cc` | `#e8d6cc` | `#f4c5ad` | `#9c644b` (Brown) | `#9c644b` (Brown) | `#9c644b` (Brown) |
| **0** | `#dcdddd` *(Centre)* | `#dcdddd` *(Centre)* | `#dcdddd` *(Centre)*| `#ffff00` (Yellow) | `#ffff00` (Yellow) | `#ffff00` (Yellow) |
| **1** | `#d1dae9` | `#d1dae9` | `#b7cff9` | `#ffcc00` (Gold) | `#ffcc00` (Gold) | `#ffcc00` (Gold) |
| **2** | `#c3d5f4` | `#c3d5f4` | `#8caffe` | `#ffa200` | `#ffa200` | `#ffa200` |
| **3** | `#b5cdfa` | `#b5cdfa` | `#6180e9` | `#7DCCFF` | `#7DCCFF` | `#7DCCFF` |
| **4** | `#a5c3fe` | `#a5c3fe` | `#3b4cc0` | `#0042ff` | `#0042ff` | `#0042ff` |
| **...**| ... | ... | ... | ... | ... | ... |
| **11**| `#3b4cc0` | `#3b4cc0` | `#3b4cc0` (Clamp)| `#97ce2f` | `#97ce2f` | `#97ce2f` |
| **12**| `#219f11` (Green Sent.)| `#219f11` (Green Sent.)| `#219f11` (Green Sent.)| `#219f11` (Sent.)| `#219f11` (Sent.) | `#219f11` (Sent.) |

### ListedColormap (Discrete) Strategy
The discrete path natively utilizes a **"Cut Out" / Slicing** approach explicitly. The categorical `amino_acid_changes` map conceptually maintains its central `[-19, 19]` anchor internally identically. Natively, we extract cleanly the exact subset array span demanded by the dataset empirical bounds while 100% physically duplicating the integer gap skips of the original Matplotlib arrays. This seamlessly renders a physically shorter colorbar legend reflecting perfectly the identical hex colors without mathematically rescaling/interpolating categorical color groupings inherently!
