# Fault-Length Distribution Fitter (K-S Method)

Automatic detection of the valid power-law fitting range in fault-length frequency distributions, implemented in MATLAB.

---

## Background

Fault-length distributions in nature follow a **power-law** (linear in log-log space):

$$\log_{10}(N) = -b \cdot \log_{10}(L) + c$$

where $N$ is frequency, $L$ is fault length, and $b$ is the scaling exponent (always positive). However, raw field data deviate from this ideal at both ends:

| End | Observation | Physical cause |
|-----|-------------|----------------|
| **Low** | Curve flattens (slope ≈ 0) | Under-sampling — small faults are below detection resolution |
| **High** | Curve steepens | Truncation — large faults are cut by the edges of the study area |

Manually picking the "valid" fitting range is subjective and irreproducible. This tool automates that decision using the **Kolmogorov-Smirnov (K-S) test**.

---

## Algorithm

For every candidate window `[i, j]` of data points (subject to a minimum length):

1. Fit a straight line by least squares
2. Standardise the residuals to zero mean and unit variance
3. Apply a one-sample K-S test against N(0, 1)
   - A **high p-value** means residuals are consistent with Gaussian scatter → the window genuinely follows a straight line
4. Record the window index, p-value, K-S statistic, and length

**Selection rule:**
- Among all windows where `p-value ≥ α` (null not rejected), select the **longest** one to maximise the data used
- If no window passes at the chosen `α`, fall back to the window with the highest p-value and issue a warning

---

## Files

| File | Description |
|------|-------------|
| `fit_fault_ks.m` | Core function — automatic range detection and line fitting |
| `test_fault_ks.m` | Test script — generates a synthetic dataset and produces diagnostic plots |

---

## Usage

```matlab
[b, c, x_fit, y_fit, idx_low, idx_high] = fit_fault_ks(x, y)
[b, c, x_fit, y_fit, idx_low, idx_high] = fit_fault_ks(x, y, min_frac, alpha)
```

### Inputs

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `x` | vector | — | log₁₀(fault length), sorted ascending |
| `y` | vector | — | log₁₀(frequency) |
| `min_frac` | scalar | `0.25` | Minimum fraction of total points a window must contain. Prevents trivially short windows from passing the K-S test with insufficient data. |
| `alpha` | scalar | `0.05` | K-S test significance level. Windows with `p-value ≥ alpha` are considered linear. |

### Outputs

| Variable | Description |
|----------|-------------|
| `b` | Fitted slope (positive); model is `y = -b·x + c` |
| `c` | Fitted intercept |
| `x_fit` | x values in the accepted window |
| `y_fit` | y values in the accepted window |
| `idx_low` | Start index of accepted window in the original vectors |
| `idx_high` | End index of accepted window in the original vectors |

### Quick start

```matlab
% Run the built-in test with a synthetic dataset
run test_fault_ks.m
```

---

## Parameter Guide

### `min_frac`
Controls the minimum window size as a fraction of total data:
- **Increase** (e.g. `0.40`) to force broader fits on large, well-sampled datasets
- **Decrease** (e.g. `0.15`) if the valid power-law range is genuinely narrow
- Default `0.25` ensures at least a quarter of the data is always used, giving the K-S test sufficient statistical power

### `alpha`
Sets how strictly "linear" a window must be to be accepted:
- `0.05` — standard threshold; rejects windows with clear non-linearity *(default)*
- `0.10` — stricter; tends to select a narrower but cleaner range
- `0.01` — very lenient; almost any window passes, selection dominated by the longest-window rule

> **Note:** The K-S test is used *inversely* here. A **high** p-value is desirable — it means there is no detectable departure from linearity, so the straight-line model is retained with confidence.

---

## Synthetic Test Dataset

`test_fault_ks.m` generates 59 points across three zones:

| Zone | x range | Slope | Represents |
|------|---------|-------|------------|
| Low (14 pts) | 0.10 – 0.95 | ≈ −0.4 | Under-sampled small faults |
| Mid (35 pts) | 1.00 – 3.00 | −2.0 (true) | Valid power-law range |
| High (10 pts) | 3.05 – 3.80 | ≈ −5.5 | Truncated large faults |

The script produces a three-panel diagnostic figure:
- **Top:** log-log scatter with shaded excluded zones, the true model, and the fitted line
- **Bottom-left:** residual scatter for the accepted window
- **Bottom-right:** Q-Q plot confirming Gaussian residuals within the window

---

## Requirements

- MATLAB R2019b or later
- **Statistics and Machine Learning Toolbox** (for `kstest` and `qqplot`)

---

## License

MIT License. See `LICENSE` for details.
