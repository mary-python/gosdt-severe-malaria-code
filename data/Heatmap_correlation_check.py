from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

# Load dataset
csv_path = Path("Imputed_dataset_long.csv")
df = pd.read_csv(csv_path, low_memory = False)

# Variables used in imputation workflow
cols = ["HCT", "Hb", "Age", "Respiratory_distress", "para_ul", "Temperature_C",
        "Coma", "Lactate", "Glucose_mmol_L", "Base_Excess_mmol_L", "Bicarbonate_mEq_L",
        "Creatinine_umol_L", "BUN", "Died_28D"]

# Keep only selected columns and make sure values are numeric where possible
sub = df[cols].apply(pd.to_numeric, errors = "coerce")

# Correlation of the actual variables
corr = sub.corr(method = "pearson")

# Mask the diagonal so self-correlations do not dominate the figure
plot_corr = corr.copy()
np.fill_diagonal(plot_corr.values, np.nan)

# Use only the visible off-diagonal values to define the scale
offdiag_values = plot_corr.to_numpy()
offdiag_values = offdiag_values[~np.isnan(offdiag_values)]
vmin = offdiag_values.min()
vmax = offdiag_values.max()

# Blue -> yellow colormap
cmap = LinearSegmentedColormap.from_list("blue_yellow",
       ["#081d58", "#225ea8", "#41b6c4", "#ffff00"])
cmap.set_bad(color = "lightgrey")

# Slightly emphasise stronger correlations
norm = PowerNorm(gamma = 0.8, vmin = vmin, vmax = vmax)

# Plot heatmap
plt.figure(figsize = (12, 10))
im = plt.imshow(plot_corr, aspect = "auto", cmap = cmap, norm = norm)
plt.colorbar(im, label = "Correlation")

plt.xticks(range(len(plot_corr.columns)), plot_corr.columns, rotation = 90)
plt.yticks(range(len(plot_corr.columns)), plot_corr.columns)
plt.title("Correlation Heatmap")
plt.tight_layout()
plt.savefig("Heatmap_correlation_check.png", dpi = 300, bbox_inches = "tight")
plt.show()