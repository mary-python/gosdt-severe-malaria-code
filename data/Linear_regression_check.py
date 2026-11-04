from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext, FuncFormatter, NullFormatter

# Load dataset
csv_path = Path("Imputed_dataset_long.csv")
df = pd.read_csv(csv_path, low_memory = False)

# Choose how to handle the imputations
imputation_mode = "first"   # "first" or "average"

# Column names for the imputation structure
imputation_col = ".imp"
id_col = ".id"

# Variables used in plotting
pairs = [("Hb", "HCT", "Hb vs HCT", False),
         ("Bicarbonate_mEq_L", "Base_Excess_mmol_L", "Bicarbonate vs Base Excess", False),
         ("BUN", "Creatinine_umol_L", "BUN vs Creatinine", True)]

# Display labels for axes
display_name = {"Hb": "Hb [g/dL]", "HCT": "HCT [L/L]",
                "Bicarbonate_mEq_L": "Bicarbonate [mEq/L]",
                "Base_Excess_mmol_L": "Base Excess [mmol/L]", "BUN": "BUN [mmol/L]",
                "Creatinine_umol_L": "Creatinine [umol/L]"}

# Label only exact powers of 10 as plain numbers
def plain_decade_label(val, pos = None):
    if val <= 0:
        return ""
    exp = np.log10(val)
    if np.isclose(exp, np.round(exp)):
        return f"{val:g}"
    return ""

# Automatically choose plain decades of powers-of-10-notation
def choose_log_tick_style(x, y):
    positive = np.concatenate([x[x > 0], y[y > 0]])
    if positive.size == 0:
        return "plain"

    decade_min = np.floor(np.log10(positive.min()))
    decade_max = np.ceil(np.log10(positive.max()))
    span = decade_max - decade_min

    if span <= 4:
        return "plain"
    return "scientific"

# Plot each regression in a separate figure
for xcol, ycol, title, use_log in pairs:

    # Keep only selected columns and make sure values are numeric where possible
    if imputation_mode == "average":

        # Average values across imputations for each participant
        sub = df[[id_col, xcol, ycol]].copy()
        sub[xcol] = pd.to_numeric(sub[xcol], errors = "coerce")
        sub[ycol] = pd.to_numeric(sub[ycol], errors = "coerce")
        sub = sub.dropna()
        sub = sub.groupby(id_col, as_index = False)[[xcol, ycol]].mean()

    else:

        # Keep only the first imputation
        sub = df[df[imputation_col] == 1][[xcol, ycol]]
        sub = sub.apply(pd.to_numeric, errors = "coerce").dropna()

    # Pull out x and y values
    x = sub[xcol].to_numpy()
    y = sub[ycol].to_numpy()

    # Display labels
    xlab = display_name.get(xcol, xcol)
    ylab = display_name.get(ycol, ycol)

    plt.figure(figsize = (8, 6))

    if use_log:

        # Keep only positive values for log-scale plots
        mask = (x > 0) & (y > 0)
        x = x[mask]
        y = y[mask]

        # Fit a line in log10 space
        logx = np.log10(x)
        logy = np.log10(y)
        slope, intercept = np.polyfit(logx, logy, 1)

        # Line for plotting on log axes
        x_line = np.logspace(np.log10(x.min()), np.log10(x.max()), 200)
        y_line = 10 ** (slope * np.log10(x_line) + intercept)

        # Calculate the correlation coefficient and R-squared in log space
        r = np.corrcoef(logx, logy)[0, 1]
        r2 = r ** 2

        # Plot scatter and regression line
        plt.scatter(x, y, alpha = 0.35, s = 18)
        plt.plot(x_line, y_line, linewidth = 2)

        ax = plt.gca()
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Automatically choose log tick formatting
        tick_style = choose_log_tick_style(x, y)

        ax.xaxis.set_major_locator(LogLocator(base = 10.0))
        ax.yaxis.set_major_locator(LogLocator(base = 10.0))
        ax.xaxis.set_minor_locator(LogLocator(base = 10.0, subs = np.arange(2, 10) * 0.1))
        ax.yaxis.set_minor_locator(LogLocator(base = 10.0, subs = np.arange(2, 10) * 0.1))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())

        if tick_style == "plain":
            ax.xaxis.set_major_formatter(FuncFormatter(plain_decade_label))
            ax.yaxis.set_major_formatter(FuncFormatter(plain_decade_label))
        else:
            ax.xaxis.set_major_formatter(LogFormatterMathtext(base = 10.0))
            ax.yaxis.set_major_formatter(LogFormatterMathtext(base = 10.0))

        # Add labels and annotation
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(f"After MI: {title}")
        plt.text(0.05, 0.95, f"n = {len(x)}\nR^2 = {r2:.3f} (log10 space)",
                 transform = plt.gca().transAxes, va = "top", ha = "left")

    else:

        # Fit a linear regression line
        slope, intercept = np.polyfit(x, y, 1)
        x_line = np.linspace(x.min(), x.max(), 200)
        y_line = slope * x_line + intercept

        # Calculate the correlation coefficient and R-squared
        r = np.corrcoef(x, y)[0, 1]
        r2 = r ** 2

        # Plot scatter and regression line
        plt.scatter(x, y, alpha = 0.35, s = 18)
        plt.plot(x_line, y_line, linewidth = 2)

        # Add labels and annotation
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(f"After MI: {title}")
        plt.text(0.05, 0.95, f"n = {len(x)}\nR^2 = {r2:.3f}",
                 transform = plt.gca().transAxes, va = "top", ha = "left")

    plt.tight_layout()
    plt.savefig(f"LR_{xcol}_vs_{ycol}_check.png", dpi = 300, bbox_inches = "tight")
    plt.show()