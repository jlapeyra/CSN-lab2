import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


# --- 1. Load table ---
# Read fixed-width formatted text, skipping the separator line
df = pd.read_csv("plot/tables/aic2.csv")

# --- 3. Set Language (Best Model) as index ---
df = df.set_index('Language (Best Model)')

# --- 4. Convert all numeric columns to floats (safe conversion) ---
df = df.apply(pd.to_numeric, errors='coerce')

# --- 5. Plot heatmap ---
plt.figure(figsize=(5, 5))

#df_log = np.log10(df.replace(0, np.nan))
sns.heatmap(df[['AIC']], vmin=np.float64(42402.2), vmax=np.float64(77715.3), annot=True, fmt=".1f", cbar_kws={'label': 'AIC'})
plt.title("AIC Heatmap")
plt.ylabel("Language (Best Model)")
plt.tight_layout()
plt.show() 
